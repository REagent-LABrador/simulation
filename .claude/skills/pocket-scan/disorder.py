"""Intrinsic-disorder prediction for the druggability pipeline.

Fills the `disorder_fraction` output field, which nothing previously populated.

WHY THIS FILE EXISTS SEPARATELY FROM modal_app.py
--------------------------------------------------
Disorder prediction has a fragile dependency (metapredict needs a C toolchain,
see IMAGE NOTES below). Keeping it in its own module means an import failure
degrades to a documented fallback instead of taking the whole Modal app down,
which is exactly what happened when a bare `metapredict` pin was added to the
image.

THE CARDINAL RULE
-----------------
A folded protein and a failed prediction must NEVER look identical.

`disorder_fraction == 0.0` is a real, meaningful answer: CDK2 and KRAS both
genuinely score 0.000. So failure can never be encoded as 0.0. Every public
entry point here returns **None** on failure and logs at ERROR level. Callers
must branch on `result is None`, not on a falsy fraction.

    res = predict_disorder(accession="P24941")
    if res is None:
        ...            # prediction FAILED - propagate null, do not score
    elif res["disorder_fraction"] < 0.1:
        ...            # protein is genuinely FOLDED

METHOD SELECTION
----------------
Tried in descending order of quality; the first that works wins, and the
`method` / `source` fields always record which one actually produced the number.

  1. metapredict  (MIT, CPU-only, ~155k residues/s, local)
     Best separation by a wide margin - MYC 0.828 vs folded controls 0.000-0.015.

  2. MobiDB REST API  (`prediction-disorder-th_50`)
     Network, accession-only, no install. Separation is real but ~3x weaker
     than metapredict (MYC 0.471 vs CA2 0.150), so it is a fallback and is
     flagged as reduced confidence.

  3. AlphaFold pLDDT via MobiDB  (`prediction-plddt-alphafold`)
     Last resort. See the INVERTED FIELD warning in _mobidb_fetch.

IUPred3 and ANCHOR2 are deliberately absent: their licence forbids commercial
use, which rules them out for this pipeline regardless of accuracy.

There is intentionally NO silent sequence-composition heuristic. A composition
score is available only by explicitly passing method="composition", and it
reports confidence="low" so it can never be mistaken for a real prediction.

IMAGE NOTES (Modal)
-------------------
metapredict publishes NO manylinux wheel - PyPI has only a cp38 macOS arm64
wheel and an sdist. On Linux, pip is therefore forced to build the sdist, which
compiles a Cython extension (metapredict/backend/cython/domain_definition.pyx)
and needs a C compiler. `modal.Image.micromamba()` is debian-slim based and has
no gcc, so the build dies with:

    error: command 'gcc' failed: No such file or directory

That is the whole failure. It is a missing build toolchain - NOT a torch pin,
NOT a numpy 2.x incompatibility (numpy 2.5.2 works fine), and NOT a
python-version ceiling (metapredict declares requires-python >=3.8).

To add it to the image, BOTH lines are required:

    .apt_install("git", "curl", "openjdk-17-jre-headless", "build-essential")
    .pip_install("torch", index_url="https://download.pytorch.org/whl/cpu")
    .pip_install("metapredict==3.0.2")

The CPU torch index is not optional. metapredict depends on torch, and on
Linux x86_64 plain `pip install torch` pulls a 527 MB wheel plus
nvidia-cudnn-cu13, nvidia-nccl-cu13, nvidia-cusparselt-cu13, nvidia-nvshmem-cu13
and triton - several GB of CUDA payload for a pipeline that never touches a GPU.
The CPU wheel is 184 MB and pulls none of it.
"""

from __future__ import annotations

import json
import logging
import re
import urllib.error
import urllib.request
from typing import Any, Iterable, Sequence

logger = logging.getLogger(__name__)

__all__ = [
    "predict_disorder",
    "disorder_fraction",
    "metapredict_available",
    "DisorderError",
]

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

METHOD_METAPREDICT = "metapredict"
METHOD_MOBIDB = "mobidb-th_50"
METHOD_PLDDT = "alphafold-plddt"
METHOD_COMPOSITION = "composition"

#: Order in which methods are attempted when method="auto".
_AUTO_ORDER = (METHOD_METAPREDICT, METHOD_MOBIDB, METHOD_PLDDT)

#: Standard metapredict disorder call. Residues scoring above this are
#: disordered. 0.5 is the value metapredict itself is calibrated on.
DEFAULT_THRESHOLD = 0.5

#: Short excursions above threshold are noise, not IDRs. Spans shorter than
#: this are dropped from `disordered_regions` but still counted in the
#: fraction, which is a per-residue quantity.
MIN_REGION_LENGTH = 12

UNIPROT_FASTA_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
MOBIDB_URL = "https://mobidb.org/api/download?acc={acc}&format=json"

_HTTP_TIMEOUT = 30
_ACCESSION_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$|^[OPQ][0-9][A-Z0-9]{3}[0-9]$")
#: UniProt FASTA headers look like ">sp|P24941|CDK2_HUMAN ...". Recovering the
#: accession from the header is what lets a caller pass a bare FASTA blob and
#: still reach the network fallbacks when metapredict is not installed.
_FASTA_ACC_RE = re.compile(r"^>(?:sp|tr)\|([A-Z0-9]+)\|", re.IGNORECASE)
_VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZUO")

#: Residues enriched in intrinsically disordered regions, and those enriched in
#: folded cores. Used only by the explicitly-requested composition heuristic.
_DISORDER_PROMOTING = set("PESKQRGDAN")
_ORDER_PROMOTING = set("WFYILVCMHN")


class DisorderError(RuntimeError):
    """Raised for unrecoverable disorder-prediction problems.

    Only escapes when ``strict=True``; the default path logs and returns None.
    """


# --------------------------------------------------------------------------
# metapredict availability - probed once, never at import time in a way that
# can raise. A broken metapredict must not stop this module importing.
# --------------------------------------------------------------------------

_METAPREDICT: Any = None
_METAPREDICT_CHECKED = False
_METAPREDICT_ERROR: str | None = None


def metapredict_available() -> bool:
    """Return True if metapredict imports and is usable in this process."""
    return _load_metapredict() is not None


def _load_metapredict() -> Any:
    """Import metapredict lazily, caching both success and failure."""
    global _METAPREDICT, _METAPREDICT_CHECKED, _METAPREDICT_ERROR
    if _METAPREDICT_CHECKED:
        return _METAPREDICT
    _METAPREDICT_CHECKED = True
    try:
        import metapredict as _mp  # type: ignore

        # Touch the Cython extension explicitly. If the image built the pure
        # -Python parts but the compiled domain_definition module is missing,
        # the plain import still succeeds and only blows up later, mid-run.
        _mp.predict_disorder("MEEPQSDPSV")
        _METAPREDICT = _mp
        logger.info("metapredict %s available", getattr(_mp, "__version__", "?"))
    except Exception as exc:  # noqa: BLE001 - any failure means unusable
        _METAPREDICT_ERROR = f"{type(exc).__name__}: {exc}"
        logger.warning(
            "metapredict unavailable (%s); falling back to network methods. "
            "To enable it the Modal image needs build-essential plus a "
            "CPU-only torch - see IMAGE NOTES in disorder.py",
            _METAPREDICT_ERROR,
        )
    return _METAPREDICT


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def _looks_like_accession(text: str) -> bool:
    return bool(_ACCESSION_RE.match(text.strip().upper()))


def _accession_from_fasta(text: str) -> str | None:
    """Pull the UniProt accession out of a FASTA header, if there is one."""
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            match = _FASTA_ACC_RE.match(line)
            if match and _looks_like_accession(match.group(1)):
                return match.group(1).upper()
            break
    return None


def _clean_sequence(seq: str) -> str | None:
    """Strip FASTA headers/whitespace and validate the alphabet."""
    if not seq:
        return None
    lines = [ln.strip() for ln in seq.splitlines()]
    body = "".join(ln for ln in lines if ln and not ln.startswith(">"))
    body = re.sub(r"[\s\-\*]", "", body).upper()
    if not body:
        return None
    bad = set(body) - _VALID_AA
    if bad:
        logger.error("sequence contains non-amino-acid characters: %s", sorted(bad))
        return None
    return body


def _http_get(url: str) -> bytes | None:
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "druggability-dossier/1.0"})
        with urllib.request.urlopen(req, timeout=_HTTP_TIMEOUT) as resp:
            if resp.status != 200:
                logger.error("GET %s returned HTTP %s", url, resp.status)
                return None
            return resp.read()
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
        logger.error("GET %s failed: %s", url, exc)
        return None


def fetch_sequence(accession: str) -> str | None:
    """Fetch a canonical sequence from UniProt. Returns None on any failure."""
    acc = accession.strip().upper()
    raw = _http_get(UNIPROT_FASTA_URL.format(acc=acc))
    if raw is None:
        logger.error("could not fetch sequence for %s from UniProt", acc)
        return None
    seq = _clean_sequence(raw.decode("utf-8", "replace"))
    if not seq:
        logger.error("UniProt returned no usable sequence for %s", acc)
    return seq


def _spans_from_scores(
    scores: Sequence[float],
    threshold: float,
    min_length: int,
) -> list[dict[str, int]]:
    """Contiguous runs above `threshold`, as 1-based inclusive spans."""
    regions: list[dict[str, int]] = []
    start: int | None = None
    for i, val in enumerate(scores):
        if val >= threshold:
            if start is None:
                start = i
        elif start is not None:
            regions.append((start, i - 1))  # type: ignore[arg-type]
            start = None
    if start is not None:
        regions.append((start, len(scores) - 1))  # type: ignore[arg-type]
    return [
        {"start": s + 1, "end": e + 1, "length": e - s + 1}
        for s, e in regions  # type: ignore[misc]
        if (e - s + 1) >= min_length
    ]


def _normalise_regions(pairs: Iterable[Sequence[int]], length: int) -> list[dict[str, int]]:
    """metapredict emits 0-based half-open [start, end]; convert to 1-based
    inclusive so spans line up with UniProt residue numbering."""
    out: list[dict[str, int]] = []
    for pair in pairs:
        s, e = int(pair[0]), int(pair[1])
        s1 = max(1, s + 1)
        e1 = min(length, e)
        if e1 >= s1:
            out.append({"start": s1, "end": e1, "length": e1 - s1 + 1})
    return out


def _result(
    fraction: float,
    method: str,
    regions: list[dict[str, int]],
    source: str,
    length: int,
    confidence: str,
    accession: str | None,
    threshold: float,
    **extra: Any,
) -> dict[str, Any]:
    payload = {
        "disorder_fraction": round(float(fraction), 4),
        "method": method,
        "disordered_regions": regions,
        "source": source,
        "length": length,
        "confidence": confidence,
        "accession": accession,
        "threshold": threshold,
        "n_disordered_regions": len(regions),
        "longest_disordered_region": max((r["length"] for r in regions), default=0),
    }
    payload.update(extra)
    return payload


# --------------------------------------------------------------------------
# Method 1: metapredict (preferred)
# --------------------------------------------------------------------------


def _predict_metapredict(
    seq: str, accession: str | None, threshold: float, min_length: int
) -> dict[str, Any] | None:
    mp = _load_metapredict()
    if mp is None:
        return None
    try:
        scores = [float(x) for x in mp.predict_disorder(seq)]
    except Exception as exc:  # noqa: BLE001
        logger.error("metapredict scoring failed: %s", exc)
        return None
    if len(scores) != len(seq):
        logger.error(
            "metapredict returned %d scores for a %d-residue sequence",
            len(scores),
            len(seq),
        )
        return None

    fraction = sum(1 for s in scores if s > threshold) / len(scores)

    # Prefer metapredict's own IDR domain calls: they merge short gaps and drop
    # noise, so on MYC they correctly report a single 1-376 IDR and EXCLUDE the
    # folded bHLH-LZ at the C-terminus rather than smearing across it.
    regions: list[dict[str, int]] = []
    try:
        domains = mp.predict_disorder_domains(seq)
        regions = _normalise_regions(domains.disordered_domain_boundaries, len(seq))
    except Exception as exc:  # noqa: BLE001 - non-fatal, fall back to spans
        logger.warning("predict_disorder_domains failed (%s); deriving spans from scores", exc)
        regions = _spans_from_scores(scores, threshold, min_length)

    version = getattr(mp, "__version__", "unknown")
    return _result(
        fraction,
        METHOD_METAPREDICT,
        regions,
        f"metapredict=={version} (local, CPU)",
        len(seq),
        "high",
        accession,
        threshold,
        mean_score=round(sum(scores) / len(scores), 4),
    )


# --------------------------------------------------------------------------
# Method 2/3: MobiDB REST
# --------------------------------------------------------------------------


def _mobidb_fetch(accession: str) -> dict[str, Any] | None:
    raw = _http_get(MOBIDB_URL.format(acc=accession))
    if raw is None:
        return None
    try:
        data = json.loads(raw.decode("utf-8", "replace"))
    except json.JSONDecodeError as exc:
        logger.error("MobiDB returned non-JSON for %s: %s", accession, exc)
        return None
    # The endpoint sometimes wraps the record in a single-element list.
    if isinstance(data, list):
        data = data[0] if data else None
    if not isinstance(data, dict):
        logger.error("MobiDB payload for %s was not an object", accession)
        return None
    return data


def _regions_from_mobidb(entry: dict[str, Any], min_length: int) -> list[dict[str, int]]:
    """MobiDB `regions` are [start, end] pairs, already 1-based inclusive."""
    out: list[dict[str, int]] = []
    for pair in entry.get("regions") or []:
        try:
            s, e = int(pair[0]), int(pair[1])
        except (TypeError, ValueError, IndexError):
            continue
        if e >= s and (e - s + 1) >= min_length:
            out.append({"start": s, "end": e, "length": e - s + 1})
    return out


def _predict_mobidb(accession: str | None, length: int | None, min_length: int) -> dict[str, Any] | None:
    if not accession:
        logger.error("MobiDB lookup needs a UniProt accession; only a sequence was given")
        return None
    data = _mobidb_fetch(accession)
    if data is None:
        return None

    # TRAP: for FOLDED proteins this key is frequently ABSENT rather than
    # present-and-zero. A bare data["prediction-disorder-th_50"] raises
    # KeyError on exactly the folded controls the pipeline cares about, so the
    # two-step .get with a dict default is load-bearing, not defensive noise.
    entry = data.get("prediction-disorder-th_50", {}) or {}
    fraction = entry.get("content_fraction")
    if fraction is None:
        logger.error(
            "MobiDB has no prediction-disorder-th_50 content_fraction for %s; "
            "refusing to infer 0.0 from an absent key",
            accession,
        )
        return None

    seq_len = int(data.get("length") or length or 0)
    return _result(
        float(fraction),
        METHOD_MOBIDB,
        _regions_from_mobidb(entry, min_length),
        MOBIDB_URL.format(acc=accession),
        seq_len,
        # Real separation but ~3x weaker than metapredict, so never "high".
        "medium",
        accession,
        DEFAULT_THRESHOLD,
    )


def _predict_plddt(accession: str | None, length: int | None, min_length: int) -> dict[str, Any] | None:
    """AlphaFold pLDDT as a disorder proxy. Last resort."""
    if not accession:
        logger.error("pLDDT lookup needs a UniProt accession")
        return None
    data = _mobidb_fetch(accession)
    if data is None:
        return None

    entry = data.get("prediction-plddt-alphafold", {}) or {}
    high_conf = entry.get("content_fraction")
    if high_conf is None:
        logger.error("MobiDB has no prediction-plddt-alphafold for %s", accession)
        return None

    # INVERTED FIELD: despite the name, content_fraction here is the fraction
    # of HIGH-confidence (well-folded) residues, not the disordered fraction.
    # Verified: MYC 0.251, CA2 0.988 - i.e. the IDP scores LOW. Disorder is the
    # complement. Failing to invert would rank every folded protein as maximally
    # disordered and vice versa.
    fraction = 1.0 - float(high_conf)

    seq_len = int(data.get("length") or length or 0)
    return _result(
        fraction,
        METHOD_PLDDT,
        _regions_from_mobidb(entry, min_length),
        MOBIDB_URL.format(acc=accession) + " [prediction-plddt-alphafold, inverted]",
        seq_len,
        "low",
        accession,
        DEFAULT_THRESHOLD,
        plddt_high_confidence_fraction=round(float(high_conf), 4),
    )


# --------------------------------------------------------------------------
# Method 4: composition heuristic - opt-in only, never reached by "auto"
# --------------------------------------------------------------------------


def _predict_composition(
    seq: str, accession: str | None, min_length: int
) -> dict[str, Any] | None:
    """Documented last resort. Windowed disorder-promoting residue fraction.

    This is NOT a trained predictor. It is reported with confidence="low" and a
    distinct method name so a downstream reader can always tell it apart from a
    real prediction. It is never selected automatically.
    """
    if len(seq) < 5:
        logger.error("sequence too short (%d) for composition heuristic", len(seq))
        return None
    window = min(25, len(seq))
    half = window // 2
    scores: list[float] = []
    for i in range(len(seq)):
        lo, hi = max(0, i - half), min(len(seq), i + half + 1)
        chunk = seq[lo:hi]
        d = sum(1 for c in chunk if c in _DISORDER_PROMOTING)
        o = sum(1 for c in chunk if c in _ORDER_PROMOTING)
        scores.append((d - o) / len(chunk))
    # Empirical cut; this heuristic is coarse by construction.
    thr = 0.25
    fraction = sum(1 for s in scores if s > thr) / len(scores)
    return _result(
        fraction,
        METHOD_COMPOSITION,
        _spans_from_scores(scores, thr, min_length),
        "sequence-composition heuristic (NOT a trained predictor)",
        len(seq),
        "low",
        accession,
        thr,
    )


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------


def predict_disorder(
    sequence: str | None = None,
    accession: str | None = None,
    method: str = "auto",
    threshold: float = DEFAULT_THRESHOLD,
    min_region_length: int = MIN_REGION_LENGTH,
    strict: bool = False,
) -> dict[str, Any] | None:
    """Predict intrinsic disorder for a protein.

    Accepts a raw sequence, a FASTA blob, or a UniProt accession. A single
    positional string is auto-detected: anything matching the UniProt accession
    pattern is treated as an accession, otherwise as a sequence.

    Args:
        sequence: Protein sequence or FASTA text. May also be an accession.
        accession: UniProt accession. Required for the network methods.
        method: "auto" (default), or one of "metapredict", "mobidb-th_50",
            "alphafold-plddt", "composition".
        threshold: Per-residue disorder cutoff (metapredict path).
        min_region_length: Spans shorter than this are omitted from
            `disordered_regions`; they still count toward the fraction.
        strict: Raise DisorderError instead of returning None on failure.

    Returns:
        A dict with at least ``disorder_fraction`` (float 0-1), ``method``,
        ``disordered_regions`` (1-based inclusive start/end spans) and
        ``source``; plus ``length``, ``confidence``, ``accession``,
        ``threshold``, ``n_disordered_regions`` and
        ``longest_disordered_region``.

        **None if prediction failed.** 0.0 always means genuinely folded, never
        "something went wrong" - see THE CARDINAL RULE at the top of this file.
    """

    def _fail(msg: str) -> None:
        logger.error("disorder prediction failed: %s", msg)
        if strict:
            raise DisorderError(msg)
        return None

    # --- resolve inputs ---------------------------------------------------
    if sequence and not accession and _looks_like_accession(sequence):
        accession, sequence = sequence.strip().upper(), None
    if sequence and not accession:
        # A FASTA blob carries its own accession; recovering it keeps the
        # network fallbacks reachable when metapredict is not installed.
        accession = _accession_from_fasta(sequence)
    if accession:
        accession = accession.strip().upper()
    if not sequence and not accession:
        return _fail("neither a sequence nor an accession was supplied")

    seq: str | None = None
    if sequence:
        seq = _clean_sequence(sequence)
        if seq is None:
            return _fail("supplied sequence could not be parsed")

    needs_seq = method in (METHOD_METAPREDICT, METHOD_COMPOSITION) or (
        method == "auto" and metapredict_available()
    )
    if seq is None and needs_seq and accession:
        seq = fetch_sequence(accession)
        if seq is None and method in (METHOD_METAPREDICT, METHOD_COMPOSITION):
            return _fail(f"could not obtain a sequence for {accession}")

    # --- dispatch ---------------------------------------------------------
    if method == METHOD_METAPREDICT:
        res = _predict_metapredict(seq or "", accession, threshold, min_region_length)
        return res if res is not None else _fail("metapredict unavailable or errored")
    if method == METHOD_MOBIDB:
        res = _predict_mobidb(accession, len(seq) if seq else None, min_region_length)
        return res if res is not None else _fail(f"MobiDB lookup failed for {accession}")
    if method == METHOD_PLDDT:
        res = _predict_plddt(accession, len(seq) if seq else None, min_region_length)
        return res if res is not None else _fail(f"pLDDT lookup failed for {accession}")
    if method == METHOD_COMPOSITION:
        res = _predict_composition(seq or "", accession, min_region_length)
        return res if res is not None else _fail("composition heuristic failed")
    if method != "auto":
        return _fail(f"unknown method {method!r}")

    # --- auto: best available, recording which one actually answered ------
    attempted: list[str] = []
    for candidate in _AUTO_ORDER:
        if candidate == METHOD_METAPREDICT:
            if seq is None:
                continue
            res = _predict_metapredict(seq, accession, threshold, min_region_length)
        elif candidate == METHOD_MOBIDB:
            res = _predict_mobidb(accession, len(seq) if seq else None, min_region_length)
        else:
            res = _predict_plddt(accession, len(seq) if seq else None, min_region_length)
        attempted.append(candidate)
        if res is not None:
            if len(attempted) > 1:
                res["fallback_from"] = attempted[:-1]
            return res

    return _fail(
        "every method failed (tried: "
        + ", ".join(attempted or ["none"])
        + "). Install metapredict or check network access to mobidb.org"
    )


def disorder_fraction(
    sequence: str | None = None,
    accession: str | None = None,
    **kwargs: Any,
) -> float | None:
    """Convenience wrapper returning just the fraction, or None on failure.

    Prefer :func:`predict_disorder` - the full dict carries the provenance that
    makes a number auditable.
    """
    res = predict_disorder(sequence=sequence, accession=accession, **kwargs)
    return None if res is None else res["disorder_fraction"]


if __name__ == "__main__":  # pragma: no cover - manual smoke test
    import sys

    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    targets = sys.argv[1:] or ["P01106", "P24941", "P00918", "P01116", "P01375"]
    for target in targets:
        out = predict_disorder(target)
        if out is None:
            print(f"{target}: FAILED (None)")
        else:
            print(
                f"{target}: frac={out['disorder_fraction']:.4f} "
                f"method={out['method']} conf={out['confidence']} "
                f"regions={[(r['start'], r['end']) for r in out['disordered_regions']]}"
            )
