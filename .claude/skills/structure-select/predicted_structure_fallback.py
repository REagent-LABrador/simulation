"""predicted_structure_fallback — fold a target that has NO experimental
structure and NO usable homolog, so the computed axis still returns a result.

WHEN THIS RUNS, AND WHEN IT MUST NOT
------------------------------------
This is the LAST resort of structure selection, reached only when every cheaper
source is exhausted: no experimental structure of the target, and no homolog
that `homolog_transfer.py` could build a `transferred_homolog_site` from. In that
state the alternative is nulling the whole computed axis, and the contract is
that a target must still get an answer. So we fold it with ESMFold and hand the
predicted structure to the ordinary pocket scan — with warnings, never silently.

It MUST NOT run when an experimental structure exists. A predicted structure is
strictly worse than a deposited one for pocket geometry (the whole tractability
axis is built on real structures), and scoring pockets on a model when a crystal
is available is a downgrade dressed as coverage. Structure selection calls this
only after returning empty from the experimental and homolog routes.

THE GPU DEPENDENCY IS ISOLATED HERE
-----------------------------------
This is the ONLY place the dossier pipeline touches a GPU proto-tool. Everything
else — fpocket, mdpocket, the geometry, disorder, cryptic, the ligand and
modality classifiers — is CPU or stdlib. `esmfold_predict` lives in
`cofold-check/predict.py`, runs in the `proto-env` Modal environment, and is
imported lazily so importing THIS module costs nothing and never needs a GPU.
See CLAUDE.md rule 4c for why the GPU tools sit off the default path.

THE GATE, AND THE CONTRACT THAT OVERRIDES IT
--------------------------------------------
ESMFold returns its own pTM. Its benchmark (14 complexes / 28 runs, recorded in
`predict.OBSERVATIONS`) found pTM >= 0.80 was 5 of 5 with zero false alarms — so
that is the confidence band, NOT a pass/fail switch. This function NEVER refuses
to return a structure on a low pTM. A low pTM means a louder warning and a lower
`fold_confidence`, not a null. The reason: "we could not crystallise it and the
model is unsure" is itself the finding a structure-less orphan is entitled to,
and nulling it would hide exactly the target that most needs flagging.

The one honest refusal is a tool FAILURE — ESMFold erroring, timing out, or
proto-tools being absent from the environment. That is `not_run`, distinct from
a low-confidence result, and it is the only case where `structure_cif` is None.

THE ORDERED-CORE LESSON
-----------------------
The ESMFold interface benchmark reproduced a catastrophic failure that turned
out to be an INPUT artifact: folding a full UniProt mature chain with disordered
termini gave 1 inter-chain contact where the ordered core gave 55. So the
default input is the UniProt canonical sequence, and `residue_range` lets the
caller pass the crystallographically/annotated ordered core when one is known.
When no range is given we fold the whole chain and say so in the warnings — a
disordered tail lowers pTM and drags the pocket scan onto flexible loops, and
the reader has to know that is why.
"""

from __future__ import annotations

import urllib.error
import urllib.request
from typing import Any

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{}.fasta"

# The confidence bands, from predict.OBSERVATIONS["il17a_esmfold_dimer"] and the
# 28-run benchmark. These set the WARNING strength and `fold_confidence`; they
# never gate the return. Declared here, not fitted to this pipeline's targets.
PTM_USABLE = 0.80      # 5/5 with zero false alarms in the benchmark
PTM_MARGINAL = 0.55    # below this the fold topology itself is in doubt


def _fetch_uniprot_sequence(accession: str, *, timeout: int = 60) -> str:
    """The canonical sequence for a structure-less target, from UniProt.

    A structure-less target has, by definition, no PDB entry to lift a sequence
    out of, so the sequence comes from UniProt directly. Raises rather than
    returning empty: a fold with no input is not a low-confidence result, it is
    a failed run, and the caller must be able to tell them apart.
    """
    url = UNIPROT_FASTA.format(accession.strip())
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:  # noqa: S310
            text = r.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        raise RuntimeError(
            f"{accession}: UniProt returned HTTP {exc.code} for the canonical "
            f"sequence; cannot fold a target with no input sequence."
        ) from None
    except Exception as exc:  # noqa: BLE001
        raise RuntimeError(
            f"{accession}: could not reach UniProt for the canonical sequence "
            f"({type(exc).__name__}: {exc})."
        ) from None
    seq = "".join(
        line.strip() for line in text.splitlines() if line and not line.startswith(">")
    )
    if not seq:
        raise RuntimeError(
            f"{accession}: UniProt returned an empty sequence body — the "
            f"accession may be a deleted or demerged entry."
        )
    return seq


def _band(ptm: float | None) -> tuple[str, list[str]]:
    """Map pTM to a confidence label and the warnings that band mandates."""
    warnings: list[str] = [
        "PREDICTED STRUCTURE — this pocket scan ran on an ESMFold model, not a "
        "deposited experimental structure. Every pocket is a hypothesis about a "
        "computational model, one level less trustworthy than the rest of the "
        "computed axis.",
    ]
    if ptm is None:
        return "unknown", warnings + [
            "ESMFold did not report a pTM for this run; treat the fold as "
            "unvalidated."
        ]
    if ptm >= PTM_USABLE:
        return "usable_low_confidence", warnings + [
            f"ESMFold pTM {ptm:.3f} >= {PTM_USABLE:.2f} — the fold is in the "
            f"band that was 5/5 with zero false alarms in the benchmark. The "
            f"pockets are still model pockets; report them as such."
        ]
    if ptm >= PTM_MARGINAL:
        return "marginal", warnings + [
            f"ESMFold pTM {ptm:.3f} is below the {PTM_USABLE:.2f} usable band — "
            f"the fold is MARGINAL. Global topology may be right but domain "
            f"placement and loop geometry are unreliable, and pockets that "
            f"depend on inter-domain packing may be artifacts. Returned so the "
            f"target is not nulled; do not treat a pocket here as evidence of "
            f"tractability on its own."
        ]
    return "unreliable", warnings + [
        f"ESMFold pTM {ptm:.3f} is below {PTM_MARGINAL:.2f} — the fold itself is "
        f"UNRELIABLE and the topology may be wrong. This result is returned ONLY "
        f"to honour the contract that a structure-less target still gets an "
        f"answer; the correct reading is 'no usable structure by any route', "
        f"and the pocket numbers should carry no weight."
    ]


def predict_structure_for_pocket_scan(
    accession: str,
    *,
    sequence: str | None = None,
    residue_range: tuple[int, int] | None = None,
    seed: int = 42,
    num_recycles: int = 4,
    timeout: int | None = None,
) -> dict[str, Any]:
    """Fold a structure-less target and package it for `pocket_scan`.

    Returns a dict that ALWAYS carries `warnings` and `fold_confidence`, and
    carries `structure_cif` unless the tool itself failed. On success the caller
    passes the returned `predicted_structures` mapping straight into
    `pocket_scan(predicted_structures=..., uniprot_accession=accession)`.

    Never raises for a low-confidence fold. Raises only for a missing input
    sequence (no UniProt) — a genuinely failed run, reported as `status:
    not_run` by the caller wrapper below, never as a zero.

    Keys:
      status            "ok" | "not_run"
      structure_cif     the ESMFold mmCIF text, or None if not_run
      label             synthetic entry id used in the pocket-scan ensemble
      predicted_structures  {label: structure_cif} — the pocket_scan handoff
      ptm / avg_plddt / avg_pae   ESMFold's own confidence, unmodified
      fold_confidence   "usable_low_confidence" | "marginal" | "unreliable"
                        | "unknown"
      warnings          list[str], never empty on an ok result
      input_note        what sequence was folded and whether it was trimmed
      provenance        pass-through of the esmfold_predict provenance
    """
    # Lazy import: this is the only GPU touchpoint, and importing this module
    # must not require proto-tools or a Modal GPU. `predict` lives one skill
    # over, in cofold-check; the caller's sys.path must include it, exactly as
    # the deployed bundle lays the skills out side by side.
    from predict import esmfold_predict  # noqa: PLC0415

    seq = sequence or _fetch_uniprot_sequence(accession, timeout=timeout or 60)
    input_note = f"UniProt canonical sequence for {accession}, {len(seq)} residues"
    if residue_range is not None:
        lo, hi = residue_range
        seq = seq[lo - 1 : hi]  # 1-indexed inclusive, the UniProt convention
        input_note = (
            f"ordered core {lo}-{hi} of {accession} ({len(seq)} residues) — "
            f"folded the core, not the full chain, per the ESMFold ordered-core "
            f"lesson"
        )
    else:
        input_note += (
            " — FULL chain folded (no ordered-core range supplied). A "
            "disordered tail lowers pTM and can drag pockets onto flexible "
            "loops; pass residue_range when the ordered core is known."
        )

    label = f"ESMFOLD_{accession.strip().upper()}"

    try:
        pred = esmfold_predict(
            seq, seed=seed, num_recycles=num_recycles, timeout=timeout
        )
    except Exception as exc:  # noqa: BLE001
        # A tool failure is the ONE honest null: distinct from a low-confidence
        # fold, and the only case with no structure_cif. The computed axis is
        # then reported `not_established` with this reason, never as a zero.
        return {
            "status": "not_run",
            "structure_cif": None,
            "label": label,
            "predicted_structures": {},
            "ptm": None,
            "avg_plddt": None,
            "avg_pae": None,
            "fold_confidence": "not_run",
            "warnings": [
                f"ESMFold did not run: {type(exc).__name__}: {exc}. The computed "
                f"axis has no structure by any route; report site_hypothesis_"
                f"basis as not_established with this reason — do NOT return a "
                f"zero druggability."
            ],
            "input_note": input_note,
            "provenance": None,
        }

    self_report = pred.get("self_report") or {}
    ptm = self_report.get("ptm")
    confidence, warnings = _band(ptm)
    return {
        "status": "ok",
        "structure_cif": pred.get("structure_cif"),
        "label": label,
        "predicted_structures": {label: pred.get("structure_cif")},
        "ptm": ptm,
        "avg_plddt": self_report.get("avg_plddt"),
        "avg_pae": self_report.get("avg_pae"),
        "fold_confidence": confidence,
        "warnings": warnings,
        "input_note": input_note,
        "provenance": pred.get("provenance"),
    }
