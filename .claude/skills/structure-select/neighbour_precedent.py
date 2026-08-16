"""Structural-neighbour precedent: Foldseek fold neighbours -> holo/apo precedent.

Fills the dossier's `structural_neighbour_precedent` axis.

The question this answers is *not* "what family is my target in" — that is
`family_precedent`, the weakest axis. It is "what other folds look like mine,
and has anyone ever put a drug-like small molecule into one of them". A jump
along fold similarity carries a testable hypothesis; a jump along family
membership carries almost nothing (TNF-alpha and IL-17A are both cytokines and
share nothing mechanically).

Run it under the proto-tools python (needs `proto_tools`), with `paperclip` on
PATH and `PAPERCLIP_API_KEY` reachable — either exported, or in a `.env` this
module can find via `env_file=`.

    from neighbour_precedent import neighbour_precedent
    r = neighbour_precedent("/path/6OIM_A.pdb", "P01116")

DOWNSTREAM: `homolog_transfer.py`, in this directory. This module says a fold
neighbour exists and carries a drug-like ligand; that one turns the pair into
COORDINATES — a `transferred_homolog_site` anchor for a target with no ligand
of its own — under three guards, none of which this module applies. It is the
only ligand-free anchor with a positive record (2 of 3 constructible against
1 of 6 for interface, 1 of 3 for symmetry and 0 of 4 for annotated function),
and it is also where the ways a fold neighbour can be the *wrong* neighbour
show up. Two of them are visible from here and neither is caught here:

  * a neighbour whose only ligand is a **cofactor** (NLRP3's sole constructible
    donor was NOD2 with ADP), which this module reports as a drug-like-holo
    negative but which still gets picked when it is the only candidate; and
  * a neighbour whose ligand is on the **right chain and the wrong domain**
    (7KRZ's bortezomib at auth 768-898 against a NACHT match at auth 506-721).
    The chain-attribution fix documented in SKILL.md is necessary and does not
    reach this — a chain is not a domain.

Nothing here needs to change for that; read it as: a positive from this module
is a candidate, not a site.

Four verified Foldseek gotchas are handled here; see SKILL.md for the
measurements. In short: `hit.evalue` is really the probability and
`hit.bit_score` is really the E-value (remote mode mislabels columns), a
TM-score is only obtainable via `mode='tmalign'` where it lands in `bit_score`,
`target_id` is a filename-plus-title blob, not an ID, and **`foldseek-search`
does not assemble a multi-chain input into a complex** — which is why an
assembly with more than one chain is routed to `foldseek-multimer-search`.

CORRECTION, and it matters because the earlier version of this note was wrong.
`foldseek-search` was believed to search only ONE chain of a multi-chain file,
on the evidence that IL-17A 8DYG assembly1 (188 residues) returned no
`query_end` above 95. That inference does not hold. The raw m8 shows **both**
query chains present (`job_A` 144 rows, `job_B` 139) — query residues are simply
numbered per protomer, in both searches. What `foldseek-search` actually does is
align every chain INDEPENDENTLY and concatenate the results with no complex
assignment, no complex TM-score, and (because the wrapper drops m8 column 1) no
way to tell the chains apart. So it yields a union of per-protomer
neighbourhoods, which is still not evidence about an interface — the practical
conclusion survives, the stated mechanism did not. See `SEARCH_PATH_CAVEAT`.

The load-bearing caveat is `_ENTRY_LEVEL_CAVEAT` below and it is *implemented*,
not just documented: an entry's ligands are attributed to the entry, so in a
multi-protein complex a ligand bound to the partner counts toward every protein
in it. Every holo count therefore comes back twice — an entry-level upper bound
and a single-protein-entry lower bound — with titles attached so a reader can
adjudicate the gap.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any

# --------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------

_ENTRY_LEVEL_CAVEAT = (
    "pdb_v.entry_ligands is keyed at ENTRY level and its entity_id is the "
    "LIGAND's nonpolymer entity, not the protein chain it touches. In a "
    "multi-protein complex a ligand bound to the PARTNER still counts toward "
    "the entry. Verified damage: RAN (P62826) showed 36/139 holo, but the "
    "leptomycin-class ligands bind exportin, not RAN. Read "
    "n_holo_entry_level as an UPPER BOUND and n_holo_single_protein_entries "
    "as the defensible floor; check holo_titles before believing either."
)

SEARCH_PATH_CAVEAT = (
    "foldseek-search aligns each chain of a multi-chain file INDEPENDENTLY and "
    "returns the union, with no complex assignment and no complex TM-score. "
    "foldseek-multimer-search additionally pairs query chains to target chains "
    "into complex matches and scores the complex. Measured on IL-17A 8DYG "
    "assembly1: single-chain 283 rows / 125 entries, multimer 863 rows / 174 "
    "entries with 135 of 137 surviving entries matching BOTH query chains and "
    "ranked by complex TM-score. On TNF-alpha (3 chains) multimer returned "
    "6,891 rows and 1,206 complex assignments of size 3 -- genuine trimer "
    "matches. NOTE: max(query_end) is NOT a test of which path ran; query "
    "residues are numbered per protomer in both (IL-17A: 95 either way). Check "
    "meta['query_chains'] instead. On BOTH calibration targets the multimer "
    "path did NOT change the precedent conclusion -- IL-17A returns the "
    "cystine-knot superfamily and TNF-alpha the TNF superfamily, with zero "
    "defensible small-molecule holo either way. It is used because on the next "
    "oligomer that agreement cannot be checked, not because it changed these."
)

#: The multimer m8 has **26 columns**, not the 12 the shared parser reads, and
#: the extra ones are the only thing that makes it a *multimer* result.
#: Verified on the shipped 1HSG fixture and on IL-17A 8DYG assembly1.
#: 1-indexed, so index i-1 into a tab-split row:
_MULTIMER_M8_COLUMNS: dict[int, str] = {
    1: "query_chain, e.g. 'job_A' -- DROPPED by the shared 12-column parser",
    2: "target_id (filename_chain + free title)",
    3: "sequence identity, PERCENT",
    4: "alignment_length",
    5: "mismatches",
    6: "gap_openings",
    7: "query_start",
    8: "query_end",
    9: "target_start",
    10: "target_end",
    11: "PROBABILITY -- the parser names this `evalue` (gotcha 1, unchanged)",
    12: "E-VALUE -- the parser names this `bit_score` (gotcha 1, unchanged)",
    13: "true bit score",
    14: "query length",
    15: "target length",
    16: "query aligned sequence",
    17: "target aligned sequence",
    18: "target C-alpha coordinates",
    19: "target sequence",
    20: "complexassignid -- groups the chain-pair rows of ONE complex match",
    21: "complex TM-score normalised by QUERY -- DROPPED by the parser",
    22: "complex TM-score normalised by TARGET -- DROPPED by the parser",
    23: "rotation matrix (superposition)",
    24: "translation vector (superposition)",
    25: "taxonomy id",
    26: "species name",
}

# Ligand identity is decided by `ligand_filter.classify_ligands`, not by a
# comp_id denylist plus a size window. That pairing cannot work and both of this
# module's false positives proved it: 4PHH's `2UK` (a GppNHp analog) and 4EC7's
# `L44` (a diacylglycerol) both sailed through a 250-1200 Da window because
# their comp_ids were not on any list. The old `EXCLUDED_LIGANDS` / `MW_MIN` /
# `MW_MAX` / `_druglike_pred()` are DELETED, not deprecated.
#
# The structural consequence, and the reason this is not a one-line swap: **the
# SQL predicate cannot express chemistry**. Candidate ligands are selected
# without any exclusion clause and classified afterwards, in Python.

#: Fold-not-sequence filter. Below this identity the hit cannot be explained by
#: sequence homology; above this alignment length it is a domain, not a motif.
MAX_SEQUENCE_IDENTITY = 0.30
MIN_ALIGNMENT_LENGTH = 120

#: Bump whenever the parsed-hit shape or parse_target_id changes.
#: 3: hits carry query_chain / complex_assign_id / tm_score_kind, and the cache
#: blob records which search path produced them.
_CACHE_VERSION = 3


# --------------------------------------------------------------------------
# Paperclip plumbing
# --------------------------------------------------------------------------


def _load_env(env_file: str | os.PathLike[str] | None) -> dict[str, str]:
    """Return an environment with the .env's KEY=VALUE pairs merged in.

    Equivalent to `set -a; . <env_file>; set +a` for the simple assignments
    paperclip needs. Does not evaluate shell syntax.
    """
    env = dict(os.environ)
    if env_file is None:
        return env
    path = Path(env_file)
    if not path.is_file():
        raise FileNotFoundError(f"env_file not found: {path}")
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        if key.startswith("export "):
            key = key[len("export "):].strip()
        value = value.strip().strip("'").strip('"')
        env[key] = value
    return env


def _parse_table(stdout: str) -> list[dict[str, str]]:
    """Parse paperclip's aligned pipe table into dicts.

    Column boundaries come from the `---+---` rule, so a `|` inside a value
    cannot split a row. Returns [] for an empty result set.
    """
    lines = [ln.rstrip("\n") for ln in stdout.splitlines()]
    rule_idx = next(
        (
            i
            for i, ln in enumerate(lines)
            if ln.strip() and set(ln.strip()) <= {"-", "+"} and "-" in ln
        ),
        None,
    )
    if rule_idx is None or rule_idx == 0:
        return []
    rule = lines[rule_idx]
    header = lines[rule_idx - 1]

    # Slice spans from the rule's '+' positions.
    bounds: list[tuple[int, int | None]] = []
    start = 0
    for i, ch in enumerate(rule):
        if ch == "+":
            bounds.append((start, i))
            start = i + 1
    bounds.append((start, None))

    def cut(line: str) -> list[str]:
        return [line[a:b].strip() for a, b in bounds]

    names = cut(header)
    rows: list[dict[str, str]] = []
    for ln in lines[rule_idx + 1:]:
        s = ln.strip()
        if not s:
            continue
        if s.startswith("(") and ("row" in s):  # "(12 rows, 8ms)" trailer
            break
        if s.startswith("[") and s.endswith("]"):  # "[30ms]" timing line
            continue
        rows.append(dict(zip(names, cut(ln), strict=False)))
    return rows


def _run_sql(
    sql: str,
    *,
    env: dict[str, str],
    paperclip: str = "paperclip",
    timeout: float = 120.0,
) -> list[dict[str, str]]:
    """Run one read-only query against Paperclip's protein database."""
    exe = shutil.which(paperclip, path=env.get("PATH", os.defpath)) or paperclip
    proc = subprocess.run(  # noqa: S603
        [exe, "sql", "-s", "proteins", " ".join(sql.split())],
        capture_output=True,
        text=True,
        env=env,
        timeout=timeout,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"paperclip sql failed (rc={proc.returncode}): "
            f"{(proc.stderr or proc.stdout).strip()[:800]}"
        )
    out = proc.stdout
    if "..." in out and "|" not in out.split("\n")[0]:
        # Long single-column values get truncated with a literal ellipsis.
        raise RuntimeError("paperclip truncated the result; narrow the columns")

    # PAPERCLIP FAILS WITH rc=0. Three of its failure signatures come back on
    # stdout with a zero exit status, and `_parse_table` finds no `---+---`
    # rule in any of them, so it returns [] — a failed query and an empty
    # result set become the same value. Measured: `3LKJ` returned
    # "vsh: cd: /papers/: Permission denied" three times in a row, rc=0, 17 ms,
    # and a TNF-alpha run turned 25 unresolved entries into a confident
    # "25 apo / 0 holo". Raise instead, so the caller's retry sees it.
    low = out.lower()
    for sig in (
        "permission denied",       # a shell error, for a SQL query
        "something went wrong",    # undocumented, no code, no detail
        "request timed out",       # not a statement-cost signal
    ):
        if sig in low:
            raise RuntimeError(
                f"paperclip sql failed with rc=0 and a failure signature "
                f"({sig!r}): {out.strip()[:300]}"
            )

    rows = _parse_table(out)
    if not rows and "row" not in low and "|" not in out:
        # No rule line, no "(N rows)" trailer, no table at all. Whatever this
        # is, it is not an empty result set.
        raise RuntimeError(
            f"paperclip sql returned no parseable table (rc=0): {out.strip()[:300]!r}"
        )
    return rows


def _run_sql_retry(
    sql: str,
    *,
    env: dict[str, str],
    paperclip: str = "paperclip",
    budgets: tuple[float, ...] = (90.0, 240.0),
) -> list[dict[str, str]] | None:
    """`_run_sql` with short-then-long retries. `None` means it never returned.

    Short first, long second — retrying a deterministically slow statement at
    the same generous budget just doubles the bill. `None` is deliberately a
    distinct value from `[]`: the caller must be able to tell "the query said
    there is nothing" from "the query never ran", which is the whole point of
    the rc=0 signature checks in `_run_sql`.
    """
    for budget in budgets:
        try:
            return _run_sql(sql, env=env, paperclip=paperclip, timeout=budget)
        except (RuntimeError, subprocess.TimeoutExpired):
            continue
    return None


def _sql_list(values: object) -> str:
    """Render a Python iterable of strings as a SQL literal list."""
    items = [str(v).replace("'", "''") for v in values]  # type: ignore[union-attr]
    return ", ".join(f"'{v}'" for v in items) if items else "''"


# --------------------------------------------------------------------------
# Foldseek
# --------------------------------------------------------------------------


def parse_target_id(target_id: str) -> tuple[str, str | None, str]:
    """Split a Foldseek `target_id` into (pdb_id, chain, title).

    GOTCHA 3: `target_id` is not an ID. It arrives as
    ``"7r0n-assembly1.cif.gz_A KRasG12C in complex with GDP and compound 2"``.
    The PDB code is the first four characters; the chain is the `_X` at the end
    of the filename token; everything after the first space is a free title.
    """
    head, _, title = target_id.partition(" ")
    pdb_id = head[:4].upper()
    chain = None
    if "_" in head:
        tail = head.rsplit("_", 1)[1]
        if tail and len(tail) <= 6:
            # Foldseek disambiguates repeated auth chain ids in an assembly with
            # a "-N" suffix ("A-2"); the auth id itself may be multi-character
            # ("CCC"). Strip the suffix to get back to the deposited chain id.
            chain = tail.split("-", 1)[0]
    return pdb_id, chain, title.strip()


def count_chains(structure_path: str | os.PathLike[str]) -> int:
    """Number of distinct polymer chains Foldseek will see in this file.

    Counted off the normalised PDB that proto-tools actually POSTs, not off the
    CIF's `_struct_asym` — an assembly CIF renames chains on the way through
    `Structure`, and the number that matters is the one the server receives.
    """
    from proto_tools.entities import Structure  # noqa: PLC0415

    # ATOM only. HETATM chains are ligands and waters, and counting them would
    # route a single-chain holo structure to the multimer search on the
    # strength of its bound ligand sitting in its own chain id.
    text = Structure.from_file(str(structure_path)).structure_pdb
    return len({ln[21] for ln in text.splitlines() if ln.startswith("ATOM")})


def _foldseek(
    structure_path: str | os.PathLike[str],
    *,
    mode: str,
    database: str,
    timeout_seconds: float,
) -> Any:
    from proto_tools.tools.structure_alignment.foldseek.foldseek_search import (  # noqa: PLC0415
        FoldseekSearchConfig,
        FoldseekSearchInput,
        run_foldseek_search,
    )

    # NOTE: `structure` takes a file path or raw PDB/CIF text. It does NOT
    # accept a PDB ID or a URL — resolve those to a file first.
    return run_foldseek_search(
        FoldseekSearchInput(structure=str(structure_path)),
        FoldseekSearchConfig(
            search_mode="remote",
            databases=[database],
            mode=mode,
            timeout_seconds=timeout_seconds,
        ),
    )


def _foldseek_multimer(
    structure_path: str | os.PathLike[str],
    *,
    mode: str,
    database: str,
    timeout_seconds: float,
) -> Any:
    """Run `foldseek-multimer-search`, remote.

    No Modal, no GPU, no MODAL_PROFILE: the tool declares `uses_gpu=False` and
    its own `local_execution_reason` says a remote worker "would only add a
    network hop". It is an HTTPS POST to search.foldseek.com/api/ticket with
    `mode='complex-{mode}'` — the same ticket API the single-chain path uses.
    """
    from proto_tools.tools.structure_alignment.foldseek.foldseek_multimer_search import (  # noqa: PLC0415
        FoldseekMultimerSearchConfig,
        FoldseekMultimerSearchInput,
        run_foldseek_multimer_search,
    )

    return run_foldseek_multimer_search(
        FoldseekMultimerSearchInput(structure=str(structure_path)),
        FoldseekMultimerSearchConfig(
            search_mode="remote",
            databases=[database],
            mode=mode,
            timeout_seconds=timeout_seconds,
        ),
    )


def _fetch_raw_m8(result_url: str, *, timeout: float = 300.0) -> list[list[str]]:
    """Re-download a finished Foldseek result and return its RAW m8 rows.

    This is not a second search — the ticket is already computed and the
    archive is a static download, so it costs one HTTP GET and no queue time.

    It exists because `FoldseekMultimerHit = FoldseekHit` and the shared
    `_parse_m8_text` reads 12 fixed positions, so the wrapper silently discards
    columns 1, 20, 21 and 22 — the query chain, the complex-assignment id and
    both complex TM-scores. Without them a multimer result is indistinguishable
    from a single-chain one that happened to return each entry twice.
    """
    import io  # noqa: PLC0415
    import tarfile  # noqa: PLC0415

    import requests  # noqa: PLC0415

    resp = requests.get(result_url, timeout=timeout)
    resp.raise_for_status()
    rows: list[list[str]] = []
    with tarfile.open(fileobj=io.BytesIO(resp.content), mode="r:gz") as tf:
        for member in tf.getmembers():
            if not member.name.endswith(".m8"):
                continue
            handle = tf.extractfile(member)
            if handle is None:
                continue
            for line in handle.read().decode().splitlines():
                cols = line.split("\t")
                if len(cols) >= 22:
                    rows.append(cols)
    return rows


def foldseek_multimer_neighbours(
    structure_path: str | os.PathLike[str],
    *,
    database: str = "pdb100",
    timeout_seconds: float = 900.0,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run `foldseek-multimer-search` and return (hits, meta), same hit shape.

    A multimer result is **one row per (query chain, target chain) pair**, not
    one row per neighbour. Rows sharing a `complexassignid` are the chain
    correspondences of a single complex match and carry a single complex
    TM-score. Verified on the 1HSG fixture: 1,830 rows, 915 groups, every group
    of size 2, no group spanning two target entries, one TM-score per group.

    Two consequences the single-chain path does not have:

    * `query_end` is still **per protomer** — IL-17A's multimer search tops out
      at 95 exactly as the single-chain one did. The measurement that exposed
      the single-chain limit is therefore *not* a test for it. Count the
      distinct values of `query_chain` instead.
    * no second `mode='tmalign'` search is needed. Column 21 is the complex
      TM-score, which is a better number than the chain TM-score anyway when
      the site is an interface — and it halves exposure to a server whose
      latency spans two orders of magnitude.
    """
    path = Path(structure_path)
    if not path.is_file():
        raise FileNotFoundError(f"query structure not found: {path}")

    res = _foldseek_multimer(
        path, mode="3diaa", database=database, timeout_seconds=timeout_seconds
    )
    meta: dict[str, Any] = {
        "database": database,
        "mode": "complex-3diaa",
        "search_path": "multimer",
        "num_rows": res.num_hits,
        "execution_time_s": res.execution_time,
        "result_url": getattr(res, "result_url", None),
        "ticket_id": getattr(res, "ticket_id", None),
        "warnings": list(res.warnings or []),
    }

    raw = _fetch_raw_m8(str(res.result_url)) if res.result_url else []
    if not raw:
        raise RuntimeError(
            "multimer search returned no parseable 26-column m8; refusing to "
            "fall through to the 12-column parse, which would silently drop "
            "the query chain and the complex TM-score"
        )

    hits: list[dict[str, Any]] = []
    for rank, c in enumerate(raw, start=1):
        pdb_id, chain, title = parse_target_id(c[1])
        hits.append(
            {
                "rank": rank,
                "target_id": c[1],
                "pdb_id": pdb_id,
                "chain": chain,
                "target_title": title,
                # Column 1. The whole point: 'job_A' vs 'job_B' is the proof
                # that more than one query chain reached the search.
                "query_chain": c[0],
                "sequence_identity": float(c[2]) / 100.0,
                "alignment_length": int(c[3]),
                # Gotcha 1 holds unchanged on /foldmulti: col 11 is the
                # probability, col 12 the E-value. Verified on both fixtures.
                "probability": float(c[10]),
                "evalue": float(c[11]),
                "complex_assign_id": int(c[19]),
                "tm_score": float(c[20]),
                "tm_score_kind": "complex_qtm",
                "complex_ttm": float(c[21]),
                "query_range": [int(c[6]), int(c[7])],
                "target_range": [int(c[8]), int(c[9])],
            }
        )
    meta["num_hits"] = len(hits)
    meta["query_chains"] = sorted({h["query_chain"] for h in hits})
    meta["n_complex_assignments"] = len({h["complex_assign_id"] for h in hits})
    return hits, meta


def foldseek_neighbours(
    structure_path: str | os.PathLike[str],
    *,
    database: str = "pdb100",
    with_tm_score: bool = True,
    timeout_seconds: float = 900.0,
    cache_path: str | os.PathLike[str] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run Foldseek and return (hits, meta) with the mislabelled fields renamed.

    GOTCHA 1: in remote mode the server emits a 17-column m8 while the parser
    reads the standard 12-column layout, so `hit.evalue` actually holds the
    Foldseek PROBABILITY (higher is better) and `hit.bit_score` holds the true
    E-VALUE (lower is better). They are re-labelled here as `probability` and
    `evalue`; the raw attributes are never sorted or thresholded on. Hits
    arrive best-first, so list order is the safe ranking.

    GOTCHA 2: there is no TM-score field. With `mode='tmalign'`, `bit_score`
    carries the TM-score — the only route to one from this tool, and it works
    only because of gotcha 1. That is a second search, joined on target_id.
    """
    path = Path(structure_path)
    if not path.is_file():
        raise FileNotFoundError(f"query structure not found: {path}")

    # The public server's latency is wildly variable — 4.4 s for a 169-residue
    # KRAS chain, 323 s for a 186-residue IL-17A dimer on the same day. Cache,
    # so re-filtering never costs another search.
    cache = Path(cache_path) if cache_path else None
    if cache is not None and cache.is_file():
        import json  # noqa: PLC0415

        blob = json.loads(cache.read_text())
        # The cache holds PARSED hits, so a change to parse_target_id silently
        # invalidates it. This bit us once: chains cached as "A-2" stopped
        # matching after suffix-stripping was added, and five neighbours
        # quietly fell back to entry-level accessions. Version the cache.
        if blob.get("cache_version") == _CACHE_VERSION:
            blob["meta"]["from_cache"] = str(cache)
            return blob["hits"], blob["meta"]

    res = _foldseek(
        path, mode="3diaa", database=database, timeout_seconds=timeout_seconds
    )
    meta: dict[str, Any] = {
        "database": database,
        "mode": "3diaa",
        "search_path": "single_chain",
        "num_hits": res.num_hits,
        "execution_time_s": res.execution_time,
        "result_url": getattr(res, "result_url", None),
        "warnings": list(res.warnings or []),
    }

    tm_by_target: dict[str, float] = {}
    if with_tm_score:
        try:
            tm = _foldseek(
                path,
                mode="tmalign",
                database=database,
                timeout_seconds=timeout_seconds,
            )
            # bit_score is the TM-score here — see gotcha 2.
            tm_by_target = {h.target_id: float(h.bit_score) for h in tm.hits}
            meta["tmalign_num_hits"] = tm.num_hits
            meta["tmalign_execution_time_s"] = tm.execution_time
        except Exception as exc:  # noqa: BLE001 - TM-score is optional enrichment
            meta["tm_score_error"] = f"{type(exc).__name__}: {exc}"

    hits: list[dict[str, Any]] = []
    for rank, h in enumerate(res.hits, start=1):
        pdb_id, chain, title = parse_target_id(h.target_id)
        hits.append(
            {
                "rank": rank,
                "target_id": h.target_id,
                "pdb_id": pdb_id,
                "chain": chain,
                "target_title": title,
                "sequence_identity": h.sequence_identity,
                "alignment_length": h.alignment_length,
                # Renamed, deliberately. See gotcha 1.
                "probability": h.evalue,
                "evalue": h.bit_score,
                "tm_score": tm_by_target.get(h.target_id),
                "tm_score_kind": "chain_tmalign",
                # One chain went in, so every row belongs to that chain. Named
                # to match the multimer path so downstream code needs no branch.
                "query_chain": None,
                "query_range": [h.query_start, h.query_end],
                "target_range": [h.target_start, h.target_end],
            }
        )
    if cache is not None:
        import json  # noqa: PLC0415

        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(
            json.dumps(
                {"cache_version": _CACHE_VERSION, "hits": hits, "meta": meta}
            )
        )
    return hits, meta


# --------------------------------------------------------------------------
# Precedent lookup
# --------------------------------------------------------------------------

#: EVERY component of the entry is a candidate. There is no SQL-side filter at
#: all — not even `comp_type='non-polymer'`, and that omission is deliberate on
#: two independent grounds.
#:
#: Correctness: `comp_type` does not separate cofactors from drugs anyway (GDP
#: and UDP are `'RNA linking'` while ATP, GTP, NAD, FAD, HEM are all
#: `'non-polymer'`), so filtering on it in SQL both misses components and
#: pretends to a discrimination it does not have. `ligand_filter` reads the CCD
#: `type` itself and buckets polymers, ions and solvents on chemistry.
#:
#: Performance, measured: on the same 25-entry IL-17A neighbourhood,
#: `... WHERE l.entry_id IN (...) GROUP BY 1` returns in **6 ms**, and adding
#: `AND l.comp_type='non-polymer'` makes the identical query **time out**
#: (>120 s, `[error] Request timed out`). The column is not usefully indexed.
#: Do not "tidy up" by adding the predicate back.
_CANDIDATE = "TRUE"


def _classify(comp_ids: list[str]) -> dict[str, Any]:
    """comp_id -> LigandVerdict, in one batched Paperclip round trip per 40."""
    from ligand_filter import classify_ligands  # noqa: PLC0415

    return classify_ligands(sorted({c for c in comp_ids if c})) if comp_ids else {}


def _druglike_from(
    verdicts: dict[str, Any], comp_ids: list[str]
) -> tuple[list[str], list[str], list[str], list[str]]:
    """Split a ligand list into (druglike, names, reasons, undetermined).

    `undetermined` carries the comp_ids whose CCD lookup FAILED, which is not
    the same as "not drug-like" and must never be counted as apo. A Paperclip
    timeout returning `unknown` for two components that classify fine is a
    measured failure, not a hypothetical one.

    This is `ligand_filter.holo_call`'s contract implemented inline — `is_holo`
    is `bool(druglike)`, `determined` is `not undetermined` — because this call
    site also needs the per-ligand rejection reason, which `holo_call` returns
    only inside `verdicts`. If the two ever disagree, `holo_call` is right.
    """
    dl, names, reasons, undet = [], [], [], []
    for c in comp_ids:
        v = verdicts.get(c.upper())
        if v is None:
            continue
        if "lookup_failed" in getattr(v, "flags", ()):
            undet.append(c)
            continue
        if v.is_druglike:
            dl.append(c)
            names.append((v.name or c)[:44])
        else:
            # Why it was dropped, not just that it was. On a 25-neighbour
            # sweep with one or two holo hits, the reason IS the finding.
            reasons.append(f"{c}={v.verdict}: {v.reason[:70]}")
    return dl, names, reasons, undet


def _entry_facts(
    hits: list[tuple[str, str | None]], *, env: dict[str, str], paperclip: str
) -> dict[str, dict[str, Any]]:
    """Per-entry facts, plus the accession of the chain Foldseek actually aligned.

    The chain match matters and is the same class of error as the ligand
    attribution caveat, one level up. 2XAC is "VEGFR1 in complex with VEGF-B";
    Foldseek matched IL-17A to **chain A = P49765 (VEGF-B)**, but the entry also
    contains chains C/X = P17948 (VEGFR1, a receptor tyrosine kinase). Taking
    entry accessions wholesale imports VEGFR1 — and with it 3HNG's genuine
    kinase inhibitor — into a cystine-knot cytokine's fold neighbourhood, which
    is exactly the wrong answer.
    """
    entry_ids = sorted({e for e, _ in hits})
    pairs = ", ".join(
        f"({_sql_list([e])}, {_sql_list([c])})" for e, c in hits if c
    )
    chain_cte = (
        f"""q AS (SELECT * FROM (VALUES {pairs}) AS v(entry_id, chain)),
    cm AS (
      SELECT q.entry_id,
             STRING_AGG(DISTINCT pe.uniprot_accession, ',') AS chain_acc
      FROM q JOIN pdb_v.polymer_entities pe
        ON pe.entry_id = q.entry_id
       AND pe.auth_asym_ids @> to_jsonb(q.chain::text)
      WHERE pe.uniprot_accession IS NOT NULL
      GROUP BY 1),"""
        if pairs
        else "cm AS (SELECT NULL::text AS entry_id, NULL::text AS chain_acc),"
    )
    # NOTE: candidate ligands come back UNFILTERED — comp_ids only, aggregated
    # so 25 entries stay well inside the 200-row cap. Chemistry is decided
    # below, in Python; the SQL cannot do it.
    sql = f"""
    WITH {chain_cte}
    e AS (SELECT unnest(ARRAY[{_sql_list(entry_ids)}]) AS entry_id),
    np AS (
      SELECT pe.entry_id,
             COUNT(DISTINCT pe.uniprot_accession) AS n_acc,
             COUNT(*) FILTER (WHERE pe.polymer_type='polypeptide(L)') AS n_poly,
             STRING_AGG(DISTINCT pe.uniprot_accession, ',') AS accs
      FROM pdb_v.polymer_entities pe
      WHERE pe.entry_id IN (SELECT entry_id FROM e)
      GROUP BY 1),
    lg AS (
      SELECT e.entry_id,
             COALESCE(LEFT(STRING_AGG(DISTINCT CASE WHEN {_CANDIDATE}
                      THEN l.comp_id END, ','), 400), '-') AS cands
      FROM e LEFT JOIN pdb_v.entry_ligands l ON l.entry_id = e.entry_id
      GROUP BY 1)
    SELECT e.entry_id,
           COALESCE(np.n_acc, 0) AS n_acc,
           COALESCE(np.n_poly, 0) AS n_poly,
           COALESCE(np.accs, '-') AS accs,
           COALESCE(lg.cands, '-') AS cands,
           COALESCE(cm.chain_acc, '-') AS chain_acc,
           COALESCE(LEFT(en.title, 78), '-') AS title
    FROM e
    LEFT JOIN np ON np.entry_id = e.entry_id
    LEFT JOIN lg ON lg.entry_id = e.entry_id
    LEFT JOIN cm ON cm.entry_id = e.entry_id
    LEFT JOIN pdb_v.entries en ON en.entry_id = e.entry_id
    """
    rows = _run_sql_retry(sql, env=env, paperclip=paperclip)
    if rows is None:
        # Never ran. Return {} so the caller's missing-entry reconciliation
        # marks every entry UNDETERMINED, rather than raising and losing the
        # Foldseek result entirely. The empty dict is safe ONLY because that
        # reconciliation exists — without it this is the false-negative path.
        return {}

    def _cands(cell: str) -> list[str]:
        # LEFT(..., 400) is OUR truncation and it leaves no ellipsis, so a
        # comp_id can be cut in half and become a plausible-looking fake id.
        # 400 chars is ~66 comp_ids and no PDB entry comes close, but drop the
        # last token when the cell is exactly at the limit rather than trust it.
        ids = _split(cell)
        if len(cell) >= 400 and ids:
            ids = ids[:-1]
        return ids

    per_entry = {r["entry_id"]: _cands(r.get("cands", "-")) for r in rows}
    verdicts = _classify([c for ids in per_entry.values() for c in ids])

    out: dict[str, dict[str, Any]] = {}
    for r in rows:
        cands = per_entry.get(r["entry_id"], [])
        dl, names, reasons, undet = _druglike_from(verdicts, cands)
        out[r["entry_id"]] = {
            "n_protein_accessions": int(r["n_acc"] or 0),
            "n_polypeptide_entities": int(r["n_poly"] or 0),
            "accessions": [] if r["accs"] == "-" else r["accs"].split(","),
            "chain_accessions": (
                [] if r.get("chain_acc", "-") == "-" else r["chain_acc"].split(",")
            ),
            "has_druglike_holo": bool(dl),
            "druglike_ligands": dl,
            # Names, not just comp_ids: a comp_id tells a reader nothing, and
            # `"5'-O-[(R)-hydroxy...]guanosine"` tells them it is a nucleotide.
            "druglike_ligand_names": " / ".join(names) if names else "-",
            # Why each candidate was REJECTED. Keeps the negative auditable.
            "rejected_ligands": reasons,
            # A failed CCD lookup is not an apo call. If this is non-empty the
            # entry's state is undetermined and must not be reported as apo.
            "undetermined_ligands": undet,
            "holo_determined": not undet,
            "title": r["title"],
        }
    return out


def _druglike_comp_ids_for_accessions(
    accessions: list[str], *, env: dict[str, str], paperclip: str, page: int = 190
) -> tuple[list[str], list[str]]:
    """(drug-like comp_ids, accessions whose ligand sweep FAILED).

    Paged with LIMIT/OFFSET because Paperclip caps at 200 rows and a
    `STRING_AGG` of a few hundred comp_ids would hit the ~880-character cell
    truncation instead — the aggregate trick does not defeat the cap.

    The second return value is load-bearing. A statement timeout here is a
    *lookup failure*, and an accession whose ligands were never read is not an
    accession with no drug-like ligands. Reporting the two the same way is the
    original false-negative bug wearing new clothes, so they are kept apart and
    the failures surface as `lookup_failed` on the accession block.
    """
    # ONE ACCESSION AT A TIME. Unioning 17 accessions into a single DISTINCT
    # over every entry's ligands times out (measured: the IL-17A single-chain
    # neighbourhood, which includes VEGF-A with 75+ entries). Per accession it
    # is a handful of fast statements, and the union is done here.
    seen: set[str] = set()
    failed: list[str] = []
    for acc in accessions:
        offset = 0
        while True:
            # `IN (SELECT ...)` on entry_ligands is the fast plan; a direct
            # JOIN of structures_by_accession to entry_ligands is NOT, and
            # times out on the same accession (P15692: 9 ms vs >120 s). Keep
            # the subquery.
            sql = f"""
            WITH s AS (SELECT DISTINCT st.entry_id
                       FROM pdb_v.structures_by_accession st
                       WHERE st.accession = {_sql_list([acc])})
            SELECT DISTINCT l.comp_id FROM pdb_v.entry_ligands l
            WHERE l.entry_id IN (SELECT entry_id FROM s) AND {_CANDIDATE}
            ORDER BY 1 LIMIT {page} OFFSET {offset}
            """
            # P67861 times out reproducibly, so short-then-long matters here:
            # at 2 x 300 s one bad accession would cost ten minutes of a run
            # that has sixteen good ones left to do.
            rows = _run_sql_retry(sql, env=env, paperclip=paperclip)
            if rows is None:
                failed.append(acc)
                break
            seen.update(r["comp_id"] for r in rows if r.get("comp_id"))
            if len(rows) < page:
                break
            offset += page
    verdicts = _classify(sorted(seen))
    dl = sorted(c for c in seen if (v := verdicts.get(c.upper())) and v.is_druglike)
    return dl, sorted(set(failed))


def _accession_precedent(
    accessions: list[str], *, env: dict[str, str], paperclip: str
) -> tuple[dict[str, dict[str, Any]], dict[str, list[dict[str, Any]]]]:
    """Aggregate holo precedent across every PDB entry of each accession.

    Returns (per-accession summary, per-accession sample holo titles).

    Both halves exist because of the entry-level attribution caveat: the summary
    reports an upper bound and a floor, and the titles are what lets a reader
    decide which is right.

    This one aggregates over EVERY entry of each accession — hundreds — so it
    cannot pull ligand rows client-side and classify them; the 200-row cap
    would bite immediately. Instead it makes two passes: page out the DISTINCT
    non-polymer comp_ids (a few hundred at most, deduplicated across all those
    entries), classify that set once, then run the aggregation with the
    resulting drug-like comp_ids as an explicit `IN` list. Same chemistry, one
    server-side aggregation, no row-cap exposure.
    """
    druglike_ids, failed_accs = _druglike_comp_ids_for_accessions(
        accessions, env=env, paperclip=paperclip
    )

    def _mark(block: dict[str, Any], acc: str) -> dict[str, Any]:
        if acc in failed_accs:
            block["lookup_failed"] = True
            block["holo_determined"] = False
            block["basis"] = (
                "the ligand sweep for this accession FAILED (statement "
                "timeout). Its holo count is UNDETERMINED and a zero here "
                "must not be read as 'no drug-like ligands'."
            )
        else:
            block.setdefault("holo_determined", True)
        return block

    # FAST PATH, and it is the common one: nothing in the neighbourhood
    # classified drug-like, so every holo count is zero by construction and the
    # only remaining fact is how many structures each accession has. The full
    # aggregation below joins entry_ligands across every entry of every
    # accession and TIMES OUT on the server when there is nothing to find —
    # measured on the 10-accession IL-17A multimer neighbourhood. Do not pay
    # for a query whose answer is already known.
    if not druglike_ids:
        sql = f"""
        WITH a AS (SELECT unnest(ARRAY[{_sql_list(accessions)}]) AS acc)
        SELECT a.acc, COUNT(DISTINCT st.entry_id) AS n_struct
        FROM a LEFT JOIN pdb_v.structures_by_accession st ON st.accession = a.acc
        GROUP BY 1 ORDER BY 1
        """
        # Retry, and tolerate a total failure: when every accession's ligand
        # sweep already failed, this structure-count query is the one thing
        # still worth having, but its absence must not crash the axis either.
        rows_n = _run_sql_retry(sql, env=env, paperclip=paperclip) or []
        summary = {
            r["acc"]: _mark(
                {
                    "n_structures": int(r["n_struct"] or 0),
                    "n_holo_entry_level": 0,
                    "n_holo_single_protein_entries": 0,
                    "attribution_ambiguous_holo": 0,
                    "ligands_single_protein_entries": [],
                    "ligands_complex_entries_UNATTRIBUTED": [],
                    "basis": (
                        "no component across any structure of this accession "
                        "classified as drug-like by ligand_filter"
                    ),
                },
                r["acc"],
            )
            for r in rows_n
        }
        # An accession whose sweep failed still needs a block saying so, even
        # when the count query returned nothing to hang it on.
        for acc in accessions:
            if acc not in summary:
                summary[acc] = _mark({"n_structures": None}, acc)
        return summary, {}

    dl = f"l.comp_id IN ({_sql_list(druglike_ids)})"

    def _base(acc: str) -> str:
        # ONE ACCESSION PER STATEMENT. Unioning the accessions into a single
        # `unnest` array makes `lg` scan every entry of every accession at
        # once, and that TIMES OUT on the server (measured: the 8-accession
        # TNF-superfamily neighbourhood, even at a 300 s client budget). Per
        # accession the same aggregation is bounded and the union is done here.
        return f"""
        WITH s AS (
          SELECT DISTINCT {_sql_list([acc])} AS acc, st.entry_id
          FROM pdb_v.structures_by_accession st
          WHERE st.accession = {_sql_list([acc])}),
        np AS (
          SELECT pe.entry_id,
                 COUNT(DISTINCT pe.uniprot_accession) AS n_acc,
                 COUNT(*) FILTER (WHERE pe.polymer_type='polypeptide(L)') AS n_poly
          FROM pdb_v.polymer_entities pe
          WHERE pe.entry_id IN (SELECT entry_id FROM s)
          GROUP BY 1),
        lg AS (
          SELECT s.acc, s.entry_id,
                 COALESCE(np.n_acc, 1) AS n_acc,
                 COALESCE(np.n_poly, 1) AS n_poly,
                 COUNT(l.comp_id) FILTER (WHERE {dl}) AS n_dl,
                 COALESCE(STRING_AGG(DISTINCT CASE WHEN {dl} THEN l.comp_id END, ','), '') AS dl
          FROM s
          LEFT JOIN pdb_v.entry_ligands l ON l.entry_id = s.entry_id
          LEFT JOIN np ON np.entry_id = s.entry_id
          GROUP BY 1, 2, 3, 4)
        """

    summary_tail = """
    SELECT acc,
           COUNT(*) AS n_struct,
           COUNT(*) FILTER (WHERE n_dl > 0) AS n_holo_entry,
           COUNT(*) FILTER (WHERE n_dl > 0 AND n_acc <= 1 AND n_poly <= 1) AS n_holo_single,
           COALESCE(LEFT(STRING_AGG(DISTINCT dl, ' ')
                    FILTER (WHERE n_dl > 0 AND n_acc <= 1 AND n_poly <= 1), 110), '-') AS lig_single,
           COALESCE(LEFT(STRING_AGG(DISTINCT dl, ' ')
                    FILTER (WHERE n_dl > 0 AND (n_acc > 1 OR n_poly > 1)), 110), '-') AS lig_complex
    FROM lg GROUP BY 1
    """

    titles_tail = """
    , ranked AS (
      SELECT lg.*, ROW_NUMBER() OVER (
        PARTITION BY acc ORDER BY (n_acc <= 1 AND n_poly <= 1) DESC, entry_id) AS rn
      FROM lg WHERE n_dl > 0)
    SELECT r.acc, r.entry_id, r.n_acc, r.n_poly, LEFT(r.dl, 24) AS dl,
           COALESCE(LEFT(e.title, 66), '-') AS title
    FROM ranked r LEFT JOIN pdb_v.entries e ON e.entry_id = r.entry_id
    WHERE r.rn <= 3 ORDER BY r.rn
    """

    def _try(sql: str) -> list[dict[str, str]] | None:
        return _run_sql_retry(sql, env=env, paperclip=paperclip)

    summary: dict[str, dict[str, Any]] = {}
    titles: dict[str, list[dict[str, Any]]] = {}
    for acc in accessions:
        rows_s = _try(_base(acc) + summary_tail)
        if rows_s is None:
            failed_accs.append(acc)
            summary[acc] = _mark({"n_structures": None}, acc)
            continue
        for r in rows_s:
            n_entry = int(r["n_holo_entry"] or 0)
            n_single = int(r["n_holo_single"] or 0)
            summary[r["acc"]] = _mark(
                {
                    "n_structures": int(r["n_struct"] or 0),
                    "n_holo_entry_level": n_entry,
                    "n_holo_single_protein_entries": n_single,
                    "attribution_ambiguous_holo": n_entry - n_single,
                    "ligands_single_protein_entries": _split(r["lig_single"]),
                    "ligands_complex_entries_UNATTRIBUTED": _split(r["lig_complex"]),
                },
                r["acc"],
            )
        # Titles only exist where something was found; skip the second, slower
        # statement for the (common) accession with no drug-like holo at all.
        if not summary.get(acc, {}).get("n_holo_entry_level"):
            continue
        for r in _try(_base(acc) + titles_tail) or []:
            n_acc, n_poly = int(r["n_acc"] or 1), int(r["n_poly"] or 1)
            titles.setdefault(r["acc"], []).append(
                {
                    "pdb_id": r["entry_id"],
                    "ligands": _split(r["dl"]),
                    "n_protein_accessions": n_acc,
                    "n_polypeptide_entities": n_poly,
                    "attribution": (
                        "unambiguous" if (n_acc <= 1 and n_poly <= 1)
                        else "ambiguous_multiprotein"
                    ),
                    "title": r["title"],
                }
            )
    return summary, titles


def _split(cell: str) -> list[str]:
    if not cell or cell == "-":
        return []
    parts: list[str] = []
    for chunk in cell.replace(",", " ").split():
        if chunk and chunk not in parts and not chunk.endswith("..."):
            parts.append(chunk)
    return parts


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def neighbour_precedent(
    structure_path: str | os.PathLike[str],
    accession: str,
    *,
    max_neighbours: int = 25,
    min_alignment_length: int | None = None,
    relax_if_fewer_than: int = 5,
    database: str = "pdb100",
    with_tm_score: bool = True,
    multimer: bool | str = "auto",
    cache_path: str | os.PathLike[str] | None = None,
    env_file: str | os.PathLike[str] | None = "/Users/bb/repos/claude-agent-starter/.env",
    paperclip: str = "paperclip",
    timeout_seconds: float = 900.0,
) -> dict[str, Any]:
    """Structural-neighbour precedent for one target.

    Args:
        structure_path: PDB/CIF file for the query chain or assembly. A path or
            raw text — **never** a PDB ID or URL.
        accession: the query's UniProt accession, excluded from its own results.
        max_neighbours: cap on distinct neighbour PDB entries carried forward.
        min_alignment_length: override the length floor. Leave `None` for the
            verified default of 120 with automatic relaxation on short queries —
            see below.
        relax_if_fewer_than: with an auto floor, relax when the strict filter
            leaves fewer than this many neighbours.
        multimer: `'auto'` (default) routes any input with more than one chain
            to `foldseek-multimer-search`, because `foldseek-search` reads only
            one chain and an oligomeric site is not a property of one chain.
            `True` forces it, `False` pins the single-chain path. On a multimer
            failure the single-chain path runs as a fallback and both the
            attempt and the error are recorded — a degraded answer that says so
            beats no answer, but it must never be reported as a multimer one.
            NOTE: `cache_path` applies to the SINGLE-CHAIN path only. A multimer
            search is not cached here; re-parse its `result_url` archive
            instead, which is a static download and costs no queue time.

    Returns a dict shaped for the dossier's `structural_neighbour_precedent`
    block, with two holo counts per neighbour accession (see
    `caveats.entry_level_attribution`) and titles to adjudicate between them.
    `search_path` says which of the two searches produced it — check it before
    writing anything about an interface.
    """
    env = _load_env(env_file)

    n_chains = count_chains(structure_path)
    want_multimer = n_chains > 1 if multimer == "auto" else bool(multimer)
    multimer_error: str | None = None
    hits: list[dict[str, Any]] = []
    meta: dict[str, Any] = {}

    if want_multimer:
        try:
            hits, meta = foldseek_multimer_neighbours(
                structure_path, database=database, timeout_seconds=timeout_seconds
            )
        except Exception as exc:  # noqa: BLE001 - fall back, but loudly
            multimer_error = f"{type(exc).__name__}: {exc}"
    if not hits:
        hits, meta = foldseek_neighbours(
            structure_path,
            database=database,
            with_tm_score=with_tm_score,
            timeout_seconds=timeout_seconds,
            cache_path=cache_path,
        )
    meta["query_chains_in_file"] = n_chains
    if multimer_error:
        meta["multimer_attempted_and_failed"] = multimer_error
    search_path = meta.get("search_path", "single_chain")

    def apply(min_len: int) -> list[dict[str, Any]]:
        # Fold, not sequence family: drop anything explicable by sequence
        # homology, and anything too short to be a domain-level match.
        #
        # One entry, one neighbour — but a multimer result is one row per
        # (query chain, target chain) pair, so the surviving rows of an entry
        # must be MERGED rather than deduplicated away. `chains` therefore
        # holds every target chain Foldseek actually aligned, and only those:
        # the 2XAC/VEGFR1 fix ("only the aligned chain is a fold neighbour")
        # is unchanged, it just now admits a set instead of a singleton.
        out: list[dict[str, Any]] = []
        by_entry: dict[str, dict[str, Any]] = {}
        for h in hits:
            if h["sequence_identity"] >= MAX_SEQUENCE_IDENTITY:
                continue
            if h["alignment_length"] < min_len:
                continue
            rec = by_entry.get(h["pdb_id"])
            if rec is None:
                rec = dict(h)
                rec["chains"] = []
                rec["query_chains"] = []
                rec["n_chain_pairs"] = 0
                by_entry[h["pdb_id"]] = rec
                out.append(rec)
            if h["chain"] and h["chain"] not in rec["chains"]:
                rec["chains"].append(h["chain"])
            if h["query_chain"] and h["query_chain"] not in rec["query_chains"]:
                rec["query_chains"].append(h["query_chain"])
            rec["n_chain_pairs"] += 1
            rec["alignment_length"] = max(
                rec["alignment_length"], h["alignment_length"]
            )
            if h.get("tm_score") is not None:
                rec["tm_score"] = max(rec.get("tm_score") or 0.0, h["tm_score"])
        return out

    strict_len = (
        MIN_ALIGNMENT_LENGTH if min_alignment_length is None else min_alignment_length
    )
    kept = apply(strict_len)
    used_len, relaxed_note = strict_len, None

    # The 120-residue floor is calibrated on a ~170-residue GTPase domain. A
    # small cytokine cannot reach it: IL-17A's 8DYG assembly resolves 93
    # residues per protomer, and the strict filter left 2 neighbours out of 283
    # hits. Relax to 70% coverage of the query rather than report a starved set
    # as if it were the whole neighbourhood. Both counts are reported.
    #
    # The multimer path does NOT break this, and the reason is worth knowing:
    # multimer `query_end` is numbered within a protomer, not across the
    # assembly, so query_span stays ~95 for IL-17A whether one chain or two
    # were searched, and the relaxed floor comes out at the same 67. A longer
    # query span would have raised the floor and starved the result; it does
    # not happen. (It is also why max(query_end) is NOT a test for whether the
    # multimer path ran — count meta['query_chains'] instead.)
    if min_alignment_length is None and len(kept) < relax_if_fewer_than and hits:
        query_span = max(h["query_range"][1] for h in hits)
        relaxed_len = max(60, round(0.7 * query_span))
        if relaxed_len < strict_len:
            relaxed = apply(relaxed_len)
            relaxed_note = {
                "reason": (
                    f"strict floor {strict_len} left {len(kept)} neighbours "
                    f"(< relax_if_fewer_than={relax_if_fewer_than}); query spans "
                    f"only {query_span} residues"
                ),
                "strict_min_alignment_length": strict_len,
                "n_passing_strict": len(kept),
                "query_span_residues": query_span,
                "relaxed_min_alignment_length": relaxed_len,
                "n_passing_relaxed": len(relaxed),
                "basis": "70% coverage of the aligned query span, floored at 60",
            }
            kept, used_len = relaxed, relaxed_len

    filtered = kept[:max_neighbours]
    if search_path == "multimer":
        method = (
            f"foldseek-multimer-search (Proto, remote, {database}, "
            f"mode=complex-3diaa; TM-scores are COMPLEX TM-scores read from "
            f"raw m8 column 21, not a second tmalign search)"
        )
    else:
        method = (
            f"foldseek-search (Proto, remote, {database}, mode=3diaa"
            f"{'; TM-scores from a second mode=tmalign search' if with_tm_score else ''})"
        )
    result: dict[str, Any] = {
        "method": method,
        # The single most important field for a reader deciding whether this
        # block is evidence about an interface. Never omit it.
        "search_path": search_path,
        "n_query_chains_in_file": n_chains,
        # NOT "how many chains were searched" — foldseek-search searches all of
        # them and the wrapper drops the column that would prove it, so on the
        # single-chain path this is genuinely unknowable and is reported as
        # None rather than guessed at 1. On the multimer path it is counted
        # from raw m8 column 1.
        "n_query_chains_observed_in_result": (
            len(meta.get("query_chains", [])) if search_path == "multimer" else None
        ),
        "chains_assembled_into_complexes": search_path == "multimer",
        "query_structure": str(structure_path),
        "query_accession": accession,
        "foldseek": meta,
        "filter": {
            "sequence_identity_lt": MAX_SEQUENCE_IDENTITY,
            "alignment_length_gte": used_len,
            "auto_relaxed": relaxed_note,
            "rationale": (
                "keeps structural-but-not-sequence neighbours so this is a FOLD "
                "signal, not a sequence family in disguise"
            ),
            "n_hits_total": len(hits),
            "n_hits_passing": len(kept),
            "n_carried": len(filtered),
        },
        "neighbours": [],
        "neighbour_accessions": {},
        "caveats": {
            "entry_level_attribution": _ENTRY_LEVEL_CAVEAT,
            "ligand_identity": (
                "holo is decided by ligand_filter.classify_ligands, from "
                "chemistry, not by a comp_id denylist plus a size window. That "
                "pairing produced BOTH of this module's historical false "
                "positives (4PHH 2UK, a GppNHp analog; 4EC7 L44, a "
                "diacylglycerol) and cannot be repaired by extending the list. "
                "Each neighbour carries rejected_ligands saying why a "
                "candidate was dropped."
            ),
            "undetermined_is_not_apo": (
                "a FAILED chemical-component lookup returns undetermined, not "
                "apo. Check neighbour_entry_summary.n_undetermined before "
                "reading a zero holo count as a finding — 'no drug-like holo "
                "among the neighbours' is a real result, and a lookup failure "
                "wearing that costume is the worst confusion on this axis."
            ),
            "foldseek_field_mislabelling": (
                "remote-mode hit.evalue is the probability and hit.bit_score is "
                "the E-value; renamed here, never sorted on raw. Verified to "
                "hold identically on the multimer /foldmulti endpoint."
            ),
            "search_path": SEARCH_PATH_CAVEAT,
            "single_chain_on_an_oligomer": (
                None
                if search_path == "multimer" or n_chains <= 1
                else f"WARNING: the query file has {n_chains} chains but the "
                "SINGLE-CHAIN search ran. Each chain was aligned "
                "independently and the results concatenated, so this is a "
                "union of per-protomer neighbourhoods with no complex "
                "assignment. It is NOT evidence about an interface site."
            ),
        },
        "not_found": [],
    }
    if not filtered:
        result["not_found"].append(
            "no structural-but-not-sequence neighbours passed the filter"
        )
        return result

    # EVERY aligned chain, not just the first. On a multimer hit the match is a
    # set of chain correspondences, and resolving only one of them would drop
    # half the neighbour. `_entry_facts` STRING_AGGs DISTINCT per entry, so a
    # list of pairs yields the union of the matched chains' accessions.
    entry_facts = _entry_facts(
        [(h["pdb_id"], c) for h in filtered for c in (h["chains"] or [h["chain"]]) if c],
        env=env,
        paperclip=paperclip,
    )

    # A MISSING ENTRY IS UNDETERMINED, NOT APO. Measured: a TNF-alpha run under
    # heavy concurrent load got zero rows back from _entry_facts with rc=0, and
    # every downstream `.get()` then defaulted so cleanly that the module
    # reported "25 apo / 0 holo, 0 undetermined, 0 accessions, 0 rejected
    # ligands" — a confident clean negative produced by a retrieval that never
    # returned. Reconcile the entries asked for against the entries answered,
    # and refuse to let the gap default to a finding.
    missing = [h["pdb_id"] for h in filtered if h["pdb_id"] not in entry_facts]
    if missing:
        result["not_found"].append(
            f"entry facts missing for {len(missing)} of {len(filtered)} "
            f"neighbour entries ({', '.join(missing[:8])}"
            f"{' ...' if len(missing) > 8 else ''}); their holo/apo state is "
            "UNDETERMINED and is not reported as apo"
        )

    acc_set: list[str] = []
    for h in filtered:
        f = entry_facts.get(h["pdb_id"], {})
        entry_accs = [a for a in f.get("accessions", []) if a != accession]
        # Only the ALIGNED chain is a fold neighbour. The rest of the entry are
        # its crystallisation partners and must not inherit the precedent.
        chain_accs = [a for a in f.get("chain_accessions", []) if a != accession]
        neighbour_accs = chain_accs or entry_accs
        result["neighbours"].append(
            {
                "pdb_id": h["pdb_id"],
                "chain": h["chain"],
                # Multimer: every target chain that matched, and how many query
                # chains reached it. A neighbour with one query chain in a
                # multimer run matched a protomer, not the assembly — that is a
                # weaker claim and the field is here so it can be seen.
                "chains": h["chains"],
                "n_query_chains_matched": len(h["query_chains"]) or 1,
                "n_chain_pairs": h["n_chain_pairs"],
                "rank": h["rank"],
                "tm_score": h["tm_score"],
                "tm_score_kind": h.get("tm_score_kind"),
                "probability": h["probability"],
                "evalue": h["evalue"],
                "sequence_identity": h["sequence_identity"],
                "alignment_length": h["alignment_length"],
                "accessions": neighbour_accs,
                "chain_accession_resolved": bool(chain_accs),
                "other_accessions_in_entry": [
                    a for a in entry_accs if a not in neighbour_accs
                ],
                # Entry-level: true for the ENTRY, not necessarily for this chain.
                "has_druglike_holo": f.get("has_druglike_holo"),
                "ligands": f.get("druglike_ligands", []),
                "ligand_names": f.get("druglike_ligand_names"),
                # Why each candidate ligand was rejected, from ligand_filter.
                "rejected_ligands": f.get("rejected_ligands", []),
                # `f` empty => the entry never came back => UNDETERMINED.
                # The default here must be False, not True: `.get(k, True)` on
                # an absent record is how a failed retrieval becomes a clean
                # apo call.
                "holo_determined": f.get("holo_determined", False) if f else False,
                "undetermined_ligands": f.get("undetermined_ligands", []),
                "attribution": (
                    "unambiguous"
                    if f.get("n_protein_accessions", 0) <= 1
                    and f.get("n_polypeptide_entities", 0) <= 1
                    else "ambiguous_multiprotein"
                ),
                "n_protein_accessions": f.get("n_protein_accessions"),
                "title": f.get("title") or h["target_title"][:78],
            }
        )
        for a in neighbour_accs:
            if a not in acc_set:
                acc_set.append(a)

    n_holo = sum(1 for n in result["neighbours"] if n["has_druglike_holo"])
    result["neighbour_entry_summary"] = {
        "n_neighbour_entries": len(result["neighbours"]),
        "n_holo_entry_level": n_holo,
        "n_apo": len(result["neighbours"]) - n_holo,
        "holo_pdb_ids": [
            n["pdb_id"] for n in result["neighbours"] if n["has_druglike_holo"]
        ],
        # Multimer only. On the single-chain path there are no complex
        # assignments, so "matched every query chain" has no meaning and the
        # honest value is None, not 25.
        "n_matching_all_query_chains": (
            sum(
                1
                for n in result["neighbours"]
                if n["n_query_chains_matched"] >= len(meta.get("query_chains") or [])
            )
            if search_path == "multimer"
            else None
        ),
        # A zero here is only a finding if n_undetermined is also zero.
        "n_undetermined": sum(
            1 for n in result["neighbours"] if not n["holo_determined"]
        ),
        "undetermined_pdb_ids": [
            n["pdb_id"] for n in result["neighbours"] if not n["holo_determined"]
        ],
        "holo_is_a_flag_not_a_finding": (
            "n_holo_entry_level counts ENTRIES, and ligand_filter decides "
            "chemistry but cannot decide WHICH protein a ligand touches. The "
            "surviving failure is attribution: 8I2G's O6F is a genuine 468 Da "
            "FSHR allosteric agonist, correctly classified druglike, bound to "
            "the RECEPTOR — while the chain Foldseek matched was the FSH "
            "cystine-knot hormone. Read the title and `attribution` before "
            "reporting any holo hit as precedent."
        ),
    }

    if acc_set:
        summary, titles = _accession_precedent(
            acc_set, env=env, paperclip=paperclip
        )
        for acc in acc_set:
            block = summary.get(acc, {})
            block["holo_titles"] = titles.get(acc, [])
            result["neighbour_accessions"][acc] = block
        result["accession_summary"] = {
            "n_accessions": len(acc_set),
            "n_with_holo_entry_level": sum(
                1
                for b in result["neighbour_accessions"].values()
                if b.get("n_holo_entry_level", 0) > 0
            ),
            "n_with_holo_single_protein": sum(
                1
                for b in result["neighbour_accessions"].values()
                if b.get("n_holo_single_protein_entries", 0) > 0
            ),
        }
    else:
        result["not_found"].append(
            "no UniProt accessions mapped to the neighbour PDB entries"
        )

    tool_name = (
        "foldseek-multimer-search" if search_path == "multimer" else "foldseek-search"
    )
    result["sources"] = [
        f"{tool_name} remote {database} ({meta['num_hits']} hits, "
        f"{meta.get('execution_time_s')} s)",
        "paperclip pdb_v.polymer_entities / pdb_v.entry_ligands / pdb_v.entries "
        "/ pdb_v.structures_by_accession",
    ]
    if meta.get("result_url"):
        result["sources"].append(str(meta["result_url"]))
    return result


if __name__ == "__main__":
    import argparse
    import json

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("structure")
    ap.add_argument("accession")
    ap.add_argument("--max-neighbours", type=int, default=25)
    ap.add_argument("--min-alignment-length", type=int, default=None)
    ap.add_argument("--cache", default=None)
    ap.add_argument("--no-tm", action="store_true")
    ap.add_argument(
        "--multimer",
        choices=["auto", "yes", "no"],
        default="auto",
        help="auto (default) = multimer search when the file has >1 chain",
    )
    ap.add_argument("--env-file", default="/Users/bb/repos/claude-agent-starter/.env")
    a = ap.parse_args()
    print(
        json.dumps(
            neighbour_precedent(
                a.structure,
                a.accession,
                max_neighbours=a.max_neighbours,
                min_alignment_length=a.min_alignment_length,
                with_tm_score=not a.no_tm,
                multimer={"auto": "auto", "yes": True, "no": False}[a.multimer],
                cache_path=a.cache,
                env_file=a.env_file,
            ),
            indent=2,
        )
    )
