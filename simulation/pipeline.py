"""Produce a druggability dossier LOCALLY, dependency-light, from a bundled cache.

DEFAULT RUN — no cloud deploy, no Modal, no Paperclip, no managed agent, no API
key. The default ``run_pipeline`` needs only Python's stdlib + ``jsonschema``. It
resolves the request against a BUNDLED CACHE of REAL dossiers shipped inside this
module at ``simulation/cache/``, keyed on the documented tuple

    (uniprot_accession, mechanism_hypothesis, as_of_date)

Every number in a cache hit is a real value the deployed agent produced and that
was copied into the cache verbatim; the resolver NEVER recomputes or fabricates a
scientific value. On a cache MISS it returns a schema-valid but honestly EMPTY
dossier (``verdict: insufficient_evidence``, ``verdict_basis: none``, every axis
null-with-reason in ``not_found``, and a ``NOT_COMPUTED_LOCALLY`` limitation in
the interpretability block). It never invents numbers to fill the space.

Cache-key match:
  * EXACT — accession, (normalised) mechanism_hypothesis and as_of_date all agree.
  * ACCESSION-ONLY FALLBACK — only the accession agrees; the returned dossier is
    stamped with a ``CACHE_KEY_FALLBACK`` limitation because its mechanism/date
    differ from what was asked. A null/absent ``mechanism_hypothesis`` normalises
    to ``"unknown"`` on both sides (input.schema.json: absent == unknown).

Cache hits keep their own authored ``interpretability`` block; the resolver only
STAMPS ``interpretability.extensions`` with ``runtime_maturity: "LOCAL"`` and a
``CACHED_DOSSIER`` qualifier (the values the orchestrator/UI reads — see
labrador-demo-orchestrator module-lock.json and store.py). A cache miss builds
its interpretability via the deterministic ``build_interpretability``.

THE MANAGED-AGENT PATH IS RETAINED BUT NON-DEFAULT. The original cloud route
(vendored bun runner against the deployed Claude Managed Agent) is preserved as
``_run_via_agent`` and reached ONLY when the environment flag
``SIMULATION_USE_AGENT=1`` is set. The default run does not import, launch or
require any of it — no bun, no ``simulation/runtime`` install, no
``ANTHROPIC_API_KEY``, no network. The vendored runtime is left in place; the
default path simply does not depend on it.
"""

from __future__ import annotations

import copy
import json
import logging
import os
import shutil
import subprocess
from functools import lru_cache
from pathlib import Path

from .interpretability import build_interpretability

log = logging.getLogger("simulation.pipeline")

# Station root is the parent of this package; schemas and the bundled cache are
# resolved relative to this module so the whole thing runs from a copy of ONLY
# the module directory, with no reference to the parent repo.
_STATION_ROOT = Path(__file__).resolve().parent.parent
_SCHEMA_DIR = _STATION_ROOT / "schemas"
_CACHE_DIR = Path(__file__).resolve().parent / "cache"

# Env flag that opts INTO the non-default managed-agent path. Absent/anything
# other than "1" keeps the dependency-light local resolver.
_AGENT_FLAG = "SIMULATION_USE_AGENT"

# The vendored agent runner lives beside this file, under runtime/. Only touched
# on the guarded agent path.
_RUNTIME_DIR = Path(__file__).resolve().parent / "runtime"
_RUN_ONCE = _RUNTIME_DIR / "run-once.ts"
_BUN_FALLBACK = Path("/opt/homebrew/bin/bun")


class PipelineError(RuntimeError):
    """Base class for pipeline failures."""


class PipelineUnavailableError(PipelineError):
    """The pipeline cannot run: invalid request, or (agent path) missing tooling."""


class PipelineInvocationError(PipelineError):
    """The agent was invoked but did not return a parseable dossier."""


# --------------------------------------------------------------------------- #
# Request validation (defensive backstop — __main__ validates first).
# --------------------------------------------------------------------------- #
@lru_cache(maxsize=1)
def _input_validator():
    # jsonschema is OPTIONAL — return None when absent so the module runs on a
    # bare python3. _validate_request then falls back to a minimal manual check.
    try:
        from jsonschema import Draft202012Validator
    except ImportError:
        return None
    schema = json.loads((_SCHEMA_DIR / "input.schema.json").read_text())
    return Draft202012Validator(schema)


def _validate_request(request: dict) -> None:
    if not isinstance(request, dict):
        raise PipelineUnavailableError("request must be a JSON object")
    validator = _input_validator()
    if validator is None:
        acc = request.get("uniprot_accession")
        if not (isinstance(acc, str) and acc.strip()):
            raise PipelineUnavailableError(
                "request missing a non-empty 'uniprot_accession' (jsonschema not installed)"
            )
        return
    errors = sorted(
        validator.iter_errors(request), key=lambda e: list(e.absolute_path)
    )
    if errors:
        joined = "; ".join(
            f"{'/'.join(str(p) for p in e.absolute_path) or '<root>'}: {e.message}"
            for e in errors
        )
        raise PipelineUnavailableError(
            f"request does not satisfy schemas/input.schema.json: {joined}"
        )


# --------------------------------------------------------------------------- #
# The bundled cache and its key.
# --------------------------------------------------------------------------- #
def _norm_mech(value) -> str:
    """Absent/null/empty mechanism_hypothesis == 'unknown' (input.schema.json)."""
    if value is None or value == "":
        return "unknown"
    return value


def _dossier_key(dossier: dict) -> tuple | None:
    """The documented cache key of a stored dossier, from its `input` block."""
    inp = dossier.get("input") or {}
    acc = inp.get("uniprot_accession")
    if not acc:
        # Fall back to the resolved accession if the echoed one is missing.
        acc = (dossier.get("target") or {}).get("uniprot_accession")
    if not acc:
        return None
    return (acc, _norm_mech(inp.get("mechanism_hypothesis")), inp.get("as_of_date"))


@lru_cache(maxsize=1)
def _load_cache() -> tuple:
    """Load the bundled cache once. Returns a tuple of (filename, key, dossier),
    sorted by filename so lookup is deterministic across duplicate keys."""
    entries = []
    if not _CACHE_DIR.is_dir():
        log.warning("cache directory not found at %s — every target will miss", _CACHE_DIR)
        return tuple(entries)
    for path in sorted(_CACHE_DIR.glob("*.json")):
        try:
            dossier = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            log.warning("skipping unreadable cache file %s: %s", path.name, exc)
            continue
        key = _dossier_key(dossier)
        if key is None:
            log.warning("skipping cache file %s: no resolvable accession", path.name)
            continue
        entries.append((path.name, key, dossier))
    log.info("loaded %d cached dossiers from %s", len(entries), _CACHE_DIR)
    return tuple(entries)


def _request_key(request: dict) -> tuple:
    return (
        request.get("uniprot_accession"),
        _norm_mech(request.get("mechanism_hypothesis")),
        request.get("as_of_date"),
    )


def _lookup(request: dict):
    """Return (dossier, match_kind, filename, matched_key) or None.

    match_kind is 'exact' when the full tuple agrees, 'accession_only' when only
    the accession agrees. Deterministic: cache is iterated in sorted filename
    order and the first match at each precedence wins.
    """
    acc, mech, date = _request_key(request)
    entries = _load_cache()
    for fn, key, dossier in entries:
        if key == (acc, mech, date):
            return dossier, "exact", fn, key
    for fn, key, dossier in entries:
        if key[0] == acc:
            return dossier, "accession_only", fn, key
    return None


# --------------------------------------------------------------------------- #
# Stamping: mark local runtime maturity + cached-dossier qualifier.
# --------------------------------------------------------------------------- #
def _stamp_extensions(interp: dict, *, runtime_maturity: str, qualifiers: list[str],
                      output_origin: str, extra: dict | None = None) -> None:
    ext = interp.setdefault("extensions", {})
    ext["runtime_maturity"] = runtime_maturity
    quals = ext.setdefault("qualifiers", [])
    if not isinstance(quals, list):
        quals = []
        ext["qualifiers"] = quals
    for q in qualifiers:
        if q not in quals:
            quals.append(q)
    ext["output_origin"] = output_origin
    if extra:
        ext.update(extra)


def _mark_cache_hit(dossier: dict, request: dict, match_kind: str, filename: str,
                    matched_key: tuple) -> dict:
    """Return a stamped deep copy of a cached dossier. Keeps its own authored
    interpretability; only ADDS the LOCAL/CACHED_DOSSIER extension markers, and
    on an accession-only fallback appends a CACHE_KEY_FALLBACK limitation. Never
    alters any scientific value."""
    out = copy.deepcopy(dossier)
    interp = out.get("interpretability")
    if not isinstance(interp, dict):
        # Cached dossiers ship with interpretability; rebuild only if one is
        # somehow absent, from the dossier's own real fields (no fabrication).
        interp = build_interpretability(out)
        out["interpretability"] = interp

    quals = ["CACHED_DOSSIER"]
    if match_kind == "accession_only":
        quals.append("CACHE_KEY_FALLBACK_ACCESSION_ONLY")
    _stamp_extensions(
        interp,
        runtime_maturity="LOCAL",
        qualifiers=quals,
        output_origin="cached_dossier",
        extra={
            "cache_hit": {
                "kind": match_kind,
                "cache_file": filename,
                "matched_key": {
                    "uniprot_accession": matched_key[0],
                    "mechanism_hypothesis": matched_key[1],
                    "as_of_date": matched_key[2],
                },
            }
        },
    )

    if match_kind == "accession_only":
        req_acc, req_mech, req_date = _request_key(request)
        limitations = interp.setdefault("limitations", [])
        if isinstance(limitations, list):
            limitations.append({
                "code": "CACHE_KEY_FALLBACK",
                "severity": "WARNING",
                "message": (
                    "No exact cached dossier for the requested "
                    f"(uniprot_accession, mechanism_hypothesis, as_of_date) = "
                    f"({req_acc}, {req_mech}, {req_date}); returned the "
                    f"accession-only match '{filename}' keyed on "
                    f"(mechanism_hypothesis={matched_key[1]}, as_of_date={matched_key[2]}). "
                    "Its mechanism-specific answer may not correspond to the "
                    "requested mechanism, which changes the dossier."
                ),
                "field_path": "input.mechanism_hypothesis",
            })
    log.info("cache %s hit on %s for %s", match_kind, filename, request.get("uniprot_accession"))
    return out


# --------------------------------------------------------------------------- #
# Cache miss: an honest, schema-valid, empty dossier. No fabricated numbers.
# --------------------------------------------------------------------------- #
def _honest_empty_dossier(request: dict) -> dict:
    acc = request.get("uniprot_accession")
    reason = (
        "This target was not precomputed in the dependency-light cached/local "
        "build. The local runner replays a fixed bundled cache of dossiers and "
        "does not run the comp-chem stack, so no evidence was retrieved and no "
        "value is reported for this (accession, mechanism_hypothesis, as_of_date)."
    )

    def nf(field: str) -> dict:
        return {"field": field, "reason": reason, "code": "NOT_COMPUTED_LOCALLY"}

    dossier = {
        "input": {
            "uniprot_accession": request.get("uniprot_accession"),
            "as_of_date": request.get("as_of_date"),
            "disease_context": request.get("disease_context"),
            "interaction_to_disrupt": request.get("interaction_to_disrupt"),
            "mechanism_hypothesis": request.get("mechanism_hypothesis"),
        },
        # target.uniprot_accession echoes the caller-supplied, schema-valid
        # accession. It is what was asked, not a resolved/looked-up value.
        "target": {
            "uniprot_accession": acc,
            "gene_symbol": None,
            "protein_name": None,
            "organism": None,
            "sequence_length": None,
            "sources": [],
        },
        "as_of_date": request.get("as_of_date"),
        "verdict": "insufficient_evidence",
        "verdict_basis": "none",
        "axis_conflict": None,
        "target_precedent": {
            "chembl_target_id": None,
            "distinct_actives": None,
            "modality_unknown_count": None,
            "best_potency_nm": None,
            "best_potency_modality": None,
            "best_potency_assay": None,
            "best_potency_characterised": None,
            "approved_small_molecules_count": None,
            "approved_small_molecules": [],
            "clinical_stage_small_molecules": [],
            "patents": {"count": None, "source": None},
            "terminated_programs": [],
            "as_of_leakage": [],
            "sources": [],
        },
        "biologic_precedent": {
            "approved_biologics": [],
            "note": (
                "Not computed in this local cached build; no biologic precedent "
                "was retrieved. Presence of an approved biologic is target "
                "validation, NOT small-molecule tractability."
            ),
        },
        "family_precedent": {
            "pfam": None,
            "family_actives": None,
            "best_family_potency_nm": None,
            "best_family_potency_modality": None,
            "best_family_target": None,
            "sources": [],
        },
        "structural_neighbour_precedent": {
            "query_structure": None,
            "neighbours": [],
            "sources": [],
        },
        "pocket_neighbour_precedent": {
            "candidates": [],
            "sources": [],
        },
        "structure": {
            "tier": "none",
            "pdb_id": None,
            "total_pdb_structures": None,
            "holo_count": None,
            "apo_count": None,
            "ensemble_used": [],
            "sources": [],
        },
        "tractability": {
            "cryptic_pocket_risk": "undetermined",
            "cryptic_mechanism": "undetermined",
            "mdpocket_site_definition_used": "none",
            "site_hypothesis_basis": None,
            "caveat": reason,
            "pocket_residues": [],
            "sources": [],
        },
        "affinity": {
            "positive_control_ligand": None,
            "positive_control_measured_nm": None,
            "positive_control_predicted_nm": None,
            "reliable": None,
            "predictions": [],
            "sources": [],
        },
        "falsification": {
            "checks_run": [
                "No precedent, structure or tractability evidence was retrieved "
                "in this local cached build, so there was no claim to attempt to "
                "falsify."
            ],
            "findings": [],
            "survived": None,
        },
        "next_experiment": {
            "description": (
                "Run this target through the full druggability pipeline "
                "(structure-select, pocket-scan, precedent-lookup, falsification, "
                "assemble-dossier) to produce a real dossier; the local cached "
                "build cannot compute one."
            ),
            "rationale": (
                "The dependency-light local runner only replays a fixed bundled "
                "cache of precomputed dossiers and does not run the comp-chem stack."
            ),
            "resolves": (
                "Whether the target has retrievable small-molecule precedent and a "
                "computationally tractable pocket."
            ),
        },
        "not_found": [
            nf("target_precedent"),
            nf("structure"),
            nf("tractability"),
            nf("family_precedent"),
            nf("structural_neighbour_precedent"),
            nf("pocket_neighbour_precedent"),
            nf("affinity"),
        ],
    }

    # Build interpretability deterministically from these (all-null) fields, then
    # stamp the local runtime maturity. The NOT_COMPUTED_LOCALLY limitations flow
    # in from not_found via build_interpretability.
    interp = build_interpretability(dossier)
    _stamp_extensions(
        interp,
        runtime_maturity="LOCAL",
        qualifiers=["NOT_COMPUTED_LOCALLY"],
        output_origin="not_computed_locally",
    )
    dossier["interpretability"] = interp
    log.info("cache miss for %s — returning honest insufficient_evidence dossier", acc)
    return dossier


def _resolve_local(request: dict) -> dict:
    hit = _lookup(request)
    if hit is not None:
        dossier, match_kind, filename, matched_key = hit
        return _mark_cache_hit(dossier, request, match_kind, filename, matched_key)
    return _honest_empty_dossier(request)


# --------------------------------------------------------------------------- #
# Public entrypoint.
# --------------------------------------------------------------------------- #
def run_pipeline(request: dict) -> dict:
    """Return a schema-valid druggability dossier for one request.

    DEFAULT (dependency-light, offline): resolve against the bundled cache; on a
    miss return an honest insufficient-evidence dossier. Never fabricates a
    scientific value.

    OPT-IN (``SIMULATION_USE_AGENT=1``): drive the vendored managed-agent runner
    instead. That path needs bun, the vendored runtime and ANTHROPIC_API_KEY and
    is NOT used by the default run.
    """
    _validate_request(request)

    if os.environ.get(_AGENT_FLAG) == "1":
        log.info("%s=1 set — using the non-default managed-agent path", _AGENT_FLAG)
        return _run_via_agent(request)

    return _resolve_local(request)


# --------------------------------------------------------------------------- #
# NON-DEFAULT managed-agent path (guarded by SIMULATION_USE_AGENT=1).
# The default run never reaches this code.
# --------------------------------------------------------------------------- #
def _resolve_bun() -> str | None:
    bun = shutil.which("bun")
    if bun:
        return bun
    if _BUN_FALLBACK.is_file():
        return str(_BUN_FALLBACK)
    return None


def _build_task_prose(request: dict) -> str:
    acc = request.get("uniprot_accession")
    lines = [
        "Assemble a druggability dossier for the following target.",
        f"uniprot_accession: {acc}",
    ]
    for field in ("as_of_date", "disease_context", "interaction_to_disrupt", "mechanism_hypothesis"):
        value = request.get(field)
        if value is not None:
            lines.append(f"{field}: {value}")
    return "\n".join(lines)


def _extract_dossier_json(stdout: str) -> dict:
    text = stdout.strip()
    if not text:
        raise PipelineInvocationError("agent produced no stdout to parse a dossier from")
    try:
        parsed = json.loads(text)
        if isinstance(parsed, dict):
            return parsed
    except json.JSONDecodeError:
        pass
    obj = _last_json_object(text)
    if obj is None:
        raise PipelineInvocationError(
            "could not find a JSON object in the agent reply; the agent must paste the "
            "complete dossier JSON into its final reply (see CLAUDE.md)"
        )
    return obj


def _last_json_object(text: str) -> dict | None:
    result: dict | None = None
    depth = 0
    start = -1
    in_str = False
    escape = False
    for i, ch in enumerate(text):
        if in_str:
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == '"':
                in_str = False
            continue
        if ch == '"':
            in_str = True
        elif ch == "{":
            if depth == 0:
                start = i
            depth += 1
        elif ch == "}":
            if depth > 0:
                depth -= 1
                if depth == 0 and start >= 0:
                    candidate = text[start : i + 1]
                    try:
                        parsed = json.loads(candidate)
                    except json.JSONDecodeError:
                        parsed = None
                    if isinstance(parsed, dict):
                        result = parsed
    return result


def _run_via_agent(request: dict) -> dict:
    """Invoke the deployed Claude Managed Agent via the vendored bun runner.

    Reached only when SIMULATION_USE_AGENT=1. Requires bun on PATH, the vendored
    runner present, ANTHROPIC_API_KEY, network access, and a deployed agent.
    Raises PipelineUnavailableError / PipelineInvocationError; never fabricates.
    """
    if not _RUN_ONCE.is_file():
        raise PipelineUnavailableError(
            f"vendored runner not found at {_RUN_ONCE} — the agent path is unavailable."
        )
    bun = _resolve_bun()
    if bun is None:
        raise PipelineUnavailableError(
            "`bun` is not on PATH — the vendored runner cannot be launched. Install "
            "bun, then `cd simulation/runtime && bun install`."
        )

    task = _build_task_prose(request)
    cmd = [bun, str(_RUN_ONCE), "--once", task]
    log.info("invoking vendored runner (task for %s)", request.get("uniprot_accession"))

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(_RUNTIME_DIR),
            capture_output=True,
            text=True,
            check=False,
            env=os.environ.copy(),
        )
    except OSError as exc:
        raise PipelineInvocationError(f"failed to launch the vendored runner: {exc}")

    if proc.stderr:
        for line in proc.stderr.splitlines():
            log.info("[runner] %s", line)

    if proc.returncode != 0:
        tail = proc.stderr[-2000:].strip()
        if proc.returncode == 3 or "deployment.agent_id" in proc.stderr:
            raise PipelineUnavailableError(
                "the managed agent has not been deployed (no deployment.agent_id in "
                "manifest.json). Deployment is the integrator's job — deploy it "
                "(run: bun run deploy), then re-run. Runner stderr:\n" + tail
            )
        raise PipelineInvocationError(
            f"vendored runner exited {proc.returncode}. Runner stderr:\n{tail}"
        )

    return _extract_dossier_json(proc.stdout)
