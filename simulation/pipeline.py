"""Produce a druggability dossier in one explicit execution mode.

``replay`` is dependency-light and resolves the request against a BUNDLED CACHE
of REAL dossiers shipped inside this module at ``simulation/cache/``, keyed on
the documented tuple

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

``live`` calls the already-configured managed-agent checkout named by
``LABRADOR_RUNTIME_ROOT``.  The split station never deploys or mutates that
checkout.  It only runs its supported ``scripts/console.ts`` entrypoint and
fails with a stable reason code when the runtime, deployment, credential, or
provider is unavailable.  There is deliberately no fallback from ``live`` to
``replay``.
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

# The live path uses the full LABrador checkout because the split repository
# contains the station contract but not the managed-agent client/runtime.
_RUNTIME_ROOT_ENV = "LABRADOR_RUNTIME_ROOT"
_LIVE_TIMEOUT_SECONDS_ENV = "SIMULATION_LIVE_TIMEOUT_SECONDS"
_DEFAULT_LIVE_TIMEOUT_SECONDS = 90 * 60
_AGENT_NAMES = (
    "small-molecule-tractability-review",
    "druggability-dossier",
)
_BUN_FALLBACK = Path("/opt/homebrew/bin/bun")


class PipelineError(RuntimeError):
    """Base class for pipeline failures."""

    code = "PIPELINE_FAILED"

    def __init__(self, message: str, *, code: str | None = None) -> None:
        super().__init__(message)
        if code is not None:
            self.code = code


class PipelineUnavailableError(PipelineError):
    """The pipeline cannot run: invalid request, or (agent path) missing tooling."""

    code = "LIVE_UNAVAILABLE"


class PipelineInvocationError(PipelineError):
    """The agent was invoked but did not return a parseable dossier."""

    code = "PROVIDER_FAILED"


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
def run_pipeline(request: dict, *, mode: str) -> dict:
    """Return a schema-valid dossier from exactly one requested execution mode.

    ``replay`` resolves the bundled cache and never makes a provider call.
    ``live`` invokes the configured managed-agent runtime and never falls back.
    """
    _validate_request(request)

    if mode == "live":
        log.info("live mode selected — invoking the configured managed-agent runtime")
        return _run_via_agent(request)

    if mode == "replay":
        log.info("replay mode selected — resolving the bundled dossier cache")
        return _resolve_local(request)

    raise PipelineUnavailableError(
        f"unsupported execution mode {mode!r}; expected 'live' or 'replay'",
        code="INVALID_MODE",
    )


# --------------------------------------------------------------------------- #
# Live managed-agent bridge. It never deploys or modifies the runtime checkout.
# --------------------------------------------------------------------------- #
def _resolve_bun() -> str | None:
    bun = shutil.which("bun")
    if bun:
        return bun
    if _BUN_FALLBACK.is_file():
        return str(_BUN_FALLBACK)
    return None


def _runtime_root() -> Path:
    configured = os.environ.get(_RUNTIME_ROOT_ENV, "").strip()
    if not configured:
        raise PipelineUnavailableError(
            f"{_RUNTIME_ROOT_ENV} is not set; point it at a LABrador checkout "
            "that already contains the deployed tractability managed agent",
            code="RUNTIME_NOT_CONFIGURED",
        )
    root = Path(configured).expanduser().resolve()
    if not root.is_dir():
        raise PipelineUnavailableError(
            f"{_RUNTIME_ROOT_ENV} does not name a directory: {root}",
            code="RUNTIME_NOT_FOUND",
        )
    if not (root / "scripts" / "console.ts").is_file():
        raise PipelineUnavailableError(
            f"LABrador headless runner not found at {root / 'scripts' / 'console.ts'}",
            code="RUNTIME_INCOMPLETE",
        )
    return root


def _managed_agent(root: Path) -> tuple[str, dict]:
    for name in _AGENT_NAMES:
        path = root / "managed" / name / "manifest.json"
        if not path.is_file():
            continue
        try:
            manifest = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError) as exc:
            raise PipelineUnavailableError(
                f"managed-agent manifest is unreadable at {path}: {exc}",
                code="MANIFEST_INVALID",
            ) from exc
        if not isinstance(manifest, dict):
            raise PipelineUnavailableError(
                f"managed-agent manifest is not a JSON object: {path}",
                code="MANIFEST_INVALID",
            )
        deployment = manifest.get("deployment")
        if not isinstance(deployment, dict) or not deployment.get("agent_id"):
            raise PipelineUnavailableError(
                f"managed agent {name!r} has no existing deployment.agent_id; "
                "this runner will not deploy it",
                code="DEPLOYMENT_NOT_CONFIGURED",
            )
        return name, manifest
    expected = ", ".join(str(root / "managed" / name / "manifest.json") for name in _AGENT_NAMES)
    raise PipelineUnavailableError(
        f"no tractability managed-agent manifest found; checked {expected}",
        code="MANAGED_AGENT_NOT_INSTALLED",
    )


def _live_timeout_seconds() -> int:
    raw = os.environ.get(_LIVE_TIMEOUT_SECONDS_ENV, "").strip()
    if not raw:
        return _DEFAULT_LIVE_TIMEOUT_SECONDS
    try:
        value = int(raw)
    except ValueError as exc:
        raise PipelineUnavailableError(
            f"{_LIVE_TIMEOUT_SECONDS_ENV} must be a positive integer",
            code="RUNTIME_CONFIGURATION_INVALID",
        ) from exc
    if value <= 0:
        raise PipelineUnavailableError(
            f"{_LIVE_TIMEOUT_SECONDS_ENV} must be a positive integer",
            code="RUNTIME_CONFIGURATION_INVALID",
        )
    return value


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
        raise PipelineInvocationError(
            "agent produced no stdout to parse a dossier from",
            code="INVALID_OUTPUT",
        )
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
            "complete dossier JSON into its final reply (see CLAUDE.md)",
            code="INVALID_OUTPUT",
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
    """Invoke the existing managed agent through its supported headless CLI."""

    root = _runtime_root()
    agent_name, _manifest = _managed_agent(root)
    bun = _resolve_bun()
    if bun is None:
        raise PipelineUnavailableError(
            "`bun` is not on PATH, so the LABrador managed-agent runner cannot start",
            code="BINARY_MISSING",
        )
    if not os.environ.get("ANTHROPIC_API_KEY", "").strip():
        raise PipelineUnavailableError(
            "ANTHROPIC_API_KEY is not set for live managed-agent execution",
            code="CREDENTIAL_MISSING",
        )

    task = _build_task_prose(request)
    timeout_seconds = _live_timeout_seconds()
    cmd = [
        bun,
        "scripts/console.ts",
        agent_name,
        "--",
        "--once",
        task,
        "--quiet",
        "--timeout",
        str(timeout_seconds),
    ]
    log.info(
        "invoking managed agent %s for %s",
        agent_name,
        request.get("uniprot_accession"),
    )

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(root),
            capture_output=True,
            text=True,
            check=False,
            env=os.environ.copy(),
            timeout=timeout_seconds + 30,
        )
    except subprocess.TimeoutExpired as exc:
        raise PipelineInvocationError(
            f"managed-agent execution exceeded {timeout_seconds} seconds",
            code="PROVIDER_TIMEOUT",
        ) from exc
    except OSError as exc:
        raise PipelineInvocationError(
            f"failed to launch the managed-agent runner: {exc}",
            code="RUNTIME_LAUNCH_FAILED",
        ) from exc

    if proc.stderr:
        for line in proc.stderr.splitlines():
            log.info("[runner] %s", line)

    if proc.returncode != 0:
        tail = proc.stderr[-2000:].strip()
        stderr_lower = proc.stderr.lower()
        if "deployment.agent_id" in proc.stderr or "not deployed" in stderr_lower:
            raise PipelineUnavailableError(
                "the configured managed agent has no usable deployment; this runner "
                "will not deploy it. Runner stderr:\n" + tail,
                code="DEPLOYMENT_NOT_CONFIGURED",
            )
        if "api key" in stderr_lower or "credential" in stderr_lower or "unauthorized" in stderr_lower:
            raise PipelineUnavailableError(
                "the live provider rejected or could not find a required credential. "
                "Runner stderr:\n" + tail,
                code="CREDENTIAL_MISSING",
            )
        if "timed out" in stderr_lower or "timeout" in stderr_lower:
            raise PipelineInvocationError(
                "the live provider timed out. Runner stderr:\n" + tail,
                code="PROVIDER_TIMEOUT",
            )
        if any(
            marker in stderr_lower
            for marker in (
                "not on path",
                "command not found",
                "no such file or directory",
                "enoent",
                "missing binary",
                "missing dependency",
            )
        ):
            raise PipelineUnavailableError(
                "the managed-agent runtime is missing a required dependency. "
                "Runner stderr:\n" + tail,
                code="DEPENDENCY_MISSING",
            )
        raise PipelineInvocationError(
            f"managed-agent runner exited {proc.returncode}. Runner stderr:\n{tail}",
            code="PROVIDER_FAILED",
        )

    dossier = _extract_dossier_json(proc.stdout)
    interp = dossier.get("interpretability")
    if isinstance(interp, dict):
        _stamp_extensions(
            interp,
            runtime_maturity="MANAGED_AGENT",
            qualifiers=["LIVE"],
            output_origin="live_provider",
        )
    return dossier
