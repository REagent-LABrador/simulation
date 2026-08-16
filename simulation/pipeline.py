"""Produce a druggability dossier by invoking the Claude Managed Agent.

This station's real computation is a Claude Managed Agent (see the repo README:
"prototype in Claude Code, ship on the Claude Developer Platform"). There is NO
pure-Python re-implementation of the science -- the dossier is produced by the
deployed agent running its skills (structure-select, pocket-scan, precedent-lookup,
falsification, assemble-dossier) in its cloud sandbox. So ``run_pipeline`` drives
the documented headless invocation route and parses the dossier JSON the agent
emits. It does NOT and MUST NOT fabricate a dossier: when the agent cannot be
invoked, it raises loudly.

The documented route (README step 3, scripts/console.ts) is:

    bun run console <name> -- --once "<task prose>"

which starts a session against the deployed Managed Agent, answers any custom-tool
round-trips in this process, prints the agent's final reply (the dossier JSON) to
STDOUT, and logs progress to STDERR. The agent's CLAUDE.md instructs it to paste
the complete dossier JSON into its final reply, which is what we parse here.

A live run requires ANTHROPIC_API_KEY, network access, ``bun`` on PATH, and the
agent to have been deployed (``manifest.deployment.agent_id`` set). Those are
expected for a real run; when any is absent this raises a typed error rather than
hanging or inventing data.
"""

from __future__ import annotations

import json
import logging
import os
import shutil
import subprocess
from pathlib import Path

log = logging.getLogger("simulation.pipeline")

# The managed agent lives at managed/<AGENT_DIR_NAME>/ in the repo. `bun run
# console <name>` resolves the agent by this directory name (lib/claude-managed-agent.ts
# loadManagedAgent reads managed/<name>/).
AGENT_DIR_NAME = "druggability-dossier"


class PipelineError(RuntimeError):
    """Base class for pipeline failures."""


class PipelineUnavailableError(PipelineError):
    """The managed agent cannot be invoked: missing credentials, tools, or deployment."""


class PipelineInvocationError(PipelineError):
    """The agent was invoked but did not return a parseable dossier."""


def _repo_root() -> Path:
    """Walk up from this file to the repo root (the dir holding package.json + managed/)."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "package.json").is_file() and (parent / "managed").is_dir():
            return parent
    raise PipelineUnavailableError(
        "could not locate the repo root (a directory containing package.json and managed/) "
        f"above {here}"
    )


def _build_task_prose(request: dict) -> str:
    """Render the five contract fields as the prose the agent parses from {task}.

    The agent receives the request as prose inside a single {task} argument and
    parses the five fields back out (see CLAUDE.md '## Contract'). We echo exactly
    what the caller supplied; we do NOT infer or back-fill absent fields.
    """
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
    """Pull the dossier JSON object out of the agent's final reply text."""
    text = stdout.strip()
    if not text:
        raise PipelineInvocationError("agent produced no stdout to parse a dossier from")
    # Fast path: the whole reply is the JSON object.
    try:
        parsed = json.loads(text)
        if isinstance(parsed, dict):
            return parsed
    except json.JSONDecodeError:
        pass
    # Otherwise, find the last balanced top-level {...} block in the reply.
    obj = _last_json_object(text)
    if obj is None:
        raise PipelineInvocationError(
            "could not find a JSON object in the agent reply; the agent must paste the "
            "complete dossier JSON into its final reply (see CLAUDE.md)"
        )
    return obj


def _last_json_object(text: str) -> dict | None:
    """Scan for balanced {...} objects and return the last one that parses to a dict."""
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


def run_pipeline(request: dict) -> dict:
    """Invoke the deployed Claude Managed Agent and return the dossier it produces.

    Raises ``PipelineUnavailableError`` when the agent cannot be invoked (no
    ANTHROPIC_API_KEY, no ``bun``, or the agent has not been deployed) and
    ``PipelineInvocationError`` when the agent runs but does not return a
    parseable dossier. Never fabricates a dossier.
    """
    if not isinstance(request, dict) or not request.get("uniprot_accession"):
        raise PipelineUnavailableError("request must be a dict carrying a uniprot_accession")

    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise PipelineUnavailableError(
            "ANTHROPIC_API_KEY is not set — the managed agent cannot be invoked. This is a "
            "live-run requirement (see .env.example); the module fails loudly rather than "
            "fabricating a dossier."
        )

    bun = shutil.which("bun")
    if bun is None:
        raise PipelineUnavailableError(
            "`bun` is not on PATH — the documented invocation route "
            "(`bun run console <name> -- --once ...`) is unavailable."
        )

    repo_root = _repo_root()
    manifest_path = repo_root / "managed" / AGENT_DIR_NAME / "manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise PipelineUnavailableError(f"cannot read agent manifest at {manifest_path}: {exc}")
    agent_id = (manifest.get("deployment") or {}).get("agent_id")
    if not agent_id:
        raise PipelineUnavailableError(
            f"managed agent '{AGENT_DIR_NAME}' has no deployment.agent_id in {manifest_path} — "
            f"it has not been deployed yet. Run: bun run deploy {AGENT_DIR_NAME}"
        )

    task = _build_task_prose(request)
    cmd = [bun, "run", "console", AGENT_DIR_NAME, "--", "--once", task]
    log.info("invoking managed agent: %s (task for %s)", " ".join(cmd[:5]), request.get("uniprot_accession"))

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(repo_root),
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise PipelineInvocationError(f"failed to launch the agent invocation: {exc}")

    if proc.stderr:
        # The agent's progress/logging goes to its stderr; surface it through ours.
        for line in proc.stderr.splitlines():
            log.info("[agent] %s", line)

    if proc.returncode != 0:
        raise PipelineInvocationError(
            f"agent invocation exited {proc.returncode}. Last stderr:\n{proc.stderr[-2000:]}"
        )

    return _extract_dossier_json(proc.stdout)
