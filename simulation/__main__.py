"""`python -m simulation` entrypoint.

The ONE documented command:

    python -m simulation run --input input.json --output output.json

Exit code 0 for success, nonzero for failure. All machine-readable results are
written to the --output file. Human-facing logging goes to STDERR; STDOUT stays
empty.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

try:  # jsonschema is OPTIONAL at runtime — the cached dossiers are pre-validated,
    # so the module stays runnable on a bare python3 with no third-party deps.
    from jsonschema import Draft202012Validator
except ImportError:  # pragma: no cover
    Draft202012Validator = None  # type: ignore[assignment]

from .interpretability import build_interpretability
from .pipeline import PipelineError, run_pipeline

# Exit codes (0 = success; everything else is a distinct failure mode).
EXIT_OK = 0
EXIT_USAGE = 1
EXIT_INVALID_INPUT = 2
EXIT_PIPELINE_FAILED = 3
EXIT_VALIDATION_FAILED = 4

# Station root is the parent of this package; schemas live beside it.
STATION_ROOT = Path(__file__).resolve().parent.parent
SCHEMA_DIR = STATION_ROOT / "schemas"

log = logging.getLogger("simulation")


def _load_schema(name: str) -> dict:
    return json.loads((SCHEMA_DIR / name).read_text())


def _format_errors(validator, instance: object) -> list[str]:
    """Return sorted, human-readable violation strings for a schema validation."""
    out = []
    for err in sorted(validator.iter_errors(instance), key=lambda e: list(e.absolute_path)):
        loc = "/".join(str(p) for p in err.absolute_path) or "<root>"
        out.append(f"  at {loc}: {err.message}")
    return out


def _run(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    output_path = Path(args.output)

    # --- (a) Read and validate the request. On invalid input: write nothing, exit nonzero.
    try:
        request = json.loads(input_path.read_text())
    except FileNotFoundError:
        log.error("input file not found: %s", input_path)
        return EXIT_INVALID_INPUT
    except json.JSONDecodeError as exc:
        log.error("input file is not valid JSON (%s): %s", input_path, exc)
        return EXIT_INVALID_INPUT

    if Draft202012Validator is not None:
        input_errors = _format_errors(Draft202012Validator(_load_schema("input.schema.json")), request)
        if input_errors:
            log.error("request does not satisfy schemas/input.schema.json:")
            for line in input_errors:
                log.error("%s", line)
            log.error("nothing written to %s", output_path)
            return EXIT_INVALID_INPUT
    else:
        log.warning("jsonschema not installed — skipping input schema validation; minimal check only")
        acc = request.get("uniprot_accession") if isinstance(request, dict) else None
        if not (isinstance(acc, str) and acc.strip()):
            log.error("request missing a non-empty 'uniprot_accession'; nothing written to %s", output_path)
            return EXIT_INVALID_INPUT

    log.info("request valid; target %s", request.get("uniprot_accession"))

    # --- (b) Produce the dossier via the managed agent. Fail loudly, never fabricate.
    try:
        dossier = run_pipeline(request)
    except PipelineError as exc:
        log.error("dossier production failed: %s", exc)
        log.error("nothing written to %s", output_path)
        return EXIT_PIPELINE_FAILED
    except Exception as exc:  # noqa: BLE001 - surface any unexpected pipeline error to stderr
        log.error("unexpected error producing dossier: %s: %s", type(exc).__name__, exc)
        log.error("nothing written to %s", output_path)
        return EXIT_PIPELINE_FAILED

    if not isinstance(dossier, dict):
        log.error("run_pipeline did not return a dossier object; got %s", type(dossier).__name__)
        return EXIT_PIPELINE_FAILED

    # --- (c) Ensure the interpretability object is present.
    # The local cached path returns dossiers that already carry a real,
    # authored (or locally-built + stamped) interpretability block — KEEP those.
    # Only build one here when a dossier arrives without it (e.g. the non-default
    # managed-agent path returns a raw dossier).
    if not isinstance(dossier.get("interpretability"), dict):
        try:
            dossier["interpretability"] = build_interpretability(dossier)
        except Exception as exc:  # noqa: BLE001
            log.error("failed to build interpretability object: %s: %s", type(exc).__name__, exc)
            # Still write the dossier below so it can be inspected; then fail.
            _write_output(dossier, output_path)
            return EXIT_VALIDATION_FAILED

    # --- (d) Validate the dossier and its interpretability object; write regardless.
    if Draft202012Validator is not None:
        output_errors = _format_errors(Draft202012Validator(_load_schema("output.schema.json")), dossier)
        interp_errors = _format_errors(
            Draft202012Validator(_load_schema("interpretability.schema.json")), dossier.get("interpretability")
        )
    else:
        log.warning(
            "jsonschema not installed — skipping output/interpretability validation "
            "(cached dossiers are pre-validated at build time)"
        )
        output_errors = []
        interp_errors = []

    _write_output(dossier, output_path)

    ok = True
    if output_errors:
        ok = False
        log.error("dossier does not satisfy schemas/output.schema.json:")
        for line in output_errors:
            log.error("%s", line)
    if interp_errors:
        ok = False
        log.error("interpretability object does not satisfy schemas/interpretability.schema.json:")
        for line in interp_errors:
            log.error("%s", line)

    if not ok:
        log.error("dossier written to %s for inspection, but validation FAILED", output_path)
        return EXIT_VALIDATION_FAILED

    # --- (e) Success.
    log.info("dossier written to %s; output and interpretability both valid", output_path)
    return EXIT_OK


def _write_output(dossier: dict, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(dossier, indent=2, ensure_ascii=False))


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        stream=sys.stderr,  # logs to STDERR only; STDOUT stays clean.
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    parser = argparse.ArgumentParser(
        prog="python -m simulation",
        description="Druggability-dossier simulation module (standardized runnable interface).",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    run_parser = sub.add_parser("run", help="produce a dossier for one target")
    run_parser.add_argument("--input", required=True, help="path to the request JSON (input.schema.json)")
    run_parser.add_argument("--output", required=True, help="path to write the dossier JSON")
    run_parser.set_defaults(func=_run)

    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help(sys.stderr)
        return EXIT_USAGE
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
