"""Offline contract tests for the simulation module. No paid calls, no network.

Run standalone:
    /Users/bb/micromamba/envs/druggability/bin/python simulation/test_module.py
or:
    micromamba run -n druggability python -m unittest simulation.test_module

The DEFAULT run is now the dependency-light LOCAL resolver over the bundled
cache in simulation/cache/ (Python + jsonschema only — no cloud, no Modal, no
Paperclip, no managed agent, no API key). So the end-to-end tests exercise the
REAL default path with NO monkeypatch: a cache hit on examples/input.json (the
real IRAK4 dossier) and a cache miss on an unknown accession (an honest,
schema-valid insufficient-evidence dossier). The managed-agent behaviours are
still covered by monkeypatching run_pipeline, since that path is non-default and
requires tooling absent offline.
"""

from __future__ import annotations

import io
import json
import logging
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from unittest import mock

from jsonschema import Draft202012Validator

# Allow `python simulation/test_module.py` from the station root by ensuring the
# station root (which contains the `simulation` package) is importable.
STATION_ROOT = Path(__file__).resolve().parent.parent
if str(STATION_ROOT) not in sys.path:
    sys.path.insert(0, str(STATION_ROOT))

import simulation.__main__ as cli  # noqa: E402
from simulation.interpretability import build_interpretability  # noqa: E402
from simulation.pipeline import PipelineUnavailableError  # noqa: E402

SCHEMA_DIR = STATION_ROOT / "schemas"
EXAMPLES = STATION_ROOT / "examples"
FIXTURES = STATION_ROOT / ".claude" / "skills" / "assemble-dossier" / "examples" / "integration-fixtures"


def _schema(name: str) -> dict:
    return json.loads((SCHEMA_DIR / name).read_text())


def _run_cli(argv: list[str]) -> tuple[int, str, str]:
    """Invoke cli.main(argv), capturing stdout and any logging emitted to stderr."""
    log_buf = io.StringIO()
    handler = logging.StreamHandler(log_buf)
    handler.setLevel(logging.DEBUG)
    root = logging.getLogger()
    root.addHandler(handler)
    prev_level = root.level
    root.setLevel(logging.DEBUG)
    out_buf = io.StringIO()
    try:
        with redirect_stdout(out_buf):
            code = cli.main(argv)
    finally:
        root.removeHandler(handler)
        root.setLevel(prev_level)
    return code, out_buf.getvalue(), log_buf.getvalue()


class EndToEndDefaultPathCacheHit(unittest.TestCase):
    """The REAL default path: no monkeypatch, offline, resolves from the cache."""

    def test_examples_input_is_a_real_cache_hit(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "o.json"
            code, stdout, stderr = _run_cli(
                ["run", "--input", str(EXAMPLES / "input.json"), "--output", str(out_path)]
            )
            self.assertEqual(code, cli.EXIT_OK, msg=stderr)
            self.assertEqual(stdout, "", "stdout must stay empty; results go to --output")
            self.assertTrue(out_path.exists(), "output file must be written")

            written = json.loads(out_path.read_text())
            # Schema-valid dossier + schema-valid interpretability.
            Draft202012Validator(_schema("output.schema.json")).validate(written)
            self.assertIn("interpretability", written)
            Draft202012Validator(_schema("interpretability.schema.json")).validate(
                written["interpretability"]
            )
            self.assertEqual(written["interpretability"]["schema_version"], "1.0.0")

            # It is the real IRAK4 dossier, carrying its real numbers.
            self.assertEqual(written["target"]["uniprot_accession"], "Q9NWZ3")
            self.assertEqual(written["verdict"], "small_molecule_tractable")

            # Marked as a LOCAL cached dossier for the orchestrator/UI.
            ext = written["interpretability"]["extensions"]
            self.assertEqual(ext["runtime_maturity"], "LOCAL")
            self.assertIn("CACHED_DOSSIER", ext["qualifiers"])
            self.assertEqual(ext["cache_hit"]["kind"], "exact")


class EndToEndUnknownTargetHonest(unittest.TestCase):
    """A target absent from the cache: valid dossier, honestly empty, no fabrication."""

    def test_unknown_accession_returns_insufficient_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            inp = Path(tmp) / "in.json"
            inp.write_text(json.dumps({
                "uniprot_accession": "P99999",
                "as_of_date": None,
                "disease_context": None,
                "interaction_to_disrupt": None,
                "mechanism_hypothesis": None,
            }))
            out_path = Path(tmp) / "o.json"
            code, stdout, stderr = _run_cli(
                ["run", "--input", str(inp), "--output", str(out_path)]
            )
            self.assertEqual(code, cli.EXIT_OK, msg=stderr)
            self.assertEqual(stdout, "", "stdout must stay empty")

            d = json.loads(out_path.read_text())
            Draft202012Validator(_schema("output.schema.json")).validate(d)
            Draft202012Validator(_schema("interpretability.schema.json")).validate(
                d["interpretability"]
            )

            self.assertEqual(d["verdict"], "insufficient_evidence")
            self.assertEqual(d["verdict_basis"], "none")

            # A NOT_COMPUTED_LOCALLY limitation, and no fabricated numbers.
            codes = [lim["code"] for lim in d["interpretability"]["limitations"]]
            self.assertIn("NOT_COMPUTED_LOCALLY", codes)
            self.assertIsNone(d["target_precedent"]["distinct_actives"])
            self.assertIsNone(d["target_precedent"]["best_potency_nm"])
            self.assertIsNone(d["target_precedent"]["approved_small_molecules_count"])
            self.assertEqual(d["structure"]["tier"], "none")
            # Every axis is recorded as null-with-reason, not as 0/[] elsewhere.
            self.assertTrue(len(d["not_found"]) >= 1)


class MalformedInput(unittest.TestCase):
    def test_missing_uniprot_accession_fails_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad_input = Path(tmp) / "bad.json"
            bad_input.write_text(json.dumps({"as_of_date": None}))  # no uniprot_accession
            out_path = Path(tmp) / "o.json"
            # run_pipeline must never be reached on invalid input.
            with mock.patch.object(
                cli, "run_pipeline", side_effect=AssertionError("pipeline reached on bad input")
            ):
                code, stdout, stderr = _run_cli(
                    ["run", "--input", str(bad_input), "--output", str(out_path)]
                )
            self.assertNotEqual(code, cli.EXIT_OK)
            self.assertFalse(out_path.exists(), "nothing must be written on invalid input")
            self.assertIn("input.schema.json", stderr)


class PipelineRaises(unittest.TestCase):
    """The non-default agent path fails loudly and writes nothing (monkeypatched)."""

    def test_pipeline_failure_is_loud_and_writes_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "o.json"
            with mock.patch.object(
                cli,
                "run_pipeline",
                side_effect=PipelineUnavailableError("ANTHROPIC_API_KEY is not set"),
            ):
                code, stdout, stderr = _run_cli(
                    ["run", "--input", str(EXAMPLES / "input.json"), "--output", str(out_path)]
                )
            self.assertNotEqual(code, cli.EXIT_OK)
            self.assertFalse(out_path.exists(), "no dossier means nothing to write")
            self.assertIn("failed", stderr.lower())
            self.assertIn("ANTHROPIC_API_KEY", stderr)


class InterpretabilityValidates(unittest.TestCase):
    def _assert_valid(self, dossier_path: Path):
        dossier = json.loads(dossier_path.read_text())
        dossier.pop("interpretability", None)  # strip any pre-existing block first
        interp = build_interpretability(dossier)
        Draft202012Validator(_schema("interpretability.schema.json")).validate(interp)
        # Sanity: the two evidence axes are preserved separately under extensions,
        # never merged into one scalar (LABrador contract, section E).
        axes = interp["extensions"]["axes"]
        self.assertIn("retrieved_precedent", axes)
        self.assertIn("computed_tractability", axes)

    def test_example_output(self):
        self._assert_valid(EXAMPLES / "output.json")

    def test_integration_fixture_orthosteric(self):
        self._assert_valid(FIXTURES / "irak4_Q9NWZ3_orthosteric.json")

    def test_integration_fixture_oligomer_destabilisation(self):
        self._assert_valid(FIXTURES / "irak4_Q9NWZ3_oligomer_destabilisation.json")


if __name__ == "__main__":
    unittest.main(verbosity=2)
