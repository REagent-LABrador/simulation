"""Offline contract tests for the simulation module. No paid calls, no network.

Run standalone:
    /Users/bb/micromamba/envs/druggability/bin/python simulation/test_module.py
or:
    micromamba run -n druggability python -m unittest simulation.test_module

Every test that needs a dossier uses a RECORDED REAL RESULT as a test double:
examples/output.json is a real dossier the deployed agent produced, checked into
the repo. We monkeypatch run_pipeline to return it (or to raise) so the contract
mechanics -- input validation, output validation, exit codes, stdout cleanliness,
interpretability attachment -- are exercised WITHOUT invoking the paid agent.
"""

from __future__ import annotations

import copy
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


def _recorded_dossier() -> dict:
    """The recorded real result used as a test double, with any pre-existing
    interpretability block stripped so the runner regenerates it."""
    dossier = json.loads((EXAMPLES / "output.json").read_text())
    dossier.pop("interpretability", None)
    return dossier


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


class EndToEndSuccess(unittest.TestCase):
    def test_success_writes_valid_output_and_clean_stdout(self):
        dossier = _recorded_dossier()
        with tempfile.TemporaryDirectory() as tmp:
            out_path = Path(tmp) / "o.json"
            with mock.patch.object(cli, "run_pipeline", return_value=copy.deepcopy(dossier)):
                code, stdout, stderr = _run_cli(
                    ["run", "--input", str(EXAMPLES / "input.json"), "--output", str(out_path)]
                )
            self.assertEqual(code, cli.EXIT_OK, msg=stderr)
            self.assertEqual(stdout, "", "stdout must stay empty; results go to --output")
            self.assertTrue(out_path.exists(), "output file must be written")

            written = json.loads(out_path.read_text())
            Draft202012Validator(_schema("output.schema.json")).validate(written)
            self.assertIn("interpretability", written)
            self.assertEqual(written["interpretability"]["module"], "simulation")


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
        # Sanity: the two axes are present and never merged.
        self.assertEqual([a["id"] for a in interp["axes"]], ["retrieved_precedent", "computed_tractability"])

    def test_example_output(self):
        self._assert_valid(EXAMPLES / "output.json")

    def test_integration_fixture_orthosteric(self):
        self._assert_valid(FIXTURES / "irak4_Q9NWZ3_orthosteric.json")

    def test_integration_fixture_oligomer_destabilisation(self):
        self._assert_valid(FIXTURES / "irak4_Q9NWZ3_oligomer_destabilisation.json")


if __name__ == "__main__":
    unittest.main(verbosity=2)
