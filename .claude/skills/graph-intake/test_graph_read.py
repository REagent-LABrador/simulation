#!/usr/bin/env python3
"""Tests for graph_read.py. Stdlib only, no dependencies.

    python3 test_graph_read.py

Two of these were written to FAIL against the code as it stood on 2026-08-15 and
are the reason the ask machinery was touched at all:

- `AlreadyAskedTargetless` -- `already_asked()` matched on (verb, target), and
  `new_question` carries `target: null` by design. The first new_question ever
  issued against a graph therefore retired the verb for the life of that graph.
- `PostResolutionGate2` -- the post-resolution contradiction ask was exempted
  from gate 3 and not from gate 2, so the one ask type that has demonstrably
  worked could not fire against the `primary` row it corrects.

Fixtures built here carry `_fixture: true`, and `FixtureGuard` checks that
graph_read.py still refuses them without --allow-fixture.
"""

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE / "graph_read.py"

sys.path.insert(0, str(HERE))
import graph_read  # noqa: E402

# Fixtures owned by other directories. Read-only, and skipped rather than failed
# if this skill is unbundled from them.
DOSSIER = HERE.parent.parent.parent
ASKBACK = DOSSIER / "fixtures" / "upstream_graph_askback.json"
WORKED_ASK = HERE.parent / "ppi-hypothesis" / "fixtures" / "worked_ask.json"

# g_1a4f -- the only real upstream graph. Not a fixture we wrote, which is what
# makes it the one that can surprise us.
REAL = DOSSIER / "fixtures" / "upstream_graph_real.json"
RUN_INTAKE = HERE / "run_intake.py"
CALIBRATE = HERE / "calibrate.py"


def minimal_graph(rounds, stop_reason="max_papers"):
    """The smallest graph the reader accepts, plus whatever `rounds` we need."""
    return {
        "_fixture": True,
        "_fixture_note": "SYNTHETIC. Built in-process by test_graph_read.py.",
        "schema_version": "1.1",
        "graph_id": "g_test",
        "round": 1,
        "status": "ok",
        "things": [{"id": "t1", "name": "TL1A", "kind": "protein"}],
        "papers": [],
        "findings": [],
        "links": [],
        "gaps": [],
        "rounds": rounds,
        "coverage": {"depth": "deep", "truncated": True, "stop_reason": stop_reason},
    }


def gate(gates, name):
    return next(g for g in gates if g["gate"] == name)


def check(graph, ask):
    """check_ask(), tolerant of the pre-fix 2-tuple return.

    So that running this file against the code as it stood before the fix fails
    on the BEHAVIOUR under test rather than on the changed arity.
    """
    result = graph_read.check_ask(graph, ask)
    if len(result) == 2:
        return result[0], result[1], []
    return result


# Two DISTINCT new_question asks. Distinct in the only sense that matters to a
# routing team: they name different sources and would be answered by different
# searches.
Q1 = ("Does any primary report measure a direct binding constant for a small "
      "molecule against TL1A? PMC10762860 asserts it only as background and "
      "CHEMBL25 holds no such assay.")
Q2 = ("Is there a deposited co-structure of the TL1A ectodomain with any "
      "receptor other than DcR3? 3K51 and 3MI8 are the only ones we find, and "
      "PMC11642585 implies a third.")


class AlreadyAskedTargetless(unittest.TestCase):
    """FAILS BEFORE THE FIX. One new_question retired the verb permanently."""

    def test_second_distinct_new_question_is_not_already_asked(self):
        graph = minimal_graph([
            {"n": 1, "ask": "new_question", "target": None, "depth": "standard",
             "question": Q1, "outcome": "new_evidence"},
        ])
        gates, _, _ = check(
            graph, {"ask": "new_question", "target": None, "question": Q2})
        g = gate(gates, "NOT_ALREADY_ASKED")
        self.assertTrue(
            g["ok"],
            "a second, unrelated new_question was reported as already asked "
            "purely because an earlier one existed: " + g["detail"])

    def test_the_same_new_question_is_still_caught(self):
        """The gate must still do its job -- the fix is not 'always pass'."""
        graph = minimal_graph([
            {"n": 1, "ask": "new_question", "target": None, "question": Q1},
        ])
        gates, _, _ = check(
            graph, {"ask": "new_question", "target": None, "question": Q1})
        self.assertFalse(gate(gates, "NOT_ALREADY_ASKED")["ok"])

    def test_a_rephrased_question_still_counts_as_asked(self):
        """Upstream rewords a question when it services it. Identity is the set
        of source identifiers, not the wording -- a text hash would let this
        through as a new ask."""
        graph = minimal_graph([
            {"n": 1, "ask": "new_question", "target": None, "question": Q1},
        ])
        reworded = ("Has anyone measured a Kd or IC50 for any small molecule on "
                    "TL1A in a primary paper? We have CHEMBL25 (silent) and "
                    "PMC10762860 (background only).")
        self.assertNotEqual(reworded, Q1)
        gates, _, _ = check(
            graph, {"ask": "new_question", "target": None, "question": reworded})
        self.assertFalse(
            gate(gates, "NOT_ALREADY_ASKED")["ok"],
            "a rephrasing of an already-issued question was treated as new")

    def test_a_round_with_no_question_text_is_reported_not_ignored(self):
        """`rounds` is not required to record the question, and the real
        ask-back fixture's round 1 does not. That prior ask cannot be compared
        against, so the gate must SAY so rather than return a silent all-clear."""
        graph = minimal_graph([
            {"n": 1, "ask": "new_question", "target": None, "outcome": "new_evidence"},
        ])
        gates, _, _ = check(
            graph, {"ask": "new_question", "target": None, "question": Q2})
        g = gate(gates, "NOT_ALREADY_ASKED")
        self.assertTrue(g["ok"])
        self.assertIn("no `question` text", g["detail"])
        self.assertIn("[1]", g["detail"])

    def test_targeted_verbs_still_match_on_target(self):
        graph = minimal_graph([
            {"n": 1, "ask": "expand_node", "target": "t1", "depth": "deep"},
        ])
        hits, unmatchable = graph_read.already_asked(graph, "expand_node", "t1")
        self.assertEqual(len(hits), 1)
        self.assertEqual(unmatchable, [])
        hits, _ = graph_read.already_asked(graph, "expand_node", "t9")
        self.assertEqual(hits, [])


class PostResolutionGate2(unittest.TestCase):
    """FAILS BEFORE THE FIX. Gate 2 blocked the correction ask.

    Uses the worked TL1A case: `resolve_link` on L4, whose basis is `primary`,
    carrying a PDB census that contradicts the row.
    """

    @classmethod
    def setUpClass(cls):
        if not (ASKBACK.exists() and WORKED_ASK.exists()):
            raise unittest.SkipTest(f"need {ASKBACK} and {WORKED_ASK}")
        cls.graph = json.loads(ASKBACK.read_text())
        cls.ask = json.loads(WORKED_ASK.read_text())["ask"]

    def test_worked_ask_still_passes_all_five_mechanical_gates(self):
        gates, _, _ = check(self.graph, self.ask)
        self.assertEqual([g["gate"] for g in gates if not g["ok"]], [])
        self.assertEqual(len(gates), 5)

    def test_gate_2_is_exempted_for_the_post_resolution_ask(self):
        gates, unchecked, exempt = check(self.graph, self.ask)
        del gates
        joined_unchecked = " ".join(unchecked)
        joined_exempt = " ".join(exempt)
        self.assertIn("SUPPORT_IS_SECONDARY_ONLY", joined_exempt)
        self.assertIn("WE_TRIED_AND_FAILED", joined_exempt)
        self.assertNotIn(
            "SUPPORT_IS_SECONDARY_ONLY", joined_unchecked,
            "gate 2 is still reported as required for an ask that carries our "
            "own answer; that is what blocks the correction")
        self.assertIn("AFFECTS_THE_DOSSIER", joined_unchecked)

    def test_an_ordinary_ask_is_not_exempted(self):
        """The exemption must not leak. An ask that does not declare itself
        post-resolution still faces both judgment gates."""
        ordinary = dict(self.ask, question=self.ask["question"].replace(
            "THIS ASK IS POST-RESOLUTION AND NOT BLOCKING", "WE WOULD LIKE TO KNOW"))
        gates, unchecked, exempt = check(self.graph, ordinary)
        del gates
        self.assertEqual(exempt, [])
        self.assertIn("SUPPORT_IS_SECONDARY_ONLY", " ".join(unchecked))

    def test_ask_context_separates_L4_from_L2(self):
        """Both are `basis: primary` and both were `mechanical_gates_clear:
        false`. L2 is the ask that must never fire; L4 is the one that worked.
        The output has to tell them apart."""
        ctx = graph_read.ask_context(self.graph)
        rows = {r["link"]: r for r in ctx["links"]}
        self.assertFalse(rows["L4"]["mechanical_gates_clear"])
        self.assertTrue(rows["L4"]["clear_if_post_resolution_contradiction"])
        # L2 is primary too, so it is also correctable-in-principle; what keeps
        # it out is that we have not contradicted it. L7 is the useful contrast:
        # already asked in round 3, so nothing about it is clear.
        self.assertFalse(rows["L7"]["mechanical_gates_clear"])
        self.assertFalse(rows["L7"]["clear_if_post_resolution_contradiction"])
        # L1/L3 are the ordinary secondary-only asks and are unaffected.
        self.assertTrue(rows["L1"]["mechanical_gates_clear"])
        self.assertTrue(rows["L3"]["mechanical_gates_clear"])
        self.assertFalse(rows["L1"]["clear_if_post_resolution_contradiction"])


class FixtureGuard(unittest.TestCase):
    """`_fixture: true` must keep being refused without --allow-fixture."""

    def run_cli(self, graph, *args):
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as fh:
            json.dump(graph, fh)
            path = fh.name
        try:
            return subprocess.run([sys.executable, str(SCRIPT), path, *args],
                                  capture_output=True, text=True)
        finally:
            Path(path).unlink()

    def test_fixture_is_refused_without_the_flag(self):
        r = self.run_cli(minimal_graph([]))
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("_fixture: true", r.stderr)

    def test_fixture_is_accepted_with_the_flag(self):
        r = self.run_cli(minimal_graph([]), "--allow-fixture")
        self.assertEqual(r.returncode, 0, r.stderr)
        json.loads(r.stdout)

    def test_guard_fires_before_check_ask(self):
        """--check-ask must not be a way around the guard."""
        ask = json.dumps({"ask": "new_question", "target": None, "question": Q2})
        r = self.run_cli(minimal_graph([]), "--check-ask", ask)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("_fixture: true", r.stderr)

    def test_a_real_graph_carrying_no_fixture_key_passes(self):
        graph = minimal_graph([])
        del graph["_fixture"]
        r = self.run_cli(graph)
        self.assertEqual(r.returncode, 0, r.stderr)


class CliShape(unittest.TestCase):
    """The CLI contract the ppi-hypothesis skill calls."""

    @classmethod
    def setUpClass(cls):
        if not (ASKBACK.exists() and WORKED_ASK.exists()):
            raise unittest.SkipTest("need the ask-back fixtures")

    def check_ask_cli(self, ask):
        return subprocess.run(
            [sys.executable, str(SCRIPT), str(ASKBACK), "--allow-fixture",
             "--check-ask", json.dumps(ask)],
            capture_output=True, text=True)

    def test_worked_ask_exits_zero(self):
        ask = json.loads(WORKED_ASK.read_text())["ask"]
        r = self.check_ask_cli(ask)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        out = json.loads(r.stdout)
        self.assertTrue(all(g["ok"] for g in out["gates"]))
        self.assertTrue(out.get("exempt_for_this_ask"))

    def test_unrelated_new_question_exits_zero(self):
        """The regression, end to end: round 1 of this fixture is a
        new_question, and it used to fail every later one."""
        r = self.check_ask_cli({"ask": "new_question", "target": None,
                                "depth": "deep", "question": Q2})
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_a_fifth_verb_is_still_refused(self):
        r = self.check_ask_cli({"ask": "ask_nicely", "target": None,
                                "depth": "deep", "question": Q2})
        self.assertNotEqual(r.returncode, 0)


class BuildStillWorks(unittest.TestCase):
    """The non-ask half of the file, unchanged by this work but exercised so a
    regression in it is not silent."""

    def test_build_on_the_askback_fixture(self):
        if not ASKBACK.exists():
            self.skipTest("fixture absent")
        out = graph_read.build(json.loads(ASKBACK.read_text()))
        self.assertEqual(out["status"], "ok")
        self.assertEqual(out["integrity"]["dangling_ids"], [])
        self.assertEqual(out["coverage"]["literature_exhausted"], False)

    def test_absent_list_is_refused(self):
        graph = minimal_graph([])
        del graph["links"]
        with self.assertRaises(graph_read.GraphShapeError):
            graph_read.build(graph)

    def test_null_list_is_refused(self):
        graph = minimal_graph([])
        graph["gaps"] = None
        with self.assertRaises(graph_read.GraphShapeError):
            graph_read.build(graph)


class SymbolCandidatesOnTheRealGraph(unittest.TestCase):
    """The merge regression these tests exist to prevent.

    `graph_read.build()` stopped emitting `symbol_candidates`. `run_intake.py`
    subscripts that key directly, so the one-command graph-to-dossier path died
    on every input; `ask_emit.py` and `calibrate.py` read it with `.get()` and
    reported zero instead of failing, which is worse. Nothing caught it because
    nothing ran the pipeline end to end.

    g_1a4f is the case that matters. It has zero `protein` and zero `gene`
    nodes -- SKILL.md failure mode 12b -- so the kind-based route CANNOT
    nominate, on any round, for any target. IRAK4 and MyD88 exist in it only as
    substrings of `small_molecule` and `process` node names. Recovering them is
    exactly and only what this route does, so a zero here is the whole real-graph
    path being dead, not a graph with nothing in it.
    """

    @classmethod
    def setUpClass(cls):
        if not REAL.exists():
            raise unittest.SkipTest("need fixtures/upstream_graph_real.json")
        cls.graph = json.loads(REAL.read_text())
        cls.out = graph_read.build(cls.graph)

    def test_build_still_emits_the_key(self):
        """A hard subscript, because that is how run_intake.py reads it."""
        self.assertIn("symbol_candidates", self.out)
        self.assertIn("candidates", self.out["symbol_candidates"])

    def test_the_real_graph_has_no_entity_node_to_nominate(self):
        """Failure mode 12b, asserted rather than trusted."""
        kinds = {t.get("kind") for t in self.graph["things"]}
        self.assertEqual(kinds & graph_read.TARGET_KINDS, set())
        self.assertEqual(self.out["nominations"], [])

    def test_zero_nominations_and_non_zero_symbol_candidates(self):
        """Both halves of the apparent contradiction, in one assertion.

        'the answer is IRAK4/MYD88' and 'this graph yields zero nominations' are
        both true: they describe different routes.
        """
        self.assertEqual(len(self.out["nominations"]), 0)
        self.assertGreater(len(self.out["symbol_candidates"]["candidates"]), 0)

    def test_irak4_and_myd88_are_recovered_from_node_names(self):
        cands = self.out["symbol_candidates"]["candidates"]
        by_symbol = {}
        for c in cands:
            by_symbol.setdefault(graph_read.symbol_key(c["symbol"]), []).append(c)

        self.assertIn("IRAK4", by_symbol)
        irak4 = by_symbol["IRAK4"][0]
        # Typed small_molecule and named for the intervention. The kind is
        # precisely what this route must not trust.
        self.assertEqual(irak4["thing"], "t1")
        self.assertEqual(irak4["thing_kind"], "small_molecule")
        self.assertEqual(irak4["field"], "name")
        self.assertEqual(irak4["action"], "inhibition")

        self.assertIn("MYD88", by_symbol)
        # Reached from both t4 ("MyD88 dimerization inhibition") and t5.
        self.assertEqual({c["thing"] for c in by_symbol["MYD88"]}, {"t4", "t5"})

    def test_myd88_query_forms_carry_the_spelling_that_resolves(self):
        """`gene_name` in uniprot_v.proteins is case-sensitive: 'MyD88' returns
        no row and 'MYD88' returns Q99836. The graph spells it 'MyD88'. Without
        the upper-case form in query_forms the MYD88 half of the known-good
        answer is unreachable, so this is the assertion that pins it."""
        cands = self.out["symbol_candidates"]["candidates"]
        myd88 = next(c for c in cands if c["symbol"] == "MyD88")
        self.assertEqual(myd88["query_forms"], ["MyD88", "MYD88"])

    def test_the_three_way_phrase_stays_ambiguous(self):
        """'TLR/MyD88/NF-kB signalling axis' is three candidates, not one.
        Collapsing it is the failure this route exists to avoid."""
        sc = self.out["symbol_candidates"]
        self.assertIn("t5", sc["ambiguous_things"])
        t5 = {c["symbol"] for c in sc["candidates"] if c["thing"] == "t5"}
        self.assertEqual(t5, {"TLR", "MyD88", "NF-kB"})

    def test_compound_codes_are_not_proposed_as_genes(self):
        """ST2825, PF-06650833 and KIC-0101 are symbol-SHAPED. A run of four
        digits is what separates a compound code from a gene symbol."""
        proposed = {c["symbol"] for c in self.out["symbol_candidates"]["candidates"]}
        for code in ("ST2825", "PF-06650833", "KIC-0101"):
            self.assertNotIn(code, proposed)

    def test_nothing_here_leaks_into_nominations(self):
        """This route PROPOSES. An accession is never asserted by the reader."""
        for c in self.out["symbol_candidates"]["candidates"]:
            self.assertIsNone(c["verified"])
            self.assertIsNone(c["uniprot_accession"])


class IntakeEndToEnd(unittest.TestCase):
    """Runs the actual entry point against the actual real graph.

    The absence of this is why the regression survived: every script still
    IMPORTED cleanly, so an import-level check saw nothing. Only running the
    pipeline touches the bundle keys it reads.

    `--dry-run` on purpose -- it makes no network call, so this runs offline and
    deterministically, and it still exercises every bundle-shape read. The
    KeyError fired regardless of dry-run.
    """

    @classmethod
    def setUpClass(cls):
        if not (REAL.exists() and RUN_INTAKE.exists()):
            raise unittest.SkipTest("need run_intake.py and the real fixture")
        cls.proc = subprocess.run(
            [sys.executable, str(RUN_INTAKE), "--graph", str(REAL),
             "--dry-run", "--json"],
            capture_output=True, text=True)

    def test_pipeline_exits_zero(self):
        self.assertEqual(self.proc.returncode, 0,
                         self.proc.stdout[-3000:] + self.proc.stderr[-3000:])

    def test_pipeline_did_not_raise_on_a_bundle_key(self):
        """The exact failure: `bundle["symbol_candidates"]["candidates"]`."""
        self.assertNotIn("KeyError", self.proc.stderr)
        self.assertNotIn("Traceback", self.proc.stderr)

    def test_traversal_reports_the_symbols_it_recovered(self):
        report = json.loads(self.proc.stdout)
        tr = report["traversal"]
        # Zero nominations is correct and expected on this graph.
        self.assertEqual(tr["nominations_from_entity_nodes"], [])
        # A zero here is the regression, not a property of the graph.
        self.assertGreater(tr["symbol_candidates_proposed"], 0)

    def test_orphan_findings_survived_being_moved(self):
        """`orphan_findings` moved under `integrity` as
        `findings_referenced_by_no_link`. g_1a4f ships one (f6), so a zero here
        means a consumer is reading a key that no longer exists."""
        report = json.loads(self.proc.stdout)
        orphans = report["traversal"]["orphan_findings"]
        self.assertEqual([o["finding"] for o in orphans], ["f6"])

    def test_every_recovered_symbol_reaches_verification(self):
        report = json.loads(self.proc.stdout)
        symbols = {s["symbol"] for s in report["verification"]["symbols"]}
        self.assertIn("IRAK4", symbols)
        self.assertIn("MyD88", symbols)


class CalibrateReadsTheCurrentShape(unittest.TestCase):
    """calibrate.py exists to measure drift. Reading a key that no longer exists
    with .get() made it report a confident zero for every graph."""

    @classmethod
    def setUpClass(cls):
        if not (REAL.exists() and CALIBRATE.exists()):
            raise unittest.SkipTest("need calibrate.py and the real fixture")
        sys.path.insert(0, str(HERE))

    def test_summary_counts_are_not_silently_zero(self):
        import calibrate  # noqa: E402
        graph = json.loads(REAL.read_text())
        row = calibrate.survey(graph)
        self.assertEqual(row["graph_id"], "g_1a4f")
        self.assertEqual(row["nominations"], 0)
        self.assertGreater(row["symbol_candidates"], 0)
        self.assertIn("IRAK4", row["symbols_proposed"])
        self.assertEqual(row["orphan_findings"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
