"""Exercise the input/output JSON Schemas to death.

Four kinds of check, and the middle two are the ones that matter most:

1. VALIDITY      — both files are legal Draft 2020-12 schemas.
2. REALITY       — every real, validator-passing dossier in examples/ validates,
                   and a corpus of diverse hand-built valid dossiers validates.
3. NON-DRIFT     — the schema's enums and required-key set are IDENTICAL to the
                   validator's (`validate_dossier.py`) and consistent with the
                   CLAUDE.md output template's pipe-string placeholders. A schema
                   that quietly disagrees with the semantic gate is worse than no
                   schema, so this is asserted, not assumed.
4. MUTATION      — for every constraint the schema claims to enforce, a fixture
                   that violates exactly that constraint is REJECTED. A schema
                   that never says no proves nothing.

Offline, stdlib + jsonschema only. Runs standalone (unittest) so it needs no
pytest:  python schema/test_schema.py
"""

from __future__ import annotations

import copy
import json
import re
import sys
import unittest
from pathlib import Path

from jsonschema import Draft202012Validator

HERE = Path(__file__).resolve().parent
STATION = HERE.parent
EXAMPLES = STATION / ".claude" / "skills" / "assemble-dossier" / "examples"
VALIDATOR_DIR = STATION / ".claude" / "skills" / "assemble-dossier"
CLAUDE_MD = STATION / "CLAUDE.md"

INPUT_SCHEMA = json.loads((HERE / "input.schema.json").read_text())
OUTPUT_SCHEMA = json.loads((HERE / "output.schema.json").read_text())
INTERP_SCHEMA = json.loads((HERE / "interpretability.schema.json").read_text())
IN = Draft202012Validator(INPUT_SCHEMA)
OUT = Draft202012Validator(OUTPUT_SCHEMA)
INTERP = Draft202012Validator(INTERP_SCHEMA)
CONTRACT_EXAMPLES = STATION / "examples"  # module-contract examples/input.json + output.json

# The deterministic interpretability builder — the same one run_pipeline uses.
sys.path.insert(0, str(STATION / "simulation"))
from interpretability import build_interpretability  # noqa: E402


def _errs(validator, doc):
    return sorted(validator.iter_errors(doc), key=lambda e: list(e.absolute_path))


# --------------------------------------------------------------------------- #
# A minimal dossier that is valid by construction. Negative fixtures mutate a
# COPY of this, so each mutation isolates exactly one violated constraint.
# --------------------------------------------------------------------------- #
def minimal_valid_dossier() -> dict:
    d = {
        "input": {
            "uniprot_accession": "P01375",
            "as_of_date": None,
            "disease_context": None,
            "interaction_to_disrupt": None,
            "mechanism_hypothesis": None,
        },
        "target": {"uniprot_accession": "P01375"},
        "as_of_date": None,
        "verdict": "insufficient_evidence",
        "verdict_basis": "none",
        "axis_conflict": None,
        "target_precedent": {},
        "biologic_precedent": {"note": "Presence of a biologic is validation, not tractability."},
        "family_precedent": {},
        "structural_neighbour_precedent": {},
        "pocket_neighbour_precedent": {},
        "structure": {"tier": "none"},
        "tractability": {
            "cryptic_pocket_risk": "undetermined",
            "cryptic_mechanism": "undetermined",
            "cryptic_potency_prior": {"expected_ceiling": "unknown"},
            "pocket_vs_interface": {"classification": "no_partner_structure"},
            "mdpocket_site_definition_used": "none",
        },
        "affinity": {},
        "falsification": {"checks_run": ["at least one check"], "findings": [], "survived": None},
        "next_experiment": {"description": "do a thing", "rationale": "because", "resolves": "the question"},
        "not_found": [],
    }
    # interpretability is now REQUIRED in the output schema; attach the real
    # builder's output so every consumer of this helper stays schema-valid.
    d["interpretability"] = build_interpretability(d)
    return d


class Validity(unittest.TestCase):
    def test_both_schemas_are_valid_draft_2020_12(self):
        Draft202012Validator.check_schema(INPUT_SCHEMA)
        Draft202012Validator.check_schema(OUTPUT_SCHEMA)

    def test_the_minimal_dossier_is_valid(self):
        self.assertEqual(_errs(OUT, minimal_valid_dossier()), [])


class Reality(unittest.TestCase):
    """Every real dossier validates. This is the check that catches an
    over-strict schema — the one that rejects correct work."""

    def test_every_example_validates(self):
        # Both the curated examples AND the dated integration fixtures — the
        # latter are real multi-structure runs (IRAK4, two mechanism modes) and
        # are what exposed chains_used being an object, not a string/array.
        files = sorted(EXAMPLES.glob("*.json")) + sorted(EXAMPLES.glob("integration-fixtures/*.json"))
        if not files:
            self.skipTest("no examples/ in this checkout")
        for f in files:
            with self.subTest(example=f.name):
                self.assertEqual(_errs(OUT, json.loads(f.read_text())), [])

    def test_diverse_valid_variants(self):
        # Each variant flexes a different legitimately-varying region.
        base = minimal_valid_dossier()
        variants = []

        v = copy.deepcopy(base)  # a positive small-molecule verdict on both axes
        v["verdict"], v["verdict_basis"] = "small_molecule_tractable", "both"
        v["target_precedent"] = {
            "distinct_actives": 900, "best_potency_nm": 0.4,
            "approved_small_molecules": [{"name": "tofacitinib", "year": 2012, "modality": "small_molecule", "source": "ChEMBL"}],
            "sources": ["chembl_v"],
        }
        v["structure"] = {"tier": "holo_experimental", "pdb_id": "6OIM", "holo_count": 300, "apo_count": 221}
        variants.append(("positive_both_axes", v))

        v = copy.deepcopy(base)  # clustering_d as a LIST, selection as a LIST
        v["tractability"]["pocket_volume_a3"] = {
            "primary_d1_6_a3": 320.5, "clustering_d": [1.6, 2.4],
            "site_pocket_selected_by": ["ligand_site_jaccard", "site_signature_overlap"],
        }
        variants.append(("list_valued_fields", v))

        v = copy.deepcopy(base)  # a PREDICTED structure (the GPU fallback path)
        v["structure"] = {"tier": "predicted", "predicted_plddt": 71.2, "pdb_id": None}
        v["tractability"]["caveat"] = "pocket scan ran on an ESMFold model; see warnings"
        variants.append(("predicted_structure_fallback", v))

        v = copy.deepcopy(base)  # a not_found carrying a failure signature
        v["not_found"] = [
            "target_precedent.distinct_actives: [error] Request timed out",
            {"field": "structure", "reason": "no experimental structure", "signature": "none"},
        ]
        variants.append(("not_found_shapes", v))

        v = copy.deepcopy(base)  # the optional follow-up channel populated
        v["follow_up_questions"] = [
            {"ask": "resolve_link", "target": "L4", "question": "no such structure exists; verify", "depth": "deep"}
        ]
        variants.append(("optional_follow_up_questions", v))

        for name, doc in variants:
            with self.subTest(variant=name):
                self.assertEqual(_errs(OUT, doc), [])


class NonDrift(unittest.TestCase):
    """The schema must not disagree with the semantic gate or the template."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(VALIDATOR_DIR))
        try:
            import validate_dossier as vd  # noqa: PLC0415
        except Exception as exc:  # noqa: BLE001
            cls.vd = None
            cls.reason = f"validate_dossier not importable: {exc}"
        else:
            cls.vd = vd

    def _need_vd(self):
        if self.vd is None:
            self.skipTest(self.reason)

    def test_required_top_level_matches_validator(self):
        self._need_vd()
        schema_required = set(OUTPUT_SCHEMA["required"])
        # `interpretability` is deterministically builder-attached at output time
        # and required by the PUBLISHED OUTPUT SCHEMA (LABrador contract), but it
        # is deliberately NOT in validate_dossier's authored-body
        # REQUIRED_TOP_LEVEL — that set also governs the hand-filled CLAUDE.md
        # template, which the builder, not the author, populates. Every other
        # required top-level key must still match the validator exactly.
        self.assertIn("interpretability", schema_required)
        self.assertEqual(schema_required - {"interpretability"}, set(self.vd.REQUIRED_TOP_LEVEL))

    def _enum_at(self, path):
        node = OUTPUT_SCHEMA["properties"]
        parts = path.split(".")
        for i, p in enumerate(parts):
            node = node[p]
            if i < len(parts) - 1:
                node = node["properties"]
        return set(node["enum"])

    def test_every_validator_enum_matches_the_schema(self):
        self._need_vd()
        # map validator ENUMS paths to the schema location
        for path, allowed in self.vd.ENUMS.items():
            with self.subTest(enum=path):
                self.assertEqual(
                    self._enum_at(path), set(allowed),
                    f"schema enum at {path} disagrees with validator",
                )

    def test_selection_basis_seven_values(self):
        vals = set(OUTPUT_SCHEMA["$defs"]["selectionValue"]["enum"])
        self.assertEqual(len(vals), 7)
        self.assertIn("ligand_site_jaccard", vals)

    def test_template_pipe_strings_match_schema_enums(self):
        # The CLAUDE.md template writes an unfilled enum as "a | b | c". Those
        # placeholders must be exactly the schema's enum members.
        if not CLAUDE_MD.exists():
            self.skipTest("no CLAUDE.md")
        text = CLAUDE_MD.read_text()
        m = re.search(r"## Output template.*?```json\n(.*?)\n```", text, re.S)
        self.assertIsNotNone(m, "could not find the json template fence")
        raw = m.group(1)
        for key, schema_enum in (
            ("verdict", set(OUTPUT_SCHEMA["properties"]["verdict"]["enum"])),
            ("verdict_basis", set(OUTPUT_SCHEMA["properties"]["verdict_basis"]["enum"])),
        ):
            hit = re.search(rf'"{key}":\s*"([^"]+)"', raw)
            self.assertIsNotNone(hit, f"{key} not found in template")
            placeholder = {p.strip() for p in hit.group(1).split("|")}
            self.assertEqual(placeholder, schema_enum, f"{key} template/schema drift")


class Mutation(unittest.TestCase):
    """Every claimed constraint rejects a fixture that violates exactly it."""

    def assertRejected(self, doc, *, where=None):
        errs = _errs(OUT, doc)
        self.assertTrue(errs, "expected the schema to REJECT this fixture")
        if where is not None:
            paths = {tuple(e.absolute_path) for e in errs}
            self.assertTrue(
                any(where == p or (len(p) >= len(where) and p[: len(where)] == where) for p in paths),
                f"expected a violation at {where}, got {sorted(paths)}",
            )

    def test_missing_each_required_top_level_key_is_rejected(self):
        for key in OUTPUT_SCHEMA["required"]:
            d = minimal_valid_dossier()
            d.pop(key)
            with self.subTest(missing=key):
                self.assertRejected(d)

    def test_unknown_top_level_key_is_rejected(self):
        # additionalProperties:false reports at the parent (root) path, not at
        # the offending key — so assert rejection and that the message names it.
        d = minimal_valid_dossier()
        d["surprise"] = 1
        errs = _errs(OUT, d)
        self.assertTrue(errs)
        self.assertTrue(any("surprise" in e.message for e in errs))

    def test_null_enum_is_rejected(self):
        d = minimal_valid_dossier()
        d["verdict"] = None
        self.assertRejected(d, where=("verdict",))

    def test_bad_enum_value_is_rejected(self):
        for path, bad in (
            (("verdict",), "maybe_tractable"),
            (("verdict_basis",), "averaged"),
            (("structure", "tier"), "guessed"),
            (("tractability", "cryptic_pocket_risk"), "extreme"),
            (("tractability", "pocket_vs_interface", "classification"), "orthosteric"),
            (("tractability", "mdpocket_site_definition_used"), "site_from_vibes"),
        ):
            d = minimal_valid_dossier()
            node = d
            for p in path[:-1]:
                node = node[p]
            node[path[-1]] = bad
            with self.subTest(path=path):
                self.assertRejected(d, where=path)

    def test_wrong_type_is_rejected(self):
        d = minimal_valid_dossier()
        d["target_precedent"]["distinct_actives"] = "lots"  # want int|null
        self.assertRejected(d, where=("target_precedent", "distinct_actives"))

    def test_negative_count_is_rejected(self):
        d = minimal_valid_dossier()
        d["target_precedent"]["distinct_actives"] = -3
        self.assertRejected(d, where=("target_precedent", "distinct_actives"))

    def test_fraction_out_of_range_is_rejected(self):
        d = minimal_valid_dossier()
        d["tractability"]["ligand_site_jaccard"] = 1.7  # 0..1
        self.assertRejected(d, where=("tractability", "ligand_site_jaccard"))

    def test_bad_selection_basis_in_list_is_rejected(self):
        d = minimal_valid_dossier()
        d["tractability"]["pocket_druggability"] = {"site_pocket_selected_by": ["ligand_site_jaccard", "nonsense"]}
        self.assertRejected(d, where=("tractability", "pocket_druggability", "site_pocket_selected_by"))

    def test_empty_falsification_checks_is_rejected(self):
        d = minimal_valid_dossier()
        d["falsification"]["checks_run"] = []  # minItems 1
        self.assertRejected(d, where=("falsification", "checks_run"))

    def test_empty_next_experiment_description_is_rejected(self):
        d = minimal_valid_dossier()
        d["next_experiment"]["description"] = ""  # minLength 1
        self.assertRejected(d, where=("next_experiment", "description"))

    def test_bad_follow_up_verb_is_rejected(self):
        d = minimal_valid_dossier()
        d["follow_up_questions"] = [{"ask": "delete_node", "target": "n1"}]
        self.assertRejected(d, where=("follow_up_questions", 0, "ask"))


class InputSchema(unittest.TestCase):
    def test_minimal_request_is_valid(self):
        self.assertEqual(_errs(IN, {"uniprot_accession": "Q9NWZ3"}), [])

    def test_full_request_is_valid(self):
        self.assertEqual(_errs(IN, {
            "uniprot_accession": "P01375",
            "as_of_date": "2018-01-01",
            "disease_context": "IBD",
            "interaction_to_disrupt": "TNFR2 engagement",
            "mechanism_hypothesis": "oligomer_destabilisation",
        }), [])

    def test_isoform_accession_is_valid(self):
        self.assertEqual(_errs(IN, {"uniprot_accession": "P23458-2"}), [])

    def test_missing_accession_is_rejected(self):
        self.assertTrue(_errs(IN, {"disease_context": "cancer"}))

    def test_gene_symbol_as_accession_is_rejected(self):
        # A gene symbol is exactly the mistake the pattern is here to catch.
        self.assertTrue(_errs(IN, {"uniprot_accession": "IRAK4"}))

    def test_unknown_request_field_is_rejected(self):
        self.assertTrue(_errs(IN, {"uniprot_accession": "P01375", "molecule_type": "small_molecule"}))

    def test_bad_mechanism_is_rejected(self):
        self.assertTrue(_errs(IN, {"uniprot_accession": "P01375", "mechanism_hypothesis": "covalent"}))

    def test_bad_date_is_rejected(self):
        self.assertTrue(_errs(IN, {"uniprot_accession": "P01375", "as_of_date": "August 2018"}))


class ModuleContract(unittest.TestCase):
    """The module-interface contract: examples/input.json + output.json exist,
    validate against their schemas, and the output carries a valid
    interpretability object."""

    def test_interpretability_schema_is_valid(self):
        Draft202012Validator.check_schema(INTERP_SCHEMA)

    def test_example_input_validates(self):
        f = CONTRACT_EXAMPLES / "input.json"
        if not f.exists():
            self.skipTest("no examples/input.json")
        self.assertEqual(_errs(IN, json.loads(f.read_text())), [])

    def test_example_output_validates(self):
        f = CONTRACT_EXAMPLES / "output.json"
        if not f.exists():
            self.skipTest("no examples/output.json")
        self.assertEqual(_errs(OUT, json.loads(f.read_text())), [])

    def test_example_output_carries_a_valid_interpretability_object(self):
        f = CONTRACT_EXAMPLES / "output.json"
        if not f.exists():
            self.skipTest("no examples/output.json")
        doc = json.loads(f.read_text())
        self.assertIn("interpretability", doc, "the example output should demonstrate the panel object")
        self.assertEqual(_errs(INTERP, doc["interpretability"]), [])


class InterpretabilityContract(unittest.TestCase):
    """The LABrador shared interpretability contract (message.txt, section E).
    These are the red tests the spec asks for."""

    def _built(self) -> dict:
        """Interpretability built by the real builder from the example dossier."""
        doc = json.loads((CONTRACT_EXAMPLES / "output.json").read_text())
        base = {k: v for k, v in doc.items() if k != "interpretability"}
        return build_interpretability(base)

    def test_interpretability_is_required_in_output_schema(self):
        self.assertIn("interpretability", OUTPUT_SCHEMA["required"])

    def test_deleting_interpretability_fails_output_validation(self):
        d = minimal_valid_dossier()
        d.pop("interpretability")
        self.assertTrue(_errs(OUT, d), "output schema must reject a dossier with no interpretability")

    def test_example_carries_a_deeply_valid_interpretability(self):
        doc = json.loads((CONTRACT_EXAMPLES / "output.json").read_text())
        self.assertIn("interpretability", doc)
        self.assertEqual(_errs(INTERP, doc["interpretability"]), [])

    def test_builder_output_is_deeply_valid(self):
        self.assertEqual(_errs(INTERP, self._built()), [])

    def test_builder_valid_for_every_fixture_including_abstaining(self):
        files = sorted(EXAMPLES.glob("integration-fixtures/*.json"))
        d = minimal_valid_dossier()
        d.pop("interpretability", None)
        cases = [("minimal_abstaining", d)] + [
            (f.name, {k: v for k, v in json.loads(f.read_text()).items() if k != "interpretability"})
            for f in files
        ]
        for name, dossier in cases:
            with self.subTest(case=name):
                self.assertEqual(_errs(INTERP, build_interpretability(dossier)), [])

    def test_ids_unique_and_all_references_resolve(self):
        i = self._built()
        for coll in ("metrics", "steps", "evidence", "assumptions"):
            ids = [x["id"] for x in i[coll]]
            self.assertEqual(len(ids), len(set(ids)), f"{coll} ids not unique")
        ev = {e["id"] for e in i["evidence"]}
        asm = {a["id"] for a in i["assumptions"]}
        met = {m["id"] for m in i["metrics"]}
        for m in i["metrics"]:
            for r in m["evidence_ids"]:
                self.assertIn(r, ev, f"metric {m['id']} references missing evidence {r}")
            for r in m["assumption_ids"]:
                self.assertIn(r, asm, f"metric {m['id']} references missing assumption {r}")
        for s in i["steps"]:
            for r in s["evidence_ids"]:
                self.assertIn(r, ev)
            for r in s["assumption_ids"]:
                self.assertIn(r, asm)
        for iv in i["uncertainty"]["intervals"]:
            self.assertIn(iv["metric_id"], met, "interval references a missing metric")

    def test_numeric_metrics_have_units(self):
        for m in self._built()["metrics"]:
            v = m["value"]
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                self.assertTrue(m.get("unit"), f"numeric metric {m['id']} has no unit")

    def test_unknown_actives_stays_null_and_earns_a_limitation(self):
        d = minimal_valid_dossier()
        d.pop("interpretability", None)
        d["target_precedent"] = {"distinct_actives": None}
        i = build_interpretability(d)
        codes = {limn["code"] for limn in i["limitations"]}
        self.assertIn("ACTIVES_NOT_RETRIEVED", codes)
        # never silently converted to 0 or an empty metric
        self.assertFalse(any(m["id"] == "metric.distinct_actives" for m in i["metrics"]))

    def test_no_ambiguous_none_actives_phrase(self):
        self.assertNotIn("None actives", json.dumps(self._built()))

    def test_two_axes_preserved_separately_in_extensions(self):
        ax = self._built()["extensions"]["axes"]
        self.assertIn("retrieved_precedent", ax)
        self.assertIn("computed_tractability", ax)
        # never collapsed into one invented scalar
        self.assertNotIn("combined_score", ax)

    def test_every_not_found_maps_to_a_limitation_with_field_path(self):
        d = minimal_valid_dossier()
        d.pop("interpretability", None)
        d["not_found"] = [{"field": "target_precedent.distinct_actives", "reason": "[error] Request timed out"}]
        i = build_interpretability(d)
        hits = [limn for limn in i["limitations"]
                if limn.get("field_path") == "output.target_precedent.distinct_actives"]
        self.assertTrue(hits, "not_found item did not produce a limitation with a field_path")

    def test_figure_is_packaged_or_flagged_absent(self):
        fig = self._built()["extensions"]["figure"]
        self.assertIn("present", fig)
        self.assertIn("note", fig)
        self.assertFalse(fig["present"], "figure is not packaged, so present must be false with a note")

    def test_cache_key_uses_full_input_not_just_accession(self):
        d = minimal_valid_dossier()
        d.pop("interpretability", None)
        d["input"] = {"uniprot_accession": "P00533", "as_of_date": None, "disease_context": None,
                      "interaction_to_disrupt": None, "mechanism_hypothesis": "orthosteric"}
        k1 = build_interpretability(d)["extensions"]["cache_key"]
        d2 = copy.deepcopy(d)
        d2["input"]["mechanism_hypothesis"] = "allosteric"  # same accession, different input
        k2 = build_interpretability(d2)["extensions"]["cache_key"]
        self.assertTrue(k1.startswith("sha256:"))
        self.assertNotEqual(k1, k2, "cache key must change when the input changes, not only the accession")

    def test_outputs_are_strict_json_without_nan_or_infinity(self):
        def reject(token):
            raise ValueError(f"non-finite JSON constant: {token}")
        files = [CONTRACT_EXAMPLES / "output.json"] + sorted(EXAMPLES.glob("integration-fixtures/*.json"))
        for f in files:
            with self.subTest(f=f.name):
                json.loads(f.read_text(), parse_constant=reject)  # raises on NaN/Infinity


if __name__ == "__main__":
    unittest.main(verbosity=2)
