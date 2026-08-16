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
IN = Draft202012Validator(INPUT_SCHEMA)
OUT = Draft202012Validator(OUTPUT_SCHEMA)


def _errs(validator, doc):
    return sorted(validator.iter_errors(doc), key=lambda e: list(e.absolute_path))


# --------------------------------------------------------------------------- #
# A minimal dossier that is valid by construction. Negative fixtures mutate a
# COPY of this, so each mutation isolates exactly one violated constraint.
# --------------------------------------------------------------------------- #
def minimal_valid_dossier() -> dict:
    return {
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
        files = sorted(EXAMPLES.glob("*.json"))
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
        self.assertEqual(schema_required, set(self.vd.REQUIRED_TOP_LEVEL))

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


if __name__ == "__main__":
    unittest.main(verbosity=2)
