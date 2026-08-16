#!/usr/bin/env python3
"""The template and the validator must not drift apart. This is the check.

`CLAUDE.md`'s output template is the shape an agent fills. `validate_dossier.py`
is the gate that shape has to pass. They are two files, edited by two different
kinds of change, and **they have now disagreed three times in one day** — most
recently when the rule 4 volume re-prioritisation added
`pocket_volume_a3.clustering_d`, `pocket_volume_a3.primary_d1_6_a3`,
`pocket_druggability.load_bearing` and `pocket_druggability._false_negative_rate`
to the validator's requirements and not to the template. An agent filling the
template exactly as written produced a dossier the gate rejected, and the only
signal was a live run failing.

Nothing detected any of those. Every existing test in `test_validate_dossier.py`
builds its dossiers from `examples/*.json`, which are hand-maintained real runs —
so they get updated whenever the validator does, and the template never enters
the loop at all. The template is the one artifact the suite never read.

So this file reads it. Two independent checks, because they fail differently:

1. `test_the_filled_template_passes_the_validator` — parse the literal JSON out
   of `CLAUDE.md`, fill every leaf with a plausible value, run the gate. Zero
   violations, or the template does not describe a dossier that can pass. This
   is the behavioural check and it catches requirements that live in rule
   *bodies* (`drug.get("load_bearing")`) rather than in a module constant.

2. `test_every_validator_path_exists_in_the_template` — a static check. Pull
   every dotted path the validator names, from its constants and from the
   `_get(d, "...")` literals in its source, and assert the template has it. This
   one names the missing key directly instead of making you read a violation.

Both are hostile to silence. A key that genuinely must exist on only one side is
allowed, but it has to be written into `EXEMPTIONS` with a reason.

**On the filler.** `FILL` below is a path-keyed table, and it encodes *semantics
only* — that `distinct_actives` must sit above the insufficient-evidence
threshold, that a cryptic census must be internally consistent, that a site
basis must be one the tool actually returns. It deliberately encodes no
structure: a key added to the template that needs no special value is filled by
`_default_for` and needs no entry here. That is what keeps this test from
becoming a second copy of the template that can itself drift.

Pure stdlib. Runs with the rest of the suite:

    python3 -m unittest discover -s . -p 'test_*.py' -v
"""

from __future__ import annotations

import copy
import json
import re
import unittest
from pathlib import Path

import validate_dossier as m
from validate_dossier import validate_dossier

HERE = Path(__file__).resolve().parent
CLAUDE_MD = HERE.parents[2] / "CLAUDE.md"

# A leaf that is deliberately absent on one side, FOREVER, and why. Empty is the
# goal. An entry here is a standing decision, not a to-do.
EXEMPTIONS: dict[str, str] = {
    "drug_name": (
        "validator-only alias. `_names()` accepts `name`, `drug_name` or "
        "`program` so it can read a hand-written or fixture-derived entry; the "
        "template offers exactly one of them, `name`, on purpose. Being liberal "
        "in what a checker accepts is not a promise to document all of it."
    ),
}

# A leaf that is absent on one side RIGHT NOW because the patch closing it has
# been written and not yet applied. Distinct from EXEMPTIONS on purpose: an
# entry here is a to-do with an owner, and
# `test_every_exemption_is_still_needed` fails the moment the patch lands, which
# is what forces the entry back out again. A permanent home for a temporary
# thing is how the last three drifts survived.
PENDING: dict[str, str] = {
    "classifications_seen": (
        "`check_mixed_interface_is_resolvable` requires it whenever "
        "`pocket_vs_interface.classification` is `mixed`, and `modal_app.py` "
        "emits it. CLAUDE.md is owned by another agent this session, so the "
        "template patch adding "
        "`tractability.pocket_vs_interface.classifications_seen: []` is written "
        "and routed rather than applied. Delete this entry when it lands."
    ),
}

# Strings that must never be produced by the filler: `check_null_is_not_zero`
# reads them as a failed measurement smuggled in as text.
_SENTINELS = {"n/a", "na", "none", "null", "unknown", "not_found", "-", "?"}

# A placeholder chosen so it contains no substring of any MEASURED_FIELDS path or
# leaf name. `_mentions_field` matches by substring, so a careless placeholder in
# `not_found` would silently excuse a null somewhere else.
_STR = "<supplied by test_template_drift>"

_BOOL_PREFIXES = ("is_", "has_", "measures_", "meets_", "placed_", "reproduces_", "matches_")
_BOOL_LEAVES = frozenset(
    {
        "trusted",
        "reliable",
        "survived",
        "leakage_risk",
        "load_bearing",
        "best_potency_characterised",
        "ligand_anchored",
    }
)
_NUM_SUFFIXES = (
    "_a3", "_nm", "_pct", "_a", "_count", "_score", "_fraction", "_length",
    "_plddt", "_overlap", "_similarity", "_jaccard", "_density", "_range",
    "_atoms", "_structures",
)
_NUM_LEAVES = frozenset(
    {
        "min", "max", "count", "distinct_actives", "family_actives", "year",
        "phase", "evalue", "confidence", "fpocket", "prank", "n_pockets",
        "binding_B", "functional_F", "n_structures", "n_measurements",
        "n_apo_examined", "n_apo_site_absent",
    }
)

# ---------------------------------------------------------------------------
# The semantic layer. Every entry answers "what makes this value COHERENT",
# never "what shape is this key" — shape comes from the template itself.
# ---------------------------------------------------------------------------

FILL: dict[str, object] = {
    # -- the verdict has to be self-consistent or half the rules fire --------
    # `retrieved_precedent` with a named approved small molecule below, so the
    # MODALITY_LEAK clause about biologic-only precedent stays quiet.
    "verdict": "small_molecule_tractable",
    "verdict_basis": "retrieved_precedent",
    # Left null on purpose: the cutoff path is exercised separately, by
    # `test_the_filled_template_passes_under_an_as_of_cutoff`.
    "as_of_date": None,
    # -- precedent: above the decline floor, below the conflict threshold ----
    "target_precedent.distinct_actives": 120,
    "target_precedent.assay_concentration.top_assay_share_pct": 4.8,
    "target_precedent.assay_concentration.measures_a_different_target": False,
    "target_precedent.best_potency_characterised": True,
    # No `-mab`/`-cept` stem, and no collision with the biologic below.
    "target_precedent.approved_small_molecules[0].name": "ruxolitinib",
    "target_precedent.clinical_stage_small_molecules[0].name": "itacitinib",
    "biologic_precedent.approved_biologics[0].name": "adalimumab",
    "biologic_precedent.approved_biologics[0].modality": "antibody",
    # -- structure: holo, so the cryptic-risk-from-tier clause is inert ------
    "structure.tier": "holo_experimental",
    "structure.holo_count": 42,
    "structure.apo_count": 10,
    "structure.total_pdb_structures": 52,
    # -- geometry: a two-structure, two-D pool over one identified site ------
    "tractability.method.ensemble_pdb_ids": ["3EYG", "1OPI"],
    "tractability.ensemble_consensus_fraction.n_structures": 2,
    "tractability.ensemble_consensus_fraction.n_measurements": 4,
    "tractability.ensemble_consensus_fraction.fraction_with_strong_pocket": 0.75,
    "tractability.pocket_volume_a3.min": 305.9,
    "tractability.pocket_volume_a3.max": 913.8,
    "tractability.pocket_volume_a3.spread_pct": 66.5,
    "tractability.pocket_volume_a3.clustering_d": [1.6, 2.4],
    # The primary number, and it must stay under the 1000 A^3 merge signature.
    "tractability.pocket_volume_a3.primary_d1_6_a3": 305.9,
    "tractability.pocket_volume_a3.site_pocket_selected_by": "ligand_site_jaccard",
    "tractability.pocket_druggability.min": 0.402,
    "tractability.pocket_druggability.max": 0.735,
    "tractability.pocket_druggability.fold_range": 1.8,
    "tractability.pocket_druggability.site_pocket_selected_by": "ligand_site_jaccard",
    # -- the cryptic census has to agree with itself and with the prognosis --
    "tractability.cryptic_pocket_risk": "high",
    "tractability.cryptic_mechanism": "loop_or_backbone_motion",
    "tractability.cryptic_potency_prior.expected_ceiling": "nanomolar",
    "tractability.cryptic_evidence.is_cryptic": True,
    "tractability.cryptic_evidence.n_apo_examined": 5,
    "tractability.cryptic_evidence.n_apo_site_absent": 5,
    "tractability.cryptic_evidence.site_present_in_apo_ensemble": False,
    # -- the interface class is measured, and it is not a disagreement -------
    "tractability.pocket_vs_interface.classification": "orthosteric_candidate",
    # NOTE: `classifications_seen` has no entry here yet, and deliberately so.
    # `check_mixed_interface_is_resolvable` needs it, but it is a pending
    # CLAUDE.md patch (that file is owned elsewhere) and
    # `test_every_FILL_key_is_a_real_template_path` refuses to let this table
    # reference a path the template does not have. When the patch lands the
    # generic fill covers it — a single-element list is not a disagreement, so
    # the mixed rule stays inert, which is correct for this dossier.
    "tractability.pocket_vs_interface.pocket_interface_overlap": 0.41,
    "tractability.pocket_vs_interface.partner_pdb_id": "3EYG",
    "tractability.pocket_vs_interface.matches_mechanism_hypothesis": True,
    # Quoted off the ligand site, well inside the proposed 4 A off-site bound.
    "tractability.mdpocket_site_definition_used": "site_from_ligand",
    "tractability.site_centroid_to_ligand_distance_a": 1.86,
    # Everything above is measured, so nothing is missing.
    "not_found": [],
}


# ---------------------------------------------------------------------------
# Reading the template out of CLAUDE.md
# ---------------------------------------------------------------------------


def load_template() -> dict:
    """The literal ```json block under `## Output template`.

    Parsed, not eyeballed. A template that is not valid JSON is a defect on its
    own — the agent is told to "fill this literally".
    """
    text = CLAUDE_MD.read_text(encoding="utf-8")
    head = text.index("## Output template")
    fence = re.search(r"```json\n(.*?)\n```", text[head:], re.DOTALL)
    if fence is None:
        raise AssertionError("no ```json fence found under '## Output template'")
    return json.loads(fence.group(1))


def _default_for(path: str) -> object:
    """The shape-only fallback, for any leaf `FILL` says nothing about.

    A new template key lands here and is filled without anybody editing this
    file, which is the whole point — the test must not need maintenance in
    lockstep with the thing it is checking.
    """
    leaf = path.rsplit(".", 1)[-1].split("[", 1)[0]
    low = leaf.lower()
    if leaf in _BOOL_LEAVES or low.startswith(_BOOL_PREFIXES):
        return True
    if leaf == "year":
        return 2015  # comfortably inside any cutoff the tests below apply
    if leaf == "phase":
        return 2
    if leaf in _NUM_LEAVES or low.endswith(_NUM_SUFFIXES):
        # Never 0. `check_null_is_not_zero` treats a measured zero as a RESULT,
        # so a zero here would be a claim rather than a filler value — and a
        # zero that also appears in `not_found` is claiming to be both.
        return 12.5
    return _STR


def _fill(node: object, path: str) -> object:
    if path in FILL:
        return copy.deepcopy(FILL[path])
    if isinstance(node, dict):
        return {k: _fill(v, f"{path}.{k}" if path else str(k)) for k, v in node.items()}
    if isinstance(node, list):
        if not node:
            return [_STR]
        if isinstance(node[0], (dict, list)):
            # The template shows list shape with one illustrative stub. Fill it.
            return [_fill(node[0], f"{path}[0]")]
        # A list of literals the template already committed to, e.g.
        # `clustering_d_swept: [1.6, 2.4]`. Keep every element — taking only the
        # first turned the mandatory two-value D sweep into a single D.
        return copy.deepcopy(node)
    if isinstance(node, str):
        if " | " in node:
            # An unfilled enum placeholder. Take the first alternative unless
            # FILL overrode it above — the first is a legal value by
            # construction, so an enum that GROWS needs no edit here.
            return node.split(" | ")[0].strip()
        if node.strip() == "":
            return _default_for(path)
        return node  # prose the template already supplies; keep it verbatim
    if node is None:
        return _default_for(path)
    return node  # a literal the template already committed to (false, [1.6, 2.4])


def fill_template(template: dict | None = None) -> dict:
    filled = _fill(copy.deepcopy(template or load_template()), "")
    assert isinstance(filled, dict)
    return filled


# ---------------------------------------------------------------------------
# Every dotted path the validator names
# ---------------------------------------------------------------------------

_STRING_LITERAL = re.compile(r'"([A-Za-z_][A-Za-z0-9_.]*)"')
_PATH_SEGMENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def validator_paths() -> set[str]:
    """From the module's constants and from every dotted literal in its source.

    Both halves matter. The constants carry the declarative requirements
    (`MEASURED_FIELDS`, `ENUMS`); the source scan catches paths that only a rule
    body mentions.

    The source scan is deliberately not `_get\\(d, "..."\\)`-shaped. Three of the
    four keys in the regression this file exists for are never passed to `_get`
    at all — `load_bearing` and `_false_negative_rate` are read with
    `drug.get(...)` off an already-fetched sub-dict, and `primary_d1_6_a3`
    reaches `_mentions_field`. What they DO all have is a fully-qualified dotted
    literal somewhere in the module, because every `Violation` names the path it
    is complaining about. So the scan takes any dotted string literal whose
    first segment is a known top-level key. That is a property of how this
    validator is written — a rule that reports a path a reader can look up — and
    it is a better hook than the call shape, which varies.

    Fragments produced by implicit string concatenation across lines (a literal
    ending in `.`) are dropped: they are half a path, and the other half is on
    the next source line.
    """
    paths: set[str] = set(m.REQUIRED_TOP_LEVEL)
    paths |= set(m.ENUMS)
    paths |= set(m.MEASURED_FIELDS)
    paths |= set(m.NOT_DATE_FILTERABLE)
    paths |= set(m.ALWAYS_FLAG_UNDER_CUTOFF)

    tops = set(m.REQUIRED_TOP_LEVEL)
    src = Path(m.__file__).read_text(encoding="utf-8")
    for lit in _STRING_LITERAL.findall(src):
        if "." not in lit:
            continue
        parts = lit.split(".")
        if parts[0] not in tops:
            continue
        if not all(_PATH_SEGMENT.match(p) for p in parts):
            continue  # a split fragment, not a whole path
        paths.add(lit)
    return paths


_BARE_GET = re.compile(r'\.get\(\s*"([A-Za-z_][A-Za-z0-9_]*)"')


def validator_leaf_reads() -> set[str]:
    """Leaf names the validator reads with a bare `.get("…")` on a sub-dict.

    `validator_paths()` cannot see these — there is no dotted path in the source
    to find, because the enclosing dict was fetched on a previous line. That is
    not a corner case: it is how **two of the four keys** in the regression this
    file exists for are read (`drug.get("load_bearing")`,
    `drug.get("_false_negative_rate")`), and it is how the validator's newest
    requirement is read too (`pvi.get("classifications_seen")`).

    So this half checks a weaker thing — that the leaf name exists *somewhere*
    in the template — and the weaker check is still decisive, because the names
    are specific enough that a collision is not plausible. A generic name like
    `min` or `name` matching some unrelated block is exactly the case where the
    behavioural half is the one doing the work anyway.
    """
    src = Path(m.__file__).read_text(encoding="utf-8")
    return set(_BARE_GET.findall(src))


def template_paths(node: object, path: str = "") -> set[str]:
    """Every dotted path in the template, list indices collapsed away."""
    out: set[str] = set()
    if isinstance(node, dict):
        for k, v in node.items():
            p = f"{path}.{k}" if path else str(k)
            out.add(p)
            out |= template_paths(v, p)
    elif isinstance(node, list):
        for v in node:
            out |= template_paths(v, path)
    return out


# ---------------------------------------------------------------------------
# The checks
# ---------------------------------------------------------------------------


class TestTemplateIsParseable(unittest.TestCase):
    def test_the_template_is_valid_json(self):
        """"Fill this literally" is not a coherent instruction otherwise."""
        self.assertIsInstance(load_template(), dict)

    def test_the_template_has_every_required_top_level_key(self):
        t = load_template()
        missing = [k for k in m.REQUIRED_TOP_LEVEL if k not in t]
        self.assertEqual(missing, [], f"REQUIRED_TOP_LEVEL keys not in template: {missing}")

    def test_every_FILL_key_is_a_real_template_path(self):
        """FILL shadows the template, so FILL is itself somewhere drift can hide.

        A `FILL` entry for a path the template no longer has is a semantic
        constraint being applied to nothing — the fill silently stops taking
        effect and the filled dossier quietly stops resembling a real one. Same
        failure shape as the dead `isinstance` guards in the validator: the code
        still reads as if it is doing something.
        """
        have = template_paths(load_template())
        # `a.b[0].c` in FILL addresses `a.b.c` once list indices collapse.
        stale = sorted(k for k in FILL if re.sub(r"\[\d+\]", "", k) not in have)
        self.assertEqual(stale, [], f"FILL keys with no template path: {stale}")

    def test_no_sentinel_strings_survive_the_fill(self):
        """Guards the filler itself, which is otherwise unchecked."""
        def walk(node, path=""):
            if isinstance(node, dict):
                for k, v in node.items():
                    walk(v, f"{path}.{k}" if path else str(k))
            elif isinstance(node, list):
                for v in node:
                    walk(v, path)
            elif isinstance(node, str):
                self.assertNotIn(node.strip().lower(), _SENTINELS, f"at {path}")

        walk(fill_template())


class TestTemplateMatchesValidator(unittest.TestCase):
    def test_every_validator_path_exists_in_the_template(self):
        """The static half. Names the missing key rather than a violation.

        This is the check that would have caught the rule 4 regression the
        moment it landed: `tractability.pocket_volume_a3.primary_d1_6_a3` is a
        `_get(d, "...")` literal in `check_volume_is_primary`, so it appears in
        `validator_paths()` whether or not anybody remembered the template.
        """
        have = template_paths(load_template())
        missing = sorted(p for p in validator_paths() if p not in have and p not in EXEMPTIONS)
        self.assertEqual(
            missing,
            [],
            "the validator requires paths the CLAUDE.md template does not have. "
            "Add them to the template, or add an entry to EXEMPTIONS with a "
            f"reason: {missing}",
        )

    def test_every_leaf_the_validator_reads_exists_in_the_template(self):
        """The half that closes `validator_paths()`'s blind spot.

        `check_druggability_not_load_bearing` reads `load_bearing` and
        `_false_negative_rate` off an already-fetched sub-dict, so no dotted
        literal names them and the path scan cannot see them. Before this check
        existed they were caught only behaviourally — a violation to read rather
        than a key to add.
        """
        have = template_paths(load_template())
        leaves = {p.rsplit(".", 1)[-1] for p in have}
        missing = sorted(
            n
            for n in validator_leaf_reads()
            if n not in leaves and n not in EXEMPTIONS and n not in PENDING
        )
        self.assertEqual(
            missing,
            [],
            "the validator reads these leaf names and the template has none of "
            f"them. Add them, or add an entry to EXEMPTIONS/PENDING: {missing}",
        )

    def test_every_exemption_is_still_needed(self):
        """An exemption that no longer applies is a lie the next reader inherits.

        For `PENDING` this is load-bearing rather than tidy: the entry is a
        to-do, and this assertion is what forces it out again the moment the
        patch it describes lands. Otherwise the exemption outlives the problem
        and quietly re-opens the hole it documented.
        """
        have = template_paths(load_template())
        leaves = {p.rsplit(".", 1)[-1] for p in have}
        # Recompute both scans with the tables ignored. An exemption is still
        # needed exactly when its key is something a scan would otherwise
        # report — asking "is the key in the template" instead conflates a
        # top-level PATH with a LEAF NAME of the same spelling, which
        # `as_of_leakage` is (exempted at the top level, present at
        # `target_precedent.as_of_leakage`).
        would_report = {p for p in validator_paths() if p not in have}
        would_report |= {n for n in validator_leaf_reads() if n not in leaves}
        for label, table in (("EXEMPTIONS", EXEMPTIONS), ("PENDING", PENDING)):
            for key, reason in table.items():
                with self.subTest(table=label, key=key):
                    self.assertTrue(reason.strip(), "an exemption must state a reason")
                    self.assertIn(
                        key,
                        would_report,
                        f"{key!r} no longer needs an exemption — delete its "
                        f"{label} entry",
                    )


class TestFilledTemplatePasses(unittest.TestCase):
    def test_the_filled_template_passes_the_validator(self):
        """The behavioural half, and the one that catches everything.

        An agent that fills the template exactly as written must produce a
        dossier the gate accepts. Any other outcome means the deployed prompt
        and the deployed grader disagree, and a live run is where you find out.
        """
        v = validate_dossier(fill_template())
        self.assertEqual(v, [], "\n".join(str(x) for x in v))

    def test_the_filled_template_passes_under_an_as_of_cutoff(self):
        """Rule 8's leakage registry is a template key; exercise it.

        With no cutoff the whole as-of rule is inert, so the base test above
        never touches `target_precedent.as_of_leakage` at all.
        """
        d = fill_template()
        d["as_of_date"] = "2030-01-01"
        d["target_precedent"]["as_of_leakage"] = [
            {"field": f, "leakage_risk": True, "note": "not date-filterable at source"}
            for f in (
                "clinical_stage_small_molecules",
                "distinct_actives",
                "best_potency_nm",
                "patents",
            )
        ]
        v = validate_dossier(d)
        self.assertEqual(v, [], "\n".join(str(x) for x in v))

    # -- the proof that this test is not vacuous ----------------------------

    def test_it_catches_the_rule_4_keys_going_missing(self):
        """CONFLICT 2, reconstructed. Remove the four keys; both halves fire.

        These are the keys the rule 4 volume re-prioritisation added to the
        validator. When the template lacked them an agent filling it exactly as
        written produced a dossier the gate rejected, and nothing in the suite
        noticed. Delete them from the template here and this file fails twice —
        once statically by name, once behaviourally by violation.
        """
        regressed = {
            ("tractability", "pocket_volume_a3"): ["clustering_d", "primary_d1_6_a3"],
            ("tractability", "pocket_druggability"): [
                "load_bearing",
                "_false_negative_rate",
            ],
        }
        t = load_template()
        for (block, sub), keys in regressed.items():
            for k in keys:
                self.assertIn(k, t[block][sub], f"{block}.{sub}.{k} left the template")
                del t[block][sub][k]

        # Static half: three of the four are `_get(d, "...")` literals or
        # constants, so they are named directly.
        have = template_paths(t)
        named = sorted(p for p in validator_paths() if p not in have and p not in EXEMPTIONS)
        self.assertIn("tractability.pocket_volume_a3.primary_d1_6_a3", named)
        self.assertIn("tractability.pocket_volume_a3.clustering_d", named)

        # Leaf half: `load_bearing` and `_false_negative_rate` are read with
        # `drug.get(...)` inside a rule body, so the PATH scan cannot see them.
        # This scan is what closes that hole — before it existed these two were
        # caught only behaviourally, as a violation to read rather than a key to
        # add.
        leaves = {p.rsplit(".", 1)[-1] for p in have}
        by_leaf = sorted(
            n
            for n in validator_leaf_reads()
            if n not in leaves and n not in EXEMPTIONS and n not in PENDING
        )
        self.assertIn("load_bearing", by_leaf)
        self.assertIn("_false_negative_rate", by_leaf)

        # Behavioural half: all four, as violations.
        v = validate_dossier(fill_template(t))
        fired = {x.path for x in v}
        self.assertIn("tractability.pocket_druggability.load_bearing", fired)
        self.assertIn("tractability.pocket_druggability._false_negative_rate", fired)
        self.assertIn("tractability.pocket_volume_a3.primary_d1_6_a3", fired)
        self.assertIn("tractability.pocket_volume_a3.clustering_d", fired)

    def test_it_catches_a_required_top_level_key_going_missing(self):
        """The coarsest drift, and the cheapest to introduce in a big edit."""
        t = load_template()
        del t["falsification"]
        self.assertIn("falsification", [p for p in validator_paths() if p not in template_paths(t)])
        self.assertIn("WELL_FORMED", {x.rule for x in validate_dossier(fill_template(t))})

    def test_it_catches_the_input_echo_block_going_missing(self):
        """The instance that landed WHILE this file was being written.

        The template grew a top-level `input` block echoing the five contract
        fields, and `REQUIRED_TOP_LEVEL` did not — so the block was documented
        and unenforced. Note which direction of drift each half sees:

        - **validator requires, template lacks** is the direction that breaks a
          live run, and `test_every_validator_path_exists_in_the_template`
          names the key. That is the direction conflict 2 had.
        - **template has, validator ignores** is not a failure on its own — an
          unenforced key is legal, and demanding the validator mention every
          template key would be demanding it grow a rule per field. What this
          suite gives instead is that the moment somebody DOES enforce it, the
          pair is locked together, which is what this test asserts.
        """
        t = load_template()
        self.assertIn("input", t, "the input echo block left the template")
        del t["input"]
        self.assertIn("input", [p for p in validator_paths() if p not in template_paths(t)])
        v = [
            x
            for x in validate_dossier(fill_template(t))
            if x.rule == "WELL_FORMED" and x.path == "input"
        ]
        self.assertTrue(v)

    def test_it_catches_an_enum_value_the_validator_rejects(self):
        """CONFLICT 1's shape: the template offers a value the gate refuses.

        The filler takes the first alternative of a pipe-separated placeholder,
        so a template that advertises an illegal value produces a dossier that
        fails — which is the drift, caught at the value level rather than the
        key level.
        """
        path = "tractability.pocket_vs_interface.classification"
        t = load_template()
        t["tractability"]["pocket_vs_interface"]["classification"] = (
            "consensus_candidate | orthosteric_candidate"
        )
        # FILL shadows the template for this path, so lift the override for the
        # duration — otherwise the illegal value never reaches the validator.
        saved = FILL.pop(path)
        try:
            v = [
                x
                for x in validate_dossier(fill_template(t))
                if x.rule == "WELL_FORMED" and x.path == path
            ]
        finally:
            FILL[path] = saved
        self.assertTrue(v, "an illegal template enum value was not caught")


if __name__ == "__main__":
    unittest.main(verbosity=2)
