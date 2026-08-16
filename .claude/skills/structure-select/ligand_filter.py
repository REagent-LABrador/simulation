"""
Ligand chemistry classifier — decide whether a PDB chemical component is a
DRUG-LIKE LIGAND THAT CONSTITUTES EVIDENCE OF A BINDABLE SITE.

WHY THIS EXISTS. Every "is this entry holo?" decision in the dossier pipeline
used to be made by (a) checking a hardcoded `comp_id` denylist of cofactors and
(b) a heavy-atom floor. Both halves fail, and they fail in the flattering
direction — they invent holo structures, which inflates apparent druggability.
Four measured wrong answers:

  * CD20 — 3 "holo" structures whose ligands were cholesterol hemisuccinate
    (`Y01`) and phosphatidylcholine, cryo-EM sample additives. True holo: 0.
  * KRAS fold neighbours — 4PHH's `2UK` is a GppNHp analog. A nucleotide
    cofactor whose comp_id nobody had listed. True holo: 0 of 25.
  * IL-17A fold neighbours — 4EC7's `L44` is a 625 Da diacylglycerol. It clears
    an 18-heavy-atom floor because it is a big greasy lipid. True holo: 0 of 25.
  * NLRP3 — binds ADP/ATP in the NACHT domain. ADP is 27 heavy atoms, so a pure
    size threshold calls an ADP-bound apo structure holo.

The pattern: **a hardcoded comp_id list cannot enumerate chemistry, and
molecular size does not distinguish a drug from a lipid.** Both a 625 Da
diacylglycerol and a 625 Da inhibitor clear a size gate; only chemistry
separates them. So this module classifies on chemistry — element composition,
ring topology, functional groups — read out of the PDB Chemical Component
Dictionary, which is the authority.

WHERE THE CHEMISTRY COMES FROM. `pdb_v.chemcomps` in Paperclip, which mirrors
the CCD and carries exactly the fields the old code ignored:

    comp_id, name, formula, formula_weight, type, smiles, inchikey, drugbank_id

`type` is `_chem_comp.type`. It alone settles the polymer question — `L-peptide
linking`, `D-saccharide, beta linking`, `RNA linking` are polymer residues and
are never small-molecule evidence. It is how 6OIM's GDP is caught: the CCD types
GDP as `RNA linking`, not `non-polymer`. Nothing outside Paperclip is fetched.

NO RDKIT, AND NO OTHER IMPORT OUTSIDE THE STANDARD LIBRARY. `pocket-scan`'s
Modal image has no RDKit — see its `image = modal.Image.micromamba(...)` block,
which installs fpocket, gemmi, numpy, torch and metapredict and nothing else —
and adding a chemistry toolkit to that image is not this module's call. So the
whole classifier runs on a self-contained SMILES graph parser in this file
(`SmilesGraph`): ~250 lines, deterministic, no network, no toolkit. The verdict
therefore cannot vary with the environment it is evaluated in, which for a
module three call sites depend on is worth more than RDKit's ring perception.
Everything here is 2D topology; no conformer is generated and no force field is
ever run, per the project's standing rule.

CHEMISTRY IS NOT ENOUGH, AND THE ENTRY SAYS SO. Added 2026-08-15 after the
first measured false positive against the held-out result. `LFI` in 8QFZ — the
TATA tri-electrophile that cyclises a Bicycle peptide — is drug-like chemistry
and is not a ligand: it is covalently bonded to all three cysteines of a
12-residue peptide that IS the ligand. Pass a `StructureContext` built from the
entry (or header) mmCIF and four rules re-attribute it; a covalent inhibitor
bonded to the TARGET is deliberately left alone. See `_apply_structure_context`
for the rules and for why this is a new verdict rather than a flag.

MEASURED ACCURACY, NOT A CLAIM. Every figure below re-measured 2026-08-15 with
all rules added since; the two originals are unchanged.

REPRODUCE THEM: `python3 tests/test_v2.py` from this directory prints every
figure below in one offline, stdlib-only pass. The harnesses and their cached
`pdb_v.chemcomps` rows are in `tests/` — see `tests/README.md`. Until
2026-08-15 they lived only in session scratch, which made every figure here
unverifiable from a checkout; if a number below ever disagrees with what the
harness prints, THE HARNESS WINS and the citation is what needs fixing.
  * 262-component ground-truth set (the four historical failures, every member
    of `modal_app.COFACTORS` and `NON_LIGANDS`, every member of
    `neighbour_precedent.EXCLUDED_LIGANDS`, known true-positive inhibitors,
    fragments, peptides, steroids, ions): **259/262 = 98.9%**. The three
    misses are `BTN` (biotin -> druglike), `ACE` and `NH2` (polymer capping
    groups -> additive / ion). `ACE`/`NH2` are CORRECT once a context is
    supplied — the CCD lists them in `_entity_poly_seq`, so they are residues.
  * The same set extended with 9 chemistry cases (`BEN`, `B3P`, `JEF` and the
    crosslinker/controls `LFI`, `ZBR`, `A1I4O`, `8VY`, `260`, `0WN`) and 9
    named-entry CONTEXT cases: **277/280 = 98.9%**, same three misses, no new
    ones. The context block is 9/9 and includes both covalent-inhibitor
    controls (6OIM `MOV`, 4G5J `0WN`), which must stay `druglike` and do.
  * 70-component HELD-OUT sample drawn blind from `pdb_v.chemcomps`:
    **61/70 = 87.1%**, and — the number that matters — **0 false positives**,
    unchanged by every rule added since. Nothing that was really a cofactor,
    lipid or additive was called drug-like. All 9 disagreements are the
    conservative direction. Two known boundaries on that zero: the sample
    contained neither TNF 5UUI's `MTN` spin label (still a false positive; one
    covalent bond to the target is rule C3, the covalent-inhibitor case) nor a
    peptide-conjugated crosslinker (`LFI`, now fixed by context).

KNOWN FALSE NEGATIVES, by class. Each was measured, not guessed. A false
negative costs a holo structure; a false positive invents one, which is what
produced all four bugs, so the bias is deliberate:
  * **Nucleoside and SAM analog inhibitors** (`V47`, `YB0`) -> `cofactor`. A
    purine plus a ribose is the nucleotide signature whether or not a medicinal
    chemist made it.
  * **Metallodrugs** (`U5U`, a palladacycle) -> `cofactor`, via the metal rule.
  * **Long-tailed natural-product antibiotics** (myxopyronin B) -> `lipid`.
  * **Glycosylated natural products** (abamectin) -> `sugar_or_glycan`.
  * **Peptidomimetic drugs** typed `peptide-like` (`LK0`, an HIV-protease
    inhibitor scaffold) -> `peptide_or_polymer`, correctly per the CCD type but
    wrongly as evidence.
  * **Biotin** -> `druglike`. The only member of the old `COFACTORS` set whose
    chemistry this does not recognise; catching it would require naming its
    ureido-thiophane bicycle, i.e. exactly the enumeration this module exists to
    replace.
  * **Bisphosphonate drugs and phosphate prodrugs** -> `cofactor`, by the
    phosphate-ester rule. Carries the flag
    `phosphate_rule_may_misfile_a_phosphate_drug`.
A call site that cares about any of these classes should read `evidence` and
`flags`, not just `verdict`.

VERDICTS
    druglike                 evidence of a bindable site
    cofactor                 nucleotide, flavin, heme, metal cluster, phosphate
                             metabolite — endogenous, not evidence of tractability
    lipid_or_detergent       acyl chains, sterols, phospholipids, detergents
    crystallisation_additive PEG, polyols, buffers, cryoprotectants
    sugar_or_glycan          saccharides and glycans
    ion_or_solvent           bare ions, simple inorganics, water
    peptide_or_polymer       polymer residues and peptide-like components
    polymer_conjugate        a covalent constituent of a polymer LIGAND — a
                             crosslinker, staple or macrocyclisation reagent.
                             NEEDS A `StructureContext`; `classify_record`
                             never returns it. `evidence["conjugate_of"]` names
                             the polymer so its precedent can be filed under
                             the right MODALITY (rule 1), which is the whole
                             point: the peptide is real evidence, just not
                             small-molecule evidence
    unknown                  the CCD has no record, or has no SMILES and the
                             non-structural fields do not decide it

`unknown` is a real answer and is NEVER coerced to `druglike`. A guess in this
code path is exactly what produced all four bugs.

USAGE

    from ligand_filter import classify_ligand, is_druglike_ligand, classify_ligands

    v = classify_ligand("L44")
    v.verdict        # 'lipid_or_detergent'
    v.reason         # 'longest aliphatic carbon chain is 21 (>= 8) ...'
    v.evidence       # every fact the verdict rests on

    is_druglike_ligand("2UK")                  # False
    classify_ligands(["MOV", "GDP", "MG"])     # one SQL round trip

    # Call sites that already hold CCD rows (e.g. from pdb_v.entry_ligands, or
    # from a gemmi mmCIF parse) should skip the network entirely:
    classify_record({"comp_id": "L44", "type": "non-polymer",
                     "formula": "C39 H76 O5", "formula_weight": 625.018,
                     "smiles": "CCCCCC...", "name": "..."})

    # Context-aware. THE HEADER, NOT THE ASSEMBLY — RCSB strips _struct_conn
    # from assembly files and the empty link table looks exactly like "nothing
    # is bonded to anything".
    ctx = StructureContext.from_mmcif_path("8QFZ_header.cif",
                                           target_accession="Q969D9")
    classify_ligand("LFI", context=ctx).verdict   # 'polymer_conjugate'
    holo_call(["LFI"], context=ctx)["is_holo"]    # False
    holo_call(["LFI"], context=ctx)["polymer_ligand_precedent"]
    # [{'modality': 'peptide', 'n_monomers': 12, 'sequence': 'CHWLENCWRGFC'...}]
"""

from __future__ import annotations

import json
import os
import re
import subprocess
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

__all__ = [
    "LigandVerdict",
    "VERDICTS",
    "CONTEXT_VERDICTS",
    "classify_ligand",
    "classify_ligands",
    "classify_record",
    "is_druglike_ligand",
    "filter_druglike",
    "holo_call",
    "SmilesGraph",
    "ChemCompSource",
    "StructureContext",
    "PolymerEntity",
    "CovalentLink",
    "read_mmcif_categories",
]

VERDICTS = (
    "druglike",
    "cofactor",
    "lipid_or_detergent",
    "crystallisation_additive",
    "sugar_or_glycan",
    "ion_or_solvent",
    "peptide_or_polymer",
    "polymer_conjugate",
    "unknown",
)

#: The verdicts that require a `StructureContext` and cannot be reached from
#: chemistry alone. `classify_record` never returns one of these.
CONTEXT_VERDICTS = ("polymer_conjugate",)

# --------------------------------------------------------------------------
# Tunables. Every one of these is a measured trade-off, not a preference.
# --------------------------------------------------------------------------

#: Longest unbranched sp3 carbon chain that marks an acyl/alkyl tail. 8 is the
#: shortest chain in the PDB detergent set that has to be caught (`OCT`, `C8E`,
#: octyl glucoside `BOG`); drug-like ligands essentially never carry one.
LIPID_MIN_CHAIN = 8

#: Below this, a component is bench chemistry, not a ligand. Deliberately well
#: under the 12-20 heavy atoms of a real fragment-screen hit — 5QQE's `N5S` is
#: 24 — so it cannot eat fragment precedent. It replaces nothing: the old
#: DRUGLIKE_MIN_HEAVY_ATOMS = 18 floor is gone, because size was never the
#: discriminator.
TRIVIAL_MAX_HEAVY_ATOMS = 9

#: Upper bound on a small molecule. Above this a `non-polymer` is a natural
#: product, a polymer or a macrocyclic peptide, not a small-molecule ligand.
DRUGLIKE_MAX_MW = 1200.0

#: The ubiquity prior (R14) fires only at or below this size. Additives and
#: fragment hits overlap between roughly 5 and 15 heavy atoms and nowhere else,
#: so bounding it here means the prior can never demote a real inhibitor. 15 is
#: above benzamidine's 9 and below the 24 of 5QQE's `N5S`.
UBIQUITY_MAX_HEAVY_ATOMS = 15

#: Entries a small component must appear in before ubiquity overrides
#: drug-like chemistry.
#:
#: MEASURED BLAST RADIUS, which is the number that makes this safe. Across the
#: whole 332-component union of the ground-truth and held-out sets, exactly
#: **11** components are `druglike` AND at or below
#: `UBIQUITY_MAX_HEAVY_ATOMS`, so those 11 are the only things this rule can
#: ever touch. Their counts, from RCSB (paperclip agrees where it answered):
#:
#:     BEN 361 | CFF 87 | 4NC 24 | LZ1 10 | ZBR 9 | 8VY 2
#:     260 1 | 363 1 | A1AX6 1 | KJM 1 | XLJ 1
#:
#: At 150 the rule fires on BEN and on nothing else, and the nearest thing it
#: does not touch is caffeine at 87 — 4.1x below. **UNCALIBRATED**: it is
#: fitted on one measured false positive and the margin is one case wide. It is
#: a proposal in exactly the sense rule 4a's volume guide is a proposal. Every
#: firing is flagged `ubiquity_prior_applied`, carries the count in `evidence`
#: and drops confidence to `medium`.
UBIQUITY_MIN_ENTRIES = 150

_METALS = frozenset(
    """LI BE NA MG AL K CA SC TI V CR MN FE CO NI CU ZN GA RB SR Y ZR NB MO TC
    RU RH PD AG CD IN SN SB CS BA LA CE PR ND PM SM EU GD TB DY HO ER TM YB LU
    HF TA W RE OS IR PT AU HG TL PB BI PO FR RA AC TH PA U NP PU AM CM
    AS SE TE""".split()
)

#: Advisory ONLY. Never touches the verdict. Literature-established promiscuous
#: binders / colloidal aggregators whose presence in a structure is weak
#: evidence even though the chemistry is drug-like. The falsification-sweep
#: skill owns promiscuity properly; this is a hint so a call site does not have
#: to rediscover the 2AZ5 problem. Keep it SHORT and cited.
_FREQUENT_HITTERS: dict[str, str] = {
    # Duan et al.; bis-electrophilic, widely regarded as promiscuous. The PDB
    # title says SPD304. See falsification-sweep/SKILL.md.
    "307": "2AZ5 TNF-alpha ligand SPD304: bis-electrophilic, widely reported "
           "as a promiscuous/aggregating binder",
}


# --------------------------------------------------------------------------
# Formula
# --------------------------------------------------------------------------

_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)\s*(\d*)")


def parse_formula(formula: str | None) -> dict[str, int]:
    """`'C39 H76 O5'` -> `{'C': 39, 'H': 76, 'O': 5}`. Charges are ignored.

    CCD formulae are whitespace-separated element-count pairs, sometimes with a
    trailing charge (`'O4 P 3-'`). Anything that is not an element token is
    dropped rather than guessed at.
    """
    if not formula:
        return {}
    out: dict[str, int] = {}
    for el, n in _FORMULA_TOKEN.findall(formula.replace("+", " ").replace("-", " ")):
        out[el.upper()] = out.get(el.upper(), 0) + (int(n) if n else 1)
    return out


def heavy_atom_count(formula: str | None) -> int | None:
    els = parse_formula(formula)
    if not els:
        return None
    return sum(n for el, n in els.items() if el not in ("H", "D", "T"))


# --------------------------------------------------------------------------
# SMILES -> graph. Self-contained; no chemistry toolkit.
# --------------------------------------------------------------------------

_BRACKET = re.compile(
    r"\[(\d*)([A-Za-z][a-z]?|\*)(@{0,2})(H\d*)?([+-]\d*|[+-]+)?(:\d+)?\]"
)
_TWO_CHAR = ("Cl", "Br")


@dataclass
class _Atom:
    idx: int
    element: str          # upper case, 'C', 'N', 'FE' ...
    aromatic: bool
    charge: int
    h_explicit: int | None
    in_ring: bool = False


class SmilesGraph:
    """A minimal SMILES parser producing atoms, bonds and small rings.

    Enough for the discriminations this module needs and nothing more: element
    composition, ring sizes and contents, ring fusion, longest aliphatic carbon
    chain, and a handful of functional-group queries. It does NOT do
    stereochemistry, aromaticity perception, valence checks or canonicalisation
    — the CCD SMILES is trusted as written, including its own kekulisation.

    Written rather than taken from RDKit because the Modal image that runs
    pocket-scan has no RDKit and the verdict must not depend on the environment.
    """

    def __init__(self, smiles: str):
        self.smiles = smiles
        self.atoms: list[_Atom] = []
        self.bonds: dict[tuple[int, int], float] = {}
        self.adj: dict[int, set[int]] = {}
        self.rings: list[list[int]] = []
        self.ok = False
        try:
            self._parse(smiles)
            self._find_rings()
            self.ok = True
        except Exception as exc:                      # noqa: BLE001
            self.parse_error = f"{type(exc).__name__}: {exc}"

    # -- parsing ---------------------------------------------------------

    def _add_atom(self, element: str, aromatic: bool, charge: int,
                  h_explicit: int | None) -> int:
        a = _Atom(len(self.atoms), element.upper(), aromatic, charge, h_explicit)
        self.atoms.append(a)
        self.adj[a.idx] = set()
        return a.idx

    def _add_bond(self, i: int, j: int, order: float) -> None:
        if i == j:
            return
        key = (min(i, j), max(i, j))
        self.bonds[key] = max(self.bonds.get(key, 0.0), order)
        self.adj[i].add(j)
        self.adj[j].add(i)

    def _parse(self, s: str) -> None:
        prev: list[int | None] = [None]
        pending_bond: float | None = None
        i = 0
        n = len(s)
        while i < n:
            ch = s[i]
            if ch == "[":
                m = _BRACKET.match(s, i)
                if not m:
                    raise ValueError(f"bad bracket atom at {i}: {s[i:i + 12]!r}")
                _iso, el, _chir, hs, chg, _map = m.groups()
                aromatic = el[0].islower()
                h_exp = None
                if hs:
                    h_exp = int(hs[1:]) if len(hs) > 1 else 1
                charge = 0
                if chg:
                    if chg[0] in "+-" and len(chg) > 1 and chg[1:].isdigit():
                        charge = int(chg[1:]) * (1 if chg[0] == "+" else -1)
                    else:
                        charge = len(chg) * (1 if chg[0] == "+" else -1)
                idx = self._add_atom(el, aromatic, charge, h_exp)
                if prev[-1] is not None:
                    self._add_bond(prev[-1], idx, pending_bond or 1.0)
                pending_bond = None
                prev[-1] = idx
                i = m.end()
                continue
            if ch == "(":
                prev.append(prev[-1])
                i += 1
                continue
            if ch == ")":
                prev.pop()
                i += 1
                continue
            if ch in "-=#$:/\\~":
                pending_bond = {"-": 1.0, "=": 2.0, "#": 3.0, "$": 4.0,
                                ":": 1.5, "/": 1.0, "\\": 1.0, "~": 1.0}[ch]
                i += 1
                continue
            if ch == ".":
                prev[-1] = None
                pending_bond = None
                i += 1
                continue
            if ch == "%":
                label = s[i + 1:i + 3]
                i += 3
                self._ring_closure(label, prev, pending_bond)
                pending_bond = None
                continue
            if ch.isdigit():
                self._ring_closure(ch, prev, pending_bond)
                pending_bond = None
                i += 1
                continue
            if s[i:i + 2] in _TWO_CHAR:
                idx = self._add_atom(s[i:i + 2], False, 0, None)
                if prev[-1] is not None:
                    self._add_bond(prev[-1], idx, pending_bond or 1.0)
                pending_bond = None
                prev[-1] = idx
                i += 2
                continue
            if ch.isalpha() or ch == "*":
                idx = self._add_atom(ch, ch.islower(), 0, None)
                if prev[-1] is not None:
                    self._add_bond(prev[-1], idx, pending_bond or 1.0)
                pending_bond = None
                prev[-1] = idx
                i += 1
                continue
            # Unrecognised (isotopes outside brackets, whitespace, junk).
            i += 1


    def _ring_closure(self, label: str, prev: list[int | None],
                      pending: float | None) -> None:
        if not hasattr(self, "_ring_map"):
            self._ring_map: dict[str, tuple[int, float | None]] = {}
        cur = prev[-1]
        if cur is None:
            return
        if label in self._ring_map:
            other, other_bond = self._ring_map.pop(label)
            self._add_bond(other, cur, pending or other_bond or 1.0)
        else:
            self._ring_map[label] = (cur, pending)

    # -- rings -----------------------------------------------------------

    def _find_rings(self, max_size: int = 9) -> None:
        """Smallest ring through each ring bond, by BFS with that bond removed.

        Not a formal SSSR — it can return a superset — but it recovers every
        ring needed here (purine 5+6, steroid 5+6+6+6, pyranose, porphyrin
        pyrroles) and it is deterministic.
        """
        seen: set[frozenset[int]] = set()
        for (a, b) in list(self.bonds):
            path = self._shortest_path(a, b, banned_edge=(a, b), max_len=max_size)
            if path is None:
                continue
            key = frozenset(path)
            if key in seen:
                continue
            seen.add(key)
            self.rings.append(path)
        for r in self.rings:
            for idx in r:
                self.atoms[idx].in_ring = True

    def _shortest_path(self, start: int, goal: int, *,
                       banned_edge: tuple[int, int],
                       max_len: int) -> list[int] | None:
        ban = (min(banned_edge), max(banned_edge))
        prev = {start: None}
        frontier = [start]
        depth = 0
        while frontier and depth < max_len:
            nxt = []
            for u in frontier:
                for v in self.adj[u]:
                    if (min(u, v), max(u, v)) == ban:
                        continue
                    if v in prev:
                        continue
                    prev[v] = u
                    if v == goal:
                        path, cur = [], v
                        while cur is not None:
                            path.append(cur)
                            cur = prev[cur]
                        return path
                    nxt.append(v)
            frontier = nxt
            depth += 1
        return None

    # -- queries ---------------------------------------------------------

    def element_counts(self) -> dict[str, int]:
        out: dict[str, int] = {}
        for a in self.atoms:
            out[a.element] = out.get(a.element, 0) + 1
        return out

    def neighbours(self, idx: int) -> list[_Atom]:
        return [self.atoms[j] for j in self.adj[idx]]

    def ring_elements(self, ring: Sequence[int]) -> dict[str, int]:
        out: dict[str, int] = {}
        for idx in ring:
            el = self.atoms[idx].element
            out[el] = out.get(el, 0) + 1
        return out

    def fused_groups(self) -> list[list[list[int]]]:
        """Rings grouped into fused systems (sharing >= 2 atoms)."""
        groups: list[list[list[int]]] = []
        used = [False] * len(self.rings)
        for i, r in enumerate(self.rings):
            if used[i]:
                continue
            used[i] = True
            group = [r]
            changed = True
            while changed:
                changed = False
                for j, r2 in enumerate(self.rings):
                    if used[j]:
                        continue
                    if any(len(set(r2) & set(m)) >= 2 for m in group):
                        used[j] = True
                        group.append(r2)
                        changed = True
            groups.append(group)
        return groups

    def longest_aliphatic_chain(self) -> int:
        """Longest simple path through acyclic, non-aromatic carbons.

        Each carbon on the path must be acyclic and bonded to no heteroatom
        other than as a chain terminus — i.e. a genuine hydrocarbon tail, which
        is what separates a lipid from a polar molecule that merely has a lot of
        carbons. `PEG` scores 2 because every second atom is an ether oxygen.
        """
        cand = [
            a.idx for a in self.atoms
            if a.element == "C" and not a.aromatic and not a.in_ring
            and sum(1 for nb in self.neighbours(a.idx)
                    if nb.element not in ("C", "H")) <= 1
        ]
        cset = set(cand)
        if not cset:
            return 0
        best = 0
        limit = 60000
        steps = 0

        def dfs(u: int, seen: set[int]) -> int:
            nonlocal steps, best
            steps += 1
            if steps > limit:
                return len(seen)
            local = len(seen)
            for v in self.adj[u]:
                if v in cset and v not in seen:
                    seen.add(v)
                    local = max(local, dfs(v, seen))
                    seen.discard(v)
            return local

        for s in cand:
            best = max(best, dfs(s, {s}))
            if steps > limit:
                break
        return best

    # -- functional groups ----------------------------------------------

    def phosphorus_groups(self) -> list[dict[str, Any]]:
        """Every P atom with the count of O it carries and whether it esterifies
        a carbon (directly, `C-P`, or through an oxygen, `C-O-P`)."""
        out = []
        for a in self.atoms:
            if a.element != "P":
                continue
            o_n = 0
            ester_to_c = False
            direct_c = False
            for nb in self.neighbours(a.idx):
                if nb.element in ("O", "N", "S"):
                    o_n += 1
                    if any(x.element == "C" for x in self.neighbours(nb.idx)
                           if x.idx != a.idx):
                        ester_to_c = True
                if nb.element == "C":
                    direct_c = True
            out.append({"o_or_n": o_n, "ester_to_c": ester_to_c,
                        "direct_c_p_bond": direct_c})
        return out

    def n_ester_carbonyls(self) -> int:
        """Count of `C(=O)O-C` ester linkages — the acyl attachment of a lipid."""
        n = 0
        for a in self.atoms:
            if a.element != "C":
                continue
            dbl_o = [nb for nb in self.neighbours(a.idx)
                     if nb.element == "O"
                     and self.bonds.get((min(a.idx, nb.idx), max(a.idx, nb.idx))) == 2.0]
            sng_o = [nb for nb in self.neighbours(a.idx)
                     if nb.element == "O"
                     and self.bonds.get((min(a.idx, nb.idx), max(a.idx, nb.idx))) != 2.0
                     and any(x.element == "C" and x.idx != a.idx
                             for x in self.neighbours(nb.idx))]
            if dbl_o and sng_o:
                n += 1
        return n

    def n_amide_bonds(self) -> int:
        n = 0
        for a in self.atoms:
            if a.element != "C":
                continue
            has_dbl_o = any(
                nb.element == "O"
                and self.bonds.get((min(a.idx, nb.idx), max(a.idx, nb.idx))) == 2.0
                for nb in self.neighbours(a.idx))
            has_n = any(nb.element == "N" for nb in self.neighbours(a.idx))
            if has_dbl_o and has_n:
                n += 1
        return n

    def hydroxyl_count(self) -> int:
        n = 0
        for a in self.atoms:
            if a.element != "O":
                continue
            heavy = [nb for nb in self.neighbours(a.idx) if nb.element != "H"]
            if len(heavy) == 1 and heavy[0].element == "C":
                bo = self.bonds.get((min(a.idx, heavy[0].idx),
                                     max(a.idx, heavy[0].idx)))
                if bo != 2.0:
                    n += 1
        return n

    def ether_oxygens(self) -> int:
        n = 0
        for a in self.atoms:
            if a.element != "O" or a.in_ring:
                continue
            heavy = [nb for nb in self.neighbours(a.idx) if nb.element != "H"]
            if len(heavy) == 2 and all(x.element == "C" for x in heavy):
                if not any(self._is_carbonyl_carbon(x.idx) for x in heavy):
                    n += 1
        return n

    def _is_carbonyl_carbon(self, idx: int) -> bool:
        return any(
            nb.element == "O"
            and self.bonds.get((min(idx, nb.idx), max(idx, nb.idx))) == 2.0
            for nb in self.neighbours(idx))

    def n_aromatic_rings(self) -> int:
        return sum(1 for r in self.rings
                   if len(r) in (5, 6) and all(self.atoms[i].aromatic for i in r))

    def alkyl_sulfonate_groups(self) -> int:
        """S with >= 3 oxygens AND a direct S-C bond — an alkyl/aryl sulfonate.

        The direct C-S bond is load-bearing. Without it the test also matches
        an N-O-SO3 sulfamate, and `9CP` — an avibactam-class beta-lactamase
        inhibitor — was filed as a Good's buffer. Found on a held-out sample,
        not on the tuning set.
        """
        n = 0
        for a in self.atoms:
            if a.element != "S":
                continue
            os_ = [nb for nb in self.neighbours(a.idx) if nb.element == "O"]
            cs = [nb for nb in self.neighbours(a.idx) if nb.element == "C"]
            if len(os_) >= 3 and cs:
                n += 1
        return n

    def electrophilic_halide_carbons(self) -> int:
        """sp3 carbons carrying a Cl, Br or I — alkyl-halide electrophiles.

        Counted because a component carrying SEVERAL of them is a bifunctional
        or trifunctional CROSSLINKING REAGENT, not a drug: a covalent drug
        carries exactly one warhead. `LFI` (the TATA reagent that cyclises
        Bicycle peptides) carries three, `ZBR` (TBMB) three, `A1I4O` three,
        `8VY` (bis(bromomethyl)benzene) two.

        Aryl halides are excluded — a chloro- or fluoro-phenyl is on half the
        drugs ever made and is not an electrophile. `260`
        (2-(bromomethyl)-1,3-difluorobenzene) therefore counts 1, not 3.

        THE COUNT IS ADVISORY AND THE THRESHOLD IS 3, NOT 2, ON PURPOSE.
        Nitrogen mustards — chlorambucil, melphalan, bendamustine — are
        approved drugs carrying exactly TWO alkyl chlorides on one nitrogen. A
        threshold of 2 would file them as reagents, which is the false negative
        this module refuses to trade for a false positive. Three arms is a
        symmetric crosslinker; two is a mustard.
        """
        n = 0
        for a in self.atoms:
            if a.element != "C" or a.aromatic:
                continue
            if any(nb.element in ("CL", "BR", "I") for nb in self.neighbours(a.idx)):
                n += 1
        return n

    # -- scaffold signatures ---------------------------------------------

    def purine_like_rings(self) -> bool:
        """Fused 5+6 all-C/N bicycle carrying >= 3 nitrogens.

        Deliberately NOT keyed on aromatic flags: the CCD writes guanine in 2UK
        as a kekulised `C1=NC(=O)...N1` six-ring fused to an aromatic five-ring,
        so an aromaticity requirement would miss the exact ligand that caused
        the KRAS wrong answer.
        """
        for group in self.fused_groups():
            sizes = sorted(len(r) for r in group)
            if len(group) < 2:
                continue
            for i in range(len(group)):
                for j in range(i + 1, len(group)):
                    r5, r6 = group[i], group[j]
                    if sorted((len(r5), len(r6))) != [5, 6]:
                        continue
                    if len(set(r5) & set(r6)) != 2:
                        continue
                    atoms = set(r5) | set(r6)
                    els = [self.atoms[k].element for k in atoms]
                    if any(e not in ("C", "N") for e in els):
                        continue
                    if els.count("N") >= 3:
                        return True
            del sizes
        return False

    def furanose_like(self) -> bool:
        """5-ring of 4 C + 1 O carrying >= 2 exocyclic oxygens — a ribose."""
        for r in self.rings:
            if len(r) != 5:
                continue
            els = self.ring_elements(r)
            if els.get("O") != 1 or els.get("C") != 4:
                continue
            exo = 0
            for idx in r:
                for nb in self.neighbours(idx):
                    if nb.element == "O" and nb.idx not in r:
                        exo += 1
            if exo >= 2:
                return True
        return False

    def pyranose_like_rings(self) -> int:
        """Sugar rings: 5- or 6-membered, 1 ring O, rest C, >= 2 hydroxyls on
        the ring carbons. A drug's tetrahydropyran or dioxolane does not qualify
        because it does not carry the hydroxyl belt."""
        n = 0
        for r in self.rings:
            if len(r) not in (5, 6):
                continue
            els = self.ring_elements(r)
            if els.get("O") != 1 or els.get("C") != len(r) - 1:
                continue
            if any(self.atoms[i].aromatic for i in r):
                continue
            oh = 0
            for idx in r:
                for nb in self.neighbours(idx):
                    if nb.element != "O" or nb.idx in r:
                        continue
                    heavy = [x for x in self.neighbours(nb.idx) if x.element != "H"]
                    if len(heavy) <= 2:
                        oh += 1
            if oh >= 2:
                n += 1
        return n

    def steroid_nucleus(self) -> bool:
        """Cyclopenta[a]phenanthrene: a fused system containing three fused
        6-rings and one fused 5-ring, every ring atom carbon."""
        for group in self.fused_groups():
            carbo = [r for r in group
                     if all(self.atoms[i].element == "C" for i in r)]
            six = [r for r in carbo if len(r) == 6]
            five = [r for r in carbo if len(r) == 5]
            if len(six) >= 3 and len(five) >= 1:
                # require the five-ring to actually be fused to a six-ring
                if any(len(set(f) & set(s)) >= 2 for f in five for s in six):
                    return True
        return False

    def sterol_side_chain(self) -> int:
        """Longest acyclic CARBON run hanging off a fused carbocyclic system.

        This is what separates cholesterol (isooctyl tail, 8) and the bile
        acids (cholic acid, 5, ending in a carboxylate) from a steroid HORMONE
        or steroid DRUG — testosterone 0, progesterone 2, dexamethasone 2 —
        which have no tail. Without it every steroid drug would be filed as a
        lipid, which is a wrong answer in the opposite direction.

        Carbons on the tail MAY bear oxygens (cholic acid's tail terminates in
        -COOH); what is counted is the length of the carbon run, not its
        purity. Requiring purity truncated cholic acid to 3 and let a bile salt
        through as drug-like.
        """
        ring_atoms = {i for r in self.rings for i in r}
        acyclic_c = {a.idx for a in self.atoms
                     if a.element == "C" and not a.in_ring}
        best = 0
        for start in acyclic_c:
            if not any(nb.idx in ring_atoms for nb in self.neighbours(start)):
                continue
            seen = {start}
            frontier = [start]
            depth = 1
            while frontier:
                nxt = []
                for u in frontier:
                    for v in self.adj[u]:
                        if v in acyclic_c and v not in seen:
                            seen.add(v)
                            nxt.append(v)
                if nxt:
                    depth += 1
                frontier = nxt
            best = max(best, depth)
        return best

    def pyrimidine_base(self) -> bool:
        """A lone (unfused) 6-ring of C and N with >= 2 N carrying an exocyclic
        =O or -NH2 — uracil, thymine, cytosine. Nucleobases are endogenous."""
        for r in self.rings:
            if len(r) != 6:
                continue
            if any(len(set(r) & set(o)) >= 2 for o in self.rings if o is not r):
                continue
            els = [self.atoms[i].element for i in r]
            if any(e not in ("C", "N") for e in els) or els.count("N") < 2:
                continue
            exo = 0
            for idx in r:
                for nb in self.neighbours(idx):
                    if nb.idx not in r and nb.element in ("O", "N"):
                        exo += 1
            if exo >= 1:
                return True
        return False

    def peptide_backbone_residues(self) -> int:
        """Count `N-C-C(=O)` alpha-amino-acid units linked head to tail."""
        n = 0
        for a in self.atoms:
            if a.element != "C" or a.in_ring:
                continue
            if not self._is_carbonyl_carbon(a.idx):
                continue
            for nb in self.neighbours(a.idx):
                if nb.element != "C" or nb.in_ring:
                    continue
                if any(x.element == "N" for x in self.neighbours(nb.idx)
                       if x.idx != a.idx):
                    n += 1
                    break
        return n


# --------------------------------------------------------------------------
# Structural context — the entry the component actually appears in
# --------------------------------------------------------------------------
#
# WHY THIS EXISTS. Everything above answers "what is this molecule?". That is
# the right question for ADP and for cholesteryl hemisuccinate, and it was the
# whole module until 2026-08-15, when a Foldseek/TSLP run measured the first
# false positive against the held-out set's 0-false-positive result:
#
#   8QFZ carries `LFI`, C12H18Br3N3O3 — 1,3,5-tris(3-bromopropanoyl)-1,3,5-
#   triazinane, the TATA tri-electrophile that CYCLISES a Bicycle peptide. It
#   is covalently bonded to all three cysteines of a 12-residue polypeptide
#   (`CHWLENCWRGFC`, entity `8QFZ_2`). It is a macrocyclisation reagent inside
#   a peptide ligand. As a molecule in isolation `druglike` is defensible; the
#   chemistry really is drug-like. In the entry it is not a ligand at all.
#
# The consequence was a MODALITY error, which is what makes it serious rather
# than merely wrong: `druglike` made the entry a holo small-molecule structure,
# the run set `tier: holo` with `tier_note: "drug-like ligand LFI"`, anchored
# the site on it and emitted `ligand_site_jaccard` of 0.769 and 1.000 — the
# strongest site-hypothesis basis this pipeline can produce — for what is
# actually PEPTIDE precedent. Rule 1 of the dossier exists to stop exactly that
# substitution, and this module was making it one layer down.
#
# NOTE WHAT THIS IS *NOT*. It is not "a peptide-binding site is not a site".
# A groove that binds a bicyclic peptide is a demonstrated ligandable surface
# and a perfectly good small-molecule target — MDM2/p53, protease substrate
# grooves, peptide GPCRs. Nothing here rejects a site. It attributes a
# component to the right molecule, so the peptide is reported as peptide
# precedent (see `holo_call`'s `polymer_ligand_precedent`) instead of being
# laundered into small-molecule precedent by the reagent that staples it.
#
# WHERE THE CONTEXT COMES FROM, and the trap in getting it. `_struct_conn` and
# `_struct_ref`. Both are present in the ENTRY mmCIF and in
# `files.rcsb.org/header/<ID>.cif`, and BOTH ARE STRIPPED FROM THE ASSEMBLY
# FILE. Verified on 8QFZ: `8QFZ-assembly1.cif` — the file the pocket-scan
# pipeline actually downloads, because the biological assembly is the right
# coordinate set — contains 23 categories and `_struct_conn` is not one of
# them. A caller that builds context from the coordinates it already has will
# silently get an empty link table and every verdict will fall through
# unchanged. Build it from the header, which `pocket-scan` already fetches per
# entry for `_struct_ref`, so this costs no extra network call.
#
# STILL NO NON-STDLIB IMPORT. The mmCIF reader below is ~90 lines and reads
# five categories. `gemmi` is not imported here and must not be: `pocket-scan`'s
# Modal image has it, the dossier sandbox does not, and this module's whole
# point is that the verdict does not vary with the environment.

_MMCIF_CONTEXT_CATEGORIES = (
    "_entity.",
    "_entity_poly.",
    "_entity_poly_seq.",
    "_struct_asym.",
    "_struct_conn.",
    "_struct_ref.",
)


def _mmcif_tokens(text: str):
    """(value, was_quoted) over an mmCIF document. Stdlib, no dependency.

    Handles the three quoting forms that appear in RCSB files: bare tokens,
    single/double-quoted tokens, and semicolon-delimited multi-line text
    fields. `was_quoted` matters because an unquoted token beginning with `_`
    is a TAG while a quoted one is a value, and `loop_` is a keyword only when
    unquoted.
    """
    lines = text.splitlines()
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        if line.startswith(";"):
            buf = [line[1:]]
            i += 1
            while i < n and not lines[i].startswith(";"):
                buf.append(lines[i])
                i += 1
            i += 1
            yield ("\n".join(buf), True)
            continue
        j, m = 0, len(line)
        while j < m:
            ch = line[j]
            if ch in " \t":
                j += 1
                continue
            if ch == "#":
                break
            if ch in "'\"":
                k = j + 1
                while k < m:
                    if line[k] == ch and (k + 1 >= m or line[k + 1] in " \t"):
                        break
                    k += 1
                yield (line[j + 1:k], True)
                j = k + 1
                continue
            k = j
            while k < m and line[k] not in " \t":
                k += 1
            yield (line[j:k], False)
            j = k
        i += 1


def read_mmcif_categories(
    text: str, categories: Sequence[str] = _MMCIF_CONTEXT_CATEGORIES
) -> dict[str, list[dict[str, str]]]:
    """`{'_struct_conn.': [{'conn_type_id': 'covale', ...}, ...]}`.

    Single-item categories come back as a one-row list, so a consumer does not
    have to care whether the depositor used a `loop_`. Values of `.` and `?`
    (mmCIF's not-applicable and unknown) become `''`.
    """
    wanted = tuple(categories)
    out: dict[str, list[dict[str, str]]] = {}

    def norm(v: str, quoted: bool) -> str:
        return "" if (not quoted and v in (".", "?")) else v

    toks = list(_mmcif_tokens(text))
    i, n = 0, len(toks)
    while i < n:
        val, quoted = toks[i]
        if not quoted and val == "loop_":
            i += 1
            tags: list[str] = []
            while i < n and not toks[i][1] and toks[i][0].startswith("_"):
                tags.append(toks[i][0])
                i += 1
            cat = tags[0].rsplit(".", 1)[0] + "." if tags else ""
            keys = [t.split(".", 1)[1] for t in tags]
            rows: list[list[str]] = []
            cur: list[str] = []
            while i < n and not (
                not toks[i][1]
                and (toks[i][0].startswith("_") or toks[i][0] in ("loop_",)
                     or toks[i][0].lower().startswith("data_"))
            ):
                cur.append(norm(*toks[i]))
                i += 1
                if len(cur) == len(keys):
                    rows.append(cur)
                    cur = []
            if cat in wanted:
                out.setdefault(cat, []).extend(dict(zip(keys, r)) for r in rows)
            continue
        if not quoted and val.startswith("_") and "." in val:
            cat = val.rsplit(".", 1)[0] + "."
            key = val.split(".", 1)[1]
            if i + 1 < n:
                if cat in wanted:
                    rows = out.setdefault(cat, [])
                    if not rows:
                        rows.append({})
                    rows[0][key] = norm(*toks[i + 1])
                i += 2
                continue
        i += 1
    return out


@dataclass(frozen=True)
class PolymerEntity:
    """One polymer entity of an entry — a protein chain, a peptide, a nucleic
    acid. `accessions` is what `_struct_ref` declares for it."""
    entity_id: str
    poly_type: str | None          # `_entity_poly.type`, e.g. 'polypeptide(L)'
    description: str | None        # `_entity.pdbx_description`
    n_monomers: int | None
    sequence: str | None           # one-letter, canonical
    accessions: tuple[str, ...] = ()
    strand_ids: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["accessions"] = list(self.accessions)
        d["strand_ids"] = list(self.strand_ids)
        return d


@dataclass(frozen=True)
class CovalentLink:
    """One `_struct_conn` row of type `covale` touching a component."""
    partner_comp_id: str
    partner_entity_id: str | None
    partner_is_polymer: bool
    partner_chain: str | None
    partner_seq_id: str | None
    partner_is_target: bool | None      # None == target identity unknown


@dataclass(frozen=True)
class StructureContext:
    """The facts about ONE PDB entry that change a component's attribution.

    Built from the entry or header mmCIF; never from the assembly file, which
    does not carry `_struct_conn`. `target_entity_ids` is what makes the
    covalent-inhibitor control work — see `_apply_structure_context`.
    """
    entry_id: str | None = None
    polymer_entities: Mapping[str, PolymerEntity] = field(default_factory=dict)
    entity_of_comp: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    polymer_monomers: Mapping[str, tuple[str, ...]] = field(default_factory=dict)
    covalent_links: Mapping[str, tuple[CovalentLink, ...]] = field(default_factory=dict)
    target_entity_ids: tuple[str, ...] = ()
    target_accession: str | None = None
    target_basis: str | None = None
    source: str = "mmcif"
    #: Was `_struct_conn` PRESENT AS A CATEGORY in the file this was built
    #: from? This is not the same as "there were covalent links", and the
    #: difference is the single easiest way to reintroduce the bug: an
    #: ASSEMBLY file has no `_struct_conn` at all, so a context built from one
    #: has an empty link table that looks exactly like "nothing is bonded to
    #: anything". Verified on 8QFZ — a context built from
    #: `8QFZ-assembly1.cif` puts `LFI` back to `druglike`. Every verdict
    #: reached with this False carries the flag
    #: `struct_conn_absent_from_context`.
    has_struct_conn_category: bool = False

    # -- construction ----------------------------------------------------

    @classmethod
    def from_mmcif_text(cls, text: str, *, entry_id: str | None = None,
                        target_accession: str | None = None,
                        target_entity_ids: Iterable[str] | None = None
                        ) -> "StructureContext":
        cats = read_mmcif_categories(text)
        ent_rows = cats.get("_entity.", [])
        poly_rows = cats.get("_entity_poly.", [])
        seq_rows = cats.get("_entity_poly_seq.", [])
        asym_rows = cats.get("_struct_asym.", [])
        conn_rows = cats.get("_struct_conn.", [])
        ref_rows = cats.get("_struct_ref.", [])

        ent_type = {r.get("id", ""): (r.get("type") or "").lower() for r in ent_rows}
        ent_desc = {r.get("id", ""): r.get("pdbx_description") or None
                    for r in ent_rows}
        poly_by_ent = {r.get("entity_id", ""): r for r in poly_rows}
        asym_to_ent = {r.get("id", ""): r.get("entity_id", "") for r in asym_rows}

        accs: dict[str, list[str]] = {}
        for r in ref_rows:
            if (r.get("db_name") or "").upper() not in ("UNP", "UNIPROT"):
                continue
            a = (r.get("pdbx_db_accession") or r.get("db_code") or "").strip()
            if a:
                accs.setdefault(r.get("entity_id", ""), []).append(a)

        monomers: dict[str, list[str]] = {}
        n_mon: dict[str, int] = {}
        for r in seq_rows:
            e = r.get("entity_id", "")
            monomers.setdefault(e, [])
            mid = (r.get("mon_id") or "").upper()
            if mid and mid not in monomers[e]:
                monomers[e].append(mid)
            n_mon[e] = n_mon.get(e, 0) + 1

        polymers: dict[str, PolymerEntity] = {}
        for eid, ty in ent_type.items():
            if ty != "polymer":
                continue
            p = poly_by_ent.get(eid, {})
            seq = (p.get("pdbx_seq_one_letter_code_can")
                   or p.get("pdbx_seq_one_letter_code") or "")
            seq = "".join(seq.split()) or None
            strands = tuple(
                s for s in (p.get("pdbx_strand_id") or "").replace(",", " ").split()
            )
            polymers[eid] = PolymerEntity(
                entity_id=eid,
                poly_type=(p.get("type") or None),
                description=ent_desc.get(eid),
                n_monomers=n_mon.get(eid) or (len(seq) if seq else None),
                sequence=seq,
                accessions=tuple(accs.get(eid, ())),
                strand_ids=strands,
            )

        # comp_id -> entity ids, recovered from the linkage table's own
        # asym ids. Enough to say which entity a bonded component belongs to.
        comp_ent: dict[str, list[str]] = {}
        for r in conn_rows:
            for side in ("1", "2"):
                c = (r.get(f"ptnr{side}_label_comp_id") or "").upper()
                a = r.get(f"ptnr{side}_label_asym_id") or ""
                e = asym_to_ent.get(a)
                if c and e and e not in comp_ent.setdefault(c, []):
                    comp_ent[c].append(e)

        # -- the target entity, which decides bonded-to-the-target vs
        # bonded-to-another-polymer. Explicit ids win; then an accession match;
        # then nothing, and "nothing" is a real state that the rules below
        # handle conservatively rather than guessing at.
        tgt: tuple[str, ...] = tuple(target_entity_ids or ())
        basis = "caller: target_entity_ids" if tgt else None
        if not tgt and target_accession:
            want = target_accession.strip().upper()
            tgt = tuple(e for e, p in polymers.items()
                        if want in {a.upper() for a in p.accessions})
            basis = (f"_struct_ref accession {target_accession}"
                     if tgt else None)

        links: dict[str, list[CovalentLink]] = {}
        for r in conn_rows:
            if (r.get("conn_type_id") or "").lower() != "covale":
                continue
            for side, other in (("1", "2"), ("2", "1")):
                c = (r.get(f"ptnr{side}_label_comp_id") or "").upper()
                oc = (r.get(f"ptnr{other}_label_comp_id") or "").upper()
                oa = r.get(f"ptnr{other}_label_asym_id") or ""
                oe = asym_to_ent.get(oa)
                if not c or not oc:
                    continue
                links.setdefault(c, []).append(CovalentLink(
                    partner_comp_id=oc,
                    partner_entity_id=oe,
                    partner_is_polymer=bool(oe and oe in polymers),
                    partner_chain=(r.get(f"ptnr{other}_auth_asym_id") or oa or None),
                    partner_seq_id=(r.get(f"ptnr{other}_auth_seq_id") or None),
                    partner_is_target=(None if not tgt else (oe in tgt)),
                ))

        return cls(
            entry_id=entry_id or None,
            polymer_entities=polymers,
            entity_of_comp={k: tuple(v) for k, v in comp_ent.items()},
            polymer_monomers={k: tuple(v) for k, v in monomers.items()},
            covalent_links={k: tuple(v) for k, v in links.items()},
            target_entity_ids=tgt,
            target_accession=target_accession,
            target_basis=basis,
            has_struct_conn_category=("_struct_conn." in cats),
        )

    @classmethod
    def from_mmcif_path(cls, path: str | os.PathLike[str], **kw) -> "StructureContext":
        p = Path(path)
        return cls.from_mmcif_text(
            p.read_text(errors="replace"),
            entry_id=kw.pop("entry_id", None) or p.stem.split("-")[0].split("_")[0].upper(),
            **kw,
        )

    # -- queries ---------------------------------------------------------

    def links_for(self, comp_id: str) -> tuple[CovalentLink, ...]:
        return tuple(self.covalent_links.get(comp_id.upper(), ()))

    def is_polymer_monomer(self, comp_id: str) -> str | None:
        """Entity id of the polymer this comp_id is a MONOMER of, if any.

        `_entity_poly_seq` is the authority: a component listed there is a
        residue of that chain, not a bound ligand, whatever the CCD type says.
        Catches `NH2` and `ACE` capping groups, which are two of the ground-
        truth set's three standing misses.
        """
        c = comp_id.upper()
        for eid, mons in self.polymer_monomers.items():
            if c in mons:
                return eid
        return None

    def is_available(self) -> bool:
        """Can the COVALENT rules run at all?

        False when `_struct_conn` was not in the file, which is exactly what an
        assembly file looks like. Deliberately NOT "there are polymer entities"
        — that was the first version of this method and it returned True on
        `8QFZ-assembly1.cif`, where the category is stripped, so `LFI` went
        straight back to `druglike` with no sign anything had gone wrong.
        """
        return self.has_struct_conn_category

    def to_dict(self) -> dict[str, Any]:
        return {
            "entry_id": self.entry_id,
            "n_polymer_entities": len(self.polymer_entities),
            "target_entity_ids": list(self.target_entity_ids),
            "target_accession": self.target_accession,
            "target_basis": self.target_basis,
            "has_struct_conn_category": self.has_struct_conn_category,
            "n_components_with_covalent_links": len(self.covalent_links),
            "source": self.source,
        }


# --------------------------------------------------------------------------
# Chemical component source — Paperclip `pdb_v.chemcomps`
# --------------------------------------------------------------------------

_DEFAULT_ENV_FILE = "/Users/bb/repos/claude-agent-starter/.env"


def _load_env(env_file: str | os.PathLike[str] | None) -> dict[str, str]:
    env = dict(os.environ)
    path = Path(env_file) if env_file else Path(_DEFAULT_ENV_FILE)
    if path.is_file():
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            env.setdefault(k.strip(), v.strip().strip('"').strip("'"))
    return env


def _parse_paperclip_table(out: str) -> list[dict[str, str]]:
    """Parse paperclip's ASCII table. Column spans come from the `---+---` rule,
    which is the only reliable delimiter — cell values contain `|` (SMILES do
    not, but names do)."""
    lines = [l for l in out.splitlines() if l.strip()]
    hdr = None
    for i, l in enumerate(lines[:-1]):
        if set(lines[i + 1].strip()) <= set("-+ ") and "-" in lines[i + 1]:
            hdr = i
            break
    if hdr is None:
        return []
    sep = lines[hdr + 1]
    spans, start = [], 0
    for j, ch in enumerate(sep):
        if ch == "+":
            spans.append((start, j))
            start = j + 1
    spans.append((start, len(sep) + 4096))
    cols = [lines[hdr][a:b].strip() for a, b in spans]
    rows = []
    for l in lines[hdr + 2:]:
        if l.startswith("(") or l.startswith("["):
            break
        rows.append({c: l[a:b].strip() for c, (a, b) in zip(cols, spans)})
    return rows


class ChemCompSource:
    """Fetches CCD rows from Paperclip `pdb_v.chemcomps`, with a process cache.

    Paperclip is the only source. It carries `type`, `name`, `formula`,
    `formula_weight`, `smiles`, `inchikey` and `drugbank_id` — every field the
    classification needs — so nothing is fetched from RCSB or anywhere else.

    TWO PAPERCLIP CONSTRAINTS ARE HANDLED HERE, both measured:
      * wide cells are truncated in the rendered table, so `smiles` and `name`
        are pulled as fixed-width SUBSTRING slices and rejoined. `B12`'s SMILES
        is 197 characters and comes back whole this way.
      * 200-row cap and a statement timeout, so comp_ids are batched at 40 and
        the list is inlined rather than subqueried.
    """

    def __init__(self, *, env_file: str | os.PathLike[str] | None = None,
                 paperclip: str = "paperclip", cache_path: str | os.PathLike[str] | None = None,
                 timeout: float = 60.0):
        self._env = _load_env(env_file)
        self._paperclip = paperclip
        self._timeout = timeout
        self._cache: dict[str, dict[str, Any] | None] = {}
        self._cache_path = Path(cache_path) if cache_path else None
        if self._cache_path and self._cache_path.is_file():
            try:
                self._cache.update(json.loads(self._cache_path.read_text()))
            except Exception:                          # noqa: BLE001
                pass
        self.last_error: str | None = None
        #: comp_id -> why the lookup failed. A LOOKUP FAILURE IS NOT A MISS.
        #: Paperclip's public endpoint intermittently exceeds its statement
        #: timeout; without this the failure is indistinguishable from "the CCD
        #: has no such component", and both would silently render as apo. That
        #: is the same fail-open shape as the bugs this module replaces.
        self.fetch_errors: dict[str, str] = {}

    def preload(self, records: Mapping[str, Mapping[str, Any]]) -> None:
        """Seed the cache — used by tests and by call sites that already hold
        `pdb_v.entry_ligands` rows and should not re-query."""
        for k, v in records.items():
            self._cache[k.upper()] = dict(v)

    def get(self, comp_id: str) -> dict[str, Any] | None:
        return self.get_many([comp_id]).get(comp_id.upper())

    def get_many(self, comp_ids: Iterable[str]) -> dict[str, dict[str, Any] | None]:
        want = [c.upper() for c in comp_ids if c]
        todo = sorted({c for c in want if c not in self._cache})
        for i in range(0, len(todo), 40):
            self._fetch_batch(todo[i:i + 40])
        return {c: self._cache.get(c) for c in want}

    def _fetch_batch(self, batch: list[str], *, attempts: int = 3) -> None:
        inlist = ", ".join("'" + c.replace("'", "''") + "'" for c in batch)
        q = (
            "SELECT comp_id, type, formula, formula_weight, drugbank_id, inchikey, "
            "SUBSTRING(smiles,1,70) s0, SUBSTRING(smiles,71,70) s1, "
            "SUBSTRING(smiles,141,70) s2, SUBSTRING(smiles,211,70) s3, "
            "SUBSTRING(smiles,281,70) s4, SUBSTRING(name,1,90) nm "
            f"FROM pdb_v.chemcomps WHERE comp_id IN ({inlist})"
        )
        # Retried because the endpoint is measurably flaky: identical queries
        # that return in 30 ms also intermittently exceed the statement timeout.
        # An unretried timeout costs a real holo structure.
        err = "no attempt made"
        for _ in range(max(1, attempts)):
            try:
                proc = subprocess.run(
                    [self._paperclip, "sql", "-s", "proteins", q],
                    capture_output=True, text=True, env=self._env,
                    timeout=self._timeout, stdin=subprocess.DEVNULL,
                )
            except Exception as exc:                   # noqa: BLE001
                err = f"{type(exc).__name__}: {exc}"
                continue
            if proc.returncode != 0:
                err = (proc.stderr or proc.stdout)[:400]
                continue
            rows = _parse_paperclip_table(proc.stdout)
            if not rows and not _looks_like_empty_result(proc.stdout):
                err = f"unparseable paperclip output: {proc.stdout[:200]!r}"
                continue
            break
        else:
            # Every attempt failed. Mark the WHOLE batch as errored; do NOT
            # leave it looking like a clean CCD miss.
            self.last_error = err
            for c in batch:
                self.fetch_errors[c] = err
            return

        for c in batch:
            self._cache.setdefault(c, None)
        for r in rows:
            cid = (r.get("comp_id") or "").strip()
            if not cid:
                continue
            smi = "".join(r.get(f"s{i}") or "" for i in range(5)).strip()
            self._cache[cid.upper()] = {
                "comp_id": cid,
                "type": _nn(r.get("type")),
                "formula": _nn(r.get("formula")),
                "formula_weight": _nf(r.get("formula_weight")),
                "drugbank_id": _nn(r.get("drugbank_id")),
                "inchikey": _nn(r.get("inchikey")),
                "smiles": _nn(smi),
                "name": _nn(r.get("nm")),
            }

    def with_entry_counts(self, comp_ids: Iterable[str]) -> dict[str, int]:
        """`comp_id -> number of PDB entries carrying it`, folded into the
        cached records as `n_pdb_entries` so R14 can read it.

        BEST EFFORT AND NEVER FATAL. A failure leaves `n_pdb_entries` absent,
        which R14 reads as NOT CHECKED rather than as rare — the same
        distinction `fetch_errors` makes for the CCD row itself, and for the
        same reason. Measured: this aggregate returns in ~0.6 s for 15 ids;
        the analogous join against `uniprot_v.pdb_chains`, which would give the
        better statistic (distinct proteins rather than entries), times out at
        120 s and is not usable.

        TWO SOURCES, BECAUSE ONE OF THEM FAILED UNDER LOAD ON THE DAY THIS WAS
        WRITTEN. Paperclip's aggregate is tried first; with a dozen agents
        against the same endpoint every batch of 35 came back empty over ten
        minutes while single-row selects still answered in milliseconds. RCSB's
        search API answers the identical question by `return_counts` and is
        stdlib-reachable. The two agree where both answered: BEN 361/361,
        LZ1 10/10, ZBR 9/9, LFI 8/8, MOV 7/7, GOL 26004/26117 (paperclip is a
        snapshot, RCSB is current). A prior that silently degrades to "not
        checked" whenever the cluster is busy is a prior that does nothing on
        exactly the runs that matter.
        """
        want = [c.upper() for c in comp_ids if c]
        out: dict[str, int] = {}
        for i in range(0, len(want), 40):
            batch = want[i:i + 40]
            inlist = ", ".join("'" + c.replace("'", "''") + "'" for c in batch)
            q = ("SELECT comp_id, count(distinct entry_id) n FROM "
                 f"pdb_v.entry_ligands WHERE comp_id IN ({inlist}) GROUP BY 1")
            try:
                proc = subprocess.run(
                    [self._paperclip, "sql", "-s", "proteins", q],
                    capture_output=True, text=True, env=self._env,
                    timeout=self._timeout, stdin=subprocess.DEVNULL,
                )
            except Exception as exc:                   # noqa: BLE001
                self.last_error = f"entry counts: {type(exc).__name__}: {exc}"
                continue
            if proc.returncode != 0:
                self.last_error = f"entry counts: {(proc.stderr or proc.stdout)[:200]}"
                continue
            for r in _parse_paperclip_table(proc.stdout):
                cid = (r.get("comp_id") or "").strip().upper()
                n = _ni(r.get("n"))
                if cid and n is not None:
                    out[cid] = n
        for cid in want:
            if cid not in out:
                n = _rcsb_entry_count(cid)
                if n is not None:
                    out[cid] = n
        for cid, n in out.items():
            rec = self._cache.get(cid)
            if isinstance(rec, dict):
                rec["n_pdb_entries"] = n
        return out

    def save_cache(self) -> None:
        if self._cache_path:
            self._cache_path.parent.mkdir(parents=True, exist_ok=True)
            self._cache_path.write_text(json.dumps(self._cache, sort_keys=True))


def _rcsb_entry_count(comp_id: str, *, timeout: float = 30.0) -> int | None:
    """How many PDB entries carry this component, from RCSB's search API.

    The fallback for `ChemCompSource.with_entry_counts`. Stdlib `urllib` only —
    no `requests`, no toolkit. The searchable attribute is
    `rcsb_chem_comp_container_identifiers.comp_id`; `chem_comp.id` is NOT
    indexed for search and returns HTTP 400 with
    "search is not enabled on [ chem_comp.id ] attribute".

    Returns None on any failure, which R14 reads as NOT CHECKED.
    """
    import urllib.parse
    import urllib.request

    q = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
            "operator": "exact_match", "value": comp_id.upper()}},
        "return_type": "entry",
        "request_options": {"return_counts": True},
    }
    url = ("https://search.rcsb.org/rcsbsearch/v2/query?"
           + urllib.parse.urlencode({"json": json.dumps(q)}))
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:  # noqa: S310
            return _ni(json.loads(r.read().decode()).get("total_count"))
    except Exception:                                  # noqa: BLE001
        return None


def _looks_like_empty_result(out: str) -> bool:
    """A genuine zero-row answer, as opposed to an error or a hang."""
    return "(0 rows" in out


def _nn(v: Any) -> Any:
    return None if v in (None, "", "NULL") else v


def _nf(v: Any) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _ni(v: Any) -> int | None:
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


_DEFAULT_SOURCE: ChemCompSource | None = None


def _default_source() -> ChemCompSource:
    global _DEFAULT_SOURCE
    if _DEFAULT_SOURCE is None:
        _DEFAULT_SOURCE = ChemCompSource()
    return _DEFAULT_SOURCE


def set_default_source(src: ChemCompSource) -> None:
    """Swap the module-level source. Tests and offline callers use this."""
    global _DEFAULT_SOURCE
    _DEFAULT_SOURCE = src


# --------------------------------------------------------------------------
# Verdict
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class LigandVerdict:
    comp_id: str
    name: str | None
    formula: str | None
    heavy_atoms: int | None
    mw: float | None
    verdict: str
    reason: str
    evidence: dict[str, Any] = field(default_factory=dict)
    flags: tuple[str, ...] = ()
    comp_type: str | None = None
    smiles: str | None = None
    drugbank_id: str | None = None
    confidence: str = "high"
    source: str = "paperclip:pdb_v.chemcomps"

    @property
    def is_druglike(self) -> bool:
        return self.verdict == "druglike"

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["flags"] = list(self.flags)
        d["is_druglike"] = self.is_druglike
        return d


def _v(comp_id, rec, verdict, reason, evidence, *, flags=(), confidence="high",
       source="paperclip:pdb_v.chemcomps") -> LigandVerdict:
    rec = rec or {}
    extra = list(flags)
    if comp_id.upper() in _FREQUENT_HITTERS:
        extra.append("promiscuity_advisory")
        evidence = dict(evidence)
        evidence["promiscuity_advisory"] = _FREQUENT_HITTERS[comp_id.upper()]
    return LigandVerdict(
        comp_id=comp_id.upper(),
        name=rec.get("name"),
        formula=rec.get("formula"),
        heavy_atoms=evidence.get("heavy_atoms"),
        mw=rec.get("formula_weight"),
        verdict=verdict,
        reason=reason,
        evidence=evidence,
        flags=tuple(extra),
        comp_type=rec.get("type"),
        smiles=rec.get("smiles"),
        drugbank_id=rec.get("drugbank_id"),
        confidence=confidence,
        source=source,
    )


# --------------------------------------------------------------------------
# The rules
# --------------------------------------------------------------------------


def classify_record(rec: Mapping[str, Any] | None, comp_id: str | None = None
                    ) -> LigandVerdict:
    """Classify from a CCD record already in hand. No I/O.

    `rec` needs `type`, `formula`, `formula_weight`, `smiles`, `name` — the
    exact column set of `pdb_v.chemcomps` and (bar `type` vs `comp_type`) of
    `pdb_v.entry_ligands`, so a call site that already selected either can
    classify without a second round trip.
    """
    cid = (comp_id or (rec or {}).get("comp_id") or "?").upper()
    if not rec:
        return _v(cid, None, "unknown",
                  "no record for this comp_id in pdb_v.chemcomps; refusing to "
                  "guess (an unclassified ligand is not evidence of a site)",
                  {"heavy_atoms": None, "ccd_hit": False},
                  confidence="low", source="paperclip:pdb_v.chemcomps (miss)")

    # Accept a `pdb_v.entry_ligands` row verbatim as well as a
    # `pdb_v.chemcomps` one: the two views spell the same two fields
    # differently, and a call site should not have to rename them.
    rec = dict(rec)
    if not rec.get("type") and rec.get("comp_type"):
        rec["type"] = rec["comp_type"]
    if not rec.get("name") and rec.get("ligand_name"):
        rec["name"] = rec["ligand_name"]

    ctype = (rec.get("type") or "").strip().lower()
    formula = rec.get("formula")
    mw = _nf(rec.get("formula_weight"))
    smiles = rec.get("smiles")
    els = parse_formula(formula)
    heavy = heavy_atom_count(formula)

    ev: dict[str, Any] = {
        "ccd_hit": True,
        "chem_comp_type": rec.get("type"),
        "formula": formula,
        "formula_weight": mw,
        "heavy_atoms": heavy,
        "elements": els or None,
        "has_nitrogen": bool(els.get("N")),
        "has_phosphorus": bool(els.get("P")),
        "metals": sorted(e for e in els if e in _METALS) or None,
        "smiles_present": bool(smiles),
        # Optional, supplied by the caller or by
        # `ChemCompSource.with_entry_counts()`. `None` means NOT CHECKED, which
        # is not the same as "rare" and never reads as one — see R14.
        "n_pdb_entries": _ni(rec.get("n_pdb_entries")),
    }

    # ---- R0. The CCD's own placeholder components. `UNK`, `UNL`, `UNX` and
    # friends are density that was never identified. Keyed on the CCD `name`
    # field, not on a hardcoded id list, so a new placeholder code is caught.
    name_l = (rec.get("name") or "").strip().lower()
    if name_l.startswith("unknown") or name_l.startswith("unidentified"):
        return _v(cid, rec, "unknown",
                  f"the CCD names this component {rec.get('name')!r}: "
                  "unidentified density, not a characterised ligand", ev,
                  confidence="high")

    # ---- R1. `_chem_comp.type`: polymer residues are never small-molecule
    # evidence. This is the field the old code ignored and it is decisive on
    # its own. 6OIM's GDP is typed `RNA linking`.
    if ctype:
        is_sacch = ("saccharide" in ctype)
        if is_sacch:
            # A saccharide carrying a phosphate is a sugar phosphate / metabolite;
            # a saccharide carrying an alkyl tail is an alkyl-glycoside detergent
            # (BOG, HTG, LMT). Neither is "a sugar" in the sense a call site
            # cares about, and both are still not drug-like.
            if parse_formula(formula).get("P"):
                return _v(cid, rec, "cofactor",
                          f"CCD type {rec['type']!r} plus phosphorus: a sugar "
                          "phosphate / nucleotide-sugar metabolite", ev)
            gg = SmilesGraph(smiles) if smiles else None
            if gg is not None and gg.ok:
                ch = gg.longest_aliphatic_chain()
                ev["longest_aliphatic_chain"] = ch
                if ch >= 6:
                    return _v(cid, rec, "lipid_or_detergent",
                              f"CCD type {rec['type']!r} with an unbranched alkyl "
                              f"chain of {ch} carbons: an alkyl-glycoside "
                              "detergent (octyl glucoside / dodecyl maltoside "
                              "class), not a glycan", ev)
            return _v(cid, rec, "sugar_or_glycan",
                      f"CCD _chem_comp.type is {rec['type']!r}", ev)
        if "linking" in ctype or "terminus" in ctype:
            if "peptide" in ctype:
                return _v(cid, rec, "peptide_or_polymer",
                          f"CCD _chem_comp.type is {rec['type']!r}: a peptide "
                          "polymer residue (modified/standard amino acid)", ev)
            if "dna" in ctype or "rna" in ctype:
                return _v(cid, rec, "cofactor",
                          f"CCD _chem_comp.type is {rec['type']!r}: a nucleotide "
                          "polymer residue — a nucleotide cofactor when it "
                          "appears free (this is how GDP in 6OIM is typed)", ev)
        if ctype == "peptide-like":
            return _v(cid, rec, "peptide_or_polymer",
                      f"CCD _chem_comp.type is {rec['type']!r}: a peptide-like "
                      "component. It may be a genuine binder, but it is not "
                      "small-molecule evidence", ev,
                      flags=("standalone_peptide_ligand",))

    # ---- R2. Elemental / size floor. No SMILES needed.
    n_c = els.get("C", 0)
    metals = [e for e in els if e in _METALS]
    n_metal_atoms = sum(els[e] for e in metals)

    if heavy is not None:
        # Metal clusters (Fe2S2, Fe4S4, Mo-pterin) are cofactors, not ions.
        if n_metal_atoms >= 2 and n_c == 0:
            return _v(cid, rec, "cofactor",
                      f"inorganic metal cluster ({''.join(sorted(metals))}, "
                      f"{n_metal_atoms} metal atoms, no carbon): an "
                      "iron-sulfur / metal-cluster cofactor", ev)
        if heavy <= 9 and n_c == 0:
            return _v(cid, rec, "ion_or_solvent",
                      f"no carbon and {heavy} heavy atoms: a bare ion or simple "
                      "inorganic species (sulfate, phosphate, pyrophosphate, "
                      "halide, metal)", ev)
        # One carbon and essentially no hydrogen is an inorganic oxyanion
        # (carbonate, bicarbonate, thiocyanate). The hydrogen bound keeps
        # formic acid and urea out — they are bench additives, not ions.
        if heavy <= 5 and n_c <= 1 and els.get("H", 0) <= 1:
            return _v(cid, rec, "ion_or_solvent",
                      f"{heavy} heavy atoms, {n_c} carbon and "
                      f"{els.get('H', 0)} hydrogen: a simple inorganic "
                      "oxyanion (carbonate / thiocyanate class)", ev)
        if heavy <= 2:
            return _v(cid, rec, "ion_or_solvent",
                      f"{heavy} heavy atom(s)", ev)

    # ---- R3. Metal-containing organic cofactors: heme, chlorophyll, B12.
    if metals and heavy and heavy >= 10:
        return _v(cid, rec, "cofactor",
                  f"organometallic component carrying {'/'.join(sorted(metals))} "
                  f"over {heavy} heavy atoms: a metalloporphyrin/corrin-class "
                  "cofactor, not a synthetic ligand", ev)

    # ---- Beyond here the topology matters.
    g = SmilesGraph(smiles) if smiles else None
    if g is not None and not g.ok:
        ev["smiles_parse_error"] = getattr(g, "parse_error", "unknown")
        g = None

    if g is None:
        return _classify_without_smiles(cid, rec, ev, heavy, els, ctype)

    ev.update({
        "n_rings": len(g.rings),
        "n_aromatic_rings": g.n_aromatic_rings(),
        "longest_aliphatic_chain": g.longest_aliphatic_chain(),
        "phosphorus_groups": g.phosphorus_groups(),
        "n_ester_carbonyls": g.n_ester_carbonyls(),
        "n_amide_bonds": g.n_amide_bonds(),
        "hydroxyls": g.hydroxyl_count(),
        "ether_oxygens": g.ether_oxygens(),
        "purine_like": g.purine_like_rings(),
        "furanose_like": g.furanose_like(),
        "sugar_rings": g.pyranose_like_rings(),
        "steroid_nucleus": g.steroid_nucleus(),
        "sterol_side_chain": g.sterol_side_chain(),
        "peptide_residues": g.peptide_backbone_residues(),
        "alkyl_sulfonates": g.alkyl_sulfonate_groups(),
        "pyrimidine_base": g.pyrimidine_base(),
        "electrophilic_halide_carbons": g.electrophilic_halide_carbons(),
    })

    chain = ev["longest_aliphatic_chain"]
    pgroups = ev["phosphorus_groups"]
    phospho_ester = any(p["ester_to_c"] and p["o_or_n"] >= 2 for p in pgroups)
    aromatic = ev["n_aromatic_rings"]

    # ---- R4. Free nucleobases. Endogenous, and checked BEFORE the size floor
    # because uracil is 8 heavy atoms and would otherwise be swallowed by it.
    # Note the purine test does not require aromatic flags — the CCD kekulises
    # several of these, uracil among them.
    if ev["purine_like"] and heavy is not None and heavy <= 12:
        return _v(cid, rec, "cofactor",
                  f"bare purine base ({heavy} heavy atoms, no ribose or "
                  "phosphate): adenine/guanine-class nucleobase", ev)
    if ev["pyrimidine_base"] and heavy is not None and heavy <= 10:
        return _v(cid, rec, "cofactor",
                  f"bare pyrimidine base ({heavy} heavy atoms): "
                  "uracil/thymine/cytosine-class nucleobase", ev)

    # ---- R5. Below the reporting floor. Two exemptions, both measured:
    #   * an aromatic ring — `LZ1` (1H-indazole) is 9 heavy atoms and IS a real
    #     fragment hit, whereas `DEP` (diethyl phosphonate, 8, no ring) is bench
    #     chemistry. Without it the floor eats fragment-screen precedent, which
    #     is the same class of error the old 18-heavy-atom floor made.
    #   * a long alkyl chain — `OCT` is n-octane, 8 heavy atoms and unambiguously
    #     greasy; it belongs with the lipids, not the buffers.
    if heavy is not None and heavy <= TRIVIAL_MAX_HEAVY_ATOMS \
            and not (aromatic >= 1 and heavy >= 8) \
            and chain < LIPID_MIN_CHAIN:
        return _v(cid, rec, "crystallisation_additive",
                  f"only {heavy} heavy atoms, no aromatic ring and no alkyl "
                  "chain: below anything a fragment screen reports; a buffer "
                  "component, cryoprotectant or solvent", ev)

    # ---- R6. Sterol. Cholesterol and its hemisuccinate — the CD20 bug.
    # The side chain is what separates a membrane sterol from a steroid DRUG.
    if ev["steroid_nucleus"]:
        if ev["sterol_side_chain"] >= 4:
            return _v(cid, rec, "lipid_or_detergent",
                      "cyclopenta[a]phenanthrene (steroid) nucleus carrying an "
                      f"aliphatic side chain of {ev['sterol_side_chain']} carbons: "
                      "a membrane sterol or bile salt, not a drug. This is the "
                      "cholesterol / cholesteryl-hemisuccinate class that "
                      "produced the CD20 wrong answer", ev)
        return _v(cid, rec, "druglike",
                  "steroid nucleus with no aliphatic side chain "
                  f"({ev['sterol_side_chain']} carbons): a steroid hormone or "
                  "steroid drug, which IS evidence of a bindable site. Reported "
                  "separately from sterols on purpose — filing every steroid as "
                  "a lipid is the same failure in the other direction", ev,
                  flags=("steroid_nucleus",))

    # ---- R7. Nucleotide / nucleoside. The 2UK bug, and ADP/ATP on NLRP3.
    if ev["purine_like"] and (ev["furanose_like"] or phospho_ester):
        return _v(cid, rec, "cofactor",
                  "purine (adenine/guanine-class) base fused bicycle plus "
                  + ("a ribose/deoxyribose furanose" if ev["furanose_like"]
                     else "a phosphate ester")
                  + ": a nucleotide/nucleoside cofactor. Catches ADP, ATP, GDP, "
                    "GTP, GNP and non-hydrolysable analogs such as 2UK "
                    "(GppNHp) without naming any of them", ev)
    if ev["furanose_like"] and phospho_ester:
        return _v(cid, rec, "cofactor",
                  "ribofuranose bearing a phosphate ester: a nucleotide-class "
                  "cofactor or sugar phosphate", ev)

    # ---- R6. Phospholipid: phospho head group + acyl chains. CD20's PC.
    if phospho_ester and (chain >= LIPID_MIN_CHAIN or ev["n_ester_carbonyls"] >= 2):
        return _v(cid, rec, "lipid_or_detergent",
                  f"phosphate ester head group with an aliphatic chain of "
                  f"{chain} carbons and {ev['n_ester_carbonyls']} acyl esters: a "
                  "phospholipid/detergent (phosphatidylcholine class)", ev)

    # ---- R7. Any other phosphate ester on a small molecule. Free phosphates
    # are the signature of endogenous metabolites and cofactors (FMN, FAD, PLP,
    # TPP, CoA, sugar phosphates), not of drugs.
    if phospho_ester:
        return _v(cid, rec, "cofactor",
                  "carries a phosphate/phosphonate ester. Free phosphates are "
                  "the signature of endogenous cofactors and metabolites; drugs "
                  "essentially never carry one (bisphosphonate drugs and "
                  "phosphate prodrugs are the known exception — see flags)", ev,
                  flags=("phosphate_rule_may_misfile_a_phosphate_drug",))

    # ---- R8. Lipid / detergent by acyl or alkyl chain. The L44 bug.
    #
    # The chain must DOMINATE the molecule (>= 30% of its heavy atoms) or the
    # molecule must be essentially ring-free. A bare `chain >= 8` test filed
    # abamectin, myxopyronin B and other long-tailed natural-product ANTIBIOTICS
    # as lipids — found on a held-out sample. A fatty acid, an acylglycerol and
    # a detergent are mostly chain; a macrolide is mostly ring.
    chain_frac = (chain / heavy) if (heavy and chain) else 0.0
    ev["chain_fraction_of_heavy_atoms"] = round(chain_frac, 3)
    if chain >= LIPID_MIN_CHAIN and (len(g.rings) <= 1 or chain_frac >= 0.30):
        return _v(cid, rec, "lipid_or_detergent",
                  f"longest unbranched aliphatic carbon chain is {chain} "
                  f"(>= {LIPID_MIN_CHAIN}) and accounts for {chain_frac:.0%} of "
                  f"{heavy} heavy atoms across {len(g.rings)} ring(s): a fatty "
                  "acid, acylglycerol or alkyl detergent. This is the 625 Da "
                  "diacylglycerol class (L44) that cleared the old "
                  "18-heavy-atom floor", ev)

    # ---- R9. Sugars not caught by `type`.
    if ev["sugar_rings"] >= 1 and not ev["has_nitrogen"] and not ev["n_aromatic_rings"]:
        return _v(cid, rec, "sugar_or_glycan",
                  f"{ev['sugar_rings']} pyranose/furanose ring(s) with a "
                  "hydroxyl belt and no aromatic ring or nitrogen: a sugar or "
                  "glycan", ev)
    if ev["sugar_rings"] >= 2:
        return _v(cid, rec, "sugar_or_glycan",
                  f"{ev['sugar_rings']} linked pyranose/furanose rings: an "
                  "oligosaccharide/glycan", ev)

    # ---- R10. Peptides written as one component.
    if ev["peptide_residues"] >= 3:
        return _v(cid, rec, "peptide_or_polymer",
                  f"{ev['peptide_residues']} alpha-amino-acid backbone units: a "
                  "peptide. May be a real binder but it is not small-molecule "
                  "evidence", ev, flags=("standalone_peptide_ligand",))

    # ---- R11. Bench chemistry: PEGs, polyols, buffers, cryoprotectants.
    if ev["ether_oxygens"] >= 3 and not ev["has_nitrogen"] and not g.rings \
            and set(els) <= {"C", "H", "O"}:
        return _v(cid, rec, "crystallisation_additive",
                  f"acyclic C/H/O chain with {ev['ether_oxygens']} ether "
                  "oxygens: a polyethylene glycol", ev)
    # ---- R11b. Acyclic polyol / polyether, AT ANY SIZE and with nitrogen
    # allowed. Two measured false positives, both from an IL-17A Foldseek run,
    # missed the three rules around this one by a hair:
    #   * `B3P` bis-tris propane — 19 heavy atoms, 6 hydroxyls, ring-free. The
    #     `heavy <= 14 and not rings` rule misses it by 5 atoms and the
    #     `heavy <= 18` polyol rule by ONE.
    #   * `JEF` Jeffamine — 41 heavy atoms, 9 ether oxygens, ring-free. The PEG
    #     rule above misses it only because a Jeffamine carries a terminal
    #     amine, so `set(els) <= {C,H,O}` fails.
    # Neither is a size question. A molecule with NO ring system at all, no
    # amide, and four or more hydroxyl/ether oxygens is a polyol, a polyether
    # or a Tris-family buffer, at 19 heavy atoms or at 41. Drug-like ligands
    # essentially always carry a ring; the ones that do not are peptides and
    # lipids, and both are decided well above this line.
    if not g.rings and ev["n_amide_bonds"] == 0 \
            and (ev["hydroxyls"] + ev["ether_oxygens"]) >= 4:
        return _v(cid, rec, "crystallisation_additive",
                  f"no ring system at all, no amide, and "
                  f"{ev['hydroxyls']} hydroxyl(s) + {ev['ether_oxygens']} ether "
                  "oxygen(s): an acyclic polyol, polyether or Tris-family "
                  "buffer (bis-tris propane / Jeffamine / PEG class)", ev)
    # Good's buffers. An ALKYL sulfonate (S with >=3 O, on an sp3 carbon) with
    # no aromatic ring is MES/HEPES/MOPS/PIPES/CHES chemistry. Aryl sulfonamide
    # drugs are not touched: their sulfur carries two oxygens and a nitrogen.
    if ev["alkyl_sulfonates"] >= 1 and aromatic == 0 and heavy is not None \
            and heavy <= 20:
        return _v(cid, rec, "crystallisation_additive",
                  f"{ev['alkyl_sulfonates']} alkyl sulfonate group(s) (direct "
                  f"C-S bond), no aromatic ring, {heavy} heavy atoms: a Good's "
                  "buffer (MES/HEPES/MOPS/PIPES/CHES class)", ev)
    # Nothing a screen reports as a hit is both ring-free and this small. Catches
    # polyamines (spermidine, spermine), N-oxalylglycine, small polyacids.
    if heavy is not None and heavy <= 14 and not g.rings:
        return _v(cid, rec, "crystallisation_additive",
                  f"{heavy} heavy atoms and no ring system at all: an additive, "
                  "polyamine or small metabolite, not a reported screening hit", ev)
    if heavy is not None and heavy <= 18 and ev["n_aromatic_rings"] == 0 \
            and (ev["hydroxyls"] + ev["ether_oxygens"] + ev["alkyl_sulfonates"]) >= 3 \
            and ev["n_amide_bonds"] == 0:
        return _v(cid, rec, "crystallisation_additive",
                  f"{heavy} heavy atoms, no aromatic ring, "
                  f"{ev['hydroxyls']} hydroxyls / {ev['ether_oxygens']} ethers / "
                  f"{ev['alkyl_sulfonates']} sulfonates and no amide: a polyol, "
                  "sulfonate buffer (MES/HEPES class) or cryoprotectant", ev)

    # ---- R12. Size ceiling.
    if mw is not None and mw > DRUGLIKE_MAX_MW:
        return _v(cid, rec, "unknown",
                  f"{mw:.0f} Da exceeds the {DRUGLIKE_MAX_MW:.0f} Da "
                  "small-molecule ceiling and no cofactor/lipid/peptide "
                  "signature fired; not classified rather than assumed", ev,
                  confidence="low")

    # ---- R13. Drug-like by exclusion of every endogenous signature, with a
    # positive check: drug-like ligands are overwhelmingly N-containing or at
    # least aromatic. A pure C/H/O acyclic molecule that got this far is greasy
    # bench chemistry, not a drug.
    if not ev["has_nitrogen"] and ev["n_aromatic_rings"] == 0 and set(els) <= {"C", "H", "O"}:
        return _v(cid, rec, "crystallisation_additive",
                  "pure C/H/O, no nitrogen and no aromatic ring: not drug-like "
                  "chemistry", ev, confidence="medium")

    # ---- R14. UBIQUITY. The one place a measured LIST is the right
    # instrument, and the argument for it is narrow on purpose.
    #
    # `BEN` (benzamidine, 120 Da, 9 heavy atoms) is a ubiquitous protease
    # crystallisation additive AND a bona fide fragment: its close neighbours
    # are real thrombin and trypsin inhibitors. No structural test can reject
    # it without rejecting a real ligand class, which is a different situation
    # from ADP, where chemistry separates cleanly and a list never could.
    #
    # So the discriminator is not what the molecule looks like but how it
    # BEHAVES across the PDB: a component that turns up in hundreds of
    # unrelated entries is laboratory practice, not precedent. That is a query
    # (`SELECT comp_id, count(distinct entry_id) FROM pdb_v.entry_ligands`),
    # not an opinion, and `ChemCompSource.with_entry_counts()` runs it.
    #
    # BOUNDED TWO WAYS, because an unbounded frequency prior would eat real
    # chemistry. It fires only (a) below `UBIQUITY_MAX_HEAVY_ATOMS`, the size
    # band where additives and fragments actually overlap, so it can never
    # demote a 500 Da inhibitor, and (b) at the very end, so it can only ever
    # convert `druglike` and never overrides a chemistry verdict.
    #
    # MEASURED, and note what the measurement also says: BEN 361 entries,
    # B3P 232, JEF 21, LZ1 (a genuine fragment hit) 10, ZBR 9, LFI 8, MOV 7,
    # N5S 1 — against GOL 26004 and EDO 17548. A frequency prior ALONE would
    # not have worked: JEF is a real additive at 21 entries, below any cut that
    # keeps LZ1 safe. JEF and B3P are caught by chemistry (R11b) and only BEN
    # needs this. The better statistic — spread across unrelated proteins
    # rather than entries — is not usable here: the
    # `entry_ligands x pdb_chains` join times out at 120 s while the plain
    # count returns in 617 ms.
    n_entries = ev.get("n_pdb_entries")
    if n_entries is not None and heavy is not None \
            and heavy <= UBIQUITY_MAX_HEAVY_ATOMS \
            and n_entries >= UBIQUITY_MIN_ENTRIES:
        return _v(cid, rec, "crystallisation_additive",
                  f"{heavy} heavy atoms and present in {n_entries} PDB entries "
                  f"(>= {UBIQUITY_MIN_ENTRIES}): a component this small that "
                  "appears across hundreds of unrelated entries is a "
                  "crystallisation additive, whatever its chemistry looks "
                  "like. This is the benzamidine class, which no structural "
                  "test can separate from a real fragment hit", ev,
                  flags=("ubiquity_prior_applied",), confidence="medium")

    bits = []
    if els.get("N"):
        bits.append(f"{els['N']} nitrogen(s)")
    bits.append(f"{ev['n_aromatic_rings']} aromatic ring(s)")
    bits.append(f"{heavy} heavy atoms")
    if mw is not None:
        bits.append(f"{mw:.0f} Da")

    dl_flags: list[str] = []
    conf = "high"
    # A component carrying THREE OR MORE alkyl-halide electrophiles is a
    # bifunctional/trifunctional crosslinking reagent, not a drug — see
    # `SmilesGraph.electrophilic_halide_carbons` for why the threshold is 3 and
    # not 2 (nitrogen mustards are approved drugs with exactly two). The
    # chemistry cannot settle it on its own, so this does NOT change the
    # verdict; it lowers confidence and says what to go and check. With a
    # `StructureContext` the question is answered outright — `LFI` in 8QFZ
    # becomes `polymer_conjugate`.
    if ev["electrophilic_halide_carbons"] >= 3:
        dl_flags.append("multi_electrophile_may_be_a_crosslinking_reagent")
        conf = "medium"
    if n_entries is None and heavy is not None \
            and heavy <= UBIQUITY_MAX_HEAVY_ATOMS:
        # Said out loud rather than left as a silent pass: at this size the
        # additive/fragment ambiguity is real and it was not checked.
        dl_flags.append("ubiquity_not_checked")

    return _v(cid, rec, "druglike",
              "no cofactor, lipid, sugar, peptide, polymer or additive "
              "signature fired; " + ", ".join(bits)
              + ("; NOTE: 3+ alkyl-halide electrophiles — check "
                 "_struct_conn before treating this as a ligand"
                 if "multi_electrophile_may_be_a_crosslinking_reagent" in dl_flags
                 else ""),
              ev, flags=tuple(dl_flags), confidence=conf)


# --------------------------------------------------------------------------
# The context rules. Applied ON TOP of the chemistry verdict, never instead
# of it — `classify_record` stays a pure function of the CCD row.
# --------------------------------------------------------------------------


def _apply_structure_context(v: LigandVerdict, context: "StructureContext | None"
                             ) -> LigandVerdict:
    """Re-attribute a component using the entry it appears in.

    WHY A NEW VERDICT RATHER THAN A FLAG ON `druglike`, OR RE-USING
    `peptide_or_polymer`. Three arguments, in order of weight:

    1. **A flag fails open, and this module exists because fail-open is how the
       four historical bugs happened.** Four call sites test
       `verdict == "druglike"` — `is_druglike_ligand`, `filter_druglike`,
       `holo_call`, and `pocket-scan`'s `_ligands`. A hard flag that the holo
       call "must respect" is a flag that one of them will not read, and the
       failure mode when it is missed is silent and flattering. A distinct
       verdict is respected by every one of those four without any of them
       changing, because none of them equals `"druglike"`.
    2. **`peptide_or_polymer` would be chemically false and would break the
       invariant that a verdict is a property of the component.** `LFI` is not
       a peptide residue and is not a polymer; saying so would mislead the next
       reader. Worse, the same comp_id would carry the verdict
       `peptide_or_polymer` in 8QFZ and `druglike` in some future entry where
       it sits free — the same label meaning two different things. A separate
       verdict makes it explicit that a DIFFERENT question was answered, and
       `evidence["chemistry_verdict"]` keeps the chemistry answer readable.
    3. **It names the finding and carries the thing the dossier actually
       needs.** `polymer_conjugate` means "this component is a covalent
       constituent of a polymer, and here is which one" —
       `evidence["conjugate_of"]` holds that entity's description, length and
       sequence, which is exactly what a PEPTIDE-PRECEDENT block is made of.
       The peptide in 8QFZ is real evidence about the target; it is just not
       small-molecule evidence. Downgrading to a flag on `druglike` would keep
       the modality error; downgrading to `peptide_or_polymer` would throw the
       pointer away.

    THE RULES, and the control each one has to survive.

    C0 — the component is a MONOMER of a polymer entity (`_entity_poly_seq`).
         Then it is a residue, not a ligand: `peptide_or_polymer`. Catches the
         `NH2`/`ACE` capping groups.

    C1 — covalently bonded to a polymer entity that is NOT the target.
         `polymer_conjugate`. This is `LFI` in 8QFZ (three bonds to the Bicycle
         peptide, entity 2, while the target TSLP is entity 1).

    C2 — covalently bonded to TWO OR MORE polymer residues.
         `polymer_conjugate`, and this one needs no target identity at all. A
         covalent DRUG carries one warhead and makes one bond; two or more
         bonds means the component is stapling a chain together, i.e. it is
         part of that chain's covalent constitution rather than a thing bound
         to it. Measured: `LFI` 3 bonds, `ZBR` (TBMB) 3, `A1I4O` 3, `8VY`
         (bis(bromomethyl)benzene) 2 — and 8VY's two bonds are to Cys427 and
         Cys432 of the SAME chain in 5V2P, a crosslinked protein, which C1
         cannot catch because that chain may well be the target.

    C3 — exactly one covalent bond, to the TARGET polymer. **The chemistry
         verdict is preserved untouched** and a flag is added. This is the
         control that stops the fix trading one false positive for a worse
         false negative: 6OIM's `MOV` (sotorasib, one bond to KRAS Cys12) and
         4G5J's `0WN` (afatinib, one bond to EGFR Cys797) both stay `druglike`.

    C4 — exactly one covalent bond, to a polymer, target identity UNKNOWN.
         Undecidable, so nothing is changed except `confidence` -> `medium` and
         a flag. Refusing to guess here is the difference between C3 and C1 and
         it is the only place this fix could manufacture a false negative.
    """
    if context is None:
        return v
    cid = v.comp_id
    ev = dict(v.evidence)
    ev["structure_context"] = context.to_dict()

    # ---- C0. A monomer of a polymer chain is a residue, not a ligand.
    mono = context.is_polymer_monomer(cid)
    if mono is not None and v.verdict != "peptide_or_polymer":
        pe = context.polymer_entities.get(mono)
        ev["chemistry_verdict"] = v.verdict
        ev["chemistry_reason"] = v.reason
        ev["polymer_monomer_of"] = pe.to_dict() if pe else {"entity_id": mono}
        return _v(cid, _rec_of(v), "peptide_or_polymer",
                  f"{context.entry_id or 'this entry'} lists {cid} in "
                  f"_entity_poly_seq for polymer entity {mono}"
                  + (f" ({pe.description})" if pe and pe.description else "")
                  + ": it is a RESIDUE of that chain, not a bound ligand. "
                  f"Chemistry alone called it {v.verdict!r}",
                  ev, flags=tuple(v.flags) + ("polymer_monomer",),
                  confidence="high", source=v.source)

    if not context.is_available():
        # SAID OUT LOUD. A context with no `_struct_conn` category cannot
        # answer the covalent question, and the answer it appears to give —
        # "no linkages" — is indistinguishable from the real thing. This is
        # the assembly-file trap; see `has_struct_conn_category`.
        ev["context_rule"] = "not_applied_no_struct_conn_category"
        return _v(cid, _rec_of(v), v.verdict,
                  v.reason + ". NOTE: the structural context supplied carries "
                  "no _struct_conn category (an assembly file does not), so "
                  "the covalent-attribution rules did NOT run. This is not "
                  "evidence that the component is unbonded — fetch the entry "
                  "or header mmCIF",
                  ev, flags=tuple(v.flags) + ("struct_conn_absent_from_context",),
                  confidence=("medium" if v.confidence == "high" else v.confidence),
                  source=v.source)

    links = context.links_for(cid)
    poly_links = [l for l in links if l.partner_is_polymer]
    if not poly_links:
        return v

    ev["covalent_links"] = [asdict(l) for l in poly_links]
    partners = sorted({l.partner_entity_id or "?" for l in poly_links})
    ev["covalently_bonded_to_entities"] = partners
    ev["n_covalent_bonds_to_polymer"] = len(poly_links)

    def conjugate(reason: str, rule: str, confidence: str = "high") -> LigandVerdict:
        e = dict(ev)
        e["chemistry_verdict"] = v.verdict
        e["chemistry_reason"] = v.reason
        e["context_rule"] = rule
        hosts = [context.polymer_entities[p].to_dict()
                 for p in partners if p in context.polymer_entities]
        e["conjugate_of"] = hosts
        # The thing a dossier has to record INSTEAD of a small-molecule holo
        # call. Rule 1 of the dossier: modality first, always.
        e["precedent_modality"] = _precedent_modality(hosts)
        return _v(cid, _rec_of(v), "polymer_conjugate", reason, e,
                  flags=tuple(v.flags) + ("covalent_to_polymer_ligand",),
                  confidence=confidence, source=v.source)

    non_target = [l for l in poly_links if l.partner_is_target is False]
    on_target = [l for l in poly_links if l.partner_is_target is True]

    # ---- C1. Bonded to a polymer that is not the target.
    if non_target and not on_target:
        host = ", ".join(
            f"{p}"
            + (f" ({context.polymer_entities[p].description})"
               if p in context.polymer_entities
               and context.polymer_entities[p].description else "")
            for p in sorted({l.partner_entity_id or "?" for l in non_target})
        )
        return conjugate(
            f"covalently bonded ({len(non_target)} _struct_conn covale "
            f"linkage(s)) to polymer entity {host}, which is NOT the target "
            f"({context.target_basis or 'target given by caller'}). The "
            "component is a covalent constituent of that polymer's assembly, "
            f"not an independent ligand. Chemistry alone called it "
            f"{v.verdict!r}; that chemistry is not disputed, its ATTRIBUTION "
            "is — precedent here belongs to the polymer, not to a small "
            "molecule",
            "C1_bonded_to_non_target_polymer")

    # ---- C2. Two or more bonds to a polymer: a crosslink, not a warhead.
    if len(poly_links) >= 2:
        chains = sorted({l.partner_chain or "?" for l in poly_links})
        return conjugate(
            f"{len(poly_links)} covalent linkages to polymer residues "
            f"({', '.join(f'{l.partner_comp_id} {l.partner_chain}/{l.partner_seq_id}' for l in poly_links)}"
            f") across chain(s) {', '.join(chains)}: a CROSSLINKER or "
            "macrocyclisation reagent, not a ligand. A covalent drug carries "
            "one warhead and makes one bond; two or more means the component "
            "is stapling the chain rather than binding it",
            "C2_multivalent_covalent_crosslink",
            confidence="high" if not on_target else "medium")

    # ---- C3 / C4. Exactly one bond.
    only = poly_links[0]
    if only.partner_is_target is True:
        ev["chemistry_verdict"] = v.verdict
        ev["context_rule"] = "C3_single_covalent_bond_to_target"
        return _v(cid, _rec_of(v), v.verdict,
                  v.reason + f". Covalently bonded to the TARGET polymer "
                  f"({only.partner_comp_id} {only.partner_chain}/"
                  f"{only.partner_seq_id}) by one linkage — a covalent "
                  "inhibitor, which IS evidence of a bindable site and is "
                  "deliberately left untouched",
                  ev, flags=tuple(v.flags) + ("covalent_to_target",),
                  confidence=v.confidence, source=v.source)
    ev["chemistry_verdict"] = v.verdict
    ev["context_rule"] = "C4_single_covalent_bond_target_unknown"
    return _v(cid, _rec_of(v), v.verdict,
              v.reason + f". Covalently bonded to polymer entity "
              f"{only.partner_entity_id} by ONE linkage, and the target entity "
              "could not be resolved, so bonded-to-the-target cannot be told "
              "from bonded-to-a-partner. Verdict left unchanged and confidence "
              "lowered rather than guessed; pass uniprot_accession or "
              "target_entity_ids to decide it",
              ev, flags=tuple(v.flags) + ("covalent_to_unidentified_polymer",),
              confidence="medium", source=v.source)


def _precedent_modality(hosts: list[dict[str, Any]]) -> str:
    """What modality the HOST polymer's precedent is, per dossier rule 1.

    The point of the whole context fix: 8QFZ is real evidence — a 12-residue
    bicyclic peptide binds this groove — and it belongs in the peptide block,
    not the small-molecule one.
    """
    if not hosts:
        return "unknown"
    kinds = set()
    for h in hosts:
        t = (h.get("poly_type") or "").lower()
        n = h.get("n_monomers") or 0
        if "polypeptide" in t or "peptide" in t:
            kinds.add("peptide" if n and n <= 50 else "protein")
        elif "polyribonucleotide" in t or "polydeoxyribonucleotide" in t:
            kinds.add("nucleic_acid")
        else:
            kinds.add("other")
    return "/".join(sorted(kinds))


def _rec_of(v: LigandVerdict) -> dict[str, Any]:
    """Rebuild the minimal CCD row `_v` reads, so a re-verdict keeps its
    name/formula/type/smiles without the caller passing them again."""
    return {"name": v.name, "formula": v.formula, "formula_weight": v.mw,
            "type": v.comp_type, "smiles": v.smiles,
            "drugbank_id": v.drugbank_id}


def _classify_without_smiles(cid, rec, ev, heavy, els, ctype) -> LigandVerdict:
    """Degraded path: the CCD row exists but carries no usable SMILES.

    Only the checks that need nothing but the formula and `type` are allowed to
    fire. Everything else returns `unknown` — deliberately, because this is the
    exact code path where a guess becomes a fabricated holo structure.
    """
    ev = dict(ev)
    ev["degraded_no_smiles"] = True
    if heavy is not None and heavy <= TRIVIAL_MAX_HEAVY_ATOMS and not els.get("N"):
        return _v(cid, rec, "crystallisation_additive",
                  f"no SMILES in the CCD row; {heavy} heavy atoms and no "
                  "nitrogen is bench chemistry on formula alone", ev,
                  confidence="medium")
    return _v(cid, rec, "unknown",
              "the CCD row has no SMILES, so no chemistry test can run. "
              "Reported as unknown rather than assumed drug-like — an "
              "unclassified ligand must not become a holo structure", ev,
              confidence="low")


# --------------------------------------------------------------------------
# Public API
# --------------------------------------------------------------------------


def classify_ligand(comp_id: str, *, source: str = "pdb",
                    chemcomps: ChemCompSource | None = None,
                    context: StructureContext | None = None) -> LigandVerdict:
    """Classify one chemical component.

    `comp_id` is the FULL component ID from the mmCIF — `A1JPS`, not the first
    three characters of it. Five-character codes have been issued since 2023 and
    the legacy PDB format cannot hold them; that truncation is a separate,
    already-documented wrong answer on IL-17A.

    `source` currently accepts only `'pdb'` (the PDB Chemical Component
    Dictionary via Paperclip). It exists so a future CCD mirror can be selected
    without changing call sites.
    """
    if source != "pdb":
        raise ValueError(f"unsupported source {source!r}; only 'pdb' is implemented")
    return classify_ligands([comp_id], source=source, chemcomps=chemcomps,
                            context=context)[comp_id.upper()]


def classify_ligands(comp_ids: Iterable[str], *, source: str = "pdb",
                     chemcomps: ChemCompSource | None = None,
                     context: StructureContext | None = None
                     ) -> dict[str, LigandVerdict]:
    """Batch form — ONE Paperclip round trip per 40 comp_ids.

    Call sites process whole entries, and an entry has 1-15 components. Looping
    `classify_ligand` would be one subprocess each.
    """
    if source != "pdb":
        raise ValueError(f"unsupported source {source!r}; only 'pdb' is implemented")
    ids = [c.upper() for c in comp_ids if c]
    src = chemcomps or _default_source()
    recs = src.get_many(ids)
    out: dict[str, LigandVerdict] = {}
    for c in ids:
        rec = recs.get(c)
        if rec is None and c in src.fetch_errors:
            # A LOOKUP FAILURE, NOT A CCD MISS. Distinct verdict text and a
            # distinct flag, so a call site can retry or refuse to report
            # rather than quietly counting the entry as apo.
            out[c] = _v(c, None, "unknown",
                        "the chemical-component lookup FAILED (not: the "
                        f"component is absent). {src.fetch_errors[c]}. This "
                        "entry's holo/apo state is undetermined and must not "
                        "be reported as apo",
                        {"heavy_atoms": None, "ccd_hit": None,
                         "lookup_error": src.fetch_errors[c]},
                        flags=("lookup_failed",), confidence="none",
                        source="paperclip:pdb_v.chemcomps (lookup failed)")
        else:
            out[c] = _apply_structure_context(classify_record(rec, c), context)
    return out


def is_druglike_ligand(comp_id: str, *, chemcomps: ChemCompSource | None = None,
                       context: StructureContext | None = None) -> bool:
    """True only for `druglike`. `unknown` is False — an unclassified ligand is
    not evidence of a bindable site. `polymer_conjugate` is False as well: a
    crosslinker inside a peptide ligand is that peptide's chemistry."""
    return classify_ligand(comp_id, chemcomps=chemcomps,
                           context=context).verdict == "druglike"


def filter_druglike(comp_ids: Iterable[str], *,
                    chemcomps: ChemCompSource | None = None,
                    context: StructureContext | None = None) -> list[str]:
    verdicts = classify_ligands(comp_ids, chemcomps=chemcomps, context=context)
    return [c for c, v in verdicts.items() if v.verdict == "druglike"]


def holo_call(comp_ids: Iterable[str], *,
              chemcomps: ChemCompSource | None = None,
              context: StructureContext | None = None) -> dict[str, Any]:
    """Entry-level holo/apo call with the full reasoning attached.

    Returns `is_holo`, the drug-like ligands that justify it, and every other
    ligand bucketed by verdict — so a dossier can say "apo, but carrying GDP"
    rather than "apo" full stop, and can show WHY a rejected ligand was
    rejected.
    """
    verdicts = classify_ligands(comp_ids, chemcomps=chemcomps, context=context)
    buckets: dict[str, list[str]] = {}
    for c, v in verdicts.items():
        buckets.setdefault(v.verdict, []).append(c)
    dl = sorted(buckets.get("druglike", []))
    failed = sorted(c for c, v in verdicts.items() if "lookup_failed" in v.flags)
    # WHAT THE ENTRY IS EVIDENCE OF, WHEN IT IS NOT EVIDENCE OF A SMALL
    # MOLECULE. A `polymer_conjugate` is not an absence — it says a peptide or
    # another polymer is bound here, which is real and belongs in the dossier's
    # PEPTIDE precedent block under rule 1, not in the small-molecule one.
    # 8QFZ is a demonstrated ligandable groove with a 12-residue bicyclic
    # peptide against it; what it is not is a small-molecule holo structure.
    precedent: list[dict[str, Any]] = []
    for c in sorted(buckets.get("polymer_conjugate", [])):
        ev = verdicts[c].evidence
        for host in ev.get("conjugate_of", []) or []:
            precedent.append({
                "via_comp_id": c,
                "modality": ev.get("precedent_modality"),
                "entity_id": host.get("entity_id"),
                "description": host.get("description"),
                "n_monomers": host.get("n_monomers"),
                "sequence": host.get("sequence"),
            })
    return {
        "is_holo": bool(dl),
        "druglike_ligands": dl,
        "polymer_conjugates": sorted(buckets.get("polymer_conjugate", [])),
        "polymer_ligand_precedent": precedent,
        "context_applied": context is not None and context.is_available(),
        "by_verdict": {k: sorted(v) for k, v in sorted(buckets.items())},
        "unknown_ligands": sorted(buckets.get("unknown", [])),
        # `is_holo=False` with a non-empty `undetermined` is NOT an apo call.
        # A caller that reports it as apo has reintroduced the original bug in
        # a new place.
        "undetermined": failed,
        "determined": not failed,
        "verdicts": {c: v.to_dict() for c, v in verdicts.items()},
        "flags": sorted({f for v in verdicts.values() for f in v.flags}),
    }
