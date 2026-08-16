"""Structural modality classification from SMILES.

Why this exists
---------------
`chembl.molecule_dictionary.molecule_type` is the authoritative modality field
for *approved drugs* (11/11 populated on JAK1's approved rows). It is NOT
usable as the authoritative field for *bioactivity compounds*. Measured across
twelve targets, 59.2% of compounds (41,358 of 69,824) carry no usable value,
and the abstention is concentrated on the potent end: on IL-17A the 26
confidently-typed compounds top out at pchembl 6.26 while the 91 it abstains on
reach 9.10.

Worse, the field cannot express "peptide" at all -- there is no `Peptide` value
in the enum -- so oral peptides land in `Small molecule`. ICOTROKINRA
(molregno 3283615), an oral IL-23R peptide, is typed `Small molecule`.

So modality for compounds is decided from structure, and `molecule_type` is
kept only as corroboration. This module is that decision.

2D only. No conformer embedding, no force-field optimisation, nothing needed
from geometry.

Usage
-----
    from modality import classify_smiles, classify_compound

    classify_smiles("CC(C)C[C@@H]1NC(=O)...")   # -> Verdict
    classify_compound(smiles, molecule_type, structure_type)

Run `python modality.py --selftest` for the calibration set.

Measured performance
--------------------
256 cases: 98 IL-17C macrocyclic peptides, 117 IL-17A compounds, 16 named
drugs chosen to be hard, 22 ChEMBL oligosaccharides, 3 drugs with no SMILES.

    coarse (small_molecule vs not)   256/256
    false small-molecule calls         0/256   <- the error that matters
    missed small molecules             0/256
    fine-grained label               255/256

Known limitations, named
------------------------
1. **Depsipeptides.** ROMIDEPSIN is the one fine-grained miss. Its backbone
   alternates amide and *ester*, so only 2 alpha linkages are found and
   structure alone calls it `small_molecule`. It is caught only because ChEMBL
   types it `Protein`. **A depsipeptide that ChEMBL typed `Small molecule`
   would pass as a small molecule.** No case of that is known in the fixture
   set; it is an untested hole, not a measured pass.
2. **No oligonucleotide control was retrievable.** Every ChEMBL
   `Oligonucleotide` row has a SMILES of 1,150-1,818 characters and the
   Paperclip transport truncates long text fields (see SKILL.md). The
   oligonucleotide rule is therefore written but **unverified against real
   data** -- treat an `oligonucleotide` call as untested.
3. **No SMILES means no structural call.** ICOTROKINRA is the reason
   `Small molecule` + `structure_type = NONE` returns `modality_unknown` rather
   than a modality, and that is a disclosure, not a classification.
4. Thresholds are calibrated on immunology and oncology targets. The
   alpha-linkage gap (3 -> 12) is wide enough that it is unlikely to be
   target-specific, but it has not been checked outside the fixture set.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# --------------------------------------------------------------------------
# Modality labels. These are the values that may appear in the dossier.
# --------------------------------------------------------------------------
SMALL_MOLECULE = "small_molecule"
PEPTIDE = "peptide"
MACROCYCLIC_PEPTIDE = "macrocyclic_peptide"
OLIGONUCLEOTIDE = "oligonucleotide"
OLIGOSACCHARIDE = "oligosaccharide"
PROTEIN = "protein_or_antibody"
UNKNOWN = "modality_unknown"

#: Only this label counts toward `target_precedent`. Everything else is
#: biologic/other precedent or is disclosed as unknown. Nothing is silently
#: dropped.
COUNTS_AS_SMALL_MOLECULE = frozenset({SMALL_MOLECULE})

# --------------------------------------------------------------------------
# SMARTS
# --------------------------------------------------------------------------

# An internal alpha-amino-acid linkage: N-CA-C(=O)-N.
# The alpha carbon must be sp3 and both nitrogens must be amine-type (the
# `!$(N=*)` guard rejects imines and amidines; aromatic ring N is excluded by
# NX3 with the aromatic perception RDKit applies).
# One match == one residue that both accepts and donates a backbone amide, so a
# linear peptide of R residues gives R-1 matches.
#
# Calibration: this count is bimodal with an empty gap on real data. IL-17A's
# 117 compounds give {0:94, 1:5, 2:5, 3:2, 12:6, 13:1, 14:3, 24:1} -- nothing
# between 3 and 12 -- and IL-17C's 98 macrocyclic peptides all sit at 12-15.
_ALPHA_RESIDUE = Chem.MolFromSmarts("[NX3;!$(N=*)][CX4][CX3](=[OX1])[NX3;!$(N=*)]")

# Backbone amide alone (C(=O)-N with an sp3 alpha carbon on the acyl side).
_BACKBONE_AMIDE = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[NX3;!$(N=*)]")

# Beta/gamma amino acid linkage -- peptoid and beta-peptide backbones.
_BETA_RESIDUE = Chem.MolFromSmarts("[NX3;!$(N=*)][CX4][CX4][CX3](=[OX1])[NX3;!$(N=*)]")

# Nucleic acid pieces.
_PHOSPHODIESTER = Chem.MolFromSmarts("[OX2,OX1-]P(=O)([OX2,OX1-])[OX2]")
_THIOPHOSPHATE = Chem.MolFromSmarts("[OX2,OX1-]P(=[SX1])([OX2,OX1-,NX3])[OX2]")
_NUCLEOBASE = [
    Chem.MolFromSmarts(s)
    for s in (
        "c1ncc2[nX3]cnc2n1",  # purine core
        "c1nc2c(n1)c(=O)[nX3]cn2",  # guanine-ish
        "O=c1[nX3]cccn1",  # pyrimidinone
        "O=c1[nX3]c(=O)ccn1",  # uracil/thymine
        "Nc1ccn([CX4])c(=O)n1",  # cytosine
    )
]
_FURANOSE = Chem.MolFromSmarts("[CX4;R][OX2;R][CX4;R][CX4;R][CX4;R]")

# Sugar ring: pyranose/furanose bearing multiple exocyclic oxygens.
_SUGAR_RING = [
    Chem.MolFromSmarts(
        "[OX2;R]1[CX4;R][CX4;R]([OX2])[CX4;R]([OX2])[CX4;R]([OX2])[CX4;R]1"
    ),  # pyranose
    Chem.MolFromSmarts("[OX2;R]1[CX4;R][CX4;R]([OX2])[CX4;R]([OX2])[CX4;R]1"),  # furanose
]
_GLYCOSIDIC = Chem.MolFromSmarts("[CX4;R]([OX2;R])[OX2][CX4;R]")

# --------------------------------------------------------------------------
# Thresholds. Every one of these is calibrated on the sets named in --selftest,
# not chosen by taste. See the block comment on each.
# --------------------------------------------------------------------------

#: >= this many alpha-residue linkages makes it a peptide. Calibrated on
#: IL-17C's 98 macrocyclic peptides (min 9 linkages) against IL-17A's 117
#: compounds (max 1) and BCL-2/EGFR/KRAS small molecules. The gap between the
#: two populations is 1 -> 9, so 4 sits in empty space and is not tuned to
#: either edge.
MIN_PEPTIDE_RESIDUE_LINKS = 4

#: A peptide whose largest ring is at least this size is macrocyclic. 12 is the
#: conventional macrocycle floor; IL-17C's macrocycles are 30-60+ atoms.
MACROCYCLE_RING_SIZE = 12

#: Fraction of heavy atoms belonging to the peptide backbone. Guards the
#: residue count against a small molecule that happens to string amides
#: together. IL-17C peptides sit at 0.41-0.51; small molecules with 4+ amide
#: linkages sit far below.
#:
#: Glycopeptides are the documented exception and are handled separately, not
#: by lowering this: VANCOMYCIN sits at 0.248 and ORITAVANCIN at 0.20 because
#: the appended sugars dilute the backbone by construction. Tuning the floor
#: to admit them would be fitting to two points; the glycopeptide rule below
#: keys on the sugar instead, which is the actual cause.
MIN_PEPTIDE_BACKBONE_FRACTION = 0.25

#: Small-molecule sanity bounds. A structure outside these is not called a
#: small molecule on structure alone -- it is reported unknown with the reason,
#: because it is more likely a peptide/conjugate the residue test missed than a
#: genuine 2 kDa small molecule.
MAX_SMALL_MOLECULE_MW = 1500.0
MAX_SMALL_MOLECULE_HEAVY_ATOMS = 100

#: Above this it is a polymer of some kind, whatever the residue test says.
PROTEIN_MW = 5000.0
PROTEIN_RESIDUE_LINKS = 40

#: Oligonucleotide: this many phosphodiester/thiophosphate linkages AND this
#: many nucleobases.
MIN_OLIGO_PHOSPHATES = 2
MIN_OLIGO_BASES = 2

#: Oligosaccharide: this many glycosidically-linked sugar rings.
MIN_SUGAR_RINGS = 3


@dataclass
class Verdict:
    """One modality call, with everything needed to audit it."""

    modality: str
    basis: str  # "structure" | "molecule_type" | "none"
    confidence: str  # "high" | "medium" | "low"
    reason: str
    features: dict = field(default_factory=dict)
    #: Populated when ChEMBL's molecule_type disagrees with the structure.
    #: Never silently resolved -- rule 1 requires the disagreement be reported.
    disagreement: Optional[str] = None

    def counts_as_small_molecule(self) -> bool:
        return self.modality in COUNTS_AS_SMALL_MOLECULE


# --------------------------------------------------------------------------
# Feature extraction
# --------------------------------------------------------------------------


def _n_matches(mol, patt) -> int:
    if patt is None:
        return 0
    return len(mol.GetSubstructMatches(patt, uniquify=True))


def _backbone_atoms(mol) -> set:
    """Heavy atoms participating in an alpha-amino-acid backbone linkage."""
    atoms = set()
    if _ALPHA_RESIDUE is not None:
        for m in mol.GetSubstructMatches(_ALPHA_RESIDUE, uniquify=True):
            atoms.update(m)
    return atoms


def _max_ring_size(mol) -> int:
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    return max((len(r) for r in rings), default=0)


def _n_sugar_rings(mol) -> int:
    seen = set()
    for patt in _SUGAR_RING:
        if patt is None:
            continue
        for m in mol.GetSubstructMatches(patt, uniquify=True):
            seen.add(frozenset(m))
    return len(seen)


def _n_nucleobases(mol) -> int:
    total = 0
    for patt in _NUCLEOBASE:
        if patt is not None:
            total += _n_matches(mol, patt)
    return total


def structural_features(smiles: str) -> Optional[dict]:
    """Return the 2D feature vector used by the classifier, or None if the
    SMILES will not parse."""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    heavy = mol.GetNumHeavyAtoms()
    bb = _backbone_atoms(mol)
    return {
        "heavy_atoms": heavy,
        "mw": round(Descriptors.MolWt(mol), 1),
        "alpha_residue_links": _n_matches(mol, _ALPHA_RESIDUE),
        "beta_residue_links": _n_matches(mol, _BETA_RESIDUE),
        "backbone_amides": _n_matches(mol, _BACKBONE_AMIDE),
        "amide_bonds": rdMolDescriptors.CalcNumAmideBonds(mol),
        "backbone_fraction": round(len(bb) / heavy, 3) if heavy else 0.0,
        "max_ring_size": _max_ring_size(mol),
        "n_rings": rdMolDescriptors.CalcNumRings(mol),
        "phosphates": _n_matches(mol, _PHOSPHODIESTER)
        + _n_matches(mol, _THIOPHOSPHATE),
        "nucleobases": _n_nucleobases(mol),
        "sugar_rings": _n_sugar_rings(mol),
        "glycosidic": _n_matches(mol, _GLYCOSIDIC),
        "n_fragments": len(Chem.GetMolFrags(mol)),
    }


# --------------------------------------------------------------------------
# The classifier
# --------------------------------------------------------------------------


def classify_smiles(smiles: str) -> Verdict:
    """Decide modality from structure alone.

    Order matters: nucleic acid and sugar polymers are tested before the
    peptide test, and the peptide test before small molecule, because a small
    molecule is defined here as "what is left over and also within bounds".
    """
    f = structural_features(smiles)
    if f is None:
        return Verdict(
            UNKNOWN,
            "none",
            "low",
            "no SMILES, or SMILES did not parse",
            features={},
        )

    # --- polymer of nucleotides -------------------------------------------
    if f["phosphates"] >= MIN_OLIGO_PHOSPHATES and f["nucleobases"] >= MIN_OLIGO_BASES:
        return Verdict(
            OLIGONUCLEOTIDE,
            "structure",
            "high",
            f"{f['phosphates']} phosphodiester linkages and "
            f"{f['nucleobases']} nucleobases",
            f,
        )

    # --- polymer of sugars -------------------------------------------------
    if f["sugar_rings"] >= MIN_SUGAR_RINGS and f["glycosidic"] >= 1:
        return Verdict(
            OLIGOSACCHARIDE,
            "structure",
            "high",
            f"{f['sugar_rings']} glycosidically-linked sugar rings",
            f,
        )

    # --- polymer of amino acids -------------------------------------------
    links = max(f["alpha_residue_links"], f["beta_residue_links"])

    # Glycopeptide: an amino-acid backbone carrying glycosidically-attached
    # sugar. The sugar depresses backbone_fraction below the peptide floor by
    # construction, so it gets its own gate. VANCOMYCIN (6 links, 0.248
    # backbone, 1 sugar ring) and ORITAVANCIN (6 links, 0.20, 1 sugar) are
    # both typed `Small molecule` by ChEMBL and are both glycopeptide
    # antibiotics -- they are the reason this branch exists.
    is_glycopeptide = links >= MIN_PEPTIDE_RESIDUE_LINKS and f["sugar_rings"] >= 1

    is_peptide = is_glycopeptide or (
        links >= MIN_PEPTIDE_RESIDUE_LINKS
        and f["backbone_fraction"] >= MIN_PEPTIDE_BACKBONE_FRACTION
    )

    if is_peptide or f["mw"] >= PROTEIN_MW or links >= PROTEIN_RESIDUE_LINKS:
        if f["mw"] >= PROTEIN_MW or links >= PROTEIN_RESIDUE_LINKS:
            return Verdict(
                PROTEIN,
                "structure",
                "high",
                f"{links} residue linkages, MW {f['mw']} -- a polypeptide, "
                f"not a peptide drug",
                f,
            )
        if f["max_ring_size"] >= MACROCYCLE_RING_SIZE:
            return Verdict(
                MACROCYCLIC_PEPTIDE,
                "structure",
                "high",
                f"{links} alpha-amino-acid linkages "
                f"({f['backbone_fraction']:.0%} of heavy atoms) "
                f"closed into a {f['max_ring_size']}-membered ring",
                f,
            )
        return Verdict(
            PEPTIDE,
            "structure",
            "high",
            f"{links} alpha-amino-acid linkages "
            f"({f['backbone_fraction']:.0%} of heavy atoms)",
            f,
        )

    # --- what is left, if it is within small-molecule bounds ---------------
    if f["mw"] > MAX_SMALL_MOLECULE_MW or (
        f["heavy_atoms"] > MAX_SMALL_MOLECULE_HEAVY_ATOMS
    ):
        return Verdict(
            UNKNOWN,
            "structure",
            "low",
            f"MW {f['mw']} / {f['heavy_atoms']} heavy atoms is outside "
            f"small-molecule bounds and no polymer backbone was detected",
            f,
        )

    # A borderline peptide: enough linkages to worry about, not enough
    # backbone to call. Report it rather than counting it.
    if links >= MIN_PEPTIDE_RESIDUE_LINKS:
        return Verdict(
            UNKNOWN,
            "structure",
            "low",
            f"{links} amino-acid-like linkages but only "
            f"{f['backbone_fraction']:.0%} backbone -- peptidomimetic, "
            f"not classifiable on structure",
            f,
        )

    return Verdict(
        SMALL_MOLECULE,
        "structure",
        "high",
        f"no peptide/nucleotide/sugar backbone; MW {f['mw']}, "
        f"{f['heavy_atoms']} heavy atoms",
        f,
    )


# --------------------------------------------------------------------------
# Combining structure with ChEMBL's field
# --------------------------------------------------------------------------

#: ChEMBL molecule_type -> our labels, for the values that are trustworthy on
#: their own face. Note `Small molecule` is deliberately absent: it is the one
#: value with known false positives (ICOTROKINRA), so it is only accepted when
#: structure agrees or when structure is unavailable AND structure_type is MOL.
_CHEMBL_BIOLOGIC = {
    "Antibody": PROTEIN,
    "Antibody drug conjugate": PROTEIN,
    "Protein": PROTEIN,
    "Enzyme": PROTEIN,
    "Cell": PROTEIN,
    "Vaccine component": PROTEIN,
    "Gene": OLIGONUCLEOTIDE,
    "Oligonucleotide": OLIGONUCLEOTIDE,
    "Oligosaccharide": OLIGOSACCHARIDE,
}


def classify_compound(
    smiles: Optional[str],
    molecule_type: Optional[str] = None,
    structure_type: Optional[str] = None,
) -> Verdict:
    """Modality for one compound or drug.

    Precedence:

    1. A ChEMBL molecule_type naming a biologic class wins outright. Those
       values have no measured false positives and structure is usually absent
       for them anyway.
    2. Otherwise structure decides.
    3. Where structure cannot decide, `Small molecule` is accepted ONLY if
       `structure_type` is MOL or BOTH -- i.e. ChEMBL claims to hold a chemical
       structure for it. `Small molecule` + `structure_type = NONE` is the
       ICOTROKINRA signature: an unverifiable claim, and 5,191 molecules in
       ChEMBL carry it. It returns modality_unknown.
    4. Disagreements are recorded, never silently resolved.
    """
    mt = (molecule_type or "").strip() or None
    st = (structure_type or "").strip() or None

    if mt in _CHEMBL_BIOLOGIC:
        coarse = _CHEMBL_BIOLOGIC[mt]
        # ChEMBL has no `Peptide` value, so every peptide drug it does type
        # lands in `Protein` -- CYCLOSPORINE, OCTREOTIDE, LEUPROLIDE,
        # CARFILZOMIB and ROMIDEPSIN are all `Protein`. Where a structure
        # exists, let it supply the finer label. Both labels sit outside
        # `target_precedent`, so this refines the report without moving any
        # compound into the small-molecule block.
        if coarse == PROTEIN and smiles:
            s = classify_smiles(smiles)
            if s.modality in (PEPTIDE, MACROCYCLIC_PEPTIDE):
                return Verdict(
                    s.modality,
                    "structure",
                    "high",
                    f"{s.reason}; ChEMBL molecule_type = {mt} "
                    f"(ChEMBL has no 'Peptide' value)",
                    s.features,
                )
        return Verdict(
            coarse,
            "molecule_type",
            "high",
            f"ChEMBL molecule_type = {mt}",
            structural_features(smiles) or {},
        )

    v = classify_smiles(smiles or "")

    if v.basis == "structure" and v.modality != UNKNOWN:
        if mt == "Small molecule" and v.modality != SMALL_MOLECULE:
            v.disagreement = (
                f"ChEMBL molecule_type = 'Small molecule' but structure is "
                f"{v.modality}: {v.reason}. Structure wins; report both."
            )
        elif mt and mt not in ("Unknown", "Small molecule"):
            v.disagreement = (
                f"ChEMBL molecule_type = '{mt}' but structure is {v.modality}."
            )
        return v

    # Structure could not decide. Two very different reasons for that, and
    # conflating them is what let VANCOMYCIN and ORITAVANCIN through as small
    # molecules in the first version of this module.
    #
    #   basis == "structure": we READ the structure and it came back
    #       ambiguous -- peptidomimetic, or out of small-molecule bounds.
    #       That is positive evidence against `Small molecule`, not an
    #       absence of evidence. molecule_type must NOT rescue it.
    #   basis == "none": there was no structure to read at all. Only here is
    #       molecule_type the best available witness.
    if v.basis == "structure":
        v.disagreement = (
            f"ChEMBL molecule_type = '{mt}' but the structure is ambiguous: "
            f"{v.reason}. Not counted as a small molecule."
        ) if mt == "Small molecule" else v.disagreement
        return v

    if mt == "Small molecule":
        if st in ("MOL", "BOTH"):
            return Verdict(
                SMALL_MOLECULE,
                "molecule_type",
                "medium",
                "no SMILES retrieved here, but molecule_type = Small molecule "
                f"with structure_type = {st} (ChEMBL holds a molfile). "
                "Retrieve the structure and re-check if this compound matters.",
                v.features,
            )
        return Verdict(
            UNKNOWN,
            "none",
            "low",
            f"molecule_type = 'Small molecule' but structure_type = {st or 'NULL'} "
            "and no parsable SMILES -- the claim cannot be corroborated. "
            "This is the ICOTROKINRA signature (an oral IL-23R peptide typed "
            "'Small molecule'); 5,191 ChEMBL molecules share it.",
            v.features,
        )

    return v


# --------------------------------------------------------------------------
# Aggregation for the dossier
# --------------------------------------------------------------------------


def summarise(verdicts, potency=None) -> dict:
    """Roll a list of Verdicts into the counts the dossier needs.

    `potency` is an optional parallel list of pchembl values, so the caller can
    report best potency PER MODALITY. That is the number that went wrong on
    IL-17A: a single best_potency figure over a mixed pool is not attributable
    to a modality.
    """
    potency = potency or [None] * len(verdicts)
    buckets: dict[str, dict] = {}
    for v, p in zip(verdicts, potency):
        b = buckets.setdefault(v.modality, {"n": 0, "best_pchembl": None})
        b["n"] += 1
        if p is not None and (b["best_pchembl"] is None or p > b["best_pchembl"]):
            b["best_pchembl"] = p
    total = len(verdicts)
    sm = buckets.get(SMALL_MOLECULE, {"n": 0, "best_pchembl": None})
    unk = buckets.get(UNKNOWN, {"n": 0, "best_pchembl": None})
    return {
        "n_compounds": total,
        "by_modality": buckets,
        "small_molecule_count": sm["n"],
        "small_molecule_best_pchembl": sm["best_pchembl"],
        "modality_unknown_count": unk["n"],
        "modality_unknown_best_pchembl": unk["best_pchembl"],
        "modality_unknown_fraction": round(unk["n"] / total, 3) if total else None,
        "n_disagreements": sum(1 for v in verdicts if v.disagreement),
    }


if __name__ == "__main__":
    import sys

    if "--selftest" in sys.argv:
        from modality_selftest import main  # noqa: PLC0415

        sys.exit(main())
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        v = classify_smiles(line)
        print(f"{v.modality}\t{v.confidence}\t{v.reason}")
