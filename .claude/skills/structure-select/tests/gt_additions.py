"""GROUND-TRUTH ADDITIONS, 2026-08-15. Every one is a MEASURED wrong answer or
its control, not an invented case.

Two blocks, because two different questions are being asked.

CHEM  — expected verdict from the CCD row alone (`classify_record`). This is
        where the three IL-17A Foldseek false positives go.
CTX   — expected verdict for a comp_id AS IT APPEARS IN A NAMED ENTRY, with a
        `StructureContext` built from that entry's mmCIF. A crosslinker has no
        context-free right answer, so putting one in CHEM would be inventing a
        label the classifier is not allowed to reach.
"""
D, C, L, A, S, I, P, PC, U = (
    "druglike", "cofactor", "lipid_or_detergent", "crystallisation_additive",
    "sugar_or_glycan", "ion_or_solvent", "peptide_or_polymer",
    "polymer_conjugate", "unknown")

CHEM = {
    # --- IL-17A Foldseek run, 2026-08-15: 8 of 8 holo flags were false ------
    "BEN": (A, "benzamidine: ubiquitous protease-crystallisation additive, 361 "
               "PDB entries. THE REGRESSION — the superseded MW_MIN=250 floor "
               "excluded it and the chemistry classifier does not. It was the "
               "single holo call in IL-17A's single-chain top-25 neighbourhood "
               "(2GNN). Chemistry cannot reach it: benzamidine IS a fragment "
               "and its neighbours are real thrombin/trypsin inhibitors"),
    "B3P": (A, "bis-tris propane, 282 Da buffer. Ring-free, 6 hydroxyls, missed "
               "the heavy<=18 polyol rule by ONE heavy atom"),
    "JEF": (A, "Jeffamine, 598 Da precipitant. Ring-free polyether, missed the "
               "PEG rule only because it carries a terminal amine"),
    # --- crosslinkers, CONTEXT-FREE. The honest chemistry answer is drug-like
    # and the classifier must say so with the reagent flag raised and
    # confidence lowered, NOT invent a verdict it has no basis for.
    "LFI":   (D, "TATA tri-electrophile (8QFZ). Context-free the chemistry IS "
                 "drug-like; the flag is the deliverable, see FLAG_REQUIRED"),
    "ZBR":   (D, "TBMB, the canonical Bicycle crosslinker (3QN7)"),
    "A1I4O": (D, "1,3,5-tris(2-chloroethylsulfonyl)triazinane, TATA-sulfone (9Q8N)"),
    "8VY":   (D, "1,3-bis(bromomethyl)benzene, a two-armed staple (5V2P). Two "
                 "arms, so the >=3 reagent flag deliberately does NOT fire — "
                 "that is the nitrogen-mustard guard working"),
    "260":   (D, "2-(bromomethyl)-1,3-difluorobenzene: ONE alkyl halide plus "
                 "two ARYL fluorides. The control that proves the flag counts "
                 "electrophiles and not halogens"),
    "0WN":   (D, "afatinib (4G5J): a real covalent inhibitor, one warhead"),
}

#: comp_id -> required flags on the context-free verdict.
FLAG_REQUIRED = {
    "LFI":   ("multi_electrophile_may_be_a_crosslinking_reagent",),
    "ZBR":   ("multi_electrophile_may_be_a_crosslinking_reagent",),
    "A1I4O": ("multi_electrophile_may_be_a_crosslinking_reagent",),
}
FLAG_FORBIDDEN = {
    "8VY": ("multi_electrophile_may_be_a_crosslinking_reagent",),
    "260": ("multi_electrophile_may_be_a_crosslinking_reagent",),
    "0WN": ("multi_electrophile_may_be_a_crosslinking_reagent",),
}

#: (comp_id, pdb_id, uniprot_accession or None) -> (expected verdict, why)
CTX = {
    ("LFI", "8QFZ", "Q969D9"): (
        PC, "THE MEASURED FALSE POSITIVE. TATA crosslinker on the three "
            "cysteines of the 12-residue Bicycle peptide 8QFZ_2; target is "
            "TSLP, entity 1. Was druglike/high, tier=holo, "
            "ligand_site_jaccard 0.769 and 1.000"),
    ("LFI", "8B9P", "Q9BYF1"): (
        PC, "same reagent, different target (ACE2) and TWO peptide copies — "
            "6 covale rows. Proves the rule is not one-entry tuning"),
    ("ZBR", "3QN7", "P00749"): (
        PC, "TBMB stapling a bicyclic peptide inhibitor of uPA. Different "
            "reagent chemistry (benzylic, not acyl), same attribution"),
    ("A1I4O", "9Q8N", None): (
        PC, "TATA-sulfone in a bicyclic peptide against LptD, and NO target "
            "accession supplied — so C1 cannot fire and rule C2 has to carry "
            "it alone. This is the target-identity-free path"),
    ("8VY", "5V2P", None): (
        PC, "bis(bromomethyl)benzene crosslinking Cys427 to Cys432 of ONE "
            "chain: a stapled PROTEIN, not a peptide ligand. C1 could never "
            "catch it because that chain may well be the target"),
    ("NH2", "8B9P", "Q9BYF1"): (
        P,  "C-terminal amide cap listed in _entity_poly_seq. One of the "
            "ground-truth set's three standing misses, closed by context"),
    # ---- THE CONTROLS. A covalent inhibitor is bonded to the target too, and
    # trading the false positive for a false negative here would be worse than
    # the bug.
    ("MOV", "6OIM", "P01116"): (
        D,  "sotorasib, ONE covalent bond to KRAS Cys12 — the target. MUST "
            "stay druglike"),
    ("0WN", "4G5J", "P00533"): (
        D,  "afatinib, ONE covalent bond to EGFR Cys797 — the target. MUST "
            "stay druglike"),
    ("MOV", "6OIM", None): (
        D,  "same entry with the accession withheld: rule C4, undecidable, so "
            "the verdict must NOT move and confidence drops to medium"),
}
