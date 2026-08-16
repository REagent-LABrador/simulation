"""Constructs used in the ppi-hypothesis control panel.

Every range is a UniProt range, and every range for a positive control is the
range that was actually crystallised in the reference entry (read off
`rcsb_polymer_entity_align`, not guessed). DR3 has no complex to read a range
from, so its range is the UniProt extracellular topological domain, which is
the same definition the other three ectodomains satisfy.

Using the deposited construct for a positive and a guessed one for the unknown
would make any difference between them uninterpretable, so the rule is:
crystallised range where one exists, UniProt domain boundary where none does,
and say which per entry.
"""
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, "data")


def _seq(acc: str) -> str:
    with open(os.path.join(_DATA, f"{acc}.fasta")) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))


def _frag(acc: str, start: int, end: int) -> str:
    return _seq(acc)[start - 1:end]


CONSTRUCTS = {
    # name: (accession, start, end, provenance, role)
    "IL17A":  ("Q16552", 25, 155, "7UWM entity 1 aligned region Q16552 25-155", "ligand"),
    "IL17RA": ("Q96F46", 33, 304, "7UWM entity 2 aligned region Q96F46 33-304", "receptor"),
    "TNFA":   ("P01375", 77, 233, "3ALQ entity 1 aligned region P01375 77-233 (soluble TNF)", "ligand"),
    "TNFR2":  ("P20333", 33, 205, "3ALQ entity 2 aligned region P20333 33-205 (CRD1-4)", "receptor"),
    "TL1A":   ("O95150", 72, 251, "3K51 entity 1 aligned region O95150 72-251 (secreted form)", "ligand"),
    "DCR3":   ("O95407", 30, 195, "3K51 entity 2 aligned region O95407 30-195 (CRD1-4)", "receptor"),
    "DR3":    ("Q93038", 25, 199, "UniProt Q93038 topological domain Extracellular 25-199; "
                                  "CRD1-4 are 34-192. NO complex exists to read a range from.", "receptor"),
    "CA2":    ("P00918", 2, 260, "UniProt P00918 chain 2-260, cytosolic carbonic anhydrase II", "control"),
}

STOICHIOMETRY = {"IL17A": 2, "TNFA": 3, "TL1A": 3}  # biological oligomer of the ligand


def chains_for(ligand: str, partner: str) -> tuple[list[str], list[str], list[str]]:
    """(sequence list, chain ids of the ligand group, chain ids of the partner).

    Boltz-2 names output chains A, B, C ... in input order, so the ligand
    oligomer occupies the first N and the partner the next one.
    """
    n = STOICHIOMETRY.get(ligand, 1)
    lig = _frag(*CONSTRUCTS[ligand][:3])
    par = (scrambled(partner[:-6]) if partner.endswith("_SCRAM")
           else _frag(*CONSTRUCTS[partner][:3]))
    seqs = [lig] * n + [par]
    ids = [chr(ord("A") + i) for i in range(n + 1)]
    return seqs, ids[:n], ids[n:]


# ---------------------------------------------------------------------------
# The panel. `expected` is what a human who knows the biology would say.
# ---------------------------------------------------------------------------
PANEL = [
    # id, ligand, partner, expected, class, why this pair is in the set
    ("P1", "IL17A", "IL17RA", "positive", "deposited",
     "7UWM. Real interface, deposited, in Boltz-2's training set."),
    ("P2", "TNFA", "TNFR2", "positive", "deposited",
     "3ALQ. Real interface, deposited."),
    ("P3", "TL1A", "DCR3", "positive", "deposited",
     "3K51/3MI8. Real interface, deposited. Same ligand as the unknown case."),

    ("N1", "TNFA", "DCR3", "negative", "hard_same_superfamily",
     "DcR3 binds FasL, LIGHT and TL1A, NOT TNF-alpha. Same fold, same size and "
     "same CRD architecture as TNFR2, which TNF does bind. This is the negative "
     "the gate has to survive."),
    ("N2", "TL1A", "TNFR2", "negative", "hard_same_superfamily",
     "TL1A binds DR3 and DcR3, NOT TNFR2. Mirror image of N1."),

    ("N3", "IL17A", "TNFR2", "negative", "cross_family_receptor",
     "A secreted cytokine against a receptor ectodomain of an unrelated family."),
    ("N4", "TNFA", "IL17RA", "negative", "cross_family_receptor",
     "Mirror of N3, with the larger receptor."),
    ("N5", "TL1A", "IL17RA", "negative", "cross_family_receptor",
     "Same, on the ligand the unknown case uses."),

    ("N6", "IL17A", "CA2", "negative", "cross_compartment",
     "Cytosolic enzyme against a secreted cytokine. Different compartment, no "
     "shared pathway, similar size to a receptor ectodomain."),
    ("N7", "TL1A", "CA2", "negative", "cross_compartment",
     "Same, on the unknown case's ligand."),

    ("U1", "TL1A", "DR3", "unknown", "the_real_test",
     "TL1A signals through DR3. No TL1A/DR3 complex is deposited anywhere; DR3's "
     "only PDB entries are its intracellular death domain. This is the one pair "
     "in the panel whose answer we do not have."),
]


# ---------------------------------------------------------------------------
# Second wave. Both were added AFTER the first ten cases showed that contact
# count and buried area cannot separate a cytosolic enzyme from a real
# receptor. They test the two things the first wave could not:
#   N8 — does the model know WHICH DOMAIN of the receptor binds? DR3's death
#        domain is intracellular and is the only part of DR3 in the PDB.
#   N9 — does anything about the partner's actual sequence matter, or will the
#        model dock a composition-matched random string?
# ---------------------------------------------------------------------------
import random as _random

CONSTRUCTS["DR3_DD"] = ("Q93038", 322, 415,
                        "UniProt Q93038 Death domain 332-413, extended to 322-415 to "
                        "match 5YGP/5YGS. Intracellular. Cannot contact a secreted "
                        "ligand. The ONLY part of DR3 that is deposited.", "control")


def scrambled_dr3(seed: int = 20260815) -> str:
    """DR3 ectodomain residues, shuffled. Same length, same composition, no fold."""
    acc, s, e, _p, _r = CONSTRUCTS["DR3"]
    chars = list(_frag(acc, s, e))
    _random.Random(seed).shuffle(chars)
    return "".join(chars)


PANEL_2 = [
    ("N8", "TL1A", "DR3_DD", "negative", "wrong_domain_same_protein",
     "DR3's intracellular death domain against TL1A. Same protein as the real "
     "test case, wrong domain, wrong side of the membrane. If the loop cannot "
     "tell DR3's ectodomain from DR3's death domain it is not reading the "
     "sequence, it is reading the name."),
    ("N9", "TL1A", "DR3_SCRAM", "negative", "scrambled_sequence",
     "DR3 ectodomain sequence shuffled — identical amino-acid composition, no "
     "fold, no evolutionary signal. The floor. Anything that passes here "
     "invalidates every other number in the panel."),
]

CONSTRUCTS["DR3_SCRAM"] = ("Q93038", 25, 199,
                           "DR3 ectodomain 25-199 with residues shuffled (seed 20260815). "
                           "Composition-matched random sequence.", "control")
PANEL = PANEL + PANEL_2


def scrambled(name: str, seed: int = 20260815) -> str:
    """Any construct's residues, shuffled. The in-run floor control."""
    acc, s, e = CONSTRUCTS[name][:3]
    chars = list(_frag(acc, s, e))
    _random.Random(seed).shuffle(chars)
    return "".join(chars)


for _n in ("IL17RA", "TNFR2"):
    CONSTRUCTS[_n + "_SCRAM"] = CONSTRUCTS[_n][:3] + (
        f"{_n} construct with residues shuffled (seed 20260815). "
        "Composition-matched random sequence — the in-run floor.", "control")

PANEL_3 = [
    ("N10", "IL17A", "IL17RA_SCRAM", "negative", "scrambled_sequence",
     "Scramble floor for the IL-17A candidate set. Every ligand needs its own; "
     "an ipTM has no scale without one."),
    ("N11", "TNFA", "TNFR2_SCRAM", "negative", "scrambled_sequence",
     "Scramble floor for the TNF-alpha candidate set."),
]
PANEL = PANEL + PANEL_3
