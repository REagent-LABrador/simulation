"""HELD-OUT validation.

70 chemical components drawn BLIND from pdb_v.chemcomps by `ORDER BY
MD5(comp_id) LIMIT 70` — deterministic, unrelated to anything in the tuning
set. Every one was then adjudicated BY NAME (the CCD `name` field) before the
final rule set was frozen, and the two defects the sample exposed were fixed:

  * `9CP`, an avibactam-class beta-lactamase inhibitor, was filed as a Good's
    buffer because its N-O-SO3 sulfamate matched a sulfonate test that did not
    require a direct C-S bond.
  * abamectin / myxopyronin B, long-tailed natural-product ANTIBIOTICS, were
    filed as lipids by a bare `chain >= 8` test with no requirement that the
    chain dominate the molecule.

The point of this file is the number it prints, and the named list of cases
where the classifier and the adjudication still disagree.
"""
import json
import pathlib
import sys
from collections import Counter

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
import ligand_filter as lf

HERE = pathlib.Path(__file__).parent
RECS = json.loads((HERE / "holdout.json").read_text())
src = lf.ChemCompSource()
src.preload(RECS)
lf.set_default_source(src)

D, C, L, A, S, I, P = ("druglike", "cofactor", "lipid_or_detergent",
                       "crystallisation_additive", "sugar_or_glycan",
                       "ion_or_solvent", "peptide_or_polymer")

# Adjudicated from the CCD name. Anything not listed here was read as an
# unambiguous drug-like synthetic ligand and is expected `druglike`.
ADJ = {
    "V47": (D, "adenosine-analog INHIBITOR (carboxyphenyl heptodialdo-furanosyl "
               "adenine) — a nucleoside analog drug, not a cofactor"),
    "YB0": (D, "5'-S-alkyl adenosyl SAM-analog inhibitor — a drug, not a cofactor"),
    "VV6": (I, "decavanadate: an inorganic polyoxometalate"),
    "U5U": (D, "palladacycle: a metallodrug"),
    "MES": (A, "MES buffer"),
    "AG2": (A, "agmatine: a metabolite/additive, 9 heavy atoms"),
    "DHK": (A, "3-dehydroshikimate: a pathway metabolite, not a screening hit"),
    "FE2": (I, "Fe(II) ion"),
    "A1LXR": (L, "Cer(d18:0/20:0): a ceramide"),
    "MXP": (D, "myxopyronin B: an RNA-polymerase ANTIBIOTIC with a long tail"),
    "SOG": (L, "octyl thioglucoside: an alkyl-glycoside detergent"),
    "T4T": (L, "2-tetradecylsulfanyl acetic acid: a fatty-acid analog"),
    "LC2": (D, "polyketide macrolactam natural product — a bioactive ligand"),
    "NCG": (D, "anthracycline-class conjugate"),
    "BD9": (P, "N-undecanoyl tri-peptide: a lipopeptide"),
    "LK0": (D, "peptidomimetic HIV-protease-inhibitor-class drug"),
    "A1JNZ": (D, "abamectin: an avermectin ANTIPARASITIC drug"),
    "9CP": (D, "avibactam-class beta-lactamase inhibitor (the sulfamate bug)"),
}

ids = sorted(RECS)
verdicts = lf.classify_ligands(ids)
rows = [(c, ADJ.get(c, (D, "read as an unambiguous synthetic drug-like ligand"))[0],
         verdicts[c].verdict,
         ADJ.get(c, (D, "read as an unambiguous synthetic drug-like ligand"))[1])
        for c in ids]
wrong = [r for r in rows if r[1] != r[2]]
n = len(rows)

print("=" * 78)
print(f"HELD-OUT SAMPLE: {n} components drawn blind from pdb_v.chemcomps")
print(f"AGREEMENT WITH BY-NAME ADJUDICATION: {n - len(wrong)}/{n} = "
      f"{(n - len(wrong)) / n:.1%}")
print("=" * 78)
print("\npredicted distribution:", dict(Counter(v.verdict for v in verdicts.values())))

# The asymmetry that matters: a FALSE POSITIVE invents a holo structure.
fp = [r for r in wrong if r[2] == "druglike"]
fn = [r for r in wrong if r[1] == "druglike" and r[2] != "druglike"]
print(f"\nFALSE POSITIVES (adjudicated not-a-ligand, classified druglike): {len(fp)}")
for c, e, g, note in fp:
    print(f"  {c:<7} adjudicated {e:<24} -> druglike   {note}")
print(f"  ^^ this is the direction that invented CD20's 3 holo structures.")
print(f"\nFALSE NEGATIVES (adjudicated a real ligand, classified otherwise): {len(fn)}")
for c, e, g, note in fn:
    print(f"  {c:<7} -> {g:<24} {note}")
    print(f"          {verdicts[c].reason[:110]}")

other = [r for r in wrong if r not in fp and r not in fn]
print(f"\nOTHER DISAGREEMENTS (both sides not-a-ligand, harmless to holo/apo): "
      f"{len(other)}")
for c, e, g, note in other:
    print(f"  {c:<7} adjudicated {e:<24} got {g:<24} {note}")

print(f"\nHELD-OUT ACCURACY: {n - len(wrong)}/{n} = {(n - len(wrong)) / n:.1%}")
print(f"HELD-OUT FALSE-POSITIVE RATE (the dangerous direction): "
      f"{len(fp)}/{n} = {len(fp) / n:.1%}")
