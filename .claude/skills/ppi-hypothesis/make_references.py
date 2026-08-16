"""Regenerate fixtures/reference_footprints.json from the deposited assemblies.

Downloads the biological assemblies (never the ASU), measures the interface with
the same code the predictions are measured with, and converts both sides to
UniProt numbering via an offset detected by string-matching a 25-mer -- never
assumed. 3ALQ's TNF chain A is offset +76 and its receptor chain R is +22; a
guessed offset would silently fabricate a footprint overlap.

    $PROTO_PY make_references.py
"""
import json, os, sys, urllib.request
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import footprint as F, ppi_hypothesis as P

DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)

# entry: (ligand chain group, partner chain, ligand acc, partner acc,
#         chain to detect the ligand offset on, chain for the partner offset)
REFS = {
    "7UWM": (["A", "B"], ["C"], "Q16552", "Q96F46", "A", "C"),
    "3ALQ": (["A", "B", "C"], ["R"], "P01375", "P20333", "A", "R"),
    "3K51": (["A", "A-2", "A-3"], ["B"], "O95150", "O95407", "A", "B"),
}

out = {}
for pdb, (ca, cb, acc_a, acc_b, off_a, off_b) in REFS.items():
    path = os.path.join(DATA, f"{pdb}-assembly1.cif")
    if not os.path.exists(path):
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb}-assembly1.cif", path)
    m = P.interface_metrics(open(path).read(), ca, cb, with_bsa=True)
    oa = F.detect_offset(path, off_a, acc_a)
    ob = F.detect_offset(path, off_b, acc_b)
    out[pdb] = {
        "ligand_acc": acc_a, "partner_acc": acc_b,
        "offset_ligand": oa, "offset_partner": ob,
        "ca_pairs_8a": m["ca_pairs_8a"], "bsa_total_a2": m["bsa_total_a2"],
        "ligand_footprint_uniprot": sorted(F.to_uniprot(m["interface_res_a"], None, oa)),
        "partner_footprint_uniprot": sorted(F.to_uniprot(m["interface_res_b"], None, ob)),
    }
    print(pdb, "CA pairs", m["ca_pairs_8a"], "BSA", m["bsa_total_a2"],
          "offsets", oa, ob)

json.dump(out, open(os.path.join(HERE, "fixtures/reference_footprints.json"), "w"), indent=1)

# Convention control: 8DYG must return exactly 97 CA-CA pairs. If it does not,
# the contact code is measuring something other than the station's reference.
p8 = os.path.join(DATA, "8DYG-assembly1.cif")
if not os.path.exists(p8):
    urllib.request.urlretrieve("https://files.rcsb.org/download/8DYG-assembly1.cif", p8)
n = P.interface_metrics(open(p8).read(), ["A"], ["B"], with_bsa=False)["ca_pairs_8a"]
print(f"8DYG A/B CA-CA pairs = {n} (must be 97)")
assert n == 97, "contact convention broken -- do not trust any other number here"
