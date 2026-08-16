"""Ground-truth evaluation of ligand_filter.

Runs entirely off the cached chemcomps.json (fetched from Paperclip
pdb_v.chemcomps) so it is deterministic and offline.

Ground truth is EXPECTED VERDICT per comp_id, assigned from chemistry knowledge,
NOT from the classifier's output. Every disagreement is printed by name.
"""
import json
import pathlib
import sys
from collections import Counter, defaultdict

sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
import ligand_filter as lf

HERE = pathlib.Path(__file__).parent
RECS = json.loads((HERE / "chemcomps.json").read_text())

src = lf.ChemCompSource()
src.preload(RECS)
lf.set_default_source(src)

D, C, L, A, S, I, P, U = ("druglike", "cofactor", "lipid_or_detergent",
                          "crystallisation_additive", "sugar_or_glycan",
                          "ion_or_solvent", "peptide_or_polymer", "unknown")

# ---------------------------------------------------------------- GROUND TRUTH
GT: dict[str, tuple[str, str]] = {}


def add(bucket, ids, note):
    for i in ids.split():
        GT[i] = (bucket, note)


# --- THE FOUR HISTORICAL FAILURES -------------------------------------------
add(L, "Y01 CLR", "CD20 bug: cholesterol / cholesteryl hemisuccinate, cryo-EM additive")
add(L, "PC1 PEE PGV PEF PCW CDL", "CD20 bug: phosphatidylcholine & phospholipid family")
add(C, "2UK", "KRAS neighbour bug: GppNHp analog, 4PHH")
add(L, "L44", "IL-17A neighbour bug: 625 Da diacylglycerol, 4EC7")
add(C, "ADP ATP", "NLRP3 bug: NACHT-domain nucleotide, 27 heavy atoms")

# --- THE HARDCODED COFACTORS frozenset from pocket-scan/modal_app.py --------
# Read out of modal_app.py; the classifier was never told them.
add(C, "GDP GTP GNP GSP GCP G2P GGZ AMP ANP ACP AGS APC ADX", "modal_app COFACTORS")
GT["GGL"] = (P, "CCD types GGL as 'L-gamma-peptide, C-delta linking' — a peptide residue")
add(C, "CDP CTP CMP UDP UTP UMP TTP TMP TDP IDP ITP", "modal_app COFACTORS")
add(C, "NAD NAI NAP NDP NAJ FAD FMN FDA SAM SAH SFG COA ACO MCA", "modal_app COFACTORS")
add(C, "HEM HEC HEA HDD BCL CLA PLP TPP B12 COB BTN MTA APR PRP 3PG", "modal_app COFACTORS")
add(C, "F6P G6P G1P UPG", "modal_app COFACTORS: sugar phosphates")
add(S, "NAG NDG BMA MAN GAL GLA GLC BGC FUC SIA XYS XYP", "modal_app COFACTORS: sugars")
add(L, "MYR PLM OLA STE DAO D12 LDA LMT CHD", "modal_app COFACTORS: lipids/detergents")

# --- modal_app.py NON_LIGANDS ----------------------------------------------
add(I, "SO4 PO4 NO3 AZI IOD BR CL NA K CA MG MN ZN FE FE2 CU NI CD CO CS RB SR BA HG NH4",
    "modal_app NON_LIGANDS: ions")
add(A, "GOL EDO PEG PG4 1PE P6G MPD ACT ACY CIT FLC TRS EPE MES DMS BME DTT TLA FMT",
    "modal_app NON_LIGANDS: buffers/cryo")
add(U, "UNX UNL UNK", "modal_app NON_LIGANDS: unknown ligand placeholders")
GT["IMD"] = (A, "imidazole buffer, 5 heavy atoms")

# --- neighbour_precedent.py EXCLUDED_LIGANDS not already covered ------------
add(C, "IMP 5GP SF4 FES F3S MTE MGD PAP PPS", "neighbour EXCLUDED: cofactors")
add(P, "GSH", "neighbour EXCLUDED: glutathione tripeptide")
add(C, "UD1", "neighbour EXCLUDED: UDP-GlcNAc")
add(A, "PGE P6G 2PE 7PE 12P 15P XPE P33 PE4 PE8 M2M MPO", "neighbour EXCLUDED: PEGs")
add(A, "MLI MLA TAR BTB OXL SIN DTV DTU IPA URE CAC", "neighbour EXCLUDED: buffers")
GT["BCT"] = (I, "CCD name is BICARBONATE ION — an ion, not a buffer organic")
# GT CORRECTIONS made after reading the CCD `name` field, NOT after seeing the
# classifier's answer. Each is a label I had wrong, verified against pdb_v.chemcomps.
GT["BIS"] = (C, "CCD name is '1,1,5,5-TETRAFLUOROPHOSPHOPENTYLPHOSPHONIC ACID "
                "ADENYLATE ESTER' — an AMP-PCP-class nucleotide analog, not bis-tris")
GT["D1D"] = (A, "CCD name is '(4S,5S)-1,2-DITHIANE-4,5-DIOL' — a reducing agent")
GT["GDS"] = (P, "oxidized glutathione disulfide: two tripeptides, a peptide")
GT["URA"] = (C, "uracil — a nucleobase, endogenous")
add(A, "PIN CXS P4C HEZ HEX BU1 BU3 PDO PGO PGR 1BO MRD MRY",
    "neighbour EXCLUDED: cryo/polyols")
add(L, "TWT OCT", "CCD names them DOCOSANE and N-OCTANE: hydrocarbons")
add(A, "SPD SPM PUT", "neighbour EXCLUDED: polyamines")
add(I, "CO3 SCN WO4 MOO VO4 PER PPV POP AF3 ALF BEF MGF", "neighbour EXCLUDED: inorganics")
add(L, "BOG C8E SDS HTG LMN UND", "neighbour EXCLUDED: detergents")
GT["DEP"] = (A, "CCD name is DIETHYL PHOSPHONATE, 8 heavy atoms — bench chemistry")
add(S, "XYL", "neighbour EXCLUDED: xylose")
add(A, "OGA IPH", "neighbour EXCLUDED: small organics")

# --- TRUE POSITIVES that must NOT be filtered ------------------------------
GT["307"] = (D, "2AZ5 TNF-alpha SPD304 — druglike, flagged promiscuous separately")
GT["MOV"] = (D, "6OIM KRAS G12C covalent inhibitor AMG 510 / sotorasib")
GT["A1JPS"] = (D, "9SQX IL-17A ligand — FIVE-character comp_id")
GT["N5S"] = (D, "5QQE RAC1-Kalirin fragment screen hit")
GT["N64"] = (D, "5QQG RAC1-Kalirin fragment screen hit")
GT["STU"] = (D, "staurosporine, ATP-site kinase inhibitor")
GT["1N1"] = (D, "JAK1 ATP-site inhibitor (4E4L series)")
GT["4WI"] = (D, "kinase inhibitor")
GT["LZ1"] = (D, "small-molecule inhibitor")
GT["CFF"] = (D, "caffeine — small but genuinely drug-like chemistry")

# --- EDGE CASES ------------------------------------------------------------
GT["MSE"] = (P, "selenomethionine — a polymer residue")
GT["SEP"] = (P, "phosphoserine — modified polymer residue")
GT["TPO"] = (P, "phosphothreonine — modified polymer residue")
GT["PTR"] = (P, "phosphotyrosine — modified polymer residue")
GT["MLY"] = (P, "dimethyl-lysine — modified polymer residue")
GT["MLZ"] = (P, "methyl-lysine — modified polymer residue")
GT["M3L"] = (P, "trimethyl-lysine — modified polymer residue")
GT["LYS"] = (P, "lysine — polymer residue")
GT["ALA"] = (P, "alanine — polymer residue")
GT["ABA"] = (P, "aminobutyric acid — polymer residue")
GT["CSX"] = (P, "S-oxy-cysteine — modified polymer residue")
GT["ACE"] = (P, "acetyl capping group")
GT["NH2"] = (P, "amide capping group")
GT["ADN"] = (C, "adenosine — nucleoside cofactor")
GT["ADE"] = (C, "adenine — nucleobase")
GT["GUN"] = (C, "guanine — nucleobase")
GT["THM"] = (C, "thymidine — nucleoside")

GT["HYP"] = (P, "hydroxyproline — polymer residue")
GT["D10"] = (L, "decane — alkane")

GT["LP3"] = (L, "lysophospholipid")
GT["LPP"] = (L, "phospholipid")
GT["LHG"] = (L, "phosphatidylglycerol")
GT["PSC"] = (L, "phosphatidylcholine variant")
GT["POV"] = (L, "POPC phosphatidylcholine")
GT["DGA"] = (L, "diacylglycerol")
GT["FUL"] = (S, "fucose")
GT["RAM"] = (S, "rhamnose")

# steroid DRUGS — the deliberate stress test on the sterol rule
GT["DEX"] = (D, "dexamethasone — steroid DRUG, must survive the sterol rule")
GT["EST"] = (D, "estradiol — steroid hormone ligand")
GT["TES"] = (D, "testosterone — steroid hormone ligand")
GT["AND"] = (D, "androstenedione — steroid hormone ligand")
GT["ASD"] = (D, "androstenedione variant")
GT["STR"] = (D, "progesterone (CCD name) — steroid hormone ligand")
GT["PRD"] = (D, "pyrido[2,3-d]pyrimidine-triamine antifolate (CCD name), a real inhibitor")
GT["HCY"] = (D, "hydrocortisone-class steroid")
GT["4NC"] = (D, "4-nitrocatechol — a genuine COMT-class inhibitor/probe, not an additive")

MISSING_OK = {"HOH", "DOD", "SUC", "TRE", "MAL", "MLR", "CME", "CSD", "CSO",
              "KCX", "SNN", "TDP", "X2R", "Q6D", "JAK", "KIN", "R4A"}

# ---------------------------------------------------------------- RUN
ids = [c for c in sorted(GT) if c in RECS]
absent = [c for c in sorted(GT) if c not in RECS]
verdicts = lf.classify_ligands(ids)

rows = []
for c in ids:
    exp, note = GT[c]
    got = verdicts[c]
    rows.append((c, exp, got.verdict, note, got))

wrong = [r for r in rows if r[1] != r[2]]
n = len(rows)
acc = (n - len(wrong)) / n if n else 0.0

print("=" * 78)
print(f"GROUND-TRUTH SET: {n} components with an expected verdict "
      f"({len(absent)} not present in pdb_v.chemcomps, listed below)")
print(f"ACCURACY: {n - len(wrong)}/{n} = {acc:.1%}")
print("=" * 78)

labels = list(lf.VERDICTS)
cm = defaultdict(Counter)
for _, e, g, _, _ in rows:
    cm[e][g] += 1
w = max(len(x) for x in labels) + 1
print("\nCONFUSION MATRIX  (rows = expected, cols = predicted)\n")
print(" " * (w + 2) + "".join(f"{x[:6]:>7}" for x in labels) + "   tot")
for e in labels:
    if not cm[e]:
        continue
    tot = sum(cm[e].values())
    print(f"{e:<{w}} |" + "".join(f"{cm[e][g] or '.':>7}" for g in labels) + f"{tot:>6}")
print(" " * (w + 2) + "".join(f"{sum(cm[e][g] for e in labels) or '.':>7}" for g in labels))

print("\nPER-CLASS RECALL")
for e in labels:
    tot = sum(cm[e].values())
    if tot:
        print(f"  {e:<26} {cm[e][e]}/{tot} = {cm[e][e] / tot:6.1%}")

print("\nPER-CLASS PRECISION")
for g in labels:
    tot = sum(cm[e][g] for e in labels)
    if tot:
        print(f"  {g:<26} {cm[g][g]}/{tot} = {cm[g][g] / tot:6.1%}")

print(f"\nEVERY MISCLASSIFICATION ({len(wrong)}), NAMED\n" + "-" * 78)
for c, exp, got, note, v in wrong:
    print(f"  {c:<6} expected {exp:<24} got {got:<24}")
    print(f"         {note}")
    print(f"         name    : {v.name}")
    print(f"         formula : {v.formula}  mw {v.mw}")
    print(f"         reason  : {v.reason[:150]}")
    print()

if absent:
    print("NOT IN pdb_v.chemcomps (classifier returns `unknown` by design):")
    print("  " + " ".join(absent))
    for c in absent:
        v = lf.classify_ligand(c)
        assert v.verdict == "unknown", f"{c} should be unknown, got {v.verdict}"
        assert not lf.is_druglike_ligand(c)
    print("  -> all return verdict='unknown' and is_druglike_ligand()==False. OK")

# ------------------------------------------------------- promiscuity flagging
print("\nPROMISCUITY FLAG (advisory, does not change the verdict)")
v307 = verdicts["307"]
print(f"  307: verdict={v307.verdict} flags={v307.flags}")
assert v307.verdict == "druglike"
assert "promiscuity_advisory" in v307.flags

# ------------------------------------------------------- five-char comp_id
print("\nFIVE-CHARACTER comp_id")
v = verdicts["A1JPS"]
print(f"  A1JPS: verdict={v.verdict} heavy={v.heavy_atoms} mw={v.mw}")
assert v.comp_id == "A1JPS" and v.verdict == "druglike"

print("\n" + "=" * 78)
print("THE FOUR HISTORICAL FAILURES, RE-RUN")
print("=" * 78)

CASES = [
    ("CD20 (P11836) — 3 reported holo structures",
     [("6Y92", ["Y01"]), ("6VJA", ["CLR"]), ("7PP0", ["PC1"])], 0),
    ("KRAS fold neighbours — 1 of 25 reported holo (4PHH)",
     [("4PHH", ["2UK", "MG", "CL"])], 0),
    ("IL-17A fold neighbours — 1 of 25 reported holo (4EC7)",
     [("4EC7", ["L44"])], 0),
    ("NLRP3 — ADP-bound NACHT domain called holo",
     [("7ALV", ["ADP"]), ("8EJ4", ["ATP", "MG"])], 0),
]
allok = True
for title, entries, expect_holo in CASES:
    calls = [(e, lf.holo_call(ls)) for e, ls in entries]
    n_holo = sum(1 for _, c in calls if c["is_holo"])
    ok = n_holo == expect_holo
    allok &= ok
    print(f"\n{title}")
    print(f"  -> holo count {n_holo} (expected {expect_holo})  {'PASS' if ok else 'FAIL'}")
    for e, c in calls:
        print(f"     {e}: is_holo={c['is_holo']}  {c['by_verdict']}")
        for cid, vv in c["verdicts"].items():
            print(f"        {cid:<6} {vv['verdict']:<24} {vv['reason'][:88]}")

print("\nCONTROL — a genuinely holo entry must still be holo")
for entry, ligs in [("6OIM", ["MOV", "GDP", "MG"]), ("2AZ5", ["307"]),
                    ("9SQX", ["A1JPS", "SO4"]), ("5QQE", ["N5S", "EDO"])]:
    c = lf.holo_call(ligs)
    print(f"  {entry}: is_holo={c['is_holo']} druglike={c['druglike_ligands']} "
          f"other={ {k: v for k, v in c['by_verdict'].items() if k != 'druglike'} }")
    allok &= c["is_holo"]

print("\nHISTORICAL RE-RUN:", "ALL PASS" if allok else "SOME FAILED")
print(f"\nFINAL MEASURED ACCURACY: {n - len(wrong)}/{n} = {acc:.1%}")
