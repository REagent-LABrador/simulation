"""Full re-validation: ground truth + additions + context set + held-out."""
import json, pathlib, sys, runpy, io, contextlib
from collections import Counter
sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
sys.path.insert(0, str(pathlib.Path(__file__).parent))
import ligand_filter as lf, offline, gt_additions as GA
lf.ChemCompSource = offline._Offline

HERE = pathlib.Path(__file__).parent
recs = json.loads((HERE/"chemcomps.json").read_text())
recs.update(json.loads((HERE/"extra_recs.json").read_text()))
try:
    counts = json.loads((HERE/"entry_counts.json").read_text())
except FileNotFoundError:
    counts = {}
for c, n in counts.items():
    if c in recs and isinstance(recs[c], dict):
        recs[c]["n_pdb_entries"] = n
src = offline._Offline(); src.preload(recs); lf.set_default_source(src)

# ---- the original 262, verbatim from test_ligand_filter.py -----------------
mod = {"__file__": str(HERE/"test_ligand_filter.py"), "lf": lf}
src_gt = (HERE/"test_ligand_filter.py").read_text().split("# ------------------------------------"
          "---------------------------- RUN")[0].replace("import ligand_filter as lf", "pass")
exec(compile(src_gt, "gt", "exec"), mod)
GT = {k: v for k, v in mod["GT"].items()}
GT_ONLY = {k: v for k, v in GT.items() if k in recs}

rows = []           # (comp_id, label, expected, got, ok, block)
def run(block, cid, expected, note, context=None, key=None):
    v = lf.classify_ligand(cid, chemcomps=src, context=context)
    rows.append((key or cid, block, expected, v.verdict, v.verdict == expected, note, v))
    return v

for cid, (exp, note) in sorted(GT_ONLY.items()):
    run("gt262", cid, exp, note)
for cid, (exp, note) in sorted(GA.CHEM.items()):
    run("gt_add", cid, exp, note)
for (cid, pid, acc), (exp, note) in sorted(GA.CTX.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2] or "")):
    ctx = lf.StructureContext.from_mmcif_path(HERE/"structures"/f"{pid}.cif",
                                              target_accession=acc)
    run("context", cid, exp, note, context=ctx,
        key=f"{cid}@{pid}" + ("" if acc else "[no-acc]"))

# ---- held-out, unchanged --------------------------------------------------
ho = json.loads((HERE/"holdout.json").read_text())
for c, n in counts.items():
    if c in ho: ho[c]["n_pdb_entries"] = n
hsrc = offline._Offline(); hsrc.preload(ho)
hmod = {"__file__": str(HERE/"test_holdout.py"), "lf": lf}
src_ho = (HERE/"test_holdout.py").read_text()
src_ho = src_ho[:src_ho.index("\nrows =") if "\nrows =" in src_ho else src_ho.index("ADJ = {")+src_ho[src_ho.index("ADJ = {"):].index("\n}\n")+3]
exec(compile(src_ho.replace("import ligand_filter as lf","pass"), "ho", "exec"), hmod)
ADJ = hmod["ADJ"]
hrows = []
for cid in sorted(ho):
    exp = ADJ.get(cid, ("druglike", "unlisted -> drug-like by adjudication"))[0]
    v = lf.classify_ligand(cid, chemcomps=hsrc)
    hrows.append((cid, exp, v.verdict, v.verdict == exp))

# ---- report ---------------------------------------------------------------
def block(name, rs):
    ok = sum(1 for r in rs if r[4]); n = len(rs)
    print(f"\n=== {name}: {ok}/{n} = {100*ok/n:.1f}%")
    for r in rs:
        if not r[4]:
            print(f"   MISS {r[0]:10s} expected {r[2]:24s} got {r[3]:24s} | {r[5][:70]}")

for b in ("gt262", "gt_add", "context"):
    block(b, [r for r in rows if r[1] == b])
allr = rows
ok = sum(1 for r in allr if r[4])
print(f"\n=== COMBINED GROUND TRUTH: {ok}/{len(allr)} = {100*ok/len(allr):.1f}%")

# flags
print("\n--- flag assertions")
for cid, need in GA.FLAG_REQUIRED.items():
    v = lf.classify_ligand(cid, chemcomps=src)
    print(f"   {cid:6s} flags={list(v.flags)} conf={v.confidence} "
          f"{'OK' if all(f in v.flags for f in need) else 'FAIL'}")
for cid, bad in GA.FLAG_FORBIDDEN.items():
    v = lf.classify_ligand(cid, chemcomps=src)
    print(f"   {cid:6s} flags={list(v.flags)} "
          f"{'OK' if not any(f in v.flags for f in bad) else 'FAIL'}")

hok = sum(1 for r in hrows if r[3])
fp = [r for r in hrows if not r[3] and r[2] == "druglike"]
print(f"\n=== HELD-OUT: {hok}/{len(hrows)} = {100*hok/len(hrows):.1f}%   "
      f"false positives {len(fp)}/{len(hrows)}")
for r in fp: print("   FP", r)

# confusion matrix over the combined ground truth
print("\n=== CONFUSION MATRIX (combined ground truth: rows = expected, cols = got)")
labels = sorted({r[2] for r in allr} | {r[3] for r in allr})
cm = Counter((r[2], r[3]) for r in allr)
w = max(len(l) for l in labels) + 1
print(" " * w + "".join(f"{l[:9]:>10s}" for l in labels))
for e in labels:
    print(f"{e:<{w}s}" + "".join(f"{cm.get((e,g),0):>10d}" for g in labels))
