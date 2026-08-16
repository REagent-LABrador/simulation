"""Compare two ligand_filter versions verdict-by-verdict over the cached CCD rows.

Why this exists: `fixtures/targets.json` pins the classifier that produced every
regenerated pdb_holo/pdb_apo/pdb_undetermined figure by sha256 and line count.
On 2026-08-15 that pin was found to identify a file that existed only in session
scratch, not in this checkout. This script is how you check whether a version
difference actually moves any verdict, rather than assuming either way.

    python3 version_diff.py <old_ligand_filter.py> <new_ligand_filter.py>

Prints every comp_id whose verdict differs, with the CCD name and formula weight
so the difference can be adjudicated by chemistry rather than by trust.
"""
import importlib.util
import json
import pathlib
import sys

HERE = pathlib.Path(__file__).parent


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def verdict(mod, row):
    try:
        return mod.classify_record(row).verdict
    except Exception as exc:                      # a crash IS a difference
        return f"ERROR:{type(exc).__name__}"


def main(old_path, new_path):
    recs = json.loads((HERE / "chemcomps.json").read_text())
    recs.update(json.loads((HERE / "extra_recs.json").read_text()))
    recs.update(json.loads((HERE / "holdout.json").read_text()))

    old = load("_lf_old", old_path)
    new = load("_lf_new", new_path)

    rows = {c: r for c, r in recs.items() if isinstance(r, dict)}
    diff = []
    for cid, row in rows.items():
        o, n = verdict(old, row), verdict(new, row)
        if o != n:
            diff.append((cid, o, n, str(row.get("name"))[:54],
                         row.get("formula_weight")))

    print(f"old: {old_path}")
    print(f"new: {new_path}")
    print(f"components compared: {len(rows)}")
    print(f"VERDICT DIFFERENCES: {len(diff)}\n")
    for cid, o, n, name, mw in sorted(diff):
        print(f"  {cid:<7} {o:<26} -> {n:<26} mw={mw}")
        print(f"          {name}")
    if not diff:
        print("  none — the two versions agree on every cached component, so a "
              "count taken under one reproduces under the other.")
    return 1 if diff else 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1], sys.argv[2]))
