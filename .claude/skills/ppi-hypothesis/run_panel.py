"""Run the ppi-hypothesis control panel through cofold-check.

    $PROTO_PY run_panel.py <case-id> [<case-id> ...]
    $PROTO_PY run_panel.py all

One Boltz-2 multimer cofold per case, n_seeds=3, plus one ESMFold call on the
same chains. Results land in runs/<case>.json — the raw payload plus the
per-seed interface metrics this module computes itself, because cofold-check
measures the interface on seed 0 only.
"""
import json
import os
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/Users/bb/repos/claude-agent-starter/managed/druggability-dossier/.claude/skills/cofold-check")

import ppi_hypothesis as P                                  # noqa: E402
from constructs import PANEL, chains_for, CONSTRUCTS        # noqa: E402
from predict import cofold_complex, esmfold_predict         # noqa: E402

N_SEEDS = int(os.environ.get("N_SEEDS", "3"))
RUNS = os.environ.get("PPI_RUNS", os.path.join(HERE, "runs"))
os.makedirs(RUNS, exist_ok=True)


def run_case(case):
    cid, lig, par, expected, klass, why = case
    out_path = os.path.join(RUNS, f"{cid}.json")
    if os.path.exists(out_path):
        print(f"[{cid}] already done, skipping")
        return
    seqs, ca, cb = chains_for(lig, par)
    print(f"[{cid}] {lig}x{len(ca)} + {par}  total {sum(len(s) for s in seqs)} aa  "
          f"chains {ca} vs {cb}  expected={expected}", flush=True)

    rec = {"case": cid, "ligand": lig, "partner": par, "expected": expected,
           "class": klass, "why_in_set": why,
           "constructs": {lig: CONSTRUCTS[lig], par: CONSTRUCTS[par]},
           "n_seeds": N_SEEDS, "chains_a": ca, "chains_b": cb,
           "total_residues": sum(len(s) for s in seqs)}

    t0 = time.time()
    r = cofold_complex(seqs, n_seeds=N_SEEDS)
    rec["boltz2_wall_s"] = round(time.time() - t0, 1)
    rec["structural_confidence"] = r["structural_confidence"]
    rec["seed_dispersion"] = r["seed_dispersion"]
    rec["cofold_interface_seed0"] = r["interface"]
    rec["provenance"] = r["provenance"]

    per_seed = []
    for i, cif in enumerate(r["structures_cif"]):
        with open(os.path.join(RUNS, f"{cid}_seed{i}.cif"), "w") as fh:
            fh.write(cif)
        m = P.interface_metrics(cif, ca, cb)
        m["seed"] = i
        per_seed.append(m)
        print(f"   seed{i}: ca_pairs={m.get('ca_pairs_8a')} "
              f"bsa={m.get('bsa_total_a2')} nres_a={m.get('n_interface_res_a')}", flush=True)
    rec["per_seed_interface"] = per_seed
    rec["seed_concordance_ligand_side"] = P.seed_concordance(per_seed, "a")
    rec["seed_concordance_partner_side"] = P.seed_concordance(per_seed, "b")

    try:
        t0 = time.time()
        e = esmfold_predict(seqs)
        rec["esmfold_wall_s"] = round(time.time() - t0, 1)
        rec["esmfold"] = {k: e[k] for k in ("self_report", "interface") if k in e}
        cif = e.get("structure_cif") or e.get("structures_cif", [None])[0]
        if cif:
            with open(os.path.join(RUNS, f"{cid}_esmfold.cif"), "w") as fh:
                fh.write(cif)
            rec["esmfold_interface_measured"] = P.interface_metrics(cif, ca, cb)
    except Exception as exc:
        rec["esmfold_error"] = f"{type(exc).__name__}: {exc}"
        print(f"   esmfold FAILED: {rec['esmfold_error']}", flush=True)

    with open(out_path, "w") as fh:
        json.dump(rec, fh, indent=1, default=str)
    print(f"[{cid}] done in {rec['boltz2_wall_s']}s -> {out_path}", flush=True)


if __name__ == "__main__":
    want = sys.argv[1:]
    cases = PANEL if want == ["all"] else [c for c in PANEL if c[0] in want]
    if not cases:
        sys.exit(f"no cases matched {want}; ids are {[c[0] for c in PANEL]}")
    for c in cases:
        try:
            run_case(c)
        except Exception:
            traceback.print_exc()
            print(f"[{c[0]}] FAILED, continuing", flush=True)
