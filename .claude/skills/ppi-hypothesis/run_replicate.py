"""Replicate the one comparison the whole conclusion turns on.

U1 (TL1A/DR3, the true partner) beat N2 (TL1A/TNFR2, a same-superfamily
non-partner) by 0.050 ipTM in the first run. That margin is about twice the
seed-to-seed spread within a run, which is not obviously more than noise. So
re-run both on a disjoint set of seeds and see whether the ordering survives.
"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/Users/bb/repos/claude-agent-starter/managed/druggability-dossier/.claude/skills/cofold-check")
import ppi_hypothesis as P
from constructs import chains_for
from predict import cofold_complex

for cid, lig, par, seed in [("U1r", "TL1A", "DR3", 142), ("N2r", "TL1A", "TNFR2", 142),
                            ("P3r", "TL1A", "DCR3", 142), ("N9r", "TL1A", "DR3_SCRAM", 142),
                            ("P1r", "IL17A", "IL17RA", 142), ("P2r", "TNFA", "TNFR2", 142)]:
    out = os.path.join(HERE, "runs", f"{cid}.json")
    if os.path.exists(out):
        print(cid, "done"); continue
    seqs, ca, cb = chains_for(lig, par)
    print(f"[{cid}] {lig}/{par} seeds {seed}..{seed+2}", flush=True)
    r = cofold_complex(seqs, n_seeds=3, seed=seed)
    ps = []
    for i, cif in enumerate(r["structures_cif"]):
        open(os.path.join(HERE, "runs", f"{cid}_seed{i}.cif"), "w").write(cif)
        ps.append(P.interface_metrics(cif, ca, cb))
    rec = {"case": cid, "ligand": lig, "partner": par, "base_seed": seed,
           "expected": {"U1r": "unknown", "N2r": "negative", "P3r": "positive",
                         "N9r": "negative", "P1r": "positive", "P2r": "positive"}[cid],
           "class": "seed_replicate", "why_in_set": __doc__.strip(),
           "total_residues": sum(len(s) for s in seqs), "n_seeds": 3,
           "chains_a": ca, "chains_b": cb,
           "structural_confidence": r["structural_confidence"],
           "seed_dispersion": r["seed_dispersion"],
           "per_seed_interface": ps,
           "seed_concordance_ligand_side": P.seed_concordance(ps, "a"),
           "seed_concordance_partner_side": P.seed_concordance(ps, "b"),
           "provenance": r["provenance"]}
    json.dump(rec, open(out, "w"), indent=1, default=str)
    print(f"[{cid}] ipTM {r['structural_confidence']['iptm']['values']}", flush=True)
