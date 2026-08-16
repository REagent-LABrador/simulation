"""Apply the gate to every case in the panel and report the confusion table."""
import json, os, sys, statistics
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import ppi_hypothesis as P, footprint as F
from constructs import CONSTRUCTS

REPL = {"P1": "P1r", "P2": "P2r", "P3": "P3r", "U1": "U1r", "N2": "N2r", "N9": "N9r"}
REF = json.load(open(os.path.join(HERE, "fixtures/reference_footprints.json")))
HOM = {"IL17A": "7UWM", "TNFA": "3ALQ", "TL1A": "3K51"}
DEPOSITED = {"IL17A/IL17RA": "7UWM", "TNFA/TNFR2": "3ALQ", "TL1A/DCR3": "3K51 / 3MI8"}
IN_GRAPH = {"TL1A/DR3": "L4 in g_tl1a1 (TL1A binds DR3, basis primary)"}

def load(c):
    d = json.load(open(os.path.join(os.environ.get("PPI_RUNS", os.path.join(HERE, "runs")), f"{c}.json")))
    v = d["structural_confidence"]["iptm"]["values"]
    ps = d["per_seed_interface"]
    if c in REPL and os.path.exists(os.path.join(os.environ.get("PPI_RUNS", os.path.join(HERE, "runs")), f"{REPL[c]}.json")):
        d2 = json.load(open(os.path.join(os.environ.get("PPI_RUNS", os.path.join(HERE, "runs")), f"{REPL[c]}.json")))
        v += d2["structural_confidence"]["iptm"]["values"]
        ps += d2["per_seed_interface"]
        blocks = 2
    else:
        blocks = 1
    return d, v, ps, blocks

CASES = ["P1","P2","P3","N1","N2","N3","N4","N5","N6","N7","N8","N9","N10","N11","U1"]
rows = {}
for c in CASES:
    d, v, ps, blocks = load(c)
    sc = P.seed_concordance(ps, "a")
    lig, par = d["ligand"], d["partner"]
    start = CONSTRUCTS[lig][1]
    rs = set(REF[HOM[lig]]["ligand_footprint_uniprot"])
    cov = statistics.median([F.compare(m["interface_res_a"], start, rs, "")["ref_coverage"] for m in ps])
    rows[c] = dict(case=c, pair=f"{lig}/{par}", ligand=lig, expected=d["expected"],
                   klass=d["class"], n_seeds=len(v), n_blocks=blocks,
                   iptm_median=round(statistics.median(v), 4),
                   iptm_min=round(min(v), 4), iptm_max=round(max(v), 4),
                   conc_mean=sc["mean_jaccard"], conc_min=sc["min_jaccard"],
                   ca_pairs=statistics.median([m["ca_pairs_8a"] for m in ps]),
                   bsa=round(statistics.median([m["bsa_total_a2"] for m in ps])),
                   fp_transfer_coverage=round(cov, 3),
                   consensus_uniprot=[r + start - 1 for r in sc["consensus_residues"]],
                   why_in_set=d["why_in_set"])

SCRAMBLE = {"TL1A": rows["N9"]["iptm_median"],   # in-run floor, PER LIGAND
            "IL17A": rows["N10"]["iptm_median"],
            "TNFA": rows["N11"]["iptm_median"]}
ranks = {}
for lig in ("IL17A", "TNFA", "TL1A"):
    order = sorted([r for r in rows.values() if r["ligand"] == lig],
                   key=lambda r: -r["iptm_median"])
    for i, r in enumerate(order, 1):
        ranks[r["case"]] = (i, order[0]["expected"] == "positive")

out = []
for c in CASES:
    r = rows[c]
    rank, rank1_known = ranks[c]
    g = P.gate(iptm_median=r["iptm_median"],
               iptm_scramble_median=SCRAMBLE.get(r["ligand"]),
               concordance={"mean_jaccard": r["conc_mean"], "min_jaccard": r["conc_min"]},
               n_seeds=r["n_seeds"], n_seed_blocks=r["n_blocks"],
               rank_in_candidate_set=rank, rank1_is_known_partner=rank1_known,
               already_in_graph=r["pair"] in IN_GRAPH,
               deposited_complex=DEPOSITED.get(r["pair"]))
    r["rank_in_candidate_set"] = rank
    r["gate_passed"] = g.passed
    r["gate_failed_on"] = g.failed
    r["novelty"] = g.novelty
    r["ask_issued"] = bool(g.passed and not DEPOSITED.get(r["pair"]) and r["pair"] not in IN_GRAPH)
    out.append(r)

json.dump(out, open(os.path.join(HERE, "fixtures/panel_results.json"), "w"), indent=1)

print(f"{'case':5}{'pair':17}{'exp':9}{'seeds':6}{'ipTM':8}{'conc':7}{'cmin':7}{'CApr':6}{'BSA':7}{'fpCov':7}{'rk':4}{'GATE':6} why-not")
for r in out:
    print(f"{r['case']:5}{r['pair']:17}{r['expected'][:8]:9}{r['n_seeds']:<6}{r['iptm_median']:<8}"
          f"{str(r['conc_mean']):<7}{str(r['conc_min']):<7}{r['ca_pairs']:<6.0f}{r['bsa']:<7}"
          f"{r['fp_transfer_coverage']:<7}{r['rank_in_candidate_set']:<4}"
          f"{'PASS' if r['gate_passed'] else 'fail':6} {','.join(r['gate_failed_on'])[:44]}")

pos = [r for r in out if r["expected"] == "positive"]
neg = [r for r in out if r["expected"] == "negative"]
print(f"\npositives recovered by the gate : {sum(r['gate_passed'] for r in pos)}/{len(pos)}")
print(f"negatives passing the gate (FPR): {sum(r['gate_passed'] for r in neg)}/{len(neg)}"
      f"  -> 95% upper bound by rule of three: {3/len(neg):.0%}")
hard = [r for r in neg if r["klass"] in ("hard_same_superfamily", "scrambled_sequence",
                                         "wrong_domain_same_protein")]
print(f"  of which HARD negatives        : {sum(r['gate_passed'] for r in hard)}/{len(hard)}"
      f"  -> upper bound {3/len(hard):.0%}, i.e. not measured")
print(f"asks actually issued            : {sum(r['ask_issued'] for r in out)}/{len(out)}")
for r in out:
    if r["gate_passed"]:
        print(f"  {r['case']} {r['pair']}: gate passed, ask={r['ask_issued']} -- {r['novelty'][:110]}")
