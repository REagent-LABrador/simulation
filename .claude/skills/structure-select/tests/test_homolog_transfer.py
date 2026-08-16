#!/usr/bin/env python3
"""The four measured anchor cases, reproduced through `homolog_transfer.py`.

    python3 test_homolog_transfer.py

**NETWORK-DEPENDENT, unlike every other harness in this directory.** The
`ligand_filter` harnesses are fully offline on purpose (see `README.md`); this
one is not, because it needs coordinates and there is no way to ship nine
mmCIFs beside a skill without shipping them to the Skills API — `deploy.ts`
zips `.claude/skills/<dir>/` whole. Entries are fetched from RCSB and cached
under `$STRUCTURE_SELECT_CACHE` (default `~/.cache/structure-select/mmcif`).
Chem-comp rows come from RCSB too, for the reason in `RcsbChemComps`.

WHAT IS BEING REPRODUCED
------------------------
An anchor-agreement test built sixteen ligand-free site anchors across four
targets. Four of sixteen found the ligand site, and transferred homolog was the
only anchor kind with a positive record (2 of 3 constructible). These are the
four cases that decided the three guards:

    TNF-alpha   2AZ5 <- 3LKJ  LKJ   accept  (the anchor that works)
    S1PR1       3V2Y <- 4Z34  ON7   accept  (a mediocre fit that still landed)
    NLRP3       7ALV <- 5IRN  ADP   GUARD 1 (a perfect transfer of a cofactor)
    IL-17A      9SQX <- 1M48  FRG   GUARD 3 (IL-2 forced onto IL-17A)

Plus two cases that are not in the four but without which two of the guards
would be untested:

    positive control  2AZ5 <- 2AZ5  307   the pipeline against itself
    guard 2           7ALV <- 7KRZ  BO2   right chain, wrong domain

The control matters more than it looks. Every number below is produced by a
stdlib superposition written for this module, and a transfer pipeline that is
subtly wrong will still produce confident coordinates — which is the entire
failure this file exists to guard against. The control transfers 2AZ5's own
SPD304 from one of its two TNF dimers to the other and must recover it exactly.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import homolog_transfer as ht  # noqa: E402


class Case:
    def __init__(self, label, target, donor, ligand, reference,
                 expect_accept, expect_guard=None, note="", **kw):
        self.label = label
        self.target = target
        self.donor = donor
        self.ligand = ligand
        self.reference = reference
        self.expect_accept = expect_accept
        self.expect_guard = expect_guard
        self.note = note
        self.kw = kw


CASES = [
    Case("CONTROL  2AZ5 <- 2AZ5", "2AZ5", "2AZ5", "307", "307", True,
         note="pipeline against itself: SPD304 moved between 2AZ5's two TNF "
              "dimers. Must come back on top of the other copy.",
         donor_chains=["C", "D"], target_chains=["A", "B"]),
    Case("TNF      2AZ5 <- 3LKJ", "2AZ5", "3LKJ", "LKJ", "307", True,
         note="CD40L's subunit-fracture inhibitor onto TNF-alpha. The one "
              "ligand-free anchor with a positive record."),
    Case("S1PR1    3V2Y <- 4Z34", "3V2Y", "4Z34", "ON7", "ML5", True,
         note="LPA1 onto S1PR1 — a ~45%-identical paralogue, so a mediocre "
              "fit that still lands. Do not read its RMSD as a general floor."),
    Case("NLRP3    7ALV <- 5IRN", "7ALV", "5IRN", "ADP", "ADP", False,
         expect_guard="donor_ligand_druglike",
         note="NOD2's only ligand is ADP. The transfer is EXCELLENT and the "
              "answer is the nucleotide lobe, not the drug site."),
    Case("IL-17A   9SQX <- 1M48", "9SQX", "1M48", "FRG", "A1JPS", False,
         expect_guard="alignment_and_sterics",
         note="IL-2 forced onto IL-17A. 21.59 A away and inside the protein, "
              "produced silently."),
    Case("GUARD2   7ALV <- 7KRZ", "7ALV", "7KRZ", "BO2", "RM5", False,
         expect_guard="domain_attribution",
         note="Bortezomib is on the right LONP1 chain (auth 768-898) and "
              "outside the NACHT-aligned region (auth 506-721)."),
]


def main() -> int:
    chem = ht.RcsbChemComps(
        cache_path=Path(__file__).parent / "chemcomps_transfer.json"
    )
    rows = []
    failures = []
    for c in CASES:
        try:
            r = ht.transfer_homolog_site(
                c.target, c.donor, donor_ligand=c.ligand,
                reference_ligand=c.reference, chemcomps=chem, **c.kw,
            )
        except Exception as exc:                       # noqa: BLE001
            failures.append(f"{c.label}: raised {type(exc).__name__}: {exc}")
            continue
        r.assert_reportable()
        p = r.provenance
        ok = r.accepted == c.expect_accept
        if not c.expect_accept:
            ok = ok and c.expect_guard in r.rejected_by
        if not ok:
            failures.append(
                f"{c.label}: expected accept={c.expect_accept}"
                + (f" via {c.expect_guard}" if c.expect_guard else "")
                + f", got accept={r.accepted} rejected_by={r.rejected_by}"
            )
        rows.append((c, r, p, ok))

    print("=" * 100)
    print("THE FOUR MEASURED CASES, THROUGH THE THREE GUARDS")
    print("=" * 100)
    hdr = (f"{'case':<22} {'TM':>6} {'RMSD':>6} {'nCA':>4} {'bbCl':>5} "
           f"{'scCl':>5} {'refD':>7} {'Jac':>5}  verdict")
    print(hdr)
    print("-" * 100)
    for c, r, p, ok in rows:
        v = r.validation or {}
        verdict = ("ACCEPT" if r.accepted
                   else "REJECT by " + ",".join(r.rejected_by))
        print(f"{c.label:<22} {p['tm_score']:>6.3f} {p['rmsd_a']:>6.2f} "
              f"{p['aligned_length']:>4} {p['backbone_clash_count']:>5} "
              f"{p['sidechain_clash_count']:>5} "
              f"{(v.get('centroid_distance_a') if v.get('centroid_distance_a') is not None else float('nan')):>7.2f} "
              f"{(v.get('contact_shell_jaccard') if v.get('contact_shell_jaccard') is not None else float('nan')):>5.2f}"
              f"  {'OK ' if ok else 'XX '}{verdict}")

    print()
    print("refD = centroid of the transferred ligand to the target's own "
          "reference ligand, in angstrom.")
    print("Jac  = contact-shell Jaccard against that reference ligand.")
    print("Both are VALIDATION, computed after the guards and never fed into "
          "one.")
    print()
    print("NOTE ON refD FOR THE TNF CASE. The anchor-agreement test reports "
          "TNF at 0.00 A, and")
    print("that figure is a POCKET-SELECTION distance from fpocket — the "
          "transferred anchor and the")
    print("SPD304 anchor selected the same pocket, so their centroids "
          "coincide exactly. It is not a")
    print("ligand-centroid distance and this harness cannot produce it: "
          "fpocket lives in pocket-scan.")
    print("What is reproducible here is the guard verdict and the geometry "
          "the guards read.")

    print()
    print("=" * 100)
    print("PER-CASE DETAIL")
    print("=" * 100)
    for c, r, p, ok in rows:
        print(f"\n--- {c.label}   donor ligand {p['donor_ligand']}"
              f" ({p.get('donor_ligand_copies')} cop"
              f"{'y' if p.get('donor_ligand_copies') == 1 else 'ies'})")
        print(f"    {c.note}")
        print(f"    chain map {p.get('chain_map_donor_to_target')}"
              f"   {p.get('tm_normalised_by')}")
        if p.get("chain_map_tie_broken_on"):
            print(f"    tie-break: {p['chain_map_tie_broken_on']}")
        lr = p.get("local_refit") or {}
        print(f"    local re-fit: {lr.get('n_pairs')} pairs within "
              f"{lr.get('radius_a')} A, local RMSD {lr.get('local_rmsd_a')} A")
        for g in r.guards:
            mark = "PASS" if g.passed else "FAIL"
            print(f"    [{mark}] {g.name}: {g.detail}")

    print()
    print("=" * 100)
    if failures:
        print(f"FAILURES: {len(failures)}")
        for f in failures:
            print("  " + f)
        return 1
    print(f"ALL {len(rows)} CASES REPRODUCE AS EXPECTED")
    print("=" * 100)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
