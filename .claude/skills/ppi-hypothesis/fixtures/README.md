# ppi-hypothesis fixtures

Fifteen cases, one row each, with the column that matters: **why that case is in
the set**. Six of them exist because a plausible version of this skill would not
have run them and would have shipped a broken gate.

Reproduce the whole table without touching a GPU:

```bash
PPI_RUNS=fixtures/runs $PROTO_PY gate_panel.py
$PROTO_PY make_references.py      # re-derives reference_footprints.json; asserts 8DYG == 97
```

`fixtures/runs/*.json` are the run records (metrics, per-seed confidence,
provenance). The CIFs they were computed from are not shipped — 33 MB — but
every number in `panel_results.json` is recomputable from the records.

## The panel

| case | pair | expected | class | why this case is in the set |
| --- | --- | --- | --- | --- |
| P1 | IL-17A / IL-17RA | positive | deposited (7UWM) | A real, deposited interface the loop must recover. Also the case where ESMFold fails on a positive (4 CA pairs vs 72), which is what removes two-model agreement from the signal list. |
| P2 | TNF-alpha / TNFR2 | positive | deposited (3ALQ) | Second deposited positive, different family. The only positive whose predicted interface is *smaller* than the crystal — everything else over-predicts. |
| P3 | TL1A / DcR3 | positive | deposited (3K51, 3MI8) | Same ligand as the real test case, so it fixes the scale for TL1A. Its footprint is the reference every TL1A transfer number is measured against. |
| N1 | TNF-alpha / DcR3 | negative | hard, same superfamily | DcR3 binds FasL, LIGHT and TL1A, **not** TNF. Same fold, same CRD count and same size as TNFR2, which TNF does bind. Fails everything — the easy half of the hard class. |
| N2 | TL1A / TNFR2 | negative | hard, same superfamily | **The negative that nearly broke the gate.** On three seeds it scored concordance 0.797/0.725, inside the positive band; on six it collapsed to 0.536/0.128. This case is why the gate requires two disjoint seed blocks. |
| N3 | IL-17A / TNFR2 | negative | cross-family receptor | Cytokine against a receptor ectodomain of an unrelated family — the shape a naive candidate set would generate. |
| N4 | TNF-alpha / IL-17RA | negative | cross-family receptor | Mirror of N3 with the larger receptor, to check the result is not a size effect. |
| N5 | TL1A / IL-17RA | negative | cross-family receptor | Same, on the unknown case's ligand. Produced 76 CA pairs and 4160 Å² — more than any positive. |
| N6 | IL-17A / carbonic anhydrase II | negative | cross-compartment | A cytosolic enzyme against a secreted cytokine. No pathway, no compartment, similar size to a receptor ectodomain. |
| N7 | TL1A / carbonic anhydrase II | negative | cross-compartment | Same on the unknown case's ligand. Buried 2448 Å², comparable to the deposited TL1A/DcR3 complex at 2049 Å². |
| N8 | TL1A / DR3 **death domain** | negative | wrong domain, same protein | Same protein as the real test case, wrong domain, wrong side of the membrane — and the only part of DR3 that is deposited. Scored ipTM 0.800, above TNF-alpha's scrambled floor. Tests whether the loop reads the construct or the name. |
| N9 | TL1A / DR3 **scrambled** | negative | floor | DR3 ectodomain residues shuffled. Same length, same composition, no fold. Produced **166 CA pairs and 6569 Å²**, the largest interface in the panel, and 0.84 footprint coverage. This single case removes contact count, buried area and footprint transfer from the gate. |
| N10 | IL-17A / IL-17RA scrambled | negative | floor | Per-ligand floor for IL-17A. ipTM 0.479. |
| N11 | TNF-alpha / TNFR2 scrambled | negative | floor | Per-ligand floor for TNF-alpha. ipTM 0.787 — a 0.31 spread from N10's floor in the same panel, which is why there is no absolute ipTM threshold anywhere in this skill. Also ranked **second** in TNF-alpha's candidate set, above the real DcR3 ectodomain. |
| U1 | **TL1A / DR3** | **unknown** | the real test | TL1A signals through DR3 and no TL1A/DR3 complex is deposited anywhere — DR3's only three PDB entities (5YEV, 5YGP, 5YGS) are its intracellular death domain. Passed every gate, and issued no ask, because `L4` already carries the interaction. |

Replicate blocks (`*r.json`, base seed 142) exist for P1, P2, P3, N2, N9 and U1
— the six cases where the first three seeds were load-bearing.

## Other files

| file | what it is |
| --- | --- |
| `panel_results.json` | the gated table, one object per case, including the gate's per-check values and the consensus epitope |
| `reference_footprints.json` | interfaces measured on the deposited assemblies 7UWM, 3ALQ, 3K51, in UniProt numbering, with the detected numbering offsets |
| `dcr3_dr3_alignment.json` | DcR3 -> DR3 residue map (26.0% identity) used to map the deposited DcR3 contact set onto DR3 numbering |
| `worked_ask.json` | the one ask this run produced: `resolve_link` on `L4`, post-resolution, carrying a PDB census that contradicts a row the graph asserts. Passes all five mechanical gates of `graph_read.py --check-ask` |
| `runs/` | per-case run records |

## Sequence data

`../data/*.fasta` are UniProt canonical sequences. Construct ranges are in
`../constructs.py`, one line of provenance each — crystallised range where one
exists, UniProt topological domain where none does. DR3 is the only receptor in
the second category, and that is the whole reason the case is interesting.
