---
name: ppi-hypothesis
description: >
  Cofolds a target against a written-down set of candidate partners, measures
  the predicted interface against in-run controls (a composition-scrambled
  partner, a same-superfamily non-partner, disjoint seed blocks), and decides
  whether the prediction has earned a question back to the upstream evidence
  graph in that graph's own four ask verbs. It does NOT decide that a predicted
  interface is real, does NOT gate on confidence, contact count or buried
  surface area, does NOT invent a fifth ask verb, and does NOT let a pending
  ask block a verdict or license a null. On the measured panel it issued ZERO
  asks in fifteen cases, including the one it was built for — read "The hard
  question" before running it expecting an ask.
---

# ppi-hypothesis

The station consumes hypotheses from an upstream evidence graph. This is the
piece that lets it emit one: cofold the graph's target against candidate
partners, and if a predicted interface survives its controls **and** is absent
from the graph, send it back as a question.

That is the design. What it actually does, measured, is narrower and the
narrowing is the deliverable. **Read "Where the gate came from" and "The hard
question" before you run anything.**

## The loop

```
graph nominates target + mechanism
        |
        v
enumerate candidate partners  ------> write the set down, and what you excluded
        |
        v
cofold target x each candidate  -----> >= 6 seeds, >= 2 disjoint seed blocks
        + a composition-scrambled partner (the in-run floor)
        |
        v
measure the interface yourself -------> cofold-check measures seed 0 only
        |
        v
gate  ---------------------------------> five conjunctive checks, all of them
        |                                 reproducibility or novelty tests.
        |                                 NONE of them establishes truth.
        v
is it new to the graph? ---- no ------> report it in the dossier. NO ASK.
        |
       yes
        v
never-fire list ----------- hit ------> NO ASK.
        |
      clear
        v
ask (one of the graph's four verbs) + one not_found[] entry in the dossier
```

Step 4 is the one that closes the loop and it is also the one that almost never
fires. That is the finding, not a shortfall.

## Setup

```bash
set -a; . <repo>/.env; set +a
export MODAL_PROFILE=rafwiewiora
```

Plain Python import, in process, under the proto-tools interpreter — the same
pattern as `cofold-check`, which this module imports and never modifies.

```bash
$PROTO_PY run_panel.py all      # cofold the panel
$PROTO_PY gate_panel.py         # apply the gate, print the confusion table
```

`ppi_hypothesis.py` needs `gemmi` and `numpy`. `footprint.py` and `align.py`
are stdlib+numpy. Nothing here writes to the upstream graph.

## Procedure

### 1. Pick partners honestly, and write down what you excluded

Cofolding a target against everything is not a plan; it is a multiple-testing
problem with the corrections left out. Three sources, in this order:

| source | example, from this run |
| --- | --- |
| partners the graph already names, for this target or a related one | DR3 (`t6`) is named for TL1A on `L4` |
| members of the cognate receptor or ligand family | TL1A is a TNFSF ligand, so the TNFRSF ectodomains: DcR3, TNFR2, DR3 |
| Foldseek structural neighbours of a **known** partner | not used here — DR3's ectodomain has no structure to search from, which is the whole reason the case is interesting |

**Excluded, and why:**

- the rest of the proteome — an unstated candidate set makes a false-positive
  rate meaningless, because you cannot compute a rate over a denominator you
  did not write down;
- intracellular partners of a secreted ligand, **except** as a negative control
  (DR3's own death domain is in the panel for exactly this);
- anything whose construct boundaries could not be read off a deposited entry
  or a UniProt topological domain. A guessed construct and a crystallised one
  are not comparable, and the difference between them will be attributed to
  biology.

Constructs are in `constructs.py`, one line of provenance each. Positive
controls use the range that was actually crystallised, read off
`rcsb_polymer_entity_align`. DR3 uses its UniProt extracellular topological
domain, 25-199, because there is nothing to read a range from.

### 2. Cofold, with the controls in the same run

```python
from predict import cofold_complex          # cofold-check, unmodified
seqs, chains_a, chains_b = chains_for("TL1A", "DR3")
r = cofold_complex(seqs, n_seeds=3, seed=42)     # block 1
r2 = cofold_complex(seqs, n_seeds=3, seed=142)   # block 2, disjoint
```

Two blocks, not one run of six. A single base seed advances `base_seed + i`
internally, and a contiguous block can be lucky as a block — measured, below.

**Every ligand needs its own scrambled-partner run.** Shuffle the partner
construct's residues (same length, same composition, no fold) and cofold it
against the same target. That run is the floor. Without it the ipTM number has
no scale: the scrambled floor measured **0.479 for IL-17A, 0.787 for TNF-alpha
and 0.785 for TL1A** — a 0.31 spread across three ligands in one panel. There
is no transferable absolute ipTM threshold and this skill does not define one.

### 3. Measure the interface yourself

`cofold_complex` returns `interface` computed on **seed 0 only**. Use
`interface_metrics(cif, chains_a, chains_b)` per seed instead. It returns the
CA-CA **pair** count under 8 Å (the same quantity as the 97 in the 8DYG
reference, not the 29-residue count), heavy-atom contacts under 5 Å, both
interface residue sets, and buried surface area from a Shrake-Rupley
implementation verified against a rigid-body control (translate one chain 500 Å:
BSA returns exactly 0.0, contacts 0).

Chain groups, not chains. A TNF trimer contacting one receptor is `["A","B","C"]`
against `["D"]`; passing one protomer measures a third of the epitope, which is
dossier rule 2b one level up.

### 4. The gate

Five conjunctive checks, in `ppi_hypothesis.gate()`. Conjunctive for the same
reason `graph-intake`'s ask trigger is: each one alone is satisfied by large
numbers of pairs that do not interact.

| # | check | threshold | what it actually tests |
| --- | --- | --- | --- |
| 1 | `enough_seeds` | >= 6 | 3 seeds lied — see failure mode 4 |
| 2 | `disjoint_seed_blocks` | >= 2 | same |
| 3 | `seed_concordance_mean` | mean pairwise Jaccard >= 0.70 | the model is **consistent** about where the partner goes |
| 4 | `seed_concordance_min` | min pairwise Jaccard >= 0.60 | one disagreeing pair is not averaged away |
| 5 | `scramble_control_margin` | median ipTM >= scramble floor + 0.08 | the model is more confident about this partner than about a random string of the same composition |
| 6 | `rank_in_candidate_set` | 1, or 2 behind a known partner | it prefers this partner to the others we offered |

Then two novelty tests, which are not scores:

- **deposited?** If a complex of the pair is in the PDB, the cofold is recall of
  training data. Report it as a method check. Never as a hypothesis.
- **already in the graph?** If the interaction is already a `links[]` row, only
  the *interface* is new, and an interface is a structural claim — see the
  never-fire list.

**Read what is not in that table.** Contact count, buried surface area, ipTM as
an absolute number, footprint transfer against a homologous complex, and
agreement with a second model are all measured, all reported, and none is
gated on. Each has a case in the panel where it points the wrong way.

### 5. The ask, if there is one

Use the graph's own four verbs — `expand_node`, `resolve_link`, `test_gap`,
`new_question`. `build_ask()` raises on anything else, because a fifth verb is
unconsumable upstream and the thing you want to ask with it is almost certainly
on the never-fire list.

| situation | verb |
| --- | --- |
| a predicted PPI the graph has no row for | `new_question`, `target: null` |
| a predicted PPI landing on a gap the graph already names | `test_gap` on the gap id |
| **our structure contradicts a row the graph carries** | `resolve_link` on that link id, carrying our answer |

The third is the one that fired on this run. It is the post-resolution ask
`graph-intake` already allows past gate 3 — it states our answer rather than
asking for theirs, and exists so a wrong row does not propagate.

Validate before sending:

```bash
python3 ../graph-intake/graph_read.py <graph.json> --check-ask '<ask json>'
```

Five mechanical gates. It does **not** check the three judgment gates and says
so; an all-green result is permission to *consider* an ask.

### 6. Where a generated hypothesis lives in the dossier

**Follow the existing convention exactly.** One entry in `not_found[]`, reason
prefixed `ASK[<verb>:<target id>]`. No template change, and none should be made.

```json
{
  "field": "tractability.pocket_vs_interface.partner_pdb_id",
  "reason": "ASK[resolve_link:L4] issued to graph g_tl1a1 round 3. Not blocking: the predicted interface is reported in tractability.pocket_vs_interface from our own measurement and this entry records the residual literature question. <the question text>"
}
```

The predicted interface itself is **not** a `not_found` entry. It is a
measurement and it goes where measurements go:

- `tractability.pocket_vs_interface.interface_residues` — the consensus epitope;
- `tractability.pocket_vs_interface.classification` — from `interface_analysis.classify_pocket`;
- `tractability.pocket_vs_interface.partner_pdb_id` — **null** when the partner
  complex is predicted rather than deposited, with the reason in `not_found`.
  A predicted complex is not a `partner_pdb_id`, and writing one there is how a
  prediction becomes indistinguishable from a crystal structure two steps later;
- `structure.cofold_control` and `leakage_risk: true`, always, because Boltz-2
  trained on the PDB and every TNFSF/TNFRSF prediction is downstream of a
  memorised complex.

**The rule is unchanged and absolute: an outstanding ask never blocks a verdict
and never licenses a null.** Complete the dossier as if the ask will never be
answered. If the only reason a field is null is that an ask is pending, the ask
is illegitimate.

## Where the gate came from

Fifteen cases, 21 cofold runs, 63 Boltz-2 diffusion jobs. Full rows in
`fixtures/panel_results.json`. `ipTM` is the median over all seeds; `conc` is
the mean pairwise Jaccard of the target-side interface residue set across seeds,
`cmin` the minimum pairwise; `CApr` the median CA-CA pair count under 8 Å;
`BSA` the median total buried area; `fpCov` the fraction of the deposited
homologous epitope the prediction covers.

| case | pair | expected | seeds | ipTM | conc | cmin | CApr | BSA | fpCov | rank | gate |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| P1 | IL17A/IL17RA | positive (7UWM) | 6 | 0.880 | 0.837 | 0.794 | 94 | 4835 | 0.918 | 1 | **PASS** |
| P2 | TNFA/TNFR2 | positive (3ALQ) | 6 | 0.937 | 0.771 | 0.694 | 40 | 2432 | 0.828 | 1 | **PASS** |
| P3 | TL1A/DCR3 | positive (3K51) | 6 | 0.902 | 0.875 | 0.771 | 38 | 2250 | 0.929 | 1 | **PASS** |
| N1 | TNFA/DCR3 | negative, same superfamily | 3 | 0.712 | 0.367 | 0.162 | 19 | 1597 | 0.406 | 3 | fail |
| N2 | TL1A/TNFR2 | negative, same superfamily | 6 | 0.838 | 0.536 | 0.128 | 47 | 2716 | 0.750 | 3 | fail |
| N3 | IL17A/TNFR2 | negative, cross-family | 3 | 0.527 | 0.264 | 0.143 | 57 | 3215 | 0.400 | 2 | fail |
| N4 | TNFA/IL17RA | negative, cross-family | 3 | 0.580 | 0.285 | 0.196 | 49 | 1969 | 0.531 | 4 | fail |
| N5 | TL1A/IL17RA | negative, cross-family | 3 | 0.673 | 0.307 | 0.180 | 76 | 4160 | 0.464 | 6 | fail |
| N6 | IL17A/CA2 | negative, cross-compartment | 3 | 0.411 | 0.548 | 0.522 | 39 | 1277 | 0.255 | 4 | fail |
| N7 | TL1A/CA2 | negative, cross-compartment | 3 | 0.596 | 0.234 | 0.133 | 52 | 2448 | 0.464 | 7 | fail |
| N8 | TL1A/DR3-death-domain | negative, wrong domain | 3 | 0.800 | 0.514 | 0.375 | 17 | 1493 | 0.429 | 4 | fail |
| N9 | TL1A/DR3-scrambled | negative, floor | 6 | 0.785 | 0.577 | 0.427 | **166** | **6569** | 0.840 | 5 | fail |
| N10 | IL17A/IL17RA-scrambled | negative, floor | 3 | 0.479 | 0.337 | 0.256 | 108 | 4199 | 0.473 | 3 | fail |
| N11 | TNFA/TNFR2-scrambled | negative, floor | 3 | 0.787 | 0.581 | 0.500 | 143 | 5333 | 0.688 | **2** | fail |
| U1 | **TL1A/DR3** | **unknown** | 6 | 0.886 | 0.789 | 0.708 | 48 | 2826 | 0.857 | 2 | **PASS** |

**Positives recovered: 3/3. Negatives passing: 0/11.** By the rule of three,
0/11 puts a 95% upper bound of **27%** on the false-positive rate — and the
subset that matters, the six hard negatives (same superfamily, wrong domain,
scrambled), is 0/6, upper bound **50%**. That is not a measured false-positive
rate. It is a rate consistent with anything up to one in two.

**Asks issued: 0/15.** The three positives are deposited (recall, not
prediction). The unknown is already a link in the graph. Nothing else passed.

### Reference interfaces, measured on the deposited assemblies

Everything above is calibrated against these, computed by the same code:

| entry | chains | CA pairs (8 Å) | BSA total | interface residues |
| --- | --- | --- | --- | --- |
| 8DYG | A/B | **97** | 3513 | 38 + 45 |
| 7UWM | AB/C | 72 | 4379 | 55 + 57 |
| 3ALQ | ABC/R | 42 | 2470 | 32 + 33 |
| 3K51 | A,A-2,A-3 / B | 30 | 2049 | 28 + 29 |

8DYG returns exactly 97, reproducing the documented convention, which is how we
know the contact code is measuring the same thing the rest of the station is.

### What the run proposes, on the case where we do not know the answer

TL1A/DR3 passed every gate. The consensus epitope on TL1A — the residues in the
predicted interface in **all six** seeds, UniProt O95150 numbering:

```
100 101 102 103 112 113 114 120 123 124 156 157 170 171 173 174 175 176
187 188 189 190 207 209 221 231 234 235 238 239 241
```

31 residues, against the 28-residue DcR3 epitope on TL1A measured from 3K51;
20 shared, Jaccard 0.513, covering **71%** of the DcR3 epitope. Note this is the
coverage of the six-seed **intersection**; the per-seed median in the table above
is 0.857, because a single seed predicts a larger epitope than the six seeds
agree on. Quote the one you mean.

On the DR3 side the predicted contacts map onto the DcR3 residues that contact
TL1A at coverage 0.82-0.89 across seeds, through a 26%-identity global
alignment — and TL1A/TNFR2, which does not bind, gives 0.83-0.90 on the same
measurement (failure mode 3).

**And it is not a hypothesis anyone needed.** The graph already carries
`L4: TL1A binds DR3` at `basis: primary`. Every immunology review says the same.
What the prediction adds is the interface — a structural claim, which the
never-ask list assigns to us.

## What must never become an ask

Longer than the trigger on purpose, and inherited from `graph-intake`'s list
rather than replacing it. Everything on that list still applies. These are the
additions this instrument makes possible.

**1. An interface, for an interaction the graph already has.** This is the
temptation the whole skill creates, and it fired on our one clean case. "The
graph says X binds Y; we predict *where*; ask them for literature evidence of
the contact" sounds like a literature question. It is not. Never-ask item 5 —
"anything structural ... interface classification" — is exactly this, and the
fact that we phrased it as a question about literature does not move it. If
mutagenesis or HDX data exist for a contact, we retrieve them ourselves.

**2. A pair whose complex, or a close homolog's complex, is deposited.**
Boltz-2 trained on the PDB. TL1A/DcR3 (3K51, 3MI8) is deposited and DR3 is a
26%-identity paralogue of DcR3, so our TL1A/DR3 prediction is downstream of a
memorised complex and its agreement with the DcR3 epitope is the expected result
whether or not DR3 binds there. Asking upstream to confirm something we recalled
from their own field's structure database is not a hypothesis, it is a lookup
with extra GPU.

**3. A pair the PDB settles.** If the loop nominates a partner absent from the
graph but present in the PDB — TL1A/DcR3 would have been exactly this, had the
fixture graph not happened to omit DcR3 — the answer is a `paperclip`/RCSB
query, not a round. Write the PDB id into the dossier and move on. This is
never-ask item 12 restated: an unfinished lookup is not an ask.

**4. A pair the candidate-selection step already implies.** Every partner this
loop proposed for TL1A came from the TNFRSF family list. Asking the graph
"does TL1A bind this TNFRSF member" returns the family logic we used to build
the candidate set. The ask has to carry information the candidate set did not.

**5. Anything gated on confidence, contact count or buried area.** Not a
literature rule, a measurement rule, and it belongs here because it is what an
agent will reach for. See failure modes 1 and 2.

**6. A prediction with no scrambled control in the same run.** The floor moved
0.31 in ipTM between three ligands in one panel. Without an in-run floor the
number has no scale, and "ipTM 0.85" is not a finding.

**7. A prediction on fewer than six seeds in fewer than two blocks.** Failure
mode 4 is the measurement that sets this, and it inverted a call.

**8. Anything you would not sign.** The ask goes to a team that cannot see our
run. A confident-looking `new_question` with a residue list in it reads as a
result. If it is wrong, it will be in their graph, cited, on the next round —
and our own `graph-intake` will read it back as a finding.

## Failure modes

The longest section, on purpose, because everything above is the easy half and
because five of these were found by controls that a plausible version of this
skill would not have run.

### 1. Contact count and buried area are not weak signals — they are inverted

The obvious gate is "a real interface buries a lot of surface". Measured, on
this panel, the opposite is true.

A **composition-scrambled random sequence** with DR3's amino-acid composition,
cofolded against the TL1A trimer, produced **166 CA-CA pairs and 6569 Å² buried**
— the largest interface in the whole panel. The deposited TL1A/DcR3 complex, run
through the identical pipeline, produced **38 pairs and 2250 Å²**. Every one of
the three scrambled controls landed in the top four by buried area:

| | CA pairs | BSA (Å²) |
| --- | --- | --- |
| TL1A/DR3-scrambled | 166 | 6569 |
| TNFA/TNFR2-scrambled | 143 | 5333 |
| IL17A/IL17RA-scrambled | 108 | 4199 |
| TL1A/DcR3 **(real, deposited)** | 38 | 2250 |
| TNFA/TNFR2 **(real, deposited)** | 40 | 2432 |

The mechanism is not mysterious: a sequence with no fold has nothing else to do,
so the model drapes it over the target. A real receptor ectodomain is a rigid
beta-rich stack that touches its ligand along one edge. **The more the partner
looks like a disordered string, the better it scores on every size-of-interface
metric.**

Reference values from the deposited assemblies say the same thing more quietly:
3K51's real TL1A/DcR3 interface is 30 CA pairs and 2049 Å². Any threshold set to
admit a real TNFSF/TNFRSF interface admits everything.

So `ca_pairs_8a` and `bsa_total_a2` are computed on every run, reported on every
run, and gated on never. They are in `REPORTED_NOT_GATED` in the module with the
numbers attached, because the next person to read this will want to gate on them.

### 2. ipTM does not separate at the seed level, and has no transferable scale

Two separate problems, and they compound.

**Seed level.** Across 18 positive-control seeds and 39 negative-control seeds,
the lowest positive single-seed ipTM was **0.8623** (IL-17A/IL-17RA, a deposited
complex) and the highest negative single-seed ipTM was **0.8693** — from the
**scrambled** DR3 string. There is no absolute ipTM threshold that admits every
deposited positive and excludes a random sequence. The medians separate; the
individual measurements do not, and a run reporting one seed is reporting a
number from the overlapping region.

**Scale.** The scrambled-partner floor is ligand-dependent and the dependence is
large: **0.479 for IL-17A, 0.787 for TNF-alpha, 0.785 for TL1A**. An ipTM of
0.79 is the floor for TNF-alpha and a spectacular result for IL-17A. Any fixed
threshold is a statement about which ligands happened to be in the calibration
set.

This is the same shape as the KRAS sealed-pocket result the station already
carries — confidence is a statement about the geometry the model drew, not about
whether that geometry should exist. Here it is worse, because the number also
has no fixed origin.

The gate therefore uses **margin over an in-run scrambled floor**, never an
absolute value. That margin is `0.401 / 0.150 / 0.117` for the three positives
and `0.101` for TL1A/DR3, against `0.052` for the hard negative TL1A/TNFR2. The
threshold of 0.08 sits in a gap of 0.049 between those two groups and was fitted
on exactly these 15 points. It is a line drawn after seeing the data.

### 3. Footprint transfer tests the site, not the partner — and it looks like it tests the partner

This is the signal that most deserved to work, and the one whose failure is
easiest to miss, because it is the only candidate signal that is not the model's
own opinion.

The idea: superpose nothing, just compare residue sets. Take the predicted
interface on the target, convert to UniProt numbering, and compare against the
interface measured on a *deposited* complex of the same target with a
*homologous* partner. If the model puts the new partner where the family's
receptors actually bind, that is external corroboration.

It does exactly that, and it corroborates the wrong things:

| pair | footprint coverage vs 3K51 TL1A/DcR3 epitope |
| --- | --- |
| TL1A/DcR3 (deposited, self-recall — not comparable) | 0.929 |
| **TL1A/DR3 (the true partner, unknown)** | **0.857** |
| **TL1A/DR3-scrambled (a random string)** | **0.840** |
| **TL1A/TNFR2 (a real receptor that does not bind TL1A)** | **0.750** |
| TL1A/IL17RA (cross-family) | 0.464 |
| TL1A/CA2 (cytosolic enzyme) | 0.464 |

The scrambled string covers 84% of the true epitope — because a string draped
over the whole trimer covers everything, and coverage has no penalty for
over-prediction. TNFR2 covers 75% because **every TNFRSF ectodomain binds its
ligands through the same CRD2-CRD3 face into the same groove between two
protomers**. The check confirms the family's canonical binding mode, which is
true of members that bind and members that do not.

The partner-side version fails identically. Mapping the predicted receptor
contacts onto DcR3 numbering through a global alignment: TL1A/DR3 gives Jaccard
0.535-0.571 and coverage 0.821-0.893 across seeds; TL1A/TNFR2 gives
**0.480-0.605 and 0.828-0.897** — indistinguishable, and TNFR2 is in fact the
*closer* paralogue of DcR3 (36.4% identity against DR3's 26.0%).

So footprint transfer answers "did the model put it in the family's site?" and
never "is this the right partner?". It is reported, with the negatives' values
beside it so the reader can see the ceiling, and it is not in the gate.

### 4. Three seeds lied, in both directions, on the case the conclusion turns on

The first pass ran three seeds per case, from one base seed. On those three
seeds the hard negative TL1A/TNFR2 scored **mean pairwise Jaccard 0.797 and
minimum 0.725** — squarely inside the positive band (0.746-0.886). It would have
passed a concordance gate.

Re-run on a disjoint block (seeds 142-144), the same pair gave **0.335 / 0.146**.
Pooled over all six seeds: **0.536 / 0.128**. The three-seed result was a block
that happened to agree.

Movement of the medians between blocks, same pairs:

| pair | ipTM, seeds 42-44 | ipTM, seeds 142-144 |
| --- | --- | --- |
| TL1A/DcR3 (true, deposited) | 0.9005 | 0.9036 |
| TL1A/DR3 (true, unknown) | 0.8899 | 0.8745 |
| TL1A/TNFR2 (false) | 0.8404 | 0.8347 |
| TL1A/DR3-scrambled | 0.8652 | 0.7249 |

The scrambled control moved **0.14** between blocks. On the first block alone it
sat 0.025 below the true answer; pooled, 0.101 below. A single block of three
seeds is not a measurement of either quantity.

Hence gates 1 and 2: at least six seeds, in at least two disjoint blocks. And
note *why* this was found — the margin looked thin, so the comparison was
replicated. Nothing in the first pass flagged it.

### 5. ESMFold has no interface signal, so two-model agreement is not available

The obvious independent check is a second model. It is not available here, and
the reason is that ESMFold fails the **positives**.

On IL-17A/IL-17RA — a deposited complex with 72 CA pairs in the crystal —
ESMFold returned **4 CA pairs, 376 Å² buried, centre-of-mass separation 54 Å**,
pTM 0.368, average PAE 23.4. Across the panel its CA-pair count ranged 0-12 with
no relation to whether the pair binds (TL1A/CA2, a cytosolic enzyme, scored 0;
TL1A/DR3, the true partner, scored 9; IL-17A/CA2 scored 12, the highest).

This reproduces the station's existing IL-17A observation (1 contact against 97
deposited) on a different complex. ESMFold's own confidence does flag its
failure, which makes it a usable *self*-filter — and useless as a corroborator,
because agreement with a model that always says no carries no information. The
founder's framing that "the two models disagree sharply, which is itself a
usable signal" does not survive: they disagree on the positives too.

ESMFold is still run on every case and reported, so that this claim stays
falsifiable if the model changes.

### 6. Predicted interfaces are bigger than crystallographic ones, uniformly

Every positive control over-predicted against its own deposited reference:

| pair | predicted CA pairs | deposited | predicted BSA | deposited |
| --- | --- | --- | --- | --- |
| IL17A/IL17RA | 94 | 72 | 4835 | 4379 |
| TNFA/TNFR2 | 40 | 42 | 2432 | 2470 |
| TL1A/DcR3 | 38 | 30 | 2250 | 2049 |

TNF-alpha is close; the other two are 10-30% high. Predicted models pack more
tightly than crystal structures and have no disordered loops or partial
occupancy. **Do not report a predicted BSA against a literature BSA as if they
were the same measurement**, and do not read a predicted-versus-deposited
difference as a real difference in interface size. If you must compare, compare
predicted against predicted.

This also disposes of the fixture graph's `1840 Å²` claim as a discriminator: our
deposited measurement is 2049 Å² total / 1024 per side and our prediction is
2826 Å², so the asserted figure matches neither, and the mismatch is a reason to
ask about the *provenance* of the figure rather than to conclude anything about
the structure.

### 7. Homo-oligomers make chain identity a trap in both directions

The ligands here are a homodimer and two homotrimers, so the same epitope exists
on two or three chains and different seeds place the receptor on different
protomers. Comparing interface residue **labels** across seeds would report
0 agreement for what is the same interface.

So `seed_concordance` drops chain identity and compares residue numbers only.
That is the correct choice here and it is the exact operation dossier rule 4
warns against for pocket matching on a C3-symmetric site — where discarding
chain identity makes the site unresolvable in principle.

The difference: there, the question was "is this the same pocket in two
structures", and residue numbers triplicate. Here the question is "did the
partner land on the same surface", and the surface is symmetry-equivalent by
construction. **If you reuse `_seqid_set` for anything other than cross-seed
comparison of one complex, you are reproducing the withdrawn 650-fold error.**

### 8. The candidate set is doing more work than the gate

TL1A/TNFR2 fails the rank check because DcR3 and DR3 are in the candidate set
and score above it. Remove them and TNFR2 ranks first among what remains, and
the rank check passes. The gate's discrimination on the hardest negative is
partly a property of having included the right answer.

That is not a bug that can be fixed inside the gate. It is why step 1 requires
the candidate set to be written down, and why the panel includes an
in-set scrambled control for every ligand — the scramble is the one candidate
you can always add, and it is the only one guaranteed not to be a true partner.

Note also what the ranks say on their own: **a scrambled version of TNF-alpha's
own receptor ranked second in TNF-alpha's candidate set, above the real DcR3
ectodomain**. Ranking is a usable signal within a target — it recovered the
correct partner at rank 1 for all three ligands — and it is not a clean one.

### 9. Positive controls are training data, so the positive band is a recall band

All three positives are deposited: 7UWM, 3ALQ, 3K51. Boltz-2 trained on the PDB.
Their ipTM of 0.88-0.94 and concordance of 0.77-0.88 measure how well the model
reproduces structures it has seen, and the gate's thresholds are set from that
band.

There is no held-out positive in this panel and none was available — every
TNFSF/TNFRSF and IL-17 family complex we could use is older than the training
cutoff. So the claim "3/3 positives recovered" is a claim about recall, and its
transfer to a genuinely novel complex is **unmeasured**, not merely uncertain.

This is the single largest hole in the calibration and no amount of extra
negatives fills it. What would fill it is in "The hard question".

### 10. A confident interface on the wrong domain

TL1A against DR3's intracellular **death domain** — the wrong side of the
membrane, and the only part of DR3 that is deposited — scored ipTM 0.800 and
mean concordance 0.514. It failed the gate, comfortably. But 0.800 is above the
TNF-alpha scrambled floor (0.787) and would read as a strong result to anyone
quoting an absolute number, and the death domain is exactly the construct
somebody would pick if they took "DR3" from a database without reading the
topology.

Read the construct provenance before reading the score. `constructs.py` carries
one line per entry for this reason.

### 11. The ask becomes the thing you do instead of the structural work

`graph-intake`'s failure mode 11, restated for this instrument, because the
shape here is worse. A dossier carrying a predicted interface *and* an ask about
it looks like a station that did the work and then went further. It can be a
station that produced a picture and outsourced the question of whether the
picture is right.

The test: **would you have written this ask if the cofold had not run?** If the
answer is yes — if the ask is really "does TL1A bind DR3", or "is there
structural information about this pair" — the cofold was decoration and the ask
was always a literature question you could have asked without spending GPU. If
the answer is no, the ask depends on a specific measured claim, and it should
name that claim, its controls and its floor in the `question` text.

### 12. `new_question` cannot pass gate 4 as `graph_read.py` implements it

Mechanical, and it blocks the intended verb for this whole skill.
`already_asked()` matches on `(verb, target)`. `new_question` always carries
`target: null`. So once **any** open question has been asked against a graph,
`NOT_ALREADY_ASKED` fails for every subsequent `new_question`, whatever it is
about. Verified against `upstream_graph_askback.json`, whose `rounds[0]` is a
`new_question` at `null`: a completely unrelated structural `new_question`
returns `FAIL NOT_ALREADY_ASKED` and exit 1.

We did not change `graph_read.py` — it is not this skill's directory. The patch
is proposed at the end of this file. Until it lands, a genuinely new PPI has no
route through the checker, and routing it through `expand_node` or `test_gap`
instead would be pointing the ask at a row that is not what it is about.

## Cost

Modal GPU, warm containers, this panel:

| | |
| --- | --- |
| cofold runs | 21 (15 cases, 6 of them replicated on a disjoint seed block) |
| Boltz-2 diffusion jobs | 63 (3 seeds per run) |
| ESMFold runs | 11 |
| total wall clock | 69.4 min (3503 s cofold + 659 s ESMFold) |
| estimated spend | **$3.21 - $5.27** (A100-40GB $2.78/h to H100 $4.56/h) |

Scaling notes: a 534-residue two-chain complex took 141 s for 3 seeds; an
812-residue four-chain complex took 651 s. Cold starts dominate the first call.
`n_seeds` multiplies GPU time roughly linearly, so the six-seed requirement
doubles the cost of every case — and failure mode 4 is what that buys.

A candidate set of *n* partners for one ligand costs *n+1* runs (the +1 is the
mandatory scramble), each at 6 seeds. For TL1A's seven candidates that is
roughly 25 minutes of GPU, about $1.50. **The controls are most of the cost and
all of the value.**

## The hard question: does this make the pipeline better?

Asked directly, and answered directly.

**On the evidence measured here, a predicted protein-protein interface should
not be sent upstream as a new-PPI hypothesis.** The loop should be run — it
produces a real, reportable measurement for the dossier — but step 4 should not
fire in the forward direction.

Four measured reasons:

1. **The signals that would justify it either fail or are unavailable.** Contact
   count and buried area are inverted (failure mode 1). Absolute ipTM does not
   separate at the seed level and has no transferable scale (2). Footprint
   transfer tests the family's site, not the pair (3). Two-model agreement is
   not available because the second model fails the positives (5).

2. **What survives is a reproducibility measure, not a truth measure.** Seed
   concordance over six seeds, and margin over an in-run scrambled floor. Both
   say the model is consistent and more confident than it is about noise.
   Neither is evidence the interaction happens. And the false-positive rate is
   0/11 with a 27% upper bound — 0/6 and a 50% upper bound on the hard
   negatives, which is not a measurement.

3. **The positive band is a recall band** (9). Every threshold above was set on
   complexes the model was trained on. Applying them to a novel pair is the
   extrapolation the whole calibration cannot support.

4. **On the one case where it fired cleanly, the answer was already in the
   graph.** TL1A/DR3 passed every gate — and `L4` already says TL1A binds DR3,
   at `basis: primary`. What the prediction added was the interface, which
   never-ask item 5 assigns to us. Every other partner the family logic
   nominated was either deposited (DcR3 — settled by an RCSB query in
   milliseconds) or a negative.

   This is the same shape the ask-back agent found: three of its four motivating
   cases were answerable inside our own pipeline in under half a second. Here it
   is four of four, and the instrument that answers them is the PDB.

**So the loop runs in the other direction, and that direction does work.** On
this run the structural work **falsified a claim the graph carries**: `L4` cites
a 2.4 Å crystal structure of TL1A with three DR3 protomers burying 1840 Å². A
PDB census against Q93038 returns three entities, all the intracellular death
domain, no ectodomain structure of any kind, no complex. That is a wrong row
that will propagate to every other consumer of the graph, we settled it
ourselves, and telling them is the post-resolution `resolve_link`
`graph-intake` already allows. It carries our answer rather than a request, it
cannot pollute the graph with speculation because it is a database census, and
it costs one round. The worked ask is in `fixtures/worked_ask.json` and passes
all five mechanical gates of `--check-ask`.

That is a closed active-learning loop. It is just not the direction it was
imagined running: **graph claim -> our structure -> correction back to the
graph**, rather than **our structure -> new claim -> graph confirms**.

### What would have to be true to change the answer

Stated so this can be revisited rather than re-argued:

1. **A held-out positive set.** Ten or more complexes deposited after the
   model's training cutoff, run through the identical pipeline. Until that
   exists, every threshold here is fitted on memorised structures and failure
   mode 9 stands. This is the blocking item.

2. **Hard negatives at n >= 20**, all within the target's own superfamily. Six
   gives a 50% upper bound on the FPR; twenty would give ~15%; sixty, ~5%. The
   negatives are cheap — a scrambled partner is free to construct — so this is
   a budget question, not a research question.

3. **A signal that is neither the model's own opinion nor the family's
   canonical geometry.** The two candidates worth trying, neither available in
   this deployment: paired-MSA coupling computed for the *specific* pair, and
   cross-species conservation of the predicted contact residues on both sides
   simultaneously. Both are independent of the structure model. A predicted
   interface whose contact residues co-vary across species on both sides would
   be evidence of a kind nothing above provides.

4. **A candidate that the family logic does not already imply.** Every partner
   this loop proposed for TL1A was proposed *because* it is a TNFRSF member, so
   the prediction adds nothing the candidate-selection step did not assume. A
   Foldseek-neighbour candidate from outside the family, passing the gate, would
   be a genuinely new claim — and would then need points 1-3 behind it.

If 1 and 2 land and the FPR holds under 10% on held-out positives, the forward
ask becomes defensible for pairs that clear points 3 and 4. Not before.

## Proposed patches, not applied

This skill owns only `.claude/skills/ppi-hypothesis/`. Two changes elsewhere
would help and were not made:

**`graph-intake/graph_read.py`, `already_asked()`.** Match `new_question` on a
stable hash of the `question` text when `target` is null, not on
`(verb, None)`. As written, one prior open question retires the verb for the
life of the graph (failure mode 12).

**`graph-intake/SKILL.md`, the post-resolution ask.** It is documented as "the
only ask allowed past gate 3" and is silent on gate 2. Our contradiction ask
targets a link whose `basis` is `primary`, so it fails gate 2 as written — yet
gate 2's own rationale ("a primary-supported claim is either usable or
contradicted by something we can measure") names this case and then routes it
nowhere. A wrong `primary` row is more damaging than a wrong
`background_only` one, not less. Propose making the post-resolution
contradiction ask exempt from gate 2 as well, with the same "state our answer
and its source" requirement.

**`managed/druggability-dossier/CLAUDE.md`.** No change proposed. Rule 2b's
`pocket_vs_interface` block already has the right fields for a predicted
interface, and `structure.cofold_control` plus `leakage_risk` already carry the
warnings. The convention in rule 13 — null the field, record the reason — is the
correct home for `partner_pdb_id` when the partner complex is predicted.

## What this skill does not do

- It does not decide that a predicted interface is real.
- It does not gate on confidence, contact count or buried surface area.
- It does not find pockets, score druggability, or rank targets.
- It does not invent a fifth ask verb.
- It does not write to the upstream graph, `rounds` included.
- It does not let a pending ask block a verdict or license a null.
- It does not treat a cofold of a deposited complex as a prediction.
