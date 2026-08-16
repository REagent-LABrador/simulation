---
name: pocket-scan
description: >
  Detects and measures ligand-binding pockets across an ensemble of structures,
  sweeping fpocket clustering, and quantifies whether a site is cryptic and by
  which mechanism (backbone collapse vs steric occlusion). Reports volume with
  its spread as an absolute quantity comparable across structures, and
  druggability ONLY as a within-structure rank among that structure's pockets
  with the pocket count beside it. It does NOT decide whether a target is
  druggable, does NOT rank targets, does NOT compare a druggability value across
  structures or against a threshold, and does NOT interpret a low score as
  evidence against tractability.
---

# pocket-scan

Pocket geometry over an ensemble, with the method's blind spot measured rather
than assumed.

Two calibrations drive everything below. Both were run in-repo; the numbers are
ours, not literature.

## Setup

Runs inside the Modal CPU image, which carries `fpocket` (conda-forge, 4.2.3),
**P2Rank 2.5.1** (MIT, needs JDK 17 — Java 11 dies with
`UnsupportedClassVersionError`, class file v61), and `proto-tools`. The fpocket
binary self-reports `fpocket 4.0` in its banner — cosmetic upstream mismatch,
not a wrong install.

**`mdpocket` is already installed** — it ships inside the same conda package as
`fpocket`, so there is nothing to add. It is not an optional extra here; it is
how a site gets fixed across an ensemble (see below).

## fpocket detects, PRANK ranks

Keep these jobs separate in your head. fpocket's alpha-sphere *detection* is
sound; its *ranking* is the weak link, and its druggability score is a
three-descriptor logistic regression fitted on 21 positives.

The best-recall configuration in the LIGYSIS benchmark of 13 predictors is
**fpocket detection + PRANK rescoring** — 60% top-(N+2) recall, ahead of
DeepPocket at 58% and P2Rank standalone at 52%. Note that even the winner
recovers only 60% of known sites; there is no method here that finds everything.

Measured on isolated fixtures:

| site | fpocket rank | PRANK rank |
| --- | --- | --- |
| 6OIM switch-II (sotorasib) | 9 | **2** |
| 2AZ5 SPD304 | 2 | **1** |
| SPD304 site across 4 apo TNF-alpha trimers | druggability noise | rank 2-3 in all four |

**Read the third row narrowly.** "The SPD304 site" in those apo structures was
assigned by the residue-number matching heuristic that is now **withdrawn**
(failure modes below): mdpocket places that matched pocket 7.7 A from the real
SPD304 site. So that row says PRANK ranks *a* consistently-detected pocket
highly across apo trimers — it does not say PRANK found the SPD304 site. The two
holo rows are unaffected, because there the site is defined by the bound ligand
rather than by matching.

**Settled at n = 70 ligand-anchored measurements across 8 targets: PRANK
promotes 79% and demotes 1%.** Median rank **5 → 1**; top-3 recall **37% → 91%**.
The single demotion is the known KRAS case below. Earlier versions of this file
said rescoring "has not yet helped, and once it hurt" — **that is falsified** and
is void.

Supporting per-target runs, in the order they were measured:

- **IL-17A** promoted the ligand site in all three structures: fpocket rank 6 to
  PRANK 2, rank 5 to PRANK 2, rank 11 to PRANK 1.
- **NLRP3**, the clearest single case: fpocket ranked the true site 3rd to 34th,
  **median 18**. PRANK put it at **rank 2 in 11 of 14 instances and never worse
  than 4.**

**But PRANK rank is a site finder, not a quality score, and used as the latter it
is inverted.** As a *cross-target druggability* classifier its AUC is **0.25** —
worse than chance, in the systematic direction — because on a target with no
ligand to anchor to, the top-ranked pocket is top-ranked by construction, so
"rank 1" carries no information about whether the site is good. **It finds sites.
It says nothing about their quality.** Never let a PRANK rank stand in for a
druggability judgement.

**Why this one worked when the cross-target evaluations did not.** Rescoring is a
**within-structure reordering** — the operation this whole family of quantities
supports (see "fpocket druggability is a within-structure quantity" below). So
the n=70 promotion result and the AUC 0.25 result are consistent rather than in
tension: the first was measured with a supported operation, the second was not.
**PRANK rank is a second within-structure ordering, legitimate on exactly the
same footing as fpocket's own rank** — not better, not a tiebreaker, and not a
substitute. Report both, and report a disagreement as a disagreement.

**The KRAS negative still stands and is not deleted.** A method that helps on
three targets and hurt on one is a more useful thing to know than a method that
always helps:

- Our fpocket invocation already ranks 6OIM's switch-II pocket **#1**, so there
  was nothing to promote. The rank-9 figure came from a different invocation
  producing 11 pockets against our 9. On that structure it has no work to do.
- At 6OIM D=1.6 PRANK **demotes** the true site, fpocket rank 1 to PRANK rank 3.

Read the two together: rescoring helps most where fpocket's own ranking has
buried the site (IL-17A ranks 5-11, NLRP3 median 18) and can hurt where fpocket
already has it at rank 1 (KRAS). It is not uniformly an improvement, and it is
not a tiebreaker. It is a second, independently trained opinion over the same
geometry — valuable because two methods disagreeing is information, not because
one is right.

`prank_rank` is reported **alongside** fpocket's rank, never replacing it. Read
a gap between them as a flag for manual attention, **not** as evidence that
PRANK found something fpocket missed — on our structures it points the other way
at least as often.

### Two P2Rank gotchas, both confirmed by direct test

**The `probability` column is only calibrated in `predict` mode, not `rescore`.**
In rescore mode the true SPD304 site scored 0.011 while a large decoy scored
0.783 — the ranking is usable, the probability is not. `predict` mode on the
same site gives 0.735. So use `rescore` for within-structure ranking and a
separate `predict` run if you need a cross-structure comparable score.

**`-chains A` is silently ignored by `predict`.** Passing it returned all 3,483
atoms of a trimer and an identical score. The documentation shows the flag; it
does nothing. Use the `chains` column of a dataset (`.ds`) file instead.

Also: `rescore` emits no `_residues.csv` — P2Rank only lays SAS points over the
surface in `predict` mode. And skip the `conservation_*` models, which need
HMMER and MSAs; `default` and `rescore_2024` use structure-derived features only.

## mdpocket: fix the site by construction, not by matching

Across an ensemble, the hard part is not detecting pockets — it is asserting
that two detected pockets are *the same site*. Matching them after the fact on
shared residue numbers is the single worst thing this skill has done, and the
failure modes below record exactly how it broke. mdpocket removes the step:
you define the site once, as a grid, and apply that one definition to every
superposed structure.

Two modes, both verified:

```bash
# mode 1 — exploration: pocket density over the whole set
micromamba run -n druggability mdpocket --pdb_list list.txt -o PREFIX

# mode 2 — characterization: one fixed site, measured per structure
micromamba run -n druggability mdpocket --pdb_list list.txt \
    --selected_pocket sel.pdb -o PREFIX
```

**Mode 1 (exploration)** emits `PREFIX_dens.dx` and `PREFIX_freq.dx` on a 1.0 A
grid, plus isosurface PDBs. Use it to find where a pocket is and how often it is
open.

**Mode 2 (characterization)** emits `PREFIX_descriptors.txt` with **one row per
snapshot and 41 columns**. This is the mode that makes an ensemble comparable,
because every row was measured inside the same grid.

Runtime is **0.85 s for five ~3.5k-atom trimers**. Compute is not a factor in
this decision — use it whenever there is more than one structure.

### What it bought, measured

Fixing the site by construction instead of by post-hoc matching cut the
across-ensemble CV on the TNF-alpha apo set:

| how the site was established | volume CV across the ensemble |
| --- | --- |
| post-hoc residue-number matching | ~28% (measured 28.1% at D=1.6) |
| **site fixed by construction (mode 2)** | **~10%** (measured 9.9%) |

The matching heuristic inflated the spread roughly **2.8-fold**, and essentially
all of that came from the single structure that matched a pocket 12 A away from
the others. The spread we were reporting was mostly the matcher, not the protein.

**Two significant figures, never three.** fpocket estimates pocket volume by
Monte Carlo and mdpocket inherits it. Three *identical* reruns of one 5-structure
ensemble gave volumes of 280.6 / 276.1 / 274.6 A^3 and CVs of 12.1 / 11.3 /
10.8%; the deployed run gave 9.9%. So **about 1 percentage point of any CV here
is fpocket's own volume estimator**, and a CV difference smaller than that is not
a difference between sites. An earlier version of this table read "27.8% to
10.2%" — the improvement is real, the third digit never was.

**This CV was measured on `site_from_density`, which is not the ligand site.**
See the next section. It measures reproducibility, not correctness.

### The two site definitions, and which one is the pocket

`pocket_scan` returns `mdpocket.sites` with up to two entries, and **fixing the
site by construction buys reproducibility, not correctness**. It guarantees every
structure was measured at the same grid points; it says nothing about whether
those points are the site anyone asked about.

| key | definition | is it the pocket? |
| --- | --- | --- |
| `site_from_ligand` | grid points within 3.0 A of the holo ligand, transferred by superposition | **yes**, by construction |
| `site_from_density` | largest connected cluster of grid points open in *every* structure | **not necessarily** — the most *persistent* cavity |

On the apo TNF-alpha ensemble the density site's centroid sits **7.73 A** from
the transferred SPD304 ligand. It is the on-axis cavity — real, well-formed,
reproducible, and **precisely the pocket the withdrawn residue-number matcher
reported as "the SPD304 site"** (see the failure mode below). Detecting it is not
the error. Calling it the ligand site is, and doing so reproduces the retracted
finding in a form that looks like a result.

So, before quoting any number off a site entry:

1. **Prefer `site_from_ligand` when it exists.**
2. **Read `distance_to_donor_ligand_centroid_a`.** Every entry carries it
   unconditionally, plus `ligand_anchored` and an `off_site_warning`. A number
   quoted without it is unverified.
3. **A proposed — not calibrated — threshold: 4 A.** Beyond that, treat the
   centroid as a *different pocket*. It is roughly half the single 7.73 A error
   we have measured and well above the ~1 A grid spacing, and it rests on one
   case. Label it a proposal wherever you use it. It gates a warning, never a
   refusal; no number is dropped because of it.
4. **A null distance is itself the finding.** A pure-apo ensemble with no
   transferable ligand can return `site_from_density` as the *only* site with
   `mdpocket_status: "ok"` — a confident single answer about a cavity of unknown
   identity. `distance_reason` says why. Carry it into the caveat.

### KRAS is the richer case, because the pocket only partially collapses

TNF-alpha gives a clean zero (below), which is honest but not very informative.
KRAS shows what mode 2 is actually for:

| mdpocket characterization | 6OIM holo (ligand stripped) | 4OBE apo |
| --- | --- | --- |
| volume | **1152.3 A^3** | **452.1 A^3** |
| alpha spheres | 500 | 182 |
| mean local hydrophobic density | **185.8** | **12.6** |

And mode 1 over the same two structures localises *which part* of the site goes
away: **179 grid points at frequency 1.0** — the nucleotide-adjacent shelf,
present in both — against **322 points at frequency 0.5** — the cryptic
switch-II sub-pocket, present only in the ligand-bound conformer. The site does
not vanish; a specific sub-pocket does, and mdpocket says which one.

(With N=2, "frequency 0.5" is presence-in-one-of-two, not a frequency. It is
being used here to *localise*, not to quantify how often the pocket is open.
See failure mode 3.)

## Procedure

### 1. Build an ensemble, not a structure

**One structure is not a measurement.** Query every PDB entry for the accession
(`pdb_v.structures_by_accession`), classify each as holo or apo using
`pdb_v.entry_ligands`, and take a set — all of them if few, otherwise the best
resolution across distinct crystal forms. Record which entries you used.

### 2. Prepare each

Protein only, altloc A or blank, hydrogens stripped (high-resolution entries
ship riding hydrogens that skew alpha-sphere volumes). Record missing residues.

**Chain selection is per-target, not a default.** KRAS is a monomer — one chain.
TNF-alpha's site sits on the trimer's 3-fold axis and **vanishes if you take one
chain**. Ask where the site is before deciding, and record the choice.

### 3. Sweep clustering — never pin it

```bash
fpocket -f <prepared>.pdb -D 1.6
fpocket -f <prepared>.pdb -D 2.4
```

Report both. See failure modes for why a single value is a coin flip.

### 4. Measure, in priority order

1. **Volume** at the site, with its spread across the ensemble. An **absolute**
   quantity — comparing it across structures is a legitimate operation.
2. **Druggability as a WITHIN-STRUCTURE RANK** — the site pocket's rank among
   that structure's pockets, with the pocket count and the PDB ID. "rank 1 of 30
   in 6OIM" is the claim. The value may sit beside the rank; the rank is what is
   asserted.
3. Hydrophobic density, lining residues with chain IDs.

**The gap between 1 and 2 is not that volume is more reproducible. It is that
they are different KINDS of quantity.** Volume is a number about a cavity.
Druggability is a number about a cavity *relative to the other cavities in the
same file* — see "fpocket druggability is a within-structure quantity" below,
which reads it out of the source.

**The one-protein proof, and the only demonstration you need:**

| RORgt entry | site MLHD | that structure's MLHD max | normalises to | druggability |
| --- | --- | --- | --- | --- |
| **4NB6** | 30.722 | 30.722 — *the site is the maximum* | **1.0** | **0.827** |
| **6C1P** | 19.0 | 52.767 | **0.36** | **0.009** |

Same protein, same orthosteric site, comparable absolute hydrophobic density.
**The 90-fold gap comes entirely from which other pockets happened to co-exist in
the file.** A value quoted without its structure's pocket population is not
interpretable, and a value compared against another structure's is not a
comparison.

Three consequences for this skill:

- **Never pool druggability across structures into a min/max.** A spread of
  druggability across an ensemble measures nothing. A range across the **D
  sweep within one structure** is legitimate and is a different object — say
  which one you built.
- **Never take a max over pockets.** r(n_pockets, max druggability) = **0.702**
  at D=1.6: max-over-pockets measures pocket count. That is what
  `max_druggability_no_ligand_site` computes, and it contaminated **70% of the
  hard class** of the retracted calibration set.
- **The clustering swing is the same mechanism, not a separate defect.** Median
  within-structure |D=2.4 − D=1.6| is **0.229**, maximum **0.955**. Changing D
  changes the pocket population, which changes the normalisation. Sweep and read
  both; do not average them.

Low values on holo structures with a drug physically bound are common and are
**not** findings about the pocket — JAK1's median is **0.009** across nine
approved drugs, TYK2 6NZP with deucravacitinib is **0.169**, BCL-2 6QGK is
**0.025**, NLRP3 runs **0.001–0.018** across seven holo crystals including one
carrying a clinical compound.

**Two named cases that used to sit in this list are struck.** *"EGFR 6LUD with
osimertinib bound scores 0.013"* is **off-site**: Jaccard **0.077** to the
osimertinib site, centroid spread **10.49 Å**, and at D=2.4 the pocket that
genuinely overlaps scores **0.174**. It was the same failure mode as the 651x
retraction, sitting inside the sentence used to justify the demotion — do not
quote it. *"RORgt 6C1P is 0.009 at rank 55 of 60"* is struck as a false-negative
case because **6C1P contains no RORgt** (sole entity A8EVM5, an ion transport
protein; anchor ligand `1N7` is CHAPSO). 6C1P survives above only as the
normalisation demonstration, where what it contains does not matter — what is
being shown is arithmetic on a pocket population.

**And the evaluation that produced the AUCs could not have established the
negative it was read as establishing.** Exact permutation over all 3,003 label
assignments: observed AUC 0.720 gives **p = 0.103**, and the minimum AUC this
design can call significant is **0.760**. See `falsification-sweep` check 10b.

### 5. Establish which pocket matters — apo is the normal case

**Assume no holo structure.** Most targets worth asking about are apo, cofolded,
or a bioemu ensemble; a bound drug-like ligand is the lucky exception, not the
baseline. A procedure that only works on holo structures is a validation harness,
not a working method.

Ranked by strength, use whichever routes are available and say which you used in
`tractability.site_hypothesis_basis`:

**(a) Holo ligand site — only when it exists.** Residues within 5 A of a
drug-like ligand; report Jaccard against the detected pocket. Strongest, rarely
available.

**(b) Persistence across the ensemble — the primary signal when nothing is
bound.** Run every conformer and rank by *how often a pocket appears*, not by its
best score. A pocket found in most conformers is credible; one appearing in a
single frame is noise. Persistence and volume are the reproducible quantities.

**Do not establish "the same pocket" by matching residue numbers across
independently detected pockets.** That step is withdrawn, and it took a headline
finding down with it — see the failure modes. On any homo-oligomer it cannot
work even in principle.

Instead: superpose the ensemble, define the site once, and push that one
definition through every structure with **mdpocket** characterization mode.
Exploration mode over the same set gives pocket *density*, which is the honest
form of "how often is it open". Both are covered above. Per-frame fpocket calls
leave you a pile of unaligned results to reconcile by hand, and the
reconciliation is the part that breaks.

If for some reason you must match pockets post hoc, **report the matched
centroid distance across the ensemble**, not an overlap fraction. Two pockets
sharing residue numbers can be 12 A apart and the overlap score will not tell
you.

**(c) Site transfer from a structural neighbour.** Foldseek the apo structure,
find a neighbour that *does* have a drug-like holo entry, superpose, and map its
ligand site onto your model. This manufactures a site hypothesis where the target
itself has none. Flag it clearly as transferred, with the source PDB ID and the
alignment quality — it is a hypothesis, not a measurement.

**(d) Curated annotation.** UniProt `Binding site` features give residue
positions for cofactor and substrate sites. Report the fraction recovered. Weak
on its own — it tells you where the *natural* ligand goes, which is often not
where a drug would.

If none of (a)-(d) apply, report the top pockets by persistence and say plainly
that no site hypothesis could be established. Do not silently promote the
highest-druggability pocket; on an apo structure that number is nearly
uninformative.

### 5b. Cofolded and bioemu ensembles carry an extra caveat

A pocket that appears only in predicted or sampled conformers, and in no
experimental structure, is a prediction about a pocket — two inferences deep.
Record it, mark `structure.tier` accordingly, and do not let it carry the same
weight as a site seen in a crystal.

For **cofolded** structures specifically: a cofold produced *with* a ligand has
had the pocket opened by the ligand you supplied. Finding a pocket there is
close to circular. Where a crystal structure also exists, score the cofold
against it first (`structure.cofold_control`) — if it cannot reproduce a known
structure for this target, its pockets are not evidence.

### 6. Measure cryptic risk — do not flag it from tier

Where both apo and holo exist, superpose on core C-alpha excluding the mobile
region, place the ligand in the apo frame, and compute:

- **max backbone C-alpha displacement at the site**
- **clash attribution** — backbone vs side-chain vs another subunit
- **free-volume fraction** of the ligand in the apo frame

Then classify the mechanism (table below). Where no holo exists, report the
ensemble backbone spread at the site instead and say the mechanism is
undetermined.

## Failure modes

### There is no correct fixed `-D`, and pinning one produces false negatives

`-D 1.6` was tuned on KRAS, where the default `-D 2.4` fuses the nucleotide and
switch-II sites into one 1540 A^3 mega-pocket scoring 0.886 — a cavity no
molecule occupies.

Applied to TNF-alpha, that same pin gives **druggability 0.002 at the site of a
co-crystallised 570 Da ligand in 2AZ5**. A false negative on a *holo* structure.
Diagnosis with `-i 5`: the channel fragments into alpha-sphere clusters of 15,
12 and 5, and in the apo the cluster sitting exactly on the ligand position has
12 spheres — below fpocket's default `-i 15` minimum — so it is **discarded
silently**. The same site at `-D 2.4` scores 0.346, rank 2 of 14, Jaccard 0.74.

Sweep. Report the range **within one structure across the D sweep** — that is
the range this section is about, and it is the only druggability range there is;
a range across structures is the type error above. A volume above ~1000 A^3 means
sites have merged; "not detected" at one D and present at another means
fragmentation, not absence.

### VOID: the 651-fold druggability spread across five apo TNF-alpha trimers

Earlier versions of this file reported, in two places, that five apo TNF-alpha
trimers held volume to +/-16% at "the same site" while druggability swung
**651-fold** (0.001 in 2ZJC to 0.651 in 1A8M, volumes 206.7–309.2 A^3; the same
claim appears as 650-fold in older copies — 651 is what 0.651/0.001 gives, and
it is the form used everywhere the retraction is now cited). **That figure is
WITHDRAWN — and so is the volume range printed beside it.**

**There were TWO defects, and the deeper one was named later.**

1. **The pooling.** Druggability is normalised inside each structure (see
   "fpocket druggability is a within-structure quantity" above), so **pooling it
   across five structures manufactures a meaningless spread on its own** — no
   matcher error required. **This alone was sufficient to produce 651x.** It is
   the same root cause as the whole within-structure rule, and it is why no
   improvement to pocket matching would have rescued the number.
2. **The matcher.** Separately and additionally, the pockets were matched across
   structures on shared residue *numbers*, chain-agnostically. That is a real
   defect, it is what withdrew the **volume** range printed beside the spread,
   and it is documented in full below.

Both retractions stand. The paragraphs that follow are the second defect.

mdpocket over the superposed ensemble showed what that matcher was actually
tracking:

- the matched pocket's centroid sits **7.7 A** from the SPD304 site it claimed
  to be measuring;
- it was not even self-consistent — **1TNF matched a pocket 12.2 A away** from
  where the other four matched. "The same site" spanned 12 A across five
  structures;
- the cause is structural, not a threshold to tune. A 19-residue reference on a
  homotrimer collapses to only **11 distinct residue numbers**, because the
  three protomers triplicate them. Throw away chain identity and a C3-symmetric
  site is **unresolvable in principle**.

The pocket it matched is real, just not the one claimed: an **on-axis cavity
lined symmetrically by Q61/K98/P117/I118/Y119 from all three chains, 107 A^3 at
frequency 1.0**. Well-formed, reproducible, and the wrong pocket.

Do not cite 650x, 651x, +/-16%, or 206.7–309.2 A^3. If you meet them in an older
dossier, they are void. What replaces them is below.

Note what does *not* change: **never build a verdict on a druggability score**,
and druggability remains a 3-descriptor regression fitted on 21 positives. That
claim never rested on the 651x figure — it rests on the KRAS holo/apo collapse
and on the provenance of the score itself.

Note also what this does **not** forbid: reporting the site pocket's **rank among
that structure's pockets, with the count**. That is the reportable form and it is
a within-structure statement. What is forbidden is a *value* travelling between
structures — as a spread, a threshold comparison, or a verdict.

### What replaces it: the spread was mostly the matcher

Fixing the site by construction rather than by post-hoc matching cut the
across-ensemble volume CV from **~28% to ~10%** (measured 28.1% at D=1.6 against
9.9%) — roughly a **2.8-fold** inflation, essentially all of it contributed by
the one structure that matched 12 A away. Quote two significant figures: ~1
percentage point of either number is fpocket's Monte-Carlo volume noise.

**A pocket-matching step is itself a measurement, and it needs its own
controls.** It was never treated as one, which is why a 12 A error survived to
become a headline number.

### At the true SPD304 site, the honest answer is that there is no pocket

mdpocket returns **0.00 A^3 in four of the five apo structures**. Not a low
score — nothing.

That is consistent with the physics rather than in tension with it. Place the
ligand into each *intact* apo trimer and you get **27–29 heavy-atom clashes
under 2.0 A**, minimum interatomic distance **0.28–0.53 A**, attributed
**identically in all five** to the third protomer (chain C: S60, Y119, L120,
G121) plus the Tyr119 triad. SPD304 does not bind this site as a trimer; it
binds after displacing a subunit. Delete that third chain and every apo
structure recovers the pocket immediately (~280–550 A^3 — see the
subunit-removed control below). Both measurements say the same thing: the site
is pre-formed, and a protomer is standing in it.

The part worth carrying to other targets is the *behaviour*, not the number:
**mdpocket returned 0.00 rather than silently substituting a nearby pocket.**
The residue-number matcher, handed the same structures, returned a confident
value for a cavity 7.7 A away. A refusal instead of a wrong number is the entire
defensibility gain here.

(Ensemble composition, **corrected**: **four** of the five apo entries differ
from wild type, not three — 1A8M is R31D, 2ZJC is **both** K98R and R31A, 2E7A
is K98R, 5TSW is Y56F. Only 1TNF is wild-type at all three positions. And the
mutation caveat does **not** attach to the SPD304 site: in holo 2AZ5 the nearest
Lys98 heavy atom is **8.74 A** from ligand `307`, and residue 56 is **7.82 A**.
Neither is in the 5 A shell — residue 98 does not line the SPD304 site at all.
The K98R concern is real, but it belongs to the *other* pocket, the on-axis
cavity above. Report ensemble composition either way.)

### fpocket druggability is a WITHIN-STRUCTURE quantity, read out of the source

**This is the finding the section below was an early, narrow instance of.** It
was written as "mdpocket cannot report a druggability". The general statement is
stronger and applies to **every** druggability number this skill emits.

`pocket.c`, `set_normalized_descriptors`, lines **736-756**:

    mean_loc_hyd_dens_norm = (mlhd - mlhd_min) / (mlhd_max - mlhd_min)

with `mlhd_min` and `mlhd_max` accumulated **over the current structure's own
pocket list**, taken whenever `n_pockets > 1`. `pscoring.c:325` feeds the result
into the logistic. The hardcoded PDB-wide constants at `pocket.c:780` —
`(mlhd - 8.23) / (24.20 - 8.23)` — are the **single-pocket branch, and it never
fires on anything we scan**: our structures carry **4 to 324** pockets.

So the score answers **"how does this pocket rank against the others in this
structure"**. It never answered "how druggable is this pocket in absolute terms".
The mdpocket case below is that statement at n_pockets = 1; the RORgt 4NB6 /
6C1P pair under "Measure, in priority order" is the same statement at n_pockets
in the hundreds, on one protein and one site, with a 90-fold gap.

**What this licenses and what it forbids:**

| operation | legal? |
| --- | --- |
| rank the site pocket among that structure's pockets, report rank + count | **yes — this is the reportable form** |
| compare fpocket rank against PRANK rank in the same structure | **yes** — two within-structure orderings, same footing |
| range across the D sweep **within one structure** | yes, with the caveat that D changes the population |
| compare a druggability value between two structures | **no** |
| compare a druggability value to a threshold | **no** |
| pool druggability into a min/max across an ensemble | **no** |
| take the max over a structure's pockets as a target value | **no** — r = 0.702 with pocket count |

### `mdpocket.sites.*` reports NO druggability, and that is the honest answer

`druggability_by_structure` used to be populated from mdpocket's `volume_score`
column. Observed values were **3.35 to 4.00 on IRAK4 and 4.36 to 4.57 on NLRP3** —
impossible for a score bounded at 1, and matching the `volume_score` descriptor
exactly on both. A plausible number under a field name that invites it to be
quoted as something else is the worst class of bug this pipeline produces, and it
was quoted.

**There is no right column to swap in.** mdpocket's characterisation table is
fixed at 22 descriptors plus 20 amino-acid counts (`M_MDP_OUTP_HEADER` in
fpocket's `mdpocket.h`) and none of them is a druggability score. Nor can it be
reconstructed: fpocket's shipped score (`pscoring.c`, `drug_score_pocket`) is

    sigmoid(-9.5699 + 7.4798*mean_loc_hyd_dens_norm
            + 0.3696*as_max_dst - 0.04672*surf_pol_vdw22)

and `mean_loc_hyd_dens_norm` is **min-max normalised across the other pockets of
the same structure** (`pocket.c`, `set_normalized_descriptors`). So a
druggability score is not a property of a pocket — it is a property of a pocket
*relative to the pocket population it was detected with*. A fixed grid has a
population of one, and the normalisation has no referent. Applying fpocket's
single-pocket fallback constants to the 6OIM switch-II row gives a saturated
1.000, which is not a measurement either.

**This is a result, not a workaround, and it generalises: fpocket druggability is
not a property of a pocket.** It is a property of a pocket *relative to the
population of pockets detected in the same structure*. A fixed grid has a
population of one, so the quantity is **undefined by construction** — not merely
unavailable, not missing from the output, not something a future mdpocket release
might add.

**Do not "fix" this by applying fpocket's single-pocket fallback constants.**
`set_normalized_descriptors` has a branch for structures with only one pocket
that substitutes `(mlhd - 8.23) / (24.20 - 8.23)`, fitted on a PDB-wide pocket
distribution. Applied to the 6OIM switch-II row (`mean_loc_hyd_dens` 185.78) it
gives a normalised 11.1 and a **saturated 1.000** for every structure. That is a
number, it is in range, it would pass the assertion, and it means nothing. It is
the most tempting wrong fix here, which is why it is written down.

So the field is `null` with `druggability_status: "not_available"` and a reason.
The descriptor is reported under its own name as `volume_score_by_structure`.
**Take druggability from the fpocket path only**, and note that it is a
per-structure number, not an mdpocket fixed-site one.

Every `[0,1]` field now passes a range assertion before it leaves the function.
A score named as a probability that comes back at 4.00 should never have escaped.

**And the loss is small, because the quantity mdpocket cannot report is not the
quantity anyone wanted.** What is missing here is an *absolute* druggability, and
fpocket never had one either — see the section above. What mdpocket can give you
is the site volume in a fixed grid, which is absolute, comparable across the
ensemble, and the thing the fixed grid was built for.

**Do not read this as "volume separated all 15 targets perfectly", which an
earlier version of this paragraph said.** That separation is **retracted** — see
the failure modes and `CLAUDE.md` rule 4a. Volume is the reported number because
it is the right *kind* of number, not because it classifies. It does not
classify; it also fails the clustering-sensitivity test worse than druggability
did (492 Å³ within-structure swing against a 139 Å³ between-group difference,
ratio 3.53 versus druggability's 1.49).

### mdpocket's own failure modes, all confirmed by direct test

It is the right tool and it is quiet about being wrong. Five things:

**1. Silent frame dropping — check this before reading any grid.** A missing
file in the list prints a message and then **exits 0**. The resulting `freq.dx`
is normalised over the frames that actually ran, so a dropped structure
**silently inflates every frequency in the grid** — the failure looks like a
stronger result. The only detector is that `time.txt` carries exactly one line
per processed frame. **Assert `len(time.txt) == len(list.txt)` before reading
any grid.** Non-negotiable.

**2. It does not superpose.** Unaligned input exits 0. Two *different proteins*
in one list also exit 0, with a non-fatal warning. Superposition and site
definition are the caller's job, and nothing downstream will notice if you skip
them.

**3. Frequency is quantised at 1/N.** With N=5 the only attainable values are
{0, 0.2, 0.4, 0.6, 0.8, 1.0}, so a genuine 1-in-5 signal is indistinguishable
from single-structure noise. **Require N >= 10 structures, or do not report a
frequency at all** — report presence/absence and say so.

**4. `_all_atom_pdensities.pdb` uses the first structure's topology only.** It
is meaningless whenever atom counts differ across the ensemble, which is the
normal case for crystal structures with different disordered regions.

**5. Superposing a homo-oligomer requires searching chain permutations.** For a
C3 trimer the three cyclic mappings agree within **0.03 A**, while the three
anticyclic ones give **~22 A** and must be rejected. Take the best mapping; do
not assume A→A, B→B, C→C.

**6. Deposited entries do not share a residue numbering, and a length check
cannot tell.** This was the worst bug the app carried. C-alpha were indexed by
raw *author* residue number, and TL1A's ensemble uses three conventions at once —
2O0O at offset 0, five entries at **+67**, 2RE9 at **+71**:

| numbering | core CA | best 3-chain RMSD, 2O0O vs the rest |
| --- | --- | --- |
| raw author | 67 | **18.70 – 20.06 A** |
| aligned | 138 | **0.51 – 1.45 A**, clean C3 split (1.3 vs 22.7) |

The ensemble superposes essentially perfectly. What was reported was
`2QE3: best chain mapping RMSD 14.84 A exceeds 5.0 A; not a superposition` —
which reads as a conformational problem and is not one. **The error message
misdiagnosed its own failure**, and the whole mdpocket stage was lost on a
target class where mdpocket is the mandated method. S1PR1 failed the same way
with `only 8 C-alpha positions are common to every chain of every structure`.

The old guard could not catch it: `len(core) < 20` tests a **count**, and 67
residues aligned by accident at a constant offset clear a count of 20
comfortably. Its message even said "the entries do not share a numbering" — it
just could not detect it.

**Align numbering first, then assert residue IDENTITY at every core position,
not just how many there are.** The app now recovers the offset by voting on
residue-name agreement against the reference and drops any core position whose
residue name disagrees between structures.

### One bad structure must not cost the ensemble

`6UYA: best chain mapping RMSD 23.87 A exceeds 5.0 A` — a 4-chain assembly
against a 2-chain reference — used to abort the entire IRAK4 run, and the same
shape killed TL1A on 2QE3. **The refusal is correct; aborting is not.** A frame
that will not superpose is dropped, recorded in `mdpocket.frames_dropped` with
its RMSD and reason, and the rest continue; the common core is then recomputed
over the survivors so a dropped entry does not shrink a measurement it is no
longer part of.

**This does not weaken failure mode 1.** Our deliberate drop happens *before*
submission, and `frame_count_check` carries three numbers —
`n_input_structures`, `n_submitted_to_mdpocket`, `n_processed` — with the
assertion made against `n_submitted`. A frame we dropped on purpose and a frame
mdpocket silently lost remain distinguishable, which is the whole point.

Below **3 surviving structures once a drop has occurred** (or 2 in any case) the
run refuses rather than reporting a CV over the survivors of a partial ensemble.

#### The common core is built over the TARGET's chains, not every chain in the file

`_common_core` requires a residue number to be present in **every chain of every
structure** and to name the same residue in all of them. Run across every chain
in the file, an entry carrying an antibody Fab, a receptor ectodomain or a fusion
partner poisons the intersection with chains that have nothing to do with the
target and nothing to do with each other.

Measured on BAFF, which returned:

> *only 0 C-alpha positions are shared by every chain of every structure AND name
> the same residue in all of them (131 numbers were shared, 131 of them dropped
> because the entries disagree about which residue that number is)*

**131 shared, all 131 rejected — and the whole by-construction site definition
lost on every target in that batch.** The identity gate is right; its *scope* was
wrong. Two chains that are different proteins are *supposed* to disagree about
residue 131, and asking them to agree is not a test any ensemble can pass.
Restricted to the target's own chains it is the test it was meant to be, and
`mdpocket.core_restricted_to_target_chains` records what was used per structure.
A run with no accession keeps the old behaviour, audibly.

#### A 60-mer must refuse the same way whether or not it fits in 62 chain IDs

1OQE, 1OQD and 4V46 refused at fetch — *"more chains than PDB chain IDs"*. That
is an **implementation limit at 62, not a judgement about the structure**, and
**1JH5, also a 60-mer, sailed through**: 378 pockets, a selected pocket **60.28 Å
from the protein centre**, and `homo_oligomer.n_polymer_chains` reporting **10**
for a 60-chain file because the guard counts target chains. Plausible numbers, no
error anywhere.

Now: any entry with ≥ `LARGE_ASSEMBLY_CHAINS` (12, PROPOSED and NOT CALIBRATED —
above NLRP3's octamer, far below 60) polymer chains and **no `chains`
restriction** is refused with `tier: "none"`, naming the escape hatch. Passing
`chains` is the caller asserting which protomers carry the site, which is rule 2b
being executed rather than guessed. `homo_oligomer` now also carries
`n_polymer_chains_in_file` and `n_polymer_chains_scored` beside the target count,
so no single number can be mistaken for the file's size.

**Related, and it was silent:** gemmi names assembly-expansion copies `<orig>-<n>`
and `_struct_ref_seq` only names `<orig>`, so on 1JH5 ten chains resolved to
Q9Y275 and **the other fifty — the same protein — carried `chain_accessions:
null`** and were reported as `non_target_chains_scored`. Expansion copies now
inherit their source chain's accession. This is also why 3K51's target chains
come back as `["A","A-2","A-3"]` rather than `["A"]`.

#### But which structure is the reference decides *who* gets dropped, and it used to be `ids[0]`

Every structure is fitted onto the reference and dropped if it will not go — so
picking an outlier as the reference does not fail the outlier, **it fails
everything else**, and the payload then reports the ensemble as unusable instead
of reporting the outlier.

Measured on NLRP3 (`8SWF, 7ZGU, 9HG4`). 8SWF came first, so it was the reference;
7ZGU and 9HG4 were dropped at **16.43 and 16.55 A**, 1 of 3 survived, and the
whole stage returned `mdpocket_status: "failed"`. **7ZGU and 9HG4 superpose onto
each other at 2.69 A** by this stage's own whole-assembly fit — and at **1.301 A**
in the *cryptic stage of the same payload*, which had already measured and
printed it. One structure cost an entire ensemble by being listed first.

| candidate | median RMSD to the rest | would keep |
| --- | --- | --- |
| 8SWF | **16.49** | **0 of 2** ← was the reference |
| 7ZGU | **9.559** | 1 of 2 ← is now |
| 9HG4 | 9.621 | 1 of 2 |

The reference is now chosen as the structure with the **lowest median core-Cα
RMSD to the rest** — not the first, not the largest. The median is what makes it
robust: an outlier's own median is large by construction so it cannot elect
itself, and one bad entry cannot unseat a good reference the way a mean would.
Ties go to the choice that keeps the most structures. Every candidate's median,
its per-structure RMSDs and how many structures it *would* have kept are in
`mdpocket.reference_selection.candidates`, and the selection is **re-run over the
survivors after a drop**, because a drop changes which structure is central to
what is left.

Two things this does **not** do, both worth knowing:

- **It does not rescue NLRP3.** The best reference is 7ZGU, 8SWF is then
  correctly dropped as the outlier, and 2 of 3 survive — still below the
  3-after-a-drop floor, so that run still refuses. It refuses **for the right
  reason and names the right outlier**, which is the whole difference.
- **It does not renumber twice.** `superposition.reference` is the *frame* every
  structure was fitted onto and may be re-elected after a drop;
  `superposition.numbering_reference` is the entry whose residue numbering the
  ensemble was put on and is fixed at the first choice. They differ only after a
  re-election. Read both before comparing two runs.

Above 12 structures the O(n²) search is skipped and the first entry is used, with
the fallback stated in `reference_selection.selected_by`.

### Geometric scoring is blind to cryptic pockets — measured, on KRAS

| | 6OIM (holo) | 4OBE (apo) |
| --- | --- | --- |
| switch-II druggability | 0.708 | **0.000** |
| rank | 1 of 9 | 4 of 5 |
| volume | 585 A^3 | 230 A^3 |
| coverage of true site | 17/22 | 7/22 |

Not merely under-scored — physically collapsed. Superposing apo on holo (0.86 A
RMSD over 128 core C-alpha) and placing sotorasib gives 6 heavy-atom clashes
under 2.0 A against a self-control baseline of 1; switch-II backbone moves up to
8.8 A at Glu63.

**A low score on an apo structure means the measurement was not made.** Say it
every time, not once.

### "Cryptic" is two different mechanisms and they need different escalations

| | KRAS | TNF-alpha |
| --- | --- | --- |
| max C-alpha displacement | **~8.8 A** | **~1.6 A** |
| what blocks the ligand | **side chains, carried in by a collapsing loop** — 12 clashes at 2.0 A, all side-chain (Arg68, Met72, His95). Backbone atoms first appear at 2.5 A. | 40/66 clashes from the displaced subunit; all 26 others are Tyr119 **side-chain** atoms. **Zero backbone clashes.** |

**Both columns show side-chain clashes, so clash composition does not separate
these mechanisms — classify on C-alpha displacement instead.** KRAS's switch-II
loop moves 8.8 A and carries its side chains with it, so the atoms sitting *in*
the site are side-chain even though the *cause* is backbone motion. Keying the
classification on which atoms clash would label KRAS as side-chain occlusion and
hand the canonical nanomolar target a micromolar prognosis.
| ligand free volume, apo | — | 62.1% intact trimer / 85.3% subunit removed / **99.8%** with two Tyr119 rotamers trimmed |
| mechanism | **backbone collapse** | **steric occlusion** |
| what would resolve it | dynamics — mixed-solvent MD, bioemu | rotamer sampling; for oligomers, test the subunit-removed state |

Backbone displacement separates these robustly at every D tested. Druggability
does not. Build the risk signal on geometry.

### The cryptic call needs a superposition gate, and it needed one badly

Three targets found this independently. In every case the module returned
`cryptic_status: "ok"` on top of a fit its own output block showed was broken:

| target | what it reported | what its own superposition said | truth |
| --- | --- | --- | --- |
| **NLRP3** (8SWF) | `is_cryptic: true`, **21.6 A**, `loop_or_backbone_motion`, prior **nanomolar** | `core_ca_rmsd: 16.627` over 487 CA, `n_excluded_ca: 0`, all four mappings at 16.629 | re-run against 7ZGU superposes at **1.248 A** and gives **0.95 A, mechanism `none`, not cryptic** — the site is pre-formed |
| **S1PR1** (8G94) | `subunit_occlusion`, **0.00 A**, 28 contacts "from a displaced chain" → prior **micromolar_at_best** | `chain_mapping {"R": "F"}`, `n_equivalent_ca: 5`, `n_residue_name_mismatches: 15` — it mapped the receptor onto **CD69, a 25-residue peptide** | 257 core CA at 1.03 A gives **1.33 A, zero clashes, mechanism `none`** |
| **TL1A** | — | numbering offsets (above) | — |

The NLRP3 case is the sharpest: **mdpocket refused the identical pair in the same
payload** (`8SWF: best chain mapping RMSD 16.22 A exceeds 5.0 A`) while the
cryptic stage built a confident mechanistic call on it. Two stages, one pair,
opposite verdicts, and only one of them had a gate.

The S1PR1 case is the most damaging: `subunit_occlusion` maps through rule 5 to a
**micromolar-at-best** ceiling, on a target with 600 sub-nanomolar compounds and
five approved drugs. **Four log units wrong, in the direction that kills a
program.**

**Gate on FOUR things: core RMSD, all-Cα RMSD, fitted-Cα count and residue-name
agreement** — not RMSD alone, and not the core RMSD alone. The S1PR1 fit had a
*low* RMSD precisely because it was fitted on five atoms. The deployed thresholds
are core RMSD ≤ 5.0 A **and all-Cα RMSD after the core fit ≤ 5.0 A** (the same
value mdpocket uses, deliberately, so the two stages cannot disagree), ≥ 20
**fitted** Cα, and ≤ 10% residue-name mismatches. On failure, refuse —
`cryptic_status: "failed"` with the gate's own numbers in `superposition_gate`.

#### The fourth check exists because the first three can be defeated by narrowing the fit

`core_ca_rmsd` is scored over the **fitted subset**. Exclude enough of the
structure and it goes green while the two entries stay where they were.
Demonstrated on NLRP3, 8SWF vs 9HG4 with `fit_residue_range=(130,370)`:

| field | value | what the old gate saw |
| --- | --- | --- |
| `core_ca_rmsd` | **1.472 A** over `n_fitted_ca` 202 | ✅ under 5.0 |
| `n_equivalent_ca` | 476 | ✅ over 20 |
| `n_residue_name_mismatches` | 0 | ✅ |
| `n_excluded_ca` | **274** | not read |
| `all_ca_rmsd_after_core_fit` | **25.619 A** | **not read** |

It **passed**, and emitted `41.7 A`, `is_cryptic: true`,
`loop_or_backbone_motion` and a **nanomolar** potency prior on two structures
16.5 Å out of frame — a *worse* confident answer than the 21.6 Å the gate was
built to stop, from a fit the payload itself scores at 25.6 Å. The field that
catches it was already computed and already in the output; the gate did not read
it. It does now, and the same pair refuses with
`after the core fit, RMSD over ALL 476 equivalent C-alpha is 25.611 A (the fitted
subset was 1.472 A over 202, with 274 excluded); the fit describes a fragment,
not the pair`.

Two details worth carrying:

- **The count gated is `n_fitted_ca`, not `n_equivalent_ca`.** The floor asks how
  many positions the fit stood on, not how many were equivalent before the
  exclusions. The floor itself scales with the smaller mapped chain — see
  *The Cα floor SCALES with the smaller mapped chain* below.
- **This route is latent, not live.** Reaching it needs `fit_residue_range` or
  `exclude_residues` and `pocket_scan` exposes neither. It is not hypothetical
  either: `auto_trim` writes the same `n_excluded_ca` on its own, so on a hinged
  protein it is one convergence away from doing this unasked. The measured
  core→all gap across six real pairs is small (0.58→1.42, 0.96→1.10, 0.71→1.02,
  1.30→2.72, 1.27→2.04), which is why it had not fired and why the check costs
  the controls nothing.

#### The mismatch fraction had no denominator, and could exceed 1

The refusal on the S1PR1 CD69 fit read *"**15 of 5** fitted positions name a
DIFFERENT residue"* — a ratio of **3.0**. `n_residue_name_mismatches` is counted
over the **pre-filter overlap**, `n_equivalent_ca` counts the **survivors** of
that filter, and dividing one by the other is not a fraction of anything. The
denominator is now survivors + mismatches (the overlap the mismatches were
counted over) while `match_residue_names` is on, and the survivors alone when it
is off. The same case now reads **15 of 20 = 0.75**, and
`superposition_gate.name_mismatch_denominator` and `name_mismatch_fraction`
carry it.

It erred toward *refusing*, so nothing was ever wrongly accepted by it — but a
gate that prints an impossible number is a gate a reader will stop trusting, and
the direction it errs in is luck rather than design.

#### A large core RMSD is not the same finding as a wrong pair — it may be a hinge

The refusal used to say *"These two entries are not superposed"*, which asserts a
mismatch. **The gate cannot tell a mismatch from a hinge and must not claim to.**
NLRP3 8SWF against 9HG4:

| fit | core RMSD | fitted Cα |
| --- | --- | --- |
| whole pair | **16.503 A** | 476, **zero trimmed** |
| single apo chain only | **16.507 A** | 476, zero trimmed |
| NBD, residues 130–370 | **1.472 A** | 202 |
| residues 220–370 | **1.377 A** | 151 |
| HD1, residues 371–430 | **0.922 A** | 58 |

An order of magnitude, per domain, with nothing trimmed at the whole-pair level —
so neither a bad chain mapping nor a trimmed outlier explains it. That is a
genuine **NACHT hinge rotation** between an open octamer and a closed NACHT. The
refusal is right; calling it "not superposed" misdiagnoses it exactly the way the
TL1A numbering-offset message used to. **A large core RMSD with well-superposing
subdomains is a hinge, and is reportable as one.**

**There is no CLI flag that expresses the domain-restricted fit which would
recover such a pair.** `fit_residue_range` / `exclude_residues` live on
`cryptic_analysis` and are deliberately not exposed by `pocket_scan` — exposing
them would also expose the narrowed-fit failure above. So a hinge is currently
**refused rather than measured**, and that is a known gap, not a verdict.

#### The Cα floor SCALES with the smaller mapped chain — the debt above is paid

**Superseded: this section previously said the floor was an absolute 20, that it
had never been exercised against a legitimately small target, and that moving the
gated count to `n_fitted_ca` had taken on a debt rather than paying it. The debt
is now paid and that framing is void.** Gating `n_fitted_ca` rather than
`n_equivalent_ca` stays — it is the fix for the narrowed-fit exploit, where a fit
restricted to one domain turns every check green while the pair sits 25.6 Å out
of frame — but an *absolute* 20 on the fitted count refuses valid small pairs,
because `auto_trim`'s `min_fit_fraction` is 0.5 and a 30-residue pair can present
15 fitted. Those are real cases here: TL1A's entries run 111–270 residues and
interface partners 25–63.

The floor is now **`min(20, max(8, 0.5 × the smaller MAPPED CHAIN's Cα))`**:

| smaller mapped chain | floor | effect |
| --- | --- | --- |
| ≥ 40 Cα | **20** | **identical to before** — the entire regression set (162–476 equivalent, 135–422 fitted) lives here, so no validated result moves |
| 30 | 15 | a 30-residue pair at the trim limit now **passes** |
| 25 | 12 | S1PR1's CD69 mis-mapping still **refuses** at 5 fitted |
| ≤ 16 | 8 | hard bottom: three points determine a rigid body exactly, so a handful of Cα reports a near-zero RMSD by construction |

**Scaled on the CHAIN, never on `n_equivalent_ca`** — the equivalent count is
itself a product of the mapping, and the S1PR1 failure was *5 equivalent
positions onto a 25-residue peptide*, which a self-referential floor would have
waved through at 5 of 5. Chain sizes are read from the two coordinate files with
gemmi (they are mmCIF at this point, so a fixed-column PDB reader silently
returns `{}` and the scaling never happens).

**Constructed and measured, not reasoned about.** KRAS 4OBE/6OIM truncated to a
single small segment of chain A, run through the real gate:

| construct | chain Cα | equivalent | fitted | floor | old absolute 20 | now |
| --- | --- | --- | --- | --- | --- | --- |
| residues 60–83 | 24 | 23 | 22 | 12 | pass | pass |
| **residues 60–79** | **20** | **20** | **19** | **10** | **REFUSED** | **passes**, core RMSD 2.37 Å, 9.56 Å displacement |
| **residues 64–81** | **18** | **17** | **16** | **9** | **REFUSED** | **passes**, core RMSD 0.77 Å |
| full chain (control) | 167 | 162+ | 135+ | 20 | pass | pass, unchanged |

`superposition_gate` now carries `n_smaller_mapped_chain_ca`, `min_fitted_ca`,
`min_fitted_ca_basis` (prose saying which of the three regimes applied),
`min_fitted_ca_ceiling`, `min_fitted_ca_chain_fraction` and
`min_fitted_ca_hard_bottom`. `min_equivalent_ca` is kept as an alias and now
carries the floor **that was actually applied**, which is what a reader of a
refusal needs. The 0.5 fraction and the bottom of 8 are **PROPOSED, NOT
CALIBRATED**; the bottom rests on a geometric argument, not a measurement.

**And whatever drops a structure from the call must drop it from every statistic
derived from it — at every level that reprints them.** After the NLRP3 re-run the
aggregate reported `is_cryptic: false, mechanism: "none"` (from 7ZGU) while still
carrying `max_backbone_ca_displacement_a: 21.6` — from the **rejected** 8SWF.
Those cannot both be true. Every headline field now comes from one named
`representative_apo_pdb_id`, chosen as the best-superposed apo entry, with the
per-structure values beside it and rejected entries listed separately.

**The same bug then survived one nesting level down.**
`structures.<PDB>.cryptic` exists so a reader looking at one apo structure does
not have to find the pairwise block — which makes it the block most likely to be
**quoted** — and it was carrying `cryptic_status: "failed"` *beside*
`mechanism: "loop_or_backbone_motion"`, `is_cryptic: true`,
`max_backbone_ca_displacement_a: 21.13` and
`cryptic_potency_prior: {expected_ceiling: "nanomolar"}`. All five fit-derived
keys — `mechanism`, `is_cryptic`, `max_backbone_ca_displacement_a`,
`clash_attribution`, `cryptic_potency_prior` — are now **null** whenever
`cryptic_status != "ok"`, with `_quarantined_keys` listing them and
`_quarantine_note` saying they were computed and refused rather than missing.
The refused values remain in `cryptic.per_apo_structure.<PDB>` for diagnosis
only. `core_ca_rmsd_a` is added to the block so a reader can see *why* without
leaving it.

**The displacement figures are protocol-dependent — quote what the run
measured.** 8.83 A (KRAS) and 1.62 A (TNF-alpha) are **hand-calibration**
numbers, from a protocol that disabled auto-trim and residue-name matching and
named the mobile regions by hand. The deployed default does neither — it finds
mobile regions nobody named and drops construct differences (KRAS
G12C/C51S/C80L/C118S, TNF L143D) out of the fit — and lands 0.1-0.2 A below:
**8.65 A for KRAS, ~1.55 A for TNF-alpha**. Mechanism and `is_cryptic` are
**identical** under both protocols, so no label changes; only the decimals do.
`pocket_scan` reports the default in `cryptic.max_backbone_ca_displacement_a` and
re-runs the calibration protocol into `calibration_protocol` beside it. Say which
one you are quoting, and do not present 8.83 or 1.62 as figures this pipeline
reproduces. The 5-fold separation is the finding; the decimals are not.

The subunit-removed control is cheap and decisive for oligomers: delete the third
chain from each apo TNF-alpha trimer and all five immediately recover the SPD304
pocket at ~280–550 A^3 against a holo dimer value of ~310 (raw 281.8–546.0 and
312.5; the same Monte-Carlo volume estimator that puts ~1 pp of noise into every
CV puts the fourth digit here out of reach). In the *intact*
apo trimer the same site measures **0.00 A^3** by mdpocket in four of five — so
"pre-formed" is a statement about the two-chain state. The site is there; the
third protomer is standing in it.

The clash attribution is what makes that reading a measurement rather than a
story, and it repeats across the whole ensemble: placing the ligand into each
intact apo trimer gives **27–29 heavy-atom clashes under 2.0 A** (minimum
interatomic distance 0.28–0.53 A), attributed **identically in all five** to
chain C — S60, Y119, L120, G121 — plus the Tyr119 triad. Five independent
crystals, one answer.

### The mechanism label is a threshold crossing — read `mechanism_margin` before quoting it

`loop_or_backbone_motion` and `sidechain_occlusion` are separated by **one
comparison**: max site Cα displacement against **2.0 A**. The two sides carry
opposite potency priors under rule 5 (nanomolar against micromolar-at-best), so
that single comparison is worth four log units of prognosis — and nothing in the
payload said how far from the boundary any label was.

Measured on S1PR1, inactive **3V2Y** → active **7TD4**. The pair passes the gate
and returns `is_cryptic: false` with a site Cα RMSD of **1.04 A**, which is
right. But `mechanism` comes back **`loop_or_backbone_motion`** because a max
site displacement of **2.16 A** clears 2.00 A by **0.16 A** on one residue
(`R:293`) — and via rule 5 that maps to a **nanomolar** ceiling. Right answer,
wrong reasoning: the label was generated by a 0.16 Å margin on a single residue,
not by the motion anyone would describe.

`cryptic.mechanism_margin` now carries `value_a`, `threshold_a`, `margin_a`, the
residue it was decided at, and `decided_by_a_narrow_margin`. **When that flag is
true, quote the displacement, not the label.** The warn band (0.5 Å) is
PROPOSED, NOT CALIBRATED — set from this one case — and it flags rather than
filters: no label is changed, suppressed or re-derived. It applies only to the
two displacement-decided mechanisms; `subunit_occlusion` and `none` are decided
earlier by different quantities and the block says so.

### A large global motion with a still site is a real state, and it used to be invisible

Same S1PR1 pair. **TM6's activation swing appeared nowhere in the payload.**
`result["global"]` was null, there was no protein-wide displacement field at all,
and the only trace was `all_ca_rmsd_after_core_fit: 2.035` — which reads as
"fine". A dossier built from that payload would have stated the site is
pre-formed and never mentioned that the two conformers differ by an
activation-state rearrangement.

**An RMSD is not a maximum.** Behind that 2.035 Å sits a **14.42 Å** single-residue
displacement at `R:248`. `cryptic.motion_scope` reports both:

| field | S1PR1 3V2Y→7TD4 |
| --- | --- |
| `site_ca_rmsd_a` | 1.04 |
| `site_max_ca_displacement_a` | 2.16 |
| `all_ca_rmsd_after_core_fit_a` | 2.035 |
| `global_max_ca_displacement_a` | **14.42** at `R:248` |
| `global_motion_with_still_site` | **true** |

The global maximum is recomputed from the superposition block — same chain
mapping, same name-matching setting, same excluded positions — and
`motion_scope.reconstruction.reconstruction_agrees` is a self-check against the
reported `core_ca_rmsd`. **If that is false the number is wrong and must not be
used.** It was true on all five controls.

"Still" is tested on the **site RMSD**, not the site maximum, deliberately: the
maximum is one residue and is the same knife-edge quantity `mechanism_margin`
warns about. Keying on it would have declared this site not-still by 0.16 Å —
hiding the one case the field exists for. 2.0 Å is also CryptoBench's own
pocket-residue criterion, so the number means the same thing on both sides.

**Report both.** A large global motion beside a still site is not crypticity and
not a failure; it is the most decision-relevant thing about the pair. "The site
is pre-formed" alone omits the activation-state change.

Measured across the controls: KRAS `false` (the global maximum **is** the site
motion, 8.65 Å at `A:63`), TNF-alpha `false` (4.58 Å, below the 5.0 Å notable
threshold), NLRP3 7ZGU `true` (0.63 Å site RMSD against 17.86 Å at `A:679`),
S1PR1 receptor-to-receptor `true` (0.62 Å against 6.07 Å). The 5.0 Å notable
threshold is PROPOSED, NOT CALIBRATED.

### Report the distribution. Stop electing one pocket.

**This supersedes the "one site pocket, chosen by a basis" design.** The
`site_pocket` / `site_pocket_selected_by` fields still exist and still carry what
they always did — nothing that reads them breaks — but they are **an annotation
on one row of a table, not the answer**.

**Scoring a whole protein and getting a score per pocket is a clean
measurement.** What went wrong in the calibration set was not that sites could
not be anchored; it was mundane and it was data: the reported pocket was on
**MAX** instead of MYC, on **IL-11 receptor α** instead of IL-11, inside
**tralokinumab** instead of on IL-13, anchored on **cholesterol hemisuccinate**
and on **CHAPSO**. Wrong chain, wrong entry, wrong ligand. Every one of them was
born the same way — nothing external applied, so selection fell back to
*whatever scored highest*, and a single number cannot show you that it came from
the wrong molecule.

So `by_clustering.<D>.pocket_table` is now the primary output: one compact row
per returned pocket carrying

    rank · prank_rank · volume_a3 · druggability · score ·
    n_lining_residues · chains · chain_accessions ·
    on_target_fraction · on_target · anchors · anchor_detail · centroid

**A table of thirty rows with nine of them on chain B cannot hide what a single
elected number hides.** It also removes the selection bias by construction:
taking the max-scoring pocket of N is a maximum over N draws and grows with N.

**That bias is now measured, and for druggability it is severe.**
r(n_pockets, max druggability) = **0.702** at D=1.6 — **max-over-pockets is
substantially a measurement of pocket count.** It contaminated **70% of the hard
class** of the retracted calibration set, because on the hard side nothing
anchored the site and selection fell back to
`max_druggability_no_ligand_site`. **That path must never produce a reportable
value.**

The earlier parenthetical here read "that bias is real but was overstated". It
was not overstated; it was measured on the wrong quantity. **Volume** does not
track pocket count (r = **−0.098**), which is what that sentence was true of.
Druggability does, and it does so because the score is normalised over exactly
that population — see "fpocket druggability is a within-structure quantity"
above.

**The confound that does hold is a separate thing and is not fixed by this.**
"Has a drug-like co-crystal" separates druggable from hard at **AUC 0.900 with no
structural measurement at all.** Any metric evaluated on that split inherits it.
That is about our labels, not our pockets — reporting a distribution does not
remove it, it just stops it hiding behind one number.

#### Anchoring is an annotation, and several can coexist

Each pocket carries every external label that applies, in `anchors`:

| label | what it means | where it comes from |
| --- | --- | --- |
| `ligand_site` | overlaps a drug-like co-crystallised ligand | Jaccard, in `anchor_detail.ligand_site_jaccard` |
| `interface` | overlaps the partner epitope | folded in by the interface stage; needs `partner_structures` |
| `symmetry_axis` | built from equivalent residues of ≥2 identical chains | geometry — **no ligand needed** |
| `annotated_functional_site` | overlaps a UniProt binding/active/site feature | UniProt, aligned through `_struct_ref_seq` |
| `buried_core` | the existing geometry flag | folded in by the interface stage, with enclosure |
| `transferred_homolog_site` | a Foldseek neighbour's ligand site | **not available in this module** — it lives in `structure-select` / `neighbour_precedent`, and its absence is stated rather than left to look like a negative |

A pocket may carry several or none. **`site_hypothesis_basis: "not_established"`
now means "no pocket in this structure carries any external label"** — a true and
useful statement about a protein — instead of "we fell back to whatever scored
highest", which is how all four bad anchors were born.

`symmetry_axis` is the one that earns its keep on a target with no chemistry:
BAFF's axial site is built from Gln144/Phe194/Leu282/Leu284 contributed by three
protomers with no ligand anywhere in the entry. It is detected as ≥2 residue
*numbers* contributed by ≥2 sequence-identical chains
(`SYMMETRY_AXIS_MIN_SHARED_RESIDUES`, PROPOSED and NOT CALIBRATED — one shared
number is a chain contact, not an axis).

#### Measured: on TNF-α the symmetry-axis label finds the SPD304 site with no ligand

This was the open question — do the ligand-free annotations agree with the
ligand-anchored one where both exist? **On the canonical case, yes.**

`2AZ5` at D=2.4, rank 5: `["ligand_site", "symmetry_axis"]` together — Jaccard
**0.636** against the SPD304 site, symmetry residues **57, 59, 119, 120** across
two chains, druggability 0.587, 322.4 Å³. The two labels land on **the same
pocket**.

And in the **apo** trimer `1TNF`, D=2.4 rank 5 carries `symmetry_axis` on
residues **57, 59, 119** — the same set, druggability 0.32 — with **no ligand
anywhere in the entry**. Leu57, Tyr59 and Tyr119 are the SPD304 contacts. So the
symmetry-axis annotation recovers the site on a structure that has no chemistry
to recover it from, which is the thing this station needs to work.

**`interface` correctly *disagrees* at that pocket, and the disagreement is the
mechanism finding.** The SPD304 pocket carries no `interface` label; the
interface-labelled pockets are elsewhere on the trimer at overlaps 0.36–0.73.
That is the same statement as the established "0.00 overlap with the TNFR2
epitope — destabiliser, not orthosteric", now readable per pocket instead of
per structure. **Anchors agreeing is evidence; anchors disagreeing is a
mechanism claim.** Do not collapse them.

One caution visible in the same table: at **D=1.6** the SPD304 pocket falls to
rank 13 with Jaccard 0.069, because the channel fragments below fpocket's `-i 15`
floor — the documented D=1.6 false negative, unchanged. Read both D values.

IL-17A is the remaining half of this test and has not been run here.

#### But agreeing with a ligand is not the same as FINDING the site without one

The result above is that the symmetry-axis *label* co-occurs with `ligand_site`
where the ligand exists. The harder question — **can a ligand-free definition
find the site when there is no ligand?** — was tested separately on TNF-alpha,
and the answer separates the four definitions sharply. **n = 1 target; do not
generalise the ranking.**

| ligand-free definition | centroid distance to the SPD304 site | shared residues |
| --- | --- | --- |
| **transferred homolog** (CD40LG **3LKJ**) | **0.00 Å**, Jaccard **0.615** | — |
| TNFR2 epitope | **14.1 Å** | **zero** |
| symmetry axis | **22.4 Å** | — |
| annotated function | **20.5 Å** | — |

Two things follow, and the second is the sharper one:

- **`transferred_homolog_site` is the strongest ligand-free anchor measured
  here** — and it is the one label this module does **not** produce (it lives in
  `structure-select` / `neighbour_precedent`). Its absence from the anchors table
  is a gap in this module, not a negative about the target.
- **"On the symmetry axis" is not a site.** TNF's C3 axis carries **five**
  distinct on-axis cavities and there is **no ligand-free rule to pick among
  them**. The runner-up sits **7.86 Å** from SPD304 — independently reproducing
  the **7.7 Å** figure in the withdrawn-matcher retraction below, arrived at from
  a completely different direction. The 22.4 Å above is what happens when the
  rule picks a different one of the five.

So a `symmetry_axis` anchor standing alone establishes *a* cavity on an axis, not
*the* site, and `site_hypothesis_basis` should say so. This does not retract the
2AZ5 / 1TNF co-occurrence result — a label agreeing with a ligand where the
ligand exists is a different and weaker claim than a definition locating the site
where it does not.

#### Payload size: read `pocket_table`, not `pockets`

A dry run hit `readScanPayload`'s 180,000-character cap, dropped 98 pocket
objects, was *still* over, and truncated mid-string — producing invalid JSON, and
deleting `_handler_note` first because it is the last key. **Per-pocket records
here carry no prose.** Every explanation lives once per clustering value, in
`on_target_selection` and `anchor_summary`. A consumer under a size cap should
read `pocket_table` and drop the verbose `pockets` list, which carries the full
residue lists and per-residue names.

### The chain resolver worked and the pocket selector ignored it

**This is the defect that suspended the volume metric**, and it is the twin of
the uniformly-null field below: a number that was uniformly *present* and quietly
measuring a different molecule.

The payload announced `"target_chains_basis": "chains mapping to P11836 in
_struct_ref_seq"`, with a `_why` naming the S1PR1/Gβ1 case it was built to
prevent — and then selection picked **the most druggable pocket anywhere in the
file**, with no check that a single lining residue was on those chains.

Measured. Selected pockets that were actually on the target:

| target | on-target | what the others were |
| --- | --- | --- |
| IL-13 | **1 of 8** (reproduced here; 1 of 9 as first reported) | cavities inside the Fabs of **tralokinumab** and **lebrikizumab** (3L5X, 5L6Y, 3L5W, 4PS4), the **receptor** chain in 3LB6 and 3BPO, chain B in 5E4E |
| BAFF | 2 of 5 | 5Y9J's rank-1 pocket, druggability **0.762**, lined 81% by **belimumab** |
| CD20 | 4 of 7 | 6Y90 and 6Y97 lined by **rituximab's Fab** |

Re-measured at D=1.6 over eight IL-13 entries after the fix:

| entry | old pick | on-target | new pick |
| --- | --- | --- | --- |
| 3L5X | rank 2, 261.7 Å³ | **0.00** (chain H) | rank 7, 68.1 Å³ |
| 5L6Y | rank 1, 225.3 Å³ | **0.00** (chain H) | rank 12, 82.9 Å³ |
| 3L5W | rank 1, 328.9 Å³ | **0.00** (chain L) | rank 3, 301.9 Å³ |
| 4PS4 | rank 15, 825.3 Å³ | **0.00** (H, L) | rank 17, 209.5 Å³ |
| 3G6D | rank 2, 188.9 Å³ | 0.67 | unchanged |
| 3LB6 | rank 1, 312.6 Å³ | **0.00** (chain C) | rank 3, 94.1 Å³ |
| 3BPO | rank 4, 317.3 Å³ | **0.00** (chain C) | rank 1, 361.6 Å³ |
| 5E4E | rank 3, 311.8 Å³ | 0.42 | rank 6, 102.4 Å³ |

**Median volume 312.2 → 145.7 Å³.** (The first report measured 312.3 → 106.8 on a
slightly different entry set and hand-filter; the *old* medians agree to 0.1 Å³,
so the disagreement is in how the corrected value is taken, not in the defect.
Quote the payload, not either of these.) BAFF and CD20 move the same way.

**So: `on_target_residue_fraction` is now computed for every pocket, and only
pockets above `POCKET_MIN_ON_TARGET_FRACTION` (0.5, PROPOSED and NOT CALIBRATED)
are eligible to be selected as the site.** Read
`by_clustering.<D>.on_target_selection` before quoting any volume. It carries the
selected pocket's fraction, the off-target chains lining it, and three counts:
`n_pockets_on_target`, `n_pockets_off_target`, and `n_pockets_fully_on_target`.

Three things to know about the rule:

- **A majority, not unanimity.** A genuine orthosteric pocket at a target/partner
  interface is legitimately lined by both, and requiring 1.0 would refuse exactly
  the pockets rule 2b exists to find. The failures being caught are not marginal
  — they sit at **0.00**.
- **`n_pockets_fully_on_target` is reported so you can be stricter than the
  module.** BAFF 5Y9J has **0 of 22** fully on-target and exactly one above 0.5
  (at 0.667, druggability 0.000). Under `≥0.5` that entry contributes 118.8 Å³;
  under `==1.0` it contributes nothing. Both are defensible; the module applies
  the looser one and gives you the number to apply the stricter.
- **An entry where no pocket is on-target contributes NOTHING.**
  `site_pocket_selected_by` becomes `no_on_target_pocket`, `site_pocket` is null,
  and no volume or druggability is emitted — rather than emitting a partner's.
  A sixth value therefore exists on `site_pocket_selected_by`.

**Off-target pockets are still returned**, with their fractions, in `pockets`. A
cavity inside an antibody is a real cavity. It is just not this target's site and
must never be quoted as its volume.

### The seven values of `site_pocket_selected_by`, and the four that do not identify a site

| value | what it means | poolable? |
| --- | --- | --- |
| `ligand_site_jaccard` | the pocket overlaps a drug-like co-crystallised ligand's contact shell | **yes** — same-site without qualification |
| `site_signature_overlap` | matched on shared residue NUMBERS from a donor's site | yes, if `collapsed_by` is 0 and `foreign_polymer_residues_dropped` is 0 |
| `site_signature_unreliable_homooligomer` | the numbers collapse across identical protomers | **no** |
| `site_signature_unreliable_foreign_polymer` | the numbers came from a DIFFERENT polymer | **no** |
| `no_pocket_matched_site_signature` | nothing matched | no site |
| `no_on_target_pocket` | no pocket sits on the target's own chains | no site |
| `max_druggability_no_ligand_site` | "the most druggable pocket anywhere on the target" | **no** |

**`site_signature_unreliable_foreign_polymer` is new and it counts a different
failure from the homo-oligomer one.** `collapsed_by` counts residue numbers lost
to **identical protomers**; `foreign_polymer_residues_dropped` counts numbers
**imported from a different polymer**, and only the first was ever guarded. The
site signature is a set of residue numbers with chain identity discarded, so a
second polymer in the file contributes its own numbering: on 8QFZ the contact
shell of `LFI` is **13 residues, 9 of them (11, 12, 13, 14, 17, 18, 19, 21, 22)
belonging to a 12-mer bicyclic peptide numbered from 1**, and those numbers also
exist on TSLP and mean something else entirely. Matched chain-agnostically they
land on a different part of the protein — `max_radius_difference_a` 33.52 Å,
per-structure site ranks 5/7/22/1/36 — while `collapsed_by` reported **0** and
the homo-oligomer guard correctly did not fire.

The signature is now filtered to the **target entity's chains** at both places it
is built, and the basis becomes `site_signature_unreliable_foreign_polymer` when
**≥⅓ of the donor residues were foreign or fewer than 6 survive** (8QFZ: 13 in,
9 foreign, 4 left). Chain-agnostic matching stays right for an inter-subunit site
on ONE protein — TNF-alpha's axial channel — and is wrong across two proteins.
Both counts are in `ensemble.site_signature`, beside `collapsed_by`.

### `lining_by_chain` is on every pocket, always — read it first

Every pocket carries `lining_by_chain` (`{chain: n lining residues}`) and
`lining_fraction_non_target`, unconditionally, whether or not the accession
resolver succeeded. It is the cheapest field in the payload and it is the one
that makes a wrong-protein pocket visible instead of invisible:

| case | what `lining_by_chain` shows |
| --- | --- |
| MYC 6G6J, the retracted 188 Å³ anchor | `{"B": 6}` — pocket 3 at 187.8 Å³, **zero MYC atoms**; chain B is MAX (P61244) |
| IL-13 3L5X | the largest pocket, 302.3 Å³, is `{"L": 8}` — entirely inside the antibody light chain |
| 8QFZ | the anchored pocket is `{"A": 4, "B": 6}` — **6 of 10 lining residues are the peptide** |

It is **necessary and not sufficient**, and the limit is worth stating: it
catches a pocket sitting on the *wrong chain of a mixed file*. It cannot catch a
**wrong PDB ID** whose every chain is the wrong protein — IL-11's 6O4P (sole
entity Q14626, the receptor) and RORgt's 6C1P (sole entity A8EVM5, an ion
channel) show one chain each and look ordinary. Those two are caught by the
accession check, which now **refuses** them outright (`entry does not contain
P20809`). Read both fields; they fail in different directions.

### A polymer LIGAND is stripped before scoring, and the result is reported as a PAIR

Rule 4 says strip every ligand before scoring. `_prep` kept `het_flag == 'A'`,
which is *every* polymer, so a peptide, a nanobody or a designed mini-binder was
never stripped and lined the pocket it was being scored in.

**How prep tells a ligand chain from a partner chain** — four cases, decided in
this order, and the order matters because chemistry and the caller's assertion
both outrank size:

| class | test | scored file |
| --- | --- | --- |
| `target` | its `_struct_ref` accession is the run's target | **keep** |
| `caller_asserted` | named in the `chains` argument (rule 2b) | **keep** |
| `polymer_ligand` | non-target **and** (hosts a component `ligand_filter` called `polymer_conjugate`, **or** ≤50 monomers, **or** lines ≥25% of the anchored pocket) | **strip** |
| `partner` | anything else non-target | **keep**, and its lining share is reported |

Homo-oligomers are safe **by construction** — TNF-alpha's three subunits are the
same entity with the same accession, so all three are `target` and the trimer
axis survives — and an unverified chain set strips **nothing**, because a chain
that cannot be shown to be non-target must not be deleted.

**Report the PAIR. Never one number.** `structures.<PDB>.polymer_ligand_control`
carries `volume_a3_stripped` (**the reported value**, and what fills
`pocket_volume_a3.primary_d1_6_a3`), `volume_a3_with_polymer_ligand` beside it,
`lining_fraction_from_polymer_ligand`, and `induced_fit_signal`.

**A site that collapses when its polymer ligand is removed is an INDUCED-FIT
signal, not a low number.** This project already has the calibration: KRAS
switch-II is druggability 0.708 on holo 6OIM and **0.000 on apo 4OBE at the same
site**, and reading the apo number as a verdict is the thirty-year KRAS error. On
8QFZ the pair is **283.6 Å³ with the peptide and no site at all without it** —
one pocket in the whole protein, 147.8 Å³, 15-18 Å away — and the pair is far
more informative than either number. **`induced_fit_signal: true` forces
`cryptic_pocket_risk: high`** and a `tractability.caveat` saying the geometric
number cannot be read in either direction without a holo *small-molecule*
structure. On 8QFZ the honest verdict is `insufficient_evidence`, for better
reasons than the run originally had.

### The covalent context, and the file it must come from

`ligand_filter` only applies its covalent rules when it is handed a
`StructureContext`, and `pocket_scan` now builds one per entry and passes it at
all three call sites — the site donor, **the signature-donor loop**, and the
per-structure pass. The signature-donor loop is the one that matters most: it is
where `dl[0]` picks the component whose contact shell becomes the site signature
for the whole ensemble.

**Build it from `files.rcsb.org/header/<ID>.cif`, never from the assembly file.**
RCSB strips `_struct_conn` from assembly files exactly as it strips `_struct_ref`
— verified on 8QFZ, whose `-assembly1.cif` carries ~22 categories and neither. A
context built from the coordinate file this module already holds returns an empty
link table that is **indistinguishable from "nothing is covalently bonded"**, and
`LFI` goes straight back to `druglike`: measured, the assembly-built context
gives `druglike` with flags `multi_electrophile_may_be_a_crosslinking_reagent`
and `struct_conn_absent_from_context`, while the header-built context gives
`polymer_conjugate`. The header is already fetched once per entry for the
accession, so this costs **no extra network call**. `holo_call.context_applied`
and `has_struct_conn_category` say whether the rules could run at all; treat any
verdict reached with them false as unestablished.

`holo_call.polymer_conjugates` and `holo_call.polymer_ligand_precedent` are
carried into the payload verbatim. The second is a deliverable, not a diagnostic:
on 8QFZ it names the peptide, its length and its sequence
(`CHWLENCWRGFC`, 12 monomers, modality `peptide`), and that belongs in the
dossier's **peptide** precedent block under rule 1. 8QFZ *is* evidence — of a
different modality.

### Failing open on the accession is what let the selector loose

`_target_chains` used to return **every** chain with the note *"no chain of this
entry maps to P35225; using every chain scored, which may include partners"* — a
warning in a string, downstream of nothing. It now **fails closed**: an entry
that declares UniProt references and contains none matching the target is
refused, with `tier: "none"` and the declared accessions named.

**Refusing on a string comparison is not safe, and two real cases prove it.**
Both were caught before shipping, and both would have thrown away valid data:

- **UniProt merges accessions.** TL1A's 2O0O, 2QE3, 2RJK, 2RJL and 2RE9 declare
  **Q8NFE9**, whose record reads `inactiveReason: {inactiveReasonType: "MERGED",
  mergeDemergeTo: ["O95150"]}`. It *is* O95150. A literal comparison would have
  refused **five of six** entries of the target's own ensemble.
- **One gene can have several live accessions.** IL-13's 3BPO declares
  **Q4VB50**, which is not merged into P35225 and never will be — it is an
  unreviewed TrEMBL entry named "Interleukin-13", gene `IL13`, human. Refusing
  3BPO would have been wrong in the damaging direction; what the entry needed was
  for chain A to be recognised as target so that B and C — IL-4Rα and IL-13Rα1 —
  are recognised as **partners**.

So a declared accession matches the target if it is the same string, **or**
UniProt has merged one into the other, **or** they are the same gene in the same
organism. And when UniProt cannot be reached, the entry is **not refused and not
verified** — an unanswered question is not a negative answer. Three statuses now
exist where there was one: `verified` (checked and matched), unverified
(mapping unreadable, or no accession supplied), and refused.

**`target_chains_verified: false` disables the on-target filter**, deliberately.
An unreadable accession mapping must not silently drop pockets — `on_target` is
null throughout and `on_target_selection._unverified_note` says the selected
pocket has not been shown to be on the target.

### Which chain is the target is a lookup, not the longest one

Anything that identifies the target by chain **length** is wrong the moment a
partner is longer, and on a GPCR–G-protein complex it always is. Measured on
S1PR1: G-beta-1 is 331–338 residues against the receptor's 278–290 in all four
entries. The interface stage therefore split 7TD4 into target `["B"]` and partner
`["A","G","R"]` — **chain B is G-beta-1, chain R is S1PR1** — computed the
G-beta/G-alpha–G-gamma interface, reported it as the target's epitope with
`interface_status: ok` and 93 interface residues, and warned about nothing.

The same longest-chain sequence was the disorder fallback, so without an explicit
accession disorder would have been computed on G-beta-1.

**And the length heuristic was still deciding in the interface stage after that
fix, because it decided what `target_seqs` was.** IL-13 3BPO split into
`target ["A","C"]` against `partner ["B"]` — A is IL-13 and **C is IL-13Rα1 at
314 aa** — and `side_a` came back containing `C:ASN240, C:PHE259, C:TYR276`,
receptor residues reported as the target's own epitope, with
`interface_status: "ok"` and no warning. 5E4E splits identically. It did not fire
on BAFF only because BAFF's receptor fragments are 31–63 aa, so BAFF is the
longest chain anyway: the heuristic was still guessing, it was just guessing
right.

The partner entry's own `_struct_ref_seq` is now read and the split made on it
whenever it resolves; the 5-mer sequence split is the fallback, flagged
`target_partner_split_verified: false`. Where both are available and **disagree**,
`target_partner_split_disagreement` reports the disagreement rather than
resolving it — a chain the sequence split claims and the accession does not is a
partner subunit similar enough to fool a 5-mer overlap, which is exactly the
IL-13/IL-13Rα1 failure.

**Resolve chains by UniProt accession**, which every entry declares in
`_struct_ref` / `_struct_ref_seq`. Two traps in doing so:

- **the assembly mmCIF does not carry `_struct_ref`.** RCSB strips it. Fetch
  `files.rcsb.org/header/<ID>.cif` for it.
- **a PDB entry's accession need not be the string you asked for, and refusing on
  the string is wrong.** Two live cases, both caught before shipping. UniProt
  **merges**: TL1A's 2O0O/2QE3/2RJK/2RJL/2RE9 all declare **Q8NFE9**, whose record
  reads `inactiveReason: {inactiveReasonType: "MERGED", mergeDemergeTo:
  ["O95150"]}` — a literal comparison would have refused five of six entries of
  the target's own ensemble. And **one gene can carry several live accessions**:
  IL-13's 3BPO declares **Q4VB50**, an unreviewed TrEMBL entry named
  "Interleukin-13", gene `IL13`, human, which is not merged into P35225 and never
  will be. Match on the string, **or** a UniProt merge, **or** same gene + same
  organism — and when UniProt cannot be reached, do not refuse.
- **that header file is not valid mmCIF as served.** It is the full entry with
  the coordinate loops deleted, and the deletion leaves bare `loop_` keywords
  with no tags — at the end of the file and in the middle of it (4OBE has three
  in a row). gemmi rejects the whole document. Strip any `loop_` not followed by
  a tag before parsing, or every accession lookup silently returns empty and you
  are back to longest-wins.

**The second trap is invisible from every angle except the parse.** The fetch
returns HTTP 200 and ~100 kB of plausible mmCIF. `gemmi.cif.read` raises, the
`except` returns `[]`, and the caller falls back to longest-chain **without a
word** — the accession machinery is present, wired, and doing nothing. It was
caught only by testing the parse rather than the fetch, and everything the
accession fix buys depends on that parse working: without it we analyse Gβ1
instead of the receptor, measure disorder on the wrong chain, and let a partner's
homodimer disqualify the target's site signature.

**So: assert on the value you came for, not on the transport.** A 200 is not a
parse, a parse is not a populated field, and a populated field is not the field
you needed. The SMILES trap below is the same failure in a different module, and
neither would have shown up in any status field.

The homo-oligomer guard has the same dependency: 8G94 reported
`is_homo_oligomer: true, identical_chains: ["F","G"]` — that is the **CD69
homodimer**, 25 and 27 residues, a partner — and disqualified an apo structure
whose rank-1 pocket matches the holo pockets at Jaccard 0.79/0.94/0.94. Measure
it over the target's chains only.

### `match_by: "seqid"` needs a residue-name check or it is silently wrong

Matching pocket residues to interface residues on residue **number** is only
legal if both entries number the protein the same way. Measured on TL1A, where
they do not:

- 2O0O at D=2.4 reported shared `A:HIS118`. 3K51's residue 118 is **THR**118.
- 2RE9 reported shared `A:THR34, A:PRO35, A:THR36` against 3K51's **VAL**34,
  **VAL**35, **ARG**36 — a spurious `overlap_fraction 0.227`, flagged
  `borderline`, with no numbering warning.
- The three entries sharing 3K51's convention are genuine and name-match.

Right where the numbering agrees, silently wrong where it does not, **no signal
either way**. A one-line residue-name assertion catches all of it; the geometric
`min_distance_to_interface_a` is unaffected, which is why the destabiliser call
still stood. `pocket_vs_interface.<D>.numbering_check` carries the identity
fraction and the mismatching positions.

#### Flagging it was not enough — the flagged numbers propagated upward unmarked

`numbering_agrees: false` fired on **exactly** the two corrupt structures (2O0O
at 2/69 = **0.029**, 2RE9 at 7/139 = **0.050**, against **0.993–1.000** for the
three valid ones) and then stopped. The classifications built on those illegal
matches travelled up into the field this SKILL.md tells callers to quote:

- `per_structure_consensus["2RE9"] = "allosteric_candidate"`, derived
  **entirely** from a 0.227 overlap on `A:THR34/PRO35/THR36` — VAL/VAL/ARG in the
  partner;
- `per_structure_consensus["2O0O"] = "mixed"`, from the `A:HIS118`-vs-THR118
  artifact.

`_aggregation_rule` said nothing about numbering. **Two fixes, both applied.**

**1. The offset is recovered and applied before the match, not merely flagged.**
A constant numbering offset between two depositions of one protein is the normal
case — TL1A carries three at once (0, **+67**, **+71**) and IL-17A carries +23 —
and it is recoverable by one vote over residue names. `interface_analysis`'s own
module docstring says `detect_numbering_offset` "should be run before any
cross-entry `match_by='seqid'` comparison"; the interface stage never ran one,
while the **mdpocket stage of the same payload had already recovered exactly
those offsets**. It does now: the partner epitope is renumbered onto our entry's
numbering before `classify_pocket` sees it, and
`pocket_vs_interface.<D>.numbering_offset_to_partner` records the offset, the
before-and-after identity fractions and whether it was applied.

The offset is applied **only** when it converts an illegal comparison into a
legal one over ≥ 20 shared positions *and* strictly increases the number of
name-agreeing positions — so an entry already on the partner's numbering is never
shifted, and a spurious offset over a handful of positions cannot buy agreement.

**2. What is still flagged after that is excluded from the consensus.** Anything
reaching the aggregation with `numbering_agrees: false` is a case no constant
offset fixes, and it is now dropped from `consensus`, from
`classifications_seen` and from the run-level `classification`. It is never
silently dropped: `excluded_numbering_mismatch` names the clustering values,
`classifications_excluded_numbering_mismatch` names the labels they would have
contributed, and `numbering_agrees` is carried onto `per_structure_consensus`
itself so the flag travels with the label.

A **sixth** consensus value exists for the case where every classification on a
structure was excluded: **`numbering_mismatch_not_interpretable`**. It is *not*
`no_pocket_to_classify` — pockets were found and classified; the classification
is what could not be trusted. The geometric fields
(`min_distance_to_interface_a`, `enclosure`, `subunit_enclosure_gain`) are
unaffected and remain in `per_structure`.

### A near-sealed hydrophobic pocket is a domain core, not a site

IRAK4's death domain gave the **top-ranked pocket of 134, druggability 0.890** —
and it is the hydrophobic core of the domain. Lining: nine Leu/Ile/Val/Phe, one
Arg, one Tyr. `enclosure = 0.998` (sealed, no solvent mouth),
`subunit_enclosure_gain = 0.020` (partner chains contribute nothing to the
burial, so it is buried within one chain), `interface_coverage = 0.026`.

fpocket's druggability regression rewards exactly that shape — large, sealed,
greasy — so a core will outscore a real site. The supporting fields caught it,
which is the system working, but nothing in the payload said so.

`buried_core_suspected` now fires on the **geometry**, never on the score:
enclosure ≥ 0.98, subunit gain ≤ 0.05, apolar lining fraction ≥ 0.7. **These are
a PROPOSAL, not calibrated** — one observed case, no held-out set. They gate a
flag and nothing else: no pocket is dropped, reordered or rescored, and a flagged
pocket still carries its rank and its score. Read the flag as "this druggability
value is uninterpretable", not as "this pocket is not there".

### A holo ligand may be a frequent hitter, not a drug

2AZ5's ligand — chemical component **`307`**, not "SPD304" as PDB titles suggest
— is a bis-electrophilic compound widely regarded as promiscuous and cytotoxic.
Its site scores 0.346 at best, well under KRAS's 0.708, which is consistent with
a micromolar tool compound rather than a drug. **A holo structure is evidence
that something bound, not that the site is drug-tractable.** Check what the
ligand is.

### An asymmetric unit is not a biological assembly

2AZ5's ASU is **four chains — two independent TNF dimers**, each with its own
copy of `307`. REMARK 350 calls it "tetrameric", which is crystal packing. Run on
the ASU and the top-ranked pocket (0.298) is a pure crystal-contact site at a
B–D interface. Pick the biological unit deliberately and record it.

### Apo does not mean ligand-free

4OBE is apo only with respect to *drug-like* ligands — it carries GDP and Mg.
1TNF has no HETATM records at all. "Apo" is a property relative to the site being
asked about. State which ligands are present rather than using the label.

### Holo vs apo is a chemistry question, and no size floor or list can answer it

The old test was `>= 18 heavy atoms` plus two hardcoded comp_id lists. Both are
now deleted, and both halves failed measurably:

- **ADP has 27 heavy atoms. So does `A1IPJ`, the genuine inhibitor in 9GU4.** No
  threshold separates them, ever.
- Identity filtering gave **16 holo / 8 apo** on NLRP3 where a naive
  molecular-weight window gave **19 / 5** — three false holo entries, a 19%
  overstatement.
- `CPS` (CHAPS, 615 Da) was simply missing from the list and passed straight
  through. **Adding it would have been a stopgap, not a fix**; a denylist cannot
  be complete, and the next detergent is the same bug again.
- The same shape produced wrong answers on CD20 (sterol tails on Y01/CLR,
  phospho-plus-two-acyl on PC1), KRAS neighbours (2UK read as purine + ribose +
  phosphate) and IL-17A neighbours (L44's 21-carbon chain).

`ligand_filter.classify_record` reads the component's actual structure from its
SMILES graph: **259/262 on ground truth, 61/70 on a blind held-out set with zero
false positives**, and it reproduces the deleted cofactor list without having
been shown it. Every remaining error is conservative — it calls a drug a cofactor
rather than the reverse.

Two behaviours to preserve wherever it is wired in:

- **`unknown` is not `apo`.** An unclassifiable component leaves the entry at
  tier **`undetermined`**, a third tier. Reporting it as apo is the same class of
  error as reporting a credential failure as "no data".
- **A lookup failure is not a CCD miss.** Records that could not be retrieved
  carry `lookup_failed` and land in `holo_call.undetermined`, so a flaky network
  cannot silently render holo structures apo.

Note it does **not** exclude 2AZ5's `307`: that comes back `druglike` with a
`promiscuity_advisory` flag. A frequent hitter is still a ligand; promiscuity
belongs to falsification, not to the holo call.

#### The classifier is useless without SMILES, and it fails silently

**This is the most dangerous thing on this page.** `ligand_filter` classifies on
the component's SMILES graph. Hand it records with no SMILES and it correctly
returns `unknown` for **every** component — and `unknown` is not `druglike`, so
every structure comes back apo or `undetermined`. The payload is well-formed,
every `<stage>_status` says `ok`, and the entire ensemble is silently holo-free.
Verified directly: MOV, GDP, ADP, CPS and `307` all return
`unknown — "the CCD row has no SMILES, so no chemistry test can run"` when the
record carries only type, name, formula and weight.

| record source | carries SMILES? |
| --- | --- |
| RCSB REST `data.rcsb.org/rest/v1/core/chemcomp/<ID>` | **yes** — `pdbx_chem_comp_descriptor`, type `SMILES_CANONICAL` or `SMILES` |
| Paperclip `pdb_v.chemcomps` | **yes** — the `smiles` column (the classifier's own default source) |
| CCD ligand file `files.rcsb.org/ligands/download/<ID>.cif` | **yes** |
| **the entry's own mmCIF `_chem_comp` block** | **NO — id, type, name, formula, formula_weight, and nothing else** |

The last row is the trap. It is the obvious source to reach for, because the file
is already fetched and parsed and on disk, and it is the one that does not work.

`modal_app.py` refuses it: `_assert_records_carry_smiles` raises
`LigandSourceError` when records *were* retrieved and not one carries a SMILES
string, and that exception is deliberately **not** caught by the per-structure
handler — a misconfigured record source is a run-level fault and must kill the
run rather than produce a full, clean, holo-free payload. A genuine 404 caches as
"no record" and does not trip it; a network failure has its own `lookup_failed`
path and its own `undetermined` tier.

**The general shape, which is worth carrying beyond this one case: test that
your source returns the FIELD the consumer needs, not that the fetch succeeded.**
Both of the invisible bugs on this page are that shape — see the header-file trap
above, where the fetch returns 200 and 100 kB and the parse silently yields
nothing.

Inside the Modal image the records come from
`data.rcsb.org/rest/v1/core/chemcomp/<ID>` rather than Paperclip, because the
`paperclip` binary is not in that image.

### A fusion chaperone inflates the pocket count and never fabricates an answer

Worth knowing so nobody spends a day on it. On the T4-lysozyme-fusion S1PR1
structures 3V2W and 3V2Y, fpocket puts **6 of 30 pockets entirely on the
lysozyme**, one of them inside the top five. But the **maximum druggability of
any lysozyme pocket is 0.003**, and the top-ranked pocket is on the receptor at
both clustering values. The chaperone inflates pocket count by roughly 30% and
never produces a druggable false positive.

It does leak elsewhere: **3 of 14 Foldseek neighbour accessions** on that target
were BRIL, thioredoxin and haemagglutinin.

Note also that a chain flag alone would not fix 3V2Y — the lysozyme sits *inside*
chain A at residues 1002-1161, alongside the receptor at 16-330. A residue-range
selection is needed as well as a chain selection.

### Missing residues near the site invalidate the pocket

6OIM chain A is missing 105–107, far from switch-II, so it does not matter there.
A gap adjacent to the site changes its shape, volume, and enclosure. Record
missing residues and whether they neighbour the pocket.

### Symmetry copies of one ligand can classify differently — aggregate, never take the first

Two copies of the same ligand in one structure can land either side of the
interface-overlap boundary. Measured on **8DYG, ligand U5Q**: copy A classified
`allosteric_candidate` at overlap **0.22**, copy B `orthosteric_candidate` at
**0.36**, both flagged `[borderline]` against the 0.25 boundary. The module is
being honest — the pocket genuinely sits on the boundary — but a caller that
reaches into `pocket_vs_interface.per_structure` and takes whichever copy came
first is tossing a coin between two different mechanistic claims.

**The rule for the caller:** quote `pocket_vs_interface.classification` (the
consensus over every pocket classified in the run) or
`per_structure_consensus[<pdb_id>]`, never a single per-D entry. When they
disagree the value is **`mixed`**, and `mixed` must be reported as `mixed` — say
the pocket sits on the boundary and give both overlaps. Do not collapse it to one
label, and do not pick the one that matches `mechanism_hypothesis`.

### A disorder number measured on the construct is about a different molecule

IRAK4 returned **0.0 over 284 residues** — the crystallised kinase domain — where
the full 460-residue protein is **0.1413**, with a disordered region at 101-162.
The old code only used the full sequence if the *caller* passed
`uniprot_accession`; with it omitted it silently fell back to the deposited
construct. A deposited construct is the ordered part of a protein **by
selection**, so that is not an understatement, it is an answer to a different
question — and a bare `0.0` in `tractability.disorder_fraction` reads as "no
disorder", not as "not measured".

Two changes, both binding:

- **The full-length path is the default wherever an accession exists**, and one
  usually does without the caller supplying it: every entry declares its own in
  `_struct_ref` (see the header-file traps above). When several are declared,
  the accession present in the **most entries of the ensemble** is the target —
  partners, fusion chaperones and scaffolds vary, the target does not — and a
  tie is reported as ambiguous rather than resolved by depositor ordering.
- **The construct-only path never populates `disorder_fraction`.** The number
  goes in `construct_disorder_fraction` alongside `scope`,
  `is_full_length_sequence: false` and `n_residues_measured`, so it cannot be
  read as the protein's. Quote it as "disorder *x* over *n* residues of the
  crystallised construct (*source*)", or supply the accession and re-run.

### A disorder number needs its method attached

`disorder.py` falls back when metapredict is unavailable, and the fallback is not
the same number: the deployed Modal image has metapredict and returned **0.3419**
where a local environment without it fell back to MobiDB and returned **0.277** —
a 23% difference on one target. The module behaved correctly (it warned, recorded
`method`, and never returned 0.0, per the cardinal rule). But a disorder fraction
quoted without `disorder.method` beside it is not comparable to any other
disorder fraction. Always carry the method.

### Never difference two centroids that were never superposed

`ensemble.site_centroid_control.max_pairwise_centroid_distance_a` is **removed**,
not nulled. It differenced pocket centroids across structures this module does
not superpose, so it was a real site displacement *plus* two arbitrary
rigid-body offsets — an IRAK4 run reported **103.9 A**, under the heading of a
control. It was not a measurement of anything.

Comparing pockets across structures without a common frame is the exact error
that retracted the 651-fold claim, and a caveat printed beside the number did not
stop it being quoted. Use **`max_radius_difference_a`** — each pocket's distance
from its own structure's protein centre, differenced across structures. It is
frame-invariant, it measures the same thing, and it already existed.

A related thing that is **not** a bug: a `site_pocket_centroid` of exactly
`[5.75, 5.75, 5.75]`, with `centroid_spread_across_clustering_a: 0.0`. That is an
on-axis pocket in an assembly whose 3-fold runs along the body diagonal — 2QE3's
assembly operators are literally `x,y,z` / `z,x,y` / `y,z,x`, so any C3-symmetric
cavity has equal coordinates by construction. It is the crystal frame showing
through, which is the same reason cross-entry centroid distances are meaningless.

### One pocket per structure makes rule 2b unsatisfiable

The app used to return only the selected site pocket. Rule 2b asks for **every**
detected pocket to be classified against the interface, so it could not be
satisfied from the output at all: the IRAK4 agent re-ran fpocket locally to see
the other 133 and reproduced the app's counts exactly, so the data existed and
was being thrown away.

Worse than lost data — on TL1A, 2RE9 reported `n_pockets: 31` while carrying only
rank 28. The agent could not tell whether the axial cavity was **absent** in that
structure or merely **unselected**, so it could honestly report neither a
persistence figure nor a zero. **A truncated payload does not just lose data; it
makes an honest answer unavailable.**

`by_clustering.<D>.pockets` now carries the top 30 by fpocket rank plus the
selected site pocket whatever its rank, each marked `is_site_pocket`. What was
left out is stated: `pockets_omitted`, the omitted rank range, and
`pockets_omitted_summary` bounding their maximum volume, druggability and site
overlap — so a reader can *check* that nothing large or site-overlapping was
hidden rather than take it on trust. Silent truncation reads as completeness.

Interface classification runs on the top 10 ranks plus the site pocket
(enclosure casts 512 rays per probe point per chain, so all 134 is not
affordable), with `n_pockets_not_classified` and the reason stated. The residue
lists for every returned pocket are present, so the overlap half of rule 2b can
be computed from `interface_residues` without re-running fpocket.

## Output

**Parse the CLI's JSON with care, or use `--out`.** `modal run` writes its own
progress banner to stdout *before and after* anything the entrypoint prints, so
`modal run modal_app.py ... > out.json` produces invalid JSON and stripping a
prefix is not enough — there is a trailing `Stopping app...` too. Either pass
`--out <path>` (the payload never touches stdout) or read with
`json.JSONDecoder().raw_decode(text[text.index("{"):])`. Without `--out` the
payload goes to **stderr**, which the banner does not share.

**Every parameter is reachable from the CLI**, including the two that matter most
on an oligomer:

```bash
modal run modal_app.py --pdb-ids 1TNF,2ZJC --chains '1TNF=A,B;2ZJC=A,B' \
    --mdpocket-site-donor 2AZ5 --ligand-codes 307 --out scan.json
```

`--chains` was missing, which made the **subunit-removed control unreachable**
without editing the file — and that control is the single experiment separating
"the cavity is too small" from "a protomer is standing in it". On TNF-alpha the
SPD304 site measures 0.00 A^3 intact and ~280-550 A^3 with a protomer deleted.
TL1A's axial cavity was reported at 49.5-141.1 A^3 intact and the control was
never run, because the CLI could not ask for it.

Fill the dossier's `tractability` block: volume with its ensemble spread (an
absolute quantity, comparable across structures), **druggability as
`site_pocket_rank` — fpocket rank, PRANK rank, `n_pockets` and
`structure_pdb_id`, which is a within-structure claim** — lining residues with
chain IDs, overlap with any annotated or ligand-derived site,
`cryptic_pocket_risk`, `cryptic_mechanism`, and a `caveat` naming what this run
could not see. Include the method block — tool, version, every D value swept,
and which PDB entries formed the ensemble.

**`ensemble.druggability.min`/`max`/`fold_range` NO LONGER EXIST, and the
sentence that used to stand here — "they stay populated for consumers that read
them, and `_comparability` must say which range you built" — is VOID.** It kept a
type error alive behind a label. Pooling a within-structure-normalised quantity
across structures is not a weak measurement, it is the wrong operation: fpocket
min-max normalises the score's dominant term over the *current structure's own*
pocket list (`pocket.c:736-756`; the hardcoded constants at `:780` are the
single-pocket branch and never fire at 4-324 pockets), so RORgt scores 0.827 at
its site in 4NB6 and 0.009 at the same site in 6C1P purely on which other pockets
co-existed in each file. A caveat beside the number did not stop the number being
quoted; a real run of the TNF-alpha pair emitted `fold_range: 195.7` with two
warning strings attached to it.

Two keys replace them, and both are legal by construction:

- **`ensemble.druggability.site_pocket_rank_by_structure`** — `{PDB@D: {fpocket,
  prank, n_pockets, structure_pdb_id, selected_by, druggability_score}}`. This is
  what fills the dossier's `tractability.site_pocket_rank`. The same object is
  also emitted per structure per clustering value at
  `by_clustering.<D>.site_pocket_rank`, which is where to read it when you are
  already looking at one structure.
- **`ensemble.druggability.druggability_range_within_structure_across_d`** — a
  min/max/fold_range **within one structure across the D sweep**, which is the
  only range the quantity supports.

`ensemble.volume_a3.min`/`max`/`spread_pct` are unaffected and stay: volume is an
absolute physical quantity and it does travel between structures.

When the ensemble number came through mdpocket, the method block also carries
**how the site was established** (grid definition, not residue matching), the
**number of frames actually processed** against the number submitted, and — if
any frequency is reported — **N**. A frequency from N < 10 is not reportable as
a frequency; give presence/absence instead.

Never return a druggability figure without its structure tier, its D value, **the
PDB ID it was measured in, its rank, and that structure's pocket count** beside
it. Separated from those, the number is not interpretable — the population it was
normalised against is the missing half of it.

(An earlier version of this sentence asked for "its ensemble spread" instead. **A
druggability spread across an ensemble is not a measurement** and asking for one
is asking for the 651x error. Rank plus count replaces it.)

And a volume of **0.00 A^3 is a result, not a failed run.** Report it. It is the
one output that cannot be an over-claim, and substituting the nearest pocket
that *does* have volume is how this skill got a headline finding wrong.

### A summary map that was all-`null` in every run ever made

`pocket_vs_interface.per_structure_consensus.<PDB>.overlap_by_clustering` read
the key `pocket_interface_overlap` off the raw classification dicts, where the
field is called `overlap_fraction` — `pocket_interface_overlap` is the name it is
*renamed to* one level down, when it is copied into
`structures.<PDB>.pocket_vs_interface.<D>`. So the map came back
`{"1.6": null, "2.4": null}` for every structure, at every clustering value, in
every run, while the overlaps it was summarising (0.267, 0.133, 0.111, 0.045)
sat one nesting level away. Fixed.

**Worth generalising: a field that is null everywhere reads as "not measured"
and never gets challenged.** No status flag was wrong, no exception was raised,
and no check in the regression covered it — it was found by reading a payload
the regression had already passed. When a summary field is uniformly null across
a whole run, check that it is not simply looking up the wrong key.
