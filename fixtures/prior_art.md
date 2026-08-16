# Prior art — what is standard, what is contested, what is ours

Researched 2026-08-15 via Paperclip full text and web. Every number here is
cited. The purpose is to stop us claiming as novel what the field settled years
ago, and to stop us quoting numbers whose provenance we do not know.

## Things we must stop claiming

**"Use an ensemble instead of one structure" is table stakes, not a
contribution.** It is the explicit design principle of CryptoBank, LIGYSIS and
HOTPocket as of 2024–25. Claim the *magnitude* we measured, never the idea.

**KRAS switch-II is the field's most-used illustration** — PocketMiner,
CryptoBank, Bowman's review and HOTPocket all use it. Good sanity check, weak
differentiator.

**"Cryptic pockets expand the druggable proteome"** is boilerplate in every
paper on the topic.

**"Cytokines are undruggable" is dead.** A 32-cytokine small-molecule-microarray
campaign against 65,000 compounds gave 864 chemotypes → 296 thermal-shift
validated binders (32.5% translation) → cellular inhibitors for IL-17, IL-13 and
IL-23 (Raevi et al., bioRxiv 2026, doi:10.64898/2026.04.20.719718).

## Things we were quoting without knowing what they were

**The fpocket druggability score is a 3-descriptor logistic regression fitted on
21 positives.** From `src/pscoring.c`, `drug_score_pocket()`:

```
score = 1/(1+exp(-(-9.5698768
                   + 7.479844   * mean_loc_hyd_dens_norm
                   + 0.3696134  * as_max_dst
                   - 0.04671833 * surf_pol_vdw22)))
```

In-code note: "21 druggable pockets vs 292". The **published** 2010
nested-logistic model (Schmidtke & Barril, *J Med Chem* 53:5858) is present in
the same file **commented out** — so the score in any current binary is not the
equation the paper describes. fpocket's own detection paper says the pocket
score "does not reflect drugability" and that rank-1 performance "drops" on apo
structures.

**Halgren's SiteMap benchmark contained exactly one PPI** (MDM2/p53). Applying
his thresholds to a curated PPI set classifies **46% of proteins as "difficult"
despite their having clinically viable inhibitors** (*Sci Rep* 2022,
doi:10.1038/s41598-022-12105-8). Our domain is PPI-heavy.

## Standard practice we are missing

| Gap | Evidence |
| --- | --- |
| **A rescorer on top of fpocket** | fpocket + PRANK or DeepPocket is the best-recall configuration in two independent benchmarks; best recall of *any* method is only 60% (Utgés & Barton, *J Cheminform* 2024, doi:10.1186/s13321-024-00923-z) |
| **A cryptic predictor** | **We ran it. It does not work on our targets — see the section below.** PocketMiner reproduces its published ROC-AUC (we measured 0.868 against 0.87) on its own test set, and inverts on ours |
| **Interface hot-spot detection** | PPI-hotspot ID: F1 0.71 / sensitivity 0.67 vs FTMap 0.13/0.07 (*eLife* 2024, doi:10.7554/eLife.96643). **27.6% of true hot spots make no cross-interface contact at all** |
| **A defensible negative-label treatment** | canSAR's PocketBagger uses positive–unlabelled learning precisely because "defining genuinely 'undruggable' pockets is nearly impossible" |
| **Cluster-aware and time-aware splits** | Random splits over homologous pockets inflate every reported metric |

**Do not use PocketMiner at all on our targets — measured, not inherited.**
HOTPocket's published complaint (5–9% DCA precision on ordinary pockets,
119 pockets per structure) is confirmed and is worse than a precision problem:
on our own four calibration cases the signal runs backwards. Details below.

**Do not build a naive consensus.** HOTPocket found that intersecting seven
pocket finders **underperformed the best single method**; ensembling only paid
off with a learned rescorer on top. This retires the consensus-of-scores idea we
sketched early on.

## PocketMiner — run in-repo 2026-08-15, and rejected

Installed from `github.com/Mickdub/gvp` branch `pocket_pred`, released weights
`models/pocketminer`. It installs, it runs in **under 0.5 s per structure on
CPU** (no GPU, no Modal), and it **reproduces its published accuracy**: pooled
ROC-AUC **0.868** on its own shipped 35-structure test set against a published
0.87, and 0.897 on the validation set. The model is genuine and correctly
installed. None of our targets appear anywhere in its dataset, so every number
below is out-of-sample in both directions.

**It still fails on our cases, and it fails by inverting.**

Operating threshold **0.536**, derived the authors' own way (Youden's J on the
shipped validation set, `src/optimal_threshold_protein_performance.py`).

At that threshold PocketMiner calls **32–54% of every residue in every
structure** a cryptic-pocket residue. "The site was flagged" is therefore
almost content-free — it fired on all four sites, including the three we know
are not cryptic. The honest test is a **spatial null**: is the real site better
than a random same-size contiguous patch of the *same* protein? 2,000 patches
per structure:

| site | ground truth (ours) | site mean | null mean | p |
| --- | --- | --- | --- | --- |
| **KRAS switch-II, apo 4OBE** | **cryptic** — 8.8 Å collapse at Glu63 | 0.632 | 0.495 | **0.203** |
| TNF-alpha SPD304, apo 1TNF | NOT cryptic — subunit occlusion, ~1.6 Å | 0.673 | 0.519 | **0.002** |
| NLRP3 MCC950, apo 6NPY | NOT cryptic — pre-formed, 0.95 Å | 0.821 | 0.483 | **0.006** |
| S1PR1 orthosteric, 3V2Y | NOT cryptic — pre-formed, 1.83 Å | 0.711 | 0.471 | **0.010** |

**The only site that fails to beat its own null is the only true cryptic site.**
All three known negatives pass. The one true positive scores *lowest of the
four* on site mean. This is not low sensitivity; it is anti-correlation with
the thing we need.

KRAS was given every fair chance — four site definitions, apo and holo:
all 22 sotorasib contacts (p=0.234), switch-II 57–76 (p=0.303), sotorasib
contacts within switch-II (p=0.088), and the five residues our own calibration
measured as moving furthest — Gly60/Glu62/Glu63/Tyr64/Tyr68 (p=0.132). None
reaches p<0.05.

**Correction to our own earlier characterisation.** This file previously said
PocketMiner "refuses multi-chain input, and failed on 38 multi-chain structures
out of roughly 220 in CryptoBench's test set". The second half is CryptoBench's
report and stands. **The first half is wrong about the released code and is
withdrawn.** There is no chain check anywhere in `pocket_pred`;
`process_strucs()` selects backbone atoms across *all* chains and flattens them
into a single residue list. It does not refuse, warn, or error.

What it does instead is worse than refusing, and we measured it. The
sequence-position encoding runs through the chain break as if it were a peptide
bond, so **the prediction depends on the order the chains appear in the file**.
Reordering 1TNF's chains A,B,C → C,B,A, with identical coordinates, changes
per-residue probabilities by up to **0.385** (mean 0.111, Pearson r 0.811) and
moves **38 of 152 residues across the decision threshold**. A multi-chain
PocketMiner number is not a well-defined function of the structure.

**And the single-protomer fallback does not rescue it** — this is the question
that decides the tool for us, since our targets are overwhelmingly oligomers.
Feed one TNF-alpha protomer instead of the trimer and the SPD304 site goes from
the top-ranked region of the structure to below its own background:

| 1TNF input | site mean | background | enrichment | site residues > threshold |
| --- | --- | --- | --- | --- |
| trimer A+B+C | 0.673 | 0.500 | **+0.174** | 34/42 |
| protomer A only | 0.392 | 0.442 | **−0.050** | **2/14** |

The same 14 chain-A site residues drop by 0.291 on average when the partners
are removed, against a whole-chain shift of only +0.034 — the loss is
concentrated at the interface, not a global offset. So for an interface site
there is no valid way to run it: the multi-chain answer is order-dependent, and
the single-chain answer cannot see the site.

**Three defects in the released code, each of which either produces wrong
answers silently or blocks the run entirely.** Recorded because anyone
re-running this will hit them:

1. **Under Keras 3 (any TF ≥ 2.16) the model does not load and does not say
   so.** `tf.train.Checkpoint.restore()` leaves **93 of 197 variables** at
   random initialisation — including the final (200,1) output kernel — while
   `util.load_checkpoint()` prints `CHECKPOINT RESTORED FROM ...`
   unconditionally and never calls `assert_consumed()`. Two runs of the
   identical script differ by up to **0.41** in per-residue probability. With
   `TF_USE_LEGACY_KERAS=1` and `tf-keras` installed: 177 variables, 0
   unrestored, bit-identical across processes. **Anyone who pip-installs
   current TensorFlow and runs `xtal_predict.py` gets confident noise.**
2. `process_strucs()` reshapes to (n_residues, 4, 3) and so cannot read any
   structure with an incomplete backbone. The 2AZ5 biological assembly (four
   copies of Leu157 with no backbone O) raises `ValueError: cannot reshape
   array of size 6492 into shape (542,4,3)`.
3. The shipped demo path is wrong: `xtal_predict.py` points at
   `../data/ACE2.pdb`; the file is at `../data/training-data/ACE2.pdb`.

**Decision: not wired in.** It would not change a single call we make. It
cannot lower `cryptic_pocket_risk` (rule 3 sets `high` on apo from tier, and a
predictor that misses the canonical cryptic site cannot license lowering it),
it cannot raise one (it fires on everything), and it cannot change
`next_experiment`, because the mechanism split that drives that decision —
backbone collapse vs side-chain/subunit occlusion — is precisely the
distinction it gets backwards here. Adding it would add a per-residue
probability we would report and not act on.

**What survives as useful:** the reproduction is a clean negative control for
the claim that single-structure cryptic prediction is ready to substitute for
an ensemble. It is not. This is the thirteenth tool whose isolated result was
more flattering than its integrated one.

## The evaluation problem — and why our as-of design matters

There is an accepted way to evaluate *pocket detection*, a partial one for
*cryptic pockets*, and **essentially nothing for druggability**.

- **NRDLD**, the de facto standard, is 71 druggable + 44 less-druggable = 115
  pockets, with a **test set of 35–37**. It is **62% positive by construction**;
  the deployment base rate — proteins with an approved drug — is **~3.5%**
  (704 Tclin of ~20,000). Every published accuracy figure (DrugPred 91%,
  DoGSiteScorer 88%, PockDrug MCC 0.885) was computed at a prior ~18× too
  generous.
- A one-class model scored **more than half of NRDLD's "less druggable" pockets
  as highly druggable** (*Front Pharmacol* 2022, doi:10.3389/fphar.2022.870479).
  Their framing: "any pocket is only non-druggable until a drug is found for it."
- **Nobody hindcasts.** No study of the form "run predictor X as of year Z,
  check whether it called the now-drugged targets druggable" was found.
- The one real time-split study is at the association layer and is damning for
  the incumbent: **OTRec** trained on Open Targets 2022 and evaluated on 2025
  trial outcomes found the **Open Targets association score scores ROC-AUC 0.559
  prospectively** — a coin flip plus epsilon (bioRxiv 2025,
  doi:10.64898/2025.12.21.695803).
- Blind prospective ligandability, for calibration: **CACHE #2 confirmed 0.7% of
  1,957 computationally nominated compounds** as binders by SPR.

So our retrospective `as_of_date` design is not a nice-to-have — it is the
missing evaluation the field has not built.

## Open Targets — measured by query, release 26.06

I claimed earlier that mechanically enforced modality separation was uncommon.
**That was wrong and is retracted.** Open Targets does it, and does it
correctly:

- `Target.tractability` returns **28 independent booleans** across four
  modalities (SM 8, AB 9, PROTAC 8, other-clinical 3). No score, no rank, no
  aggregate.
- The IL-17A / TNF test passes cleanly. `{"modality":"SM","label":"Approved
  Drug","value":false}` for both, alongside `{"modality":"AB","label":"Approved
  Drug","value":true}`. Three approved antibodies do not leak into the
  small-molecule row.

So the trap we built rule 1 around is already handled by the field's most-used
platform, at least on the clinical-precedence buckets. What survives as ours is
narrower and needs stating precisely.

**Where it is exploitably weak — all confirmed by query, not inferred:**

| Weakness | Evidence |
| --- | --- |
| The structural half is a single frozen source | Both `High-Quality Pocket` and `Med-Quality Pocket` derive from **DrugEBIlity alone**, a legacy EBI dataset with no ongoing releases |
| …and it fires wrongly on our domain | **TNF scores `SM:High-Quality Pocket = true` AND `Druggable Family = true`** — a secreted trimeric cytokine with zero small-molecule clinical candidates |
| `Structure with Ligand` conflates two questions | It requires a solved structure *and* a bound small molecule, so it cannot fire for a good apo pocket with no ligand |
| Silent all-false is indistinguishable from unassessed | **49% of 298 sampled targets have all 8 SM buckets false.** No abstention state exists |
| No confidence, no provenance, no evidence trail | The API exposes no provenance field at all; "conf" appears only inside bucket *names* |
| **No versioning — and it looks like there is** | `tractability(version:"24.06")` errors, but the URL param `?version=24.06` returns **HTTP 200 and is silently ignored** — `meta.dataVersion` still reports 26.06. Historical data exists only as 28 FTP parquet dumps |
| The adjacent numeric field is modality-blind | `Target.prioritisation.maxClinicalStage = 1` (max) for **both IL-17A and TNF**, driven entirely by approved antibodies with no modality qualifier |

**Base rates over 298 targets, useful for calibrating how much a `true` carries:**
`SM:High-Quality Pocket` fires at **10.4%** (the sharpest structural
discriminator), `SM:Approved Drug` 16.4%, `SM:Phase 1 Clinical` **0.0%**,
`PR:Database Ubiquitination` 60.1% (near-noise).

**Two things to use rather than compete with:**

1. `drugAndClinicalCandidates.drug.drugType` is a clean modality label
   (`Antibody` / `Protein` / `Small molecule` / `Unknown`) and is **more
   trustworthy than the tractability buckets themselves**. Use it as an
   independent cross-check on our SMILES-null test.
2. The 28 booleans across all targets are the **weak validation labels we
   lack**. Not ground truth, but the only broad labelled set available.

API note for whoever writes the integration: `Target.knownDrugs` no longer
exists in 26.06, and `Drug.isApproved`, `Drug.maximumClinicalTrialPhase` and
`Drug.hasBeenWithdrawn` are gone. Use `drugAndClinicalCandidates` and
`Drug.maximumClinicalStage`.

## Architecture — we are not alone, and that is good

**canSAR** keeps four axes separate (3D pocket, ligand-based, network
target-likeness, antibody accessibility) and **refuses to make negative calls**:
"The notion of a negative 'undruggable' pocket is scientifically intangible…
it is impossible to quantify the negative recall or the precision of our
predictions." That is the closest published position to ours.

**TargetDB** (SGC Oxford, GPL-3.0) emits **eight separate 0–1 component scores**
with a user-weighted MPO, and keeps fpocket-derived druggability distinct from
ChEMBL-derived chemistry. It also ships **dated SQLite snapshots pinned to
ChEMBL versions** — the best as-of story in any public tool.

**The theoretical argument for refusing to average**, which we should make
explicitly: ML pocket scorers are trained on liganded pockets, so **they already
encode precedent and cannot serve as an independent second axis**. Averaging a
precedent axis with a scorer trained on precedent double-counts.

## Where we may genuinely be novel

**Retracted from this list:** mechanically enforced modality separation. Open
Targets already does it correctly on clinical precedence (see above). Our rule 1
is still necessary — a dossier that got IL-17A wrong would be wrong — but it is
table stakes, not a differentiator. What differentiates is that we carry
*provenance and a potency figure* alongside the modality call, where Open
Targets carries a bare boolean.

1. **The magnitude of druggability-score irreproducibility — RETRACTED as a
   contribution, kept as a record.** Detection instability is published (~85%
   pocket identity under mere rotation, ~59% PDB-vs-AF2); score inflation on
   uncleaned holo is published as a distribution shift. This item used to claim
   that no fold-range of the score across an apo ensemble with volume held
   constant had been published, and offered **651× with ±16% volume** as ours.
   **That claim is withdrawn.** Volume was never held constant, because the
   pockets being compared were never the same pocket: they were matched across
   structures on shared residue numbers, and mdpocket showed the matcher was
   tracking a cavity 7.7 Å off-site with 12.2 Å of internal inconsistency
   between structures. A 19-residue reference on a homotrimer collapses to 11
   distinct residue numbers, so a C3-symmetric site is unresolvable in principle
   by that method. Both halves are void — the 651× spread and the 206.7–309.2 Å³
   volume range beside it — because both came out of the same matching step. See
   the VOID section of `pocket-scan/SKILL.md` and `tnf_result_VOID` in
   `pocket_calibration.json`.

   What survives is a smaller and better-supported claim: fixing the site by
   construction with mdpocket characterization mode, instead of matching after
   the fact, cut the ensemble volume CV from roughly 28% to roughly 10% — an
   inflation of about 2.8-fold, essentially all of it from one structure. Quote
   that to two significant figures: three identical reruns of one ensemble gave
   CVs of 12.1 / 11.3 / 10.8%, so about one percentage point of any CV is
   fpocket's Monte-Carlo volume estimator. And note what it was measured on —
   `site_from_density`, not the ligand site. It is a measurement of
   reproducibility, not of the SPD304 site.
2. **Mechanism as a routing decision.** Existing taxonomies (PocketMiner's
   forward/reverse + four backbone rearrangements; CryptoBank's
   buried/superficial × fragment/ligand) are descriptive. Using mechanism class
   to decide *which computation is valid* was not found.
3. **Subunit displacement as a first-class occlusion mechanism.** Absent from
   PocketMiner's four backbone-centric classes and excluded by CryptoBench's
   2 Å RMSD filter. **The claim that PocketMiner "refuses multi-chain input" is
   withdrawn — we ran it and it does not refuse.** What we measured is stronger
   and more specific: it accepts an oligomer silently and returns an
   order-dependent answer (reordering 1TNF's three chains moves 38 of 152
   residues across the decision threshold), while a single protomer cannot see
   the interface site at all (enrichment +0.174 → −0.050). CryptoBench's 38
   multi-chain failures out of ~220 stand as their report. So the mechanism is
   still structurally unhandleable by PocketMiner — for a demonstrated reason
   rather than an assumed one.
4. **The clustering-parameter sweep.** No paper treats fpocket's `-D` as a
   variance source to be swept.

## Definitions we must respect

**Vajda et al. 2018** (*Curr Opin Chem Biol* 44:1, PMC6088748): a cryptic site
"forms a pocket in a ligand-bound structure, but not in the unbound protein
structure", with the stringent form being absent in **all, or nearly all**
unbound structures. **Beglov et al. extended CryptoSite from 186 to 4,950
structures and found bound-like pockets partially formed in some unbound
structure for close to 50% of the 93 proteins**, with BACE-1 druggability
ranging 0.2–0.6 across 52 apo structures. So ensemble score variation was
published in 2018 — our contribution is magnitude on the raw score, not the
observation.

**CryptoBench** operationalises crypticity as pocket-residue RMSD > 2 Å. **Our
TNF-alpha case at roughly 1.6 Å sits below that threshold and the site is
present in all five apo structures — so it is NOT cryptic by either standard.**
Report it as occluded. Two significant figures on that displacement, and the
protocol matters: 1.62 Å is hand calibration, the deployed default returns
~1.55 Å, and both sit well below 2 Å, so the call does not turn on the decimals.

## Method lesson

`paperclip grep "<title>" /papers/` returns papers that **cite** a work, not the
work. Vajda 2018 looked absent from the corpus by grep (127 reference-list hits)
and was found instantly by `paperclip lookup doi "10.1016/j.cbpa.2018.05.003"`.
Use `lookup` for a known work, `grep` for named entities inside text.

**The recurring failure shape: machinery present, wired, reporting success,
doing nothing.** PocketMiner is the sharpest instance we have measured — under
current TensorFlow `tf.train.Checkpoint.restore()` silently leaves **93 of 197
variables**, output layer included, at random initialisation, while
`util.load_checkpoint()` prints `CHECKPOINT RESTORED` unconditionally and never
calls `assert_consumed()`; two identical runs then differ by **0.41**. The same
shape has now cost this project repeatedly: a rule that was never executable
through the deployed app and degraded silently to whole-assembly scoring; a
pocket matcher that confidently tracked a cavity 7.7 Å from the site it named;
an `--format csv` flag that is accepted and ignored. **A component that prints
success is not a component that ran.** For anything that loads weights, restores
a checkpoint, or accepts a flag, demand a positive assertion that it took effect
— `assert_consumed()`, a variable count, a re-read of the option — and treat a
success message with no assertion behind it as unverified.
