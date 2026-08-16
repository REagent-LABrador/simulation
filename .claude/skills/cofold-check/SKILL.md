---
name: cofold-check
description: >
  Runs the four GPU structure tools — Boltz-2 cofolding, Boltz-2 affinity,
  ESMFold and BioEmu — through one callable module, returning each tool's
  numbers with its provenance, its per-run controls (seed dispersion, the site
  the seeds converged on, a positive-control log error, inter-chain contact
  counts) and the sample size behind every observation carried alongside, plus
  a side-chain repacking step that turns BioEmu's backbone-only frames into
  full-atom structures a pocket finder can read. It does NOT find pockets, does
  NOT score druggability, does NOT rank targets, and does NOT apply a
  calibration or correction to any predicted value.
---

# cofold-check

One module, four functions, all in `predict.py` next to this file:

| function | proto-tools key | what it is for |
| --- | --- | --- |
| `cofold_complex(sequences, ligand_smiles=None, ...)` | `boltz2-prediction` | pose and geometry for a site you already found |
| `cofold_affinity(protein, ligand_smiles, ...)` | `boltz2-affinity` | **triaging** binders from non-binders within one target (benchmarked; ordering actives against each other is **not** supported — see failure mode 8) |
| `esmfold_predict(sequence, ...)` | `esmfold-prediction` | fast fold + the model's own confidence; usable on complexes **behind a pTM ≥ 0.80 gate** (see failure mode 4) |
| `bioemu_ensemble(sequence, n_samples, ...)` | `bioemu-sample` | backbone conformational ensemble |

Every function returns a `dict` carrying the numbers, a `provenance` block, and
the per-run controls. None of them returns a bare score.

And one CPU-only module beside it, `repack.py`:

| function | what it is for |
| --- | --- |
| `repack_pdb(pdb_text, sequence)` | rebuild side chains on one backbone(+CB) frame |
| `repack_frames(frames, sequence)` | the same over an ensemble, one bad frame not costing the rest |
| `sidechain_rmsd(pred, ref)` | the step's own error against a known full-atom structure |
| `strip_to_backbone_cb(pdb_text)` | manufacture a BioEmu-shaped frame from a real structure, which is how that error gets measured |

`repack.py` is what makes `bioemu_ensemble` reach a pocket finder at all. Read
the repacking section below before using it — it has a measured, large and
**directional** effect on pocket scores, and that effect is the main thing you
need to know about it.

## How these are invoked — read this before writing any call

**Plain Python import, in process. Not MCP. Not a CLI.** The prose elsewhere
implies a tool-server; there isn't one. This is the same pattern as
`structure-select`'s Foldseek call.

```python
from predict import cofold_complex, cofold_affinity, esmfold_predict, bioemu_ensemble

result = cofold_affinity(
    protein=JAK1_KINASE_DOMAIN,
    ligand_smiles=[candidate_smiles],
    positive_control_smiles=TOFACITINIB,
    positive_control_measured_nm=0.50,
)
```

Run it under the **proto-tools python** — the venv that has `proto_tools`
installed. Set `PROTO_PY` to that interpreter and invoke `$PROTO_PY your.py`.
Importing under any other interpreter raises with that instruction.

Execution is on Modal (`device="modal"`, the default), workspace `rafwiewiora`,
environment `proto-env`, apps `proto-tools-boltz2`, `proto-tools-esmfold`,
`proto-tools-bioemu`.

**Credentials come from the environment only.** `MODAL_TOKEN_ID` +
`MODAL_TOKEN_SECRET`, or `MODAL_PROFILE` with a `~/.modal.toml`. The module
never reads a dotenv file by path — it will run in sandboxes where no such path
exists — and it raises a named error when credentials are absent rather than
failing deep inside a dispatch.

## The verbatim signatures these wrap

Recorded from the proto-tools source, not from documentation:

```python
# proto_tools/tools/structure_prediction/boltz2/boltz2.py
run_boltz2(inputs: Boltz2Input, config: Boltz2Config, instance: Any = None) -> Boltz2Output

# proto_tools/tools/structure_prediction/boltz2/boltz2_affinity.py
run_boltz2_affinity(inputs: Boltz2AffinityInput, config: Boltz2AffinityConfig, instance: Any = None) -> Boltz2AffinityOutput

# proto_tools/tools/structure_prediction/esmfold/esmfold.py
run_esmfold(inputs: ESMFoldInput, config: ESMFoldConfig, instance: Any = None) -> ESMFoldOutput

# proto_tools/tools/structure_dynamics/bioemu/bioemu_sample.py
run_bioemu(inputs: BioEmuInput, config: BioEmuConfig, instance: Any = None) -> BioEmuOutput
```

Complexes are built as `{"chains": [{"sequence": ..., "entity_type": "protein"},
{"smiles": ..., "entity_type": "ligand"}]}`. Boltz-2 advances the seed per
complex (`base_seed + dispatch_idx`), which is how `n_seeds` is implemented:
the same complex is submitted N times in one call and comes back on N seeds
through one container warm-up.

## What is a measurement here and what is an anecdote

This is the distinction the whole module is built around, so it is worth
stating before the failure modes.

**Per-run measurements** — recomputed on every call, safe to act on:

- `seed_dispersion` — ligand-centroid and backbone spread across seeds;
- `converged_site` — which residues the ligand actually contacted, and whether
  the seeds agreed;
- `control` in `cofold_affinity` — the log error against a known binder **on
  your target, in this run**;
- `interface` — inter-chain CA contact count, closest approach, COM separation;
- `cofold_control` — CA RMSD of the cofold against a crystal structure you
  supply;
- `frame_caveats.atoms_in_first_frame_THIS_RUN` — re-verifies the BioEmu frame
  format instead of asserting it.

**Cross-target observations** — carried in the payload under
`single_target_observations` / `single_compound_observations` /
`single_complex_observations`, every one stating its sample size. The first four
are now **benchmarked**; the seed-dispersion row is still one target and still
carries `benchmarked: False`. None of them is applied to any returned number, as
a correction, a gate or a downweighting:

| observation | n | what it is |
| --- | --- | --- |
| sealed-pocket confidence | **5 targets, 3 seeds/state, 30 folds** | pLDDT family (`complex_plddt`, `complex_iplddt`, `confidence_score`) drops beyond 2x seed spread on **5/5**; `iptm`/`ligand_iptm` on only **3/5**, and `ligand_iptm` *rose* on TNF-alpha (0.864 → 0.906). Magnitudes mostly too small to act on (JAK1 0.967 → 0.948). Supersedes the earlier 1-target/2-seed result |
| affinity absolute error | **23 pairs, 3 targets** (`out/claim2_*.json`, 2026-08-15) | mean signed error **+0.32 log**, 95% CI (-0.07, +0.72), p=0.12 — no systematic offset; 16 too weak / 7 too strong; MAE **0.82**, RMSE 1.01 against a ground-truth spread of **0.76** (n=17 compounds with ≥3 measurements) — error is now indistinguishable from the noise of the data scoring it. The old 1.97 was one compound vs one literature value and sits ~5 SE outside the CI; vs the 64-measurement consensus tofacitinib is +0.96 |
| affinity ranking vs triage | **12 actives × 12 decoys = 144 pairs (JAK1)** (`out/claim2_JAK1.json`, 2026-08-15, 0 run failures) | triage works: AUC **0.958** on affinity, **1.000** on binder probability, **2.13 log** separation, Cohen's d 2.41. Within-target ranking does **not**, on any of 3 targets: JAK1 rho +0.483, 95% CI (-0.05, +0.77), p=0.11 (n=12); BCL-2 +0.600, p=0.28 (n=5); EGFR +0.314, p=0.54 (n=6, **provisional** — still being repaired at the read; JAK1/BCL-2 final). Void, all mid-repair reads of this same artifact: 12×6→2.08/0.972, 12×9→2.32/0.981, 12×10→2.36/0.983, 12×11→2.27/0.977; and older still, 2.36 from 1 active vs 2 decoys |
| ESMFold at interfaces | **14 complexes x 2 constructions (28 runs)** | bimodal — 6/14 above 50% contact recovery, 6/14 at zero. pTM vs recovery rho **+0.79** (+0.57, +0.91); at pTM >= 0.80, 5/5 succeed, zero false alarms in 28 runs. The 1-contact IL-17A case reproduces **only** with the full mature chain; the ordered core gives 42% recovery |
| seed dispersion / site convergence | 1 target, 24 runs, 8 seeds of one probe | median 0.21 Å dispersion, and 21 of 24 runs on a real site that was **not the one asked about** — still `benchmarked: False` |

**Why they are still not baked in — and what the benchmark did to them.** Each
of the first three used to be a conclusion about a tool drawn from **one**
example, and this project exists to refuse exactly that reasoning. The
cross-target benchmark has now landed, and **all three n=1 claims moved**:

- the sealed-pocket claim was **overturned as stated** — the pLDDT family does
  notice, on 5/5; the treacherous metrics are `iptm`/`ligand_iptm`;
- the 1.97-log affinity bias was **overturned** — no offset is detectable over
  23 pairs across 3 targets;
- the ESMFold interface failure was **reproduced exactly and then explained** —
  it was an input-construction artifact, not a property of the tool.

Every one of them had failed in the flattering direction: each made the
instrument sound more decisive, or a limitation sound more absolute and
therefore more quotable. That is the argument for the n, not against it.

They are still not applied to any returned number, because a benchmarked
observation is a statement about a population of targets and your run is one
draw from it. `OBSERVATIONS` in `predict.py` carries the new figures with their
sample sizes; nothing downstream needed un-picking, because no correction was
ever applied to a returned number.

The seed statistics **within** a target are sound; whether the magnitude
transfers is unknown, which is why dispersion is recomputed every call rather
than assumed.

## Failure modes

### 1. Un-superposed frames fabricate a dispersion number

Boltz-2 emits every diffusion sample in **its own arbitrary coordinate frame**.
Taking ligand centroids straight off the raw CIFs measured **15.57 Å** of
"seed dispersion" between two seeds whose ligand-contact residue sets were
**identical** — a number produced entirely by the frames not being aligned. The
same pair after protein-CA superposition: **0.045 Å**.

This is the single most dangerous thing in this module, because 15.57 Å is a
plausible-looking answer to "how much do the seeds disagree" and it is pure
artifact. `_ligand_centroids_common_frame` superposes first, always, and a
rigid-body control (rotate + translate one structure by 40 Å) returns 0.000.

If you compute any cross-seed or cross-frame geometry yourself, superpose
first, and verify with a rigid-body control before believing the number.

### 2. High seed agreement is not evidence the site is right

`converged_site` returns a `caution` string on every call for this reason. On
KRAS, 21 of 24 runs landed on SI/II-P at a median 0.21 Å dispersion when the
question was about switch-II. A real site, tight agreement, wrong question.

**So always read `converged_site.consensus_contact_residues` against the site
you intended, and never treat `seed_agreement_fraction` as validation.** In the
JAK1 test the consensus was 21 residues at agreement 1.0 covering the ATP site
— hinge Glu93/Leu95, gatekeeper Met92, catalytic Lys44 in kinase-domain
numbering (offset +864 to UniProt) — which is correct, and it is correct
because the residues were checked, not because the seeds agreed.

### 3. Cofolding runs from SEQUENCE, so it cannot see your structure

Apo and holo structures of one protein usually share a sequence, so a
sequence-only cofold cannot distinguish them. This is not a pocket finder and
it is not a way to test a collapsed pocket. Use it downstream of a site you
already have.

### 4. The "97 contacts" figure is a PAIR count, not a residue count

Re-verified here against 8DYG: **97 CA-CA pairs** within 8 Å, but only **29
residues-in-contact**. Quoting a residue count against the 97 compares two
different quantities and understates the reference by 3.3×. `_inter_chain_ca_contacts`
returns both and labels which is which; `contacts` is the comparable one.

**And the "1 contact against 97" result itself was an INPUT artifact, not a
property of ESMFold.** Benchmarked on 14 complexes × 2 linker constructions
(28 runs). The old result **reproduces exactly** — 1 contact, minimum
inter-chain CA 7.30 Å, pTM 0.399 — **when you feed IL-17A's full UniProt mature
chain.** Feed the crystallographically ordered core of the same dimer, scored
against the same reference and the same contact set, and it returns **55
contacts, 42% contact recovery, complex TM 0.861, pTM 0.684.** Same tool, same
complex; only the input sequence differs. It bites on chains with long
disordered termini — TNF-alpha is unaffected either way (78% recovery on both
constructions). Linker length is minor by comparison (paired G50−G25 mean
difference −0.03, n=14).

**The behaviour is bimodal and pTM tells you which mode you got.** Median
contact recovery 0.42; **6 of 14** complexes above 50% (HIV-1 protease dimer
0.90, SOD1 dimer 0.90, KRAS/RAF1-RBD 0.84, TNF trimer 0.78, barnase/barstar
0.72, lysozyme/Fab 0.54) and **6 of 14** at exactly zero. pTM vs contact
recovery Spearman **+0.79** (95% CI +0.57, +0.91, n=28); pTM vs complex TM-score
**+0.94**.

| pTM cut | runs kept | median recovery | zero-recovery runs |
| --- | --- | --- | --- |
| none | 28 | 0.414 | 10 |
| ≥ 0.60 | 18 | 0.708 | 2 |
| **≥ 0.80** | **5** | **0.873** | **0** |

**Zero false alarms in 28 runs** — nothing below pTM 0.60 recovered ≥50% of
contacts. The one false confidence is trypsin/BPTI (pTM 0.752, recovery 0):
both chains fold well but BPTI docks on the wrong face, and that is why the
usable gate is 0.80 rather than 0.70.

**So: trim the input to the ordered region, gate at pTM ≥ 0.80, and treat
pTM < 0.6 as "no answer" rather than as a negative result.** A
separated-monomers result on a protein with disordered termini is a prompt to
re-run on the core, not a finding about the complex.

### 5. BioEmu rejects multimers outright — the linker is the only route

`BioEmuInput` validates `comp.num_chains() != 1` and raises *"BioEmu only
supports single-chain proteins (monomers)"*. `bioemu_ensemble` therefore joins
chains with a poly-glycine linker (default 8, range 5–10) and records in
`linker`: that one was inserted, its sequence and length, the **0-indexed
residue range of every linker**, and the range of every original chain.

Verified on a 2×60 construct: linker at residues **60–67**, chains at **0–59**
and **68–127**.

**A linker changes what the ensemble means.** It is a covalently tethered
construct, not the biological assembly: inter-chain distances are constrained
by the tether and the relative-orientation distribution is not the free one.
Strip the recorded linker ranges before any pocket detection or RMSD, and never
report an inter-chain measurement off these frames as free-solution.

### 6. BioEmu's sanity filter crashes when it actually rejects a frame

Reproduced, twice. When the physical filter rejects frames, upstream writes a
`*_unphysical.xtc` and then dies parsing its own filename:

```
ValueError: Invalid suffix '_unphysical.xtc'
```

It surfaces as `TypeError: Tool 'bioemu-sample' result does not conform to
BioEmuOutput: ... ensembles Field required` — and **the `_unphysical` string
does not survive into that message**, even down the `__cause__` chain, so a
naive `except` matching on it will not fire. `bioemu_ensemble` walks the whole
exception chain and also matches the output-shape signature.

A glycine-linked construct is exactly the input most likely to produce
rejectable frames, so **the multimer path walks into this every time**. The
module retries once with `filter_samples=False` and records
`filter_fallback.triggered = True`. When it fires, **the returned frames were
not sanity-checked** — clashes and chain breaks may be present. Filter on
radius of gyration, SASA and secondary-structure sanity before scoring.

The monomer path (KRAS 169, 8 samples) did **not** trigger it.

### 7. BioEmu frames have no side chains and no confidence

Confirmed again on both test runs: **835 atoms / 169 residues** (monomer) and
**628 atoms / 128 residues** (linked multimer) — about 4.9 atoms per residue,
i.e. backbone + C-beta only. Residues are zero-indexed and all B-factors are
zero, so there is **no per-frame confidence to read**.

fpocket and mdpocket define pockets from side-chain atoms, so **these frames
must be repacked before pocket detection**. Frames arrive pre-superposed, so no
alignment step is needed. `frame_caveats.atoms_in_first_frame_THIS_RUN`
re-checks this per run — divide by `residues_folded` and confirm it is ~5.
Re-confirmed on a 30-frame run here: **834 atoms / 169 residues = 4.93**.

`repack.py` is the fix, and **the inflation is now measured, not asserted**. On
the same 30 BioEmu frames, scored identically at D = 1.6 against a
ligand-anchored KRAS switch-II site:

| | frames with a pocket at the site | largest site volume |
| --- | --- | --- |
| raw frames, unrepacked | **10 / 30 (33%)** | 416 Å³ |
| after repacking | **3 / 30 (10%)** | 221 Å³ |

So skipping the repack triples the apparent hit rate and nearly doubles the
apparent volume. **Both of those are artifacts of absent side chains**, and an
agent that scored raw frames would have reported a cryptic-pocket opening that
does not survive giving the protein its side chains back. This is the single
most likely way to manufacture a false positive out of BioEmu.

### 8. Nothing here emits a potency

`cofold_affinity`'s primary output is `ranking`. The absolute value is returned
under `absolute` as `affinity_pred_value_log10_ic50_um`, marked `is_a_kd:
False`, `is_a_potency_measurement: False`, `benchmarked_against_measured_affinities:
False` and `correction_applied: None`.

**It has now been benchmarked against measured affinities — 23 pairs across
JAK1, EGFR and BCL-2 (`out/claim2_*.json`, 2026-08-15) — and the result is still
"do not quote a potency", for a different reason.** There is no systematic
offset to correct for (mean signed error +0.32 log, 95% CI (-0.07, +0.72),
p=0.12), but MAE is **0.82 log** — a factor of 6.6. So do not compare it against a nanomolar threshold and do not
quote it as a Kd or an IC50. The `benchmarked_against_measured_affinities:
False` flag on the returned payload means *this run* was not calibrated against
measurements, not that the head is unstudied.

**And `ranking` is the primary output for a narrower reason than it used to
be.** The benchmark separates two uses that the earlier text ran together:

- **Triage — supported.** Actives vs decoys on JAK1: ROC AUC **0.958** on
  affinity and **1.000** on binder probability, **2.13 log** separation,
  Cohen's d 2.41 — **12 actives × 12 decoys = 144 pairs**, from
  `out/claim2_JAK1.json` as of **2026-08-15** with zero remaining run failures.
  Quote this with its n. Four earlier values (2.08/0.972 at 12×6, 2.32/0.981 at
  12×9, 2.36/0.983 at 12×10, 2.27/0.977 at 12×11) are **void** — each was read
  while a repair pass was still recovering decoys that had failed to run.
  **A decoy that failed to run is not a decoy that scored badly:** the six
  missing decoys were tautomer-matching failures, and dropping them shrank the
  effective n while flattering the AUC. The verdict is the same under all five
  counts.
- **Ordering actives against each other — NOT supported.** Three targets, all
  positive, **none significant**: JAK1 **+0.483**, 95% CI (-0.05, +0.77), p=0.11
  (n=12); BCL-2 **+0.600**, p=0.28 (n=5); EGFR **+0.314**, p=0.54 (n=6,
  **provisional** — that artifact was still being repaired at the 2026-08-15
  18:46 read and its n is still growing toward 12; JAK1 and BCL-2 are final, and
  the JAK1-only triage figures above are unaffected). Every
  interval includes zero. **The pooled +0.564 (p=0.005, n=23) must not be
  quoted** — it is inflated by between-target potency offsets.
- **The untested case.** All of the above is **diverse chemistry only**. No
  congeneric series could be assembled — Paperclip's statement timeout blocks
  the `GROUP BY assay_id` needed to find one — and a congeneric series is
  exactly the setting where chemists would use ranking and where it would most
  plausibly look better. Read this as *not supported and not yet tested where it
  matters*, not as *shown to fail*. **The missing series, not the missing
  compounds, is the real limitation.**

So read `ranking` as a binder/non-binder split, not as a series ordering. Ranks
are **within-target only** — a pooled ranking manufactures rank correlation out
of the offset. A ranking of one ligand is not a ranking, and `ranking_usable` is
`False` when fewer than two ligands scored.

### 9. Run the positive control, and read it as a check not a calibration

Pass `positive_control_smiles` + `positive_control_measured_nm` and the control
runs in the same call. It reports `log_error` and sets `reliable` against a
1-log threshold (dossier rule 12).

Measured on JAK1 / tofacitinib: predicted −1.333 log10(IC50 µM) against 0.50 nM
measured, **log_error 1.968, reliable: False** — while the ranking put
tofacitinib first with **2.43 log units** of separation from an ibuprofen
decoy and binder probabilities 0.754 vs 0.107.

`reliable: False` means the absolute values are uninformative **for this
target**. It does not by itself condemn the ranking. The control's `scope`
field says plainly that it is one compound and is not applied as an offset to
the other ligands.

**Two calibration notes the benchmark added to that one-log criterion.** First,
**one control compound is not a control.** MAE is 0.82 log against a
ground-truth spread of 0.76 log (n=17 compounds with >=3 ChEMBL measurements) —
the model's error is now essentially indistinguishable from the noise of the
data scoring it — so whether a single pair lands inside or
outside one log is largely a coin flip on that pair's own measurement noise —
which is exactly how the withdrawn 1.97 was produced. Against the
**64-measurement** ChEMBL consensus rather than one 0.50 nM paper value,
tofacitinib's error is **+0.96 log**, and the `log_error 1.968` above is the
same run scored against the single literature value. Pass several controls, or
say the control is a single point. Second, **`reliable: True` licenses triage
and nothing more** — see failure mode 8; a passing control does not license
ordering actives against each other.

### 10. `cofold_control` refuses rather than mis-pairs

Pass `reference_structure` to `cofold_complex` (a `Structure`, a file path, or
raw PDB/CIF text) and it returns CA RMSD of each seed against that reference,
filling the dossier's `structure.cofold_control`.

Residue matching is strict: **if the CA counts do not match, no number is
emitted** and the reason is returned instead, because a silently mis-paired
RMSD is worse than a null. Verified — self-reference 0.0, rigid-body copy 0.0,
0.5 Å isotropic noise 0.905 (expected 0.866), mismatched chain `None` with a
reason.

`reproduces_reference` and `trusted` are returned **null on purpose**. They are
judgements, and there is no calibrated RMSD threshold to make them from. The
2026-08-15 cross-target benchmark did **not** close this gap: it measured
confidence metrics, affinity and ESMFold interfaces, and it measured cofold CA
RMSD only as a wild-type-vs-sealed contrast, never against a deposited
reference. So this remains uncalibrated, and it is now a *named* gap rather than
a pending one.

## Repacking — `repack.py`

The step that lets a BioEmu ensemble reach a pocket finder at all. **CPU only,
no Modal, no GPU, no network after the first run.**

```python
from repack import repack_frames, sidechain_rmsd, strip_to_backbone_cb

ens = bioemu_ensemble(SEQ, n_samples=32)
packed = repack_frames(ens["frames_pdb"], SEQ)      # list, one dict per frame
full_atom = [p["pdb"] for p in packed if p["pdb"]]
```

`repack_pdb(pdb_text, sequence)` returns `{"pdb": ..., "repack": {...}}`.
`repack_frames` is the ensemble form and returns `{"error": ...}` for a frame it
could not pack instead of aborting the other thirty-one. Pass the sequence
whenever you have it: the length is checked against the residue count and a
mismatch **raises** rather than packing a frame-shifted sequence.

### The engine, and what was rejected

**FASPR** (Huang, Pearce & Zhang, *Bioinformatics* 2020). Backbone-dependent
Dunbrack-2010 rotamer library, discrete combinatorial search. One self-contained
C++ program, no runtime dependencies, **0.15 s for a 169-residue chain** on a
laptop CPU. `ensure_faspr()` downloads and builds it on first use and caches it
next to this file.

**Licence: the upstream repo ships an MIT `LICENSE` file. Its README separately
says "FASPR is free to academic users."** Those two statements are in tension and
the discrepancy is upstream, not ours. Treat MIT as the operative grant but get
the tension resolved before shipping it commercially.

| rejected | why |
| --- | --- |
| **PyRosetta** `PackRotamersMover` / `FastRelax` | Best-in-class and reachable — `proto_tools.tools.structure_scoring.pyrosetta.pyrosetta_relax` is installed. **Not freely redistributable**; needs a per-user licence, so it cannot ship in this repo. Also GPU/Modal cost per frame where FASPR is free. |
| **SCWRL4** | Same accuracy class as FASPR (FASPR was benchmarked against it). Requires a registered academic licence and is **not redistributable**. FASPR is the shippable member of this family. |
| **Rosetta `fixbb`** | No binary reachable; same licence problem as PyRosetta. |
| **PDBFixer + OpenMM** | Neither is installed in `druggability` and `openmm` is a heavy dependency. More importantly `addMissingAtoms` places side chains from a **template in a fixed conformation** — it is not a packer, and its output needs a minimisation to be usable. That would have been the one place a force-field step entered the pipeline. Avoided. |
| **FAMPNN** (`proto_tools...fampnn.fampnn_pack`) | Genuinely attractive — a learned packer that returns **per-residue pSCE**, i.e. the per-frame confidence BioEmu itself does not provide. Rejected as the *default* only because it is a GPU call per frame against a free CPU one. **Worth revisiting** if per-frame confidence becomes load-bearing. |
| **AttnPacker / DLPacker** | Need torch plus downloaded weights; `torch` is not in `druggability`. Install cost not justified by FASPR's measured accuracy. |
| **PULCHRA / REMO** | Backbone-to-full-atom reconstruction rather than rotamer optimisation; lower side-chain accuracy than FASPR for the same job. |
| **Modeller** | Academic licence with a key; heavier than the task needs. |

**No force-field optimisation anywhere in this step.** FASPR is a discrete
rotamer search, not a continuous optimiser, and no RDKit conformer generation is
involved. `repack.minimisation_performed` is `False` on every call.

### The step's own error — measured, not asserted

`strip_to_backbone_cb()` reduces a real full-atom structure to exactly BioEmu's
format (N/CA/C/O/CB, 0-indexed, B-factors 0.00 — verified at 4.93 atoms/residue
against BioEmu's 4.93). Repacking that and comparing to the original measures the
packer's error with no held-out set needed.

Measured over **11 independent KRAS chains** (3GFT A-F, 4DSN, 4DSO, 4EPR, 4OBE,
6OIM), heavy-atom side-chain RMSD **excluding CB** (CB is fixed by the backbone,
so including it flatters the packer) and symmetry-corrected for PHE/TYR/ASP/GLU/
ARG/VAL/LEU:

| | median | range |
| --- | --- | --- |
| all scored residues | **1.60 Å** | 1.43 – 1.79 |
| buried only (≥16 CB neighbours within 10 Å) | **1.23 Å** | 1.08 – 1.61 |
| χ1 within 40°, all | **0.838** | 0.803 – 0.884 |
| χ1 within 40°, buried | **0.900** | 0.880 – 0.942 |

These sit on top of the published FASPR benchmark (~1.6 Å, ~87% χ1 on native
backbones), which is the check that the wiring is right. **Buried side chains are
markedly better predicted than surface ones, and a pocket lining is mostly
buried** — that is the reason this step is usable at all.

Caveat that matters: this error is measured on **crystal** backbones. A BioEmu
backbone carries its own error underneath, so 1.60 Å is a **floor**, not the
error you get on a generated frame.

### The thing you must not skip: repacking MOVES POCKET SCORES, hard

This is the dominant fact about the step and it is bigger than the RMSD.

Same structure, ligands stripped, scored by the same fpocket run at D = 1.6,
site anchored on 6OIM's MOV ligand transferred by CA superposition:

| structure | native | after strip-and-repack |
| --- | --- | --- |
| **6OIM (holo)** switch-II | **585 Å³, druggability 0.708**, 2.35 Å off ligand, Jaccard 0.773 | **294 Å³, druggability 0.003**, 2.07 Å, Jaccard 0.409 |
| **4OBE (apo)** same site | 224 Å³, druggability 0.000 | **0 Å³ at D 1.6; 594 Å³, druggability 0.444 at D 2.4** |
| 1TNF dimer, SPD304 site (D 2.4) | 537 Å³ | **0 Å³** |

Two consequences, both binding:

1. **Repacking systematically CLOSES induced-fit pockets.** FASPR packs in the
   absence of the ligand, so side chains that the crystal shows splayed around a
   bound drug relax back into the site. The holo KRAS pocket loses half its
   volume and effectively all its druggability. The TNF dimer site closes
   completely.
2. **A repacked score and a crystal score are NOT on the same scale.** The
   repacking shift on KRAS (585 → 294 Å³) is *larger than the apo-versus-holo
   difference the whole cryptic-pocket question is about*. Comparing a repacked
   BioEmu frame against an apo crystal number is therefore invalid.

**So always run the repack-only control.** Take the crystal, strip it, repack it,
score it. That is the baseline a repacked ensemble must be compared against —
never the native crystal number. And repack the **holo** structure the same way:
its score is the pipeline's own ceiling for that site, which is the only
honest thing to call a positive against. For KRAS switch-II at D = 1.6 that
ceiling is **294 Å³**, not 585 Å³.

### What the full ladder actually did on KRAS — a NEGATIVE result, stated plainly

The ladder is: no usable structure → BioEmu ensemble → repack → pocket
detection → did a cryptic pocket open? It was run end to end on the case that
settles it, and **it did not work.**

Setup: 124 BioEmu frames from the wild-type KRAS 1-169 sequence (nothing
post-1983 about that input), repacked, fpocket at D = 1.6 and 2.4, site anchored
on 6OIM's MOV ligand transferred by CA superposition. Compared against the nine
pre-2013 nucleotide-only KRAS chains (3GFT A-F, 4DSN, 4DSO, 4EPR) — the honest
2012 evidence — and against the repack-only control on each.

| arm, D = 1.6 | site volume Å³ | druggability | centroid→ligand | Jaccard |
| --- | --- | --- | --- | --- |
| 6OIM holo, native — **the answer** | **595** | 0.708 | 2.35 Å | 0.773 |
| 6OIM holo, repacked — **pipeline ceiling** | **294** | 0.003 | 2.07 Å | 0.409 |
| best pre-2013 crystal (4DSN) | 264 | 0.222 | 3.88 Å | 0.208 |
| 6 of 9 pre-2013 chains | **0** | 0.000 | no pocket | — |
| BioEmu + repack, **best of 124 frames** | 306 | 0.034 | 2.16 Å | **0.273** |
| BioEmu + repack, median of the 20 frames that found anything | 212 | — | — | — |
| BioEmu + repack, **median over all 124** | **0** | 0.000 | — | — |

Read on volume, which is the metric that separates targets: the best frame
(306 Å³) is barely above the best 2012 crystal (264 Å³) and only **16% of frames
show anything at the site at all**. That is not a pocket opening.

**The detection half of the ladder demonstrably works — which is what makes the
negative attributable.** Hand it the *holo* backbone stripped to N/CA/C/O/CB
(i.e. a synthetic "BioEmu frame" of the open state, with every real side chain
thrown away) and repack+fpocket returns **294 Å³ at 2.07 Å, Jaccard 0.409**. Hand
it the *apo* backbone treated identically and it returns **0 Å³**. So the
repack-and-detect stage separates the open conformation from the collapsed one
by 294 Å³ using no crystal side-chain information at all. It was working. It
was never given an open backbone.

**And the falsification test says the hits are not the switch-II pocket.**

- Correlation between a frame's on-site pocket volume and how close its
  switch-II backbone is to the holo conformation: **r = 0.031**. Zero.
- Frames *with* a pocket average 4.45 Å switch-II RMSD to holo; frames *without*
  average 4.37 Å. **No difference.**
- The ten frames that come *closest* to the holo switch-II backbone
  (2.29–3.14 Å) score **0 Å³ in nine of ten**.
- Jaccard against the true MOV contact shell is 0.17–0.33, against 0.773 for the
  real pocket.

So the frames that "find a pocket" find a transient loop cavity that happens to
land within 4 Å of the ligand centroid, uncorrelated with adopting the open
state. **It is a different pocket in the right neighbourhood** — the exact
failure that cost this project a headline claim before.

> **Volume alone would have called this a positive, and been wrong.** The best
> frame's 306 Å³ *exceeds* the pipeline's own holo ceiling (294 Å³), the apo
> crystal (224 Å³) and the repacked apo control (0 Å³). Reported as
> "the ensemble recovers 306 Å³ at the switch-II site where the apo crystal
> gives 0" it reads as a clean win. **The only things that catch it are the
> Jaccard (0.273 against the real pocket's 0.773) and the r = 0.031 backbone
> correlation.** Volume being the metric that separates *targets* does not make
> it a metric that identifies a *site*. Never quote a site volume without the
> centroid distance AND the residue overlap beside it.

**The generator, not the repacker, is the broken rung.** BioEmu's 124 frames span
2.29–7.33 Å switch-II RMSD to holo (median 4.19). The apo crystal is already at
3.23 Å and the 2012 crystal 4EPR is at **1.96 Å** — closer to the answer than any
generated frame. The ensemble adds diversity, not direction. The repacker was
validated at 1.60 Å and was handed nothing to pack.

Note also that D = 2.4 is unusable on KRAS: at that clustering the **holo
reference itself** returns 0 Å³ on-site, because the nucleotide and switch-II
sites merge and the merged centroid moves off-site. Quote D = 1.6 for this
target and say why.

**Leakage caveat.** BioEmu's training data is modern PDB, which includes 6OIM.
A 2012 cutoff cannot be enforced on the model's weights, only on the structures.
This makes the negative result *stronger* — the model could have memorised the
answer and still did not produce it — and it would have made a positive result
uninterpretable.

### TNF-alpha: the other mechanism class — same failure, different physics

Run because subunit occlusion is a genuinely different question from backbone
collapse. 46 frames over two runs of the glycine-linked trimer (3 × 152 aa,
linkers at 0-indexed 152-159 and 312-319, **stripped before scoring**), site
anchored on SPD304 (`307`) from 2AZ5 transferred onto the 1TNF apo trimer.

- **`filter_fallback.triggered` was `False` on all four BioEmu runs here**
  (KRAS ×2, TNF ×2). The `_unphysical.xtc` crash did **not** fire, including on
  both linked multimers, so these frames *were* sanity-checked. 5 of 16 and 13 of
  48 samples were filtered out normally.
- 1TNF apo trimer, native and repacked: **0 Å³** at both D. Site closed, as
  expected for an occlusion mechanism.
- 1TNF **dimer** (protomer deleted), native, D = 2.4: 545 Å³, **Jaccard 0.45** —
  the site is pre-formed and recoverable. This is the positive control.
- BioEmu trimer + repack: **12 of 34 Rg-sane frames (35%)** show an on-site
  pocket at D = 1.6, best 375 Å³, one frame reaching druggability 0.781.

At first read that is a win over the crystal trimer's 0 Å³. **It does not
survive the falsification pass**, and it fails it the same way KRAS did:

| check | result |
| --- | --- |
| corr(site volume, radius of gyration) | r = −0.076 |
| corr(site volume, inter-chain CA contacts) | r = 0.002 |
| corr(site volume, CA density within 12 Å of the site) | r = −0.063 |
| site-local CA density, frames **with** a pocket | **42.7** |
| site-local CA density, frames **without** | **42.6** |
| site-local CA density, apo crystal trimer | **43** |
| Jaccard of the on-site pockets (n = 12) | 0.147 – **0.242**, median 0.213 |
| Jaccard, crystal dimer positive control | **0.45** |

The mechanism for this site is *protomer displacement*. In the frames that show
a pocket, **the protomers are not displaced** — local packing around the site is
indistinguishable from the frames that show nothing and from the crystal itself
(42.7 / 42.6 / 43). And no frame's pocket reaches even 0.25 Jaccard against a
control that reaches 0.45.

One frame in 35 was grossly expanded (Rg 41.3 Å against a median of 21.7) and
was excluded. **Filter on Rg before scoring a linked multimer** — the sanity
filter passed that frame.

So: a second mechanism class, a second negative, and the same signature —
sub-threshold cavities in the right neighbourhood, uncorrelated with the physical
motion that actually opens the site.

## Cost

GPU time on Modal, warm containers, from the test runs:

| call | wall clock |
| --- | --- |
| ESMFold monomer (KRAS 169 aa) | 25.4 s |
| ESMFold dimer (IL-17A 2×132) | 3.1 s |
| Boltz-2 cofold, 1 protein + ligand, 2 seeds (JAK1 290 aa) | 64.7 s |
| Boltz-2 cofold, 2 protein chains, 1 seed (IL-17A) | 15.8 s |
| Boltz-2 affinity, 2 ligands (JAK1 290 aa) | 77.9 s |
| BioEmu, 169 aa, 8 samples | 23.9 s |
| BioEmu, 128 aa linked, 4 samples, incl. filter retry | 40.7 s |
| BioEmu, 169 aa, 32 samples (KRAS) | 147 s |
| BioEmu, 169 aa, **128** samples (KRAS) | **157 s** |
| BioEmu, 472 aa linked trimer, 16 samples (TNF) | 131 s |
| BioEmu, 472 aa linked trimer, 48 samples (TNF) | 294 s |
| **repack, 169 aa, per frame** | **0.15 s, CPU, free** |
| fpocket at one D, per structure | ~0.3 s, CPU, free |

**Sample count is nearly free; container connect is the cost.** 128 KRAS samples
cost 10 s more than 32 did — the 147 s was almost all warm-up. So do not
iterate at n = 8 to save money; take n = 100+ in one call and pay the connect
once. Sequence length is the real scaling term (169 aa → 472 aa quadrupled it).

Everything downstream of BioEmu in this skill — repacking, all fpocket sweeps,
every superposition and RMSD — is **local CPU and costs nothing**. The entire
KRAS + TNF study below spent about **12 minutes of GPU container time across four
calls, on the order of $0.50**, and several hundred pocket scorings for free.

Cold starts are much worse — the first affinity call of a session took ~250 s,
most of it container connect. MSA generation (MMseqs2) is a separate remote
service and costs no GPU, but it is a large part of first-call latency.

**Keep `n_seeds` and `n_samples` low while iterating.** `n_seeds` multiplies GPU
time roughly linearly; `diffusion_samples_affinity` defaults to 5 internally
already.

## What these four tools are actually for

Stated bluntly, because the honest answer is narrower than the tool list
suggests:

- **`cofold_affinity` — the most useful of the four, for a narrower job than
  previously claimed.** Separating plausible binders from non-binders within one
  target, with a positive control run in the same call. That is a real job and
  nothing else here does it. **It is not a way to order actives against each
  other** — within-target Spearman +0.48, 95% CI (−0.05, +0.77), p=0.11 over 17
  benchmarked pairs, and "rank candidates within a target" is withdrawn as a
  recommendation.
- **`cofold_complex` — a pose and geometry step downstream of a site you already
  found.** Its per-run value is `converged_site` (where did the ligand actually
  go) and `cofold_control` (does it reproduce a crystal structure you have). Its
  confidence numbers are structural-confidence metrics and answer "how sure is
  the model about the geometry it drew", which is not the dossier's question.
- **`esmfold_predict` — a fast folder that reports its own confidence, and the
  confidence is the point.** Cheap enough to run as a triage step. On complexes
  it is bimodal (6/14 above 50% contact recovery, 6/14 at zero) but **pTM
  predicts which** (rho +0.79, n=28) — so it is usable behind a pTM ≥ 0.80 gate,
  on the ordered core rather than the full mature chain. See failure mode 4.
- **`bioemu_ensemble` — still the weakest link, and now measurably so.** The
  repacking step it needed exists (`repack.py`) and is validated at 1.60 Å
  side-chain RMSD, so the blocker is no longer missing side chains. The blocker
  is the ensemble. On the one case where the answer is known — KRAS switch-II,
  124 frames — it did not sample the open backbone: its closest frame sits
  2.29 Å from the holo conformation while a 2012 crystal already sat at 1.96 Å,
  and pocket detection across the ensemble correlates with backbone opening at
  **r = 0.031**. It has no per-frame confidence, and the literature figure for
  generative ensembles on apo input (56% cryptic-pocket recovery, against 86%
  from holo) applies to our normal case.

  **Do not use it to answer "is this cryptic site openable".** It has not been
  shown to answer that on the case we can check. Use it for backbone diversity
  where diversity itself is the point, report the ensemble spread, and route the
  cryptic question to mixed-solvent MD instead.

## The ladder, and where it breaks

The intended escalation — apo-only or no structure → fold → BioEmu ensemble →
repack → pocket detection → does a cryptic pocket open? — **now runs end to end.
Every rung is implemented and none of them errors.** The rung that fails is
scientific, not mechanical:

| rung | status |
| --- | --- |
| fold / take structure | works |
| BioEmu ensemble | works — 124 frames, sanity filter clean, real backbone spread (median pairwise 2.6 Å) |
| **repack** | **works** — 1.60 Å side-chain RMSD, 0.15 s/frame, free, MIT |
| pocket detection across the ensemble | works |
| **"did the cryptic pocket open?"** | **NO on backbone collapse (KRAS, n=124). NO on subunit occlusion (TNF, n=34).** |

Both mechanism classes fail, and they fail identically: ~16-35% of frames report
a cavity near the site, none of it correlated with the physical motion that opens
that site (r = 0.031 and r ≈ 0.00 respectively), all of it at roughly half the
residue overlap of the true pocket. **The ladder is not currently a cryptic-site
finder for either mechanism class.** Plan around triage, not detection.

Treat anything downstream of `bioemu_ensemble` as triage and say so, with the
r = 0.031 figure attached — a reader who sees "16% of frames showed a pocket"
and nothing else will read it as a positive.

## Future integration points (NOT on the default path)

**None of these four tools is wired into dossier assembly, and this section is a
map for later, not a change.** The dossier answers a TARGET-level question — can
a small molecule bind this protein at all — from retrieved precedent and pocket
geometry on REAL structures. Three of these tools answer COMPOUND-level,
prospective questions (*will this molecule bind, and where, and how hard*) that
belong to a downstream hit-finding, screening or design station, not to target
triage. The fourth, ESMFold, is the one exception already wired — through
`structure-select/predicted_structure_fallback.py`, for the single case of a
target with NO experimental structure and no usable homolog (CLAUDE.md rule 4c).
The reason the other three sit off the path is scope first and the benchmark
second: each is conditional, GPU-bound in the separate `proto-env` Modal
environment, and a single point of failure against an otherwise CPU/stdlib
pipeline.

| tool | downstream station / question it would serve | the benchmark caveat that gates it |
| --- | --- | --- |
| `cofold_affinity` | a screening/design station: compound triage and ranking against a target | TRIAGE is supported — binder vs decoy AUC **0.958** on JAK1 (12×12, 144 pairs). Within-target ORDERING of actives is **not** supported (Spearman +0.48, 95% CI (−0.05, +0.77), p=0.11, n=12) and untested on a congeneric series. The absolute value is **not** a potency (MAE **0.82 log**). So: triage yes, series-ordering no, absolute affinity no. See failure mode 8. |
| `cofold_complex` (seed dispersion / `converged_site`) | a design station: which site a putative binder converges on, downstream of a site already found | High seed agreement is NOT validation — on KRAS 21 of 24 seeds converged, tightly (0.21 Å), on a real site that was **not the one asked about** (failure mode 2). It runs from SEQUENCE and recalls where the PDB puts ligands, so it cannot FIND a pocket — only report where a named binder lands (CLAUDE.md rule 5). |
| `bioemu_ensemble` | a cryptic-pocket station: conformational spread for open-state / cryptic-pocket hypotheses | Did NOT recover the open state on the one checkable case — KRAS switch-II, 124 frames, on-site pocket volume vs backbone-opening **r = 0.031**; the same negative on TNF-alpha subunit occlusion (r ≈ 0.00). Use for backbone DIVERSITY as triage, never as cryptic-site detection, and attach the r = 0.031 figure. See "The ladder, and where it breaks". |

ESMFold's near-term use is the one already live and is documented above and in
CLAUDE.md rule 4c: a fast fold behind a **pTM ≥ 0.80** gate, on the ordered core
rather than the full mature chain, as the structure-less-target fallback that
keeps the computed axis from nulling.
