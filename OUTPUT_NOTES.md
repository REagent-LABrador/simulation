# Druggability dossier — output notes for consuming stations

For an engineer wiring this station into the evidence gauntlet, or assembling the
final asset-to-indication score. It describes what we emit, what each part means,
and the specific ways this output gets misread. Read the reading rules and the
worked example; the rest is reference.

Authority, in order: `.claude/skills/assemble-dossier/validate_dossier.py` (the
machine grader), `CLAUDE.md` (the deployed system prompt and the JSON template at
its foot), `rubric.md` (the grader's index). Where this document and those
disagree, they win and this is stale.

---

## What this station produces

One JSON object per run. It is an answer to one question:

> **Can a small molecule be made against this protein, for this mechanism?**

It is not an answer about a protein. It is an answer about a *(target,
mechanism)* pair — see reading rule 5. It carries no indication decision, no
ranking, no molecule, no biologic assessment and **no overall score**.

The object has exactly **16 top-level keys**, enumerated in
`validate_dossier.py:REQUIRED_TOP_LEVEL`. None may be absent; a value that could
not be retrieved is `null`, never omitted.

| key | what it holds |
| --- | --- |
| `target` | accession, gene symbol, protein name, organism, `sequence_length`, `sources` |
| `as_of_date` | ISO date or `null`. When set it is binding on every date in the file |
| `verdict` | `small_molecule_tractable` \| `not_tractable` \| `insufficient_evidence` |
| `verdict_basis` | `retrieved_precedent` \| `computed_tractability` \| `both` \| `none` |
| `axis_conflict` | free text, or `null`. Populated when the axes disagree |
| `target_precedent` | axis 1. Bioactivity, approved and clinical **small molecules**, patents, terminations, `as_of_leakage` |
| `biologic_precedent` | approved biologics. Never summed with the above |
| `family_precedent` | Pfam-family activity. Weakest transfer axis |
| `structural_neighbour_precedent` | Foldseek fold neighbours. Middle transfer axis |
| `pocket_neighbour_precedent` | pocket-descriptor neighbours. Strongest transfer axis, and a hypothesis |
| `structure` | tier, PDB entry, ligand, holo/apo census, ensemble, `cofold_control` |
| `tractability` | axis 2. Pocket geometry, cryptic call, interface class, caveat, method |
| `affinity` | predicted affinity plus its mandatory positive control |
| `falsification` | the checks we ran against our own claim, including the ones that found nothing |
| `next_experiment` | what would move the answer |
| `not_found` | one entry per thing we could not establish, and why |

Two shipped reference dossiers exit the validator at **0 violations** and are the
calibration for everything below:

- `.claude/skills/assemble-dossier/examples/jak1_P23458.json` — the filled case
- `.claude/skills/assemble-dossier/examples/tnf_P01375.json` — the refusal case

Verified in this checkout: both pass, and the validator's own suite is 103 tests,
all passing.

---

## Where to read the output: two channels, and you get the reply

The agent writes the dossier to a pinned sandbox path **and** pastes the complete
JSON into its final reply. `CLAUDE.md` requires both, every run.

| channel | path | who reads it |
| --- | --- | --- |
| file | `/mnt/session/outputs/druggability-dossier.json` | the grader, during the run |
| reply | the agent's final message text | **you** |

**Read the reply.** Sandbox files are not retrievable through the Files API after
the session ends, and the router wrapper
(`agent/tools/druggability-dossier.ts`) returns `result.text` — the reply — and
nothing else. The file is the grading channel and is gone by the time you see the
result. Parse the JSON out of the reply text; do not expect a file handle.

Two consequences:

- The reply may contain prose around the JSON. Extract the object, do not assume
  the whole message parses.
- `rubric.md` explicitly does **not** grade the reply ("this grader sees only the
  file"). So the channel you consume is the one no automated check has looked at.
  Validate what you receive: run `validate_dossier.py` on it yourself. It is pure
  stdlib and takes milliseconds.

### Invocation

`manifest.json`: `"invocation": "outcome"`, `"session_policy": "fresh"`,
`"max_iterations": 3`, model `claude-opus-5`. Outcome mode means the run is
graded against `rubric.md` and may re-iterate up to three times before returning.
The wrapper sets a **45-minute** timeout because one dossier fans out to a Modal
pocket scan (function timeout 1800 s plus cold start) plus several database
queries. Budget for that; a 10-minute default aborts mid-scan.

Sessions are fresh: nothing carries over between runs, so two runs on the same
accession are independent and may differ.

---

## Input contract

The router tool takes a single free-text `task` string. The five contract fields
below are what `CLAUDE.md` reads out of that string — there is no structured
input object today (see "Where this output is weak").

| field | required | what it does |
| --- | --- | --- |
| `uniprot_accession` | **yes** | e.g. `P01116`. A gene symbol must be resolved to an accession first, and both recorded |
| `as_of_date` | no | ISO date. **Binding**: every evidence item must predate it, filtered at the source. Turns on the `as_of_leakage` machinery — see reading rule 4 |
| `disease_context` | no | Free text. Selects relevant clinical precedent. **Never** adjusts a tractability number |
| `interaction_to_disrupt` | no | What the molecule must stop: a named partner, an oligomeric state, or a catalytic function. Determines which chains constitute the site |
| `mechanism_hypothesis` | no | `orthosteric` \| `allosteric` \| `oligomer_destabilisation` \| `unknown`. **Changes the answer** — see reading rule 5 |

What each optional field changes, concretely:

- **`as_of_date`** switches on `target_precedent.as_of_leakage[]`. Four fields
  cannot be date-filtered at source and must carry an explicit leakage flag:
  `distinct_actives` and `best_potency_nm` (the bioactivity table has no date
  column), `patents` (not filtered at source), and
  `clinical_stage_small_molecules` **unconditionally, including when the list is
  empty** — ChEMBL's `max_phase` is a current value with no phase history, so
  neither the presence nor the absence of a past clinical candidate is a
  retrievable statement. With no `as_of_date`, `as_of_leakage` is `[]`.
- **`interaction_to_disrupt`** and **`mechanism_hypothesis`** together select the
  chains that get scored, recorded in `tractability.method.chains_used`. Omitting
  them is legal and produces a *weaker* answer: the agent scores the biological
  assembly, says so in `tractability.caveat`, and asserts no pocket as the
  relevant one. Both shipped examples are in this state.
- **`disease_context`** touches only which clinical precedent is retrieved.

Supplying `uniprot_accession` alone is a valid call. Supplying a mechanism is how
you get a targeted answer.

---

## Reading rule 1 — There are two axes and no overall score. Do not create one.

The dossier reports two things that answer different questions and are allowed to
disagree:

| axis | question | where it lives | nature |
| --- | --- | --- | --- |
| **retrieved precedent** | has anyone actually made a molecule against this? | `target_precedent` | looked up. The stronger axis when it exists |
| **computed tractability** | does the structure suggest a molecule could bind? | `tractability` | computed, with declared blind spots |

`verdict_basis` names which axis carried the verdict. It exists precisely so that
one label over two axes is not a hidden average: a verdict with no basis and a
populated `axis_conflict` is an average with extra steps.

When the axes point different ways, `axis_conflict` is populated with prose
explaining the disagreement, and **the disagreement is the finding**. On
TNF-alpha (`tnf_P01375.json`) the retrieved axis reads zero approved small
molecules against five approved biologics, while the computed axis reads a
genuine ligandable cavity — reconciled mechanistically (the site is on the trimer
3-fold axis and is reached by displacing a subunit) and settled in direction by a
Phase 2 oral small molecule, balinatunfib. Collapse that to one number and you
lose everything useful about the target.

**The validator enforces this structurally.** `check_axes_never_averaged` rejects
any key at any depth named or containing `overall`, `composite`, `score`,
`averaged`, `axis_average`, `total_score`, `combined_score`, `weighted_score`,
`tractability_score`, `druggability_index`, `priority_score` and the rest of the
banned set — and rejects any numeric value under a key that mixes a precedent
token (`precedent`, `actives`, `potency`, `approved`, `clinical`) with a
tractability token (`druggability`, `pocket`, `tractability`, `volume`). There is
no field you can read to get a combined number, because it is impossible to emit
one.

**What this means for you.** If your station needs a scalar, derive it yourself,
downstream, where the weighting is visible and arguable — and carry
`verdict_basis` and `axis_conflict` alongside it so the derivation can be undone.
Do not average `distinct_actives` against `primary_d1_6_a3`. They are not
commensurable and the ordering they produce is not stable: `assemble-dossier`'s
own worked pair notes that averaging the two axes would rank TNF-alpha above JAK1
on tractability, which is the exact inversion this station exists to prevent.

## Reading rule 2 — `insufficient_evidence` is not `not_tractable`. It means nobody looked.

This is the single most important distinction in the output, and the easiest to
lose in a downstream `if`.

| verdict | claim |
| --- | --- |
| `not_tractable` | a measurement was made and came back negative |
| `insufficient_evidence` | **nobody has looked.** The question is open |
| `small_molecule_tractable` | evidence supports a small molecule, on the axis named in `verdict_basis` |

A consumer that folds `insufficient_evidence` into `not_tractable` penalises
unexplored targets exactly as hard as refuted ones, which inverts the value of
this output — the unexplored ones are where the opportunity is.

`CLAUDE.md` rule 11 makes declining a *correct* answer: "A confident score on an
unstudied target is the worst output you can return." The validator enforces both
directions:

- `check_insufficient_evidence_reachable` **fails** a dossier that reports fewer
  than 50 distinct actives, a literal `structure.holo_count: 0`, and no approved
  small molecules, and then claims any verdict other than `insufficient_evidence`.
  The agent cannot quietly upgrade a thin target.
- `insufficient_evidence` with an empty `next_experiment.resolves` also fails. A
  declination must name what would settle it. **Read `next_experiment` on every
  `insufficient_evidence` dossier** — it is the actionable part.

The shipped fixture for this shape is IL-11 (P20809): 15 ChEMBL activities, of
which **13 carry a null `standard_value`** — the entire numeric footprint is two
SPR Kd values, 140 nM and 2,600 nM, from one method on one chip. A 140 nM number
exists and is exactly the sort of thing that tempts a confident score. The
expected output is `insufficient_evidence`.

The related trap: TL1A / TNFSF15 (O95150) has zero ChEMBL targets and 7 PDB
entries, all apo or DcR3-bound, alongside Phase 3 antibodies (tulisokibart,
duvakitug). The antibody programmes give a confident-sounding precedent to latch
onto while no chemistry has been attempted. See reading rule 3.

**Also read `falsification.survived: false` correctly.** It means some recorded
check materially undercut the *precedent claim*. It does not mean the target is
not tractable, and it does not lower the verdict — the falsification skill
"attaches evidence. It does not lower a score, flip a verdict, or resolve
`axis_conflict`." TNF-alpha's shipped dossier has `survived: false` and
`verdict: small_molecule_tractable`, simultaneously and correctly.

## Reading rule 3 — An approved antibody is not evidence a small molecule is possible.

It is often evidence of the opposite: an antibody-drugged target is frequently one
where a small molecule was tried and could not be made.

Biologics live in `biologic_precedent.approved_biologics` and **never** in
`target_precedent`. The two blocks are never summed. `biologic_precedent.note`
carries the disclaimer as a mandatory template field.

Modality is decided per drug from `chembl.molecule_dictionary.molecule_type`, and
every entry in `approved_small_molecules` and `clinical_stage_small_molecules`
carries its own `modality`, for which the only legal value is `small_molecule`.
The validator (`check_modality_separation`) rejects:

- any `modality` other than `small_molecule` in either small-molecule list;
- any name appearing in both the small-molecule and biologic blocks;
- any name carrying a USAN `-mab` or `-cept` stem in a small-molecule list;
- `verdict: small_molecule_tractable` on `retrieved_precedent` or `both` with
  **both** small-molecule lists empty and no characterised potency — that is a
  biologic being leaned on, and it fails outright;
- approved biologics plus zero approved small molecules plus a tractable verdict,
  with `axis_conflict` left empty.

The motivating case is IL-17A (Q16552): three approved antibodies — secukinumab
2015, ixekizumab 2016, bimekizumab 2021 — and zero approved small molecules. A
consumer reading "approved drugs exist" off a merged count gets IL-17A exactly
backwards.

**`Unknown` is a returned value, not an absence.** Drugs whose `molecule_type`
comes back `Unknown` or NULL count toward **neither** block and are recorded in
`not_found`. Two TNF-alpha drugs (ABBV-3373, AZ9773) and two IL-17A drugs
(M-1095, CJM-112) are in this state. Do not infer their modality in either
direction, and do not treat their absence from both lists as a zero.

**Counts and lists may legitimately disagree.** Salt and parent forms are
distinct ChEMBL `molregno`s, so `approved_small_molecules_count` is a
salt-collapsed figure while `approved_small_molecules` is what could be named.
JAK1's shipped dossier reports count 9 against 8 named entries, with the gap
explained in `not_found` and the ninth deliberately unnamed rather than guessed.
**Use the count for arithmetic and the list for display; do not assume
`len(list) == count`.**

## Reading rule 3b — Four precedent axes, four separate blocks. Never merged, never discounted into one another.

Activity against something else is real signal and it is not activity against
this target. Each axis is its own top-level block, and no discount factor folds
one into another.

| block | similarity by | strength |
| --- | --- | --- |
| `target_precedent` | measured on this protein | direct evidence |
| `pocket_neighbour_precedent` | pocket descriptors plus cofold transfer | **strongest transfer** |
| `structural_neighbour_precedent` | Foldseek fold similarity | middle |
| `family_precedent` | Pfam sequence family | weakest |

**The pocket is the transferable unit, not the family.** TNF-alpha and IL-17A are
both cytokines, both PPI targets, both drugged with antibodies first — and their
small-molecule stories share nothing mechanically. TNF-alpha's site is a cavity
on the trimer 3-fold axis, opened by displacing a subunit. IL-17A's is a groove
at the homodimer interface, addressed by macrocycles from 2016. A jump along
"same cytokine family" transfers nothing.

Two consequences for a consumer:

- **When the axes disagree, that is the informative thing on the page, not noise
  to reconcile.** Measured: IL-17A has a strong `target_precedent` — 44
  structures, 20 holo, a real macrocycle series — and an empty
  `structural_neighbour_precedent`, on the same target. Both KRAS and IL-17A
  landed **0 of 25** defensible small-molecule holo fold-neighbours. That is a
  real and reportable finding, not a retrieval failure.
- **Everything in `pocket_neighbour_precedent` is a hypothesis, not a
  measurement.** It is labelled transferred, names its source target, and carries
  the similarity value so you can discount it. Its `cofold_transfer` sub-block —
  the sharp falsifiable test — is `null` throughout this deployment (reading rule
  4), so what remains is a descriptor similarity with no confirmation. Weight it
  accordingly. "No actives on this target; 340 actives across the Pfam family,
  best 2 nM" is an honest and useful statement. "Moderate precedent" is not.

## Reading rule 4 — `null` means not measured. `0` means measured zero. Read `not_found`.

`null` is never a stand-in for a low or absent value. It means the measurement
was not made, and the reason is in `not_found`.

The validator makes this checkable. `check_null_is_not_zero` walks 15 named
fields (`target_precedent.distinct_actives`, `best_potency_nm`,
`assay_concentration.top_assay_share_pct`, `structure.total_pdb_structures`,
`structure.holo_count`, the volume and druggability min/max,
`pocket_hydrophobic_density`, `disorder_fraction`,
`annotated_binding_site_overlap`, `ligand_site_jaccard`,
`max_backbone_ca_displacement_a`, `ensemble_consensus_fraction.fraction_with_strong_pocket`)
and requires, for each:

- a `null` to have a matching `not_found` entry naming it, **or the dossier
  fails**;
- a field named in `not_found` **not** to also be reported as `0` — a measured
  zero is a result, and you must pick one;
- no placeholder strings (`"n/a"`, `"unknown"`, `"none"`, `"-"`) in numeric
  fields.

Matching is deliberately tight: naming the enclosing *block* in `not_found` does
not excuse a null leaf, because one vague line would otherwise launder every null
in the block.

**So: `not_found` is a load-bearing field, not an appendix.** Its entries are
`{"field": "<dotted path or leaf>", "reason": "<why>"}`. Parse it. A consumer
that ignores it will read deliberate abstentions as negative findings.

### Whole axes are nulled by rule, not by failure

Four capabilities have **no tool in this deployment**. `CLAUDE.md` rule 13
requires them nulled with a stated reason, and `rubric.md` criterion 21 fails any
run that populates them — a recalled number is indistinguishable from a measured
one once it is in the JSON.

| field | why it is null | expect |
| --- | --- | --- |
| the entire `affinity` block, including rule 12's mandatory positive control | no affinity predictor | all `null`, `reliable: null`, `predictions: []` |
| `structure.cofold_control` | no cofolding model | all `null` |
| `pocket_neighbour_precedent.*.cofold_transfer` | no cofolding model | all `null` |
| `structure.tier` values `cofolded`, `predicted`, `sampled_ensemble` | no structure predictor | unreachable. Only `holo_experimental`, `apo_experimental` or `none` will appear |
| the `Unknown`-modality cross-check (rule 10b) | no Open Targets client | the drug stays modality-unknown |
| `target_precedent.patents.count` | the patent source returns "Patents sources are not available." | `null` |
| `structural_neighbour_precedent` when the neighbour tool raised `ModuleNotFoundError` | dependency absent on the operator machine | nulled and recorded as **unavailable**, never as "no neighbours found" |

**Do not treat these nulls as signal.** `affinity: null` is not "poor predicted
affinity". `cofold_control: null` is not "the cofold failed to reproduce the
crystal". These are declared abstentions; both shipped dossiers carry them
identically, on a target with nine approved drugs and on a target with none.

The patents case is the sharpest illustration of why this matters: the patent
source returns the string `"Patents sources are not available."` **with exit code
0**. In a pipeline that is indistinguishable from "no patents found". It is
recorded as unavailable.

## Reading rule 5 — A dossier is an answer about a target *and* a mechanism. Key on both.

Chain selection is not a preparation preference. It is an assertion about which
interaction you intend to break, and it silently changes the number:

- KRAS 4OBE gives druggability **0.442 at rank 1 on chain A** and **0.257 at rank
  6 on chains A+B** — same structure, same clustering, different verdict.
- TNF-alpha prepared as one chain has no site at all. Its site **is** the trimer:
  0.00 Å³ intact at the SPD304 site, ~280–550 Å³ with a protomer deleted.
- A measured run in the deployment notes: 1TNF whole assembly gives 14 pockets
  and an ensemble volume of 153.3–198.3 Å³; the same call restricted to chains
  A and B gives 7 pockets and 155.9–156.0 Å³.

Four mechanisms, all real and all in the fixture set, requiring different chains:

| mechanism | example | where the pocket sits | chains needed |
| --- | --- | --- | --- |
| orthosteric | BCL-2 + venetoclax | in the BH3 groove — the epitope itself | the partner's contact chain |
| allosteric | TYK2 + deucravacitinib | JH2 pseudokinase domain — neither ATP site nor interface | the domain, by residue range |
| oligomer destabilisation | TNF-alpha + SPD304 | inside the trimer axis; displaces a subunit | **all** subunits |
| adjacent cryptic, state-locking | KRAS switch-II | beside the effector interface; locks the inactive state | the single chain |

A system that only inspects the annotated binding site or the PPI epitope misses
three of these four. TYK2's structures split **29 entries for JH1 and 21 for
JH2**; picking the wrong domain scores the wrong pocket on an approved-drug
target. IRAK4 has two separable functions in the same protein — kinase activity
and scaffolding the MyD88 signalosome — and a kinase inhibitor stops only the
first, so a dossier answering the kinase question says nothing about the
scaffolding question.

**Therefore: cache and join on `(uniprot_accession, mechanism_hypothesis)`, not
on the accession.** Two dossiers with the same accession and different mechanism
hypotheses are different answers, not duplicates, and reconciling them by taking
the more favourable one is a category error.

**Where the mechanism was not supplied**, both shipped examples show the honest
form: `tractability.method.chains_used` records the assembly actually scored, and
`tractability.caveat` states that no mechanism was specified and no pocket is
asserted to be the relevant one. `rubric.md` criterion 18 makes a null
`chains_used` legal **only** with that caveat present. So a dossier with a null
`chains_used` and a caveat is a whole-assembly answer; treat it as
under-specified, not as an answer about your mechanism.

The related field is `tractability.pocket_vs_interface.classification`
(`orthosteric_candidate` / `allosteric_candidate` / `destabiliser_candidate` /
`no_partner_structure`). It is a **measurement against a complex containing the
partner**, not a reading of the literature: the validator requires both
`partner_pdb_id` and `pocket_interface_overlap` before any substantive class may
be asserted. `no_partner_structure` is the honest abstention and is what both
shipped examples carry — including TNF-alpha, whose destabiliser mechanism is
known from the literature and is *still* not asserted, because no partner complex
was analysed.

## Reading rule 6 — Volume at D=1.6 is the computed number. Druggability carries nothing. Threshold neither.

This changed on 2026-08-15 on the strength of an evaluation over **15 targets, 67
structures and 134 measurements**, and it is the most likely thing for a
downstream consumer to get wrong, because the obvious field to threshold on is
the one that does not work.

### fpocket's druggability score does not separate druggable from hard

| measurement | value |
| --- | --- |
| target-level AUC at D=1.6 | **0.720**, bootstrap 95% CI **0.44–0.94** — the interval includes chance, P(AUC ≤ 0.5) = 0.071 |
| target-level AUC at D=2.4 | **0.520** — chance |
| label-free test, 37 holo structures with a drug-like ligand physically bound and the pocket anchored to it | median **0.320**; **25 of 37** below 0.5; **15 of 37 (41%)** below **0.1** |

Named cases, because a rate is easy to discount: **EGFR 6LUD with osimertinib
bound scores 0.013.** JAK1's median across nine approved drugs is 0.009. RORgt
6C1P is 0.009 at rank 55 of 60 — **but strike that one: the same audit that
retracted the volume separation found 6C1P contains no RORgt (sole entity
A8EVM5, an ion transport protein) and that its anchor ligand `1N7` is CHAPSO, a
detergent, so it is neither the target nor a drug-anchored positive.** The
demotion does not depend on it; the remaining 36 holo structures carry it. TYK2
6NZP with deucravacitinib is 0.169. BCL-2
6QGK is 0.025. And the inversion is confirmed at target level: **MYC** — zero
holo structures, canonical undruggable — has a D=2.4 median of **0.75**, above
KRAS (0.54), BCL-2 (0.52), JAK1 (0.49), EGFR (0.44) and NLRP3 (0.12).

There is a mechanistic reason independent of the statistics: fpocket's
`drug_score_pocket` leans on `mean_loc_hyd_dens_norm`, which is **min-max
normalised across the other pockets of the same structure**. Druggability is a
property of a pocket *relative to the population detected beside it*, not of the
pocket. And the shipped score is a three-descriptor logistic regression fitted on
**21 positive pockets against 292 others** — the published 2010 nested model is
commented out in the source.

**Consequences you must honour:**

- `tractability.pocket_druggability.load_bearing` is fixed at literal `false`.
  The validator (`DRUGGABILITY_LOAD_BEARING`) accepts no other value — not
  `true`, not `"false"`, not `0`, not missing.
- `tractability.pocket_druggability._false_negative_rate` carries the measured
  rate as a string beside every reported range, so a reader meeting a 0.02 later
  can discount it.
- Druggability may **not** carry a `not_tractable` or `insufficient_evidence`
  verdict on its own, in any `verdict_basis`. A negative computed verdict needs
  the D=1.6 volume behind it; with no volume, the honest output is
  `insufficient_evidence` **with the unmeasured volume named as the reason**, not
  a poor pocket.
- It is reported as a **range** (`min`/`max`/`fold_range`), never a point. A
  one-sided range fails. A loose druggability number anywhere else in the dossier
  fails.

**Do not threshold on `pocket_druggability`. Do not rank on it. Do not display it
without `_false_negative_rate` beside it.**

### Volume is the primary number — and the separation behind it is RETRACTED

`tractability.pocket_volume_a3.primary_d1_6_a3` is the computed axis's number.
**A volume is a measurement of a cavity in a structure we scored. It carries no
verdict, and there is no volume at which this output calls a target druggable or
hard. Do not build a threshold, a gate, a sort order or a colour scale on it.**

**RETRACTED 2026-08-15 — read this before you wire anything up.** This section
previously stated, as a current result, that volume at D=1.6 separated all 15
evaluation targets perfectly: AUC 1.000, CI [1.000, 1.000], stable under all 15
leave-one-target-out refits, every known-hard target at or below 207 Å³ and every
known-druggable one at or above 242 Å³. **That claim is withdrawn in full.** It
is recorded here as a retraction rather than deleted, because a consumer who
meets the number in an older dossier, an older copy of these notes or a cached
summary needs to be able to find out what happened to it. The upstream
retractions are `CLAUDE.md` rule 4a, `rubric.md` and `pipeline.html` N2.

**Why it failed — the anchors.** A full audit resolved every lining residue of
all 67 calibration structures to an entity, by aligning SEQRES to the target's
UniProt sequence directly. **Four of the five hard anchors are compromised, and
two druggable ones with them:**

| anchor | what the audit found |
| --- | --- |
| **MYC 188 Å³** | 0 of 6 lining residues are MYC. The pocket is **100% MAX (P61244)** — a different protein |
| **IL-11 164 Å³** | 0 of 7 lining residues are IL-11. The pocket is **100% IL-11 receptor alpha (Q14626)** |
| **TNF-alpha 207 Å³** | on TNF, but **zero residue overlap** with the only genuinely drug-anchored pocket. TNF's defensible value is **129.6 Å³** |
| **CD20 154 Å³** | on CD20, but the anchor ligand `Y01` is **cholesterol hemisuccinate — a detergent.** A lipid site on a membrane protein |
| **TL1A 137 Å³** | the trimer 3-fold axis, with **no site anchor at all** |
| **KRAS 400 Å³** | a median over two different sites — and **there is no switch-II pocket in 4OBE at all** (no pocket carries three or more switch-II residues), because switch-II is closed in the GDP state. That is the cryptic-pocket story of rule 3, not a measurement |

**And one that breaks the assumption the fix rested on. RORgt's 6C1P contains no
RORgt.** Its sole entity is **A8EVM5, "Ion transport protein"** — a NavAb
sodium-channel HypoPP mutant, and the entry title says so. All eight lining
residues are the channel. It was selected by **`ligand_site_jaccard`**, the
selection path these notes call trustworthy, and its anchor ligand `1N7` is
**CHAPSO, a detergent**. So restricting scoring to the target's chains is
**necessary but not sufficient**: a wrong PDB ID passes straight through the
trusted path. A structure list needs accession verification of its own.

**The mechanism, confirmed in the artifacts.** `chain_accessions` is `{}` on
**every single entry**, with `target_chains_basis: "entry declares no _struct_ref
UniProt mapping"` — while the `_why` string sitting next to it asserts that
chains "are resolved by UniProt accession from the entry's own
`_struct_ref_seq`". Resolution never once succeeded, and every chain of every
assembly was scored as though it were the target.
`max_druggability_no_ligand_site` — "the most druggable pocket anywhere in the
assembly", which identifies no site — sets the headline median for **MYC, IL-11,
TL1A and TNF outright, and KRAS by half.**

**The two statistical points matter more than the anchors, and they generalise:**

1. **A bootstrap CI on a perfectly separated set is degenerate by
   construction.** Resampling cannot create an inversion that is not in the data,
   so the quoted `[1.000, 1.000]` was arithmetic, not evidence. It carried no
   information at any point, and it read as the strongest line in the table.
2. **The confound is fatal and it is simple.** The binary flag "a drug-like
   ligand was co-crystallised" separates the two groups at **AUC 0.900 using no
   structural measurement whatsoever.** The label and the measurability are the
   same variable. So the honest conclusion is not that the anchors were wrong:

   > **If the hard side can only ever be measured by "whichever pocket ranked
   > highest", then this axis is measuring structure availability, not biology.**

   Which is why a properly anchored calibration set is close to a contradiction
   in terms. It would need hard targets carrying a real drug-like co-crystal at a
   defined site, and a hard target is largely defined by not having one.

**Corrected values, for the record only — not a new boundary.** MYC 242.0
(exactly NLRP3's 242, a margin of **+2.1 Å³**, and 151.8–325.7 across
anchorings), IL-11 146.7, RORgt 428.5 with 6C1P removed, KRAS 597.9, TNF 129.6,
IRAK4 593.1, JAK1 257.7, NLRP3 244.1. **CD20 and TL1A have no valid value at all;
MYC and IL-11 have no site-anchored value.** Under one anchoring AUC is 1.000
with a +2.1 Å³ margin; under another it is **0.900, CI [0.680, 1.000]**, with a
margin of **−81.6 Å³**. The exact permutation p for AUC = 1.0 at n=10 druggable
against 1–2 hard is **0.015 to 0.091** — not significant. **The test has no
resolution left**, and re-running it on the corrected numbers would not restore
one.

The evaluation's own stated limitation also still travels with any future
attempt: n = 5 hard targets, all PPI / cytokine / membrane class, against a
druggable set enriched in kinases, nuclear receptors and GPCRs — **volume may
partly be tracking target class rather than tractability.**

So: **a target at 230 Å³ is not classified by this. It is unclassified by it —
and so is a target at 600 Å³ and a target at 130 Å³.**
`rubric.md` names "do not gate on a volume threshold" as a failure mode of the
grader itself, and the only numeric volume gates anywhere in the system are
structural, not classificatory:

- **1000 Å³** — a `primary_d1_6_a3` above this is the signature of sites merged
  with neighbouring cavities, not a big pocket. It fails validation.
- **240 Å³** — a low druggability sitting beside a volume at or above this
  *requires* a non-empty `tractability.caveat` stating the disagreement. It
  classifies nothing; it forces the disagreement to be written down.

Read `primary_d1_6_a3` specifically. It is **not** the min of the pooled
`min`/`max` spread: at D=2.4 volumes exceed 1000 Å³ and sites merge, so a spread
pooled over both clustering values is a spread over two different things.
`pocket_volume_a3.clustering_d` records which values were pooled.

### Two more computed fields that are not quality values

- **`tractability.site_pocket_rank.prank`** is a site-*finding* aid, reported
  beside `site_pocket_rank.fpocket`, never replacing it. On n = 70
  ligand-anchored measurements PRANK promotes the true site in 79% of cases and
  demotes it in 1% (median rank 5 → 1; top-3 recall 37% → 91%). **As a
  druggability classifier its rank is inverted, AUC 0.25** — worse than chance in
  the systematic direction, because on a target with no ligand to anchor to the
  top-ranked pocket is top-ranked by construction. Never read rank as quality.
- **`tractability.ensemble_consensus_fraction`** is an anti-cherry-picking
  control, not a tractability signal. Persistence is constant — the site pocket
  was detected in 100% of structures for all 15 targets, **AUC exactly 0.500** —
  and the published consensus criterion built on it gives **AUC 0.560 and ranks
  MYC first at 0.80**, above 8 of the 10 druggable targets. It stops the agent
  quoting its best conformer. It tells you nothing about the site. Both shipped
  dossiers report the denominators (`n_structures`, `n_measurements`) and leave
  `fraction_with_strong_pocket` and `meets_consensus_criterion` null, because the
  tool emits no consensus fraction.

## Reading rule 7 — `ASK[...]` in `not_found` is a question we sent upstream. It never blocks anything.

The station can ask questions back to the upstream literature knowledge graph.
When it has, a pending ask appears as a `not_found` entry whose `reason` is
prefixed:

```json
{
  "field": "target_precedent.clinical_stage_small_molecules",
  "reason": "ASK[resolve_link:L3] issued to graph g_tl1a1 round 3. Not blocking: the field is filled from what we could verify and this entry records the residual question. <question text>"
}
```

The four verbs are the upstream graph's own: **`expand_node`**,
**`resolve_link`**, **`test_gap`**, **`new_question`**. The target after the
colon is a graph row id.

**The rule is absolute and you may rely on it:**

> An outstanding ask never blocks a verdict and never licenses a null.

The dossier is completed as if the ask will never be answered. Every field the
ask would have improved is either filled from what was verifiable, or nulled with
a reason that stands on its own — a reason that would still be true if the ask
had never been written. If the only reason a field is null were a pending ask,
the ask would be illegitimate. And `verdict: insufficient_evidence` may never be
reached by counting a pending ask as missing evidence.

An ask only fires when all five gates pass, the strictest being gate 3: the local
lookups (ChEMBL, the trials registry, the structure, a grep on the exact
identifiers) were **actually run**, **all abstained**, and each null is already
in `not_found` before the ask is drafted. The mechanism is designed to fire
rarely.

**Delivery.** The station does not write to the upstream graph — not to its
`rounds` array, not anywhere. `graph-intake` emits the ask object into a
nomination's `asks: []` and stops; that object is not part of the dossier
template. So **the `ASK[...]` string in `not_found` is the only carrier of a
pending ask inside the JSON you receive**, and getting it to the upstream team is
an out-of-band, human step today. If your station wants to route asks, scan
`not_found[].reason` for the `ASK[` prefix; nothing else will surface them.

Treat an `ASK[...]` entry as an open question with a named addressee. It is not a
degraded measurement, and it is not a reason to discount the verdict beside it.

---

## What the machine grader guarantees

`validate_dossier.py` carries **16 rule functions emitting 17 violation codes**.
It is the machine grader; what it enforces *is* the promise we make. It checks
whether a dossier is *honest about what it knows*, not whether it is *right* — "a
dossier that passes is not necessarily correct; a dossier that fails is
definitely wrong in the way the violation names."

| code | promise to you |
| --- | --- |
| `WELL_FORMED` | All 16 top-level keys present. Every enum non-null and inside its legal set. No template placeholder strings (`"a \| b"`) left standing. No nameless stub entries in the six list blocks — `[]` means none. `verdict` is a label, never a number. `next_experiment.description` and `biologic_precedent.note` non-empty. `falsification.survived` is a real boolean and `checks_run` non-empty. No NaN or Infinity |
| `NUMBER_WITHOUT_PROVENANCE` | **Every numeric leaf sits inside a block naming a source.** Provenance inherits downward only |
| `AXES_AVERAGED` | No combined, composite or overall score exists at any depth. No number mixes a precedent term with a tractability term |
| `MODALITY_LEAK` | Nothing appears in both the small-molecule and biologic blocks. Every small-molecule entry declares `modality: "small_molecule"`. No `-mab`/`-cept` stem in a small-molecule list. A tractable verdict cannot rest on biologic precedent |
| `INSUFFICIENT_EVIDENCE_AVOIDED` | A thin target (< 50 actives, `holo_count: 0`, no approvals) cannot be given a confident verdict. A declination must name what would resolve it |
| `DRUGGABILITY_POINT_ESTIMATE` | Druggability and volume are objects with ranges, never scalars. A range requires ≥ 2 clustering values swept and a named ensemble. No loose druggability figure elsewhere |
| `DRUGGABILITY_LOAD_BEARING` | `load_bearing` is literally `false`. The false-negative rate travels with the number. **No negative verdict rests on druggability alone.** A low score beside a ≥ 240 Å³ volume forces a written caveat |
| `VOLUME_NOT_PRIMARY` | `primary_d1_6_a3` is present whenever any geometry is reported, or its absence is named in `not_found`. It never exceeds 1000 Å³. A spread records which clustering values it pools |
| `FRACTION_WITHOUT_N` | Every fraction carries its denominator, and a stated `n_structures` matches the named ensemble |
| `SAME_SITE_BASIS_MISSING` | Any pooled spread records **how the site was chosen**, from a closed list of five bases |
| `SAME_SITE_BASIS_INVALID` | The three bases that do not identify a site are never pooled across measurements |
| `SITE_INCONSISTENT` | Geometry is not quoted off a site definition measured to be a different pocket. An interface classification is measured against a named partner structure, not assumed |
| `CRYPTIC_MISCLAIM` | A cryptic call carries its apo census. `is_cryptic: true` requires the site absent in ≥ 80% of examined apo structures (Vajda 2018's "all or nearly all"). Site present in the apo ensemble means **occluded, not cryptic**. Mechanism and potency prior must agree |
| `NULL_IS_NOT_ZERO` | Every null in the 15 measured fields says why, in `not_found`. A measured zero is never also listed as missing. No placeholder strings in numeric fields |
| `AS_OF_LEAKAGE` | Under a cutoff, undated sources carry explicit `leakage_risk: true`, and no year or release date exceeds the cutoff |
| `AXIS_CONFLICT_UNDECLARED` | Disagreement between the axes is **stated**, not resolved. Triggers on biologics-without-small-molecules, ≥ 500 actives with zero holo structures, one assay at ≥ 30% share measuring a different protein, or an uncharacterised headline potency |
| `ASSAY_PROVENANCE_MISSING` | A non-zero actives count names its top contributing assay and that assay's share, and answers whether it measures a different protein. A reported potency says whether its assay is characterised |

`rubric.md` criteria 18–21 are file-verifiable but **not** validator-enforced:
chain selection recorded (18), the method block reconstructing the run (19), site
rank reporting both rankers or neither (20), and unavailable axes nulled rather
than fabricated (21). Criterion 21's *hard* half — a populated `affinity.*`,
`structure.cofold_control` or `cofold_transfer` — fails a run outright.

**One gap you should know about.** The validator walks structured fields. Prose
written into `axis_conflict` or `tractability.caveat` is a string and passes
untouched. `falsification-sweep` names this explicitly: writing "the pocket
scores 0.02 and is probably not real" into `falsification.findings` as prose
routes around `DRUGGABILITY_LOAD_BEARING` entirely. So the structured fields are
guaranteed; the prose is not. If you are extracting signal, extract it from the
structured fields.

---

## Provenance: how to follow any number

**Every numeric claim carries a source.** A figure without provenance cannot
appear in the dossier — `check_number_provenance` walks every numeric leaf and
fails any that does not sit inside a dict holding a non-empty provenance key.

The 24 accepted keys are: `source`, `sources`, `_provenance`, `provenance`,
`doi`, `pubmed_id`, `pmid`, `url`, `query`, `queries`, `chembl_target_id`,
`chembl_id`, `chembl`, `assay_id`, `best_potency_assay`, `pdb_id`,
`source_pdb_id`, `reference_pdb_id`, `partner_pdb_id`, `query_structure`,
`ensemble_pdb_ids`, `tool`, `basis`, `descriptor_basis`.

**Provenance inherits downward and only downward.** A `sources` list on a block
covers every number inside it; a source on one drug entry covers nothing in a
sibling block. An empty `sources` list attributes nothing. Four blocks hold
numbers no other key attributes and therefore always carry their own `sources`:
`target`, `tractability`, `structure`, `affinity`.

Sources take four forms, and each is followable:

| form | example, from `jak1_P23458.json` | how to follow it |
| --- | --- | --- |
| database identifier | `"CHEMBL2835 (chembl_v.bioactivities_by_accession, compounds_by_accession, drugs_by_accession)"` | query ChEMBL on that target id |
| structure identifier | `"3EYG (earliest JAK1 holo, released 2009-02-03)"` | RCSB entry |
| tool + version + parameters | `"pocket-scan Modal function `pocket_scan`, measured run 2026-08-15 on 3EYG + 10PI, D in {1.6, 2.4}"`, with `"fpocket 4.2.3 (conda-forge), P2Rank 2.5.1 (rescore_2024), mdpocket"` | re-run with the recorded ensemble and sweep |
| literature citation | `"Lazou, Kozakov, Joseph-McCarthy and Vajda, Drug Discov Today 2024"` | the paper |

`tractability.method` is designed to reconstruct the run: `tool`, `version`,
`clustering_d_swept`, `ensemble_pdb_ids`, `chains_used`, plus a prose note on
preparation (protein only, altloc A/blank, hydrogens stripped — so ligands are
stripped before scoring and the holo/apo scale question does not arise).

**Keys beginning with `_` are notes, not claims.** They are invisible to the
provenance walker by design and carry the reason a number is discounted.
`_false_negative_rate`, `_primary_note`, `_provenance` and `_note` are all of
this kind. Display them next to the numbers they qualify; do not parse them as
data.

---

## Worked example: JAK1 (P23458), field by field

`.claude/skills/assemble-dossier/examples/jak1_P23458.json`. Zero violations.
This is the shape of a dossier where the retrieved axis is abundant and the
computed axis agrees.

**Header.**

```json
"verdict": "small_molecule_tractable",
"verdict_basis": "retrieved_precedent",
"axis_conflict": null
```

Read this as: tractable, and **the retrieved axis carried it**. The computed axis
agrees but is not what the verdict rests on — the caveat says so explicitly:
"JAK1's verdict is carried by retrieved precedent — nine approved small molecules
— and would be unchanged if the pocket run had never happened." `axis_conflict`
is null because the axes agree. `as_of_date` is null, so `as_of_leakage` is `[]`
and no leakage machinery applies.

**Axis 1 — `target_precedent`.**

```json
"chembl_target_id": "CHEMBL2835",
"distinct_actives": 14472,
"assay_concentration": {
  "top_assay_description": "JAK1 enzyme inhibition (characterised); no single assay dominates",
  "top_assay_share_pct": 4.8,
  "measures_a_different_target": false,
  "assay_type_split": {"binding_B": null, "functional_F": null}
},
"best_potency_nm": 0.01,
"best_potency_characterised": true,
"approved_small_molecules_count": 9,
```

14,472 actives is a real number *because* `top_assay_share_pct` is 4.8% — no
single assay dominates. **Never read `distinct_actives` without reading
`top_assay_share_pct` beside it**; see the TNF contrast below. `best_potency_nm`
0.01 is usable because `best_potency_characterised: true` — it comes from a named
enzyme assay, not from a description reading "inhibition assay using X".

`approved_small_molecules_count` is 9 and the list names 8. The ninth is in
`not_found` — counted after collapsing salt/parent `molregno` pairs, unnameable
from the retrieved rows, deliberately not guessed. Every named entry carries
`"modality": "small_molecule"` and its own `source`.

`assay_type_split` is `{null, null}` with a `not_found` entry — and that entry
tells you the split would not have helped anyway, since TNF-alpha's contaminating
assay sits in the binding bucket.

`biologic_precedent.approved_biologics` is `[]`. `family_precedent` is all null
with a `not_found` entry saying it was not retrieved and would not change the
verdict. **Null, not zero.**

**`structure`.**

```json
"tier": "holo_experimental",
"pdb_id": "3EYG",
"bound_ligand": {"comp_id": "MI1", "heavy_atoms": 23, "is_druglike": true,
                 "is_known_frequent_hitter": null},
"total_pdb_structures": 52, "holo_count": 42, "apo_count": 10,
"ensemble_used": ["3EYG", "10PI"],
"cofold_control": { ... all null ... }
```

Tier `holo_experimental` is the top of the selection order (holo > apo >
predicted). 42 of 52 entries carry a drug-like ligand. `is_known_frequent_hitter`
is `null` — no PAINS screen was run, so the field is **not** assumed clean. The
whole `cofold_control` block is null by rule 13, not by failure.

**Axis 2 — `tractability`.**

```json
"pocket_volume_a3": {"min": 305.9, "max": 913.8, "spread_pct": 66.5,
                     "clustering_d": [1.6, 2.4],
                     "primary_d1_6_a3": null,
                     "site_pocket_selected_by": "ligand_site_jaccard"},
"pocket_druggability": {"min": 0.02, "max": 0.437, "fold_range": 21.8,
                        "site_pocket_selected_by": "ligand_site_jaccard",
                        "load_bearing": false}
```

Read in this order:

1. **`site_pocket_selected_by`** first. `ligand_site_jaccard` is the strongest of
   the five bases — it means every pooled measurement had real overlap with the
   MI1 ligand site (0.300 to 0.778 per measurement). This is what licenses the
   pool. If it read `site_signature_unreliable_homooligomer`,
   `max_druggability_no_ligand_site` or `no_pocket_matched_site_signature`, the
   values would not be one site and could not be pooled at all.
2. **`primary_d1_6_a3`** — here it is **`null`**, with a long `not_found` entry:
   the run pooled four measurements without carrying the per-measurement values
   through, so the D=1.6 figure cannot be recovered. That entry ends "This is a
   MISSING PRIMARY MEASUREMENT, not a finding that the JAK1 ATP site is small."
   **This is exactly the null you must not read as a low value.**
3. **`min`/`max` 305.9–913.8 Å³** is pooled over both clustering values, so it is
   not the primary number and must not be substituted for one. The 913.8 upper
   end is close enough to the ~1000 Å³ merge signature to be worth watching.
4. **`pocket_druggability` 0.020–0.437** with `load_bearing: false`. Its maximum
   sits below the 0.5 band in which 25 of 37 ligand-anchored holo pockets fall.
   Nothing turns on it.

`tractability.caveat` is the field that reconciles the two: it records that the
low druggability and the large volume disagree, that the disagreement is
**reported and not resolved**, and that the volume guide it would be resolved
against gates nothing. **Read `caveat` on every dossier.** It is prose, it is
unvalidated, and it is where the run's own reservations live.

The cryptic block:

```json
"cryptic_pocket_risk": "low",
"cryptic_mechanism": "none",
"cryptic_evidence": {"is_cryptic": false, "n_apo_examined": 10,
                     "n_apo_site_absent": null,
                     "site_present_in_apo_ensemble": null,
                     "basis": "the ATP site is pre-formed and occupied ...",
                     "definition": "Vajda et al. 2018 ..."}
```

`is_cryptic: false` is the call; `basis` states honestly that no apo/holo
superposition was run (the cryptic stage returned `not_run` on an all-holo
ensemble) so the call rests on the holo census rather than a displacement
measurement. `max_backbone_ca_displacement_a` is `null` with that reason in
`not_found`. **A null displacement here means the comparison was impossible, not
that the site does not move.**

`site_centroid_to_ligand_distance_a: 1.86` is the rule-4b off-site check
**passing**: mdpocket's ligand-free site definition sits 1.86 Å from the ligand
centroid, well inside the proposed (not calibrated) 4 Å threshold, so the two
independent site definitions agree. On TNF-alpha the same field reads 29.57 Å.
`mdpocket_site_definition_used` is `site_from_ligand`, which is the ligand site by
construction.

**`falsification`.** Nine `checks_run`, four `findings`, `survived: true`. Note
that `checks_run` lists checks that found nothing — the terminated-programs check
records "none retrieved for JAK1 in this run. Absence of retrieved terminations
is not evidence that none exist — the query was not run against a trials
registry." **That is the correct reading of an empty `terminated_programs`
list anywhere in this output: silence, not a clean record.**

**`next_experiment`** names the one gap that would close: add an apo entry to the
ensemble, which turns on the cryptic stage and produces a displacement
measurement instead of a census inference.

**`not_found`** has 15 entries. Every null above is in there with a reason, and
several reasons carry the measured values that could not be put in a scalar field
— e.g. `ligand_site_jaccard` was 0.560 and 0.778 for 3EYG and 0.300 and 0.455 for
10PI, and the field is null because the template has one slot and picking one
would be arbitrary. **Useful data lives in `not_found` reasons.** If you are
mining, read them.

### The contrast: TNF-alpha (P01375), where the computed axis is refused

Same tool, same day, opposite outcome. Four things to take from it:

1. **`axis_conflict` is populated and is the finding.** Zero approved small
   molecules against five approved biologics, with a Phase 2 oral small molecule
   (balinatunfib) settling the direction. `verdict_basis` is `both`.
2. **`distinct_actives: 2582` is not what it looks like.**
   `top_assay_share_pct: 45.0` and `measures_a_different_target: true` — the
   top assay is an "IRAK4 Monocyte TNFalpha Cell Based Assay" that measures a
   different protein and uses TNF only as a cellular readout. The `assay_type`
   split (B 5830 / F 617) is no defence: the contaminating assay is in the
   binding bucket. **This is why the validator requires those fields whenever
   `distinct_actives` is non-zero.**
3. **`best_potency_nm: 1.3`, not 0.03.** The headline 0.03 nM Ki sits behind an
   assay described only as "Inhibition assay using TNF-alpha" and was rejected as
   uncharacterised. The reported figure is the weaker, characterised SPR Kd —
   three orders of magnitude worse, and correct.
4. **`pocket_volume_a3` and `pocket_druggability` are entirely null, and that is a
   refusal, not a low score.** A six-entry, twelve-measurement run completed with
   every stage `ok`. Ten of twelve measurements came back
   `site_signature_unreliable_homooligomer` — a basis that by construction cannot
   tell one protomer's site from another's on a C3-symmetric trimer. Pooling them
   anyway yields a druggability fold-range of exactly **651.0**, which is the
   withdrawn "651-fold spread" regenerating from the identical defect. The
   `not_found` entry says so in words: "REFUSED, not missing."

**The generalisable lesson:** a null in `tractability` on a target with zero
approved small molecules is the exact place a consumer will infer "bad pocket".
The dossier's own caveat pre-empts it: "the computed axis is REFUSED here, not
weak. A refusal is not a low score, and it must not be read as one when it sits
beside TNF's zero approved small molecules."

---

## Failure and degradation: telling "not measured" from "measured negative"

This is the pipeline's own hard-won lesson, and it is the section to read before
writing any downstream conditional. **A timeout, a missing credential, an
unclassifiable ligand or a refused superposition must never read as a finding.**

### Where the distinction is recorded

| you see | it means | it does **not** mean |
| --- | --- | --- |
| a field is `null` **and** named in `not_found` | not measured, for the stated reason | a low or absent value |
| a field is `0` and **not** in `not_found` | measured zero. A result | a failed measurement |
| `not_found` reason says **"REFUSED, not missing"** | the tool ran and declined to attribute a number to a site | the number was bad |
| `structure.tier: "none"` | no usable experimental structure | no structure exists |
| `holo_count: 0` | measured zero holo entries. Triggers the insufficient-evidence rule | unknown |
| `holo_count: null` | the census was not run. Deliberately does **not** trigger that rule | zero |
| `terminated_programs: []` | no terminations were retrieved | none exist |
| `terminated_programs: null` | the sweep could not run | none exist |
| `approved_small_molecules: []` with a non-null count of 0 | **a real and common finding.** Validated target, zero approved small molecules | nothing was ever made — the drug view only lists drugs with a curated direct mechanism, so a target can carry thousands of bioactivities and no rows here |
| `pocket_volume_a3` values of `0.00` | a measured zero. mdpocket returned 0.00 at the true site rather than substituting a nearby pocket | a failed run |
| `druggability_status: "not_available"` quoted in a reason | the quantity is **undefined by construction** on a fixed grid | unavailable, missing, or something a future release might add |
| `falsification.survived: false` | a check undercut the *precedent claim* | the target is not tractable |
| `affinity.*` all null | no predictor exists in this deployment | poor predicted affinity |
| `structural_neighbour_precedent` null with an unavailability reason | the neighbour tool raised `ModuleNotFoundError` | no structural neighbours found |
| `patents.count: null` | the patent source is unavailable and returns a message with **exit 0** | no patents |

### Stage status strings you may see quoted in `sources`, `basis` or `not_found`

The pocket-scan tool emits a per-stage `<stage>_status`. Its values are `ok`,
`not_run` and `failed`, and the dossier quotes them verbatim.

- **`cryptic_status: "not_run"`** — the comparison was impossible. JAK1's
  `not_found` quotes it: "needs both an apo and a holo structure; got 1 holo
  reference and 0 apo. This comparison is a PAIRWISE measurement — there is
  nothing to superpose against with only one state." Not a finding that the site
  is rigid.
- **`cryptic_status: "failed"`** with a `superposition_gate` block — the
  superposition did not meet all three gates (core RMSD ≤ 5.0 Å, ≥ 20 equivalent
  Cα, ≤ 10% residue-name mismatches). This gate exists because three targets
  produced confident, wrong mechanistic calls on broken fits: NLRP3 reported
  `is_cryptic: true` at 21.6 Å on a fit whose own core RMSD was 16.6 Å (correct
  answer: not cryptic, 0.95 Å); S1PR1 reported `subunit_occlusion` — which maps
  to a **micromolar-at-best** potency ceiling — after mapping the receptor onto a
  25-residue peptide, on a target with 600 sub-nanomolar compounds and five
  approved drugs. **Four log units wrong, in the direction that kills a
  program.** A `failed` status is the gate working.
- **`druggability_status: "not_available"`** — mdpocket cannot report a
  druggability at all, by design: fpocket's score is normalised across the other
  pockets of the same structure and a fixed grid has a population of one, so the
  quantity is undefined. **A null here is correct; a number would be the bug.**
- **`off_site_warning`** with `distance_to_donor_ligand_centroid_a` — the site
  definition being quoted is a different pocket. TNF-alpha's reads 29.57 Å
  against a proposed 4 Å threshold. Geometry quoted off that definition is
  refused, and `site_hypothesis_basis` becomes `not_established`.
- **`buried_core_suspected`** — a near-sealed hydrophobic pocket that fpocket's
  regression rewards. IRAK4's death domain gave the top-ranked pocket of 134 at
  druggability 0.890, and it is the hydrophobic core of the domain, not a site.
  Read the flag as "this druggability value is uninterpretable", not as "this
  pocket is not there".
- **`LigandSourceError`** — a run-killing error, never a result. See the next
  subsection; it is the guard against the worst failure shape in the pipeline.
- **`frames_dropped`** and `n_processed < n_submitted_to_mdpocket` — structures
  fell out of the ensemble. Silent frame dropping inflates a persistence
  frequency, so **the failure looks like a stronger result**. The run refuses
  below three surviving structures once a drop occurred.
- **`filter.auto_relaxed`** on the structural-neighbour axis — the alignment
  floor was lowered because too few neighbours passed the verified one. The
  relaxed floor is a judgement call, not a calibrated threshold. A relaxed run is
  never reported without that block.

### A stage can be unavailable while every status field reads `ok`

This is the most dangerous pattern in the output and it has two documented
instances.

Ligands are classified from their **SMILES graph**. Hand the classifier a record
source that carries no SMILES — the entry's own mmCIF `_chem_comp` block is the
obvious one to reach for and does not carry it — and every component classifies
`unknown`, nothing is `druglike`, and **the whole ensemble comes back holo-free
while the payload is well-formed and every `<stage>_status` says `ok`.** The
`LigandSourceError` guard exists solely so that state cannot ship silently; if you
ever see it quoted, it is a misconfiguration, never a result.

The second instance: `files.rcsb.org/header/<ID>.cif` returns HTTP 200 and about
100 kB of plausible mmCIF that the parser rejects (the coordinate-loop deletion
leaves bare `loop_` keywords with no tags). The exception handler returned `[]`
and the caller fell back to longest-chain selection without a word.

The rule the pipeline drew from both: **assert on the value you came for, not on
the transport.** A 200 is not a parse, a parse is not a populated field, and a
populated field is not the field you needed. Apply the same rule to the dossier:
a populated `tractability` block is not evidence that a site was identified —
read `site_pocket_selected_by` and `site_hypothesis_basis`.

Three ligand-classification states you may see reflected, and they are not the
same:

| state | meaning |
| --- | --- |
| `druglike` | a drug-like ligand is bound. The only state that is evidence of a bindable site |
| `unknown` | classification could not be made. **`unknown` is not `apo`** |
| `lookup_failed` | the chemical-component lookup timed out. The entry is `undetermined`, not apo |

Reporting an `unknown` or `lookup_failed` entry as apo reintroduces the original
bug in a new place. Classification accuracy where it *is* determined: 259/262
(98.9%) on ground truth, 61/70 (87.1%) blind, **zero false positives**.

### Fields whose names invite the wrong reading

These are template fields you will consume directly. Each has a documented
mismatch between name and content.

| field | the trap |
| --- | --- |
| `structural_neighbour_precedent.neighbours[].evalue` | The underlying Foldseek wrapper's `hit.evalue` holds a **probability** (higher is better, 1.000 → 0.045) and `hit.bit_score` holds the true E-value. The skill's own output rule is to carry the **TM-score**, never the raw `evalue`. **Do not sort or threshold on `evalue`.** List order is best-first and is the safe ranking; use `tm_score` |
| `structural_neighbour_precedent.neighbours[].has_druglike_holo` | **Entry-level, not chain-level.** A drug-like ligand anywhere in the entry sets it, including one bound to a different protein in the same entry. Measured on a 137-entry IL-17A neighbourhood: **8 entries flagged and all 8 are false** — buffer components (bis-tris propane, Jeffamine), a diacylglycerol, a nucleotide analog. Read `ligand` and the entry title before believing it |
| `structure.holo_count` | Entry-level holo counts are an **upper bound**. RAN (P62826) reads 36 of 139 holo at entry level and **0** when attribution is unambiguous — leptomycin binds exportin, not RAN |
| `tractability.disorder_fraction` | A `0.0` here would read as "no disorder" when it is usually the **crystallised construct**, not the protein. IRAK4 measured 0.0 over 284 construct residues against 0.1413 over the full 460. Both shipped dossiers null this field for exactly this reason. Never quote a disorder number without knowing whether an accession was supplied |
| `tractability.pocket_druggability` | Discussed at length in reading rule 6. Its provenance string is in the block; read it |

### Credentials and binaries

`preflight()` runs **before the session is created** and aggregates every missing
credential and binary into one throw: `ANTHROPIC_API_KEY`, `PAPERCLIP_API_KEY`,
the `paperclip` / `micromamba` / `modal` binaries, the Modal profile, and a live
gemmi + numpy import. A failed preflight means **no dossier at all**, which is the
intended behaviour — you get an error, not a thin dossier.

The reason it is built this way is the one that matters to you: "A missing
credential must never look like a negative result. A `paperclip` call with no
valid key returning zero rows is indistinguishable from *a target with no
precedent* — and telling those two apart is this agent's entire job." A key that
is present but expired is caught separately, by checking failed runs for auth
signatures and converting them into a throw that says, in words, that this is an
authentication failure and not an empty result.

**So: if you receive a dossier at all, its zeros are measured zeros.** The
degradation modes above show up as nulls with reasons, never as empty results.

### Other degradations documented in the sources

- **Search-mode retrieval.** A semantic `search` returns confident hits that never
  name the compound. A "not found" produced by `search` is a retrieval failure,
  not an absence; identifier lookups use `grep`.
- **Grep match caps.** `"hit the per-shard match cap — more matches exist; raise
  -m N"`. A conclusion of absence before raising `-m` is invalid.
- **`why_stopped` truncation.** The trials registry caps it at 252 characters
  with no ellipsis flag; anything above ~240 characters should be assumed
  truncated.
- **`why_stopped` NULL on a `Terminated` row** is the most misleading state in the
  registry, not a clean record. And a `Withdrawn` status is not a failure —
  NCT06061523 was withdrawn because regulators agreed existing data sufficed.
- **Absence of a registry record is never evidence of absence of a trial.**
  VTP-43742's program-killing Phase 2a is not in the registry at all.
- **Non-US registries** (`ct_cn`, `ct_jpn`, `ct_eu`, `ct_global`) have no
  `why_stopped` column. They are registry-silent on *why*, which is not the same
  as a program stopping for no stated reason.
- **A pooled volume above ~1000 Å³** means sites have merged, and the
  druggability beside it is a merge artifact.
- **`max_radius_difference_a`** is the frame-independent same-site control.
  JAK1's is 2.29 Å; TNF-alpha's is 14.43–16.49 Å; the IL-17A pool that was ruled
  invalid was 16.61 Å. A large value means the pooled values are not one site.

---

## What we do not do

Stated so no consumer builds against it:

- **No indication decision.** We do not decide whether to pursue the indication.
- **No ranking.** We do not rank hypotheses or targets against each other. Two
  dossiers are not comparable on any field we emit; volume in particular is not
  calibrated for cross-target ordering.
- **No molecule design.** No structures, no chemotypes, no series suggestions.
- **No biologic assessment.** We record approved biologics as target validation
  and stop there.
- **No aggregate score.** See reading rule 1.
- **No clinical judgement.** Clinical failure is not evidence against
  tractability and never lowers a tractability number. RORgt has 152 holo
  structures, 12,900 compounds, 0.1 nM potency and zero approvals —
  VTP-43742 stopped on transaminase elevations, TAK-828F on preclinical
  teratogenicity. It is **small-molecule tractable and clinically failed**, and
  both belong in the dossier without either discounting the other. The
  terminations are in `target_precedent.terminated_programs` for you to weigh;
  weighing them is your station's job, not ours.

---

## Where this output is weak

Stated here rather than in a footer, because a consumer who knows where the
output is weak will use it better than one who trusts it uniformly.

**The uncalibrated numbers.** Three thresholds appear in the output and all three
are proposals, not calibrations, and the files say so at every appearance:

| threshold | status | basis |
| --- | --- | --- |
| volume ≥ 240 Å³ druggable / ≤ 210 Å³ hard | **RETRACTED, not merely uncalibrated.** Gates nothing and may not be revived | the calibration anchors did not measure the proteins they were attributed to — four of five hard anchors compromised, plus a RORgt entry containing no RORgt. See the volume section above |
| 4 Å off-site centroid distance | **proposal** | roughly half the one error ever measured (7.73 Å), resting on a single case |
| `buried_core_suspected` geometry gates | **proposal, not calibrated** | one observed case, no held-out set |

**The computed axis is often refused rather than reported.** One of the two
shipped reference dossiers has no pocket geometry at all. On homo-oligomers the
site-selection basis frequently degrades to something that cannot identify a
site, and the honest output is a null. Expect `tractability` to be sparse.

**The evaluation behind rule 4 is not followable from this checkout.**
`CLAUDE.md` cites `druggability_eval/RESULTS_TABLE.txt` and `all_rows.csv` as the
source of the 15-target / 67-structure / 134-measurement figures. Neither file
exists in this repository. Every downstream number that rests on that evaluation
— the 0.720 AUC, the 41% false-negative rate, and the now-retracted volume
separation — is currently a citation without a followable artifact. **The volume
separation is the case that shows what that costs:** it stood as the computed
axis's headline for part of a day, and what eventually falsified it was not the
statistics but a residue-level audit of the underlying structures, which nobody
reading this checkout could have run against the cited files.

**The precision of the reported figures is lower than it looks.** fpocket
estimates volume by Monte Carlo and the druggability score inherits that noise:
three identical reruns of one 5-structure ensemble gave CVs of 12.1 / 11.3 /
10.8%, so **about one percentage point of any reported CV is the method's own
noise**. The same structure read 0.673 on the deployed path and 0.708 locally.
Never read a difference in the third significant figure as a difference between
sites.

**Prose fields are unvalidated.** `axis_conflict`, `tractability.caveat`,
`clash_attribution` and the `not_found` reasons carry a large share of the
information and none of it is machine-checked. They are also where the most
decision-relevant caveats live. Surface them to a human; do not try to
regex-classify them.

**Extra keys appear.** Both shipped reference dossiers add keys the template does
not contain — `tractability.ensemble_consensus_fraction.source`,
`tractability.cryptic_potency_prior.source`,
`tractability.max_backbone_ca_displacement_note`,
`max_backbone_ca_displacement_protocol`,
`target_precedent.clinical_stage_small_molecules[].mechanism`. This is sanctioned
by `assemble-dossier`'s SKILL.md as an extension mechanism. **Do not use strict
schema validation** — parse permissively, read the 16 required top-level keys,
and ignore what you do not recognise.

**No worked example exists for two of the three verdicts.** Both shipped
dossiers are `small_molecule_tractable`. There is no reference dossier for
`not_tractable` or `insufficient_evidence`, and no reference dossier with an
`as_of_date` set — so the `as_of_leakage` machinery has never been exercised in a
shipped example. If your station handles those cases, test against a synthetic
one first.

**Figures that are withdrawn, and may still appear in older dossiers.** If you
are ingesting a backlog, treat every one of these as void rather than as data:

| figure | status |
| --- | --- |
| the **651-fold** (also written 650x) apo TNF-alpha druggability spread, and `±16%` beside it | withdrawn. A pocket-matching artifact — the matcher tracked a pocket 7.7 Å from the site it claimed, with 12.2 Å internal inconsistency |
| the volume range **206.7–309.2 Å³** for TNF-alpha | withdrawn, same cause |
| the CV improvement quoted as "**27.8% to 10.2%**" | withdrawn. The precision was never warranted; the measured figures are ~28% and ~10%, and about one percentage point of any CV is method noise |
| **8.83 Å** (KRAS) and **1.62 Å** (TNF-alpha) quoted as pipeline output | these are hand-calibration numbers from a protocol the deployment does not run. The deployed defaults are 8.65 Å and ~1.55 Å. The finding is the order-of-magnitude separation, not the decimals |
| a **1.97-log** systematic bias in the affinity predictor | overturned at n=17. Mean signed error is +0.23 log, CI (−0.28, +0.74) |
| "affinity prediction can be used to rank candidates within a target" | withdrawn. Within-target Spearman is +0.48, CI (−0.05, +0.77) |
| "PRANK rescoring has not yet helped, and once it hurt" | falsified at n=70 and void |
| the volume separation — **AUC 1.000, CI [1.000, 1.000]**, hard ≤ **207 Å³** / druggable ≥ **242 Å³**, stable under 15 leave-one-out refits | **retracted 2026-08-15.** The anchors did not measure the proteins they were attributed to; the CI was degenerate by construction; and a co-crystal flag alone separates the groups at AUC 0.900. See the volume section above |
| the per-target anchor volumes **MYC 188**, **IL-11 164**, **TNF 207**, **CD20 154**, **TL1A 137**, **KRAS 400 Å³** | void individually, not just as a set. MYC and IL-11 measured a different protein entirely; CD20's anchor is a detergent; TL1A had no anchor; TNF's had no overlap with its drug-anchored pocket; KRAS pooled two sites |
| **RORgt 386 Å³** and any figure sourced from **6C1P** | void. 6C1P contains no RORgt — its sole entity is A8EVM5, an ion transport protein — and its `1N7` anchor is CHAPSO, a detergent. RORgt is 428.5 Å³ with 6C1P removed |

**One internal contradiction to be aware of.** `assemble-dossier`'s worked pair
states that "TNF beats JAK1 on every pocket metric", and TNF-alpha's own
`axis_conflict` calls its pocket "the strongest pocket signal in the fixture
set". The later 15-target evaluation placed TNF-alpha at **207 Å³ in the hard
group** and JAK1 at **286 Å³ in the druggable group** — the opposite ordering on
the number that is primary. **The contradiction is now moot rather than
resolved:** those two figures are both void with the rest of the calibration set
(TNF's anchor had zero residue overlap with its drug-anchored pocket; its
defensible value is 129.6 Å³, against JAK1's 257.7 Å³), and there is no longer a
"hard group" or a "druggable group" for either to sit in. Neither the old
ordering nor the new one is a finding. **Do not build an ordering on any of
them.**

---

## Quick reference for a consumer

```
verdict            -> label, one of three. Never a number.
verdict_basis      -> which axis carried it. Read with verdict, always.
axis_conflict      -> non-null means the axes disagree. That IS the finding.
not_found[]        -> why every null is null. Parse it. Scan reasons for "ASK[".
tractability.caveat-> the run's own reservations. Prose. Show a human.
```

**Do:** key on `(accession, mechanism_hypothesis)`; read `verdict_basis` with
every verdict; check `not_found` before interpreting any null; read
`top_assay_share_pct` before believing `distinct_actives`; treat
`primary_d1_6_a3` as the computed number; run `validate_dossier.py` on what you
receive.

**Do not:** average the axes; merge the four precedent blocks; treat
`insufficient_evidence` as `not_tractable`; count a biologic as small-molecule
precedent; treat a null as a zero; threshold on `pocket_druggability`; threshold
on volume; sort neighbours on `evalue`; trust `has_druglike_holo` without reading
`ligand`; read a refusal as a low score; read an empty `terminated_programs` as a
clean safety record; block on a pending `ASK[...]`; validate strictly against the
template.
