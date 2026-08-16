# Fixtures

Targets arrive from the upstream pipeline, so this set is not a domain list —
it is a **ladder of increasing hardness**. Each rung isolates one way the agent
can be wrong, and the rungs are ordered so a failure tells you how far up the
system got before it broke.

Data for every target lives in `targets.json` with an `expected_output` block
that serves as the grading key.

**The classifier behind every `pdb_holo` / `pdb_apo` / `pdb_undetermined` figure
in this directory is checkable from this checkout**, as of 2026-08-15:

```bash
python3 ../.claude/skills/structure-select/tests/test_v2.py
```

Offline, stdlib-only, no Paperclip call. It prints the classifier's measured
accuracy — 259/262 on the original ground truth, 277/280 combined, 61/70 blind
with **0 false positives** — and `tests/version_diff.py` beside it will tell you
whether any two classifier versions actually disagree on a verdict. Until
2026-08-15 those harnesses existed only in session scratch, so the figures this
key rests on could not be checked by anyone reading the repo. Treat any fixture
number whose generating code is not in the checkout the same way.

## Files

| File | Status | Contents |
| --- | --- | --- |
| `targets.json` | retrieved and cited | Ten targets, `expected_output` grading keys |
| `pocket_calibration.json` | verified in-repo | KRAS holo vs apo — the backbone-collapse cryptic mechanism; the TNF-α mdpocket ensemble calibration; and the withdrawn 651-fold druggability spread, kept as the record of a retraction |
| `immunology_calibration.json` | found by execution | Four failure modes only real structures surface |
| `upstream_graph.json` | synthetic, marked `_fixture` | An upstream evidence graph, conformed to the real `SCHEMA.md` v1.1 |
| `upstream_graph_edgecases.json` | synthetic, hand-written | `kind: gene`, `basis: hedged_only`, a `no_effect` finding, `status: partial` |
| `upstream_graph_unknownverb.json` | synthetic, hand-written | Unknown `how` verbs — the `needs_adjudication` path, in three signal states |
| `upstream_graph_askback.json` | synthetic, two quotes verbatim-retrieved | The post-intake ask-back trigger — two links that should produce an ask, five that must not |
| `upstream_graph_expected.json` | derived, accessions retrieved | Grading key for `graph-intake` |

## The upstream graph

`upstream_graph.json` is where "targets arrive from the upstream pipeline"
stops being a sentence and becomes a file. It is a literature evidence graph in
its producer's own format — `things`, `links`, `findings`, `papers`, `gaps` —
and `graph-intake` reads it to fill this agent's input contract.

It is **synthetic** and carries `_fixture: true`, so its papers and quotes were
never retrieved from any corpus. That is acceptable here and nowhere else in this
directory: what it grades is the *extraction*, not the biology. `graph_read.py`
refuses it without `--allow-fixture` so the guard cannot be forgotten.

One nomination is correct: IRAK4 (Q9NWZ3), catalytic function, mechanism
`unknown`. The graded negative is IL-6 — `zimlovisertib reduces IL-6` has the
same shape as `zimlovisertib inhibits IRAK4` and only the verb separates a
readout from a target. That is the TNF-alpha assay-provenance failure of Rung 4,
moved one stage upstream where nothing else is looking for it.

### The fourth graph fixture — and why its negatives matter more

`upstream_graph_askback.json` grades a different thing from the other three:
not what the intake *extracts*, but what it decides to **send back**. Two links
should produce an ask; five must not, and the five are the point.

| link | outcome | why it is in the set |
| --- | --- | --- |
| `L3` | **ask** — `resolve_link` | A review asserts an oral small-molecule antagonist reached Phase 2. It would fill `clinical_stage_small_molecules`; the compound has no ChEMBL mechanism row, no registry record, and the only source is the review. All five gates pass. |
| `L1` | **ask** — `resolve_link`, post-resolution | The obefazimod/TL1A trap. We answer it ourselves from `chembl.drug_mechanism`, so it never blocks; the ask goes anyway so the wrong edge does not propagate. |
| `L2`, `L4` | no ask | `basis: primary`. A primary-supported claim is never an ask. |
| `L5` | no ask | An efficacy claim. Dossier rule 7 — it touches no tractability number, so it fails gate 1. |
| `L6` | no ask | A contested clinical status. `ctgov` settles it; fails gate 3. |
| `L7` | no ask | `rounds` already carries `resolve_link` at `L7`. Fails gate 4. |

A sixth negative is not expressible in the file, because `coverage.stop_reason`
is a single global value and this graph needs `max_papers` for `L3`. Set it to
`complete` and every ask must stop firing — verified: `--ask-context` then
reports no link clearing the gates at all.

Two of its quotes (`f1`, `f2`) are **verbatim from real papers**, PMC10762860
and PMC11642585, retrieved 2026-08-15. That is a deliberate departure from the
other three fixtures and it is recorded in the file's `_quote_provenance` block.
The rule was derived from those two exact sentences and a paraphrase would hide
why the trap works — both reviews call ABX464 "a prototype of TL1A", and they
share a senior author, so "two sources agree" is one source twice. The file is
still `_fixture: true` and nothing in it may be cited.

### Why three graph fixtures

The producer's schema gives `how` **no enum**, while every other categorical
field in it has one. So the verb that separates a target from a readout is open
vocabulary written by an upstream model, and our two verb lists can never be
complete. We do not get to ask for a closed vocabulary — this is ours to absorb.

That is why the set splits three ways. `upstream_graph.json` covers the happy
path with known verbs. `upstream_graph_edgecases.json` covers the four schema
branches the RA graph never reaches. `upstream_graph_unknownverb.json` covers
the case that has no clean answer: a verb we do not recognise, where the intake
must weigh the quote, the assay context and the graph shape — and is allowed to
refuse.

## The ladder

### Rung 1 — JAK1 (P23458). Can it do the easy thing?

Pre-formed ATP-site kinase. 14,342 compounds, best 0.010 nM, ruxolitinib 2011
plus two JAK1-selective approvals, **43 of 52 structures holo** (9 apo, 0
undetermined — re-derived from scratch 2026-08-15 under the chemistry
classifier; 42 was the superseded MW-window value and 40 reproduced under no
rule at all).

Worth knowing before you grade this rung: JAK1's *measured* fpocket
druggability is **0.009**, a median across nine approved drugs. If any criterion
is ever keyed to a druggability threshold, the easiest rung in the set fails
first.

The 14,342 is a *filtered* count — `n_target_components = 1`. Dropping that
predicate gives 14,472 and is wrong; see `_audit_2026_08_15` in `targets.json`.
Best potency was 0.032 nM here until the 2026-08-15 audit found that value came
from excluding every IC50 on a false "flagged" premise.

**Tests:** nothing subtle. If this is wrong, stop and fix the plumbing.
**Expect:** `small_molecule_tractable`, `cryptic_pocket_risk: low`.

### Rung 2 — RORγt (P51449). Does it confuse tractable with successful?

154 holo of **162** structures, 12,900 compounds, 0.017 nM potency — and **zero
approvals**. (152–154 all pass — two of the three entries the rule adds are a
DHEA sterol and an NDSB-256 crystallisation additive. 0.1 nM remains the best
cell-context IC50.) VTP-43742 stopped on transaminase elevations, TAK-828F on
preclinical teratogenicity, class-wide thymic lymphoma concern.

**Regenerated 2026-08-15: 162 total = 148 holo / 14 apo / 0 undetermined.** The
holo count **fell** from the superseded MW-window 154, only the second target
after BCL-2 to fall under the chemistry classifier — the apo bucket is almost
entirely the **hydroxycholesterol / sterol-agonist** structures, which the
classifier correctly declines to call small-molecule precedent. Note
`structures_by_accession` returns **165 rows for 162 entries** here, so
`COUNT(DISTINCT entry_id)` is load-bearing exactly as it is on KRAS.

> Two predictions this key made about the regeneration, one wrong. **4NIE's
> `DMX` (NDSB-256) was predicted to be caught by the Good's-buffer rule and is
> NOT** — it still classifies `druglike`, because that rule requires
> `aromatic == 0` and DMX's benzyl group gives it a ring. It costs nothing (4NIE
> is holo on `NBH` regardless), but the prediction reasoned from the chemical
> class without reading the rule's conjunction. **And the endogenous sterol in
> 29OB is comp_id `AND`, not `DHEA`** — a fixture keying on the literal string
> `DHEA` matches nothing. It classifies `druglike` as designed.

**Tests:** that clinical failure does not leak into the tractability number.

> **The 6C1P example is withdrawn, 2026-08-15.** This paragraph used to add "and
> now also that a near-zero pocket score does not, since RORγt's 6C1P measures
> 0.009 at rank 55 of 60." **6C1P is not a RORγt structure.** Verified against
> `pdb_v`: its sole polymer entity is `A8EVM5`, an ion transport protein — a
> NavAb sodium-channel mutant — and its only components are `PX4` (DMPC),
> `1N7` (**CHAPSO**) and phosphate. A near-zero score on a detergent site in a
> sodium channel tests nothing about RORγt.
>
> It never contaminated the *count*: `structures_by_accession` does not map it to
> P51449 (0 rows, against a confirmed 162-entry total). It reached the
> calibration set through **pocket selection** — chosen by `ligand_site_jaccard`,
> the route this repo calls trustworthy, anchored on a detergent. **Restricting
> scoring to the target's chains is necessary and not sufficient: a wrong PDB ID
> passes straight through a chain filter,** and nothing in the pocket pipeline
> announces it. Do not let 6C1P into any RORγt count or ensemble.
**Expect:** `small_molecule_tractable` **and** a populated `terminated_programs`.
Downgrading tractability because programs failed is the failure mode.

### Rung 3 — IL-17A (Q16552). Does it fall for modality?

Three approved antibodies (secukinumab 2015, ixekizumab 2016, bimekizumab 2021),
all reporting `action_type: INHIBITOR`, all with no chemical structure. **Zero
approved small molecules.** But 117 compounds exist with a best of 0.79 nM, and
LY3509754 reached Phase 1 before being halted for drug-induced liver injury.

Also the rung where structural plumbing breaks: 9SQX is CIF-only with a
five-character ligand code, and the site is a dimer-interface groove.

**Regenerated 2026-08-15: 44 total = 19 holo / 25 apo / 0 undetermined.** The
classifier lands back on the *original* 19, which the MW window had moved to 20
— but not by the same route, since 8DYG's `U5Q` (279 Da) is still called holo.
Two things a grader should know. **The five-character comp_ids were checked for
truncation and none occurred**: 11 five-character ids across this target and
RORγt, all resolved, all `druglike`, and the positive proof is that the CCD
lookup returned exactly 186 rows for a 186-component union — a truncated `A1J`
would have surfaced as verdict `unknown`, and the `unknown` count is 0. And
**8CDG sits in the apo 25 only because its sole component `UTF` is a 787 Da
peptidomimetic** typed `peptide-like`. That is the documented conservative
false-negative class, but on *this* target — the one whose whole point is that
macrocyclic-peptide chemistry is real precedent and not small-molecule
precedent — it should not be read as "nothing bound". It is also below the
1,000 Da hand-inspection threshold the BCL-2 degrader case established, so that
guard would not have caught it.

**Tests:** modality separation, and that neither "three approved drugs" nor "no
small molecules exist" is accepted as the answer.
**Expect:** zero approved small molecules stated explicitly, biologics in their
own block, real small-molecule chemistry acknowledged, the DILI termination
recorded.

### Rung 4 — TNF-α (P01375). Does it count assays or targets?

2,582 compounds — and **45% of all bioactivity comes from an IRAK4 assay
measuring a different protein**, labelled `assay_type = 'B'` so the obvious
filter does not catch it. Five approved biologics, zero small molecules. The
earliest holo ligand (`307`, 2AZ5) is a known promiscuous frequent hitter. The
site is a trimer-axis cavity that is **occluded, not cryptic** — steric
occlusion, not backbone collapse. It is pre-formed: delete the third chain and
all five apo structures recover the pocket, and the max backbone C-α
displacement at the site is ~1.6 Å. It therefore fails both community criteria
for cryptic (Vajda 2018; CryptoBench's apo-holo pocket-residue RMSD > 2 Å), and
must not be cited as a cryptic-pocket case. "Pre-formed" is a statement about the
subunit-removed state; in the intact trimer the third protomer is standing in
the site.

Holo count is **20 by the current rule, 19 defensibly** — regenerated
2026-08-15 under the chemistry classifier (52 total, 20 holo, 32 apo, 0
undetermined). The history is 15 (≥300 Da) → 17 (250–1200 Da) → 16 (that, minus
the spin label) → **20**. The jump is not drift: the classifier has **no lower
MW bound at all**, deliberately, and 16 inherited a 250 Da floor.

5UUI is a TNF-α carrying the MTSL nitroxide **spin label** on an engineered T77C
cysteine — a false holo, and the one measured **false positive** of the
classifier anywhere in this regeneration: MTN's chemistry really is drug-like,
and what disqualifies it lives in the entry title, which no chemistry test can
read. Subtract it → 19. The three sub-250 Da entries the new rule adds are not
noise: 6X81 (UTJ, 244 Da) and 6X83 (UTS, 208 Da) come from the *same
J. Med. Chem. paper and the same series* as 6X82/6X85/6X86, which every rule
counted, and UTS is the minimal benzimidazole core of the series 6OOY's A7M
belongs to. 4TWT (38A, 210 Da) is left unresolved and named.

**Accept 18, 19 or 20 with 5UUI named. Reject 15/16/17** — all artifacts of an
MW floor this rule does not have.

**Tests:** assay provenance, frequent-hitter detection, multi-chain handling, and
the occlusion mechanism — including that it is *not* reported as cryptic.
**Expect:** `axis_conflict` populated. Reporting 2,582 compounds as precedent is
the failure.

### Rung 5 — KRAS (P01116), `as_of_date = 2012-12-31`. Does it know what it cannot see?

Every pre-2013 structure is apo or GDP-bound. The switch-II pocket scores
**0.708 on holo and 0.000 on apo** — backbone collapsed ~8.8 Å at Glu63. That
figure is hand calibration: it comes from a protocol with auto-trim and
residue-name matching disabled and the mobile regions named by hand. The
deployed zero-knowledge default measures 8.65 Å max C-α displacement for KRAS
and ~1.55 Å for TNF-α, against hand figures of 8.83 Å and 1.62 Å. Mechanism and
`is_cryptic` are identical under both, so nothing downstream changes, but the
two sets are not interchangeable — quote what the run reported. The
order-of-magnitude separation (~8.8 vs ~1.6 Å) is the finding, not the decimals.

**Tests:** the cryptic blind spot, and the as-of cutoff.
**Expect:** `cryptic_pocket_risk: high`, and explicitly **not** `not_tractable`.
Run uncapped it must return sotorasib 2021 and adagrasib 2022. The agent must
not claim it would have found G12C early.

⚠️ **This rung no longer tests what the paragraph above used to claim, and the
change has not been applied.** It said to expect *low* computed tractability,
and that a 2012 run "says not tractable". Neither is derivable now: rule 4.2
makes druggability non-load-bearing and forbids it from carrying a verdict, and
the volume criterion that briefly replaced it has been **suspended** (its
anchors were measuring the wrong proteins). Taking the suspended band at face
value the answer was never "low" anyway — the apo switch-II pocket is 230 Å³ and
site-anchored KRAS is 226 Å³, both *inside* the unclassified band. A 2012 run
under current rules reports the computed axis as **unresolved**, not low. The
intent of the rung — "knows what it cannot see" — survives intact and is
arguably sharper. Proposed wording is in
`_expected_output_audit_2026_08_15` in `targets.json`; a human should rule.

### Rung 6 — MYC (P01106). Can it hold two contradictory facts?

1,079 compounds and **0 of 25 structures with any ligand above 122 Da** — holo 0
re-confirmed 2026-08-15, now under a fifth independent rule, with all 25 entries
strict apo. Intrinsically disordered. Best potency 0.2 nM from an assay described
only as "Inhibition of c-MYC (unknown origin)".

> **The 1,079 vs 1,249 disagreement is SETTLED, at 1,079.** `precedent-lookup`'s
> sweep table reported 1,249 compounds and flagged the pair as unreconciled. Both
> figures are correct measurements of *different quantities*: with the
> `n_target_components = 1` predicate MYC has **1,079 compounds / 1,249
> activities**; without it, **1,249 compounds / 1,675 activities**, because
> P01106 maps to seven ChEMBL targets including two PROTAC ternary complexes.
> The activity arithmetic closes exactly (1,249 + 426 = 1,675). This repo's
> authoritative definition mandates the predicate, so **1,079 is the compound
> count**.
>
> What made it look irreconcilable is a **numeric collision**: 1,249 is
> simultaneously the *unfiltered compound* count and the *filtered activity*
> count. Two different measurements landing on the same integer read exactly like
> one figure filed under two names. When two numbers disagree, check they are the
> same quantity before calling one wrong — and when they agree, check the same
> thing, because a coincidence looks just like a corroboration.

**Tests:** that reported actives with no holo structure read as conflict, not
precedent; that an uncharacterised assay is rejected however good the number.
**Expect:** `not_tractable`, `axis_conflict` populated, not rescued by family
precedent.

⚠️ **The computed axis now argues the opposite, on every sub-measure.** MYC's
D=2.4 druggability median is **0.75 — the highest in the set**, above KRAS 0.54,
BCL-2 0.52, JAK1 0.49, EGFR 0.44. The published consensus criterion ranks MYC
**top**. Persistence is chance-level. And volume, the last measure pointing the
right way, no longer does: the 188 Å³ that put MYC in the "hard" group was
measured on pockets lined entirely by MAX, by MAX plus DNA, and by apo OmoMYC,
and the corrected median is **325.7 Å³** — druggable.

There is a reason: **not one of MYC's 25 entries is wild-type MYC as an isolated
folded chain** (8 MYC:MAX dimers, 9 short MYC peptides on partner proteins, 6
fusion chimeras, 2 OmoMYC). Every pocket ever scored "on MYC" is on a partner, a
chimera or a designed miniprotein.

So `not_tractable` must now rest on **retrieved precedent** — 0 holo, 122 Da
ceiling, uncharacterised best assay — and `verdict_basis` must say so. The rung
now tests something harder and better than "hold two contradictory facts": it
tests whether the agent resists a computed axis actively shouting *druggable*.
An agent that reports the high druggability and top consensus rank **and still
returns `not_tractable`** has passed. Not applied to the ladder — see
`targets.json`.

### Rung 7 — IL-11 (P20809). Will it refuse?

15 compounds from two near-identical SPR assays on the same CAP chip
(CHEMBL6115567 ×11, CHEMBL6115571 ×4). **8 structures, 0 holo, 8 apo, 0
undetermined** — measured 2026-08-15, where before it was asserted without a
number behind it. No drugs at any phase. Just enough data to tempt a confident
score.

Sharper still: **only 2 of the 15 rows carry a number at all** — 140 nM and
2,600 nM. The other 13 have `standard_value` NULL. Reporting "15 activities" as
15 measurements already overstates the evidence 7.5×.

And `best_potency_nm` is now **null**, not 140. Both assay descriptions contain
"(unknown origin)" — the exact string this fixture's own definition says to
reject, and which it *did* apply to KRAS. The file was applying its own rule to
one target and not another; that is provable from the file against itself, with
no query. The numbers stay visible because an agent will find them, and finding
them is the temptation; asserting them is the failure.

**Tests:** the hardest thing to make a system do — decline.
**Expect:** `insufficient_evidence`. Any number here is a failure.

## Retained for method validation

EGFR (P00533), BCL-2 (P10415) and TYK2 (P29597) stay in `targets.json`. EGFR and
BCL-2 are unambiguous positives useful for regression. TYK2 earns its place on
dating: JH1 structures from 2010-06-02, JH2 pseudokinase from 2013-04-10,
deucravacitinib approved 2022. A 2012 cutoff must show no allosteric precedent
and a 2015 cutoff must show it — the cleanest as-of test in the set, because it
tests a *pocket* appearing rather than a drug.

TYK2's domain split was re-derived independently on 2026-08-15 and reproduces
exactly: **JH1 28, JH2 20, both 1 (4OLI), other 3 = 52**. These are *exclusive* —
4OLI is a JH1+JH2 tandem and is reported in its own bucket, not added to either
domain. The older inclusive 29/21 double-counted it and summed to 54 against a
real total of 52; it stays rejected. TYK2 also now carries a holo count for the
first time: **46 holo / 5 apo / 1 undetermined**.

BCL-2 was the one place where this key was worse than the better answer.
**Resolved 2026-08-15, and the numbers this paragraph used to carry — 40 graded,
42 better — are both void.** They were the MW-window generation. BCL-2 has since
been regenerated: **39 holo / 23 apo / 1 undetermined**, history `[41, 40, 39]`,
and 40 is now the *superseded middle value*, unreachable under the chemistry
classifier. The defensible count is **41** — 39 plus the two bivalent degrader
co-crystals 8FY1 (YF8, 1640 Da) and 8FY2 (YFH, 1682 Da), which are genuine
beyond-rule-of-5 precedent a scientist counting BCL-2 should count.

**The resolution is that the graded criterion on this target is not a number at
all.** Grade the *disposition of 8FY1 and 8FY2*. Both 39 and 41 pass when the
agent says what it did with them; a bare integer fails at every value. Grading
39 and merely tolerating 41 relegates the better science to a tolerance band;
grading 41 instead would penalise an agent for running `classify_record` exactly
as this repo specifies, which teaches hand-adjusting tool output toward an
expected number — a worse habit than an off-by-two count. Neither scalar is what
separates a good dossier from a bad one. Whether the agent *saw the two entries*
is. **Reject 40 outright**, disposition or not: reaching it means using the
retracted MW window.

And the earlier prediction that both degraders would land in `undetermined` was
**half wrong**, which is the more interesting half. Only 8FY1 does. 8FY2's YFH
never reaches the 1200 Da ceiling — rule R10 fires first, reading three
α-amino-acid backbone units in the linker and amide-rich warhead, so it is filed
`peptide_or_polymer` and counted **apo, silently**: no flag, no `unknown`,
nothing for a reader to notice. Degraders are amide-rich by construction, so
this misfires across the whole modality and it misfires quietly. That is the
generalisable reason to grade the disposition rather than the integer — a key
that grades the integer cannot tell an agent that inspected the
`peptide_or_polymer` bucket from one that did not.

## How the structure counts are defined (regenerated 2026-08-15)

Every `pdb_holo` / `pdb_apo` / `pdb_undetermined` figure comes from the
**chemistry classifier**, `.claude/skills/structure-select/ligand_filter.py`
(sha256 `ae9ef140f3f8…5216c1`, 2,430 lines), not from a molecular-weight window.
**The hash previously recorded here — `526610951ee1…89f1f3`, 1,490 lines — named
no file in this checkout**; it identified a pre-polymer-conjugate copy that
existed only in session scratch. Corrected 2026-08-15; the full finding, and the
check that the difference moves nothing in this key, is in `targets.json` under
`_classifier_pin_corrected_2026_08_15`. An entry is
holo iff ≥1 of its `pdb_v.entry_ligands` components classifies `druglike`.

Three things follow, and all three matter for grading:

- **There is no lower MW bound any more.** That is deliberate — size was never
  the discriminator, and the old floors were splitting congeneric series across
  the holo/apo line. This is why TNF-α moves 16 → 20 and JAK1 42 → 43.
- **`unknown` is not `apo`.** An entry with no drug-like component but an
  unclassifiable one is `undetermined`, and `total = holo + apo + undetermined`
  always. TYK2 5C01 is the worked case: its only component is `UNL`,
  "UNKNOWN LIGAND". An agent reporting TYK2 apo as 6 has made the error the
  classifier exists to prevent.
- **A lookup failure is not a chemistry miss.** Failures carry `lookup_failed`
  into `holo_call(...)["undetermined"]` and must never render as apo.

**All ten targets are regenerated** (plus CD20 as a control), and every one now
satisfies `pdb_total = holo + apo + undetermined`. Nothing carries a
`_NOT_REGENERATED_2026_08_15` flag any more. The record of which half was which
is kept in
`_structure_regeneration_2026_08_15`, because a half-regenerated key is only
dangerous when nobody can tell.

> **Corrected 2026-08-15.** This paragraph previously said six regenerated and
> listed BCL-2 among the five stale ones, as did the `targets.json` index block.
> Both were wrong and they contradicted the BCL-2 record itself, which carries
> `pdb_holo_measured`, a three-entry history and a *measured* degrader block.
> Only four targets ever carried the `_NOT_REGENERATED_2026_08_15` flag. The
> index disagreed with the territory, which is precisely the failure the index
> exists to prevent — so when you regenerate a target, fix the index in the same
> edit.

| target | total | holo | apo | undet | was |
| --- | ---: | ---: | ---: | ---: | ---: |
| KRAS | 522 | 300 | 221 | 1 | 301 |
| EGFR | 388 | 283 | 105 | 0 | 283 — same total, different entries |
| BCL-2 | 63 | 39 | 23 | 1 | 40 |
| MYC | 25 | 0 | 25 | 0 | 0 |
| TYK2 | 52 | 46 | 5 | 1 | — |
| JAK1 | 52 | 43 | 9 | 0 | 42 |
| IL-17A | 44 | 19 | 25 | 0 | 20 |
| RORγt | 162 | 148 | 14 | 0 | 154 |
| TNF-α | 52 | 20 | 32 | 0 | 16 |
| IL-11 | 8 | 0 | 8 | 0 | — |

**The direction prediction in the old notes was wrong.** They said to expect holo
counts to *rise*, since abolishing the 250 Da floor added entries "on every
target where it was measured". Of the four regenerated last, **three fell**
(BCL-2 −1, RORγt −6, KRAS −1) and one was flat. The floor does admit entries;
the sterol, lipid and additive rules removed more. Generalising a direction from
whichever targets were measured first is extrapolation.

**EGFR is the sharpest lesson here: it did not move, and that is not a no-op.**
283 under the MW window and 283 under the classifier — but the two rules disagree
about *which* entries, and the older ≥300 Da rule gave 282. Two different
definitions landing on one integer. It is the cleanest argument in this file for
recording a count's definition beside it.

**And "Paperclip's SQL backend failed partway through" was the wrong diagnosis.**
The backend was healthy. `cli_cwd` in the Paperclip *client* config was set to
`/papers/`, and the CLI `cd`s into it before every command — making each `sql`
call a ~15 ms no-op that prints `vsh: cd: /papers/: Permission denied` **on
stdout with exit code 0**. Nothing about it looks like a failure to a caller
checking the exit status. Isolate with `PAPERCLIP_CONFIG_DIR` before
re-measuring rather than resetting the shared config; the shared one has been
observed being re-poisoned minutes after a reset.

If you re-measure: the backend silently served a **moving row cap** — the same
query returned 200 rows at one moment and exactly 10 at another, well-formed and
without error, which recorded KRAS as 10 entries against a true 522 on the first
attempt. Reconcile every paged read against a separately issued `COUNT`, under a
canary on the `pdb_v` tables you are actually reading. `SELECT 1` proves nothing.

## Rules for the grading keys

- Every value carries a source: ChEMBL target ID, PDB ID, DOI, or line-pinned URL.
- Every count carries **the definition it was measured under and the date**.
  A grading key with no stated definition drifts; that is exactly how this
  happened.
- Unretrievable is `NOT_FOUND`, never an estimate. A query that timed out is
  unretrievable — it is never a zero.
- Approved drugs split by modality; the two columns are never summed.
- Dates matter — first-approval years and PDB release dates, not just counts.
- **No criterion may be keyed to an fpocket druggability score.** It has a
  measured 41% false-negative rate on structures with a drug physically bound,
  and it is min-max normalised across the other pockets of the same structure,
  so it is not even a property of the pocket.
