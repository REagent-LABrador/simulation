---
name: assemble-dossier
description: >
  Assembles the finished druggability dossier from a UniProt accession — walking
  identity, retrieved precedent, structure selection, pocket geometry, cryptic
  mechanism, interface classification and a falsification sweep into the JSON
  template, then gating it on validate_dossier.py before returning. It does NOT
  decide whether to pursue an indication, does NOT rank targets against each
  other, does NOT average the two axes into a score, and does NOT invent a value
  to fill a field.
---

# assemble-dossier

This is the station that produces the output. `precedent-lookup`,
`structure-select`, `pocket-scan` and `falsification-sweep` all feed it, and
none of them returns a dossier — they return blocks. This skill is where the
blocks become the object the caller receives, and it is the only place where the
verdict is written.

Two things make it different from the four skills upstream:

- **It is the only step allowed to write `verdict`.** Every upstream skill's
  frontmatter says it does not decide. This one decides *what to report*, which
  is not the same as deciding whether to pursue the target. See "Choosing the
  verdict".
- **It ends with a machine gate.** `validate_dossier.py` sits next to this file
  and must return zero violations before the JSON is returned. It does not check
  whether the dossier is right; it checks whether it is honest about what it
  knows. Failing it is not a warning.

## The contract in one paragraph

Two axes, reported separately, allowed to disagree, never averaged. Every number
carries a source. Anything not retrieved is `null` with a line in `not_found`,
never an estimate and never a zero. `insufficient_evidence` is a correct answer
and must stay reachable. The verdict is a label, not a number.

## Procedure

Run the steps in order. Later steps depend on earlier ones — chain selection in
step 4 depends on the mechanism question raised in step 3, and step 7 attacks
claims made in steps 2 through 6.

### 0. Input contract, when the input is a graph rather than an accession

**Skill:** `graph-intake` — only when the caller hands you an upstream literature
evidence graph instead of a filled input row. It extracts which nodes are protein
targets, what the molecule is meant to stop, and the accession, and emits
follow-up asks when the answer is missing. It does not decide a mechanism the
evidence has not stated, so a missing `mechanism_hypothesis` stays missing and
step 4 reports the biological assembly without asserting which pocket matters.

Skip this step when `uniprot_accession` was supplied directly.

### 1. Identity

**Skill:** `precedent-lookup`, step 1 (`uniprot_v.proteins`).

Resolve the input to an accession first. If given a gene symbol, record both —
gene symbols are ambiguous across organisms and the accession is the join key
for every other query in this procedure.

**Writes:** `target.*` and a `target.sources` entry naming UniProt. Also
`as_of_date` verbatim from the input, and `verdict_basis` is left unset until
step 8.

`sequence_length` is a number, so it needs a source like every other number.
`target.sources` is the extension key that provides it.

### 2. Precedent — modality first, then potency, then provenance

**Skill:** `precedent-lookup`, steps 2 through 6. **Then:**
`falsification-sweep` check 2b for the cross-accession test.

Four sub-steps, in this order, because each one can invalidate the next:

**2a. Modality.** Read `chembl.molecule_dictionary.molecule_type` per drug.
`Small molecule` goes to `target_precedent.approved_small_molecules`; `Antibody`
and `Protein` go to `biologic_precedent.approved_biologics`; `Unknown` and NULL
go to `not_found` and count toward neither. Collapse salt and parent
`molregno`s before counting, or state that the figure is a row count.

**2b. Assay provenance, before any actives count is written.** Group
bioactivities by assay description, take the top contributor and its share. If
one assay exceeds ~30%, the count is about that assay. Then read the description
and ask whether it measures this protein at all. `assay_type = 'B'` is not a
defence and must not be used as a filter.

**2c. Potency, with characterisation.** An uncharacterised assay description
makes a number unusable however good it is. Prefer a weaker number you can
characterise. Record `best_potency_characterised` explicitly — the validator
requires it whenever a potency is reported.

**2d. Cross-accession check.** Does the top assay appear under other
accessions? That is how a cellular readout for a different protein gets caught.

**Writes:** `target_precedent.*`, `biologic_precedent.*`, `family_precedent.*`.

### 3. Structures

**Skill:** `structure-select`.

Classify holo and apo by actual ligand chemistry using the exclusion list, not
by label. Apply the `as_of_date` at the source (`release_date <= cutoff`) and
record how many entries it removed. Pick the tier by the strict order: holo
experimental, apo experimental, predicted. Assemble an ensemble — one structure
is not a measurement.

Ask the biological-assembly question here, not in step 4: is the site on one
chain, at an interface, or on an oligomer axis? Prepare TNF-alpha as one chain
and its site does not exist.

**Writes:** `structure.*`, `structural_neighbour_precedent.*`, and the ensemble
that step 4 consumes.

### 4. Pockets

**Tool:** the Modal function `pocket_scan` in
`.claude/skills/pocket-scan/modal_app.py`, deployed as
`druggability-pocket-scan`. One invocation for the whole ensemble:

```python
import modal
fn = modal.Function.lookup("druggability-pocket-scan", "pocket_scan")
result = fn.remote(pdb_ids=[...], chains={...}, ligand_codes=[...])
```

Four calls pay four cold starts; one call pays one. Pass `ligand_codes` only as
an override. Pass `site_residues` when you have a site definition, because
without one an apo structure falls back to "the most druggable pocket anywhere
in the chain". **`chains` and `site_residues` both reach the tool now** — rule
2b's chain selection is directly expressible (`chains={"1TNF": ["A","B"]}`) and
the subunit-removed control is reachable, which is what separates "the cavity is
too small" from "a protomer is standing in it". A chain flag is not always
enough: a fusion chaperone can sit *inside* a chain (3V2Y's T4 lysozyme at
1002–1161 beside the receptor at 16–330), which needs a residue range.

**Read `site_pocket_selected_by` on every value before using it.** It is
returned per structure per clustering value, and **four of its six possible
values** mean the number does not describe a known site:
`max_druggability_no_ligand_site` ("the most druggable pocket anywhere in the
chain"), `site_signature_unreliable_homooligomer` (a residue-number match a
homo-oligomer makes ambiguous in principle), `no_pocket_matched_site_signature`
(matched nothing) and `no_pocket_overlapped_ligand_site` (a ligand site was
known and no pocket reached it — the strongest case of all, because the
measurement is anchored to nothing). Values carrying any of them must not be
pooled into one spread.

**This sentence read "two of its five" until 2026-08-15, and both numbers were
wrong.** `no_pocket_overlapped_ligand_site` was never transcribed into
`SELECTION_BASES` at all, so the validator rejected a basis the tool emits — and
it is the basis for the exact false negative rule 4 exists around (TNF-alpha
0.002 at D=1.6 on a co-crystallised 570 Da ligand, the cluster discarded below
fpocket's `-i 15` floor). `test_four_of_the_six_bases_do_not_identify_a_site`
now pins both counts.

**Writes:** `tractability.pocket_volume_a3` — including
`primary_d1_6_a3`, **the computed axis's primary number**, and `clustering_d`
saying which D values the spread pools — `tractability.pocket_druggability` (a
range across D and across structures, never a point, with `load_bearing: false`
and `_false_negative_rate`), `tractability.site_pocket_rank`,
`tractability.ensemble_consensus_fraction`, `tractability.method`, and
`site_pocket_selected_by` on each spread.

**Volume at D=1.6 first — but not because it separates anything.** It is the
computed axis's primary number because D=1.6 and D=2.4 measure different
cavities and a pooled spread is a spread over two different things, so one of
them has to be carried explicitly. **A volume is a measurement of a cavity in a
structure we scored. It carries no verdict, and there is no volume at which this
output calls a target druggable or hard.**

> **RETRACTED 2026-08-15 (`CLAUDE.md` rule 4a, `rubric.md`, `OUTPUT_NOTES.md`,
> `pipeline.html` N2).** This passage previously reported, as a current result,
> that volume at D=1.6 separated all 15 evaluation targets at **AUC 1.000** with
> a CI of [1.000, 1.000], stable under all 15 leave-one-target-out refits.
> **That claim is withdrawn in full.** It is recorded rather than deleted so a
> reader meeting the number in an older dossier can find out what happened to
> it. An audit resolved every lining residue of all 67 calibration structures to
> an entity by aligning SEQRES to the target's UniProt sequence, and found that
> the anchors did not measure the proteins they were attributed to: **MYC's
> pocket is 100% MAX (P61244)** with zero MYC lining residues; **IL-11's is 100%
> IL-11 receptor alpha (Q14626)**; **TNF's** sits on TNF but has **zero residue
> overlap** with its only drug-anchored pocket; **CD20's anchor ligand is
> cholesterol hemisuccinate, a detergent**; **TL1A had no anchor at all**; and
> **RORgt's 6C1P contains no RORgt** — its sole entity is an ion transport
> protein, and it arrived through `ligand_site_jaccard` anchored on CHAPSO,
> another detergent. `chain_accessions` was `{}` on **every** entry while the
> adjacent field asserted chains were resolved from `_struct_ref_seq`.
>
> Two points generalise past the anchors and should travel with any revival
> attempt:
>
> - **A bootstrap CI on a perfectly separated set is degenerate by
>   construction.** Resampling cannot produce an inversion that is not in the
>   data, so `[1.000, 1.000]` was arithmetic, not evidence.
> - **The label and the measurement share a cause.** The binary flag "a
>   drug-like ligand was co-crystallised" separates the same groups at **AUC
>   0.900** with no structural measurement at all.

Druggability's own numbers are reported below for the separate purpose of
demoting it, and that demotion is unaffected by the above — it does not turn on
any retracted anchor. Druggability manages **0.720 with a 95% CI of 0.44–0.94**
at D=1.6 and **0.520** at D=2.4. On holo structures with a drug-like ligand
physically bound and the pocket anchored to it, a large fraction of druggability
scores fall below 0.1 — EGFR 6LUD with osimertinib bound is 0.013 (see the
denominator note below). Druggability is a
three-descriptor logistic regression fitted on 21 positives, and
`mean_loc_hyd_dens_norm` is min-max normalised across the *other* pockets of the
same structure, so it is a property of a pocket relative to whatever else was
detected beside it. Quote it as a weak prior with its provenance attached and
**never let it carry a verdict**.

Three traps specific to the primary number:

- `primary_d1_6_a3` is **not** the min of the pooled spread. D=1.6 and D=2.4 do
  not measure the same cavity. Record the D=1.6 value or decline it.
- A site volume above ~1000 Å³ means sites have merged with neighbours. That is
  the normal D=2.4 failure and it is not a big pocket.
- **The volume guide is RETRACTED, not merely uncalibrated.** It was previously
  described here as an uncalibrated proposal resting on a real separation: ≥240 Å³
  falling entirely in the druggable group and ≤210 Å³ entirely in the hard group.
  **That separation is withdrawn** (rule 4a, 2026-08-15) — it was fitted post hoc
  on n=15 with only 5 hard targets, and four of those five anchors did not measure
  the target protein at all. **Do not quote it, do not mark it "uncalibrated" and
  quote it anyway, and do not revive it from the constants in
  `validate_dossier.py`.** `VOLUME_GUIDE_DRUGGABLE_A3` and `VOLUME_GUIDE_HARD_A3`
  still exist there deliberately, but they are read by exactly one rule as a
  *disclosure trigger* — deciding when a low druggability sitting beside a large
  volume disagree loudly enough to force a `tractability.caveat`. That survives
  the retraction precisely because it asserts nothing about which side of the
  number a target falls on. Nothing gates on them; do not add a rule that does.

**`prank_rank` goes in `site_pocket_rank.prank`, beside `site_pocket_rank.fpocket`,
never instead of it.** PRANK promotes the true site in 79% of 70 ligand-anchored
measurements and demotes it in 1% (the 6OIM KRAS case, kept visible) — median
rank 5→1, top-3 recall 37%→91%. It is a **site finder**. As a druggability
classifier its rank is inverted, **AUC 0.25**, because on an unanchored target
the top pocket is top by construction.

**Do not substitute persistence for the demoted score.** It is the obvious wrong
fix: the site pocket was detected in 100% of structures for all 15 targets, so
persistence is **AUC 0.500** — exactly chance — and the published consensus
criterion built on it gives AUC 0.560 and **ranks MYC first at 0.80**, above 8 of
the 10 druggable targets.

### 5. Cryptic mechanism

**Skill:** `pocket-scan`, step 6, plus `cryptic_analysis.py`.

Measure it; do not flag it from tier. Where both apo and holo exist, superpose
and compute max backbone C-alpha displacement at the site and clash attribution.
Classify on **displacement**, never on which atoms clash — keying on clash
composition classifies KRAS as side-chain occlusion and hands the canonical
nanomolar target a micromolar prognosis.

Then apply the field's definition before calling anything cryptic: cryptic means
the pocket is absent in **all, or nearly all**, unbound structures. A site
present in the apo ensemble is **occluded**, not cryptic.

**Writes:** `tractability.cryptic_pocket_risk`,
`tractability.cryptic_mechanism`, `tractability.cryptic_potency_prior`,
`tractability.max_backbone_ca_displacement_a`,
`tractability.clash_attribution`, and the `tractability.cryptic_evidence`
extension block that carries the apo census the claim rests on.

Mechanism is a prior on achievable potency, so the mechanism and the ceiling
must agree: loop or backbone motion means nanomolar is reachable, side-chain or
subunit occlusion means micromolar at best.

### 6. Interface classification

**Module:** `interface_analysis.py` (`interface_residues`, `classify_pocket`).

Only when a complex structure containing the partner exists. Compute interface
residues and classify each pocket: overlapping the interface is an
`orthosteric_candidate`, distal is `allosteric_candidate`, buried within the
oligomer is `destabiliser_candidate`. With no partner structure, the value is
`no_partner_structure` — the literature knowing the mechanism is not a
measurement, and rule 2b requires this one be measured.

**Writes:** `tractability.pocket_vs_interface.*`.

### 7. Falsification sweep

**Skills:** `falsification-sweep`, all nine checks. For check 7 — what was tried
in the clinic and why it stopped — use `terminated-programs`, which joins
ClinicalTrials.gov `why_stopped` to the literature account of the same event and
reports contradictory sources side by side with their dates rather than picking
one. A contested termination stays contested in `terminated_programs`.

Record checks that found nothing alongside checks that found something. "We
looked and found nothing" is information; silence is not. Absence of *retrieved*
terminations is not evidence that none exist — say which registry was queried,
as both worked examples do.

**Writes:** `falsification.checks_run`, `falsification.findings` (each with a
source), `falsification.survived` — a boolean, never null after a run.

Nothing here changes a number. A finding attaches evidence; it does not apply a
discount.

### 8. Verdict, then validate

Write `verdict`, `verdict_basis`, `axis_conflict` and `next_experiment`. Write
the dossier to its pinned path, then run the gate **on the file you wrote**:

```bash
python3 .claude/skills/assemble-dossier/validate_dossier.py \
        /mnt/session/outputs/druggability-dossier.json
```

**This command used to read `python3 validate_dossier.py dossier.json`**, which
named a file `CLAUDE.md` forbids writing (the output path is pinned and no
other) and a script path that only resolves if the working directory happens to
be this skill's own. Validate the artifact the grader reads, not a copy of it —
a copy is where the two can differ and nobody finds out.

Zero violations, or the dossier does not go out. Fix the dossier, not the
validator — every rule in it exists because the corresponding mistake is easy,
plausible and invisible in the output.

## Choosing the verdict

The verdict is one of three labels. It is about **small-molecule tractability
only** — not about whether the indication is worth pursuing, which other
stations answer.

| condition | verdict |
| --- | --- |
| approved or clinical small molecules exist, or characterised potent chemistry against a real site | `small_molecule_tractable` |
| structure or chemistry positively argues no small molecule can bind — an epitope with no cavity, an IDR with no holo structure anywhere | `not_tractable` |
| thin chemistry, no holo structures, no approvals — not enough to say either way | `insufficient_evidence` |

**A low druggability score is not a route to either negative verdict, and the
validator now blocks it.** A large fraction of pockets with a drug-like ligand
physically bound score below 0.1 — EGFR 6LUD with osimertinib in it is 0.013,
JAK1's median across nine approved drugs is 0.009, TYK2 6NZP with
deucravacitinib is 0.169, BCL-2 6QGK is 0.025. (**That fraction was reported as
41%, 15 of 37, and its denominator is now under audit.** RORgt 6C1P — cited here
until 2026-08-15 as `0.009 at rank 55 of 60` — is **not a certain positive** and
has been struck: the entry contains no RORgt, its sole entity is an ion
transport protein, and its anchor ligand is CHAPSO, a detergent. So the
denominator is **36, not 37**, and **the other 36 have never been audited at
residue level.** Quote the direction and the named cases, not the rate.) A
`not_tractable` reached on that number is a false negative in the most expensive
direction available. `not_tractable` on computed grounds needs the D=1.6 volume
behind it; if the volume was not measured, the answer is `insufficient_evidence`
**with the unmeasured volume named as the reason**, not a poor pocket.

**And every historic verdict that leaned on a low druggability score is flagged
for re-examination.** The failure is systematic across all 10 known-druggable
targets, not a handful of outliers, so this is not a matter of spotting the odd
bad call — assume any such verdict is unsupported until the volume is measured.

`verdict_basis` names which axis carries it: `retrieved_precedent`,
`computed_tractability`, `both`, or `none`. Without it, "tractability claimed on
precedent grounds" is not a checkable statement and the modality rule cannot be
enforced.

**Clinical failure never lowers a verdict.** RORgt has 152 holo structures,
12,900 compounds, 0.1 nM potency and zero approvals because of transaminase
elevations and teratogenicity. It is `small_molecule_tractable` **and**
clinically failed. Both go in the dossier; neither discounts the other.

**`axis_conflict` is a field, not a verdict value.** The fixture set writes TNF's
expected outcome as "axis_conflict" — that means the verdict is accompanied by a
populated `axis_conflict`, not that the verdict string is `axis_conflict`. The
validator rejects any other value.

Populate `axis_conflict` when the axes genuinely disagree — approved biologics
with no small molecules and a tractable finding; hundreds of actives against
zero holo structures; a dominant assay measuring a different protein; a headline
potency from an uncharacterised assay. Do not populate it out of caution. On
CD20 both axes agree there is nothing to bind, and a manufactured conflict there
is filler.

## Extensions to the template

The template in `CLAUDE.md` is the contract and no key in it may be omitted or
renamed.

**This section used to list nine extensions. All nine have since landed in the
template itself** — `verdict_basis`, `target.sources`, `tractability.sources`,
`approved_small_molecules_count`, `as_of_leakage[]`, `structure.apo_count`,
`pocket_*.site_pocket_selected_by`, `cryptic_evidence` and
`ensemble_consensus_fraction.n_measurements` are now template keys, not
additions, and they must be filled or nulled like any other. Do not treat them
as optional because this table once did.

**Six more landed on 2026-08-15 with the rule 4 re-prioritisation, and they are
template keys too:** `pocket_volume_a3.clustering_d`,
`pocket_volume_a3.primary_d1_6_a3`, `pocket_druggability.load_bearing`,
`pocket_druggability._false_negative_rate`, `pocket_volume_a3._primary_note` and
the `tractability.site_pocket_rank` block. Two of them are `_`-prefixed and so
are invisible to the provenance walker by design — they carry the reason a
number is discounted, not a claim. `load_bearing` has exactly one legal value.

**The first four of those six reached the validator before they reached the
template, and an agent filling the template exactly as written produced a
dossier the gate rejected.** That is the regression `test_template_drift.py`
exists for; see "Nobody was reading the template" below. Do not fix a
recurrence by editing the validator — the template is where the key is missing.

**And a seventeenth top-level key landed the same day: `input`**, echoing the
five contract fields verbatim. It exists because a downstream team was told to
key a cache on (accession, `mechanism_hypothesis`, `as_of_date`) and none of the
three survived into the output. It is an **echo**, never inferred and never
back-filled: `input.uniprot_accession` is what the caller said,
`target.uniprot_accession` is the resolved accession, and the top-level
`as_of_date` remains authoritative. Both worked examples carry it.

What remains genuinely added, and why each earns its place:

| key | why |
| --- | --- |
| `tractability.ensemble_consensus_fraction.source` | names the specific run and the raw `2 of 4` / `4 of 4` form. `tractability.sources` names the run but not the count, and the numerator survives nowhere else |
| `tractability.cryptic_potency_prior.source` | the CryptoSite citation (Lazou et al. 2024) is a *different* source from anything in `tractability.sources`, which names in-repo calibration and fpocket. Dropping it would delete a citation, not a duplicate |
| `target_precedent.clinical_stage_small_molecules[].mechanism` | on TNF-alpha, "stabilises a receptor-incompetent trimer rather than blocking the TNF/TNFR interface" is the fact that makes `axis_conflict` resolvable at all, and it is the rule 2b oligomer-destabilisation case in one line |
| `falsification.findings[]`, `not_found[]` element shape | the template leaves both lists as `[]`, so the element shape has to come from somewhere: `{check, result, source}` and `{field, reason}` |

**A derived number is not an extension, it is a liability.** The examples used to
carry `ensemble_consensus_fraction.n_measurements_meeting_criterion`. It is
exactly `fraction_with_strong_pocket x n_measurements`, so it adds nothing and
can silently drift out of agreement with the fraction if one is edited and the
other is not — a self-contradicting dossier that no rule catches. It has been
removed. Provenance strings are not derivable and stay; arithmetic is derivable
and goes.

The template shows list shape with one all-empty object. That is illustration,
not data: use `[]` for "none", never an entry with an empty name.

## The gate

`validate_dossier.py` — pure stdlib, importable as
`validate_dossier(dossier: dict) -> list[Violation]`, or runnable as a CLI that
exits 1 on any violation. **Seventeen rule functions, eighteen violation types**
— and both numbers are asserted by `test_the_rule_inventory_is_pinned`, so this
sentence cannot go stale the way the rubric's test count did:

| rule | catches |
| --- | --- |
| `WELL_FORMED` | missing keys, unfilled enum placeholders, nameless template stubs, `survived: null`, NaN |
| `NUMBER_WITHOUT_PROVENANCE` | any numeric leaf with no source in scope |
| `AXES_AVERAGED` | any overall/composite score, and any number combining a precedent term with a tractability term |
| `MODALITY_LEAK` | a drug in both blocks, a `-mab`/`-cept` stem among small molecules, tractability claimed on precedent that is entirely biologic |
| `INSUFFICIENT_EVIDENCE_AVOIDED` | thin chemistry plus zero holo plus zero approvals answered with a confident verdict |
| `DRUGGABILITY_POINT_ESTIMATE` | a scalar druggability, a one-sided range, a single clustering value, an unnamed ensemble, druggability without volume |
| `DRUGGABILITY_LOAD_BEARING` | **new, 2026-08-15.** `load_bearing` anything but `false`, a range with no `_false_negative_rate` beside it, a `not_tractable`/`insufficient_evidence` verdict on computed grounds with a low druggability and no volume, and a low druggability sitting beside a large volume with an empty `tractability.caveat` |
| `VOLUME_NOT_PRIMARY` | **new, 2026-08-15.** computed-axis geometry reported without `pocket_volume_a3.primary_d1_6_a3` and without a `not_found` line, a primary volume above 1000 Å³ (sites merged), or a spread with no `clustering_d` record |
| `FRACTION_WITHOUT_N` | a consensus fraction with no N, or an N disagreeing with the named ensemble |
| `SAME_SITE_BASIS_MISSING` / `_INVALID` | a pooled spread with no recorded basis, or pooled on a basis that does not identify a site |
| `SITE_INCONSISTENT` | geometry quoted off `site_from_density` when its centroid is past the 4 A proposed threshold, or an interface class asserted with no partner structure and no measured overlap |
| `INTERFACE_MIXED_UNRESOLVED` | a `mixed` interface classification that does not say what it is mixed between, carries one overlap instead of the disagreeing ones, names no partner, or claims to confirm a mechanism hypothesis — and, in the other direction, a `classifications_seen` naming two labels while `classification` reports one of them |
| `CRYPTIC_MISCLAIM` | cryptic asserted on a site present in the apo ensemble, on too few apo structures, from tier alone, with a mechanism-inconsistent potency prior, or with a mechanism asserted over an all-null census |
| `NULL_IS_NOT_ZERO` | an unexplained null, a measured zero also listed in `not_found`, a sentinel string in a numeric field |
| `AS_OF_LEAKAGE` | undated sources used under a cutoff without `leakage_risk: true`, and any date after the cutoff |
| `AXIS_CONFLICT_UNDECLARED` | the axes disagree and the dossier does not say so |
| `ASSAY_PROVENANCE_MISSING` | an actives count with no top assay named, or a dominant assay never checked for target identity |

Thresholds are named module constants, not inline numbers, so they can be argued
with: `INSUFFICIENT_ACTIVES_THRESHOLD = 50`, `SINGLE_ASSAY_DOMINANCE_PCT = 30.0`,
`CRYPTIC_APO_ABSENCE_FRACTION = 0.8`, `AXIS_CONFLICT_ACTIVES_THRESHOLD = 500`,
`DRUGGABILITY_FALSE_NEGATIVE_BAND = 0.5`, `DRUGGABILITY_FALSE_NEGATIVE_FLOOR = 0.1`,
`PRIMARY_VOLUME_CLUSTERING_D = 1.6`, `MERGED_VOLUME_A3 = 1000.0`,
`VOLUME_GUIDE_DRUGGABLE_A3 = 240.0`, `VOLUME_GUIDE_HARD_A3 = 210.0`, and
`OFF_SITE_CENTROID_DISTANCE_A = 4.0`.
`test_the_measured_constants_are_pinned` asserts **every one of them**, so
widening the false-negative band or promoting the volume guide into a classifier
fails a test that names the number that moved. That claim was false until
2026-08-15 — the test pinned six of the ten, and the four it skipped were the
policy thresholds, which are exactly the ones a reader would want to argue with.
`OFF_SITE_CENTROID_DISTANCE_A` was in neither the list nor the test, despite
being the threshold `CLAUDE.md` and `rubric.md` describe at most length.

**The two volume-guide constants classify nothing, deliberately.** They are used
in exactly one place: deciding when a low druggability and a large volume
disagree loudly enough that `tractability.caveat` must say so. No verdict, no
threshold, no gate. **They are not a boundary awaiting out-of-sample
validation** — an earlier version of this sentence said so, and that framing is
withdrawn with the rest of rule 4a. The separation they encoded is RETRACTED,
not pending: four of the five hard anchors did not measure their target. They
survive only as a *disclosure trigger*, which asserts nothing about which side
of the number a target falls on, and that is why it outlives the retraction.

### The gate's vocabulary is the TOOL's vocabulary, not a shorter one

`pocket_vs_interface.classification` has **seven** legal values, not four. Four
come from `interface_analysis.classify_pocket` — `orthosteric_candidate`,
`allosteric_candidate`, `destabiliser_candidate`, `no_partner_structure` — and
three come from the aggregation step in `modal_app.py` that runs over them:
`mixed`, `no_pocket_to_classify` and `numbering_mismatch_not_interpretable`.
The last two are abstentions and demand nothing.

The same defect had a second instance one field over: `SELECTION_BASES` was
transcribed with five of `site_pocket_selected_by`'s **six** values, so the
validator rejected `no_pocket_overlapped_ligand_site` — see step 4 above. **When
this skill and the tool disagree about a vocabulary, the tool is right**, because
the tool is what produces the value. A gate narrower than its own input does not
make the output stricter; it makes the agent launder.

**The validator used to reject `mixed`, and `pocket-scan` mandates it.** That is
not a difference of opinion about a label; it is a machine gate refusing a value
its own upstream tool is required to emit. Two symmetry copies of one ligand in
one structure can land either side of the 0.25 overlap boundary — measured on
**8DYG, ligand U5Q**: copy A `allosteric_candidate` at **0.22**, copy B
`orthosteric_candidate` at **0.36**, both flagged borderline. A caller that
reaches into `per_structure` and takes whichever copy came first is tossing a
coin between two different mechanistic claims, so the aggregation rule reports
the disagreement *as* a disagreement. A run that met the old enum had to record
the conflict in `not_found` and hide the true value in a `_consensus_note`.
**That is laundering, and it is precisely what this dossier exists to prevent.**
The enum was wrong; `mixed` is now legal.

**But admitting it is not the whole fix, because a bare `mixed` is worse than
either label it replaces** — it names no mechanism, so a reader cannot act on it
at all. `INTERFACE_MIXED_UNRESOLVED` is the price of admission. `mixed` must
carry:

| what | why |
| --- | --- |
| `classifications_seen` — at least two *distinct* classes | mixed **between what**. One class repeated is a consensus, not a disagreement |
| `pocket_interface_overlap` as the individual overlaps, not one scalar | 0.22 and 0.36 straddling 0.25 is a pocket **on the boundary**; a single 0.22 is a claim that it is not |
| `partner_pdb_id` | a classification is measured against a complex, mixed or not — rule 2b does not relax |
| `matches_mechanism_hypothesis` not `true` | a disagreement cannot confirm a hypothesis, and the copy that agrees with the prior is exactly the one that gets quoted |

So a consumer reading `mixed` acts on it like this: **both mechanistic
hypotheses are live**, the overlaps say how far apart, and the resolving
experiment is a structure that separates the copies — which is what
`next_experiment` should then name. It is not a hedge and it is not a
"medium confidence".

The rule also fires in the opposite direction. `classifications_seen` naming two
labels while `classification` reports one of them is the first-wins bug itself,
caught after the fact.

### Nobody was reading the template, so the template kept drifting

`test_template_drift.py` loads the literal JSON out of `CLAUDE.md`, fills every
leaf with a plausible value and asserts the gate returns **zero violations**.

It exists because the template and the validator disagreed **three times in one
day**, and the pair had already had this exact defect once before. Every other
test in the suite builds its dossiers from `examples/*.json` — hand-maintained
real runs, which get updated whenever the validator does. **The template was the
one artifact the suite never read**, so a validator change that did not reach it
was invisible until a live run failed.

Two halves, because they fail differently:

- **Static.** Every dotted path the validator names — from `REQUIRED_TOP_LEVEL`,
  `ENUMS`, `MEASURED_FIELDS`, `NOT_DATE_FILTERABLE`, and from every
  fully-qualified dotted string literal in the module source — must exist in the
  template. It names the missing key instead of making you read a violation. The
  source scan is not `_get(d, "…")`-shaped on purpose: three of the four keys in
  the regression it was built for are never passed to `_get` at all. What they
  all have is a `Violation` naming the path, which is a property of how this
  validator is written.
- **Behavioural.** Fill and run the gate. This is what catches requirements
  living inside a rule body (`drug.get("load_bearing")`), which no static scan of
  call shapes would see.

**A key that must genuinely live on only one side is fine — it needs an entry in
`EXEMPTIONS` with a stated reason, and a test asserts every exemption is still
needed.** Silence is not an option. There is currently exactly one.

Three things the file is careful about, all of which cost a debugging cycle:

- **The fill table encodes semantics, never structure.** A new template key that
  needs no special value is handled by `_default_for` and needs no edit. If the
  table had to name every key it would be a second copy of the template, free to
  drift on its own — so `test_every_FILL_key_is_a_real_template_path` refuses to
  let it reference a path the template does not have. That test caught a stale
  entry the first time it ran.
- **A list of literals is not a list stub.** Filling `clustering_d_swept:
  [1.6, 2.4]` by taking element zero turned the mandatory two-value D sweep into
  a single D, and the gate correctly called it a coin flip.
- **The placeholder string contains no substring of any `MEASURED_FIELDS` path.**
  `_mentions_field` matches by substring, so a careless placeholder dropped into
  `not_found` would silently excuse a null somewhere else in the dossier.

### A rule that asks "was this attempted?" must test the VALUE, never the key

Every time the template grows, any guard written as "is this key absent?" dies
silently — the template now ships the key, the branch becomes unreachable, and
the rule stops catching the thing it exists for **without a single test
failing**. It happened here: `CRYPTIC_MISCLAIM` keyed on `cryptic_evidence`
being absent, the template started always shipping that block, and a dossier
asserting `cryptic_mechanism: subunit_occlusion` with the census left entirely
null passed clean. The tests did not notice because they built the broken case
with `del`, which no template-derived dossier ever does.

So: test `ev.get("is_cryptic") is None`, or `_has_reported_value(ev)`, never
`isinstance(ev, dict)`. Keep the `isinstance` check only as a disjunct guarding
against `None.get`, because a hand-written dossier can still write
`cryptic_evidence: null`.

**And when you add a test for a sparse dossier, build the sparse case the way
the template would produce it — an all-null block — not with `del`.** A `del`
test passes against a dead guard and tells you nothing.

The inverse is safe and worth knowing: a key-absence check used to *skip* a
check (`check_null_is_not_zero`'s `__absent__` continue, `check_fraction_carries_n`'s
`isinstance` early return, the `ENUMS` absence skip) becomes unreachable in the
same way, but that makes the rule run *more* often, not less. Those are fine.

`test_validate_dossier.py` carries a deliberately broken dossier for every rule.
A validator that passes everything is worthless, so that file is part of the
deliverable, not an extra.

## Worked examples

`examples/jak1_P23458.json` and `examples/tnf_P01375.json` are both real measured
runs and both return zero violations. Read them together — the pair is the point.

|  | JAK1 P23458 | TNF-alpha P01375 |
| --- | --- | --- |
| approved small molecules | **9** (ruxolitinib 2011 to deuruxolitinib 2024) | **0** |
| approved biologics | 0 | **5** |
| distinct compounds | 14,472 | 2,582 |
| top assay share | 4.8%, clean | **45.0%, measures IRAK4** |
| best potency | 0.010 nM, characterised enzyme assay | 0.03 nM **uncharacterised** — rejected; reported 1.3 nM SPR instead |
| structures | 42 holo / 10 apo | 17 holo / 35 apo, the one holo ligand a known frequent hitter |
| pocket volume | **305.9 - 913.8 A^3**, spread 66.5% | **refused** (would pool to 126.9 - 809.5) |
| volume at D=1.6 — **the primary number** | **not separated**, `not_found` (the run pooled both D values) | **refused**, same grounds |
| druggability | **0.020 - 0.437**, 21.8-fold, `load_bearing: false` | **refused** (would pool to 0.001 - 0.651, fold-range 651.0) |
| site basis | **4 of 4 `ligand_site_jaccard`** | 2 of 12 `ligand_site_jaccard`, **10 of 12 `site_signature_unreliable_homooligomer`** |
| same-site control `max_radius_difference_a` | **2.29 / 2.16 A** | **14.43 / 16.49 A** |
| mdpocket off-site distance | **1.86 A** — check passes | **29.57 A** — `off_site_warning` raised |
| ensemble | 3EYG + 10PI, 2 structures / 4 measurements | 2AZ5 + 5 apo, 6 structures / 12 measurements |
| verdict | `small_molecule_tractable`, basis `retrieved_precedent` | `small_molecule_tractable`, basis `both`, **`axis_conflict` populated** |

**TNF has none of JAK1's precedent, and its computed axis is refused rather than
low — that is what earns `axis_conflict` here.** An earlier version of this line
read "TNF beats JAK1 on every pocket metric", and it is **withdrawn**: the table
above refuses TNF's volume, D=1.6 primary and druggability outright, and its
same-site and off-site controls (14.43/16.49 Å, 29.57 Å) are *worse* than
JAK1's, so "every metric" was false against the table it sits under. The
cross-target ranking it implied also has no source left — it rested on the
calibration retracted in rule 4a, and TNF's own 207 Å³ anchor was one of the
compromised ones. The honest statement is narrower and still enough: a target
with zero approved small molecules is not thereby untractable, and averaging the
two axes would let the retrieved zero erase a computed axis that was never
measured. That is the inversion this agent exists to prevent.

**The most instructive thing in the pair is that one fills and one refuses, from
the same tool on the same day.** Both dossiers previously reported no pocket
geometry at all, for a reason worth keeping: each run had relayed a *maximum*
druggability (JAK1 0.735, TNF 0.992) and no minimum, and a maximum without its
minimum is a point estimate wearing a range's clothes. That refusal was correct
on the data it had. **It has since been superseded on JAK1 by measurement, and
it has been re-earned on TNF for a better reason.** Note that 0.735 does not
appear anywhere in JAK1's new 0.020-0.437 range: the old figure was one end of
an unknown interval on an unnamed ensemble, which is exactly why it was not
reportable.

**JAK1 fills** because all four measurements were selected by
`ligand_site_jaccard` — real overlap with the MI1 site, 0.300 to 0.778 — on two
entries that are not oligomeric, with the frame-independent same-site control at
2.29 A. The spread is a spread over one pocket.

**TNF refuses** because ten of its twelve measurements came back
`site_signature_unreliable_homooligomer`. That is not a tool failure; it is the
tool declining to claim a site, and it is right to: the signature discards chain
identity, and 2AZ5's 22 site residues collapse to 14 distinct numbers across
four identical chains, so any pocket carrying those numbers matches. **Pool them
anyway and the druggability fold-range comes out at 651.0** — the withdrawn
"651-fold spread across five apo TNF-alpha structures of the same site"
regenerating from the identical defect, on a fresh run, five years of lessons
later. `test_the_real_tnf_pool_is_rejected` pins that exact case.

Three independent controls agree with the refusal, and this is why the rule asks
for more than one: `max_radius_difference_a` at 14.43 and 16.49 A against JAK1's
2.29 A and the 16.61 A that invalidated the IL-17A pool; per-structure centroid
drift of 18.43 A on 2AZ5 between the two clustering values alone; and mdpocket's
ligand-free `site_from_density` centroid sitting **29.57 A** from the ligand,
nearly four times the 7.73 A error that forced the original retraction. The
tool's own `off_site_warning` names it.

**A refusal is not a low score and must never be read as one.** On a target with
zero approved small molecules, an empty `pocket_druggability` block is the one
place a reader is most likely to supply a pessimistic number from imagination.
TNF's `tractability.caveat` says so in as many words.

Two smaller things the runs corrected, both left visible rather than smoothed
over: the in-repo calibration said the TNF site was recovered in **all five** apo
structures, and at the transferred ligand site 1TNF returns **0.00 A^3** — so
`n_apo_site_absent` is 1, not 0, and the call is unchanged but the census is
honest. And neither `disorder_fraction` is reported despite metapredict
returning one, because no accession was passed and it read the sequence off the
structure — 280 of JAK1's 1154 residues, 141 of TNF's 233.

**Both examples now decline the primary number, and that is the rule 4a
re-examination flag landing on our own exemplars rather than on somebody else's
dossier.** Neither run separated the D=1.6 site volume — JAK1 pooled 3EYG and
10PI across both clustering values into 305.9–913.8 Å³ and did not keep the
per-measurement figures, so `primary_d1_6_a3` is null with a `not_found` line
saying exactly that. Under the old ordering that was a complete answer. It is
not one now, and the honest record of the change is a declined field rather than
a number back-derived from a spread. Note the shape of the upper end: **913.8 Å³
is what a D=2.4 merge looks like** (rule 4: above ~1000 Å³ sites have merged), so
reading the spread as "the site is up to 914 Å³" would be reading a merge.

Neither example's verdict moves, and that is worth saying plainly: JAK1 is
carried by retrieved precedent — nine approved small molecules — and TNF by both
axes with `axis_conflict` populated. Nothing here was resting on a druggability
score, which is why the demotion cost these two dossiers a field and not a
conclusion. A dossier that *would* have moved is exactly the one the new rule
fires on.

JAK1's `not_found` also carries a smaller, sharper case: the run counted **9**
approved small molecules after collapsing salt/parent pairs, and only **8** could
be named from retrieved rows. The count and the list disagree by one,
deliberately, with the gap written down. Guessing the ninth would have made the
dossier look complete and been the worse output.

## Failure modes

### The verdict field is where averaging hides

The rule "never average the axes" is easy to keep at the level of fields — you
simply do not create an `overall_score`. It is hard to keep in the verdict,
because a single label over two disagreeing axes *is* an aggregation unless you
say which axis produced it. That is what `verdict_basis` is for. A
`small_molecule_tractable` with no basis and a populated `axis_conflict` is an
average with extra steps.

The subtler form: writing a hedged verdict *because* the axes disagree —
downgrading TNF to `insufficient_evidence` to split the difference. That is
averaging in the verdict field. Report the strong pocket, report the absent
precedent, populate `axis_conflict`, and let the reader carry both.

### "Approved drugs exist" is the single most expensive sentence here

IL-17A has three approved antibodies all reporting `action_type: INHIBITOR` and
zero approved small molecules. TNF-alpha has five approved biologics and zero.
A dossier that reports "approved drugs exist" for either is wrong in the way
that matters most, and it is wrong in a direction that looks like good news.

The validator catches the crude form (a `-mab` name in the small-molecule list,
or the same name in both blocks). It cannot catch the form where the modality
was never read and the biologics silently became the evidence for
`verdict_basis: retrieved_precedent` — except through the clause that fires when
tractability is claimed on precedent with zero approved and zero clinical small
molecules and no characterised potency. Read `molecule_type` per drug. There is
no shortcut.

### A row count is not a drug count

JAK1's 11 approved rows are 9 approved drugs, because salt and parent forms are
distinct `molregno`s. Deduplicating on `molregno` does not deduplicate drugs.
Either collapse the pairs or say the figure is a row count — and if you collapse
them and cannot itemise the result, say that too, as the JAK1 example does.

### An actives count is a claim about assays until you check

TNF-alpha's 6,447 activities include 2,901 — 45% — from an IRAK4 monocyte assay
that measures a different protein. The obvious defence fails: 90% of TNF's
activities are labelled `assay_type = 'B'` and the contaminating assay is one of
them. Only the description separates a direct binding measurement from a
cellular readout.

So the order in step 2 is not stylistic. Group by assay *before* writing
`distinct_actives`, or you will write a number you then have to retract three
fields later.

### A maximum is not a range, and this is the failure the real runs actually had

Earlier runs delivered `max druggability` and no minimum. The temptation is
overwhelming: you have a number, the field wants a number, the number is even
defensible. It is still a point estimate, and the whole reason druggability is
reported as a range is that a single value is a coin flip — pinning `-D 1.6`
gives TNF-alpha 0.002 at the site of a co-crystallised 570 Da ligand, and the
same site at `-D 2.4` scores 0.346.

If you have one end of the range, you do not have the range. Put the value in
`not_found` with its provenance and say plainly that it must not be lifted into
the body.

**The 2026-08-15 re-run reproduced that coin flip almost exactly.** 2AZ5's site
pocket scores **0.003 at D=1.6 and 0.587 at D=2.4** — the same false negative at
the same co-crystallised ligand, two hundred-fold apart, from the clustering
parameter alone. Whatever else is uncertain about the druggability score, its
sensitivity to `-D` is now measured twice on the same structure.

### The next failure after "I have both ends" is "both ends describe one site"

Having a min and a max is necessary and it is not sufficient, and this is the
trap the TNF re-run walked into and the dossier caught. A complete, successful
run with every stage `ok` produced a full range on six structures — and it is
not reportable, because `site_pocket_selected_by` says ten of the twelve values
are not attributable to a site. **Read that field before you read the numbers
beside it.** The clue that something is wrong is not in the range; the range
looks great. It is in `max_radius_difference_a`, in the per-structure centroid
drift, and in `distance_to_donor_ligand_centroid_a`.

The tell that you got this wrong is a fold-range in the hundreds. That is what a
homo-oligomer's protomers produce when residue numbers are matched across them,
and it is what the retracted 651-fold figure was.

### Reporting the number you have instead of the number that matters

Volume is the reproducible measurement; druggability is a weak prior from a
three-descriptor logistic regression fitted on 21 positives. When a run returns
druggability and no volume, the pull is to report druggability, because
something is better than nothing.

It is not. The validator requires that druggability reported without volume be
accompanied by a `not_found` line naming the gap, and that is the *minimum* — the
honest response is to re-run for volume, which is why it is the JAK1 example's
`next_experiment`.

**As of 2026-08-15 the case against druggability is measured; the case *for*
volume as a discriminator is retracted.** The two halves of this table have
different standing and must not be read as one result. Over 15 targets, 67
structures and 134 measurements:

| | |
| --- | --- |
| volume at D=1.6, target level | ~~**AUC 1.000**, CI [1.000, 1.000], stable under all 15 refits~~ — **RETRACTED 2026-08-15, rule 4a.** The anchors did not measure the proteins they were attributed to; the CI was degenerate by construction. **Volume carries no verdict.** See the retraction box above |
| druggability at D=1.6, target level | **0.720**, 95% CI **0.44–0.94** — includes chance, P(AUC≤0.5) = 0.071 |
| druggability at D=2.4 | **0.520** — chance |
| druggability on certain positives | ligand-anchored holo structures across the druggable targets: median **0.320**, and a large majority below 0.5. **Denominator under audit** — it was quoted as 37, one of which (RORgt 6C1P) is now void, and the remaining 36 have never been checked at residue level. See the denominator note below |
| persistence | site pocket detected in **100% of structures for all 15 targets** → **AUC 0.500** |
| published consensus criterion | **AUC 0.560**, and it **ranks MYC top at 0.80** |
| co-crystal flag alone, no structural measurement | **AUC 0.900** — the label and the measurement share a cause, which is why the volume separation above could not be evidence |

Two things follow that are easy to get wrong in opposite directions.

**The score is inverted at target level, not merely noisy.** MYC — zero holo
structures, canonical undruggable — has a D=2.4 median of **0.75**, above KRAS
0.54, BCL-2 0.52, JAK1 0.49, EGFR 0.44 and NLRP3 0.12. And the *clustering
choice* does about 1.5× more work than the biology: median within-structure
|D=2.4 − D=1.6| is **0.229** against a between-group difference of medians of
**0.154** at D=1.6. D=2.4 is not a fix — IRAK4 2O8Y goes 0.791 → 0.001.

**The obvious replacement is worse.** Persistence is *exactly* chance and the
consensus criterion built on it puts the canonical undruggable target first.
Reaching for it would reproduce the same inversion one rung down. Keep the
consensus fraction for what it does do — stopping you quoting your best
conformer — and read no tractability from it.

**The class caveat that used to sit here is now moot, and the mitigation that
answered it is RETRACTED.** The caveat was **n = 5 hard targets**, all
PPI/cytokine/membrane class, against a druggable set enriched in kinases,
nuclear receptors and GPCRs, so volume might have been tracking target class.
The mitigation offered was that the two least classical druggable targets —
NLRP3 at 242 Å³ and IL-17A at 250 Å³ — still landed above every hard target.
**That mitigation is void**, because "above every hard target" is a statement
about the hard anchors, and four of those five did not measure their target at
all. It was never a control, and it is now not even a mitigation. The class
confound is subsumed by a larger one: a co-crystal flag alone separates the
groups at **AUC 0.900**, so the label and the measurability were the same
variable throughout.

**The denominator note.** The label-free test above was reported on **37**
ligand-anchored holo structures — "certain positives by construction", holo
structures with a drug-like ligand physically bound and the pocket anchored to
it. **RORgt 6C1P, quoted elsewhere as "0.009 at rank 55 of 60", is neither**: the
entry contains no RORgt (its sole entity is an ion transport protein) and its
"ligand" is CHAPSO, a detergent. So the denominator is **36, not 37** — and **the
other 36 have never been audited at residue level**, which is the audit that
caught this one. Do not quote "15 of 37", "25 of 37" or "41%": each rests on a
denominator now known to contain at least one error and unverified in its
remaining 36 entries. **Quote the direction and the named cases, not the rate,
pending a residue-level audit of all 36.** The demotion of druggability itself is
unaffected — it does not turn on that case, and MYC's inversion was reproduced
independently.

### Laundering nulls through a vague `not_found` line

`not_found: [{"field": "tractability", "reason": "we did not run it"}]` would, in
a naive matcher, excuse every null in the block. The validator matches on the
full dotted path or the leaf name and deliberately does **not** accept the
enclosing block name, so one vague line covers nothing.

Write one line per missing field, naming the field, saying why it is missing and
what its absence does to the conclusion. "Not computed — no disorder predictor
was run. This is a missing measurement, NOT a finding that JAK1 is ordered" is
the shape.

### A zero is a result and must not be hidden

The other direction of the same rule. `disorder_fraction: 0.0` and
`disorder_fraction: null` are different claims, and a field that is both reported
as zero and listed in `not_found` is claiming both. A volume of **0.00 cubic
angstroms is a result, not a failed run** — it is the one output that cannot be
an over-claim, and substituting the nearest pocket that does have volume is how
`pocket-scan` got a headline finding wrong.

### Cryptic and occluded are different, and the difference is checkable

TNF-alpha's axial site is recovered in **four of the five** apo structures once
the third subunit is removed, at roughly 280 to 550 cubic angstroms, with
roughly 1.6 A of backbone C-alpha displacement. (**Not all five** — this said
"all five" until 2026-08-15, inheriting the in-repo calibration; measured at the
transferred ligand site, 1TNF returns **0.00 A^3**, so `n_apo_site_absent` is 1.
The call is unchanged, because one of five is far below the "all or nearly all"
the definition requires, but the census is 4 of 5 and is written down rather
than smoothed over. The same correction is recorded above and in
`tnf_P01375.json`.) It is **occluded, not cryptic**, and citing it as
a cryptic-pocket case is an error a reviewer finds immediately. KRAS switch-II —
absent on apo, druggability 0.000, backbone displaced roughly 8.8 A — is cryptic.

Two significant figures on both, deliberately. The raw volumes are 281.8-546.0
A^3 but fpocket estimates volume by Monte Carlo, so the fourth digit is noise.
And 8.83 A / 1.62 A are **hand-calibration** figures from a protocol that
disabled auto-trim and residue-name matching and named the mobile regions by
hand; the deployed zero-knowledge default returns **8.65 A** and **~1.55 A**.
Mechanism and `is_cryptic` are identical under both, so nothing downstream of
the label changes — but quote what the run reported and say which protocol it
was, and never present 8.83 or 1.62 as pipeline output.

The 2026-08-15 deployed run bears that out and sharpens it: over the five apo
entries it returned per-structure displacements of 1.55, 1.61, 1.75, 1.83 and
1.82 A. The **1.55 A is the rule 5 deployed figure reproduced exactly**, and the
1.83 A the TNF example reports is the *maximum over a five-structure ensemble* —
a different quantity from either single-pair number, which is why that field now
carries its protocol beside it.

The definition is not ours: cryptic means absent in all or nearly all unbound
structures. A site missing from one apo structure but present in others is
low-scoring, not cryptic, which is also the argument for running an ensemble at
all — a single apo structure cannot distinguish "absent" from "low-scoring in
this crystal form".

### Cryptic risk flagged from tier carries no information

Setting `cryptic_pocket_risk: high` because the tier is apo fires on every apo
target equally. It is a restatement of the tier, not a measurement. Measure the
displacement, or say `undetermined` and put the reason in `not_found`.

And classify on displacement, not on clash composition. KRAS's switch-II loop
moves roughly 8.8 A (hand calibration; 8.65 A on the deployed default) and
**zero** of the 12 clashing atoms at 2.0 A are backbone — a loop that swings
that far carries its side chains with it. Keying on clashes inverts the
answer on the canonical case.

### Pooling a spread over a site that cannot be resolved

A spread is only a measurement if every value describes the same site. On an apo
structure with no ligand to anchor to, `pocket_scan` falls back to "the most
druggable pocket anywhere in the chain" and reports
`site_pocket_selected_by: max_druggability_no_ligand_site`. Pooling those across
an ensemble compares different pockets on different parts of the surface.

On a homo-oligomer it is worse than a fallback, it is impossible in principle:
the signature is a set of residue numbers with chain identity discarded, and a
homotrimer's protomers triplicate every number. A 19-residue reference collapses
to 11 distinct numbers. That is how the **withdrawn 651-fold TNF-alpha
druggability spread** was produced — the matcher was tracking a pocket 7.7 A from
the site it claimed, with 12.2 A of internal inconsistency between structures.
Do not cite that figure. Fixing the site by construction with mdpocket instead of
by post-hoc matching cut the measured CV from roughly 28% to roughly 10%, an
inflation of about 2.8-fold.

Two significant figures, never three. fpocket estimates volume by Monte Carlo
and mdpocket inherits it: three identical reruns of one five-structure ensemble
gave CVs of 12.1 / 11.3 / 10.8% against the deployed run's 9.9%, so about **one
percentage point of any CV is the estimator's own noise**. The improvement is
real and survives that noise; the third digit does not exist, and a CV difference
smaller than ~1 pp is not a difference between sites. An earlier version of this
paragraph read "27.8% to 10.2%".

Note also what that CV was measured on: `site_from_density`, whose centroid sits
7.7 A from the transferred SPD304 ligand. It is a real measurement of
reproducibility across the ensemble. It is not a measurement of the SPD304 site.

`651`, not `650`: `fixtures/pocket_calibration.json` records `withdrawn_min`
0.001 (2ZJC) and `withdrawn_max` 0.651 (1A8M), and 0.651/0.001 is exactly 651.
The `650` in older prose was a rounding with no separate measurement behind it.
Cite 651 so the figure resolves to the record that documents its retraction.

A pocket-matching step is itself a measurement and needs its own control: report
the matched centroid distance across the ensemble, not an overlap fraction. Two
pockets sharing residue numbers can be 12 A apart and no overlap score will say
so.

### Silent as-of contamination is worse than a gap

Under a cutoff, a source that cannot be date-filtered must be omitted or carried
with `leakage_risk: true` and a note naming the source. The one that catches
people is clinical candidates: ChEMBL's `max_phase` is a *current* value with no
phase history, so neither the presence nor the absence of a clinical candidate at
a past date is retrievable. It always carries the flag under a cutoff, even when
the list is empty — an empty list under a cutoff is also an unverifiable claim.
Bioactivity counts have the same problem for a different reason: no date column.

The reason to be strict is that the whole point of an as-of run is retrospective
evaluation, and a retrospective evaluation contaminated by future data is
worthless. RA at 2010-12-31 is a biologic-only disease across every target in the
fixture set; a run that quietly returns ruxolitinib has not answered the question
asked.

### `insufficient_evidence` feels like failure and is not

IL-11 is shaped to tempt: 15 compounds, a 140 nM Kd that looks real, 8
structures. All 15 come from one assay and none of the structures is holo. Any
number here is a failure, and the hardest thing to make a system do is decline.

The validator enforces the floor mechanically — thin chemistry plus zero holo
plus zero approvals must return `insufficient_evidence` — but the floor is low
and most declines will be judgement calls above it. A confident score on an
unstudied target is the worst output this agent can return.

Note the deliberate asymmetry: `holo_count: 0` triggers the rule and
`holo_count: null` does not, because a failed structure query is not a finding of
zero structures. That is the null-is-not-zero rule doing real work.

### Clinical failure read as evidence against tractability

They are different questions and other stations answer the second. RORgt is
tractable and clinically failed. Never lower a tractability number because
programs failed; record the terminations in `terminated_programs` with their
stated reasons and let the reader weigh them. The corresponding validator
behaviour is deliberate: a `small_molecule_tractable` verdict with zero
approvals and a full `terminated_programs` block raises nothing.

### Recall dressed as retrieval

The most likely way this dossier acquires a false fact is a plausible identifier
written from memory — a PDB ID that "looks like" a JAK1 entry, a Pfam accession,
an approval year. Everything in the fixtures was retrieved and the ones that
were not are marked `NOT_FOUND` or carry an `OPEN_CAVEAT`.

The JAK1 example leaves `ensemble_pdb_ids` empty rather than name two plausible
entries, and leaves the ninth approval unnamed. An empty list with a reason is
worth more than a filled one that is 90% right, because the reader cannot tell
which 10% is wrong.

Related: accession mapping is not sufficient identification. Several
P10415-mapped PDB entries are actually Bcl-xL constructs. Check the entry title.

### Prose is a laundering channel the validator cannot see

`NUMBER_WITHOUT_PROVENANCE` walks numeric leaves. A number written into
`axis_conflict` or `caveat` as text is a string and passes untouched. That is by
design — prose has to be able to quote figures — but it means the discipline in
narrative fields is yours alone. When a number appears in prose, attach where it
came from in the same sentence, and never use prose to state a figure the
structured fields refused.

### `axis_conflict` as filler

The mirror of under-reporting. Populating it on every target with any tension
makes it noise, and the validator deliberately does **not** require it for CD20 —
three approved antibodies, a four-pass transmembrane protein whose epitope is a
small extracellular loop, nothing for a small molecule to bind. Both axes agree.
Silence is the correct output.

### Passing the validator is not being right

It checks structure and self-consistency. It cannot tell you that
`distinct_actives` is wrong, that the top assay was misread, or that the pocket
you scanned is not the site anyone cares about. A dossier that passes may still
be wrong; a dossier that fails is wrong in the way the violation names. Do not
read a clean run as verification of the content.

## Output

A single JSON object matching the template in `CLAUDE.md`, with the extension
keys above. It goes to **two places, every run**:

1. **Write it to `/mnt/session/outputs/druggability-dossier.json`** — that exact
   path and no other, creating the directory first if it does not exist.
2. **Paste the complete JSON into your final reply.**

**Neither substitutes for the other, and this file used to say "Return the JSON
and nothing else", which was wrong.** The file is the only channel the dossier
reaches the *grader* by — a dossier that exists only in the reply is ungraded and
fails every criterion. The reply is the only channel it reaches a *human* by,
because sandbox files are not retrievable through the Files API once the session
ends. A short wrap-up message in place of the JSON loses the deliverable. The
deployed prompt in `CLAUDE.md` requires both; that is the authority here.

Validate the file you wrote, not a copy of it:

```bash
python3 .claude/skills/assemble-dossier/validate_dossier.py \
        /mnt/session/outputs/druggability-dossier.json
```

Before returning: `validate_dossier.py` returns zero violations. If it does not,
fix the dossier. If a rule seems wrong for a genuine case, say so in the run's
notes rather than editing the rule — the rules encode findings that cost
measurements to establish, and one of them was already withdrawn once for being
wrong in exactly that way.
