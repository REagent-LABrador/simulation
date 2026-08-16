# Interpretability Panel — the workflow reasoning trace over a dossier

A specification for a panel that, when a user clicks one prediction (one
druggability dossier), shows **how that verdict was reached from the evidence**:
which pipeline stages ran, in what order, what each computed, every intermediate
number, where each number came from, and the caveat that number cannot be shown
without.

The panel renders exactly one artifact — a dossier JSON matching
`schema/output.schema.json` — and claims nothing that is not traceable to a field
in it.

---

## 0. What "reasoning trace" means here

**It is the WORKFLOW trace, not an LLM chain-of-thought.** This station is a
pipeline of retrieval and computation stages (structure selection, pocket
scanning, precedent lookup, falsification). The "reasoning" the panel exposes is
the flow of evidence through those stages and the arithmetic that turned it into
a verdict — not the internal deliberation of any model. The panel never shows,
infers, or invents model tokens, hidden thoughts, or a probability the pipeline
did not emit.

The dossier is the trace. Every top-level and nested field is an **intermediate
result** the panel can surface and attribute. The panel's job is to make that
trace legible and honest to a human reader.

---

## 1. Purpose & non-goals

### Purpose

- **Explain the workflow.** Show the ordered/branching stages that produced this
  dossier and what each one established.
- **Surface every intermediate number with attribution.** Volume, druggability
  rank, potency, actives counts, displacement, distances — each shown with its
  stage, its raw value, and its provenance (`source`/`pdb_id`/`chembl_target_id`/
  query).
- **Carry the caveat with the number, always.** A number shown without its caveat
  is the exact failure the dossier design fights. The panel enforces this at
  render time: a field that has a required caveat cannot be displayed without it.
- **Make refusals and nulls first-class.** A `not_found` reason is displayed
  verbatim. A null is never rendered as a zero.

### Non-goals — what the panel explicitly does NOT claim

- It is **not** an LLM introspection / model-confidence panel. No chain-of-thought.
- It does **not** produce a probability of druggability. The fpocket
  druggability score is a within-structure rank, not a probability (rule 4.0 /
  4b provenance); the panel says so.
- It does **not** average, merge, blend, or reconcile the two axes into a single
  score. The retrieved-precedent axis and the computed-tractability axis are
  rendered **side by side and never merged** (CLAUDE.md "The two axes"; rules 9,
  9b). There is no overall number to show, and the panel must not manufacture one
  (no colour-scale, no combined bar, no weighted blend).
- It does **not** decide, rank targets, design molecules, or assess biologics as
  small-molecule evidence. It reports what the dossier reports.
- It does **not** re-run or re-compute anything. It is a read-only renderer of a
  finished dossier.

---

## 2. The trace model — the stages the panel walks

The panel presents the pipeline as an ordered, branching sequence of stages. Each
stage maps to a skill under `.claude/skills/` and to specific dossier fields.
For each stage the panel shows: **inputs → dossier fields produced → the one-line
"what this establishes"**, and a status (`ran` / `refused` / `not_run` / `null —
no tool`).

### Stage order

```
intake
  → structure selection  ── experimental (holo → apo)
                          ── homolog transfer (fallback)
                          ── ESMFold predicted-structure fallback (last resort)
  → pocket scan          ── fpocket/mdpocket · disorder · cryptic · pocket-vs-interface
  → precedent lookup     ── target · family · structural-neighbour · pocket-neighbour
  → falsification
  → assembly + validation
```

### Stage table

| # | Stage (skill) | Inputs | Dossier fields produced | What this establishes |
|---|---|---|---|---|
| 1 | **intake** (`graph-intake`) | caller `{task}` prose: accession, `as_of_date`, `disease_context`, `interaction_to_disrupt`, `mechanism_hypothesis` | `input.*` (verbatim echo), `target.*`, `as_of_date`, optionally `follow_up_questions` | What was actually asked. The `input` block is echoed exactly, never inferred; `target.uniprot_accession` is the *resolved* accession, `input.uniprot_accession` is what the caller said. The answer is about a **target and a mechanism**, keyed on `(input.uniprot_accession, input.mechanism_hypothesis, input.as_of_date)`. |
| 2a | **structure selection — experimental** (`structure-select`) | resolved accession, `as_of_date`, `mechanism_hypothesis` | `structure.tier` (`holo_experimental`/`apo_experimental`), `structure.pdb_id`, `bound_ligand.*`, `total_pdb_structures`, `holo_count`, `apo_count`, `ensemble_used`, `structural_neighbour_precedent.*` | Which real structures the scan runs on. Selection order is strict: holo (drug-like ligand bound) → apo → predicted. Predicting a structure that exists in the PDB is a defect (rule 2). Classifies holo vs apo by ligand chemistry, not label. |
| 2b | **structure selection — homolog transfer** (`structure-select/homolog_transfer.py`) | fold neighbours carrying a drug-like ligand | `structure.transferred_homolog_site` (donor PDB, ligand, TM, RMSD, aligned length, clash + backbone-clash counts) | A site basis for a structure-less-but-not-orphan target — only past three guards (donor ligand passes `ligand_filter`; contact shell overlaps the aligned region; TM ≥ 0.5, RMSD ≤ 5 Å, zero backbone clashes). A perfect transfer of the wrong ligand is a wrong answer with high confidence (rule 4b). Unreachable if nothing is liganded in the neighbourhood. |
| 2c | **structure selection — predicted fallback** (`structure-select/predicted_structure_fallback.py`, ESMFold) | accession sequence, no experimental structure **and** no usable homolog | `structure.tier = "predicted"`, `structure.predicted_plddt` | The one GPU exception. **ALWAYS returns a result WITH WARNINGS, never null and never zero** (rule 4c). Low pTM ⇒ louder warnings + lower confidence, not refusal. The only null is a *tool failure* (`status: not_run`). *Note: in the current deployment this and all cofold/affinity tiers are nulled — see §6 and rule 13.* |
| 3 | **pocket scan** (`pocket-scan`) | selected ensemble, `chains`/`site_residues`, clustering swept over {1.6, 2.4} | `tractability.pocket_volume_a3.*`, `pocket_druggability.*`, `site_pocket_rank.*`, `ensemble_consensus_fraction.*`, `disorder_fraction`, `cryptic_*`, `pocket_vs_interface.*`, `mdpocket_site_definition_used`, `site_centroid_to_ligand_distance_a`, `site_hypothesis_basis`, `tractability.method.*` | The computed axis. Whether a small molecule *could* bind, from geometry — with its declared blind spots (cryptic collapse, within-structure normalisation, clustering sensitivity). Also runs disorder, cryptic-mechanism classification, and interface classification. Emits refusals when the site cannot be identified across the ensemble. |
| 4 | **precedent lookup** (`precedent-lookup`) | resolved accession, `as_of_date` cutoff, `disease_context` | `target_precedent.*`, `biologic_precedent.*`, `family_precedent.*`, `pocket_neighbour_precedent.*` (`structural_neighbour_precedent` from stage 2a) | The retrieved axis — the stronger axis when it exists. Four separate sub-axes, never merged: `target_precedent` (measured on this protein), `pocket_neighbour_precedent` (strongest transfer), `structural_neighbour_precedent` (fold), `family_precedent` (weakest). Modality is classified first; only small molecules count toward `target_precedent`. |
| 4b | **terminated programs** (`terminated-programs`) | drug codes, ClinicalTrials.gov | `target_precedent.terminated_programs[]` | Why programs stopped — registry reason beside the literature account, contradictions shown side by side. Clinical failure is **not** evidence against tractability (rule 7). |
| 5 | **falsification** (`falsification-sweep`) | the assembled precedent + pocket claims | `falsification.checks_run[]`, `falsification.findings[]`, `falsification.survived` | Attacks the claim before it ships: single-assay dominance, frequent-hitter ligand, one-crystal-form pocket, terminated programs. Records checks that found nothing as well as checks that found something. Does not produce a verdict; attaches evidence. |
| 6 | **assembly + validation** (`assemble-dossier`) | all of the above | `verdict`, `verdict_basis`, `axis_conflict`, `next_experiment`, `not_found[]`, and the full validated JSON | The verdict and *which axis carried it*. Never averages the axes; when they disagree, populates `axis_conflict`. Gated on `validate_dossier.py`. |

**Branch semantics the panel must render.** Stages 2a→2b→2c are a fallback
chain: exactly one tier wins, and the panel shows which branch fired and why the
earlier ones did not. Stage 3's sub-branches (cryptic, interface) may each be
`not_run` with a stated reason (e.g. cryptic needs both an apo and a holo frame —
JAK1's cryptic stage returned `not_run` on an all-holo ensemble). Stage 4's four
sub-axes are independent and any may be `null`/not-run.

---

## 3. Field-to-explanation map

For each field: plain-language gloss, the caveat that **must** be shown alongside,
and which axis it belongs to — **P** = retrieved-precedent axis, **C** =
computed-tractability axis, **N** = neither (verdict/provenance/meta).

### Verdict layer

| Field | Gloss | Caveat that must travel with it | Axis |
|---|---|---|---|
| `verdict` | One of `small_molecule_tractable` / `not_tractable` / `insufficient_evidence`. | This is small-molecule tractability only — not whether to pursue the indication, not clinical success. `insufficient_evidence` is a correct answer, not a failure (rule 11). | N |
| `verdict_basis` | Which axis carried the verdict: `retrieved_precedent` / `computed_tractability` / `both` / `none`. | A verdict with no basis over two axes allowed to disagree *is* an average with extra steps. Render the basis next to the verdict, always. A `computed_tractability`/`both` verdict can never be carried by druggability alone (`load_bearing: false`, rule 4.0/11). | N |
| `axis_conflict` | Free text stating that the two axes point in different directions and why. | Show it **prominently** and near the top when non-null. The disagreement is usually the most informative thing on the page (rule 9). Do not resolve it in the UI; display it. | N |

### Retrieved-precedent axis (P)

| Field | Gloss | Caveat that must travel with it | Axis |
|---|---|---|---|
| `target_precedent.distinct_actives` | Count of **small-molecule** actives measured on this protein. | Small molecules only (modality classified from structure, not `molecule_type`, which abstains on ~59% of compounds — rule 1b/1c). A raw count measures assays, not the target: pair with `assay_concentration.top_assay_share_pct` and `measures_a_different_target` (rule 6). Every count is a lower bound unless reconciled against a `COUNT` (rule 14). | P |
| `target_precedent.compound_modality_split` | Full per-modality breakdown of the compound set. | The pool is often mixed; a potency over a mixed pool is not a small-molecule claim (rule 1d). | P |
| `target_precedent.modality_unknown_count` | Compounds with no parsable structure. | Counted and disclosed, never folded into small molecules (rule 1c). | P |
| `target_precedent.best_potency_nm` + `best_potency_modality` + `best_potency_characterised` | Strongest potency, its modality, and whether the assay was characterised. | A potency without a modality is not attributable; an uncharacterised assay is unusable however good the number (rule 6). Show all three together. | P |
| `approved_small_molecules_count` vs `approved_small_molecules[]` | Collapsed drug count vs the drugs actually named. | The two may legitimately disagree (salt/parent collapse, rule 1e); when they do, the gap is in `not_found` and the panel shows both plus the `not_found` reason — never silently reconciles. Only modality `small_molecule` is legal in this list. | P |
| `clinical_stage_small_molecules[]` | Small molecules in the clinic, with phase. | Under an `as_of_date`, this needs an `as_of_leakage` entry *unconditionally*, even when empty — `max_phase` is a current value with no history (rule 8). | P |
| `terminated_programs[]` | Programs that stopped, with stated reason. | Clinical failure is **not** evidence against tractability (rule 7). Display as context, never as a negative on the axis. | P |
| `as_of_leakage[]` | Per-field flags where a source could not be date-filtered under the cutoff. | Show each flag on the field it concerns. Silent contamination is worse than a gap (rule 8). Empty list with no `as_of_date` is correct; empty list *under* a cutoff for `clinical_stage_small_molecules` is a bug the panel can surface. | P |
| `biologic_precedent.approved_biologics[]` + `note` | Approved antibodies/proteins/fusions. | **Target validation, NOT small-molecule tractability** — often evidence of the opposite. Render in its own block, visually separated from `target_precedent`, carrying the `note` verbatim (rule 1, TNF example). | P (biologic sub-block) |
| `family_precedent.*` | Pfam-family activity, with `best_family_potency_modality`. | The weakest transfer axis; never merged into `target_precedent` (rule 9b). Family best potency may be a peptide — show the modality. | P |
| `structural_neighbour_precedent.*` | Foldseek fold neighbours and whether they carry drug-like holo ligands. | Fold similarity, not family; can disagree with `family_precedent` — show both, merge neither. If the tool was unavailable, this is `null` with a `not_found` reason, never "no neighbours found" (rule 13). | P |
| `pocket_neighbour_precedent.*` | Nearest pockets on other targets + cofold transfer of their ligand into our pocket. | **The strongest transfer axis, but a hypothesis, not a measurement** — label it transferred, name the source target, carry the similarity and cofold result so a reader can discount it (rule 9). `cofold_transfer.*` is nulled in this deployment (rule 13). | P |

### Two-axis separation (N)

The panel renders the P block and the C block as two columns/panes that never
combine. There is no field that spans them; `verdict_basis` and `axis_conflict`
are the only fields that *reference* both, and they name the relationship rather
than blending the numbers.

### Structure (C, provenance for the computed axis)

| Field | Gloss | Caveat that must travel with it | Axis |
|---|---|---|---|
| `structure.tier` | `holo_experimental` / `apo_experimental` / `cofolded` / `predicted` / `sampled_ensemble` / `none`. | The tier gates everything downstream. `predicted` requires the **warnings contract**: render a predicted result with `predicted_plddt` and a warnings banner; a low pLDDT is louder warnings, not a null (rule 4c). In this deployment `cofolded`/`predicted`/`sampled_ensemble` are unreachable (rule 13). | C |
| `structure.predicted_plddt` | Model confidence for a predicted structure. | Only meaningful with `tier = predicted`; shown *inside* the warnings banner, never as a bare quality score. | C |
| `bound_ligand.is_known_frequent_hitter` | Whether the holo ligand is a promiscuous binder. | A holo structure can be exactly the misleading kind — the pocket is real, the evidence it is a productive drug site is not (TNF `307`/SPD304). Show beside `is_druglike`. | C |
| `holo_count` / `apo_count` / `total_pdb_structures` | Structure census. | The population the cryptic census is drawn from; each is a count subject to rule 14 reconciliation. | C |

### Tractability — the computed axis (C)

| Field | Gloss | Caveat that must travel with it | Axis |
|---|---|---|---|
| `pocket_volume_a3` (`min`/`max`/`spread_pct`/`clustering_d`/`primary_d1_6_a3`) | Cavity size, with its spread across the clustering sweep. **The reported computed-axis number.** | Absolute physical quantity — *may* be compared across structures — **but carries no verdict**. The 210/240 Å³ guide is **withdrawn and may not be revived** (rule 4a). Its clustering sensitivity travels with it every time (~492 Å³ swing). `primary_d1_6_a3` is the D=1.6 site volume only, not the pooled min/max; above ~1000 Å³ sites have merged. Render the spread, never a bare midpoint, never a comparison to 210/240. | C |
| `pocket_druggability` (`min`/`max`/`fold_range`/`load_bearing`) | fpocket druggability score range. | **WITHIN-STRUCTURE quantity** — a 3-descriptor logistic on 21 positives, dominant term normalised over *this structure's own* pocket list. Never compared across structures/targets/thresholds; a spread pooled across structures **measures nothing**. `load_bearing` is fixed `false` — it may not carry a `not_tractable`/`insufficient_evidence` verdict. **Do not render as a probability or a bare value**; render its reportable form instead: | C |
| `site_pocket_rank` (`fpocket`/`prank`/`n_pockets`/`structure_pdb_id`) | **The reportable form of druggability.** | Render as **"rank k of n in `<pdb>`"**, never a bare score. `fpocket` and `prank` are two within-structure orderings on equal footing — show both; a disagreement is shown as a disagreement. PRANK is a site-finding aid, inverted (AUC 0.25) as a quality classifier (rule 4d). | C |
| `pocket_volume_a3.site_pocket_selected_by` / `pocket_druggability.site_pocket_selected_by` | How the scored pocket was chosen (single basis or list). | **Only `ligand_site_jaccard` is same-site without qualification.** The other five values (`site_signature_overlap`, `site_signature_unreliable_homooligomer`, `max_druggability_no_ligand_site`, `no_pocket_matched_site_signature`, `no_pocket_overlapped_ligand_site`, `site_signature_unreliable_foreign_polymer`) do **not** identify a site and must not be pooled as one. `max_druggability_no_ligand_site` must never produce a reportable value. The panel must show the basis and mark values on a non-identifying basis as **not-poolable / refused** (TNF example). | C |
| `cryptic_pocket_risk` / `cryptic_mechanism` / `cryptic_evidence.*` | Whether the site is cryptic and by what mechanism, with the apo census behind it. | Cryptic is a geometric measurement, not a flag on apo (rule 5). Mechanism is a **prior on achievable potency** (`loop_or_backbone_motion` → nanomolar-plausible; `sidechain_occlusion`/`subunit_occlusion` → micromolar-at-best). "Cryptic" alone is not actionable. Show `n_apo_examined`/`n_apo_site_absent` (the denominator/numerator) and `site_present_in_apo_ensemble` (occluded-not-cryptic test). A geometric low score under apo-only is an *absence of measurement*, not evidence of poor tractability. | C |
| `cryptic_potency_prior.expected_ceiling` | Best achievable potency implied by the mechanism. | A prior, not a prediction about any specific compound (TNF: balinatunfib is a live counter-example carried explicitly). | C |
| `pocket_vs_interface.classification` | `orthosteric_candidate` / `allosteric_candidate` / `destabiliser_candidate` / `no_partner_structure` / `mixed` / `no_pocket_to_classify` / `numbering_mismatch_not_interpretable`. | Measured against a real partner complex, not assumed. `no_partner_structure` means the classification was not made, not that there is no interface. A pocket claimed orthosteric that does not touch the interface is a mislabelled hypothesis (rule 2b). | C |
| `ensemble_consensus_fraction` (`n_structures`/`n_measurements`/`fraction_with_strong_pocket`/`meets_consensus_criterion`) | Fraction of the ensemble showing a strong pocket. | An **anti-cherry-picking control, not a tractability signal** (AUC 0.560, ranks MYC top — rule 4c). A fraction with no N is not a measurement: `n_measurements` (D-sweep) ≠ `n_structures`. Show the denominator. | C |
| `disorder_fraction` | Predicted disordered fraction. | If measured off a truncated crystal sequence rather than the full accession, it is a scoping error and belongs in `not_found` as a missing measurement — a `0.0` here is not "the protein is ordered" (JAK1 example). | C |
| `site_centroid_to_ligand_distance_a` + `mdpocket_site_definition_used` | Distance from the density-defined site to the ligand-defined site; which definition was used. | Read this before quoting any geometry. `site_from_ligand` is the site by construction; `site_from_density` may be a *different* pocket — a centroid >4 Å (a **proposed, not calibrated** threshold) from the ligand is a distinct cavity (rule 4b). JAK1: 1.86 Å (check passes). TNF: 29.57 Å (check fires). | C |
| `tractability.caveat` | Free-text caveat for the whole computed axis. | Render verbatim, prominently. This is where refusals, chain-selection notes, and poolability reasoning live. | C |
| `tractability.method.*` | Tool, version, `clustering_d_swept`, `ensemble_pdb_ids`, `chains_used`. | Provenance for every computed number. `chains_used` records the rule-2b chain assertion — the site you block depends on it. | C |

### Affinity (C)

| Field | Gloss | Caveat that must travel with it | Axis |
|---|---|---|---|
| `affinity.positive_control_ligand` / `positive_control_measured_nm` / `positive_control_predicted_nm` / `reliable` | The predictor run on a known binder before any novel prediction is trusted. | **A prediction without its control is not a measurement** (rule 12). `reliable: true` licenses triage only — it never orders compounds within a target. In this deployment the whole block is `null` (no predictor, rule 13); render as "no predictor — nulled", with the `not_found` reason. | C |

### Falsification & meta (N)

| Field | Gloss | Caveat that must travel with it | Axis |
|---|---|---|---|
| `falsification.checks_run[]` / `findings[]` / `survived` | What was checked, what was found (including nothing), and whether the claim survived. | A check that found nothing is shown, not hidden — a claim that survived a real attack is worth more than an untested one. `survived: false` is not a verdict flip; it is attached evidence for the reader (TNF survived: false, verdict still tractable). | N |
| `not_found[]` | One entry per axis/field that could not be retrieved, quoting the failure signature verbatim. | **A null is not a zero.** Render each reason verbatim (rule 15's Kind A signatures, tool unavailability, scoping refusals). This block is the panel's proof that "no evidence" is distinguished from "we failed to retrieve it." | N |
| `next_experiment.*` | The single experiment that would most move the dossier. | Context, not a claim about the current verdict. | N |
| `follow_up_questions[]` | OPTIONAL, future capability (see §6). | Not required; a consumer need not act on it. | N |

---

## 4. How the panel renders a click

**The interaction, top to bottom:**

1. **Click a prediction (a dossier).** The panel opens on the **verdict header**:
   `verdict` + `verdict_basis` shown together (never the verdict alone), and
   `axis_conflict` rendered as a prominent banner immediately below when non-null.

2. **Expand the two axes side by side.** Two panes — **Retrieved precedent (P)**
   and **Computed tractability (C)** — laid out so they are visually never merged.
   No combined bar, no averaged score, no shared colour scale. If one axis is
   refused or nulled, its pane says so with the reason; the other pane is
   unaffected.

3. **Drill into any number.** Clicking a value expands its attribution card:
   - **stage** — which pipeline stage produced it (from §2);
   - **raw value** — exactly as it appears in the JSON, to the significant figures
     the dossier used (never re-rounded up; volumes/CVs quoted to 2 sig figs per
     rule 4a);
   - **provenance** — the owning block's `sources[]`, plus the specific
     `pdb_id` / `chembl_target_id` / assay ID / query where the field records one;
   - **caveat** — the required caveat from §3, shown in the same card. If a field
     has a required caveat and the dossier provides none, the panel shows a
     visible "caveat missing" marker rather than a clean number.

**Honest-defaults the renderer must enforce:**

- **`not_found` reasons are shown verbatim.** No paraphrase, no truncation, no
  "data unavailable" substitution.
- **A null is never a zero, an empty list, or a blank.** A field absent because it
  was not retrieved renders as "null — see not_found: `<reason>`", visually
  distinct from a genuine measured `0`.
- **A predicted structure shows its warnings banner** with `predicted_plddt`
  inside it. A low-confidence fold is a flagged result, never a hidden or nulled
  one.
- **A within-structure rank renders as "rank k of n in `<pdb>`"**, never a bare
  druggability score. If only the score is present, show it explicitly labelled
  "within-structure only, not comparable" with the provenance caveat.
- **`site_pocket_selected_by` on a non-identifying basis flags the value as
  not-poolable / refused** and blocks any range display built from it.
- **`axis_conflict` is flagged prominently** — top of the panel, not buried in the
  computed pane.
- **The biologic block is visually fenced off** from `target_precedent` with its
  "validation, not tractability" note.

---

## 5. Two worked examples

### JAK1 (P23458) — axes agree, computed axis informative

- **Verdict header:** `small_molecule_tractable`, `verdict_basis:
  retrieved_precedent`, `axis_conflict: null` — no conflict banner.
- **Retrieved pane (P):** 14,472 distinct actives; top assay only 4.8% of
  bioactivity (single-assay-dominance check *clean*); best potency 0.01 nM from a
  *characterised* enzyme assay; **9 approved small molecules** (ruxolitinib →
  deuruxolitinib). The count/list drill-down surfaces the deliberate one-drug gap:
  count is 9, only 8 named, and the `not_found` entry ("the 9th could not be named
  … and is NOT guessed") renders verbatim. This axis carries the verdict.
- **Computed pane (C):** ensemble `3EYG` + `10PI`, swept D {1.6, 2.4}, **all 4 of
  4 site pockets selected by `ligand_site_jaccard`** → the panel shows the volume
  range 305.9–913.8 Å³ (spread 66.5%) as poolable, with the "no 210/240 guide"
  caveat and the note that 913.8 Å³ approaches the ~1000 Å³ merge threshold.
  Druggability 0.020–0.437 rendered as within-structure only, non-load-bearing.
  `site_centroid_to_ligand_distance_a: 1.86` → the mdpocket off-site check
  **passes** (both site definitions agree).
- **Caveats that fire:** no `mechanism_hypothesis` supplied (no pocket asserted as
  *the* site); cryptic stage `not_run` (all-holo ensemble, nothing to superpose) →
  `max_backbone_ca_displacement_a` is null with a verbatim reason;
  `disorder_fraction` null because measured on the 280-residue kinase domain, not
  the 1154-residue protein — shown as a missing measurement, not "ordered".
- **What the panel demonstrates:** both axes point the same way; the computed axis
  is informative *and* would not have changed the verdict — the drill-down makes
  that explicit.

### TNF-alpha (P01375) — axes disagree; computed axis refused

- **Verdict header:** `small_molecule_tractable`, `verdict_basis: both`,
  **`axis_conflict` populated** → prominent banner: zero approved small molecules
  against five approved biologics, an oral small molecule (balinatunfib /
  SAR441566) in phase 2, and a computed axis that is **refused**. The banner tells
  the reader not to average, and states the mechanistic reconciliation (the site
  is a cavity on the trimer 3-fold axis, reached by displacing a subunit).
- **Retrieved pane (P):** `distinct_actives: 2582` shown with its
  single-assay-dominance finding — **45.0% of bioactivity is an IRAK4 monocyte
  assay measuring a different protein**; the panel surfaces
  `measures_a_different_target: true`. Headline 0.03 nM Ki is rejected as
  uncharacterised; reported potency is the characterised SPR Kd of 1.3 nM. Biologic
  block (5 approvals) fenced off with "validation, not tractability."
- **Computed pane (C): refused, not weak.** `pocket_volume_a3` and
  `pocket_druggability` are `null` in the body; `site_pocket_selected_by` is a list
  `[ligand_site_jaccard, site_signature_unreliable_homooligomer]`, and the panel
  marks the 10-of-12 apo measurements as **not-poolable** (homo-oligomer signature
  cannot tell one protomer's site from another's on a C3 trimer). The
  `tractability.caveat` renders verbatim: pooling anyway regenerates the withdrawn
  **651.0-fold** druggability range. `site_centroid_to_ligand_distance_a: 29.57`
  → the off-site check **fires hard** (~4× the previously recorded error), so no
  `site_from_density` number appears anywhere. The refusal is displayed as a
  refusal, explicitly *not* a low score sitting next to "zero approved small
  molecules".
- **Caveats that fire:** frequent-hitter holo ligand (`307`/SPD304,
  `is_known_frequent_hitter: true`); cryptic mechanism `subunit_occlusion` with
  `is_cryptic: false` and the census `n_apo_examined: 5`, `n_apo_site_absent: 1`
  (occluded, not collapsed); `cryptic_potency_prior.expected_ceiling:
  micromolar_at_best` shown with balinatunfib as a live counter-example;
  `falsification.survived: false` shown as attached evidence, with the verdict
  still `small_molecule_tractable`.
- **What the panel demonstrates:** the two axes disagree and **neither is wrong**;
  the panel's job is to show the disagreement and the refusal honestly, not to
  collapse them into a number.

---

## 6. Data contract

- **The panel consumes exactly `schema/output.schema.json`.** All 17 required
  top-level keys are always present (null when unretrieved); every enum field is
  present and non-null (each enum carries its own unknown member —
  `none`/`undetermined`/`unknown`). Doc keys beginning with `_` may appear
  anywhere and the panel may render them as inline annotations.
- **Degrade gracefully on null/absent optional fields.** Optional blocks
  (`family_precedent`, `structural_neighbour_precedent`, `pocket_neighbour_precedent`,
  `affinity`, most of `tractability`) may be null or empty; the panel renders each
  as "not run / null" with the corresponding `not_found` reason where one exists,
  and never as a substantive zero.
- **Never fabricate a number the dossier left null.** No interpolation, no
  back-filling, no "estimated" values, no averaging of the two axes. If a value is
  null, the panel shows null + reason.
- **Honour the deployment nulls (rule 13).** `affinity.*`, `structure.cofold_control`,
  `pocket_neighbour_precedent.*.cofold_transfer`, the `cofolded`/`predicted`/
  `sampled_ensemble` tiers, the rule-10b modality cross-check, and
  `target_precedent.patents` are all nulled in the current deployment. The panel
  renders these as "no tool in this deployment", with the `not_found` reason, and
  never recalls a value.
- **`follow_up_questions` (OPTIONAL, future affordance).** The schema defines an
  optional array of structured asks (`expand_node` / `resolve_link` / `test_gap`,
  each with a target node/link id and a human-readable question) that the
  graph-intake stage can already emit back to the upstream knowledge graph. It is
  deliberately not required and a consumer is not obliged to act on it. When
  present, the panel may render it as a **"the pipeline could ask upstream X"**
  affordance — a place to surface where more evidence would be sought — clearly
  marked as a future channel, not a current claim. When absent, the panel shows
  nothing there.

---

## 7. Implementation notes

- **Static renderer over the dossier JSON — no server.** Build it the way
  `pipeline.html`, `graph-viewer.html`, and `hypothesis-viewer.html` already work
  in this directory: a single self-contained HTML/CSS/JS file that loads one
  dossier JSON and renders it client-side. Match the existing house style
  (monospace `--mono`, the muted `--bg`/`--ink` palette, `min-width` layout) so the
  panel reads as part of the same station.
- **Every rendered claim traces to a dossier field.** The renderer must not
  compute derived quantities the dossier did not emit — no ratios, no averages, no
  re-normalised scores. If a display value is not a verbatim field (or a verbatim
  field plus its `sources`/caveat lookup), it does not belong on the panel.
- **Encode the caveat contract in code.** Maintain a field→caveat map (from §3) so
  that rendering a caveat-bearing field without its caveat is impossible by
  construction; a missing caveat renders as a visible marker, not a clean number.
- **Encode the poolability rule.** A `site_pocket_selected_by` value outside
  `ligand_site_jaccard` disables any range/spread widget built from that field and
  shows the per-structure values with a "not one site" note instead.
- **Distinguish the four states** every field can be in, with distinct styling:
  measured value · `null` with `not_found` reason · refused (present-but-not-
  reportable, e.g. TNF geometry) · not-run/no-tool (rule 13). Collapsing any of
  these into a blank or a zero defeats the panel's purpose.
- **Validate the input.** Reuse the shape guarantees of
  `.claude/skills/assemble-dossier/validate_dossier.py` as the contract; the panel
  should assume a validator-passing dossier and fail loudly (not silently render a
  zero) if a required key is missing.
- **Test against the two shipped dossiers.** `jak1_P23458.json` (axes agree,
  computed axis filled) and `tnf_P01375.json` (axes disagree, computed axis
  refused, off-site check fires) exercise the agree/disagree, poolable/refused, and
  pass/fire branches; both should render correctly before the panel is considered
  done.
