# As-of date filtering — verified queries

How to restrict each evidence type to what existed before `as_of_date`. Every
query below was executed against the live Paperclip database; the numbers in
"Tested" blocks are retrieved, not estimated.

Verified 2026-08-15.

---

## 1. Bioactivities — THE JOIN (verified)

`chembl_v.bioactivities_by_accession` has no date column. The route to a year is:

```
bioactivities_by_accession.assay_id
  -> chembl.assays.assay_id (.doc_id, .tid)
  -> chembl.docs.doc_id
  -> chembl.docs.year   (integer, 1974-2025, 100911/101100 populated)
```

`chembl.activities` exposes `doc_id` too, but see "Which doc_id" below.

### Recommended query (raw route, activity-level year)

```sql
SELECT COUNT(*)                          AS n_bioactivities,
       COUNT(a.pchembl_value)            AS n_with_pchembl,
       MAX(a.pchembl_value)              AS best_pchembl,
       COUNT(*) FILTER (WHERE d.year IS NULL) AS n_undatable_excluded
FROM chembl.activities a
JOIN chembl.assays s ON s.assay_id = a.assay_id
LEFT JOIN chembl.docs d ON d.doc_id = a.doc_id
WHERE s.tid IN (SELECT tid FROM chembl_v.target_proteins WHERE accession = :acc)
  AND d.year < :cutoff_year;
```

Drop the `AND d.year < :cutoff_year` line and read `n_undatable_excluded` to
report how much evidence the cutoff could not date.

**This raw route reproduces `bioactivities_by_accession` exactly.** Verified row
counts, view vs raw route: P23458 48836/48836, P01375 6447/6447, P01116
18245/18245.

### Equivalent view-based query

If you must start from the view, join through `assays`, never `activities`:

```sql
SELECT COUNT(*) AS n_bioactivities, MAX(b.pchembl_value) AS best_pchembl
FROM chembl_v.bioactivities_by_accession b
JOIN chembl.assays s ON s.assay_id = b.assay_id
JOIN chembl.docs  d ON d.doc_id   = s.doc_id
WHERE b.accession = :acc AND d.year < :cutoff_year;
```

Produced identical numbers to the raw route on all three test targets.

### DO NOT join the view to `chembl.activities` on (assay_id, molregno)

It fans out ~3x. Measured on P23458: 48,836 true rows -> **144,980** rows after
that join. The view does not expose `activity_id`, so there is no non-fanning
key from the view into `activities`. Use `assays.doc_id` (1:1 on `assay_id`) or
go to the raw route.

### Which doc_id — assay-level vs activity-level

`activities.doc_id` and `assays.doc_id` disagree on **643,154** activity rows
corpus-wide. Of those, **480,620** have the assay-level year *earlier* than the
activity-level year (delta up to 12 years). That is the leaky direction:
filtering on assay year would admit a measurement actually published later.

Prefer `activities.doc_id` (the document the measurement was extracted from).
The raw route above does this. On the three test targets both routes gave
identical results, so the divergence is rare per-target — but it is not zero
corpus-wide, so use the activity-level year by default.

### Undatable rows

Some activities hang off `doc_type = 'DATASET'` docs with `year IS NULL` —
deposited screening sets with no publication date. Identified:

| doc_id | title | affected target | n_activities |
|---|---|---|---|
| 61374 | Compounds: GSK PKIS; Assays: Nanosyn kinase panel | P23458 | 738 |
| 51887 | PubChem BioAssay data set | P01375 | 399 |

They cannot be dated at all. A cutoff query excludes them (NULL fails `<`).
That is the correct conservative behaviour, but **report the excluded count** —
do not let it vanish silently.

### Tested — bioactivities

| Accession | Cutoff | n_bioactivities | best pChEMBL |
|---|---|---|---|
| P23458 (JAK1) | < 2010 | 308 | 8.77 |
| P23458 (JAK1) | < 2013 | 1,034 | 10.05 |
| P23458 (JAK1) | no cutoff | 48,836 (738 undatable) | 11.00 |
| P01375 (TNF) | < 2005 | 239 | 8.05 |
| P01375 (TNF) | < 2026 | 6,048 | 10.52 |
| P01375 (TNF) | no cutoff | 6,447 (399 undatable) | 10.52 |
| **P01116 (KRAS)** | **< 2012** | **2** | **5.97** |
| **P01116 (KRAS)** | **< 2026** | **18,245** | **10.70** |

KRAS is the control: a 2012 cutoff collapses 18,245 activities to 2 and drops
best potency by 4.73 log units (10.70 -> 5.97, i.e. 20 pM -> ~1 uM). The
G12C series is correctly invisible before 2013. Earliest KRAS doc year is 2011.
**If a 2012 KRAS cutoff ever returns sub-nanomolar potency, the filter is
broken.**

---

## 2. Approved drugs — works, with a documented lossy edge

`chembl_v.drugs_by_accession.first_approval` is an integer year.

```sql
SELECT
  COUNT(*) FILTER (WHERE max_phase = 4 AND first_approval <= :cutoff_year) AS approved_by_cutoff,
  COUNT(*) FILTER (WHERE max_phase = 4 AND first_approval >  :cutoff_year) AS approved_after_cutoff_excluded,
  COUNT(*) FILTER (WHERE max_phase = 4 AND first_approval IS NULL)         AS approved_date_unknown,
  COUNT(*) FILTER (WHERE max_phase < 4)                                    AS clinical_stage_undatable
FROM chembl_v.drugs_by_accession
WHERE accession = :acc;
```

### NULL handling — the decision

A plain `WHERE first_approval <= :year` is **not** safe as a lone filter.
Corpus-wide breakdown of `drugs_by_accession`:

| max_phase | n_rows | n with NULL first_approval |
|---|---|---|
| 4.0 (approved) | 10,583 | **632** |
| 3.0 | 2,191 | 2,191 |
| 2.0 | 4,545 | 4,545 |
| 1.0 | 969 | 969 |
| -1.0 | 191 | 191 |

So NULL does **not** mean "not approved" — 632 genuinely approved rows (6.0% of
approved) carry no approval year. `WHERE first_approval <= :year` silently
discards all 632 at every cutoff.

**Correct behaviour: four buckets, none silent.**

1. `max_phase=4 AND first_approval <= cutoff` — count as approved precedent.
2. `max_phase=4 AND first_approval > cutoff` — exclude. This is the leakage
   prevention actually doing its job.
3. `max_phase=4 AND first_approval IS NULL` — approved, date unknown. Do **not**
   include (cannot prove it predates the cutoff) and do **not** drop silently.
   Report as an explicit uncertainty count.
4. `max_phase < 4` — clinical-stage. Not datable at all; see section 3.

Rationale for excluding bucket 3 from the headline count: a retrospective run
must not assert precedent it cannot date. Rationale for surfacing it: an
undated approved drug is evidence the target is drugged, and hiding it
overstates the cleanliness of the cutoff.

### Tested — drugs

| Accession | Cutoff | approved by cutoff | approved after (excluded) | approved date unknown | clinical-stage undatable |
|---|---|---|---|---|---|
| P23458 | 2010 | 0 | 11 | 0 | 12 |
| P23458 | 2013 | 2 | 9 | 0 | 12 |
| P01375 | 2005 | 3 | 2 | 0 | 10 |
| P01375 | 2026 | 5 | 0 | 0 | 10 |
| P01116 | 2012 | 0 | 2 | 0 | 1 |
| P01116 | 2026 | 2 | 0 | 0 | 1 |

JAK1 at 2013 returns exactly ruxolitinib (2011) and tofacitinib (2012); at 2010
it returns zero. Matches the known history.

Optional finer granularity: `chembl.formulations.molregno -> chembl.products.approval_date`
(timestamp) gives a real date rather than a year, but covers US products only.
Not required — `first_approval` is sufficient at year resolution.

---

## 2b. Modality — structure is stable under a cutoff, `molecule_type` is not

Modality classification (SKILL.md steps 2 and 2c) interacts with `as_of_date` in
one direction only, and it favours the structural test.

**A molecule's structure does not change.** Classifying a compound from its
`canonical_smiles` gives the same answer whatever the cutoff, so the structural
modality call carries **no leakage risk of its own**. It inherits the leakage
of whatever set you fed it.

**`chembl.molecule_dictionary.molecule_type` is a current-state annotation with
no history**, exactly like `max_phase`. There is no column recording when a
molecule was typed, and 404,621 molecules are still `Unknown` today — some of
which were typed *later* than the cutoff you are asking about. So a
`molecule_type` read under a cutoff is a present-day fact about a past compound.

Practical consequence:

- **Filter the compound set by date first, then classify.** Do not classify the
  full set and then filter — the modality split must be over the same rows as
  `distinct_actives`, or the two numbers describe different populations.
- The date filter for compounds is the bioactivity route in section 1
  (`activities.doc_id -> docs.year`); `compounds_by_accession` has no date
  column any more than `bioactivities_by_accession` does.
- `target_precedent.compound_modality_split`,
  `target_precedent.modality_unknown_count` and
  `target_precedent.best_potency_modality` inherit `distinct_actives`'s
  `as_of_leakage` entry. They do **not** need their own — the structural call
  adds nothing undatable. Say this explicitly rather than leaving them
  unflagged, because an unflagged field next to a flagged one reads as
  "checked and clean".
- **`molecule_type`-derived drug modality does need its own flag under a
  cutoff**, for the same reason `clinical_stage_small_molecules` does: the
  annotation is current-state.

## 3. Clinical candidates — NOT datable

`max_phase` is a **current-state** field. ChEMBL stores no phase history. A
compound at `max_phase=3` today may have been preclinical at your cutoff, and
nothing in the schema distinguishes the two.

Exhaustive sweep of every date-like column in the `chembl` schema (all ~80 raw
tables, `information_schema.columns` on name or type):

| table | column | usable for as-of? |
|---|---|---|
| docs | year | **yes** — bioactivity route |
| molecule_dictionary | first_approval | **yes** — approved drugs |
| molecule_dictionary | usan_year | no — see below |
| products | approval_date | partial — US products only |
| product_patents | patent_expire_date, submission_date | Orange Book only |
| drug_warning | warning_year | not a precedent date |
| chembl_release, version | creation_date | ingestion metadata |

There is no phase-transition or first-in-human table.

`usan_year` is not a usable substitute. On P23458's 23 drug rows it is NULL for
6, and it dates name assignment, not target disclosure — TRICETAMIDE carries
`usan_year = 1963` while sitting at `max_phase = 1` against JAK1.

**Any dated run that reports clinical candidates must flag them as
leakage-prone.** Report the count (bucket 4 above) and the flag; do not report
them as as-of precedent.

---

## 4. Structures — enforceable at source

`pdb_v.structures_by_accession.release_date` is date-typed and fully populated.

```sql
SELECT COUNT(*) AS n_structures, MIN(resolution) AS best_resolution
FROM pdb_v.structures_by_accession
WHERE accession = :acc AND release_date < :as_of_date;
```

Use `release_date`, not `deposit_date` — release is when the structure became
public knowledge.

---

## 5. Trials — enforceable at source

`paperclip skill trials` does not exist (`Error: Failed to fetch skill "trials"`),
and `/trials/{NCT}/meta.json` carries no trial dates — only ingestion metadata.
But the `-s trials` SQL source exposes a full AACT-style `ctgov` schema.

`ctgov.studies` date columns include `start_date`, `study_first_submitted_date`,
`study_first_posted_date`, `primary_completion_date`, `completion_date`,
`results_first_posted_date` — all `date`-typed.

Coverage (`ctgov.studies`, measured):

| metric | value |
|---|---|
| n_studies | 595,652 |
| start_date populated | 590,300 |
| study_first_posted_date populated | **595,652 (100%)** |
| results_first_posted_date populated | 75,260 |
| study_first_posted_date range | 1999-09-20 to 2026-07-24 |

```sql
SELECT COUNT(*) FROM ctgov.studies
WHERE study_first_posted_date < :as_of_date;
```

Use `study_first_posted_date` — the date the trial became publicly visible.
`start_date` may be prospective/estimated and is 5,352 rows short of complete.

Other registries also carry registration dates: `ct_global.ictrp_studies.date_registration`,
`ct_cn.chictr_studies.date_registration`, `ct_jpn.jrct_studies.first_published_date`,
`ct_eu.ctis_studies.start_date`.

**Caveat:** filtering by posted date correctly dates *trial existence*, but
`overall_status`, `phase` and results columns on the same row are current-state.
A trial registered in 2010 will show its 2026 status. Date the existence; do not
report the outcome as as-of.

---

## 6. Patents — unavailable in this deployment

`paperclip skill patents` documents `--since` / `--before` OpenSearch filters and
a `documents.publication_date` column, so the schema supports date filtering in
principle. It is not reachable here:

```
$ paperclip search -s patents "KRAS G12C inhibitor" --since 2020-01-01
Patents sources are not available.

$ paperclip sql -s patents "SELECT publication_number, publication_date FROM documents WHERE ..."
Patents sources are not available.
```

Patent evidence cannot be retrieved at all, dated or otherwise. This is an
availability gap, not a date-filtering gap.

---

## 7. Literature — enforceable at source (year granularity)

`-s pmc` `documents` exposes `pub_year` and `pub_date`. Verified populated on
sampled rows (e.g. PMC9850818 -> 2023). Note the CLI prints `pub_date` under a
`pub_year` header, so treat resolution as **year**, not day.

```sql
SELECT ... FROM documents WHERE source = 'pmc' AND pub_year < :cutoff_year;
```

Aggregate queries over the full `documents` table time out at 15 s. Always
filter to a candidate set first.

---

## Summary — enforceability by evidence type

| Evidence type | Enforceable? | Mechanism |
|---|---|---|
| Structures | (a) at source | `structures_by_accession.release_date` |
| Approved drugs | (a) at source, lossy | `first_approval`; 632 approved rows corpus-wide have NULL |
| Bioactivities | (b) verified join | `assays.assay_id -> docs.year`; excludes DATASET docs with NULL year |
| Trials | (a) at source | `ctgov.studies.study_first_posted_date`, 100% populated |
| Literature | (a) at source | `pmc documents.pub_year` (year granularity) |
| **Clinical candidates** | **(c) NOT enforceable** | `max_phase` is current-state; no phase history anywhere in `chembl` |
| **Patents** | **(c) not retrievable** | source unavailable in this deployment |

A dated run must set a leakage-prone flag whenever it reports clinical
candidates, and must report the undatable-bioactivity and undated-approved
counts alongside the filtered numbers.
