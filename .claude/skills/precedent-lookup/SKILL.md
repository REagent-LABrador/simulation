---
name: precedent-lookup
description: >
  Retrieves what has actually been made against a protein target — measured
  bioactivities, approved and clinical drugs split by modality, structures, and
  family activity — from the Paperclip protein database, joined on UniProt
  accession. Fills the retrieved-precedent axis of the dossier. It does NOT
  compute tractability, does NOT score druggability, and does NOT merge target
  activity with family activity.
---

# precedent-lookup

Everything on this axis comes from one place: Paperclip's protein database,
three schemas joined on UniProt accession. One tool, one join key, dated fields
throughout — which is what makes the `as_of_date` rule enforceable.

## Setup

`PAPERCLIP_API_KEY` in the environment. Then, **before writing any query**:

```bash
paperclip skill proteins
```

This is mandatory, not advisory. The schema is not guessable and wrong column
names are the most common failure. Run it once per task and read it.

**Modality classification needs RDKit.** `modality.py` in this skill directory
does the compound-level call (step 2c). It is 2D only — no conformer embedding,
no force-field optimisation, nothing from geometry.

```bash
micromamba run -n <env-with-rdkit> python modality.py --selftest
```

The selftest is 32 hard cases with ground truth from known chemistry rather than
from ChEMBL, and it must print `PASS` with **0 false small-molecule calls**
before you trust a modality split. Never `conda`/`mamba`; use
`/Users/bb/.local/bin/micromamba`.

## The three schemas

| Schema | Scale | What you take from it |
| --- | --- | --- |
| `uniprot_v` | 574K proteins | identity, sequence, `features` (binding sites with positions), Pfam via `cross_references` |
| `pdb_v` | 177K structures | `structures_by_accession` (with `release_date`), `entry_ligands` (**holo detection**), `entries` (resolution, deposit date) |
| `chembl_v` | ~24M bioactivities | `bioactivities_by_accession`, `compounds_by_accession`, `drugs_by_accession` (with `first_approval`) |

A fourth is reachable and you need it: the **raw** `chembl.*` tables, same
`sql -s proteins` connection. `chembl.molecule_dictionary` is the one that
matters here — it carries `molecule_type` and `structure_type`, the fields no
`chembl_v` view exposes, and modality classification depends on them entirely.

## Procedure

### 1. Identity

```sql
SELECT accession, gene_name, protein_name, organism, sequence_length
FROM uniprot_v.proteins WHERE accession = '<ACC>'
```

### 2. Drugs with modality — `molecule_type`, and ONLY for drugs

`chembl.molecule_dictionary.molecule_type` is the modality field **for approved
and clinical drugs**. It is *not* the field for bioactivity compounds — see
step 2c, which is not optional, and the failure mode "molecule_type abstains on
59.2% of compounds and does it worst where it matters".

`sql -s proteins` reaches the raw `chembl.*` tables as well as the `chembl_v.*`
views, so this is one query, not two:

```sql
SELECT DISTINCT d.drug_name, d.max_phase, d.first_approval,
       md.molecule_type, md.structure_type
FROM chembl_v.drugs_by_accession d
JOIN chembl.molecule_dictionary md ON md.molregno = d.molregno
WHERE d.accession = '<ACC>'
ORDER BY d.max_phase DESC, d.first_approval NULLS LAST
```

`max_phase` 4.0 = approved. On drugs `molecule_type` is the discriminator; the
tests it superseded are in failure modes and must not be reinstated.

The full enum, with counts over the whole `molecule_dictionary` (12 values):

| `molecule_type` | rows | modality | where it goes |
| --- | --- | --- | --- |
| `Small molecule` | 1,920,259 | `small_molecule` **only if corroborated** — see below | `target_precedent` |
| `Protein` | 22,799 | `protein_or_antibody`, often really a **peptide** — read the structure | `biologic_precedent` |
| `Antibody` | 1,032 | `antibody` | `biologic_precedent` |
| `Oligonucleotide` / `Gene` / `Enzyme` / `Antibody drug conjugate` / `Vaccine component` / `Cell` / `Oligosaccharide` | 260 / 191 / 129 / 109 / 90 / 85 / 81 | `other` | `biologic_precedent` |
| **`Unknown`** | 404,621 | resolve from structure (step 2c); `modality_unknown` only if that fails | reported, never dropped |
| **NULL** | 571,492 | resolve from structure (step 2c); `modality_unknown` only if that fails | reported, never dropped |

The last two rows are 976,113 molecules — **40.4% of the table**, and the reason
step 2c exists. Note that **570,221 of the 571,492 NULL rows carry
`structure_type = MOL`**: a NULL `molecule_type` is a missing annotation, not a
missing molecule, and the structure is almost always there to read.

Verified on the three calibration accessions:

| accession | result |
| --- | --- |
| P23458 (JAK1) | **11 of 11** approved rows `Small molecule` / `MOL`; 23 rows total, 21 `Small molecule`, 2 `Unknown`/`NONE` (INCB-047986, GLPG-0555) |
| P01375 (TNF-alpha) | 5 approved: 4 `Antibody`/`SEQ` (infliximab, adalimumab, certolizumab pegol, golimumab) + etanercept `Protein`/`SEQ`. 15 rows total, 2 `Unknown` (ABBV-3373, AZ9773) |
| Q16552 (IL-17A) | 3 approved, **all `Antibody`/`SEQ`** (secukinumab, ixekizumab, bimekizumab); izokibep `Protein`/`SEQ`; 11 rows total, 2 `Unknown` (M-1095, CJM-112) |

On **drugs** the classes separate cleanly. This is a local field — no external
API call is needed for the common case.

**But `Small molecule` is the one value with measured false positives.** ChEMBL
has **no `Peptide` value in the enum at all**, so every peptide it types lands
in `Small molecule` or `Protein`. Four confirmed:

| drug | `molecule_type` | `structure_type` | what it actually is |
| --- | --- | --- | --- |
| ICOTROKINRA (molregno 3283615) | `Small molecule` | `NONE` | oral **IL-23R peptide**, max_phase 2 |
| VANCOMYCIN (70140) | `Small molecule` | `MOL` | glycopeptide antibiotic, MW 1449 |
| ORITAVANCIN (1076140) | `Small molecule` | `MOL` | glycopeptide antibiotic, MW 1793 |
| DAPTOMYCIN (374037) | `Small molecule` | `MOL` | cyclic lipopeptide, MW 1621 |

Three of the four are caught by reading the SMILES (step 2c). **Icotrokinra is
not** — it has no SMILES at all, and `structure_type = NONE` is precisely why.
So: **`Small molecule` + `structure_type NONE` + no SMILES is an unverifiable
claim and must not be counted.** 5,191 ChEMBL molecules carry that combination.
Only the accession join kept icotrokinra out of the IL-17A dossier.

Conversely the biologic values are trustworthy on their face — `Antibody`,
`Protein`, `Enzyme`, `Cell`, `Oligonucleotide`, `Oligosaccharide`,
`Vaccine component`, `Gene`, `Antibody drug conjugate`. Note that `Protein`
absorbs peptide drugs (CYCLOSPORINE, OCTREOTIDE, LEUPROLIDE, CARFILZOMIB,
ROMIDEPSIN are all `Protein`), so read the structure to get the finer label —
it changes nothing about which block they go in, but "peptide" and "protein" are
different findings.

**`Unknown` is a real value and must not be guessed.** Two TNF-alpha drugs and
two IL-17A drugs return it. For a **drug**, map `Unknown` (and NULL) to
modality-unknown, put it in `not_found` with the drug name, and corroborate
against an independent source (CLAUDE.md rule 10b). For a **compound**, do not
stop there — go to step 2c, which resolves it from structure in 99.6% of cases.

### 2c. Compound modality — read the SMILES, not the field

**This step is mandatory before any `distinct_actives` or `best_potency_nm`
claim.** `molecule_type` abstains on 59.2% of bioactivity compounds and abstains
preferentially on the potent ones. The measurements are in the failure mode
below; the consequence is that an agent applying step 2 to compounds reports
**zero small-molecule precedent for IL-17A**.

Use `precedent-lookup/modality.py` (RDKit, 2D only — no conformers, no force
field):

```python
from modality import classify_compound, summarise
v = classify_compound(canonical_smiles, molecule_type, structure_type)
v.modality  # small_molecule | peptide | macrocyclic_peptide | oligonucleotide
            # | oligosaccharide | protein_or_antibody | modality_unknown
```

The rules, and what each rests on:

| modality | test | evidence |
| --- | --- | --- |
| `oligonucleotide` | ≥2 phosphodiester/thiophosphate linkages **and** ≥2 nucleobases | ChEMBL's own `Oligonucleotide` rows |
| `oligosaccharide` | ≥3 glycosidically linked sugar rings | 22/22 on ChEMBL `Oligosaccharide` |
| `macrocyclic_peptide` | ≥4 alpha N-CA-C(=O)-N linkages, backbone ≥25% of heavy atoms, largest ring ≥12 | 98/98 IL-17C, 112/113 overall |
| `peptide` | same, largest ring <12 | 11/11 |
| `protein_or_antibody` | MW ≥5,000 or ≥40 linkages | ChEMBL biologic types |
| `small_molecule` | none of the above, MW ≤1,500, ≤100 heavy atoms | 107/107 |
| `modality_unknown` | anything else, incl. unparsable/absent SMILES | — |

**The alpha-linkage threshold is not tuned.** The count is bimodal with an empty
gap: IL-17A's 117 compounds give `{0:94, 1:5, 2:5, 3:2, 12:6, 13:1, 14:3, 24:1}`.
Any value from 4 to 11 gives the same partition.

**Glycopeptides get their own gate.** Appended sugar dilutes backbone fraction
below the floor by construction — vancomycin 0.248, oritavancin 0.20 — so
`≥4 linkages AND ≥1 sugar ring` is a peptide whatever the backbone fraction.
Do not instead lower the floor; that is fitting to two points.

**Structure outranks `molecule_type`, and an ambiguous structure is NOT rescued
by it.** If the structure parsed and the classifier still returned
`modality_unknown`, that is evidence *against* `Small molecule`, not an absence
of evidence. Only a compound with no structure at all falls back to the field.
Getting this backwards is what let vancomycin and oritavancin through as small
molecules during testing.

**Report `modality_unknown` — never drop it.** Three figures, always:
the per-modality split, the small-molecule count, and the unknown count with a
line in `not_found`.

**A potency figure without a modality is not attributable.** IL-17A's compound
set is mixed: 106 small molecules and 11 peptides, and all 11 peptides sit in
the bucket `molecule_type` abstains on. Carry the modality of the compound that
achieved `best_potency_nm`, and do the same for `family_precedent` — IL-17C's
family best of 1.4 nM is a **macrocyclic peptide**.

### 2d. `drugs_by_accession` UNDER-reports drugs — always cross-check

The view only contains molecules with a `drug_mechanism` row. Real drugs that
are in ChEMBL are missing from it. Confirmed on four targets:

| target | missing from the view | in ChEMBL as |
| --- | --- | --- |
| IRAK4 Q9NWZ3 | KT-474 | `CHEMBL5569030`, `Small molecule`, max_phase 1 |
| IRAK4 Q9NWZ3 | emavusertib | phase 2 |
| TL1A | all three Phase 3 antibodies | — |
| IL-11 P20809 | **bazedoxifene** (approved) | — but it binds **gp130**, not IL-11 |
| MYC P01106 | 7 compounds at phase ≥2 | — |

**It is not uniform** — on CD20, BAFF and IL-13 every approved antibody does
appear. So the conclusion is "always cross-check", not "the view is broken".

The second path, which needs no `drug_mechanism` row:

```sql
SELECT c.compound_chembl_id, c.compound_name, c.max_phase, c.best_pchembl_value,
       (SELECT molecule_type FROM chembl.molecule_dictionary md
         WHERE md.molregno = c.molregno) AS molecule_type
FROM chembl_v.compounds_by_accession c
WHERE c.accession = '<ACC>' AND c.max_phase >= 1
ORDER BY c.max_phase DESC
```

**Neither number is a drug count, and you must say which you are quoting.** The
view is high-precision/low-recall; this query is high-recall/low-precision
because it returns any clinical-stage molecule with *measured activity* on the
target, including off-target screening. Measured:

| target | `drugs_by_accession` | `compounds` max_phase ≥1 |
| --- | --- | --- |
| IRAK4 Q9NWZ3 | 1 (zimlovisertib) | 271 — incl. ruxolitinib, alectinib, axitinib as off-target hits |
| EGFR P00533 | 80 | 1,135 |
| JAK1 P23458 | 23 | 304 |
| TYK2 P29597 | 14 | 295 |
| RORgt P51449 | 1 | 137 |
| KRAS P01116 | 3 | 13 |
| BCL-2 P10415 | 4 | 27 |
| TNF P01375 | 15 | 28 |
| IL-1B P01584 | 5 | 5 |
| MYC P01106 | 0 | 9 |
| IL-17A Q16552 | 11 | **0** |
| IL-17C Q9P0M4 | 0 | 0 |

IL-17A is the case that shows the second path is not a cure: it returns **zero**,
while four oral IL-17A small molecules are in the clinic. A compound with no
ChEMBL activity row against the accession is invisible to both paths. Say so
rather than reporting an absence.

**A degrader's DC50 is not a binding affinity.** KT-474's potent numbers are all
DC50, 0.46–2.0 nM, CRBN-mediated degradation; its only binding number is an
**IC50 of 41 nM** against full-length protein (pchembl 7.39). Degrader precedent
counts as small-molecule precedent; it does **not** establish that the site is
bindable. Quoting the DC50 as affinity overstates the chemistry by more than a
log.

**`structure_type` is a hint, not the test.** It is `MOL` for small molecules and
`SEQ` for sequence-based entities, but it goes `NONE` for entries with no
structure of either kind — and `NONE` appears on genuine antibodies
(VUNAKIZUMAB, REMTOLUMAB on Q16552 are `Antibody`/`NONE`) as well as on
`Unknown` rows. Read `molecule_type`; use `structure_type` only to describe why
a record is thin. AZ9773 is `Unknown`/`SEQ` — sequence-based, so almost
certainly a biologic, but ChEMBL declines to say and so do you.

### 2e. Collapse salt and parent forms before counting

Salt, hydrate and parent forms are **distinct molregnos**, so deduplicating on
`molregno` does not deduplicate drugs. Verified on JAK1 (P23458): the 11 approved
rows carry 11 distinct molregnos but represent **9 real drugs**. Two parents
appear alongside their own salts —

- FILGOTINIB `1763569` and FILGOTINIB MALEATE `2336138`
- MOMELOTINIB `617563` and MOMELOTINIB DIHYDROCHLORIDE MONOHYDRATE `3283827`

— while four others appear *only* in salt form (ruxolitinib phosphate,
tofacitinib citrate, upadacitinib hemihydrate, deuruxolitinib phosphate), so you
cannot simply drop rows whose name contains a counter-ion.

Collapse by stripping trailing salt/hydrate tokens from `drug_name` —
`PHOSPHATE`, `CITRATE`, `MALEATE`, `MESYLATE`, `SUCCINATE`, `HEMIHYDRATE`,
`MONOHYDRATE`, `DIHYDROCHLORIDE`, `HYDROCHLORIDE`, `TOSYLATE`, `FUMARATE`,
`SODIUM`, `POTASSIUM` — and group on the remaining stem, keeping the earliest
`first_approval` for the group. If you do not collapse, report the raw count as
what it is: **11 approved rows for 9 approved drugs on JAK1, an inflation of 2**.
Never present a row count as a drug count without saying which you did.

### 3. Compound-level potency

```sql
SELECT COUNT(*) AS n_compounds, MAX(best_pchembl_value) AS best_pchembl
FROM chembl_v.compounds_by_accession WHERE accession = '<ACC>'
```

`pchembl` is −log10(molar): 9.0 = 1 nM, 6.0 = 1 µM. Convert to nM for the
dossier. Pull the top compounds with SMILES to confirm they are small molecules.

### 4. Assay provenance — run this BEFORE reporting any actives count

```sql
SELECT LEFT(assay_description, 55) AS assay, COUNT(*) AS n,
       ROUND(100.0 * COUNT(*) / SUM(COUNT(*)) OVER (), 1) AS pct
FROM chembl_v.bioactivities_by_accession
WHERE accession = '<ACC>'
GROUP BY assay_description ORDER BY n DESC LIMIT 5
```

If one assay exceeds ~30% of all activity, the count describes that assay, not
the target. Put the assay name and share in
`target_precedent.assay_concentration` and read the description carefully — it
may be measuring a different protein entirely.

Also take the type split, but see the failure mode below before trusting it:

```sql
SELECT assay_type, COUNT(*) AS n FROM chembl_v.bioactivities_by_accession
WHERE accession = '<ACC>' GROUP BY assay_type ORDER BY n DESC
```

### 5. Assay-level detail when potency looks surprising

```sql
SELECT standard_type, standard_relation, standard_value, standard_units,
       pchembl_value, confidence_score, assay_description
FROM chembl_v.bioactivities_by_accession
WHERE accession = '<ACC>' AND pchembl_value IS NOT NULL
ORDER BY pchembl_value DESC LIMIT 20
```

`confidence_score` is ChEMBL's target-assignment confidence, 9 high to 0 low.
`standard_relation` of `>` is a **non-result** — see failure modes.

### 5. Structures, with holo detection

```sql
SELECT s.entry_id, s.resolution, s.exptl_method, s.release_date,
       l.comp_id, l.ligand_name, l.formula_weight, l.drugbank_id
FROM pdb_v.structures_by_accession s
LEFT JOIN pdb_v.entry_ligands l ON l.entry_id = s.entry_id
WHERE s.accession = '<ACC>'
ORDER BY s.resolution NULLS LAST
```

`entry_ligands` carries `comp_id`, `smiles`, `inchikey`, `formula_weight` and
`drugbank_id` per entry — this is how you find a holo structure without leaving
Paperclip. `release_date` is what the `as_of_date` filter keys on.

### 6. Family precedent — separately, never merged

Get the Pfam from `uniprot_v.cross_references` (`database = 'Pfam'`), find sibling
accessions, and query their activity. Report as its own object with its own
sources. Never fold it into target precedent, never apply a discount.

## Failure modes

### `molecule_type` abstains on 59.2% of COMPOUNDS, and worst where it matters

This is the largest failure mode in this skill. The field is sound on drugs
(11/11 on JAK1's approved rows) and does not transfer to bioactivity compounds.

Measured on all twelve fixture targets, one query per accession, each bracketed
by a canary and each `n` reconciled against an independent `COUNT(*)`. Three
states are distinguished, because two is what causes this bug: `molecule_type`
populated, `molecule_type` NULL/`Unknown`, and query-did-not-complete. **`no
molecule_dictionary row` was 0 for all twelve** — every NULL is a real NULL, not
a missing join partner.

| target | acc | compounds | typed | **typed %** | `Unknown` | NULL | best pchembl typed | best pchembl abstained |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| EGFR | P00533 | 19,490 | 12,783 | **65.6%** | 1,554 | 5,153 | 11.00 | 11.00 |
| BCL-2 | P10415 | 4,350 | 2,385 | **54.8%** | 376 | 1,589 | 11.00 | 11.00 |
| TNF-alpha | P01375 | 2,582 | 1,077 | **41.7%** | 92 | 1,413 | 10.52 | 9.28 |
| JAK1 | P23458 | 14,472 | 5,819 | **40.2%** | 988 | 7,665 | 11.00 | 11.00 |
| TYK2 | P29597 | 10,603 | 2,970 | **28.0%** | 706 | 6,927 | 10.70 | **11.00** |
| RORgt | P51449 | 12,900 | 2,892 | **22.4%** | 1,194 | 8,814 | 10.77 | 10.72 |
| IL-17A | Q16552 | 117 | 26 | **22.2%** | 34 | 57 | 6.26 | **9.10** |
| MYC | P01106 | 1,249 | 183 | **14.7%** | 979 | 87 | 9.22 | **9.70** |
| IL-11 | P20809 | 15 | 2 | **13.3%** | 0 | 13 | — (no pchembl) | 6.85 |
| KRAS | P01116 | 3,597 | 163 | **4.5%** | 448 | 2,986 | 9.40 | **10.70** |
| IL-1B | P01584 | 351 | 16 | **4.6%** | 1 | 334 | 9.34 | **9.40** |
| IL-17C | Q9P0M4 | 98 | 0 | **0.0%** | 0 | 98 | — | 8.85 |
| **total** | | **69,824** | **28,466** | **40.8%** | 6,372 | 34,986 | | |

**It abstains on the majority of compounds on every target without exception**,
and the abstained set carries the better best-potency on **6 of 12** and ties on
4. It is never usefully enriched for potency. On IL-17A it is actively
anti-correlated: the 26 typed compounds top out at 6.26 — the RORgt secretion
assay ceiling, a contaminated readout — while the 91 abstained reach 9.10.

The abstained compounds almost all have structure: **99.6% carry SMILES**
(IL-17A 91/91, IL-17C 98/98, JAK1 8,611/8,653, RORgt 10,000/10,008, KRAS
3,431/3,434, MYC 1,066/1,066). So step 2c resolves them.

**And the field is wrong in both directions, not merely sparse.** Confirmed
false positives: ICOTROKINRA, VANCOMYCIN, ORITAVANCIN, DAPTOMYCIN (step 2), plus
four linear TSLP peptides from a 2017 paper returning `Small molecule`/`MOL`.
Related targets where SMILES inspection does all the work: IL-11 — 13 of 15
compounds are cyclic peptides and its recorded 140 nM best belongs to a
**1,938 Da macrocycle**; TSLP — the five most potent are bicyclic peptides,
all NULL.

**Do not use `molecule_type` to classify a compound. Use
`precedent-lookup/modality.py`.**

*Provenance of the table above, stated so a reader can discount it.* Each row is
a **single-row server-side aggregate**, which the moving row cap cannot truncate
— that is the structural argument, and it is the reason the sweep was written
as aggregates rather than as row pulls. Supporting checks that were actually
run: the four counts partition `n` exactly on all twelve rows; the whole sweep
was executed twice about ten minutes apart with identical results; canaries on
`compounds_by_accession`, `molecule_dictionary` and `drugs_by_accession` passed
immediately before and immediately after the batch; and four figures reconcile
against numbers recorded independently of this run — IL-17A 117, IL-17C 98,
RORgt 12,900, JAK1 23 drug rows. **The per-target reconciliation against a
separately issued `COUNT(*)` did not complete** — the SQL backend went down
mid-task and stayed down — so treat the eight unreconciled rows as provisional.
**One figure does not reconcile: MYC 1,249 compounds here against 1,079 recorded
in `CLAUDE.md` rule 6.** Both may be right (different predicates) and neither has
been re-derived; do not quote either as settled until it is.

### The `modality` column is empty — do not use it

`chembl_v.bioactivities_by_accession` has a `modality` column. It is NULL.
Verified on IL-17A: all 305 rows return `modality = NULL`, one group. It exists
in the schema and carries nothing.

### `action_type` does not distinguish antibodies from small molecules

This is the trap this skill exists to prevent. IL-17A (Q16552) returns eleven
drugs, three approved:

| drug | max_phase | first_approval | action_type | molecule_type |
| --- | --- | --- | --- | --- |
| SECUKINUMAB | 4.0 | 2015 | INHIBITOR | **Antibody** |
| IXEKIZUMAB | 4.0 | 2016 | INHIBITOR | **Antibody** |
| BIMEKIZUMAB | 4.0 | 2021 | INHIBITOR | **Antibody** |

All three are monoclonal antibodies. All three say `INHIBITOR`. Nothing in
`drugs_by_accession` itself marks them as biologics. An agent that reports "three
approved inhibitors" for IL-17A has produced the wrong answer to the only
question this dossier asks.

**The discriminator is `chembl.molecule_dictionary.molecule_type`** (step 2). On
Q16552 it returns `Antibody` for all three, and `Small molecule` for none of the
eleven. Apply that test to every drug, always, before it enters
`target_precedent`.

A name ending in `-mab` is a useful cross-check but not the test — `IZOKIBEP`
(`Protein`) and `M-1095` (`Unknown`) are also biologics and neither ends in
`-mab`.

### SUPERSEDED — the NULL-SMILES test and its cross-accession confirmation

**Both of these are void. Do not reinstate either.** They are recorded here
because they were the documented procedure, they look plausible, and a reader who
does not know they were tried will invent them again.

*The earlier procedure was:* (1) treat a NULL `canonical_smiles` in
`chembl_v.compounds_by_accession` as a candidate biologic; (2) confirm it with a
cross-accession query asking whether that molregno has SMILES under **any**
accession, on the theory that salt forms of small molecules carry no bioactivity
against the target in hand and would otherwise read as biologics.

Step 1 is a real signal but it over-fires, exactly as previously documented: on
EGFR, nine drugs returned no SMILES and only four were real biologics
(cetuximab, panitumumab, necitumumab, amivantamab) — the rest were salt forms
(osimertinib mesylate, neratinib maleate, mobocertinib succinate, lazertinib
mesylate). That part of the finding stands.

**Step 2 does not work, and it was the part that was supposed to fix step 1.**
Verified by execution: the confirmation query returns **0 rows for both
classes**.

| molregnos queried | rows returned |
| --- | --- |
| JAK1 salt forms — upadacitinib hemihydrate `2832770`, filgotinib maleate `2336138`, deuruxolitinib phosphate `2464813`, momelotinib dihydrochloride monohydrate `3283827` | **0** |
| TNF-alpha biologics — etanercept `675371`, adalimumab `675482`, infliximab `675617`, certolizumab pegol `675782`, golimumab `675784` | **0** |

Salt forms are absent from `compounds_by_accession` entirely — no bioactivity
under any accession — and so are antibodies. The output is *identical* for
approved small molecules and approved antibodies, so `has_smiles_anywhere` can
never be `1` for the cases it was written to rescue, and the check cannot
discriminate. It produces false biologic calls on every JAK1 salt form it sees.

The claim that `SALIRASIB` "has SMILES under 33 other accessions" is not a
counter-example: a compound with broad bioactivity is a case the check never
needed to rescue.

Use step 2's `molecule_type` instead. It calls all four JAK1 salt forms
`Small molecule` and all four TNF-alpha antibodies `Antibody`.

**Also superseded:** the previous instruction to read `structure_type: NONE` as
"unresolvable" (APG-2575, JTE-151). `structure_type` is not the modality field —
`NONE` appears on drugs whose `molecule_type` is a confident `Antibody`. Decide
on `molecule_type`, and reserve modality-unknown for `molecule_type = 'Unknown'`.

### `drugs_by_accession` returns one row per mechanism, not per drug

LAZERTINIB appears three times in the EGFR small-molecule bucket and twice more
with NULL SMILES. Counting rows overstates drug counts. Deduplicate on
`molregno` before reporting any total.

**And `molregno` deduplication is not enough** — salt and parent forms are
*distinct* molregnos, so JAK1 still returns 11 approved rows for 9 approved
drugs after deduplicating. Collapse salt/parent pairs as well; step 2e says how.

### Approved biologics and tractable small molecules can coexist

Same target, both true: IL-17A has three approved antibodies **and** 117
compounds in ChEMBL with a best pchembl of 9.10 — 0.79 nM, genuine fluorinated
med-chem, real SMILES. So IL-17A is *not* "no small molecules exist". It is
"the approved drugs are biologics, and potent small molecules exist but none
approved."

Report both facts. Collapsing either direction is wrong: "druggable, three
approved drugs" ignores modality; "not small-molecule druggable" ignores 0.79 nM
compounds. This is why `biologic_precedent` is its own block.

### `standard_relation` of `>` is a failed measurement

`EC50 > 10000 nM` means the compound did **not** work up to 10 µM. Filtering on
`standard_value` alone silently turns non-results into weak actives and inflates
the count. Always read `standard_relation`; only `=` is a measurement.

### An actives count can be dominated by an assay for a different protein

TNF-alpha (P01375) has 6,447 activities. The single largest contributor:

| assay | n | pct |
| --- | --- | --- |
| IRAK4 Monocyte TNFalpha Cell Based Assay: Cryopreserved | 2901 | **45.0** |
| Inhibition Assay: Inhibition assay using TNF-alpha. | 577 | 8.9 |
| TNF-alpha Secretion Assay: Monocytic THP-1 cells | 321 | 5.0 |

**45% of TNF-alpha's bioactivity is an IRAK4 assay** — a different target, using
TNF only as a cellular readout. Report the count without this check and TNF-alpha
looks heavily precedented. It has zero approved small molecules.

### `assay_type = 'B'` does NOT mean the assay is clean

The obvious defence — filter to binding assays — does not work. Verified on
TNF-alpha: the split is **B = 5,830 / F = 617**, so ~90% are labelled binding,
**and the IRAK4 cellular assay is among them**. The type field is too coarse to
separate a direct-binding measurement from a cellular readout.

Report the split, but do not treat `B` as a filter. The assay *description* is
the only reliable signal, which is why step 4 is mandatory rather than optional.

### `confidence_score = 9` does NOT mean the assay is attributed to your protein

The obvious second defence — filter on ChEMBL's target-assignment confidence —
does not catch assay misattribution either.

The IRAK4 monocyte TNF-alpha cell-based assay that makes up 45.0% of TNF-alpha's
bioactivity is `assay_id 2591534` with 2,901 rows. The *same assay* appears on
IRAK4, correctly attributed, as `assay_id 2591535` with 2,895 rows. The two
copies are not identical — four compounds exist only under the TNF-alpha one.
**Both carry `confidence_score = 9`**, the maximum.

So ChEMBL's confidence field does not flag this class of error, and it cannot be
used as a filter for it. Neither can `assay_type`. Only the assay *description*
separates a direct measurement from a cellular readout of a different protein,
which is why step 4 is mandatory.

A second instance, same shape: **IL-6**'s best ChEMBL potency, 1.1 nM at
`confidence_score` 8 with `standard_relation = '='`, comes from a
**glucocorticoid receptor** fluorescence-polarisation assay. A different protein
contributes 7 of IL-6's top 12 values.

### `n_target_components > 1` means the hit is inherited

In `target_proteins`, `n_target_components = 1` is a clean single-protein target.
Greater than 1 means a complex or family, and activity attributed there is not
necessarily activity against your protein. Check it before counting actives.

### `drugs_by_accession` empty does not mean no chemistry

The view only includes drugs with an annotated **direct mechanism of action**. A
target can have thousands of bioactivities and no rows here. Empty means "no
drug with a curated direct mechanism", not "nothing has ever been made". Say
which you mean.

### THE ROW CAP MOVES. Reconcile every count against an independent `COUNT`.

**This is the most dangerous behaviour on this tool and nothing below it in this
file comes close.** The cap is not a fixed 200. Measured 2026-08-15 while
regenerating fixture counts: **the same query returned 200 rows one moment and
exactly 10 the next** — same SQL, same source, well-formed table, correct
columns, **no error, no warning, no truncation marker**. The first run recorded
**KRAS as 10 PDB entries against a true 522**. Only a separately issued `COUNT`
caught it.

A silently 5%-complete result set is the exact failure this skill exists to
prevent: it is a retrieval failure wearing the costume of a finding, and it
reads as "this target has almost no precedent". It is also the leading
explanation for inconsistencies we had filed elsewhere — a table that looked
degraded and then returned in 7 ms, latencies varying by two orders of
magnitude, and two agents reporting different figures off the same query.

**So, operating rule, not advice:**

- **Every count that reaches the dossier is reconciled** against an
  independently issued `SELECT COUNT(*)` / `COUNT(DISTINCT …)` over the same
  predicate, issued as its own call. This covers `distinct_actives`,
  `approved_small_molecules_count`, `total_pdb_structures`, `holo_count`,
  `family_actives`, patent counts, and **any list whose length becomes a
  number** — drugs, structures, trials, terminated programs.
- **A mismatch is a hard failure.** The field is `null` and `not_found` names
  both figures and both queries. Do not report the bigger one. Do not report the
  aggregate with a footnote. See rule 14 in `CLAUDE.md`, including where to
  record the reconciliation (the block's `sources` list, until the proposed
  per-block `count_reconciliation` field exists).
- **Best of all, never count rows.** Aggregate server-side with `COUNT`, `MAX`,
  `STRING_AGG ... GROUP BY`, or use `paperclip export` for large result sets. A
  one-row result cannot be capped. Reconciliation is for when you needed the
  rows themselves.
- **Exactly 200 or exactly 10 rows is capped until proven otherwise** — but a
  result of 47 rows is **not** thereby safe. We do not know what moves the cap
  or what other values it takes. Only the aggregate clears a count.

### Four Paperclip failure signatures, and all four mean the query did not run

**11 of 30 SQL calls in one dry run failed**, across four distinct signatures,
three of them undocumented:

| signature | what it actually is |
| --- | --- |
| `[error] Request timed out` | seen at 120 s on a **tableless `SELECT 1`**. Not a statement-cost signal; it says nothing about your query. |
| `[error] Something went wrong. Please try again.` | undocumented. No code, no detail, transient and permanent indistinguishable. |
| `vsh: cd: /papers/: Permission denied` | returned by `paperclip sql` **for a SQL query** — a shell error from another subsystem, naming a path you never asked for. |
| a silently capped row set | the section above. **No error text at all** — the only one of the four that does not announce itself. |

**Any of these means the query did not run.** The value is `null`, the reason
goes in `not_found` quoting the signature verbatim, and the retry is
short-then-long. Never `0`, never `[]`, never "no approved small molecules",
"no bioactivities" or "no precedent found". `drugs_by_accession` returning
nothing is a finding; `drugs_by_accession` erroring is not.

**Only auth failures are guarded, and that guard catches none of the four.** The
tool layer throws on `401`/`403`/`unauthorized`/`forbidden`/`invalid api key`
and only on a non-zero exit. Timeouts, "Something went wrong", `Permission
denied` and a capped table match none of those patterns, and the last one is not
even a failed run. Nothing upstream of you will stop these.

**Keep an accession predicate on `compounds_by_accession`, always.** A query
against it whose only filter was `molregno IN (SELECT ... FROM
chembl.molecule_dictionary WHERE pref_name IN (...))` was **cancelled by the
statement timeout at 85s**. Rewriting the same query with the molregnos as
literals returned in **13 ms**. Resolve molregnos in a separate query and inline
them; never make the planner scan that view unfiltered.

### A JOIN to `chembl.molecule_dictionary` times out; the same thing as a correlated subquery is instant

The documented advice — inline literals, avoid subqueries — **inverts** on this
join. Measured, same predicate, same session:

| form | result |
| --- | --- |
| `JOIN chembl.molecule_dictionary md ON md.molregno = c.molregno` | **timed out at 120 s** |
| `(SELECT md.molecule_type FROM chembl.molecule_dictionary md WHERE md.molregno = c.molregno)` as a scalar subquery | **40 ms** |

Both over `compounds_by_accession WHERE accession = 'Q16552'`, 117 rows, and
both return identical values. That is a factor of 3,000. The same inversion has
been seen on a Pfam cross-reference join (85.1 s with inline literals, 2.2 s as
a subquery).

**So try both forms when one is slow, and do not assume the documented shape is
the fast one.** The scalar-subquery form is what the twelve-target abstention
sweep runs on; the JOIN form cannot complete on any of the large targets.

Use `COUNT(*) FROM chembl.molecule_dictionary md WHERE md.molregno = c.molregno`
as a second scalar subquery when you need to tell *"the row is missing"* from
*"the field is NULL"*. They are different findings and the join collapses them.

### Retrieving SMILES through `paperclip sql` corrupts them silently, twice over

Two independent truncations, neither of which announces itself, and RDKit will
happily parse the wreckage.

1. **The table renderer truncates long text fields with a literal `...`.** A
   250-char SMILES prints whole; a 326-char one comes back as 264 characters
   ending in `...`. There is no error. Verified against IL-17C, whose SMILES run
   258–326 characters — **93 of 98 came back truncated** on the first attempt.
2. **The renderer right-pads every field to the column width**, so a
   reassembled value picks up trailing spaces.

And then: **`Chem.MolFromSmiles` treats whitespace as the start of the name
field and silently returns the truncated molecule.** No exception, no warning.
A truncated macrocyclic peptide parses as a smaller, still-valid peptide, and a
truncated one can parse as a small molecule.

**So retrieve SMILES in fixed-width chunks and reassemble:**

```sql
WITH pool AS (
  SELECT c.molregno, c.canonical_smiles s FROM chembl_v.compounds_by_accession c
  WHERE c.accession = '<ACC>' ORDER BY c.molregno LIMIT 25 OFFSET <n>
)
SELECT p.molregno::text || '~' || g.i::text || '~' ||
       SUBSTRING(p.s FROM g.i*120+1 FOR 120) AS r
FROM pool p CROSS JOIN generate_series(0,3) g(i)
WHERE LENGTH(p.s) > g.i*120 ORDER BY p.molregno, g.i
```

`rstrip()` every chunk, require the chunk indices to be gapless `0..n`, and
**verify the reassembled length against `MAX(LENGTH(canonical_smiles))` from a
separate query**. That check is what caught this; 16 of 16 named drugs matched
exactly once it was applied. Keep the page small — a result over roughly 100
rows switches to a 5-row preview and you lose the rest with no error.

### The per-protein document is not the database

`paperclip cat /proteins/<ACC>/content.lines` returns a pre-generated 8-line
summary with PDB and DrugBank cross-references. It contains **no ChEMBL data**.
Concluding from it that Paperclip lacks bioactivity is wrong — the data is in
`chembl_v` via `sql -s proteins`. The document is a summary; the SQL views are
the source.

### `scan` re-dumps the whole document per pattern

`paperclip scan <file> "A" "B" "C"` prints the entire document once per pattern.
Four patterns on an 8-line document produced ~200 lines. Use `grep` on the
section you want.

## As-of filtering

When `as_of_date` is set, filter at the source:

- structures — `WHERE s.release_date <= '<DATE>'`
- approvals — `WHERE d.first_approval <= EXTRACT(YEAR FROM DATE '<DATE>')`
- bioactivities — `bioactivities_by_accession` has no date column; join through
  `chembl.activities` / the document year, or mark the count `leakage_risk: true`
  and say so. **Do not silently report a current count under a past date.**

## Output

Fill `target_precedent`, `biologic_precedent`, `family_precedent`, and the
`structure` block. Every number carries its source — ChEMBL target ID, PDB ID, or
the query itself. Anything not retrieved is `null` with a line in `not_found`.
