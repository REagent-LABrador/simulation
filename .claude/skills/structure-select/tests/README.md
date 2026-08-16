# `ligand_filter` accuracy harness

**Every accuracy figure quoted for `ligand_filter.py` anywhere in this repo is
produced by running a file in this directory.** Before 2026-08-15 they were not:
the harnesses and their data lived only in session scratch under `/private/tmp`,
so the figures were unverifiable the moment that session ended. That is the same
"not followable from this checkout" failure that let the retracted volume
separation stand as long as it did. It is closed here.

## Run it

```bash
cd .claude/skills/structure-select/tests
python3 test_v2.py            # the master run — every figure, one command
python3 test_ligand_filter.py # ground-truth set alone, with the historical re-run
python3 test_holdout.py       # blind held-out sample alone
```

Pure stdlib, Python 3.13, **fully offline**. `offline.py` forces
`ligand_filter`'s default source offline so no Paperclip call is made and no
network is touched. Nothing here depends on the row cap, the `cli_cwd` bug, or
any live backend — which is the point.

### One harness in here is NOT offline: `test_homolog_transfer.py`

```bash
python3 test_homolog_transfer.py   # the four measured anchor cases, ~7 s
```

It reproduces the four anchor-agreement cases through `homolog_transfer.py`'s
three guards, plus a self-transfer positive control and the 7KRZ guard-2 case.
It is pure stdlib and still **needs the network**, because it needs
coordinates: nine mmCIFs cannot be parked beside a skill without shipping them
to the Skills API, since `deploy.ts` zips `.claude/skills/<dir>/` whole. So
entries are fetched from RCSB and cached under `$STRUCTURE_SELECT_CACHE`
(default `~/.cache/structure-select/mmcif`, outside the skill directory on
purpose), and chem-comp rows come from RCSB's REST API rather than Paperclip.

**That last choice is not a preference.** On 2026-08-15 `paperclip sql -s
proteins "SELECT comp_id, type FROM pdb_v.chemcomps WHERE comp_id = 'ADP'"`
returned `[error] Request timed out` after 122 s **with exit code 0**, from a
config carrying `{"cli_cwd": "/"}` — two of the four failure signatures at
once. Guard 1 must not be silently skippable, so it was given a source with no
credentials and no shared mutable state. The verdicts still come from
`ligand_filter.classify_record`; only the row does not.

Measured, 2026-08-15: **6 of 6 cases reproduce as expected.** The table is in
`../SKILL.md` under "Transferred-homolog site anchoring".

## Measured, 2026-08-15, from this checkout

| figure | value | produced by |
| --- | ---: | --- |
| original ground-truth set | **259/262 = 98.9%** | `test_ligand_filter.py`, `test_v2.py` block `gt262` |
| + 9 chemistry cases | 9/9 = 100% | `test_v2.py` block `gt_add` |
| + 9 named-entry context cases | 9/9 = 100% | `test_v2.py` block `context` |
| **combined ground truth** | **277/280 = 98.9%** | `test_v2.py` |
| blind held-out sample | **61/70 = 87.1%** | `test_holdout.py`, `test_v2.py` |
| **held-out false positives** | **0/70 = 0.0%** | `test_holdout.py`, `test_v2.py` |

The three standing misses are `BTN` (biotin → `druglike`), `ACE` and `NH2`
(polymer capping groups → additive / ion). `ACE` and `NH2` are **correct once a
`StructureContext` is supplied** — the CCD lists them in `_entity_poly_seq`, so
they are residues, not ligands. The context block covers `NH2@8B9P` for exactly
this reason and passes.

All 9 held-out disagreements run in the conservative direction: nothing that was
really a cofactor, lipid or additive was called drug-like. That asymmetry is the
deliberate bias — a false negative costs a holo structure, a false positive
*invents* one, and inventing one is what produced all four historical bugs.

## Two named boundaries on the zero

The zero is a statement about **that 70-component sample**, not about the
classifier. Two counterexamples exist and are not in it:

- **TNF 5UUI's `MTN`** spin label — a genuine false positive that no chemistry
  rule can fix, because the disqualifying fact is *why the component is there*,
  not what it is made of.
- **`LFI`**, a peptide-conjugated crosslinker — the first false positive actually
  measured against the held-out result. Now closed by `StructureContext`, and
  the fix is regression-tested in the `context` block.

Quote the zero with both boundaries attached.

## What is in here

| file | role |
| --- | --- |
| `test_v2.py` | **master runner.** Ground truth + additions + context + held-out in one pass, with the confusion matrix and the flag assertions. |
| `test_ligand_filter.py` | the 262-component ground-truth set, plus the four historical failures re-run as holo/apo calls and the genuinely-holo controls. |
| `test_holdout.py` | the 70-component blind sample and its by-name adjudication. |
| `gt_additions.py` | the 9 chemistry + 9 context cases added by the polymer-conjugate work, with `FLAG_REQUIRED` / `FLAG_FORBIDDEN`. |
| `offline.py` | forces the chem-comp source offline. |
| `chemcomps.json` | cached `pdb_v.chemcomps` rows for the ground-truth set. |
| `extra_recs.json` | cached rows for the added components. |
| `holdout.json` | cached rows for the blind sample. |
| `entry_counts.json` | `n_pdb_entries` per comp_id, for the ubiquity rule. **Two sources exist — see below.** |
| `structures/*.cif` | the 7 mmCIFs the context cases build a `StructureContext` from: 8QFZ, 8B9P, 3QN7, 9Q8N, 5V2P, 6OIM, 4G5J. **Header-trimmed** — see below. |
| `trim_cif.py` | the trimmer, kept so the trim is reproducible rather than a one-off. |
| `version_diff.py` | compares two `ligand_filter.py` versions verdict-by-verdict over all 346 cached rows. Written because the classifier sha256 pinned in `fixtures/targets.json` turned out to name a file that existed only in session scratch; this is how you check whether a version difference actually moves anything. |
| `test_homolog_transfer.py` | the four anchor-agreement cases through `homolog_transfer.py`'s three guards, plus a self-transfer positive control and the 7KRZ guard-2 case. **Network-dependent** — see above. |
| `chemcomps_transfer.json` | RCSB chem-comp rows cached by that harness (6 components). |

### The mmCIFs are trimmed, deliberately

`scripts/deploy.ts` zips each `.claude/skills/<dir>/` **whole, with no
exclusions**, so anything parked beside a skill ships to the Skills API. Full
entry files for the seven context cases are **3.6 MB**; `StructureContext` reads
only six header categories (`_entity`, `_entity_poly`, `_entity_poly_seq`,
`_struct_asym`, `_struct_conn`, `_struct_ref`) and `_atom_site` is not one of
them. Trimmed to those, they are **88 KB — 41x smaller — and every figure is
identical**, context block included at 9/9.

Regenerate with `python3 trim_cif.py <full.cif> <out.cif>`, inspect with
`python3 trim_cif.py --check <file.cif>`. The retained `_struct_conn` rows are
the covalent-linkage evidence the polymer-conjugate rules turn on, so a verdict
is still auditable from the file — nothing load-bearing was dropped.

**Do not swap in an `-assembly1.cif`.** RCSB strips `_struct_conn` from assembly
files, and a context built from one silently cannot see a covalent link. The
classifier guards this with `has_struct_conn_category`, but the guard turns the
case into "not applied", not into a correct answer.

### `entry_counts.json` disagrees with `../SKILL.md`, and both are right

Three comp_ids differ between this cache and the counts quoted in the skill:

| comp_id | `entry_counts.json` | `pdb_v.entry_ligands` (live, 2026-08-15) |
| --- | ---: | ---: |
| `B3P` | 235 | **232** |
| `GOL` | 26,117 | **26,004** |
| `EDO` | 17,718 | **17,548** |

`BEN` 361, `JEF` 21, `LZ1` 10, `ZBR` 9, `LFI` 8, `MOV` 7, `N5S` 1 agree exactly.
The higher figures are RCSB-sourced (26,117 and 17,718 also appear in
`CLAUDE.md`'s probe-library argument); the lower are Paperclip's
`pdb_v.entry_ligands`, where `COUNT(*)` and `COUNT(DISTINCT entry_id)` return
the same value, so the view is one row per (entry, component).

**It changes nothing.** Re-running `test_v2.py` with the live Paperclip values
substituted gives byte-identical results — 259/262, 9/9, 9/9, 277/280, 61/70,
0 FP. The ubiquity rule's threshold is ~150 entries and all three components sit
far above it under either source, so no verdict moves. Recorded rather than
silently reconciled, because **a count needs its source named beside it**: these
two numbers are not a discrepancy to resolve, they are two different
measurements, and picking one without saying which would be the error.

## How the sets were drawn — read this before adding a case

**Ground truth is the expected verdict assigned from chemistry knowledge, never
from the classifier's output.** The 262 set is built from sources the classifier
was never shown: every member of `modal_app.COFACTORS` and `NON_LIGANDS`, every
member of `neighbour_precedent.EXCLUDED_LIGANDS`, the four historical failures,
and known true-positive inhibitors, fragments, peptides, steroids and ions. A
handful of labels were corrected after reading the CCD `name` field — each is
annotated in place with the CCD name that justifies it, and each was corrected
*before* the classifier's answer was consulted.

**The held-out 70 were drawn blind** by `ORDER BY MD5(comp_id) LIMIT 70` over
`pdb_v.chemcomps` — deterministic and unrelated to anything in the tuning set —
then adjudicated by name before the rule set was frozen.

> **That draw was re-executed against the live database on 2026-08-15 and
> returns the identical 70 comp_ids** — set-equal to `holdout.json`, nothing
> added, nothing missing. So the sample really is the deterministic blind draw
> it claims to be, and anyone can re-derive it:
>
> ```sql
> SELECT comp_id FROM pdb_v.chemcomps ORDER BY MD5(comp_id) LIMIT 70
> ```
>
> This is the check worth copying: a claim that a set was drawn blind is only
> as good as the ability to re-draw it. Note the draw is stable *because* MD5 is
> deterministic — if `pdb_v.chemcomps` gains rows the sample can change, so
> re-verify rather than assuming. Two defects it exposed
were fixed (`9CP`'s sulfamate misread as a Good's buffer; abamectin and
myxopyronin B misread as lipids by a bare chain-length test).

The 9 context cases each name a real entry and a real accession, because a
crosslinker has **no context-free right answer** — putting one in the chemistry
block would be inventing a label the classifier is not allowed to reach. The two
covalent-inhibitor controls (`MOV`@6OIM, `0WN`@4G5J) are in that block precisely
so the fix cannot buy its false-positive reduction with a false negative on a
real covalent drug.

## One deploy caveat this directory inherits

Running these harnesses creates `tests/__pycache__/`. That is gitignored, but
`managed/druggability-dossier/.gitignore` records that **the ignore does not stop
the deploy zip** — `scripts/deploy.ts` zips each skill directory whole, and a
`.pyc` embeds the compiling machine's absolute source path in `co_filename`.
Fourteen of them once shipped paths like `/Users/<name>/repos/...`, and in one
case a path to a `.env`.

So: **delete `tests/__pycache__/` before a deploy**, until `deploy.ts` gains
`-x "*__pycache__*" -x "*.pyc"` on its `execFileSync` call. That is shared
starter infrastructure and not this skill's to change, but adding a test
directory beside a skill adds one more place for bytecode to accumulate, and
saying so is cheaper than rediscovering it.

## If a number here disagrees with a number quoted elsewhere

**The measured value wins.** Re-run `test_v2.py`, quote what it prints, and fix
the citation. The figures are cited in six places outside this directory; the
list is in `../SKILL.md` under "Where these figures are quoted".
