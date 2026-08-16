---
name: structure-select
description: >
  Picks the structures a pocket scan should run on — classifying holo vs apo by
  actual ligand chemistry rather than by label, assembling an ensemble, applying
  an as-of date cutoff, and finding structural neighbours with Foldseek to
  establish structural-homolog precedent. It does NOT score pockets, does NOT
  predict structures, and does NOT decide tractability.
---

# structure-select

Everything here runs against Paperclip's `pdb_v` views plus Proto's Foldseek.
All queries below were tested; the controls are stated so a regression is
visible.

## Ligand identity — `ligand_filter.py`, not a denylist

**There is no `$EXCL` list any more.** It was a hardcoded set of ~160 comp_ids
paired with a 250-1200 Da window, and that pairing cannot decide holo — see the
failure mode "A comp_id denylist plus a size floor CANNOT decide holo" below,
which is the systemic one on this page.

Holo is now decided from chemistry:

```python
from ligand_filter import classify_ligands, holo_call
v = classify_ligands(["MOV", "GNP", "L44"])   # one Paperclip trip per 40
v["L44"].verdict        # 'lipid_or_detergent'
v["L44"].reason         # '...unbranched aliphatic carbon chain is 21...'
holo_call(["GNP", "GOL"])["is_holo"]          # False
holo_call(["GNP", "GOL"])["determined"]       # True — and check this
```

**Chemistry alone is not enough, and one call changes that.** Pass the entry's
structural context and the classifier can tell a covalent inhibitor from a
crosslinker stapling a peptide ligand:

```python
from ligand_filter import StructureContext, holo_call
# THE HEADER, NOT THE ASSEMBLY — assembly files have no _struct_conn.
ctx = StructureContext.from_mmcif_path("8QFZ_header.cif",
                                       target_accession="Q969D9")
holo_call(["LFI"], context=ctx)["is_holo"]                 # False
holo_call(["LFI"], context=ctx)["polymer_ligand_precedent"]
# [{'via_comp_id': 'LFI', 'modality': 'peptide', 'n_monomers': 12,
#   'sequence': 'CHWLENCWRGFC', ...}]
```

Without the context `LFI` is `druglike` — see "A component covalently bonded
inside a PEPTIDE ligand is not a small molecule" below, which is the measured
false positive that put this here.

`neighbour_precedent.py` uses it. Two rules that come with it:

- **`unknown` is not apo.** A failed CCD lookup returns `unknown` with a
  `lookup_failed` flag. Read `holo_call(...)["determined"]`, and on this axis
  especially — "no drug-like holo among the neighbours" is a real finding, and
  a lookup failure wearing that costume is the worst confusion available.
- **Select candidates in SQL, classify in Python.** The SQL predicate cannot
  express chemistry. Pull the comp_ids unfiltered and decide afterwards.

**Record the entry's OTHER polymer chains when you select it.** Structure
selection is where the information exists and pocket scoring is where it is
needed: `pocket-scan`'s prep keeps every polymer (`het_flag == 'A'`) and strips
only non-polymer HETATM, so a peptide, nanobody, Fab or G-protein partner is
present in the file fpocket scores. Measured on 8QFZ at D=1.6: with the
12-residue Bicycle peptide kept, the headline pocket is 283.6 Å³ and **six of
its ten lining residues are the peptide**; with the peptide chain dropped the
site does not exist and the whole target has one pocket at 147.8 Å³, 15.4 Å
away. Those two numbers straddle rule 4a's proposed 207–242 Å³ band. A census of
this project's 15-target volume evaluation found **14 of 67 structures across 6
of 15 targets** carrying a distinct non-target polymer entity — so emit, per
selected entry, the list of polymer entities with their accessions and whether
each is the target. See the failure mode below and
`polymerfix/PATCHES_modal_app.md`.

## 1. Candidate structures with holo/apo classification

Pull the candidates, then classify them. **The SQL does not decide holo.**

```sql
SELECT s.entry_id, s.resolution, s.exptl_method, s.release_date,
  COALESCE(STRING_AGG(DISTINCT l.comp_id, ' '), '-') AS candidate_ligands,
  COALESCE(STRING_AGG(DISTINCT l.comp_id||':'||COALESCE(l.drugbank_id,'-'), ' '), '-') AS all_ligands
FROM pdb_v.structures_by_accession s
LEFT JOIN pdb_v.entry_ligands l ON l.entry_id = s.entry_id
WHERE s.accession = '<ACC>'
GROUP BY s.entry_id, s.resolution, s.exptl_method, s.release_date
ORDER BY s.resolution NULLS LAST, s.release_date DESC
LIMIT 25
```

Then `holo_call(candidate_ligands.split())` per entry.

**Do NOT add `AND l.comp_type='non-polymer'` to that WHERE clause**, however
tempting. Measured on a 25-entry list: without it the query returns in **6 ms**,
with it the identical query **times out** (>120 s, `[error] Request timed out`).
The column is not usefully indexed, and it does not separate cofactors from
drugs anyway — `ligand_filter` reads the CCD type itself.

**Do not order by the state column.** An earlier version of this query used
`ORDER BY 5 DESC`, where column 5 is the `'HOLO'|'APO'` CASE string. Since
`'HOLO' > 'APO'` alphabetically, every apo entry falls past the row limit and
becomes invisible — TNF-alpha looked like it had no apo structures at all when
it has plenty. Order by resolution, and get the counts from a separate
`GROUP BY state` query rather than by reading the first page.

> **The split quoted here used to be "17 holo / 35 apo" and that is superseded,
> 2026-08-15.** It was a `>=300 Da`-era figure. Under the chemistry classifier
> TNF-alpha is **52 total = 20 holo / 32 apo / 0 undetermined**
> (`fixtures/targets.json`, `_structure_regeneration_2026_08_15`). `targets.json`
> had this logged as `known_inconsistency_left_in_place` pointing at this line;
> it is now fixed here, so that flag can be cleared. **17 + 35 = 52 and
> 20 + 32 = 52 — both sum correctly, which is exactly why the stale pair
> survived review.** An internal consistency check does not detect a superseded
> rule; only re-measurement does.

**Controls: 4OBE must return APO, 6OIM must return HOLO (MOV).** If either
flips, `ligand_filter` has regressed.

Verified on KRAS: 4LYH HOLO 21F(516), 6OIM HOLO MOV(563), 8AZX HOLO OFU(467);
4OBE APO, 6P0Z APO. On TNF-alpha: 9OJO A1CB1(384) 1.36 A, 2AZ5 307(548) 2.1 A.

## 2. As-of filtering

`release_date` is `date`-typed and **100% populated** — 508,687 of 508,687 rows,
spanning 1976-05-19 to 2026-07-01. Future-dated entries exist, so clamp to today
if that matters.

```sql
... WHERE s.accession='<ACC>' AND s.release_date < DATE '<CUTOFF>'
```

Composes with the query above by adding the predicate. KRAS before 2013-01-01
returns 15 entries, earliest 1D8D (2000-02-09), latest 2012-05-23 — and none of
them holo with a drug-like ligand, which is the whole retrospective story.

## 3. Domain-restricted selection

When a target has multiple domains and the drug binds one of them, select on it.
Use **overlap fraction, never containment**:

```sql
WITH ranges AS (
  SELECT r.entry_id,
    GREATEST(0, LEAST(r.max_uniprot_pos,<D1_HI>) - GREATEST(r.min_uniprot_pos,<D1_LO>) + 1)::float
      / (<D1_HI>-<D1_LO>+1) AS d1_frac,
    GREATEST(0, LEAST(r.max_uniprot_pos,<D2_HI>) - GREATEST(r.min_uniprot_pos,<D2_LO>) + 1)::float
      / (<D2_HI>-<D2_LO>+1) AS d2_frac
  FROM pdb_v.uniprot_alignment_ranges r WHERE r.accession='<ACC>'
), dom AS (
  SELECT entry_id,
    CASE WHEN MAX(d1_frac)>=0.7 AND MAX(d2_frac)>=0.7 THEN 'BOTH'
         WHEN MAX(d1_frac)>=0.7 THEN 'D1'
         WHEN MAX(d2_frac)>=0.7 THEN 'D2'
         ELSE 'other' END AS domain
  FROM ranges GROUP BY entry_id
)
SELECT d.domain, COUNT(DISTINCT d.entry_id) n, MIN(s.release_date) earliest
FROM dom d JOIN pdb_v.structures_by_accession s
  ON s.entry_id=d.entry_id AND s.accession='<ACC>'
GROUP BY d.domain ORDER BY 1
```

Domain boundaries come from `uniprot_v.features WHERE feature_type='Domain'`.

Verified on TYK2 (JH1 897-1176, JH2 575-869): JH1 28 entries earliest
**2010-06-02**, JH2 20 entries earliest **2013-04-10**, plus one JH1+JH2 tandem
(4OLI) and three "other" (FERM-SH2 23-583) correctly separated out.

## 4. Foldseek — structural neighbours

```python
from proto_tools.tools.structure_alignment.foldseek.foldseek_search import (
    FoldseekSearchConfig, FoldseekSearchInput, run_foldseek_search)

result = run_foldseek_search(
    FoldseekSearchInput(structure="/path/to/query.pdb"),
    FoldseekSearchConfig(search_mode="remote", databases=["pdb100"],
                         mode="3diaa", timeout_seconds=900.0))
```

One required field, `structure`: a `Structure` object, **a file path**, or raw
PDB/CIF text. **It does not accept a PDB ID or a URL** — resolve those to a file
first.

Verified: 6OIM chain A against `pdb100` returned **992 hits in 13.3 s**, no
database download, no local disk. `mode='tmalign'`: 8.0 s, 989 hits. Remote mode
POSTs to `search.foldseek.com`; no Modal, no GPU.

Databases: `pdb100, afdb50, afdb-swissprot, afdb-proteome, mgnify_esm30,
gmgcl_id, BFVD, cath50, bfmd`. Use `pdb100` when you need PDB IDs back —
afdb/BFVD return AlphaFold/UniProt accessions and need a different parse.

### 4b. Foldseek MULTIMER — the right tool for any assembly

Use this whenever the input has more than one chain. Signature, verbatim from
`proto_tools/tools/structure_alignment/foldseek/foldseek_multimer_search.py`:

```python
@tool(
    key="foldseek-multimer-search",
    label="Foldseek Multimer Search",
    category="structure_alignment",
    input_class=FoldseekMultimerSearchInput,
    config_class=FoldseekMultimerSearchConfig,
    output_class=FoldseekMultimerSearchOutput,
    description="Search Foldseek multimer (complex) structural homology — remote (server) or local (CLI)",
    uses_gpu=False,
    example_input=example_input,
    cacheable=True,
)
def run_foldseek_multimer_search(
    inputs: FoldseekMultimerSearchInput,
    config: FoldseekMultimerSearchConfig,
    instance: Any = None,
) -> FoldseekMultimerSearchOutput:
```

`FoldseekMultimerSearchInput` has one field, `structure: Structure` — "Multi-chain
query complex (Structure object, file path, or raw PDB/CIF string)". Config is
`search_mode` (`'remote'` default | `'local'`), `databases` (default
`['pdb100']`), `mode` (`'3diaa'` default, **wire-encoded as `complex-{mode}`**),
`poll_interval_seconds`, `timeout_seconds`, plus a local-only block (`local_db`
required, `evalue`, `sensitivity`, `max_seqs`, `alignment_type`,
`tmscore_threshold`, `lddt_threshold`, `num_threads`, `use_gpu`). Output is
`ticket_id`, `hits`, `num_hits`, `databases_queried`, `result_url`.

```python
from proto_tools.tools.structure_alignment.foldseek.foldseek_multimer_search import (
    FoldseekMultimerSearchConfig, FoldseekMultimerSearchInput,
    run_foldseek_multimer_search)

result = run_foldseek_multimer_search(
    FoldseekMultimerSearchInput(structure="/path/8DYG-assembly1.cif"),
    FoldseekMultimerSearchConfig(search_mode="remote", databases=["pdb100"],
                                 mode="3diaa", timeout_seconds=900.0))
```

**No Modal, no GPU, no `MODAL_PROFILE`.** `uses_gpu=False`, and the tool's own
`local_execution_reason` says remote mode "queries the public Foldseek server
over HTTP, so device=... would only add a network hop". It POSTs to the same
`search.foldseek.com/api/ticket` endpoint as the single-chain tool with
`mode='complex-3diaa'`; `/foldmulti` is the web UI path, not the API path.

Verified: IL-17A 8DYG assembly1 (2 chains) 863 rows in **405 s**; TNF-alpha
1A8M assembly1 (3 chains) 6,891 rows. **Read `result_url` and parse the raw m8
yourself** — the wrapper's parser drops the columns that make it a multimer
result. See the failure modes.

## 5. Neighbour precedent — `neighbour_precedent.py`

This is wired up. Run the module in this directory rather than reassembling the
steps by hand:

```bash
$PROTO_PY neighbour_precedent.py <structure.pdb> <ACCESSION> \
    [--max-neighbours 25] [--min-alignment-length N] [--cache hits.json] \
    [--multimer auto|yes|no]
```

or `from neighbour_precedent import neighbour_precedent` under the proto-tools
python. It needs `paperclip` on PATH and `PAPERCLIP_API_KEY`; the module reads
`/Users/bb/repos/claude-agent-starter/.env` by default (`env_file=`).

The procedure it implements, in order:

0. **Count the chains and pick the search.** `>1` polymer chain routes to
   `foldseek-multimer-search` (`mode='complex-3diaa'`), one chain to
   `foldseek-search`. This is not a preference: `foldseek-search` aligns each
   chain of a multi-chain file *independently* and returns the union, with no
   complex assignment and no complex TM-score, so on an oligomer it answers a
   different question. `--multimer no` pins the old path; a multimer failure
   falls back to it and records the error in
   `foldseek.multimer_attempted_and_failed`. Chains are counted from ATOM
   records only — a ligand in its own chain id must not trigger the multimer
   path.
1. **Foldseek against `pdb100`**, remote. Single-chain path: `mode='3diaa'`
   plus a second `mode='tmalign'` search joined on `target_id` for TM-scores.
   Multimer path: **no second search** — the complex TM-score is already in the
   raw m8 (column 21), which the module re-downloads from `result_url`.
2. **Filter to fold-not-sequence neighbours**: `sequence_identity < 0.30 AND
   alignment_length >= 120`. Both halves matter — without the identity ceiling
   this is `family_precedent` wearing a Foldseek hat, and without the length
   floor it is a motif match.
3. **Resolve the aligned CHAIN(s) to accessions**, not the entry's accessions —
   see the failure mode below, it changes answers. On the multimer path a
   neighbour is a *set* of matched chains and all of them are resolved.
4. **Ask whether those proteins have drug-like holo entries**, classifying
   every candidate ligand with `ligand_filter`, reported as an entry-level
   upper bound *and* a single-protein-entry floor, with `pdb_v.entries.title`
   attached and a `rejected_ligands` line per neighbour saying *why* each
   candidate was dropped.

Output keys worth knowing: **`search_path`** (`multimer` | `single_chain` — read
this first; it decides whether the block is evidence about an interface),
`n_query_chains_in_file`, `n_query_chains_observed_in_result` (None on the
single-chain path — the wrapper drops the column that would tell you) and
`chains_assembled_into_complexes`, `neighbours` (per entry —
`tm_score` with `tm_score_kind`, `chains`, `n_query_chains_matched`,
`probability`, `evalue`, `has_druglike_holo`, `ligands`, `ligand_names`,
`attribution`, `title`), `neighbour_accessions` (per accession across *all* its
PDB entries, with `holo_titles`), `filter.auto_relaxed`, and `caveats`.

Fills `structural_neighbour_precedent`. Keep it separate from
`family_precedent` — fold neighbours and sequence family are different signals
and are allowed to disagree.

### Measured on both calibration targets

| | KRAS 6OIM_A (P01116) | IL-17A 8DYG asm1 (Q16552) |
| --- | --- | --- |
| Foldseek hits | 992 | 283 |
| passing the filter | 285 | **2** (relaxed to 81) |
| neighbourhood | Rab / Ran / Rac / Ypt GTPases, TM 0.73-0.89 | cystine-knot growth factors: IL-25, VEGF-A/B/C/F, NGF/NT-3, BMP-2, PDGF-B, TGF-beta2, sclerostin, coagulogen; TM 0.34-0.78 |
| entries apo / holo, old MW+denylist | **24 / 1** (4PHH) | **24 / 1** (4EC7) |
| the one "holo" | 2UK = a GppNHp analog | L44 = a 625 Da diacylglycerol |
| with `ligand_filter` | — | **25 / 0** on the multimer path |
| honest read | **no small-molecule precedent** | **no small-molecule precedent** |

Both of those "holo" hits were denylist leaks and `ligand_filter` now rejects
both, so the defensible count on both targets — **0 of 25** — is what the tool
returns rather than what a reader has to work out. Keep reading `ligand_names`
anyway: the classifier's remaining gap is crystallisation additives.

A third target, and the one the single-chain limitation was most likely to
break:

| | TNF-alpha, single-chain (1TNF asm1) | TNF-alpha, multimer (1A8M asm1) |
| --- | --- | --- |
| chains in | 3 | 3 |
| rows | 1,010 | **6,891** |
| query chains present | job_A/B/C (349/320/341) | job_A/B/C (2,119/2,113/2,659) |
| complex assignments | none — no such column | **4,244; 1,206 of size 3** |
| distinct entries | 302 | 361 |
| passing strict 120 | 236 entries | 279 entries |
| entries apo / holo | **25 / 0** | **25 / 0** |
| neighbourhood | TL1A, LIGHT, RANKL, TRAIL/DR4/DR5, CD40L, FasL, LT-alpha | TL1A, LIGHT, RANKL, TRAIL/DR4/DR5, CD40L, FasL, zebrafish TNF |

TNF-alpha needs no relaxation — the protomer is ~157 residues, so the strict
120 floor holds and 236 of 302 entries pass. Both paths return the **TNF
superfamily**, and **zero holo among the 25 carried neighbour entries**.

**But do NOT read that as "the TNF superfamily has no small-molecule
precedent" — an earlier draft of this section said exactly that and it is
wrong.** The 25 carried entries are apo; the *accession* sweep behind them is
not. Of the 8 neighbour accessions, **CD40L (P29965) returns 1 holo of 8
structures: `3LKJ`, ligand `LKJ`, attribution `unambiguous`** (one protein, one
polypeptide entity — no partner to misattribute to), titled:

> Small Molecule Inhibition of the TNF Family Cytokine CD40L

`LKJ` is a biphenyl-containing drug-like molecule, and this is the one hit on
either calibration target that survives every filter we have: right chemistry,
right protein, right entry, and a title that states the intent outright. **It is
genuine small-molecule precedent on a TNF-superfamily cytokine trimer**, which
is the single most decision-relevant thing this axis has produced.

Two honest limits on it. Re-querying `3LKJ`'s ligand row to confirm the exact
MW hit failure signature #3 (`vsh: cd: /papers/: Permission denied`) on three
consecutive attempts, so the molecular weight is **unconfirmed**; the comp_id,
the accession, the attribution flag and the title all came back cleanly from two
independent retrievals. And it is **accession-level, not neighbour-level** — it
did not appear among the 25 carried entries, so a reader looking only at
`neighbours` would never see it. Read `neighbour_accessions`, not just
`neighbours`.

The multimer path's 1,206 complex assignments **of size 3** are the thing the
single-chain path cannot produce — genuine trimer-to-trimer matches, ranked by
complex TM-score (TL1A 0.878 at the top, which matters because TL1A is itself a
target in this pipeline).

### Multimer versus single chain on the same IL-17A input

`8DYG-assembly1.cif`, Q16552, both paths, same filter:

| | single-chain | multimer |
| --- | --- | --- |
| rows | 283 | 863 |
| query chains present | job_A 144 + job_B 139 | job_A 433 + job_B 430 |
| target chains of the SELF hit | `8dyg_A` only | `8dyg_A` **and** `8dyg_B` |
| complex assignments | none | 570 (293 of size 2) |
| distinct entries | 125 | 174 |
| passing strict 120 | 2 | 2 rows / **1** entry |
| relaxed floor | 67 | **67 — same** |
| entries after relaxation | 81 | 137 |
| carried (top 25) | 25 | 25 |
| apo / holo, old MW+denylist | 24 / 1 | 23 / 2 |
| those "holo" | 4EC7 `L44` diacylglycerol | 4EC7 `L44` + 4XPJ `LPY` lysophospholipid |
| apo / holo, `ligand_filter` | 24 / 1 (2GNN `BEN`, benzamidine) | **25 / 0** |
| defensible small-molecule holo | **0 of 25** | **0 of 25** |
| entries shared with the other path | 105 | 105 (69 multimer-only, 20 single-chain-only) |
| neighbourhood | cystine-knot superfamily | **cystine-knot superfamily** |

The 42 chain accessions across the whole 137-entry multimer neighbourhood are
VEGF-A/B/C/D and PlGF, NGF / NT-3 / NT-4 / BDNF, PDGF-A/B, TGF-beta1/2, GDF-15,
BMP-2, GDF-5, sclerostin, noggin, norrin, the glycoprotein hormones, von
Willebrand factor, AMH, and IL-25 — the cystine-knot superfamily, which is the
same answer the single-chain search gave. Every one is an antibody-drugged PPI
target with no small-molecule holo structure.

**So on IL-17A the single-chain limit was harmless, and the precedent call does
not change.** Say that plainly rather than manufacturing a difference. What
does change is what the result *ranks*: multimer orders by complex TM-score, so
true homodimers (IL-25 at 0.74, then the neurotrophins) come first and
single-protomer matches sink. That is the right ordering for an interface site,
and it is a better-justified neighbourhood, not a different conclusion.

The reason to keep the multimer path is therefore **not** that it changed this
answer. It is that on this target we can *check* that it did not, and on the
next oligomer we cannot.

IL-17A is the interesting one and the prediction held. Its fold neighbourhood is
the cystine-knot superfamily — VEGF, NGF, BMP, PDGF, TGF-beta — every one of
them a PPI target approached with antibodies, none with a small-molecule holo
structure. Meanwhile IL-17A *itself* has 44 structures and 20 holo, the
macrocycle series. So `target_precedent` is strong and
`structural_neighbour_precedent` is empty **on the same target**. Report the
disagreement; it is the informative thing on the page.

The underlying family-level query, if you need it standalone:

```sql
WITH hits AS (
  SELECT DISTINCT pe.uniprot_accession AS acc FROM pdb_v.polymer_entities pe
  WHERE pe.entry_id IN (<foldseek hit pdb ids>) AND pe.uniprot_accession IS NOT NULL),
all_s AS (
  SELECT DISTINCT h.acc, s.entry_id FROM hits h
  JOIN pdb_v.structures_by_accession s ON s.accession = h.acc)
SELECT a.acc, COUNT(DISTINCT a.entry_id) n_struct,
       COUNT(DISTINCT a.entry_id) FILTER (WHERE l.comp_id IS NOT NULL) n_holo,
       COALESCE(STRING_AGG(DISTINCT l.comp_id,' '),'-') druglike
FROM all_s a
LEFT JOIN pdb_v.entry_ligands l ON l.entry_id = a.entry_id
     AND l.comp_id IN (<the comp_ids ligand_filter classified `druglike`>)
GROUP BY a.acc HAVING COUNT(DISTINCT a.entry_id) FILTER (WHERE l.comp_id IS NOT NULL) > 0
ORDER BY 3 DESC LIMIT 20
```

**This bare form is an upper bound only** — it has no chain attribution and no
title check. Use the module, or read the two failure modes below first.

## 6. Transferred-homolog site anchoring — `homolog_transfer.py`

`neighbour_precedent` tells you a fold neighbour exists and carries a drug-like
ligand. This turns that into **coordinates**: superpose the donor onto the
target and carry its ligand across, producing a ligand-free site anchor for a
target that has no ligand of its own.

```bash
python3 homolog_transfer.py <TARGET_PDB> <DONOR_PDB> \
        [--donor-ligand COMP] [--reference-ligand COMP]
```

Pure stdlib — no `gemmi`, no `numpy`, no `fpocket`, no `paperclip`. It fetches
mmCIFs from RCSB (cached under `$STRUCTURE_SELECT_CACHE`, default
`~/.cache/structure-select/mmcif`) and chem-comp rows from RCSB's REST API.
`--reference-ligand` names a ligand **in the target** and is used only to
validate afterwards; it is never an input to a guard.

### Why this is the only ligand-free anchor worth building on

An anchor-agreement test built **sixteen** ligand-free anchors across
TNF-alpha, IL-17A, NLRP3 and S1PR1. Four could not be built. **Four of sixteen
found the ligand site:**

| ligand-free anchor | found the site |
| --- | --- |
| **transferred homolog** | **2 of 3 where constructible** |
| interface | 1 of 6 |
| symmetry axis | 1 of 3 |
| annotated function | **0 of 4** |

Four targets is not a rate. What it supports is exactly two statements:
transferred homolog is the one worth building on, and **none of the four is
safe unaided**. The other three are handled by reporting rather than selection
— see "The other three anchors" below.

### The three guards

Every one closes a transfer that ran, produced a confident answer, and was
wrong.

| guard | what it requires | the failure it closes |
| --- | --- | --- |
| 1 `donor_ligand_druglike` | the donor's ligand passes `ligand_filter` | NLRP3's only constructible donor was NOD2 with **ADP**. The transfer was *excellent* and it selected the nucleotide lobe. |
| 2 `domain_attribution` | ≥50% of the donor ligand's 4.5 Å contact shell lies in the **aligned region** | 7KRZ's bortezomib is on the correct LONP1 chain at auth 768-898 while the NACHT-aligned region is auth 506-721. Right chain, wrong domain. |
| 3 `alignment_and_sterics` | TM ≥ 0.5, RMSD ≤ 5.0 Å, **zero backbone clashes** | IL-2 forced onto IL-17A at TM 0.254-0.274 put the ligand 21.59 Å away and inside the protein, silently. |

**Guard 1 rejects everything that is not `druglike`** — `cofactor`,
`ion_or_solvent`, `lipid_or_detergent`, `sugar_or_glycan`,
`crystallisation_additive`, `peptide_or_polymer`, `polymer_conjugate`, and
`unknown`. A donor we cannot classify is a donor we cannot vouch for, and a
CCD lookup that *failed* is reported as `lookup_failed`, never as a pass.

**Guard 2 is a third level of attribution.** This file already documents that
`entry_ligands` cannot attribute a ligand to a chain, and that the same bug
exists one level up on the protein side. This is the level below both: **a
chain is not a domain.** Chain-level attribution passes 7KRZ.

**Guard 3's floor is the literature fold cut, not a number read off our
cases.** TM ≥ 0.5 is where two structures are held to share a fold; our data
confirms it discriminates (IL-2/IL-17A 0.203 measured here) but n=2 successes
cannot calibrate a threshold. The RMSD ceiling is deliberately **loose at
5.0 Å**: TNF's donor fit at 1.35 Å and S1PR1's at a mediocre 3.04 Å and still
landed — but S1PR1's donor is a ~45%-identical **paralogue**, so 3.04 Å is not
a safe general floor for a distant neighbour. Anything above 2.0 Å is accepted
and flagged `alignment_marginal`.

### The steric check is the strong guard — and it must be counted on BACKBONE

Measured, and this is the finding that changed the design:

| case | backbone clashes | side-chain clashes | ref. distance |
| --- | ---: | ---: | ---: |
| control 2AZ5 ← 2AZ5 `307` | **0** | 0 | **0.06 Å**, Jaccard 1.000 |
| TNF 2AZ5 ← 3LKJ `LKJ` | **0** (min 2.40 Å) | 4 (min 1.15 Å) | 8.73 Å, Jaccard 0.50 |
| S1PR1 3V2Y ← 4Z34 `ON7` | **0** (min 2.50 Å) | 11 (min 0.16 Å) | **1.79 Å**, Jaccard 0.50 |
| NLRP3 7ALV ← 5IRN `ADP` | **0** | 0 | **0.46 Å**, Jaccard 1.000 |
| IL-17A 9SQX ← 1M48 `FRG` | **7** (min 0.65 Å) | 5 | **21.41 Å**, Jaccard 0.00 |

**On a total-heavy-atom count the S1PR1 success is worse than the IL-2
failure** — 11 clashing atoms of 38 at 0.16 Å against 12 of 33 at 0.65 Å —
while landing 1.79 Å from the reference ligand rather than 21.41 Å away. Any
total-atom threshold that rejects the failure also rejects the second success.
On backbone the two classes separate completely: 0/0/0 against 7.

The reason is physical and is already in `CLAUDE.md` rule 5. A ligand
overlapping a **side chain** is a rotamer problem — two crystal structures of
the same pocket disagree about rotamers, so a transferred ligand nearly always
overlaps some — and repacking resolves it. A ligand overlapping N, CA, C or O
is not resolvable by anything. So: **backbone contact voids the transfer;
side-chain contact is reported and flagged `sidechain_occlusion`**, which is a
finding rather than noise, because rule 5 puts that mechanism's prognosis at
micromolar-at-best.

### Nothing fires silently

`TransferResult.assert_reportable()` raises unless all seven of
`donor_pdb_id`, `donor_ligand`, `tm_score`, `rmsd_a`, `aligned_length`,
`clash_count` and `backbone_clash_count` are present — **on rejection as well
as acceptance**. A transferred-homolog anchor without its provenance is not
weak evidence, it is unusable, because a reader has no way to discount it. Call
it before writing an anchor into a dossier.

The block also carries `chain_map_donor_to_target`, `local_refit`,
`sidechain_clash_count`, `clash_by_chain` and `alignment_marginal`.

### Three things the reproductions forced, which are not obvious

**An asymmetric unit is not a biological unit, and a chain map must not span
two.** 2AZ5 is **two TNF dimers** (A+B and C+D; interfaces 24 and 23 residues,
against 9 for the packing contact between them), not one tetramer — SPD304
occupies the third protomer's place, which *is* the mechanism. Allowing a map
to cross them produced a transfer that passed all three guards and sat 8.79 Å
from SPD304 with a contact-shell Jaccard of 0.171. Chains are grouped into
units at ≥15 contacting residues; measured, every biological interface in
2AZ5/3LKJ/9SQX is ≥18 and every packing contact ≤12.

**Sequence identity cannot rank chain maps for a homo-oligomer** — every
pairing scores identically by construction, which is the same defect
`CLAUDE.md` names `site_signature_unreliable_homooligomer`. Maps are ranked on
the TM-score **numerator** (the number of structurally equivalent residues),
which carries no normalisation and so compares maps of different size, and
rewards coverage rather than collapsing to one chain. Ties within 90% are
broken on donor-ligand shell coverage: on TNF, ranking on mass alone left the
ligand's principal contact chain (donor C, 238 contacting atoms against 85 and
63) out of the fit entirely.

**The transform that moves the ligand is re-fitted locally.** A global
superposition minimises error everywhere and therefore nowhere in particular.
`tm_score` and `rmsd_a` stay global — that is the fold match, and guard 3's TM
floor applies to it — but the coordinates come from a re-fit over the aligned
pairs within 15 Å of the donor ligand, whose RMSD is reported separately as
`local_refit.local_rmsd_a`. Reporting one and transferring on the other is the
same mismatch rule 5 records for the two C-alpha displacement protocols.

### What "TNF landed at 0.00 Å" actually means

The anchor-agreement test reports TNF at **0.00 Å, Jaccard 0.615**. That is a
**pocket-selection** distance: the transferred anchor and the SPD304 anchor
selected the *same fpocket pocket*, so their centroids coincide exactly. It is
not a ligand-centroid distance, and this module cannot produce it — fpocket
lives in `pocket-scan`. What is reproducible here is the guard verdict and the
geometry the guards read. The ligand-centroid figure this module measures for
that case is 8.73 Å with a contact-shell Jaccard of 0.50, which is consistent
with selecting one large axial cavity and is **not** the same claim. Quote them
apart.

### The other three anchors: report, do not select

These belong to `pocket-scan`'s `_annotate_pocket_labels`, which already
treats anchoring as an annotation rather than an election. The
anchor-agreement test says what each label is worth, and none of the three may
choose a pocket:

- **Interface — 1 of 6.** Its one success depended on symmetrising the receptor
  epitope across the 2-fold; **as deposited the same anchor is 14.67 Å off**.
  And two defensible pocket-selection rules on the *same* epitope pick pockets
  4.11 Å and 14.19 Å from the site. So **the selection rule must be declared
  before measurement**, and when the interface anchor disagrees with another
  anchor, **report both pockets rather than choosing**.
- **Symmetry — 1 of 3.** TNF's C3 axis carries **five** buried on-axis
  cavities and no ligand-free rule picks among them: the pre-declared "widest
  cavity" rule returns one **22.37 Å** away while the runner-up at **7.86 Å**
  independently reproduces the 7.7 Å figure from the withdrawn residue-number
  matcher. On IL-17A the correct cavity wins by **0.18 Å against 2.06 Å** and
  the decoy is *the other real ligand site*. So it may fire only when the
  nearest-axis buried cavity is **unique by a stated margin**; otherwise emit
  `ambiguous` and list every candidate.
- **Annotated function — 0 of 4**, and on two of the four the canonical form
  does not exist so a **surrogate** was scored. Do not build on it.

### Reproducing it

`tests/test_homolog_transfer.py` runs all four measured cases plus the
self-transfer control and the 7KRZ guard-2 case. **It is network-dependent,
unlike every other harness in that directory** — see `tests/README.md`.

## Failure modes

### A comp_id denylist plus a size floor CANNOT decide holo — use `ligand_filter.py`

This is the systemic one. Several failure modes below are symptoms of it.

Every "is this entry holo?" decision in this pipeline was made two ways, and
both are wrong:

1. **`comp_id` against a hardcoded exclusion set** — the former §1 `$EXCL`,
   `neighbour_precedent.EXCLUDED_LIGANDS`, `modal_app.COFACTORS`.
2. **A heavy-atom or MW floor** — `modal_app.DRUGLIKE_MIN_HEAVY_ATOMS = 18`,
   `neighbour_precedent.MW_MIN/MW_MAX`.

Both fail, and **both fail in the flattering direction**: they invent holo
structures that do not exist, which inflates apparent druggability. Four
measured wrong answers, on four different targets:

| target | reported | the "ligand" | truth |
| --- | --- | --- | --- |
| CD20 | 3 holo | `Y01` cholesteryl hemisuccinate, phosphatidylcholine — cryo-EM sample additives | **0 holo** |
| KRAS fold neighbours | 1/25 holo (4PHH) | `2UK`, a GppNHp analog — a nucleotide cofactor with a comp_id nobody listed | **0 of 25** |
| IL-17A fold neighbours | 1/25 holo (4EC7) | `L44`, a 625 Da diacylglycerol — clears an 18-heavy-atom floor because it is a big greasy lipid | **0 of 25** |
| NLRP3 | ADP entry called holo | `ADP`, ~27 heavy atoms in the NACHT domain | **apo** |

The pattern: **a hardcoded comp_id list cannot enumerate chemistry, and
molecular size does not distinguish a drug from a lipid.** Both a 625 Da
diacylglycerol and a 625 Da inhibitor clear a size gate; only chemistry
separates them. Extending the list is not a fix — it is the bug, applied again.

**Use `ligand_filter.py`** (same directory). It classifies on chemistry read out
of `pdb_v.chemcomps` — the CCD, via Paperclip, nothing external:

```python
from ligand_filter import classify_ligand, is_druglike_ligand, classify_ligands, holo_call

classify_ligand("L44").verdict     # 'lipid_or_detergent'
is_druglike_ligand("2UK")          # False
classify_ligands([...])            # batch: ONE round trip per 40 comp_ids
holo_call(["MOV", "GDP", "MG"])    # {'is_holo': True, 'druglike_ligands': ['MOV'], ...}
```

Verdicts: `druglike`, `cofactor`, `lipid_or_detergent`,
`crystallisation_additive`, `sugar_or_glycan`, `ion_or_solvent`,
`peptide_or_polymer`, **`polymer_conjugate`**, `unknown`. Only `druglike` is
evidence of a bindable site by a small molecule. Every verdict carries a
`reason` string and the `evidence` it rests on.

`polymer_conjugate` is the only verdict that needs a `StructureContext` and the
only one `classify_record` can never return — it means "this component is a
covalent constituent of a polymer, and here is which one", and
`evidence["conjugate_of"]` names that polymer so its precedent can be filed
under the right modality.

What it keys on, none of which the old code read:

- **`_chem_comp.type`.** Decisive on its own for polymers. It is how 6OIM's GDP
  is caught — the CCD types GDP as `RNA linking`, not `non-polymer`. Peptide-
  and saccharide-linking components are never small-molecule evidence.
- **Element composition.** Drug-like ligands are overwhelmingly N-containing; a
  pure C/H/O molecule with a long aliphatic run is a lipid.
- **Nucleotide signature** — purine base + ribose + phosphate. Catches ADP, ATP,
  GDP, GTP, GNP and analogs like `2UK` without naming any of them.
- **Sterol signature** — steroid nucleus *plus an aliphatic side chain*. The
  side chain is what separates cholesterol (tail of 8) and bile salts from a
  steroid DRUG like dexamethasone (tail of 2), which stays `druglike`.
- **Phospho-headgroup plus acyl chains** — phosphatidylcholines and detergents.
- **Any free phosphate ester** — the signature of endogenous metabolites.
- **Longest unbranched alkyl chain**, as a fraction of the molecule.

- **The structural context of the entry**, when one is supplied: `_struct_conn`
  covalent linkages and `_entity_poly_seq` membership. Four rules, and the one
  that matters most is the one that does *nothing* — a component making ONE
  covalent bond to the TARGET polymer is a covalent inhibitor and keeps its
  chemistry verdict untouched.

Measured, not claimed. **Re-measured 2026-08-15 after the first false positive
against the held-out result; both original figures are unchanged and the new
cases are additive.**

**Reproduce every figure in this section with one command:**

```bash
python3 .claude/skills/structure-select/tests/test_v2.py
```

Offline, pure stdlib, no Paperclip call and no network — so it does not depend
on the row cap or the `cli_cwd` bug. The harnesses, the cached
`pdb_v.chemcomps` rows and the 7 context mmCIFs are in
`.claude/skills/structure-select/tests/`; read `tests/README.md` for how the
sets were drawn. **Before 2026-08-15 these harnesses existed only in session
scratch under `/private/tmp`, which meant none of the figures below traced to
anything in this checkout** — the same failure mode that let the retracted
volume separation survive. If a quoted number ever disagrees with what the
harness prints, the harness wins.

- **259/262 = 98.9%** on the original ground-truth set, which includes the four
  failures above, every member of `modal_app.COFACTORS` and `NON_LIGANDS`, and
  every member of `neighbour_precedent.EXCLUDED_LIGANDS` — none of which the
  classifier was told. Misses: `BTN` (biotin → `druglike`), `ACE`/`NH2` (capping
  groups). **`NH2` and `ACE` are now correct when a context is supplied** — the
  CCD lists them in `_entity_poly_seq`, so they are residues, not ligands.
- **277/280 = 98.9%** on the ground-truth set extended with 9 chemistry cases
  and 9 named-entry context cases (see below). The same three misses; no new
  ones.
- **61/70 = 87.1%** on a blind held-out sample from `pdb_v.chemcomps`, with
  **0 false positives** — unchanged by every rule added since. Nothing that was
  really a cofactor, lipid or additive was called drug-like. **That sample did
  not contain TNF 5UUI's `MTN` spin label, which is a genuine false positive no
  chemistry can fix** — see "A ligand can be disqualified by why it is there"
  below. Quote the zero with that boundary attached, and with the second one
  now: it also did not contain a peptide-conjugated crosslinker, and `LFI` was
  the first measured false positive against it.

Known false negatives — deliberate, since a false positive is what caused all
four bugs: nucleoside/SAM-analog inhibitors and bisphosphonates → `cofactor`;
metallodrugs → `cofactor`; long-tailed natural-product antibiotics → `lipid`;
glycosylated natural products → `sugar_or_glycan`; peptidomimetic drugs typed
`peptide-like` → `peptide_or_polymer`. Read `evidence` and `flags`, not just
`verdict`, if any of those classes matter to the target in hand.

Two behaviours to respect:

- **`unknown` is not `apo`.** An unclassified ligand is not evidence of a site,
  so `is_druglike_ligand` returns False — but do not write "apo" on the strength
  of it. Check `holo_call(...)["determined"]`.
- **A lookup failure is not a CCD miss.** Paperclip's endpoint intermittently
  exceeds its statement timeout. Those verdicts carry the flag `lookup_failed`
  and appear in `holo_call(...)["undetermined"]`. Reporting such an entry as apo
  reintroduces the original bug in a new place.

It has no dependencies outside the standard library — no RDKit, which
`pocket-scan`'s Modal image does not have — so the verdict cannot vary with the
environment it is evaluated in.

### `comp_type` does not separate cofactors from drugs

**Still true, and now a reason not to filter on `comp_type` at all** — see the
section above; `ligand_filter` reads the CCD type itself and decides on
chemistry.

GDP and UDP are `'RNA linking'`, but **ATP (507 Da), GTP (523), GNP/GppNHp (522),
NAD (663), FAD (786), COA (768), HEM (617) are all `'non-polymer'`** and sit
inside any sensible MW window. GNP is the KRAS *active-state* analog, so a
`comp_type` filter both admits every GppNHp KRAS structure as HOLO *and* drops
the `'RNA linking'` nucleotides it should have caught. It discriminates in
neither direction.

There is a second, unrelated reason: putting `comp_type` in a WHERE clause is
**slow enough to fail**. Measured on the same 25-entry list, `SELECT DISTINCT
l.comp_id ... WHERE l.entry_id IN (...)` returns in 6 ms and the identical query
with `AND l.comp_type='non-polymer'` times out past 120 s.

### `drugbank_id` is not a druglikeness signal

Glycerol is `DB09462`. Sulfate is `DB14546`. GDP is `DB04315`. **Never gate HOLO
on `drugbank_id IS NOT NULL`.** Report it; do not filter on it.

### Containment tests silently miss domain structures

Selecting JH2 entries with `min_uniprot_pos>=560 AND max_uniprot_pos<=880` gives
earliest 2015-03-18 — wrong by two years. The true earliest, **3ZON
(2013-04-10)**, spans 541-873 and overhangs the domain boundary, so containment
drops it. Use overlap fraction.

### `entry_ligands` cannot attribute a ligand to a chain

It is keyed at **entry** level, and its `entity_id` is the ligand's own
nonpolymer entity, not the protein chain it touches (verified: 2AZ5's ligand is
entity `2AZ5_2`, the protein is `2AZ5_1`). In a multi-protein complex a ligand
bound to the *partner* still counts toward the entry.

This inflates neighbour precedent badly — RAN (P62826) showed 36/139 holo, but
the leptomycin-class ligands bind exportin, not RAN. **Treat family-level holo
counts as an upper bound** and check `pdb_v.entries.title` before believing them.
Single-protein targets like KRAS and TNF-alpha are unaffected.

`neighbour_precedent.py` implements this rather than merely warning about it.
Every holo count comes back twice:

- `n_holo_entry_level` — the naive count, the **upper bound**;
- `n_holo_single_protein_entries` — only entries with one distinct UniProt
  accession *and* one polypeptide entity, where attribution is unambiguous;
- `holo_titles` — up to three PDB IDs with ligands, an `attribution` flag and
  the entry title, so a reader can adjudicate the gap.

Re-measured on the KRAS neighbourhood, the gap is enormous:

| accession | structures | holo entry-level | holo single-protein |
| --- | --- | --- | --- |
| P62826 RAN | 139 | **36** | **0** |
| P63000 RAC1 | 80 | 11 | 0 |
| P32939 Ypt7 | 5 | 3 | 3 |

**But `ambiguous` does not mean `wrong`, and the flag must not be used as a
filter.** Two cases from the same run make the point:

- P62826 RAN's 36 are `4GMX / 4GPT / 4HAT` — "KPT185 / KPT251 / Leptomycin B in
  complex with CRM1-Ran-RanBP1". The ligand is a CRM1 inhibitor. RAN is a
  bystander and the count is **spurious**.
- P63000 RAC1's `5QQE / 5QQG` are a PanDDA fragment screen on a **RAC1-Kalirin
  complex**. Same flag, and the fragments are **genuine RAC1 precedent**.

Identical `ambiguous_multiprotein` label, opposite truth. The flag exists to
route a human to the title, not to decide.

### The same attribution bug exists one level up, on the protein side

A Foldseek hit names a **chain**. Taking the entry's accessions wholesale
imports that chain's crystallisation partners as if they were fold neighbours.

Verified: searching IL-17A returned 2XAC, "Structural Insights into the Binding
of VEGF-B by VEGFR1". Foldseek matched **chain A = P49765, VEGF-B** — a
cystine-knot cytokine, correctly. The entry also contains chains C/X = **P17948,
VEGFR1**, a receptor tyrosine kinase. Reading entry accessions pulled VEGFR1
into a cytokine's fold neighbourhood, and with it **3HNG's genuine kinase
inhibitor** — the only real drug-like small molecule in the whole IL-17A result,
and it did not belong there.

Resolve the chain instead:

```sql
JOIN pdb_v.polymer_entities pe
  ON pe.entry_id = <hit entry> AND pe.auth_asym_ids @> to_jsonb('<chain>'::text)
```

`auth_asym_ids` is `jsonb`. Use `@> to_jsonb(x::text)` rather than the `?`
operator — `?` is a placeholder in several drivers.

Fixing this dropped KRAS from 29 neighbour accessions to 19 (7 holo accessions
to 4, correctly losing Rabphilin P47709, the *effector* in 1ZBD) and IL-17A from
23 to 17 (6 holo to 4). The module reports `chain_accession_resolved` and lists
what it excluded in `other_accessions_in_entry`; when the chain cannot be
matched it falls back to entry accessions and says so.

### Foldseek's `evalue` and `bit_score` are mislabeled in remote mode

The public server emits a wide m8 — **21 columns single-chain, 26 multimer**
(an earlier note here said 17; re-measured on IL-17A 8DYG, both paths) — but the
parser reads columns 10 and 11 per the standard 12-column layout. Confirmed
against the raw m8 on both endpoints:

| field | actually holds | best hit | worst hit |
| --- | --- | --- | --- |
| `hit.evalue` | Foldseek **probability** (higher = better) | 1.000 | 0.045 |
| `hit.bit_score` | the **true E-value** (lower = better) | 5.6e-36 | 3.542 |

**Do not sort or threshold on `hit.evalue`.** Hits arrive best-first, so ranking
by list order is safe; if you must threshold, use `hit.bit_score` as the E-value.
The real bit score is in raw column 13 and is dropped entirely.

The single-chain layout, 1-indexed, for anyone parsing the raw archive:
1 query chain, 2 target, 3 identity %, 4 alignment length, 5 mismatches,
6 gap openings, 7-8 query start/end, 9-10 target start/end, **11 probability**,
**12 E-value**, 13 bit score, 14-15 query/target length, 16-17 aligned
sequences, 18 target C-alpha coords, 19 target sequence, 20 taxid, 21 species.
The multimer layout matches through column 19 and then diverges — see
"`foldseek-multimer-search` returns 26 columns".

### There is no TM-score field — but `mode='tmalign'` hides it in `bit_score`

With `mode='tmalign'`, `hit.bit_score` carries the TM-score (verified: 0.9899 for
the self-hit, then 0.965, 0.9639). This is the only route to a TM-score from this
tool, and it works only *because* of the column mislabeling above.

### `target_id` is not an ID

It is `"7r0n-assembly1.cif.gz_A KRasG12C in complex with GDP and compound 2"`.
Parse `pdb_id = target_id[:4].upper()`, chain from the `_X` before the first
space, and treat the remainder as a free title.

Two things the obvious parse gets wrong. The chain token can be **multi-character
(`6T9D_CCC`)**, so do not assume one letter. And because targets are
`-assembly1` files, a repeated auth chain id is disambiguated with a **`-N`
suffix — `7AG0_A-2`, `1BTG_B-3`** — which is not a deposited chain id and will
not match `auth_asym_ids`. Split on `-` and keep the head. Missing this silently
cost five of twenty-five IL-17A neighbours their chain resolution; they fell
back to entry accessions without erroring.

### `alignment_length` counts gap columns, so the length filter is a *target* filter

It is not query coverage. Measured on the IL-17A run: hit 6YW8 has
`query_range [9, 93]` — 85 query residues — `target_range [1, 118]`, and
`alignment_length 123`. The 38 extra columns are gaps.

So `alignment_length >= 120` on a short query is mostly asking *the target* to
be long. It is a reasonable domain-level filter when the query is a ~170-residue
domain like KRAS (285 of 992 hits pass). It is close to useless on a small
protein: IL-17A's `8DYG` assembly resolves 93-95 residues per protomer, and the
strict filter left **2 hits out of 283** — a neighbourhood of exactly two NGF
structures, which would have been reported as "IL-17A has almost no fold
neighbours" when in fact it has 81.

The module relaxes automatically: when fewer than `relax_if_fewer_than` (5)
neighbours pass, it drops to `max(60, 0.7 * query_span)` and records the whole
decision — both thresholds, both counts, the query span — in
`filter.auto_relaxed`. **Never report a relaxed run without that block**, and
pass `--min-alignment-length` explicitly if you want the behaviour pinned.

### CORRECTED: `foldseek-search` DOES search every chain — it just never assembles them

**The earlier version of this section was wrong about the mechanism and it is
worth knowing why, because the wrong evidence was very convincing.**

The claim was: "`foldseek-search` takes the file, but only one chain reaches the
query", on the evidence that `8DYG` assembly 1 went in as two chains and 188
residues, and no `query_end` across all 283 hits exceeded **95** — one protomer.

That inference does not hold. The raw m8 has a query-chain column the wrapper
discards, and it shows **both chains present**: `job_A` 144 rows, `job_B` 139.
Query residues are simply numbered **per protomer** in both searches, so 95 is
what a fully-searched dimer looks like too — the multimer search on the same
file also maxes out at 95.

What is true, and what the original observation actually caught, is the
**target** side: the self-hit comes back as `8dyg_A` from both query chains and
`8dyg_B` never appears. `foldseek-search` aligns each chain independently
against each target chain independently and concatenates; it never pairs them.
There is no complex assignment and no complex TM-score, so nothing in the result
asserts that a target's chains form the same assembly as the query's.

**So the practical conclusion survives intact and the reasoning changes.** A
single-chain search on an oligomer yields a *union of per-protomer
neighbourhoods*. IL-17A's site is a groove at the **homodimer interface** and
TNF-alpha's is a cavity on the **trimer 3-fold axis**; a union of protomer
neighbourhoods is not evidence about either.

`foldseek-multimer-search` supplies the missing part, has now been run on both,
and the module routes to it automatically when the file has more than one chain.
Same input, 8DYG assembly 1:

| | `foldseek-search` | `foldseek-multimer-search` |
| --- | --- | --- |
| wire mode | `3diaa` | `complex-3diaa` |
| m8 columns | 21 | **26** |
| rows returned | 283 | **863** |
| query chains present | `job_A` 144 + `job_B` 139 | `job_A` 433 + `job_B` 430 |
| target chains of the self-hit | `8dyg_A` only | `8dyg_A` **and** `8dyg_B` |
| complex assignment | **none** | 570 groups, 293 of size 2 |
| complex TM-score | **none** | column 21 |
| max `query_end` | 95 | **95 — identical, see below** |
| distinct entries | 125 | 174 |
| entries after the relaxed filter | 81 | 137 |
| entries matching BOTH query chains | n/a | **135 of 137** |
| TM-score route | second `mode='tmalign'` search | already in the result |
| wall clock | 2.6 s – 323 s | **405 s** |

On TNF-alpha (3 chains) the difference is starker still: 1,010 rows against
**6,891**, and **1,206 complex assignments of size 3** — trimer-to-trimer
matches, which is exactly the object TNF-alpha's site belongs to.

### The complex TM-score catches a wrong biological assembly

Unexpected and useful. `2AZ5` — the canonical TNF-alpha holo entry, SPD304
bound — has **four chains in `2AZ5-assembly1.cif`**, two ligand-bound dimers,
not the TNF trimer. TNF-alpha's biological unit is a trimer and its site is on
the 3-fold axis, so that file is the wrong object to ask an interface question
about.

The multimer search says so without being told:

| query | chains | best NON-SELF complex qTM |
| --- | --- | --- |
| 1A8M assembly1 (trimer) | 3 | **0.99** (1TNF), then 0.983, 0.98 |
| 2AZ5 assembly1 | 4 | **0.495** (5MU8), then 0.495, 0.494 |

A trimer query matches the TNF superfamily at 0.85-0.99. The 4-chain query
matches nothing above 0.50, including other TNF-alpha entries — because there is
no 4-chain TNF in the PDB to match. The single-chain search cannot produce this
signal at all; every protomer matches every protomer regardless.

So **a uniformly low best complex TM-score against a family you know your target
belongs to is a signal that the assembly is wrong**, not that the fold is
unusual. Check the chain count against the known oligomeric state before
concluding anything about the neighbourhood. `assembly1` is the right default
over the asymmetric unit, but it is not a guarantee of the biological oligomer.

### `query_end` is NOT a test of which search ran

`query_end` is numbered **within a protomer** on both paths. IL-17A maxes out at
95 either way; TNF-alpha at 152 either way. The diagnostic that found the
original bug reads identically on a search that used every chain, which is how
the wrong mechanism above survived.

Test `meta['query_chains']` (raw m8 column 1: `job_A`, `job_B`, …) instead — the
distinct query chains, and for the multimer path `n_complex_assignments`. The
shared parser throws column 1 away, so this needs the raw m8; see the next
failure mode.

The upside is that the auto-relaxed length floor is unaffected: `query_span`
stays ~95, so the relaxed floor comes out the same on both paths, and a multimer
query does *not* inflate it into a starved result. That was the worry; it is
measured and it does not happen.

### `foldseek-multimer-search` returns 26 columns and the wrapper parses 12

Its own docstring says the opposite, verbatim:

```python
FoldseekMultimerHit = FoldseekHit
"""Same shape as FoldseekHit — multimer search returns the standard 12-column M8."""
```

It does not. Verified on two inputs (the shipped 1HSG fixture, 1,830 rows; IL-17A
8DYG assembly1, 863 rows): **every row has 26 columns**, and `_parse_m8_text` —
shared with the single-chain path — reads twelve fixed positions and drops the
rest. What gets dropped is precisely the multimer information:

| column | content | fate |
| --- | --- | --- |
| 1 | **query chain** (`job_A` / `job_B`) | **dropped** |
| 2-12 | target, identity, lengths, coords, prob, E-value | parsed |
| 20 | **`complexassignid`** — groups the rows of one complex match | **dropped** |
| 21 | **complex TM-score, query-normalised** | **dropped** |
| 22 | **complex TM-score, target-normalised** | **dropped** |
| 23-24 | rotation matrix, translation vector | dropped |

So `run_foldseek_multimer_search(...).hits` is a flat list of **chain-pair
rows** that is shape-identical to a single-chain result. Two rows for one entry
look exactly like one entry hit twice. **Do not consume `.hits` directly.**
`neighbour_precedent.py` re-downloads `result_url` — a static archive, already
computed, one GET, no queue time — and parses the 26 columns itself.

The grouping semantics, verified on the 1HSG fixture: 1,830 rows → **915
`complexassignid` groups, every one of size 2, none spanning two target
entries, one TM-score per group**. On 8DYG: 863 rows → 570 groups, 293 of size
2 and 277 of size 1 (a query chain that matched a protomer with no partner
assignment). Rows arrive sorted by complex TM-score, so **list order is still
the safe ranking** — the same rule as the single-chain path.

Gotcha 1 is unchanged on `/foldmulti`: column 11 is the probability, column 12
the E-value, so `hit.evalue` and `hit.bit_score` are mislabeled identically.
Confirmed on both fixtures.

### A multimer hit names several chains — the chain-attribution fix still holds

This was the thing most likely to break, and it does not. **Each m8 row still
names exactly one target chain**, so `auth_asym_ids @> to_jsonb(chain)` is
unchanged. What changes is that one entry now appears in several rows, so a
neighbour is a *set* of matched chains and the module resolves **all** of them
rather than deduplicating down to the first.

That is the correct generalisation and it does not reintroduce the VEGFR1 bug:
P17948 (VEGFR1) is **absent** from all 42 chain accessions of IL-17A's
137-entry multimer neighbourhood, exactly as it should be.

But the entry-level *ligand* bug is untouched and multimer search makes it
easier to trip over, because a complex match pulls in bigger entries. Measured
on the full 137-entry IL-17A multimer neighbourhood, 8 entries flagged
`has_druglike_holo` and **all 8 are false**:

| entry | ligand | what it actually is |
| --- | --- | --- |
| 1RV6 | `B3P` | bis-tris propane — a **buffer**; `ligand_filter` still calls it druglike |
| 4MQW | `JEF` | Jeffamine — a **crystallisation additive**; `ligand_filter` still calls it druglike |
| 4QAF | `OMA` | a cyclopropane fatty acid |
| 4EC7 | `L44` | diacylglycerol (already known) |
| 4XPJ | `LPY` | lysophospholipid (already known) |
| 7W9M / 7W9P | `9SR` / `9SL` | guanidinium channel toxins bound to **Nav1.7**, not to the matched chain |
| 8I2G | `O6F` | a **genuine 468 Da drug** — an FSHR allosteric agonist. Foldseek matched chains **X/Y = the FSH cystine-knot hormone**; the compound binds **FSHR**, the receptor. |

8I2G is the 2XAC/VEGFR1 pattern exactly, one entry later: a real small molecule,
correctly retrieved, bound to the wrong protein in the entry. The chain fix
keeps FSHR out of `accessions`, but `has_druglike_holo` is entry-level and will
still read `true` with `attribution: ambiguous_multiprotein`. **Read
`ligand_names` and the title.** `ligand_filter` removes the cofactor and lipid
rows of that table (`L44`, `LPY`, `OMA`, and the `9SR`/`9SL` lipids alongside
them) but not `B3P` or `JEF` — see "What `ligand_filter` fixed here, and what it
did NOT" below.

### What `ligand_filter` fixed here, and what it did NOT

Superseded: this section used to argue that `$EXCL` leaks and that the fix was
to read `ligand_names` by eye. The list is gone and the classifier decides. Both
of this axis's historical false positives now classify correctly, **without the
classifier having been shown either case**:

| | ligand | old verdict | `ligand_filter` verdict |
| --- | --- | --- | --- |
| KRAS → 4PHH | `2UK` | holo (635 Da, comp_id not listed) | `cofactor` — purine + ribose + phosphate |
| IL-17A → 4EC7 | `L44` | holo (625 Da, comp_id not listed) | `lipid_or_detergent` — 21-carbon chain, 48% of the molecule |
| IL-17A → 4XPJ | `LPY` | holo | `lipid_or_detergent` — phosphate head, 12-carbon chain |

IL-17A multimer accordingly went from 23 apo / 2 holo to **25 apo / 0 holo, 0
undetermined**, which is the answer the previous run argued for by hand.

**It opened a new gap, and it was a regression in one place. FIXED 2026-08-15;
the three cases below are now in the ground-truth set and all three pass.**
Recorded rather than deleted, because the two halves were fixed by two
*different* instruments and which one applies to a new case is the reusable
part.

| comp_id | what it is | MW / heavy | was | now | fixed by |
| --- | --- | --- | --- | --- | --- |
| `BEN` | **benzamidine** — protease crystallisation additive | 120 / 9 | `druglike` | `crystallisation_additive` | ubiquity prior (R14) |
| `B3P` | bis-tris propane — buffer | 282 / 19 | `druglike` | `crystallisation_additive` | chemistry (R11b) |
| `JEF` | Jeffamine — precipitant | 598 / 41 | `druglike` | `crystallisation_additive` | chemistry (R11b) |

**`B3P` and `JEF` were never a chemistry gap — they were an off-by-one.** Both
are ring-free polyols/polyethers and both missed the existing rules by a hair:
`B3P` clears the `heavy <= 18` polyol rule by ONE atom, and `JEF` fails the PEG
rule only because a Jeffamine carries a terminal amine. R11b drops the size cap
and the no-nitrogen requirement: **no ring system at all, no amide, and four or
more hydroxyl/ether oxygens** is a polyol at 19 heavy atoms and at 41.

**`BEN` is the regression, and it is the one case where a chemistry rule was
the wrong instrument.** At 120 Da the superseded `MW_MIN = 250` floor excluded
it; it was the single holo call in the IL-17A **single-chain** top 25 (2GNN, an
Orf virus VEGF variant), where the defensible count is **0 of 25**. No
structural test can reject benzamidine without rejecting a real ligand class —
its close neighbours *are* thrombin and trypsin inhibitors, and it is a bona
fide fragment. So the discriminator is not what it looks like but how it
behaves across the PDB, which is a query and not an opinion:

```sql
SELECT comp_id, count(distinct entry_id) FROM pdb_v.entry_ligands ... GROUP BY 1
```

**A frequency prior alone would NOT have worked, and the measurement says so.**
Entry counts: `BEN` 361, `B3P` 232, **`JEF` 21**, `LZ1` (a genuine fragment hit)
10, `ZBR` 9, `LFI` 8, `MOV` 7, `N5S` 1 — against `GOL` 26004 and `EDO` 17548.
**Re-confirmed against `pdb_v.entry_ligands` on 2026-08-15**; on that view
`COUNT(*)` and `COUNT(DISTINCT entry_id)` agree, so it is one row per
(entry, component). Note `tests/entry_counts.json` carries the RCSB-sourced
figures for three of these — `B3P` 235, `GOL` 26117, `EDO` 17718 — and both
sources are correct for what they measure. The harness gives identical results
under either, so name the source rather than reconciling them.
JEF is a real additive at 21 entries, below any cut that keeps LZ1 safe. The
chemistry fix and the frequency fix each catch what the other cannot.

**Why a measured list here does not contradict "a denylist cannot decide holo".**
That conclusion is about **cofactors**, and it holds: ADP/GDP/heme form an open,
growing family, `2UK` was an analog nobody had listed, and chemistry generalises
to molecules not yet made. **Crystallisation additives are not a chemical family
— they are a laboratory practice**, a small slow-growing set that changes on the
timescale of new screen kits rather than new medicinal chemistry, and one member
of it is chemically indistinguishable from a real hit. That is the narrow case
where measuring behaviour beats describing structure.

**Bounded two ways so it cannot eat real chemistry**, and the blast radius is
measured, not argued: across the whole 332-component union of the ground-truth
and held-out sets, exactly **11** components are `druglike` *and* at or below the
15-heavy-atom bar, so those 11 are the only things R14 can ever touch — `BEN`
361, `CFF` 87, `4NC` 24, `LZ1` 10, `ZBR` 9, `8VY` 2, and five singletons. At the
150-entry threshold it fires on `BEN` and nothing else; the nearest thing it does
not touch is caffeine, 4.1x below. **The threshold is UNCALIBRATED** — fitted on
one measured false positive, margin one case wide — so treat it exactly as rule
4a's volume guide is treated. Every firing is flagged `ubiquity_prior_applied`,
carries the count in `evidence` and drops confidence to `medium`.

`ChemCompSource.with_entry_counts()` runs the query and falls back to RCSB's
search API. **Use the fallback, do not remove it**: on the day this was written
Paperclip's aggregate returned nothing for ten minutes under concurrent load
while single-row selects answered in milliseconds, and a prior that silently
degrades to "not checked" whenever the cluster is busy does nothing on exactly
the runs that matter. The two sources agree where both answered (BEN 361/361,
LZ1 10/10, LFI 8/8; GOL 26004/26117, paperclip being a snapshot).

### A component covalently bonded inside a PEPTIDE ligand is not a small molecule — 8QFZ `LFI`

**This is the first measured false positive against the "0 false positives on 70
blind components" result, and it is not a gap in a list — it is the
chemistry-based method getting it wrong.** It is fixed, and how it was fixed is
the part worth carrying.

`LFI` is `C12H18Br3N3O3`, 1,3,5-tris(3-bromopropanoyl)-1,3,5-triazinane: the
**TATA tri-electrophile that cyclises Bicycle peptides**. In 8QFZ it is
covalently bonded to all three cysteines of a 12-residue polypeptide
(`CHWLENCWRGFC`, entity `8QFZ_2`) which is the actual ligand; the target is TSLP,
entity 1. Context-free it classified `druglike`, confidence `high`, reason "no
cofactor, lipid, sugar, peptide, polymer or additive signature fired", and that
chemistry read is **defensible** — the molecule really is drug-like.

**The damage was a MODALITY error, which is why it matters more than a
misclassification.** `druglike` made 8QFZ a holo small-molecule structure: the
run set `tier: holo`, `tier_note: "drug-like ligand LFI"`, auto-derived `LFI` as
the site anchor and emitted `ligand_site_jaccard` of 0.769 and 1.000 — the
strongest site-hypothesis basis this pipeline can produce — for what is
**peptide** precedent. Dossier rule 1 exists to stop exactly that substitution.

**What this is NOT.** It is not "a peptide-binding site is not a site". A groove
that binds a bicyclic peptide is a demonstrated ligandable surface and a
perfectly reasonable small-molecule target — MDM2/p53, protease substrate
grooves, peptide GPCRs. A peptide-bound groove is a **lead**, not a
disqualification. Nothing here rejects a site; it attributes a component to the
right molecule so the peptide is reported as peptide precedent instead of being
laundered into small-molecule precedent by the reagent that staples it.
`holo_call(...)["polymer_ligand_precedent"]` carries the peptide's description,
length and sequence for exactly that purpose.

**The fix: pass a `StructureContext`.** Four rules, applied on top of the
chemistry verdict, never instead of it — `classify_record` stays a pure function
of the CCD row.

| rule | test | outcome |
| --- | --- | --- |
| C0 | comp_id is in `_entity_poly_seq` for some polymer entity | `peptide_or_polymer` — it is a residue. Closes `NH2`/`ACE` |
| C1 | covalently bonded to a polymer entity that is **not** the target | `polymer_conjugate` |
| C2 | **two or more** covalent bonds to polymer residues | `polymer_conjugate`, and this needs no target identity at all |
| C3 | **exactly one** bond, to the **target** polymer | **chemistry verdict untouched** + flag `covalent_to_target` |
| C4 | exactly one bond, target identity unknown | verdict untouched, confidence → `medium`, flag |

**C2 is the rule that does not need to know who the target is, and it is the
general one.** A covalent drug carries one warhead and makes one bond; two or
more means the component is stapling a chain rather than binding it. Measured:
`LFI` 3 bonds, `ZBR` (TBMB) 3, `A1I4O` 3, `8VY` (bis(bromomethyl)benzene) 2 —
and 8VY's two bonds are to Cys427 and Cys432 of the **same chain** in 5V2P, a
crosslinked protein, which C1 could never catch because that chain may well be
the target.

**C3 is the control, and it is the reason this is not a worse bug in the other
direction.** A genuine covalent inhibitor is bonded to the target too. Verified:
6OIM's `MOV` (sotorasib, one bond to KRAS Cys12) and 4G5J's `0WN` (afatinib, one
bond to EGFR Cys797) both stay `druglike`, and `MOV` stays `druglike` even with
the accession withheld (C4). The distinction is bonded-to-a-**polymer-ligand**
versus bonded-to-the-**target**, and C1/C2/C3 encode it three different ways so
no single missing input collapses it.

**THE TRAP, and it is silent: `_struct_conn` IS NOT IN THE ASSEMBLY FILE.**
RCSB strips it exactly as it strips `_struct_ref`. Verified on 8QFZ:
`8QFZ-assembly1.cif` carries 23 categories and neither of those two, so a
context built from the coordinate file the pipeline already holds comes back
with an **empty link table that is indistinguishable from "nothing is
bonded"** — and `LFI` goes straight back to `druglike`. Build the context from
`files.rcsb.org/header/<ID>.cif`, which `pocket-scan` already fetches once per
entry for `_struct_ref`, so it costs no extra call. `StructureContext` detects
this itself: `has_struct_conn_category` is False, `is_available()` returns False,
every verdict reached that way is flagged `struct_conn_absent_from_context` and
its confidence is lowered. **Do not read the absence of links as the absence of
bonds.**

**Context-free, the chemistry still says something useful.** A component with
three or more **alkyl-halide** electrophiles is a bi/trifunctional crosslinking
reagent, not a drug. That does **not** change the verdict — the chemistry cannot
settle it — but it raises
`multi_electrophile_may_be_a_crosslinking_reagent` and drops confidence to
`medium`, so an isolated `LFI` is no longer a confident `druglike`. **The
threshold is 3 and not 2 on purpose:** nitrogen mustards — chlorambucil,
melphalan, bendamustine — are approved drugs with exactly two alkyl chlorides on
one nitrogen, and a threshold of 2 would trade this false positive for a worse
false negative. Aryl halides do not count, so `260`
(2-(bromomethyl)-1,3-difluorobenzene) scores 1, not 3.

**And it does NOT fix `MTN` (below).** MTSL makes exactly ONE covalent bond to
the target's engineered cysteine, which is rule C3 — the covalent-inhibitor
case — so it stays `druglike` by design. What context adds is the flag
`covalent_to_target`, which narrows the set of entries whose **titles** must be
read from "all of them" to "the ones with a covalent linkage". That is a
smaller job, not a solved one.

### A ligand can be disqualified by why it is there, not by what it is — TNF 5UUI `MTN`

**This is the one false positive no chemistry can fix, and it marks the boundary
of what a structural classifier can decide.**

`MTN` in TNF-alpha entry **5UUI** classifies `druglike`, and the classifier is
not wrong about the molecule: its chemistry really is drug-like. The
disqualifying fact is not in the molecule at all — it is in the **entry title**,
"Crystal Structure of Spin-Labeled T77C TNFa". `MTN` is an **MTSL spin label
covalently attached to an engineered cysteine** for EPR. It is a reagent for
measuring the protein, not a ligand bound to a site, and an entry carrying it is
**apo** for every purpose this skill serves.

Two consequences:

- **The "0 false positives on 70 blind components" figure is true of a sample
  that did not contain this case.** It is not a claim about the classifier's
  ceiling. Quote it with this caveat wherever it appears — including
  `pocket-scan/SKILL.md` and `OUTPUT_NOTES.md`, which quote the same number.
  **A second boundary was measured on 2026-08-15**: the sample also contained no
  peptide-conjugated crosslinker, and `LFI` was the first false positive
  actually observed against it. That one IS fixed, by context rather than by
  chemistry; `MTN` is not, and the two together mark the edge — a molecule's
  identity is decidable from the CCD, its *role* often is not.
### Where these figures are quoted

The accuracy figures are cited in **six places outside the harness**, and the
sole authority for all six is `tests/test_v2.py`. If you change a rule, re-run
it and update every row — a figure that disagrees with the harness is the
citation's defect, not the harness's.

| file | what it quotes |
| --- | --- |
| `structure-select/ligand_filter.py` (module docstring) | 259/262, 277/280, 61/70, 0 FP |
| `structure-select/SKILL.md` (this file, above) | 259/262, 277/280, 61/70, 0 FP |
| `pocket-scan/modal_app.py` (header comment) | 259/262, 61/70, zero FP |
| `pocket-scan/modal_app.py` (`method.ligand_classification.accuracy`) | 259/262, 61/70, zero FP — **this one ships in the output JSON** |
| `pocket-scan/SKILL.md` | 259/262, 61/70, zero FP |
| `OUTPUT_NOTES.md` | 259/262 (98.9%), 61/70 (87.1%), zero FP |

All six were **verified against a live run on 2026-08-15 and every figure is
correct as written**. Two of the six quote only the pair `259/262` + `61/70` and
omit the combined `277/280`; that is not wrong, but the combined figure is the
one that covers the crosslinker and context rules, so prefer it.

Also note `fixtures/targets.json` and `fixtures/README.md` cite the zero-FP
claim in order to *rebut* it with the `MTN` counterexample. Those are correct
and should stay.

- **Read the entry title before promoting a structure to holo.** Titles naming a
  spin label, a crosslinker, a fluorophore, a photo-affinity probe, a chaperone
  or a fusion describe why the chemistry is present, and that is information the
  molecule does not carry. `ligand_filter` decides *what a ligand is*; only the
  title decides *why it is there*. Do not add `MTN` to a denylist — the denylist
  is the defect this replaced, and the next case will be a different comp_id
  with the same shape.

### The public server's latency varies by two orders of magnitude

Same query, same day: KRAS 6OIM_A took **4.4 s**, then 8.4 s. IL-17A's 283-hit
search took **323 s** on one attempt and **2.6 s** on the next, with a
`ReadTimeoutError` on the ticket poll in between (the client retried and
recovered). Multimer is slower again: IL-17A 8DYG multimer took **405 s**. Keep
`timeout_seconds=900`, expect transient connection warnings on stderr, and do
not treat a slow run as a hung one.

### Run ONE search at a time, or you will wedge the queue

**Do not parallelise Foldseek searches.** Six were launched at once (3 multimer,
3 single-chain). Exactly one completed; the other five sat in `PENDING` and
never scheduled, failing after 15 minutes with:

```
ERROR: Tool foldseek-multimer-search: failed with TimeoutError: Timeout after 900.0s
polling https://search.foldseek.com/api/ticket/PQVvYFbyRzyvpxr_nFidnRHkFdEtXfAPHnxGjg;
last status='PENDING'
```

They were still `PENDING` **85 minutes later**, long after every client had
exited — so they were not slow, they were abandoned. A freshly submitted job
went `RUNNING` within seconds at the same moment, which is how you tell the two
apart.

Two operational consequences:

- **`PENDING` for more than a few minutes is a wedged ticket, not a slow one.**
  `RUNNING` is the state that means progress. Diagnose by submitting something
  small and fresh; if that schedules instantly, the old ticket is dead.
- **Resubmitting does not retry.** The server keys tickets on content, so the
  identical file rejoins the same wedged ticket (verified: a curl submit of
  8DYG returned the ticket a running client already held). To force a fresh
  ticket, change the bytes without changing the structure — prepend a
  `REMARK` line, which Foldseek ignores.

`--cache` writes the parsed hits to JSON so re-filtering never costs another
search. **Version any such cache.** Ours stores parsed hits, so adding the
`-N` chain-suffix strip silently invalidated every existing file — the module now
carries a `cache_version` and re-runs on mismatch.

### THE ROW CAP MOVES, AND IT IS SILENT. Reconcile every structure count against a `COUNT`.

**Read this before any other failure mode in this file.** The cap is not a fixed
200. Measured 2026-08-15 while regenerating fixture counts: **the same query
returned 200 rows one moment and exactly 10 the next** — same SQL, same source,
well-formed table, correct columns, **no error, no warning, no truncation
marker**. The first run recorded **KRAS as 10 PDB entries against a true 522**.
Only a separately issued `COUNT` caught it.

This skill is where that does the most damage, because almost everything it
produces is a *length*: `total_pdb_structures`, `holo_count`, `apo_count`, the
ensemble, the apo census behind a cryptic call, the Foldseek neighbourhood.
A capped structure list does not look broken — it looks like a target nobody has
crystallised, which is precisely the conclusion this station must never reach by
accident. It is also the leading explanation for things previously filed
elsewhere in this section: a table that looked degraded and then returned in
7 ms, and latencies spanning two orders of magnitude.

**Binding, per `CLAUDE.md` rule 14:**

- **Reconcile every count against an independently issued `SELECT COUNT(*)` /
  `COUNT(DISTINCT …)`** over the same predicate, as its own call — entry counts,
  holo/apo splits, neighbour counts, and any list whose length becomes a number.
- **A mismatch is a hard failure**: the field is `null`, and `not_found` names
  both figures and both queries. Do not report the larger one.
- **Aggregate server-side wherever the rows are not themselves needed.** A
  one-row result cannot be capped.
- **Exactly 200 or exactly 10 rows is capped until an aggregate says otherwise**
  — and 47 rows is not thereby safe. We do not know what moves the cap.
- The **per-accession loop** below is not a substitute for this. Sixteen
  per-accession queries that each came back capped still produce a wrong total,
  and each one needs its own `COUNT` before its length means anything.

### Four Paperclip failure signatures, and all four mean the query did not run

**11 of 30 SQL calls in one dry run failed**, across four distinct signatures,
three of them undocumented:

| signature | what it actually is |
| --- | --- |
| `[error] Request timed out` | seen at 120 s on a **tableless `SELECT 1`**. Not a statement-cost signal — it says nothing about your query, so do not read it as "this accession is expensive". |
| `[error] Something went wrong. Please try again.` | undocumented. No code, no detail. |
| `vsh: cd: /papers/: Permission denied` | returned by `paperclip sql` **for a SQL query** — a shell error from another subsystem, naming a path you never queried. |
| a silently capped row set | the section above. **No error text at all.** |

**Any of these means the query did not run.** The value is `null`, the reason
goes in `not_found` quoting the signature verbatim, and the retry is
short-then-long. Never `0`, never `[]`, never "no structures", never "apo". This
is the same rule already written here as `lookup_failed` /
`holo_determined: false` and "**a timeout is not a zero**" — it now covers four
signatures instead of one, and the fourth one is invisible.

**Only auth failures are guarded, and that guard catches none of the four.** The
tool layer throws on `401`/`403`/`unauthorized`/`forbidden`/`invalid api key`,
and only on a non-zero exit. Timeouts, "Something went wrong", `Permission
denied` and a capped table match none of those patterns, and the last one is not
even a failed run.

### Paperclip truncates wide cells

`sql -s proteins` truncates a long value with a literal `...`. `json_agg(t)::text` looks like a way around
the row cap and is not — a 200-row aggregate came back cut off mid-array at
~880 characters. Aggregate server-side into short columns instead
(`LEFT(title, 78)`), and rank with a window function when you need N per group.

Parse the output by the `---+---` rule's `+` positions, not by splitting on
`|` — titles contain pipes.

### Paperclip fails with rc=0, and an empty list is the disguise

**The most dangerous signature, because nothing announces it.** Three of
Paperclip's failure modes come back on stdout with a **zero exit status** and no
`---+---` rule line, so a table parser returns `[]` — and a failed query becomes
indistinguishable from "there is nothing there". On this axis that is the worst
possible confusion: the whole point is telling *no precedent exists* from *we
failed to retrieve it*.

| signature | seen |
| --- | --- |
| `vsh: cd: /papers/: Permission denied` | a **shell** error, returned for a **SQL** query, naming a path never queried. Three consecutive attempts on `3LKJ`, rc=0, 17 ms each. |
| `[error] Something went wrong. Please try again.` | undocumented; no code, no way to tell transient from permanent |
| `[error] Request timed out` | observed on queries that cost nothing; not a statement-cost signal |

Measured damage before the guard existed: a TNF-alpha run got zero rows from
`_entry_facts` with rc=0 and reported **"25 apo / 0 holo, 0 undetermined, 0
accessions, 0 rejected ligands"** — a clean, confident, entirely fabricated
negative. With the guard the same input reports **"25 apo / 0 holo,
undetermined: 25"** plus a `not_found` entry naming all 25 entries.

`_run_sql` now raises on all three signatures, and separately on any rc=0 output
with no table and no `(N rows)` trailer. A genuine empty result set — which does
have the rule line and the trailer — still returns `[]`. **Never let a bare
`len(rows) == 0` mean "none exist".**

`_run_sql_retry` wraps it with short-then-long budgets and returns **`None` when
the query never ran**, which is deliberately a different value from `[]`. Every
call site branches on that distinction:

| call site | on `None` |
| --- | --- |
| `_entry_facts` | returns `{}`; the caller's missing-entry reconciliation then marks every neighbour `holo_determined: false` and names them in `not_found` |
| `_druglike_comp_ids_for_accessions` | that accession goes in `failed_accs` → `lookup_failed: true` on its block |
| `_accession_precedent` summary | `n_structures: null`, block marked undetermined |

Verified live while writing this: `_run_sql_retry` on `3LKJ`'s ligand row
returned `None` and printed `QUERY NEVER RAN` — not `[]`, and not an empty
ligand list that would have made a holo entry look apo.

Three more query shapes that fail, all measured on `pdb_v.entry_ligands`:

- **`comp_type` in a WHERE clause times out.** 6 ms without it, >120 s with it,
  same 25 entries. See the `comp_type` section.
- **`IN (SELECT ...)` is the fast plan; a direct JOIN is not.** For P15692
  (VEGF-A, 75+ entries), `WHERE l.entry_id IN (SELECT entry_id FROM s)` returns
  in **9 ms**, while `FROM structures_by_accession st JOIN entry_ligands l ON
  l.entry_id = st.entry_id WHERE st.accession = ...` **times out**. Same rows,
  same accession.
- **One accession at a time — the brief's rule, and it is not optional.**
  Sweeping the 17 accessions of the IL-17A single-chain neighbourhood as one
  `unnest` array timed out; per accession, 16 of 17 returned in ~1.4 s and
  **P67861 timed out reproducibly at 120 s**. The same is true of the holo
  aggregation: unioned over the 8 TNF-superfamily accessions it timed out even
  at a 300 s client budget, and per accession it is bounded. Both sweeps now
  loop. On persistent failure the accession is marked `lookup_failed` /
  `holo_determined: false` rather than recorded as having no drug-like ligands.
  **A timeout is not a zero.**
- **Retry short-then-long, not long-twice.** A deterministically slow statement
  retried at the same generous budget just doubles the bill — P67861 at
  2 x 300 s would cost ten minutes of a run with sixteen good accessions left.
  90 s then 240 s.

**Attribution caveat on all of the above.** Some of these timeouts were measured
while a parallel session was running **16 concurrent `paperclip` processes**
against the same backend, so the absolute thresholds are load-dependent and
should not be read as fixed properties of the queries. What is NOT
load-dependent, because it was measured back-to-back in the same seconds, are
the *relative* results: `comp_type` present vs absent (6 ms vs timeout) and
`IN (SELECT ...)` vs JOIN (9 ms vs timeout). Trust the comparisons, treat the
timeout values as a floor on how bad it gets, and keep the retries — a shared
backend under someone else's load is the normal condition here, not an
anomaly.

`neighbour_precedent.py` also short-circuits the accession aggregation entirely
when nothing classified drug-like: the `n_dl > 0` filter is then unsatisfiable,
and running it anyway costs 120+ s to prove an already-known empty result.

### Untested, do not assume

- **Foldseek local mode** — no binary present; remote is the only tested path.
  `search_mode='local'` also requires `local_db`, which we do not have.
- **`mode='complex-tmalign'`** — never run. The multimer path takes its
  TM-scores from column 21 of the `complex-3diaa` result instead, so there has
  been no reason to try it.
- Databases other than `pdb100`.
- **The relaxed alignment-length floor** (`0.7 * query_span`, min 60) is a
  judgement call, not a calibrated threshold. It was chosen so IL-17A returned a
  neighbourhood at all. It has one sanity check — the 81 hits it admits are the
  cystine-knot superfamily, which is the correct fold answer for IL-17A — and no
  other validation. The strict 120 floor is the verified one. It has now been
  checked on the multimer path and behaves identically (`query_span` is
  per-protomer, so the floor does not move), but it is still uncalibrated.
- **The multimer path's ranking has one target's worth of evidence.** Rows come
  back sorted by complex TM-score and that ordering looked right on IL-17A
  (IL-25 first, then the neurotrophins). It has not been checked against a case
  where the correct answer is known to be a low-TM complex.
- **`complexassignid` groups of size 1** — 277 of 8DYG's 570 groups. Read as "a
  query chain matched a target chain with no partner correspondence", but the
  server's assignment rules were not confirmed against Foldseek's source. The
  module reports `n_query_chains_matched` per neighbour so a size-1 match is
  visible rather than assumed away.

## Output

Fill the dossier's `structure` block: chosen tier, entry ID, resolution,
biological unit, the ligand with its heavy-atom count and whether it is
drug-like, total and holo counts, and the ensemble actually used. Fill
`structural_neighbour_precedent` from step 5. Record the as-of cutoff if one was
applied and how many entries it removed.

For `structural_neighbour_precedent` specifically, carry through: **the
`search_path`** (`multimer` or `single_chain`) and `chains_assembled_into_complexes`, the
TM-score with its `tm_score_kind` (a complex TM-score and a chain TM-score are
not the same number — never the raw `evalue`/`bit_score`), **both** holo counts
per neighbour with the entry-level one named as an upper bound, the ligand
*names* and `rejected_ligands` not just comp_ids, `n_undetermined` alongside any
zero holo count, and `filter.auto_relaxed` whenever it fired.

If the query was an oligomer and `search_path` is `single_chain`, say so: the
neighbourhood is then a union of per-protomer matches with no complex
assignment, and it is not evidence about an interface site.

A neighbourhood with no small-molecule precedent is a real and reportable
finding — it was the answer on both calibration targets — not a retrieval
failure. Do not pad it.
