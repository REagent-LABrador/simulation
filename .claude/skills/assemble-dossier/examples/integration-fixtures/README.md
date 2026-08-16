# Integration fixtures — not calibration examples

The two files in this directory are **inputs for station 4's end-to-end run**, not
reference artifacts for the validator or the rubric.

```
irak4_Q9NWZ3_orthosteric.json                 verdict small_molecule_tractable
irak4_Q9NWZ3_oligomer_destabilisation.json    verdict insufficient_evidence
```

## Why they exist

Station 4 (`therapeutic-program-economics`) checks a dossier's subject against the
programme's target. With only `jak1_P23458.json` and `tnf_P01375.json` committed,
an IRAK4 programme quarantines **every** tractability row as `SUBJECT_MISMATCH`
and the station is a no-op in the end-to-end run. One IRAK4 dossier makes it
load-bearing.

Two, not one, because **the pair is the point**. Same accession, two
`mechanism_hypothesis` values, two different verdicts, and the two runs
**share zero residues** — the kinase-site pocket lines Q9NWZ3 residues 192–329,
the Myddosome death-domain cavity lines residues 11–100. A single mode would
quietly hide the input-contract problem the integration test exists to surface:
`ProgramInput.target` is free text with no accession and no mechanism field, so
both programmes collapse to the same input and can produce the same
`input_digest`, and `insufficient_evidence` has no landing place in the station's
schema at all.

## What they are NOT

**Do not calibrate the validator or the rubric against these files.** The two
files in the parent directory are the calibration examples; these are not. A
calibration example that goes stale makes the *grader* wrong, not just the file,
and these are a snapshot of a pipeline that changed repeatedly on the day they
were produced.

Each file carries a `_fixture` key saying the same thing, and a
`_pipeline_state` key recording which rules were in force when it was produced.
If this directory and a fresh run disagree, **the fresh run is right**. Re-run
rather than trust these; they are here so an integration test has something to
read, not to settle a question about IRAK4.

## Provenance

Produced 2026-08-15 by a full pipeline run:

- `pocket_scan` (Modal, fpocket 4.2.3 + P2Rank 2.5.1, D swept over {1.6, 2.4}),
  chains asserted per rule 2b — `6EGE=A;2OIB=A;2O8Y=A` for the orthosteric run,
  all 14 chains of `3MOP` for the oligomer-destabilisation run
- Paperclip `proteins` (`uniprot_v`, `pdb_v`, `chembl_v`, raw `chembl.*`), every
  count reconciled against an independently issued aggregate
- `precedent-lookup/modality.py` for compound modality from structure
- `neighbour_precedent.py` (Foldseek multimer, pdb100) for the structural-neighbour axis

Both validate at **0 violations** against
`../../validate_dossier.py`.
