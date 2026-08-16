# simulation

A **druggability dossier** station: given a protein target, it assembles evidence
on whether that target can be drugged **with a small molecule**. It reports
evidence; it does not decide. It is one specialist station in a larger evidence
gauntlet that scores asset-to-indication hypotheses — other stations handle
genetics, expression, perturbation, PK/PD, safety, and clinical precedent. This
one handles small-molecule tractability, and nothing else.

## What it does

Given a `uniprot_accession` (plus optional `as_of_date`, `disease_context`,
`interaction_to_disrupt`, and `mechanism_hypothesis`), the station produces a
single JSON dossier that evaluates the target along **two independent axes** and
reports them separately.

## The two-axis architecture

The core design principle: the two axes answer different questions, they are
allowed to disagree, and **they are never averaged into a single overall
score**. There is no combined number. When the axes disagree, the disagreement
is recorded and explained rather than resolved.

- **Axis 1 — retrieved precedent.** What has actually been made against this
  target: measured bioactivity, approved drugs, patents, terminated programs.
  This is *looked up*, not computed. It is the stronger axis when it exists.
  Only genuine small-molecule entries count toward small-molecule precedent;
  biologics and peptides are tracked separately so a reader can see when a target
  is *validated* but not *small-molecule tractable*.

- **Axis 2 — computed tractability.** What the structure says about whether a
  small molecule could bind: pocket geometry, disorder, affinity prediction.
  This is *computed*, and it has known blind spots that the dossier declares.

Every dossier records which axis carried the verdict (`retrieved_precedent`,
`computed_tractability`, `both`, or `none`) so that a single label over two
axes is never mistaken for an average.

## The skill pipeline

The station runs a sequence of skills, each with a narrow contract:

1. **graph-intake** — parse the incoming request and upstream hypothesis graph.
2. **structure-select** — pick the structures a pocket scan should run on;
   classify holo vs apo by actual ligand chemistry rather than by label,
   assemble an ensemble, apply the as-of date cutoff, and find structural
   neighbours to establish structural-homolog precedent.
3. **pocket-scan** — detect and measure ligand-binding pockets across the
   structure ensemble, sweep clustering, and quantify whether a site is cryptic.
   Reports volume with its spread as the primary number and druggability only as
   a range.
4. **precedent-lookup** — retrieve measured bioactivities, approved and clinical
   drugs split by modality, structures, and family activity, joined on UniProt
   accession. Fills the retrieved-precedent axis.
5. **assemble-dossier** — combine the two axes into the final JSON dossier and
   validate it.

Additional skills support this core pipeline, including `falsification-sweep`
(attaching disconfirming evidence for the reader to weigh), `terminated-programs`,
`ppi-hypothesis`, and `cofold-check`.

## Input / output contract

The input and output are JSON. The formal JSON Schemas live under `schema/`:

- `schema/input.schema.json` — the request contract.
- `schema/output.schema.json` — the dossier contract.

The `input` block is echoed back verbatim on every run and is never inferred.
See `CLAUDE.md` for the full field-by-field contract and operating rules.

## Not a substitute for experiment

This station reports computational and retrieved evidence about small-molecule
tractability. It does **not** decide whether to pursue an indication, does not
rank hypotheses, does not design molecules, and does not assess biologics. A
computed tractability signal is a prediction with declared blind spots, not a
measurement. Nothing here substitutes for experimental validation.
