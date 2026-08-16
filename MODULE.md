# Simulation module — druggability-dossier station

This station exposes the standardized runnable-module interface. It produces one
druggability dossier for one protein target and writes it as a single JSON file
that validates against `schemas/output.schema.json`, with a side-panel
interpretability object at `output.interpretability`.

## The one command

```bash
python -m simulation run --input input.json --output output.json
```

Copy-pasteable example against the shipped example request:

```bash
# from the station root (managed/druggability-dossier/), with the env from below active
micromamba run -n druggability-simulation \
  python -m simulation run --input examples/input.json --output /tmp/dossier.json
```

`--input` and `--output` are both required paths. `--input` must satisfy
`schemas/input.schema.json` (one required field, `uniprot_accession`).

## Exit codes

| code | meaning |
| --- | --- |
| `0` | success — dossier written, and BOTH the output schema and the interpretability schema validated |
| `2` | invalid input — request missing/failing `schemas/input.schema.json`; **nothing is written** |
| `3` | dossier production failed — the managed agent could not be invoked (see below); **nothing is written** |
| `4` | validation failed — the dossier was written to `--output` for inspection, but it (or its interpretability object) did not validate |
| `1` | usage error |

Any nonzero code is a failure. Exit `0` is returned only when the dossier and its
interpretability object both validate.

## Where results and logs go

- **All machine-readable results go to the `--output` file.** That JSON is the
  only data artifact; nothing structured is printed to stdout.
- **STDOUT stays empty.** A caller can safely capture stdout and expect nothing.
- **All human-facing progress and errors go to STDERR** (Python `logging`,
  configured to stderr). Schema violations are printed there, one per line.

## Setup — locked environment / one command

The runner needs only Python + `jsonschema`. The heavy scientific stages
(fpocket, mdpocket, gemmi, GPU folding) run remotely in the Claude Managed
Agent's sandbox, so none of them are needed locally. Dependencies are pinned and
shared with no other repo.

Either of these single commands (pick one):

```bash
# a) create the locked micromamba environment:
micromamba env create -f simulation/environment.yml      # env name: druggability-simulation

# b) or install the pinned runner deps into an existing env:
micromamba run -n druggability pip install -r simulation/requirements.txt
```

Both `simulation/requirements.txt` (exact versions) and
`simulation/environment.yml` are provided; versions match what is installed in the
`druggability` env (`jsonschema==4.26.0`).

## The interpretability object

- **Lives at `output.interpretability`** in the dossier the module writes — a
  consumer reads exactly one file.
- **Schema:** `schemas/interpretability.schema.json`. The full display heuristic
  is documented there and in `INTERPRETABILITY_PANEL.md`.
- **Display heuristic, in one paragraph:** it is a WORKFLOW reasoning trace (which
  pipeline stages ran and what each computed), never LLM chain-of-thought. The
  panel renders, top to bottom: a HEADLINE (verdict chip + one honest sentence +
  which axis carried it); an optional FIGURE beside it; a BANNER for an
  axis-conflict or predicted-structure warning; the TWO AXES — retrieved precedent
  and computed tractability — as side-by-side columns that are **never merged into
  one score**; a collapsible per-stage TRACE; a muted "Not retrieved" list that
  makes null-is-not-zero visible; and a provenance footer. Every value is pulled
  verbatim from a dossier field; a required caveat always travels with its number
  (e.g. druggability is a within-structure rank, never compared across
  structures). `simulation/interpretability.py` builds it as a pure function of
  the dossier.

## How the dossier is produced (and what a live run needs)

`run_pipeline` in `simulation/pipeline.py` invokes the **Claude Managed Agent**
that is this station — there is no pure-Python re-run of the science. It drives
the documented headless route (`bun run console druggability-dossier -- --once
"<task>"`, README step 3 / `scripts/console.ts`), which runs the deployed agent,
answers any custom-tool round-trips in-process, prints the agent's final reply
(the dossier JSON) to stdout, and logs to stderr; the module then parses that
JSON. A live run therefore needs `ANTHROPIC_API_KEY`, network access, `bun` on
PATH, and the agent to have been deployed (`manifest.deployment.agent_id` set).
When any of those is absent, `run_pipeline` raises a typed error and the command
fails loudly (exit `3`) rather than hanging or fabricating a dossier.

The **schemas, examples, and interpretability logic are self-contained** and need
none of that: they run offline against the checked-in `examples/` and fixtures,
which is what `simulation/test_module.py` exercises.

## Tests

```bash
micromamba run -n druggability python simulation/test_module.py
```

Offline, no paid calls: it monkeypatches `run_pipeline` with a recorded real
dossier (`examples/output.json`) to check the end-to-end contract (exit 0, output
written, stdout empty, output validates), checks that malformed input and a
raising pipeline both fail loudly with nothing written, and checks that
`build_interpretability` validates against `schemas/interpretability.schema.json`
for the example dossier and both integration fixtures.
