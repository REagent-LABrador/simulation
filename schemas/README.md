# Simulation-station JSON Schemas

The machine-readable contract for the simulation (druggability-dossier) station.
Three files, all [JSON Schema draft 2020-12](https://json-schema.org/):

| file | what it governs |
| --- | --- |
| `input.schema.json` | the **request** a caller sends the station |
| `output.schema.json` | the **dossier** the station returns |
| `interpretability.schema.json` | the **`output.interpretability`** object — the LABrador shared interpretability contract (v1.0.0) |

## The interpretability contract (required)

Every successful or scientifically degraded/abstaining run carries a top-level
`interpretability` object, and `output.schema.json` **requires** it (an
infrastructure failure that produces no result file is the orchestrator's
concern). It is the **shared LABrador contract**, identical in shape across
modules: `schema_version`, `headline{title,result,plain_language,status,basis[]}`,
`metrics[]`, `steps[]`, `evidence[]`, `assumptions[]`,
`uncertainty{method,intervals[],seed,draws,limitations[]}`, `limitations[]`,
`counterfactuals[]`, `lineage[]`, `extensions{}`. It **supersedes** this station's
earlier `module/axes/trace/figure` object; that content now lives under
`extensions` (the two evidence axes kept separate — never averaged — plus the
ordered trace, identifiers, figure reference, and an input-hash `cache_key`).
`../simulation/interpretability.py` builds it as a pure, deterministic map from the
dossier — no recomputation, no fabrication; unknowns stay `null` and earn a
`limitation`. Enums: `status` SUPPORTED|QUALIFIED|INCONCLUSIVE|FAILED|NOT_APPLICABLE;
`basis` OBSERVED|INFERRED|MODELED|SYNTHETIC; `direction` positive|negative|neutral|mixed|unknown;
`grade` HIGH|MODERATE|LOW|UNSUPPORTED; `severity` INFO|WARNING|ERROR.

## What these are, and what they are not

These schemas enforce **shape, types, and vocabularies** — the structural
contract two teams can code against. They are deliberately *not* the semantic
gate. The semantic gate is `../.claude/skills/assemble-dossier/validate_dossier.py`,
which enforces the cross-field rules a JSON Schema cannot express (a verdict must
be supported by the axis it cites, a druggability range may not be pooled across
structures, a failed retrieval must be `null`-with-reason and never a `0`, …).

Read the split this way:

- **input.schema.json** is closed and strict. Exactly one field is required —
  `uniprot_accession`, the resolved accession the whole dossier is keyed on — and
  no unknown fields are allowed, because a typo in a request field is a caller
  bug worth catching at the door.
- **output.schema.json** is strict at the **top level** (all 17 keys present,
  every enum a known value, correct types) but **open inside** each block. Real
  dossiers legitimately extend the minimal template with provenance and
  annotation keys (`source`, `mechanism`, per-field notes), so nested objects
  accept extra properties. Enforcing more than a validator-passing dossier
  actually contains would make the schema reject correct work — the opposite of
  useful.

## The three authorities, kept in sync

The output schema is derived from, and cross-checked against, two other
artifacts. `test_schema.py` fails if they drift:

1. **`../CLAUDE.md` → "## Output template"** — the canonical shape an agent fills.
2. **`../.claude/skills/assemble-dossier/validate_dossier.py`** — the authoritative
   enum vocabularies (`ENUMS`, `VERDICTS`, `VERDICT_BASES`, `INTERFACE_CLASSES`)
   and the 17 `REQUIRED_TOP_LEVEL` keys. Where the template abbreviates an enum
   (it lists four interface classes; the validator accepts seven), **the
   validator wins** and this schema matches the validator.
3. **`../.claude/skills/assemble-dossier/examples/*.json`** — real, validator-passing
   dossiers. Every one must validate against `output.schema.json`.

## Enum vocabularies (authoritative here, mirrored from the validator)

- `verdict`: `small_molecule_tractable` · `not_tractable` · `insufficient_evidence`
- `verdict_basis`: `retrieved_precedent` · `computed_tractability` · `both` · `none`
- `input.mechanism_hypothesis`: `orthosteric` · `allosteric` · `oligomer_destabilisation` · `unknown`
- `structure.tier`: `holo_experimental` · `apo_experimental` · `cofolded` · `predicted` · `sampled_ensemble` · `none`
- `tractability.cryptic_pocket_risk`: `low` · `medium` · `high` · `undetermined`
- `tractability.cryptic_mechanism`: `loop_or_backbone_motion` · `sidechain_occlusion` · `subunit_occlusion` · `none` · `undetermined`
- `tractability.cryptic_potency_prior.expected_ceiling`: `nanomolar` · `micromolar_at_best` · `unknown`
- `tractability.pocket_vs_interface.classification`: `orthosteric_candidate` · `allosteric_candidate` · `destabiliser_candidate` · `no_partner_structure` · `mixed` · `no_pocket_to_classify` · `numbering_mismatch_not_interpretable`
- `tractability.mdpocket_site_definition_used`: `site_from_ligand` · `site_from_density` · `none`
- `site_pocket_selected_by` (a single value **or a list**): `ligand_site_jaccard` · `site_signature_overlap` · `site_signature_unreliable_homooligomer` · `max_druggability_no_ligand_site` · `no_pocket_matched_site_signature` · `no_pocket_overlapped_ligand_site` · `site_signature_unreliable_foreign_polymer`. Only `ligand_site_jaccard` identifies a site without qualification; the last five do **not** identify a site, and their values must not be pooled as one.

Every enum field is **present and non-null** in a valid dossier. Each enum
carries an explicit unknown-member (`none` / `undetermined` / `unknown`) that is
used instead of `null` — the validator flags a null enum as a defect.

## For the pipeline assembler: the follow-up-question channel

`output.schema.json` carries an **optional** top-level `follow_up_questions`
array. It is **not required**, and no producer or consumer has to support it yet
— that is deliberate, so this schema does not force the other stations' schemas
to take it on.

It exists because the capability is already real: the `graph-intake` stage can
emit structured asks back to the upstream knowledge graph — verbs
`expand_node` · `resolve_link` · `test_gap`, each with a target node/link id and
a human-readable question (`ppi-hypothesis/fixtures/worked_ask.json` is a worked
example). Today those asks are not surfaced in the dossier. `follow_up_questions`
is the reserved, schema-blessed place to surface them **when the orchestrator is
ready to route them back upstream and re-run with more evidence**. Assembler:
know the channel is here and wire it when you want it; ignore it until then.

## Running the tests

```bash
/Users/bb/.local/bin/micromamba run -n druggability python schema/test_schema.py
```

Requires `jsonschema`. Offline; reads only files in this repo.
