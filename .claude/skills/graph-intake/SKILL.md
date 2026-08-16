---
name: graph-intake
description: >
  Reads an upstream literature evidence graph and extracts the dossier's input
  contract from it — which things are protein targets, what the molecule is
  meant to stop, and which UniProt accession that is. Emits follow-up asks back
  to the graph using the graph's own four ask verbs — at intake, and later in
  the build when a falsification finding turns out to be a literature-provenance
  question. It does NOT assess tractability, does NOT rank targets against each
  other, does NOT decide a mechanism the evidence has not stated, does NOT treat
  a computed result as literature evidence, and does NOT ask upstream anything a
  local table, a registry, a structure or arithmetic can settle. A pending ask
  never blocks a verdict and never justifies a null.
---

# graph-intake

The dossier's required input is a UniProt accession. Two optional inputs —
`interaction_to_disrupt` and `mechanism_hypothesis` — decide which chains get
scored, and therefore change the druggability number (rule 2b: KRAS 4OBE is
0.442 at rank 1 on chain A and 0.257 at rank 6 on chains A+B). This skill gets
all three out of an upstream evidence graph, or asks the graph for what is
missing.

It reads the graph. It does not ask the graph to change format.

## Setup

`PAPERCLIP_API_KEY` in the environment, for accession resolution:

```bash
set -a; . <repo>/.env; set +a
```

The upstream graph is a JSON file with `things`, `links`, `findings`, `papers`
and `gaps`. Its `rounds` array records the asks already issued against it.

## What the real graph actually looks like

The producer is `managed/research-evidence-mapper/` on the LABrador repo's
`labrador/main` branch. Read it without merging:

```bash
git ls-tree -r labrador/main managed/research-evidence-mapper/
git show labrador/main:managed/research-evidence-mapper/SCHEMA.md
git show labrador/main:managed/research-evidence-mapper/runs/g_1a4f.json
```

One real run exists: `g_1a4f`, round 2, `status: ok`, 5 things / 6 links /
12 findings / 4 papers / 3 gaps. Diffed against what `graph_read.py` expects:

**Fields we read that behave as documented** — `schema_version` `1.1`,
`graph_id`, `question`, `round`, `status`, `coverage.{depth,truncated,
stop_reason,no_quote_discarded,limits,found,read,used}`, `things.{id,name,kind,
aliases,mentions}`, `links.{id,from,how,to,yes,no,no_effect,state,why,basis,
confidence,changed_in_round}`, `papers.{id,title,year,journal,doi,first_author,
study_type,is_preprint,retracted}`, `gaps.{id,missing,implied_by,note,
confidence,searched_in_round}`, `rounds[]`. No breaking rename anywhere.

**Additive fields, none of which we read** — `papers[].pmid`;
`coverage.{queries_run, figures_read, figures_skipped_reason}`;
`findings[].claim`, a model-written paraphrase sitting beside the verbatim
`quote`. Harmless, but `claim` is *not* quote-guaranteed and must never be used
where a quote is required.

**Fields SCHEMA.md defines that the real run omits** — `findings[].flags` and
`findings[].round` are simply absent. `flags` matters: `is_secondary()` reads
`flags: ["background"]` as one of its three signals, so on this graph that
signal is permanently silent and only `is_own_result` and `study_type` are
doing the work.

**Guarantees that are documented and not enforced** — see failure mode 13.
`rounds[1].target` is `"g3"`, and `gaps` holds `g1`, `g2`, `g4`. Finding `f6`
is referenced by no link.

**The difference that actually matters is not in the field list.** It is that
the real graph contains no `protein` or `gene` node at all — see failure mode
12b. The format is fine. The granularity is not.

## Procedure

### 0. Check `status` before anything else

`status` is `ok | empty | partial | failed`, and SCHEMA.md note 7 is explicit
that failure is still a graph — the lists are simply empty and `coverage` is
still real. So a `failed` graph parses cleanly and returns zero nominations,
which reads exactly like "no targets in this literature."

Never report zero nominations without quoting `status` alongside it.

### 1. Traverse — mechanical, use the helper

```bash
python3 graph_read.py <graph.json>
```

Stdlib only, no dependencies. It returns seven blocks and **all seven are
load-bearing**:

| block | why you read it |
| --- | --- |
| `integrity` | dangling ids, duplicate ids, findings no link references. Any of these means a row was silently lost or overwritten. See failure mode 13. |
| `coverage` | whether an absence in this graph may be read as an absence at all. `warnings` empty is the only clean state. |
| `selection` | how many candidates, and the explicit statement that the graph offers no basis to choose between them. |
| `nominations` | the candidates, each with `evidence_class`, `support`, `quote_reading` and `direction_warning`. |
| `contested_nominations` | things that are **both** nominated and rejected. Almost always the gap door — failure mode 12. |
| `rejected` | where a wrongly-dropped target would be hiding. |
| `needs_adjudication` | a decision waiting for you: an unknown verb, or a verb whose own quotes contradict it. |

Neither `rejected` nor `needs_adjudication` is a log to skim.

**Zero nominations is never a result on its own.** It is what a granularity
mismatch, a dangling id, an unknown verb and a genuinely empty literature all
look like. Read `selection.note`, `status`, `coverage` and `integrity` before
saying the graph contains no targets. Failure mode 12b is what this looks like
on real data.

The helper now **refuses** rather than degrades on: an absent or null top-level
list, a non-object `coverage`, a `NaN`, and a `--thing` id that was never
nominated. Each of those used to produce a clean exit 0 with an empty or wrong
answer, which reads as "no targets in this literature."

Nomination rule, both halves required:

- **`kind` is `protein` or `gene`**, and
- either the object of a direct-action edge from a `small_molecule` (`inhibits`,
  `binds`, `blocks`, `degrades`, …), **or** named in a `gaps[].missing` pair.

The gap half is what carries undrugged candidates. Without it the intake can
only ever return targets somebody has already made a molecule against, which
inverts the point of the pipeline.

Both kinds, because the same target is typed either way depending on how a
paper phrased it — "IRAK4 knockdown" reads as a gene, "IRAK4 kinase activity"
reads as a protein. Keying on `protein` alone drops half of them into
`rejected`, where nobody looks.

The helper classifies an edge by the **subject's kind**, not by the verb alone.
`activates` from a small molecule is an agonist; `activates` from a receptor is
pathway biology.

**`how` has no enum in SCHEMA.md** — every other categorical field there carries
an explicit `a|b|c` comment and `how` does not. But the upstream `claim-extraction`
SKILL.md **does** pin it, to a closed list of eleven:

> `binds` · `inhibits` · `activates` · `increases` · `decreases` · `causes` ·
> `treats` · `associated_with` · `expressed_in` · `measured_by` ·
> `no_relationship_stated`

That list barely overlaps ours. Three of the eleven are in `DIRECT_ACTION`
(`binds`, `inhibits`, `activates`), exactly one is in `DOWNSTREAM_EFFECT`
(`increases`), and the remaining seven — including `decreases`, the single most
likely verb for a cytokine readout — match neither set.

And the extractor does not obey its own list. The one real graph we have
(`g_1a4f`) emits `blocks`, `suppresses` and `drives`, none of which appear in the
eleven. So the vocabulary is unaligned in both directions: our sets do not match
their spec, and their output does not match their spec either.

Two consequences. `needs_adjudication` is the **normal** path, not a corner case —
against the documented vocabulary it fires on 7 of 11 verbs. And a verb landing
in one of our sets proves nothing about whether it was the right set, which is
why the verb no longer decides alone (failure mode 12c).

### 2. Read the mechanism out of the quotes — judgment, yours

`how` is too coarse. `inhibits` does not say what is inhibited. The mechanism is
in the `quote` text, and it is why the helper hands you quotes rather than a
verdict.

`interaction_to_disrupt` accepts exactly three shapes. Match the quote to one:

| shape | quote looks like | example |
| --- | --- | --- |
| catalytic function | "…inhibited IRAK4 **kinase activity**…" | kinase, protease, ATPase |
| named partner | "…**nucleates assembly** of the MyD88 signalosome…" | a PPI |
| oligomeric state | "…displaces a subunit of the **trimer**…" | TNF-alpha |

A target can have more than one. Report every shape the evidence supports —
see failure mode 6.

### 3. Apply the three-state basis rule

This is the rule that keeps a weak answer from driving a hard output. The tier
comes from the **link's `basis`**, never from the finding's own confidence. A
0.88-confidence quote from a review is still background.

| tier | what to do |
| --- | --- |
| `primary` | usable. Set `interaction_to_disrupt`. |
| `mixed` | usable, but carry the disagreement into `notes`. |
| `hedged_only` | **record it, do not act on it.** Every supporting finding said "may" or "suggests". Same treatment as background. |
| `background_only` | **record it, do not act on it.** Treat as unstated for chain selection, and issue a `resolve_link` ask. |
| absent | unstated. Issue a `new_question` ask. |

`background_only` is more dangerous than absence. Absence trips the dossier's
existing refusal to guess. A review citation does not — it produces a confident
looking answer from a single secondary source.

### 4. Resolve the accession — from the quote, not the name

```sql
SELECT accession, gene_name, protein_name, organism, sequence_length
FROM uniprot_v.proteins
WHERE gene_name IN ('<SYM>', …) AND organism = 'Homo sapiens'
```

Resolve using what the **quotes** say, not the node's `name` string. See failure
mode 3 — this is the one that silently assesses the wrong protein.

If two accessions both fit, populate `ambiguity` with both and leave
`uniprot_accession` null. An unresolved target is a correct output. A confidently
wrong accession poisons every number downstream of it.

### 5. Emit

One object per nomination, matching the dossier's Contract table:

```json
{
  "uniprot_accession": "Q9NWZ3",
  "gene_symbol": "IRAK4",
  "disease_context": "rheumatoid arthritis",
  "interaction_to_disrupt": "kinase activity (catalytic function)",
  "mechanism_hypothesis": "unknown",
  "provenance": {
    "graph_id": "g_ra4k", "round": 2, "thing": "t1",
    "primary_findings": ["f1", "f5"],
    "recorded_not_acted_on": ["f4", "f3"]
  },
  "asks": []
}
```

Leave `mechanism_hypothesis` as `unknown` unless a quote states the site. See
failure mode 5.

## Asking back

Use the graph's own verbs. SCHEMA.md defines **four**: `expand_node`,
`resolve_link`, `test_gap`, `new_question`. Point at a row by id, never in prose.
One ask per request, one round per request.

| situation | ask |
| --- | --- |
| we have a target and want everything touching it | `{"ask": "expand_node", "target": "<thing id>", "depth": "deep"}` |
| mechanism only in `background_only` | `{"ask": "resolve_link", "target": "<link id>", "depth": "deep"}` |
| mechanism absent entirely | `{"ask": "new_question", "target": null, "depth": "deep", "question": "<specific>"}` |
| a gap names our target | `{"ask": "test_gap", "target": "<gap id>", "depth": "deep"}` |

Before reporting anything as "not stated", check `coverage.truncated` and
`coverage.stop_reason`. Of the five stop reasons — `max_papers`,
`queries_exhausted`, `no_new_results`, `time_limit`, `complete` — **only
`complete` means the literature was exhausted.** The other four mean the run
ran out of budget, so absence proves nothing. The helper surfaces this as
`coverage_warning`.

At `depth: "quick"` absence never means anything at all: SCHEMA.md note 2 says
`quick` reads page 1, and page 1 lies.

## Asking back after intake

Intake is not the only moment a literature question arises. A falsification
sweep three steps later can turn up a claim that changes the dossier and rests
on nothing but review text — and until now the pipeline either burned its own
effort on it, reported it unresolved, or dropped it. It can now issue the same
four asks, mid-build, under a much narrower rule than the one intake uses.

**Same four verbs, same id-pointing rule, same one-ask-one-round rule. No fifth
ask type and no new schema.** The upstream team consumes what it already
consumes. The only field reused beyond intake's usage is `question`, which
SCHEMA.md already defines on `new_question`; populate it on every ask type,
because an ask that is only a verb and an id is unanswerable (see "Make the ask
carry its own evidence" below). That reuse is the entire schema footprint. If
you find yourself wanting a fifth verb, the thing you want to ask is almost
certainly on the never-ask list.

### The trigger — all five gates, no exceptions

A finding becomes an ask only when **every one** of these holds. Four is not
enough; the gates are conjunctive because each one alone is satisfied by large
numbers of claims we must handle ourselves.

| # | gate | fails when |
| --- | --- | --- |
| 1 | **It changes a dossier value.** If true, the claim fills or alters a field in the template — a precedent list entry, a count, a mechanism, a structure choice. | it is an efficacy, outcome or clinical-success argument. Dossier rule 7: clinical failure is not evidence against tractability, so an efficacy claim touches no tractability number and is out at gate 1. |
| 2 | **Its support is secondary-only.** The link's `basis` is `background_only` or `hedged_only`. | `basis` is `primary` or `mixed`. A primary-supported claim is either usable or contradicted by something we can measure; neither is a *request*. **Exempt: the post-resolution contradiction ask**, which is the second branch of that same sentence — see below. |
| 3 | **We tried, we failed, and the failures are written down.** ChEMBL, the registry, the structure and a `grep` on the exact identifiers were all actually run, all abstained, and each null is in `not_found` **before** the ask is drafted. | any of them was not run. This is the gate that does all the work — see the never-ask list. **Exempt: the post-resolution contradiction ask.** |
| 4 | **There is a row to point at, and nobody has asked yet.** A `links[]`, `things[]` or `gaps[]` id, or a genuinely new question; and no matching entry in `rounds`. | `rounds` already carries the same verb at the same target — or, for a targetless `new_question`, the same **source identifiers**. Re-asking costs a round and returns the same evidence. See failure mode 18 for why `new_question` is matched differently. |
| 5 | **There is literature left to find.** `coverage.stop_reason != "complete"`, or `coverage.truncated` is true. | `stop_reason` is `complete`. The graph has exhausted the literature; another round returns what you already have. |

Gate 5 is the one that is easy to skip and cheap to check. It is also the only
gate that can make a *perfectly reasoned* ask worthless, so check it first — it
costs one field read and it retires the whole question.

**Gate 5's `or coverage.truncated` half is NOT implemented in
`graph_read.py`.** `LITERATURE_NOT_EXHAUSTED` tests `stop_reason != "complete"`
and nothing else, so a graph that is `complete` **and** `truncated` fails the
gate in the script while passing it in this table. That combination is a real
state — `coverage_notes()` reads `truncated: true` as "this is a sample, not the
literature" — and it is the one case where the script blocks an ask this rule
allows. Not fixed here because changing it moves `literature_not_exhausted` and
`mechanical_gates_clear` for every consumer of `--ask-context`; read the two
coverage fields yourself when they disagree.

```bash
python3 graph_read.py <graph.json> --ask-context
python3 graph_read.py <graph.json> --check-ask '<ask json>'
```

`--ask-context` reports gates 2, 4 and 5 per link. `--check-ask` gates one
proposed ask and exits 1 on any failure. **Neither checks gates 1 or 3, and it
says so in its own output.** Those are judgment, they are the two that stop this
becoming a way to avoid work, and no green result from a script substitutes for
them. Treat an all-green `--check-ask` as permission to *consider* an ask.

Two fields carry the exemption above, and both are mechanical facts, not
permission:

- `--check-ask` returns `exempt_for_this_ask`, populated only when `question`
  declares the ask post-resolution. When it is populated, gates 2 and 3 move out
  of `not_checked_here` and into it, each with the reason. Gate 1 never moves.
- `--ask-context` returns `clear_if_post_resolution_contradiction` per link:
  true where every mechanical gate **except gate 2** is clear, i.e. where the
  only thing blocking an ask is that the support is primary or mixed. Without
  it, `L4` — the ask that worked — and `L2` — the ask that must never fire —
  were both just `mechanical_gates_clear: false` and indistinguishable.

### The one ask that skips gates 2 and 3

When we **did** settle a claim ourselves and the settlement **contradicts a row
in the graph**, issue `resolve_link` on that row anyway, carrying our answer.

This is the only ask allowed past gate 3, and only because it inverts it: gate 3
exists to stop us outsourcing work we did not do, and here the work is done. The
ask has no bearing on our dossier — that field is already filled — and exists so
the wrong row does not propagate to every other consumer of the graph on the
next round. Mark it plainly as post-resolution so nobody upstream reads it as a
blocker. An ask of this kind must state our answer and its source, not just the
doubt.

**It skips gate 2 as well, and the reason is the same one.** An ask that carries
our answer is **not a request for work — it is a correction, and the gates were
written for requests.** Every one of them asks some version of "should somebody
else spend a round on this?", and that question does not apply to a statement of
what we already measured.

Gate 2 was the last thing blocking it. The exemption was written for gate 3
only, so a contradiction ask still had to clear "the link's `basis` is
`background_only` or `hedged_only`" — and the row you are most likely to be
correcting is a `primary` one, because a `primary` basis is what makes a wrong
row look authoritative to every downstream consumer. **A wrong `primary` row is
more damaging than a wrong `background_only` one, not less.** Gate 2's own
rationale already says so — "a primary-supported claim is either usable or
**contradicted by something we can measure**" — and then routed the second half
of its own sentence nowhere.

The worked case: `L4` in `fixtures/upstream_graph_askback.json` asserts a 2.4 Å
TL1A/DR3 crystal structure burying 1840 Å², at `basis: primary`. A PDB census
against DR3's accession returns three entities, all intracellular death domain,
no ectodomain, no complex of any kind. The claim is false and we can show it.
Under gate 2 as written, that ask could not fire. The worked ask is
`.claude/skills/ppi-hypothesis/fixtures/worked_ask.json`.

**What replaces the two exempted gates is stricter, not looser.** All four
required, and failure mode 19 is what happens when they are not:

1. `question` says **post-resolution**, in those words, so nobody upstream reads
   it as a blocker and so `--check-ask` can see it;
2. it states **our answer, its source and its date** — not the doubt;
3. it blocks **no** dossier field: every field the row touches is already filled
   from our own measurement, and the ask goes in `not_found[]` as a record;
4. every unexempted gate still applies — an id to point at, no prior round, a
   non-exhausted graph, two source identifiers.

Nothing here licenses a *forward* ask on a primary row. "Our prediction disagrees
with this row, please check" is not a correction, it is a hypothesis wearing
one — the PPI panel measured that direction and issued **zero asks in fifteen
cases**, because a predicted interface is never a thing to send upstream. The
exemption covers a row we **falsified with a measurement**, and nothing else.

### Make the ask carry its own evidence

*"Is obefazimod a TL1A agent?"* is not an ask, it is a shrug with an id attached.
The answering team cannot see our run, so an ask that does not carry the dispute
makes them redo the search that produced it.

Every ask states four things in `question`:

1. **the claim**, as the graph has it;
2. **the sources asserting it**, by id;
3. **the sources against it**, by id, and what we measured ourselves;
4. **what would settle it** — named, and phrased as something primary literature
   could actually contain.

`--check-ask` enforces a weak floor on this: at least two source identifiers
(`PMC…`, `NCT…`, a DOI, a `CHEMBL…` id, a PDB code) must appear in `question`.
Two identifiers is a floor, not a target, and it cannot see whether they sit on
opposite sides of the dispute. Clearing it is not the same as writing an
answerable ask.

**And the floor is looser and blinder than it reads.** Measured on the worked
TL1A ask, which clears it 11 tokens to 2:

- the PDB-code pattern is `[0-9][A-Z0-9]{3}`, which matches **any bare
  four-digit number**. Seven of that ask's eleven "identifiers" are `1840`,
  `2049`, `2826`, `1024`, `8863` — buried-surface measurements — plus `2026`,
  a date. Only five are real PDB ids. Any ask quoting two numbers clears the
  gate.
- **UniProt accessions are not recognised at all**, nor is anything lowercase.
  `Q93038` is the accession the entire DR3 census is keyed on and the single
  most load-bearing identifier in that ask, and it counts for nothing; so do
  `O95150` and `O95407`, and so does `3k51` written in lower case.

So the gate admits an ask carrying two measurements and no sources, and can
refuse one whose evidence is an accession census. It is a syntax check on a
string. It is not a check that the ask names its sources, and it was never able
to be.

### Where a pending ask lives, and why it can never block

**In the dossier: one entry in `not_found[]`.** No template change is required
and none should be made — `not_found[]` already exists, already carries
`{field, reason}`, and is already the place a reader looks for what the dossier
could not establish. Prefix the reason so it is machine-findable:

```json
{
  "field": "target_precedent.clinical_stage_small_molecules",
  "reason": "ASK[resolve_link:L3] issued to graph g_tl1a1 round 3. Not blocking: the field is filled from what we could verify and this entry records the residual question. <the question text>"
}
```

**In the graph: nothing.** This skill does not write to the upstream graph, and
that has not changed — `rounds` is the upstream team's record of asks they have
serviced, and they append to it when they answer. We emit the ask object in the
nomination's existing `asks: []` and stop there. Writing our own entry into
`rounds` would claim a round that was never run.

**The rule, and it is absolute:**

> An outstanding ask never blocks a verdict and never licenses a null.

Complete the dossier as if the ask will never be answered. Every field the ask
would have improved is filled from what we do have, or nulled with a reason that
stands on its own — a reason that would still be true if the ask had never been
written. If the *only* reason a field is null is that an ask is pending, the ask
is illegitimate: it means gate 3 was skipped and the resolution work was not
done. And `verdict: insufficient_evidence` may never be reached by counting a
pending ask as missing evidence. Rule 11 makes insufficient evidence a correct
answer; it does not make it a waiting room.

## What must never become an ask

Longer than the trigger on purpose, and the more important half. The trigger
describes a rare event. This list describes what an agent holding a new button
will actually reach for.

The failure is not that a bad ask gets sent. It is that **the moment an agent
can ask, the cheapest way to finish any hard sub-problem is to ask** — and the
result is a pipeline that has kept its verdict and given away its judgement. We
are the structural and chemical instrument. Everything below is ours.

**1. What a drug actually hits.** `chembl.drug_mechanism`, joined on `molregno`.
This is the single most tempting ask and it is a local table. Measured: the
obefazimod/TL1A trap — the case that motivated this whole section — is settled
by one query returning one row in 6 ms:

```sql
SELECT molregno, tid, mechanism_of_action, action_type
FROM chembl.drug_mechanism WHERE molregno = 2335315
-- 2335315 | 120082 | Cap binding complex modulator | MODULATOR
```

Not TL1A. Never ask "what does compound X target".

**2. Modality.** `chembl.molecule_dictionary.molecule_type`, read per drug
(dossier rule 1). Authoritative for `Small molecule`, `Antibody`, `Protein`. A
returned `Unknown` is *still* not an ask — it maps to modality-unknown, goes in
`not_found`, counts toward neither block, and the optional Open Targets
corroboration (rule 10b) is also ours to run.

**3. Trial status, phase, and why a program stopped.** `ctgov.studies`, and the
`terminated-programs` skill exists entirely for this. Measured on the two cases
that looked like they needed help: LY3509754's "conflicting" trial IDs are two
distinct trials, both `Terminated`, both carrying a `why_stopped`, returned by
one query in 418 ms; zimlovisertib's Phase 2 record is three `Completed` studies
plus one terminated for lack of enrollment, in 437 ms. Neither was ever a
literature question.

**4. Assay provenance.** Whether an actives count collapses to one assay, and
whether that assay measures a different protein. Falsification check 1, one
`GROUP BY`. The IRAK4-assay-inside-TNF-alpha finding came from that query, not
from a review.

**5. Anything structural.** Pocket volume, druggability range, cryptic mechanism
and its apo census, interface classification, chain selection, ligand identity,
whether a holo ligand is a frequent hitter. This is the axis we exist to
compute. An ask here is not delegation, it is abdication — and the graph team
has no structures anyway.

**6. Anything computable.** A fraction, a CV, a spread, a distance, a count, a
ratio. Arithmetic is never a literature question. If the inputs are in hand the
answer is too.

**7. Accession resolution.** `uniprot_v.proteins`. Where two accessions both fit
the quotes, step 4 already gives the correct output: populate `ambiguity`, leave
`uniprot_accession` null. That is a finished answer, not a pending one. Note the
narrow legitimate case that *sounds* the same: asking `resolve_link` because the
**evidence** on a link is review-only is fine; asking "which accession is it" is
not. The first is about what the literature says, the second is a lookup.

**8. A verdict, a score, a ranking, or an interpretation.** They map literature
relationships; they do not assess tractability, and they have never claimed to.
"Is this target druggable?" is a category error in both directions — it asks
them for our output using their input. This also rules out the softer forms:
"how strong is this evidence", "should we count this", "which of these matters
more".

**9. Efficacy and clinical-outcome arguments.** Gate 1. PMC12325316's claim that
IRAK4 kinase inhibition "cannot completely block TLR signalling" is a real and
interesting argument that bears on whether the drug works, and bears on no
number in this dossier. Reporting it unresolved is the correct output and always
was.

**10. Anything already in `rounds`.** Gate 4. Check before drafting, not after.
For `new_question` this check is weaker than it looks and you must read its
detail string rather than its boolean — failure mode 18.

**11. Anything against a `complete` graph.** Gate 5.

**12. A question whose failed attempts are not written down.** Gate 3, stated as
a prohibition because it is the one that gets rationalised away. If `not_found`
does not already carry the null results from ChEMBL, the registry and the grep,
there is no ask — there is an unfinished lookup. Write the nulls first; often
the act of writing them produces the answer.

**13. Retraction status.** `papers[].retracted` is in the graph already and
`graph_read.py` surfaces it in `retracted_papers`. Read the field.

**14. An unknown verb at intake.** "Adjudicating an unknown verb" already routes
this to `resolve_link` at intake. Do not issue a second, build-time ask for the
same link — that is one question consuming two rounds, and `rounds` will show it.

**15. Anything you have not yet finished reading.** The obefazimod trap was
caught by reading further, in our own run, in the corpus we already had. An
agent with an ask button would have pressed it at precisely the moment before
the reading that solved it. Before any ask: is this unresolvable, or merely
unread?

## Adjudicating an unknown verb

`how` has no enum. Every other categorical field in SCHEMA.md carries an explicit
`a|b|c` comment; `how` does not, because it is open vocabulary written by the
upstream extraction model. `DIRECT_ACTION` holds seventeen verbs and
`DOWNSTREAM_EFFECT` fourteen — several of them only `-ise`/`-ize` spellings of
each other — against an unbounded space, so the two sets can never be complete.
We cannot ask the upstream team to close the vocabulary. This is ours to solve.

So an unmatched verb is not dropped and is not a rare edge. It is a target the
intake could not classify, and the helper hands it back in `needs_adjudication`
with a `signals` block and a `decide` field pointing here.

The decision is one thing: **is this edge a direct action on a target, or a
downstream effect on a readout?** Getting it wrong sends the dossier to score
pockets on a secreted cytokine — failure mode 2, arriving through a verb the
helper had never seen.

`upstream_graph_unknownverb.json` exercises this path directly, with three
invented verbs chosen to land in three different signal states: `clamps` on a
kinase with an IC50 and a biochemical `where` (direct), `quenches` on IL-6 with
ACR50 and serum levels (downstream), and `perturbs` on MRTFA with quotes that
match nothing either way (refuse). All three objects are `protein` or `gene`, so
`eligible_kind` is true for all three and cannot break the tie.

The table below uses the RA graph instead, because `L1` (zimlovisertib → IRAK4)
and `L3` (zimlovisertib → IL-6) are the shape-identical pair an adjudication has
to separate, and every value in it is real output from `signals()` rather than
an illustration.

| field | evidence for | how it misleads alone |
| --- | --- | --- |
| `object_kind` | only `protein` or `gene` can be a target at all | IRAK4 and IL-6 are **both** `protein`. Kind never separates a target from a readout; it only rules out `process` and `disease` — see failure mode 7. |
| `object_has_edge_to_disease` | a readout usually carries the chain onward. IL-6 has `L6` (IL-6 drives RA); IRAK4 has none. | a well-studied target has disease edges too. IRAK4's is empty only because no paper here wrote "IRAK4 drives RA". One more round and it would not be. |
| `assay_contexts` | the `where` string of every finding on the edge | descriptive free text, not a category. `L1` carries `human whole blood` and is a direct action; `L3` carries `LPS-stimulated whole blood` and is not. Cellular does not mean downstream. |
| `direct_context` | a `where` naming a biochemical, cell-free or purified system. `L1` has `biochemical assay`. | empty means the field was empty. `f3` and `f4` have `where: null` outright, so a blank here is a missing field far more often than an absent experiment. |
| `direct_terms_in_quotes` | the vocabulary of a binding measurement. `L1` matches `ic50`, `target engagement`, `kinase activity`. | `L3` **also** matches `target engagement` — from `f10`, *"serum IL-6 trended lower in the treatment arm, which may reflect target engagement"*. The engagement is against IRAK4, and the sentence is hedged. Term matching cannot see whose engagement, or whether anything was asserted. |
| `downstream_terms_in_quotes` | outcome vocabulary. `L3` matches `release`, `output`, `serum`. | `levels` and `expression` sit happily in a direct quote — "receptor occupancy levels" matches. One hit is not a readout. |

These are evidence, not a verdict, and no single one decides. `L3` carries a
direct term and `L1` carries a cellular assay context; either read alone gets the
call backwards.

The procedure weights the quote above all of it:

1. **Read the quotes, all of them, first.** Every signal is an index computed
   over the quotes. The quote is the one thing SCHEMA.md guarantees is
   verbatim — its guarantees section says every quote is string-matched against
   the fetched abstract before a finding is written, and claims that fail to
   match are dropped and counted in `coverage.no_quote_discarded`. Nothing else
   in the packet carries that guarantee.
2. **Ask what was measured, and on what.** *"inhibited IRAK4 kinase activity
   with an IC50 of 0.2 nM"* measures the object of the edge. *"suppressed
   LPS-induced IL-6 release"* measures what happened downstream of something
   else that was measured. That distinction survives any verb the extractor
   invents.
3. **Then read the signals, as corroboration for the reading you already have.**
   If they contradict the quote, re-read the quote — do not switch on the
   signals. They are keyword matches; the quote is the sentence.
4. **Check `says` and `hedged` before acting either way.** `f10` is hedged
   (*"may reflect"*), and a hedged quote settles nothing in either direction.
   `no_effect` on a direct-action edge is real tractability evidence; on a
   downstream edge it is someone else's question — failure mode 10.

Refusing is a legitimate outcome. When the quotes do not settle it, leave the
edge unresolved, report it as unresolved, and issue `{"ask": "resolve_link",
"target": "<link id>", "depth": "deep"}` on that link. Do not promote it to a
nomination and do not file it under `rejected`; both are guesses wearing an
answer's clothes. This is the same rule as accession ambiguity in step 4 — an
unresolved target is a correct output, and a confidently wrong one poisons every
number downstream of it. A refused edge costs one round. A wrong one costs the
whole dossier.

## Failure modes

Longest section on purpose. The procedure above is the easy half.

### 1. Reading `how` instead of the quote

`L1` is `zimlovisertib inhibits IRAK4`. That verb supports a nomination and
nothing else. Every mechanism claim in this graph lives in quote text. An intake
built on the verb vocabulary returns targets with no `interaction_to_disrupt`
and silently hands the dossier its weakest input.

### 2. Readouts look exactly like targets

In the RA fixture, `L3` is `zimlovisertib reduces IL-6`. A small molecule acting
on a protein — the same shape as `L1`. But `reduces` is an outcome measurement
and IL-6 is a downstream cytokine, not the drug's target.

This is the same error the dossier already defends against one stage lower.
`targets.json` records that 45% of TNF-alpha's bioactivity comes from an "IRAK4
Monocyte TNFalpha Cell Based Assay" measuring a different protein, with TNF as
the cellular readout. Same conflation, one stage earlier, and here nothing is
looking for it. The verb split in the helper is the whole defence — keep the
two verb sets separate and send unrecognised verbs to `needs_adjudication`,
where they get decided rather than dropped.

### 3. The name resolves to a different protein than the evidence

The worst failure in this skill, because it is silent and everything downstream
still runs.

Node `t5` is named `IL-6`, alias `interleukin-6`. That string-matches **P05231**
(Interleukin-6), the ligand. Its supporting quote, `f7`, reads *"IL-6 **receptor**
blockade reduced ACR20 non-response across 14 randomized trials"* — tocilizumab
and sarilumab, which target **P08887** (IL-6 receptor subunit alpha). Those are
different proteins with different tractability.

`rheumatoid_arthritis.json` lists P08887 in `biologic_only` with an expected
verdict of `not_tractable OR insufficient_evidence`. Resolve by name and the
dossier assesses a protein the evidence never supported, while the target the
evidence does support is never assessed at all.

The same trap appears twice more in that fixture: anakinra targets IL1RN, not
IL1B; brodalumab targets IL-17RA, not IL-17A. Always read the quote before
resolving.

### 4. A review citation that looks like a finding

IRAK4's two PPI statements — `f4` (nucleates the MyD88 signalosome) and `f3`
(recruited to the receptor complex) — carry finding confidences of 0.88 and
0.85. Both come from one review, both are flagged `background`, and both sit on
links with `basis: background_only` and link confidence 0.38 and 0.35.

Taking `interaction_to_disrupt: "MyD88 signalosome assembly"` from that is a
single secondary source deciding chain selection, and therefore the druggability
number. Record it, ask `resolve_link` on L2, and do not act on it until a
primary result comes back.

### 5. Catalytic function does not imply orthosteric

`interaction_to_disrupt: "kinase activity"` says what to stop. It does not say
where to bind.

TYK2 is in your own fixture set for exactly this reason: deucravacitinib is an
approved kinase inhibitor that binds the **JH2 pseudokinase domain**
allosterically, not the ATP site. Structures split 29 entries for JH1 and 21 for
JH2, and picking the wrong domain scores the wrong pocket.

So set `mechanism_hypothesis: "unknown"` unless a quote states ATP-competitive,
allosteric, or a residue range. The dossier already handles `unknown`: it reports
pockets for the biological assembly and records in `tractability.caveat` that no
mechanism was specified. That is a correct output. A guess is not.

### 6. One target, two functions, and a drug that only hits one

IRAK4's four links say two different things about it. `f1` and `f5` are about
kinase activity. `f4` is about scaffolding — nucleating the signalosome. These
are separable functions, and a kinase inhibitor stops only the first.

No single link states this. It appears only when a target's whole neighbourhood
is read at once, which is why the helper returns the neighbourhood rather than
one edge. Report every function shape the evidence supports, and never let the
drug's mechanism stand in for the target's biology.

For this graph it also matters for the question being asked. The graph asks
whether the effect is confined to the myeloid compartment. A second explanation
sits in the same data: the molecule may be blocking half of what IRAK4 does, in
every compartment. That belongs in `notes`, not in a verdict.

### 7. Complexes and pathways are not proteins

`t3` is "MyD88 signalosome" and `t7` is "TLR/IL-1R signaling". Both are typed
`process`. Neither has a sequence, so neither can be handed to a structural step.

`t3` resolves to a component list (MYD88 Q99836, IRAK4 Q9NWZ3, IRAK1 P51617,
IRAK2 O43187). `t7` is a pathway name and does not resolve at all — TLR4
(O00206) and IL1R1 (P14778) are both plausible and the graph gives no way to
choose. Return a list with `ambiguity` populated; do not pick.

### 8. Synthetic fixtures

A graph carrying `_fixture: true` has papers, DOIs and quotes that were never
retrieved from any corpus. The helper refuses these unless `--allow-fixture` is
passed. Never lift that guard for a real run, and never cite a fixture quote.

### 9. Computed results are not literature evidence

If a downstream method (pocket scan, protein-protein cofold, interface
computation) produces a result that gets written back into a graph, it must not
re-enter this intake as a `finding`. Their `findings` require a verbatim quote —
`coverage.no_quote_discarded` shows the upstream pipeline drops entries without
one — and a computed result has provenance but no quote.

The specific hazard is laundering: a cofold contaminated by PDB training data
becomes an ordinary-looking finding on the next round, and the `leakage_risk`
flag the dossier template already carries is gone. Keep computed results in
their own type, with `leakage_risk` required rather than optional.

### 10. `no_effect` is not `no`

`says` has three values — `yes | no | no_effect` — and `links` carries three
arrays to match. A `no` finding is evidence against a claim. A `no_effect`
finding is a measurement that came back null. They are different, and folding
them together loses the more useful one.

For this node the difference is load-bearing. A `no_effect` finding on a
direct-action edge says the compound did not engage the target, which is real
tractability evidence. On a downstream edge it says the pathway did not move,
which is a biology result and someone else's question.

Two consequences. Read all three arrays — an intake that reads `yes` and `no`
silently drops every null result. And when the shared `Evidence.direction`
mapping is settled, argue for a third value rather than mapping `no_effect` onto
`contradicts`.

### 11. The ask becomes the thing you do instead of the work

This is the failure mode the ask-back mechanism introduces, and it is worse than
the gap it closes, because it is invisible in the output. A dossier with a
well-formed pending ask and a null beside it looks *more* rigorous than one with
a null and no ask. It reads as diligence. It can be the opposite.

Three symptoms, all measured against the real cases this section was built from:

- **An ask on a question a local table answers.** Every one of the four cases
  that motivated this mechanism turned out to be answerable inside the pipeline
  — obefazimod by `chembl.drug_mechanism`, both trial cases by `ctgov.studies`,
  the IRAK4 efficacy claim by not being a tractability question at all. A rule
  written from those cases and *not tested against them* would have fired on all
  four and offloaded four questions we can answer in under a second each.
- **`not_found` growing while `sources` does not.** If a run's asks outnumber
  its newly-sourced numbers, the sweep stopped retrieving and started
  forwarding.
- **An ask drafted before the nulls are written.** Gate 3 requires the failed
  attempts in `not_found` first. Drafting the ask first and backfilling the
  nulls to justify it produces an identical-looking dossier and inverts the
  reasoning.

The defence is gate 3 and it is deliberately expensive: you must have run the
lookups and recorded their nulls before you may write the ask. If that feels
like it defeats the purpose, that is the correct feeling. The mechanism is meant
to fire rarely — a handful of times a year, on claims that genuinely have no
answer inside our instruments — not to be a routing layer for hard questions.

### 12. The gap door — a readout nominated as a target, with its rejection deleted

The nomination rule has two halves and the second one has no evidence behind it.

**Upstream `gaps` are not a curated list of undrugged candidates.** They are
structural open triangles, computed from graph shape alone:
`assemble.py find_gaps()` emits a gap for every pair of things that share a
neighbour and have no edge to each other. Nothing about a gap says the pair is
biologically interesting, and nothing says either half is a target.

So run **their** generator over **our** RA fixture, unchanged, and it emits:

```
g5  missing: [t1, t5]   = [IRAK4, IL-6]   confidence 0.16
```

— because both connect to zimlovisertib. IL-6 is then nominated as a target,
through the gap half, on no evidence at all. And the old code then ran
`rejected.pop(tid)`, which **deleted** IL-6's `"reached only by downstream-effect
edge L3 (reduces) -- readout, not target"` reason. Nomination present, rejection
gone, exit 0.

That is failure mode 2 arriving through the door built to catch undrugged
candidates, and our shipped fixture never caught it because the gaps in it
(`g1 [t2,t4]`, `g2 [t1,t4]`) were **hand-authored by us** and happen to contain no
readout. We wrote both sides.

A gap nomination now carries `evidence_class: "structural_gap_only"`, a rejection
is never deleted, and a thing that is both appears in `contested_nominations`.
**A gap-only nomination is a question, not a candidate.** Treat it as
`test_gap` material, never as an input to the dossier on its own.
`upstream_graph_gapdoor.json` is the fixture.

### 12b. Zero nominations from a perfectly healthy graph

Run this intake against `g_1a4f` — the only real graph the knowledge-graph team
has produced — and it returns **zero nominations**, `status: "ok"`, exit 0.

Not because the literature is empty. Because that graph models things at a
coarser grain than we do. Its five `things` are:

| id | name | kind |
| --- | --- | --- |
| t1 | IRAK4 inhibition | `small_molecule` |
| t2 | myeloid inflammatory signalling | `process` |
| t3 | synovial fibroblast driven inflammation | `process` |
| t4 | MyD88 dimerization inhibition | `small_molecule` |
| t5 | TLR/MyD88/NF-kB signalling axis | `process` |

**There is not one `protein` or `gene` node in it.** The extractor nominalises:
it makes nodes out of *interventions* ("IRAK4 inhibition", typed
`small_molecule`, with the actual compounds PF-06650833 and KIC-0101 demoted to
`aliases`) and out of *phenotypes*, never out of the protein. IRAK4, MyD88 and
IRAK1 exist in that file only as substrings inside other nodes' names.

Our nomination rule requires `kind in {protein, gene}`. Against this graph it can
never fire, on any round, for any target.

**This is the single most likely thing to go wrong when a real graph arrives, and
before the fix we would not have noticed** — zero nominations with `status: ok`
reads exactly like "no druggable targets in this literature." It now reports
`selection.n_candidates: 0` with a note naming granularity mismatch as one of the
things that lands there. When you see it: read `rejected`, check the `kind`
distribution yourself, and if there are no `protein`/`gene` nodes at all, that is
an `expand_node` ask on the intervention node — not a finding about the target.

### 12c. The verb decided, and the quote was never read

`classify()` was a pure function of (subject `kind`, verb). The quote never
reached it. So the graded IL-6 negative passed for one reason only: `reduces`
happens to sit in `DOWNSTREAM_EFFECT`.

Measured across 18 shape-matched pairs — same grammar, same subject, only the
verb and the readout/target distance varying — the boundary was a straight line
with no gradient:

| what happened | how often | class |
| --- | --- | --- |
| a `DIRECT_ACTION` verb on a readout → **nominated as a target** | 6 of 6 | silent wrong |
| a `DOWNSTREAM_EFFECT` verb on a real target → **rejected as a readout** | 2 of 2 | silent wrong |
| a verb in neither set → `needs_adjudication` | 2 of 2 | loud, correct |

`blocks IL-6 secretion`, `inhibits IL-6 production`, `modulates IL-6 output`,
`engages IL-6 (inferred because serum levels fell)` all nominated a secreted
cytokine as a druggable target. `suppresses IRAK4 catalytic activity, Ki 1.1 nM`
and `attenuates IRAK4 enzymatic activity in a purified system` both rejected a
kinase as a readout.

The distance between readout and target made no difference at all, because the
distance lives in the quote and the quote was not consulted. **The only thing
that mattered was which set the verb fell into** — and per the vocabulary note in
step 1, that is close to a coin flip on real upstream output. The graded negative
is one verb away from failing, and the verb it needs (`inhibits` for "inhibited
IL-6 production") is the one upstream's own extraction skill tells the model to
emit.

`signals()` — the whole quote-reading apparatus — was only ever computed for
*unclassified* verbs. On a recognised verb it never ran.

Now it always runs. Where the verb and its own quotes disagree the edge goes to
`needs_adjudication` with `why_contested`, and it is neither nominated nor
rejected. Where they agree the reading travels with the nomination as
`quote_reading`.

**Read `quote_reading: "silent"` carefully.** It is the residual boundary: the
quote carries no vocabulary either way, so nothing but the verb supports that
nomination. `targets IL-6 driven inflammation in synovium` reads silent and still
nominates, because there is genuinely nothing in the sentence to catch it. That
one is yours. `upstream_graph_verbquote.json` is the fixture.

### 13. Every id resolves — except it doesn't

SCHEMA.md's guarantees section says "Every id resolves. `from`, `to`, `paper`,
and the ids inside `yes`, `no`, `no_effect`, `implied_by` all point to a row in
the same file."

It is not enforced. The one real graph we have breaks it on arrival:
`rounds[1].target` is `"g3"` against a `gaps` list of `g1`, `g2`, `g4`. It also
carries `f6`, a finding no link references, invisible to this intake entirely.

Before the fix, each of these was silent:

- a link whose `to` did not resolve **vanished** — no nomination, no rejection,
  no adjudication, nothing anywhere in the output;
- a link citing a nonexistent finding skipped it without comment;
- a finding citing a missing paper came back with `paper_ref: null` **and
  `retracted: null`**, so a retracted source reads as clean;
- a duplicate id **overwrote last-wins**, which is the worst of them: declare
  `t1` twice and the second row silently relabels the whole neighbourhood, so
  `zimlovisertib inhibits IRAK4` renders as `zimlovisertib inhibits IL-6` and the
  agent resolves an accession for the wrong protein off the right evidence. That
  is failure mode 3 produced mechanically.

All of it now lands in `integrity`. None of it stops the run — a graph with one
dangling id is still mostly usable — but **a nomination that touches one of these
rows is not trustworthy**, and the block tells you which rows those are.
`upstream_graph_integrity.json` is the fixture.

### 14. `coverage` fields that were never read

Three of them, each one a way to read an absence as a result:

- **`no_quote_discarded`** — the upstream pipeline's *only* documented removal.
  A claim whose quote will not string-match the fetched text is dropped, never
  repaired. It was surfaced nowhere. A graph that discarded 41 claims looked
  identical to one that discarded none, and anything this intake reports as
  "not stated" may be among them.
- **`depth: "quick"`** — SCHEMA.md note 2: quick reads page 1, and page 1 lies.
  At `quick`, absence means unknown *whatever* `stop_reason` says.
- **an absent `coverage` block** — `graph.get("coverage", {})` returned `{}`, so
  `stop_reason` was `None`, so no warning fired, so a graph with no coverage
  block at all read exactly like a `complete` one.

The trap is the combination: `stop_reason: "complete"` with `truncated: false`
made the old `coverage_warning` null, so a page-1 skim that threw away 41 claims
presented as an exhausted search. Read `coverage.warnings`; empty is the only
clean state, and `literature_exhausted` requires all three of untruncated,
`complete`, and not-quick. `upstream_graph_coverage.json` is the fixture.

### 15. Three candidates and no way to choose

The dossier Contract takes exactly **one** `uniprot_accession`. The graph can
nominate any number and offers nothing to rank them by — and this skill does not
rank, deliberately.

That was previously a plain list sorted by `thing` id, with no field saying a
choice was required. `thing` id order is upstream insertion order and carries no
meaning whatsoever, so an agent taking the first took an arbitrary one, silently,
with no record that three others existed.

`selection` now reports `n_candidates`, `basis_to_choose: null`, and splits
`with_direct_action_evidence` from `gap_only_no_evidence` — because a gap-only
nomination must never sit in a list beside targets with measured IC50s as though
they were the same kind of thing.

**What to do:** run the dossier once per candidate if the question allows it;
otherwise state the ambiguity, populate `ambiguity`, leave `uniprot_accession`
null, and issue one ask. Same rule as accession ambiguity in step 4 — an
unresolved target is a correct output. Never take the first.
`upstream_graph_multitarget.json` is the fixture.

### 16. A negative result nominated as a target

SCHEMA.md is explicit: a negation lives in `says`, never in the verb — there is
no `does_not_bind`. So a link reading `zimlovisertib inhibits IRAK4` whose
findings **all** say `no` is the literature reporting that it does *not* inhibit
IRAK4, and it produced a nomination indistinguishable from a positive one. The
reason string said `"object of direct-action edge from small_molecule t2
(inhibits)"` and nothing about direction.

Nominations now carry `support: {yes, no, no_effect}` and a `direction_warning`
when `yes` is empty. Keep the two empty-`yes` cases apart, per failure mode 10:
all-`no` is evidence against the relation itself; all-`no_effect` on a
direct-action edge is a measured null and **is** real tractability evidence.
`upstream_graph_negative.json` is the fixture.

### 17. The fixture guard only fires on a flag the producer never sets

`--allow-fixture` refuses a graph carrying `_fixture: true`. That guard works,
and it has now been exercised: a truthy non-boolean (`"true"`) is correctly
refused, and `_fixture: false` or an absent key correctly passes.

But note what it means. `g_1a4f`, the real graph, carries **no `_fixture` key at
all** — a real graph is distinguished from a synthetic one only by the *absence*
of a flag that we invented and upstream has never heard of. The guard protects
against our own fixtures leaking into a real run. It provides no assurance
whatsoever in the other direction, and nothing in the pipeline verifies that an
unflagged graph came from a real corpus.

### 18. One `new_question` retired the verb for the life of the graph

`already_asked()` matched a proposed ask against `rounds` on **(verb, target)**.
Three of the four verbs have a target. `new_question` does not — SCHEMA.md gives
it `target: null` by design, because it is the ask for a claim the graph has no
row for.

So every `new_question` ever issued is `("new_question", None)`, and the first
one asked against a graph made **every** later one fail `NOT_ALREADY_ASKED`,
exit 1, forever. Measured on `upstream_graph_askback.json`, whose `rounds[0]` is
a `new_question` from round 1: an unrelated structural `new_question` came back
`FAIL NOT_ALREADY_ASKED`, exit 1, with the detail string
`rounds already carries 1 ask(s) of new_question at None` — pointing at a
question that had nothing to do with it.

The verb that carries every claim the graph has no row for was the one verb that
could be used once.

**The fix is not a hash of the question text**, which was the obvious repair and
is wrong in a way that would have looked right. `rounds` is not required to
record the question — the real ask-back fixture's round 1 records `ask`,
`target`, `depth`, `papers_added` and `outcome`, and no question at all — so a
text hash would match nothing and the gate would silently become a no-op that
always passes. And where the text *is* recorded, the upstream team rewords a
question when it services it, so a hash calls a rephrasing a different ask and
lets the same question through twice. Both failures are invisible.

A targetless ask is now matched on **the set of source identifiers `question`
names** (`question_identity`). Those are what make an ask routable at all —
`QUESTION_IS_ACTIONABLE` already refuses a question naming fewer than two — and
they survive rewording, reordering and translation, which the wording does not.

**And where a prior round records no question, the gate says so instead of
passing.** Those rounds come back in the detail string as rounds this gate is
structurally unable to compare against. That is an honest "cannot tell", and it
is the state the one real ask-back fixture is in. Read the detail, not the
boolean.

### 19. "Post-resolution" as a password

Failure mode 11 is the ask becoming the thing you do instead of the work. The
gate-2/gate-3 exemption above is the same hazard with a specific key: the two
gates that stop an ask being a way to avoid work are exactly the two an agent
can now switch off, and it switches them off **by putting a word in the
`question` text**. `--check-ask` cannot tell a correction from a claim to be
one. Nothing mechanical can.

What separates them is not the wording, it is whether an artifact exists:

- a real one names **what we measured, with what, on what date, and what it
  returned** — a PDB census, a ChEMBL row, a registry record — and the number it
  contradicts;
- a fake one names our *disagreement*: a prediction, a confidence, an
  expectation, an "inconsistent with", a "we would expect". A cofold is not a
  measurement of what is deposited. A prediction contradicting a graph row is
  failure mode 11 with a new label, and the PPI panel is the precedent: fifteen
  cases, a predicted interface every time, **zero asks**, because a prediction is
  never a correction.

Two tells, both cheap. **The ask must be writable before the exemption is
claimed** — if you cannot state our answer without the phrase "post-resolution"
doing the work, there is no answer. And **it must block nothing**: a
post-resolution ask beside a null field is a contradiction in terms, because the
whole premise is that we already filled that field ourselves. If the field is
null, gate 3 was skipped and the exemption was not earned.

## What this skill does not do

- It does not assess tractability. That is the dossier.
- It does not rank nominations against each other.
- It does not choose between two plausible accessions.
- It does not assign a mechanism the evidence has not stated.
- It does not write anything back into the upstream graph, `rounds` included.
- It does not ask upstream anything a local table, the registry, a structure or
  arithmetic can answer, and it does not let a pending ask block a verdict or
  justify a null.
