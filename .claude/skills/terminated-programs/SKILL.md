---
name: terminated-programs
description: >
  Retrieves why clinical programs against a target stopped — the registry
  termination reason from ClinicalTrials.gov (ctgov.studies.why_stopped) joined
  to the literature account of the same event, recovered by exact-string grep on
  drug codes rather than semantic search. Reports contradictory sources side by
  side with their dates instead of picking one. It does NOT score druggability,
  does NOT decide whether a target is dead, and does NOT resolve a contested
  status into a single answer.
---

# terminated-programs

A target can be perfectly tractable and still fail. RORγt has 152 holo
structures, ~12,900 ChEMBL compounds and sub-nanomolar leads, and zero
approvals — because VTP-43742 stopped on reversible transaminase elevations and
TAK-828F stopped on preclinical toxicology and clinical teratogenicity. A
dossier that reports "tractable" and stops there is true and useless.

Termination reasons are almost never in ChEMBL. They live in two places, and you
need both:

| Source | What it gives | What it misses |
|---|---|---|
| `ctgov.studies.why_stopped` | Sponsor's own words, dated, per-NCT, structured | Only ~8% of studies; only for Terminated/Withdrawn/Suspended; truncated at 250 chars; mostly business/enrollment noise |
| `/papers/` full text | The mechanistic reason, the number of affected patients, class-wide context | Undated relative to the event; reviews recycle stale text for years |

This skill fills `terminated_programs[]`. It does not touch any score.

## Setup

```bash
# Credentials come from the ENVIRONMENT, never from a hardcoded path. This
# skill bundle is uploaded and executed in a sandbox where no developer
# checkout exists, so `. /some/repo/.env` is not portable and must not appear
# here. If you are running locally, source your own .env before invoking the
# agent (`set -a; . "$REPO_ROOT/.env"; set +a`); when deployed, the key is
# injected as a secret.
: "${PAPERCLIP_API_KEY:?PAPERCLIP_API_KEY is not set — STOP. Do not run this skill without it.}"
paperclip skill            # run once per session
```

**If `PAPERCLIP_API_KEY` is absent, stop and report it as a blocking failure.**
Do not fall back to a path, a cached file, or "proceed and see what comes back".
Every query below returns zero rows without credentials, and zero rows here is
indistinguishable from *a target with no terminated programs* — a silent
false-negative that lands in the dossier as a positive finding. The correct
output when the key is missing is an error naming the missing variable, and
`terminated_programs` left `null` with the reason recorded in `not_found`.

`-s` is **required** on every `paperclip search`. It is **not** used by `grep`,
`sql -s trials`, `cat`, or `ls` — those take a path or a source flag of their own.

---

## Rule 1 (the important one): grep drug codes, do not search them

`paperclip grep` is a **full-text regex over every paper body in the corpus**. A
drug code is an exact string. Semantic search embeds it into a topic vector and
returns papers about the topic, not about the compound.

Measured, same session, same corpus:

```bash
paperclip search -s pmc "LY3509754 IL-17A small molecule hepatotoxicity" -n 5
#   -> 5 hits: macrocyclic IL-17a modifiers, alisertib/doxorubicin hepatotoxicity
#      in mice, BMS-986094 mitochondrial toxicity, ...
#      ZERO of them mention LY3509754.

paperclip grep "LY3509754" /papers/
#   -> 13 papers. First snippet: "a small molecule inhibitor of IL-17, LY3509754,
#      progressed to clinical trials but was halted during Phase 1 due to
#      unfavorable hepatotoxicity"
```

```bash
paperclip grep "VTP-43742" /papers/     # 39 papers
paperclip grep "TAK-828"   /papers/     # 25 papers
paperclip grep "SPD304"    /papers/     # 76 papers
```

**Use `search -s pmc` only for topics** ("oral RORγt inhibitor safety in
psoriasis"), never for identifiers. Use `grep` for: drug codes, INNs, NCT
numbers, gene symbols, accessions, patent numbers, author surnames.

### Regex grep is how you get the reason, not just the mention

A bare grep on a code returns dozens of papers that merely cite it. Bind the code
to termination vocabulary in one regex and the corpus collapses to the papers
that say *why*:

```bash
paperclip grep -i "TAK-828.{0,300}(discontinu|terminat|halt|toxic|teratogen)|((discontinu|terminat|halt|toxic|teratogen)[^.]{0,300}TAK-828)" /papers/
# 25 papers -> 4 papers, one of which states the reason outright
```

The two-sided alternation matters: reviews write it both ways ("TAK-828F was
discontinued because…" and "…toxicity prompted discontinuation of TAK-828F").
`[^.]{0,N}` keeps the match inside one sentence; `.{0,N}` crosses sentence
boundaries and is looser — use `.` when the code and the reason sit in a table
row, `[^.]` when you want a real sentence.

`grep` reports `hit the per-shard match cap — more matches exist; raise -m N`.
Raise it (`-m 40`) before concluding a compound is absent.

### Then read the sentence with line numbers

```bash
paperclip grep -n "TAK-828F was discontinued" /papers/PMC8080595/content.lines
paperclip cat /papers/PMC8080595/meta.json     # doi, journal, pub_date
```

Every claim in the output block needs `PMC<id>` + line, a DOI, or an NCT number.
No exceptions, and no claim survives without one.

---

## Rule 2: query ctgov first — it is structured, dated, and free

`paperclip skill trials` **fails** (`Error: Failed to fetch skill "trials"`), and
trial VFS documents are empty shells:

```bash
paperclip cat /trials/us/NCT02817516/meta.json
# {"total_lines": 0, "total_blocks": 0, ... }   no dates, no text
```

Ignore both. `paperclip sql -s trials` exposes a full AACT-style `ctgov` schema
(71 tables) and that is the real interface. Sibling schemas exist for other
registries — `ct_cn` (ChiCTR), `ct_jpn` (UMIN/JRCT), `ct_eu` (EudraCT/CTIS/
ISRCTN), `ct_global` (WHO ICTRP) — but **none of them has a why_stopped / reason
column**. Verified:

```sql
SELECT table_schema, table_name, column_name FROM information_schema.columns
WHERE table_schema IN ('ct_cn','ct_jpn','ct_eu','ct_global')
  AND column_name ~* '(stop|reason|termin|status)'
-- returns only recruitment_status / overall_status / status / status_en
```

So a program that ran only in China, Japan or the EU is registry-silent on
*why*, and you are on literature alone (e.g. balinatunfib's CTR20241078, known
only from PMC11734553).

### `ctgov.studies` — the fields that matter

Verified population, whole table:

| Column | Populated | Notes |
|---|---|---|
| `nct_id` | 595,652 / 595,652 | |
| `study_first_posted_date` | **595,652 / 595,652 (100%)** | 1999-09-20 → 2026-07-24. This is your as-of date. |
| `overall_status` | 595,652 | see enum below |
| `why_stopped` | 46,557 (7.8%) | free text, ≤252 chars |
| `phase` | 454,428 | `Phase 1`, `Phase 1/Phase 2`, … , `Not Applicable` |
| `enrollment` | 588,526 | with `enrollment_type` = Actual / Estimated |
| `start_date`, `completion_date`, `primary_completion_date`, `last_update_posted_date` | high | `last_update_posted_date` is the freshness signal |
| `source`, `source_class` | 595,652 | sponsor name / industry-vs-academic |
| `last_known_status` | | what a study was doing before it went `Unknown` |
| `official_title`, `brief_title`, `acronym` | | |

Other tables you will actually use: `interventions` (+ `intervention_other_names`),
`sponsors` (`lead_or_collaborator='lead'`), `conditions`, `detailed_descriptions`,
`reported_events`, `result_groups`, `drop_withdrawals`, `outcome_measurements`.

### `overall_status` — real distribution, and `why_stopped` coverage inside it

```
overall_status              n        with why_stopped
Completed                   325,239        0
Unknown                      94,294        0
Recruiting                   65,410        0
Terminated                   34,011   30,489  (89.6%)
Not Yet Recruiting           29,087        0
Active Not Recruiting        21,968        0
Withdrawn                    16,588   14,476  (87.3%)
Enrolling By Invitation       5,248        0
Suspended                     1,748    1,591  (91.0%)
Withheld / No Longer Available / Available / Approved For Marketing / …
NOT_YET_RECRUITING (7), RECRUITING (6), ACTIVE_NOT_RECRUITING (4),
WITHDRAWN (1), AVAILABLE (1)          <-- upstream SCREAMING_CASE stragglers
```

Two consequences you must encode in every query:

1. `why_stopped` is populated **only** for Terminated / Withdrawn / Suspended.
   Never filter on `why_stopped IS NOT NULL` to *find* stopped programs if you
   also care about programs that quietly completed and were then shelved — those
   are `Completed` with a dead pipeline, and ctgov will never tell you.
2. `overall_status` is **Title Case with SCREAMING_CASE stragglers**. Use
   `overall_status ILIKE 'terminated'` or `~* '^terminat'`, never `= 'Terminated'`.

### What `why_stopped` actually contains

Free text, sponsor-authored, min 1 char, **mean 56, max 252** (the field is
truncated at ~250 by ClinicalTrials.gov). 3,602 entries are under 15 characters.

Most common verbatim values across all 46,557:

```
Sponsor decision (335) · Lack of funding (309) · Slow accrual (271) ·
Lack of enrollment (181) · Sponsor Decision (172) · Low accrual (171) ·
Business decision (141) · slow accrual (141) · low accrual (138) ·
See termination reason in detailed description. (132) · Low enrollment (129) ·
Slow enrollment (120) · No participants enrolled (119) · Poor accrual (110)
```

Bucketed by regex over all 46,557:

| Bucket | n | share |
|---|---|---|
| operational (enrollment / accrual / PI / COVID) | 17,766 | 38% |
| other / uncategorised | 15,260 | 33% |
| business / strategic / funding | 9,252 | 20% |
| **safety / toxicity / AE / hepatic** | **2,480** | **5.3%** |
| efficacy / futility / interim | 1,799 | 3.9% |

**So: `why_stopped` is usable, but it is a 5% signal in a 95% noise field.** It
is the fastest first pass and it is authoritative when it fires — but the
majority of target-relevant terminations are *not* in it, and a target-level
conclusion built only on `why_stopped` will be wrong. Always run the literature
grep too.

`See termination reason in detailed description.` (132 studies, plus variants)
is a pointer, not an answer:

```sql
SELECT s.nct_id, s.why_stopped, LEFT(d.description, 400)
FROM ctgov.studies s JOIN ctgov.detailed_descriptions d ON d.nct_id = s.nct_id
WHERE s.why_stopped ILIKE '%detailed description%'
```

### The workhorse query

```sql
SELECT DISTINCT s.nct_id, i.name, s.overall_status, s.phase, s.enrollment,
       s.start_date, s.completion_date, s.study_first_posted_date,
       s.last_update_posted_date, s.why_stopped, sp.name AS sponsor
FROM ctgov.interventions i
JOIN ctgov.studies s   ON s.nct_id = i.nct_id
LEFT JOIN ctgov.sponsors sp
       ON sp.nct_id = s.nct_id AND sp.lead_or_collaborator = 'lead'
WHERE i.name ILIKE '%<CODE>%'
ORDER BY s.study_first_posted_date
```

Also union in `intervention_other_names.name` — a compound registered as
`SAR441566` in 2022 appears as `balinatunfib` in 2025, and title-only search
misses both when the sponsor titles the study by mechanism.

### `reported_events` gives you the actual adverse-event counts

This is the highest-value table nobody uses. It turns "hepatotoxicity" into a
number, with a denominator, per dose arm:

```sql
SELECT r.nct_id, g.title AS arm, r.adverse_event_term,
       r.subjects_affected, r.subjects_at_risk
FROM ctgov.reported_events r
JOIN ctgov.result_groups g ON g.id = r.result_group_id
WHERE r.nct_id = 'NCT04586920'
  AND r.adverse_event_term ~* '(hepat|transaminase|liver|alanine|aspartate)'
  AND r.subjects_affected > 0
ORDER BY r.subjects_affected DESC
```

Returns, for LY3509754: `Hepatitis acute 3/6` in *Part D (Japanese) — 1000 mg
LY3509754 QD*, `Hepatitis acute 1/6`, `Transaminases increased 1/6` and
`Hepatic steatosis 1/6` in *Part D (Japanese) — 400 mg QD*. Four acute-hepatitis
cases, dose-concentrated at 400–1000 mg — which independently confirms the
published account without ever reading the paper.

Filter `subjects_affected > 0`: the table stores an explicit zero row for every
term × every arm, so an unfiltered query returns hundreds of `0/n` rows.

---

## Rule 3: never silently pick the convenient citation

Papers disagree. Report both, with dates, and say which registry facts each is
consistent with. **A finding that a claim is time-dependent is more useful than a
confident pick.**

The trap is not "trust the newer paper". Verified counter-example:

- **PMC10334362** (*J Pharm Anal*, pub 2023-05-20): "the trial was subsequently
  terminated as reversible transaminase elevation was found in several patients"
- **PMC10487560** (*IJMS*, pub **2023-08-28**, i.e. three months **later**):
  "VTP-43742 … is currently being evaluated in a phase III clinical trial"

The later paper is the wrong one. Reviews recycle boilerplate written years
before submission. **Publication date is a weak proxy for information currency.**
The tiebreak is the registry: ctgov shows exactly two VTP-43742 studies, both
Phase 1, no Phase 3 anywhere, latest posted 2018 — so the phase-III sentence is
unsupported by any trial record and the termination account is not.

Same pattern on LY3509754: PMC13149041 lists it as "currently in clinical
development" while NCT04152382 and NCT04586920 have both read `Terminated` since
2022.

Output format for a contested claim:

```yaml
status: CONTESTED
claims:
  - assertion: "..."
    source: PMC12795581
    source_date: 2025-12-27
    registry_consistency: "..."
  - assertion: "..."
    source: PMC13206447
    source_date: 2026-02-11
    registry_consistency: "..."
resolution: NOT_RESOLVED   # and say why
```

---

## Rule 4: patents are unavailable in this deployment

Confirmed, three ways:

```bash
paperclip search -s patents "…"   # -> "Patents sources are not available."
paperclip sql    -s patents "…"   # -> "Patents sources are not available."
paperclip ls     /patents/        # -> "Patents sources are not available."
```

Do not build against them. This is a real gap for this skill specifically:
`Expert Opin Ther Pat` patent-evaluation articles are frequently the *only*
public record of a discontinued backup series (three separate RORγt reviews cite
Gege 2017, "RORγt inhibitors as potential back-ups for the phase II candidate
VTP-43742", and the corpus has the citation but not the text). Record affected
programs as `NOT_FOUND (patents unavailable)`, never as "no backup series".

---

## Procedure

1. **Enumerate candidate programs.** From the precedent-lookup drug list, plus
   `grep` on the target's class name in `/papers/` to catch codes ChEMBL never
   had. Class-review tables are the richest source — one grep on `VTP-43742`
   returned a review sentence naming six halted RORγt programs at once.
2. **Registry pass.** Workhorse query per code (+ other-names). Record
   `nct_id / status / phase / n / dates / sponsor / why_stopped` verbatim.
3. **Escalate the `Terminated` rows with NULL `why_stopped`.** These are the
   interesting ones. Go to `detailed_descriptions`, then `reported_events`, then
   the literature grep.
4. **Literature pass.** Bare `grep "<CODE>" /papers/`, then the bound regex,
   then `grep -n` the winning paper for line numbers and `cat meta.json` for the
   date and DOI.
5. **Cross-check every literature claim against the registry.** Confirmed /
   contradicted / registry-silent. Say which.
6. **Classify** each termination (taxonomy below) and **date** it.
7. **Write NOT_FOUND** where you found nothing. Do not estimate, do not infer
   from a compound's disappearance, do not convert a Withdrawn record into a
   safety story.

### Termination taxonomy — classify, do not summarise

| Class | Bears on target druggability? | Example |
|---|---|---|
| `TARGET_MECHANISM` | **Yes, strongly** | on-target class effect; thymic aberration across RORγt series |
| `SAFETY_CLINICAL` | Yes, until proven off-target | VTP-43742 transaminase elevations |
| `SAFETY_PRECLINICAL` | Yes | TAK-828F preclinical tox + teratogenicity |
| `SAFETY_OFF_TARGET` | **No** — chemotype problem, not target problem | LY3509754 DILI, hypothesised off-target |
| `EFFICACY_FUTILITY` | Yes | sotorasib NCT04933695, `why_stopped = "Futility"` |
| `PK_FORMULATION` | Weakly — molecule problem | MRTX1133 "Formulation challenges" |
| `BUSINESS_STRATEGIC` | **No** | JDQ443; "business objectives have changed" |
| `OPERATIONAL` | **No** | "Lack of enrollment", "PI left NIH" |
| `UNKNOWN` | — | `Terminated`, `why_stopped` NULL, silent literature |

The `SAFETY_OFF_TARGET` vs `TARGET_MECHANISM` distinction is the whole point of
this skill. Getting it backwards writes off a good target or greenlights a bad
one.

---

# FAILURE MODES

Longest section deliberately. Every one of these was hit or nearly hit while
building this skill.

### F1. Semantic search for a drug code returns confident, wholly irrelevant hits

`search -s pmc "LY3509754 IL-17A small molecule hepatotoxicity"` returns five
well-ranked papers, none of which contains the string `LY3509754`. Nothing in the
output signals this. If you had summarised those five you would have written a
paragraph about doxorubicin hepatotoxicity in mice and attributed it to Lilly.
**Symptom:** results are topically perfect and never name the compound.
**Fix:** grep. Then, if you must, use those hits only as background.

### F2. Terminated with NULL `why_stopped` — the most misleading state

TAK-828 NCT02817516: `overall_status = Terminated`, Phase 1, n=24, ran
2016-06-30 → 2016-08-22 (54 days), `why_stopped = NULL`. The registry tells you a
Phase 1 died after eight weeks and refuses to say why. 3,522 Terminated studies
are in this state. Recording "Terminated, reason not stated" and stopping is a
retrieval failure — the reason was in `/papers/PMC8080595` the whole time.
**Fix:** every Terminated-with-NULL row is a mandatory literature escalation.

### F3. The stale-review inversion — newer paper, older facts

Covered above. PMC10487560 (2023-08-28) says Phase 3; PMC10334362 (2023-05-20)
says terminated. Sorting by `pub_date` and taking the top row gives the wrong
answer. **Fix:** date every claim, corroborate against ctgov, report both.

### F4. Named-target false positives — the icotrokinra trap

Grepping IL-17A literature surfaces **icotrokinra / JNJ-77242113** constantly,
often in the same sentence as oral IL-17 small molecules. It is an oral
**peptide** against **IL-23R** (PMC11289455: "a highly potent, selective peptide
targeting the IL-23 receptor"; PMC11279831: "an oral IL-23R antagonist peptide"),
now in four active Phase 2/3 trials. Filing it under IL-17A would put a thriving
program in the wrong target's dossier and make IL-17A look better than it is.
**Fix:** before adding any program, confirm *molecular target* and *modality*
from a primary or mechanistic source, not from co-occurrence in a review's
"emerging oral therapies" paragraph. Adjacent-pathway agents (IL-23, IL-23R,
TYK2, RORγt) contaminate IL-17A searches in both directions.

### F5. `why_stopped` is 95% noise — do not build a target verdict on it

38% operational, 20% business, 5.3% safety. On KRAS, 28 stopped trials of
approved-or-clinical agents yield "Lack of enrollment", "Business Reasons",
"0 participant accrual", "PI left NIH", "Funder Decision". A naive aggregate
("KRAS: 28 terminated trials!") describes academic site logistics, not KRAS.
**Fix:** classify before counting; report the class breakdown, never a raw count.

### F6. `overall_status` string matching

Enum is Title Case (`Terminated`, `Active Not Recruiting`, `Enrolling By
Invitation`) with a handful of upstream SCREAMING_CASE rows (`NOT_YET_RECRUITING`
×7, `RECRUITING` ×6, `WITHDRAWN` ×1 — and that WITHDRAWN row *does* carry a
`why_stopped`). `= 'Terminated'` silently drops rows. Use `ILIKE` / `~*`.

### F7. Duplicate rows from the interventions join

Placebo comparators register as separate interventions: `SAR441566` and
`SAR441566 matching Placebo` both point at NCT06637631, doubling every row.
Radiolabels do it too (`DC-806` + `[14C]-DC-806`). Multi-arm dosing does it
worst: GSK2981278 returns four identical rows for NCT02548052 (0.03% / 0.1% /
0.8% / 4%). **Fix:** `SELECT DISTINCT` on nct_id-level columns, or aggregate
intervention names, before you count programs.

### F8. Withdrawn ≠ failed

16,588 Withdrawn studies, 87% with a reason, and the reasons are overwhelmingly
"no participants enrolled" — the study never started. NCT06061523 (sotorasib) was
withdrawn because "FDA and EMA agreed existing data are appropriate for
approval". That is a *success* wearing a failure's status code. **Fix:** read
every `Withdrawn` reason before classifying; never bucket Withdrawn with
Terminated.

### F9. The trial that was never registered

VTP-43742's Phase 2a in psoriasis — the study that produced the transaminase
signal and killed the program — **is not in ctgov**. Only two Phase 1 records
exist (NCT02555709, NCT03724292). A registry-only sweep concludes "VTP-43742:
two completed Phase 1s, no problems". The termination exists solely in the
literature. **Fix:** absence of a registry record is never evidence of absence of
a trial, particularly pre-2018, ex-US, and for small biotechs.

### F10. The primary trial report is often outside the corpus

The definitive LY3509754 paper (Datta-Mannan et al., *Clin Pharmacol Ther*
2024;115:1152–61, doi 10.1002/cpt.3185, PMID 38294091) is **not in PMC**:

```bash
paperclip sql "SELECT * FROM documents WHERE doi='10.1002/cpt.3185'"   # 0 rows
```

You reach it only through papers that cite it, and you inherit their paraphrase.
**Fix:** cite the citing paper for the claim *and* name the primary by DOI so the
reader can go get it. Do not present a secondary paraphrase as if you read the
primary.

### F11. `paperclip skill trials` fails and the trial VFS is empty

`Error: Failed to fetch skill "trials"`. `/trials/<NCT>/meta.json` has
`total_lines: 0` and no dates; `grep` on `/trials/<NCT>/content.lines` returns
nothing; `paperclip grep … /trials/NCT…/` errors with `Cannot read path`. None of
this means the data is missing — it is all in `sql -s trials`. **Fix:** go
straight to SQL and do not spend calls proving the VFS is empty.

### F12. Patents unavailable, quietly

Both `search -s patents` and `sql -s patents` return a one-line message and exit
0. In a pipeline that is indistinguishable from "no patents found". **Fix:**
assert the string, and mark affected fields `NOT_FOUND (patents unavailable)`.

### F13. Confusing the failure of a *molecule* with the failure of a *target*

MRTX1133 died of "Formulation challenges" (registry) / "high pharmacokinetic
variability and failed to meet thresholds for advancement" (PMC12352898). KRAS
G12D is fine; that molecule was not orally viable. Similarly SPD304 was never
developed because its 3-alkylindole is CYP-metabolised to toxic electrophilic
intermediates — a chemotype liability, not a TNF liability. Conversely,
BMS-986251 stopped on thymic lymphoma, which is on-mechanism for RORγt and
generalises across the class. **Fix:** the taxonomy column above is mandatory,
per program, with a one-line justification.

### F14. Over-reading a single class-review sentence

PMC10334362's sentence naming six halted RORγt programs is enormously useful and
also a secondary source with no dates and no per-program detail. It says
GSK2981278 was "halted or put on hold", while the same paper's Table 4 records
GSK2981278 as **Completed** with a null efficacy result, and ctgov shows all
three GSK2981278 studies `Completed`. **Fix:** treat class-review lists as
*leads*, resolve each one individually against the registry, and note internal
disagreement inside a single paper when it occurs.

### F15. Silent truncation at 250 characters

`MAX(LENGTH(why_stopped)) = 252`. Sponsors with a real explanation get cut
mid-sentence: `"Trial terminated due to business decision, not based on a..."`.
There is no ellipsis flag. **Fix:** if `LENGTH(why_stopped) > 240`, assume
truncation and pull `detailed_descriptions`.

### F16. Counting programs when you mean trials

One program spawns 5–15 NCT records (SAD, MAD, DDI, hepatic-impairment, QT,
Japanese bridging, expansion). Sanofi's balinatunfib has ≥9. Reporting "9 trials"
implies nine shots on goal. **Fix:** `terminated_programs[]` is keyed by
*program* (compound × sponsor); NCT numbers are evidence inside it.

---

# WORKED RESULTS — four targets

Every claim below carries a PMC id, DOI, or NCT number. Everything not retrieved
is written `NOT_FOUND`.

## RORγt / RORC — UniProt P51449

```yaml
target: {uniprot: P51449, symbol: RORC, name: "RORγt"}
terminated_programs:

- program: VTP-43742 (vimirogant)
  sponsor: Vitae Pharmaceuticals Inc., an Allergan affiliate
  modality: small molecule, orthosteric RORγt inverse agonist
  highest_phase_reached: Phase 2a (psoriasis)
  outcome: TERMINATED
  class: SAFETY_CLINICAL
  reason: >
    Reversible transaminase elevations in several patients in the Phase 2a
    psoriasis study, despite a clear efficacy signal (PASI −29% at 700 mg and
    −23% at 350 mg at 4 weeks, with 50–75% reductions in plasma IL-17A/IL-17F).
  evidence:
    - "PMC10334362 L93: 'A phase II study in patients with psoriasis also showed
       a signal of efficacy; however, the trial was subsequently terminated as
       reversible transaminase elevation was found in several patients.'
       J Pharm Anal, 2023-05-20, doi:10.1016/j.jpha.2023.05.009"
    - "PMC10334362 L90 (Table 4): 'VTP-43742 (Vimirogant) … Phase II Psoriasis
       Terminated: due to elevated reverse transaminase' [sic]"
    - "PMC8080595 L44: 'elevations of reversible transaminase were observed and
       prompted the termination of development.'
       Sci Rep, 2021-04-28, doi:10.1038/s41598-021-88492-1"
    - "PMC12657003 L81: AEs 'included headache, flushing, nausea and transaminase
       elevation.' Drugs in Context 2025, doi:10.7573/dic.2025-8-4"
  registry:
    - NCT02555709  Phase 1  Completed  n=74   2015-08-01 → 2016-03-01  why_stopped: NULL
    - NCT03724292  Phase 1  Completed  n=40   2015-08-01 → 2016-03-08  why_stopped: NULL
  registry_note: >
    The Phase 2a psoriasis study is NOT in ClinicalTrials.gov. Only two Phase 1
    records exist and both read Completed with no why_stopped. Registry-only
    retrieval reports this program as clean. See failure mode F9.
  contradicting_source:
    - "PMC10487560 L46 (IJMS, pub 2023-08-28 — LATER than the termination
       source): 'VTP-43742 … is currently being evaluated in a phase III clinical
       trial for the treatment of plaque psoriasis.' doi:10.3390/ijms241713313
       No Phase 3 VTP-43742 record exists in ctgov. Treated as stale review text;
       NOT used to overturn the termination finding, and reported here per Rule 3."

- program: TAK-828F (TAK-828)
  sponsor: Takeda
  modality: small molecule, tetrahydronaphthyridine RORγt inverse agonist
  highest_phase_reached: Phase 1
  outcome: TERMINATED
  class: SAFETY_PRECLINICAL
  reason: >
    Discontinued on the basis of preclinical toxicology and clinical
    teratogenicity study results.
  evidence:
    - "PMC8080595 L44: 'the phase I trial of oral RORγ antagonist TAK-828F was
       discontinued based on results of preclinical toxicology and clinical
       teratogenicity studies.' Sci Rep, 2021-04-28,
       doi:10.1038/s41598-021-88492-1"
    - "PMC8080595 abstract: 'oral RORγ antagonists (VTP43742, TAK828) with high
       systemic exposure showed toxicity in phase I/II clinical trials and
       terminated development.'"
  registry:
    - NCT02817516  Phase 1  TERMINATED  n=24  2016-06-30 → 2016-08-22  why_stopped: NULL
    - NCT02706834  Phase 1  Completed   n=36  2016-03-01 → 2016-06-17  why_stopped: NULL
  registry_note: >
    Textbook F2: Terminated after 54 days with why_stopped NULL. The registry
    knows it died and will not say why; the literature says why.

- program: GSK2981278
  sponsor: GlaxoSmithKline
  modality: small molecule RORγ inverse agonist, TOPICAL ointment
  highest_phase_reached: Phase 1/2 (topical, plaque psoriasis)
  outcome: DISCONTINUED_AFTER_COMPLETION
  class: EFFICACY_FUTILITY
  reason: >
    Not a safety termination. All registered studies read Completed; the program
    stopped because topical GSK2981278 did not improve psoriatic lesions at any
    concentration tested.
  evidence:
    - "PMC10334362 L93: 'The results showed that topical use of 0.03%, 0.1%, 0.8%,
       and 4% of GSK2981278 ointment did not improve psoriatic lesions.'"
    - "PMC10334362 L90 (Table 4): 'Completed: study of safety and efficacy …
       psoriasis did not improve … scores were unchanged in most test areas
       (≥87%) at day 19.'"
    - "PMC10334362 L90: listed among agents that 'were halted or put on hold to
       allow for more research' — internally inconsistent with the same paper's
       Table 4 and with ctgov. Recorded, not resolved (F14)."
    - "PMC7080699 L52 (Front Immunol, 2020-03-12, doi:10.3389/fimmu.2020.00348):
       GSK-2981278 among agents 'either discontinued or suspended for further
       development'"
  registry:
    - NCT02548052  Phase 1  Completed  n=15  2015-10-22 → 2016-02-19  why_stopped: NULL
    - NCT03004846  Phase 1  Completed  n=8   2017-02-13 → 2017-05-05  why_stopped: NULL
  discontinuation_reason_status: >
    RECOVERED. A prior structured-database sweep could not find a discontinuation
    reason because there is no `why_stopped` to find — the studies completed. The
    reason is lack of topical efficacy, and it was recovered by grepping
    "VTP-43742", not "GSK2981278": the class-review sentence at PMC10334362 L90
    named the whole halted cohort at once ("GSK2981278, PF-06763809,
    JNJ-61803534, VTP-43742 (Vimirogant), TAK-828F, and AZD0284, were halted or
    put on hold to allow for more research"), and L90/L93 supplied the per-agent
    detail. Grepping a sibling compound to find a class table is a reusable move.

class_wide_findings:   # same target, other sponsors — bears on the target itself
- "JNJ-61803534 (Janssen): Phase 1 Terminated 'due to rabbit embryo toxicity'
   (PMC10334362 L90 Table 4); 'terminated due to the toxicity of the compound in
   rabbit embryos' (L93). class: SAFETY_PRECLINICAL"
- "BMS-986251 (BMS): 'development … was terminated as thymic lymphoma was
   observed' (PMC10334362 L94). class: TARGET_MECHANISM — thymic aberration
   recurs across the RORγt series and 50% of embryonic RORγ-deficient mice
   develop T-cell lymphoma (PMC8080595 L44)."
- "AZD0284 (AstraZeneca): Phase 1 'Terminated: under preclinical evaluation'
   (PMC10334362 L90 Table 4); 'discontinued in 2019' (L93)."
- "ABBV-157 / cedirogant (AbbVie): 'Phase Ⅱ trials … were terminated in 2022'
   (PMC10334362 L92); Table 4 records 'Terminated: due to unknown reason'.
   class: UNKNOWN"
- "JTE-451 / retezorogant (Japan Tobacco): oral development 'terminated in
   February 2021', redirected to topical Phase 1 (PMC10334362 L91).
   class: NOT_FOUND (reason not stated)"
- "PF-06763809 (Pfizer): topical, completed, no efficacy — 'did not reduce skin
   infiltration thickness or disease biomarkers' (PMC10334362 L90).
   class: EFFICACY_FUTILITY"

target_level_read: >
  Systemic RORγt inhibition has failed on toxicity across at least four sponsors
  and three distinct mechanisms (hepatic transaminase elevation, embryo/
  teratogenic toxicity, thymic lymphoma). Two independent sources attribute this
  to high systemic exposure and RORγ1/RORγt cross-reactivity, since compounds are
  screened against the shared LBD and are therefore pan-RORγ (PMC8080595 L44).
  Topical routes fail on efficacy instead. This is a TARGET_MECHANISM pattern,
  not a chemotype pattern.
```

## IL-17A — UniProt Q16552

```yaml
target: {uniprot: Q16552, symbol: IL17A}
terminated_programs:

- program: LY3509754
  sponsor: Eli Lilly and Company
  modality: oral small molecule, IL-17A/IL-17RA PPI inhibitor
  highest_phase_reached: Phase 1
  outcome: TERMINATED
  class: SAFETY_OFF_TARGET     # note: off-target, per the sources
  reason: >
    Drug-induced liver injury. Acute/lymphocytic hepatitis and raised
    transaminases in high-dose cohorts, despite strong target engagement. The
    hepatotoxicity is attributed by two independent reviews to an off-target
    effect unrelated to IL-17 biology.
  evidence:
    - "NCT04152382  Phase 1  TERMINATED  n=30   2019-11-20 → 2022-02-09
       why_stopped: 'Terminated due to liver findings'"
    - "NCT04586920  Phase 1  TERMINATED  n=104  2020-10-20 → 2022-10-14
       why_stopped: 'Terminated due to safety findings'"
    - "ctgov.reported_events, NCT04586920: 'Hepatitis acute' 3/6 in arm
       'Part D (Japanese) — 1000 mg LY3509754 QD'; 'Hepatitis acute' 1/6,
       'Transaminases increased' 1/6, 'Hepatic steatosis' 1/6 in
       'Part D (Japanese) — 400 mg LY3509754 QD'. Four acute-hepatitis cases,
       dose-concentrated at 400–1000 mg."
    - "PMC12433869 L62: 'the oral IL-17A small-molecule inhibitor LY3509754
       significantly reduced IL-17A activity but induced lymphocytic hepatitis
       and drug-induced liver injury (DILI) in high-dose cohorts (400–1,000 mg)'"
    - "PMC12182492 L19: 'the development of LY3509754 … was terminated in 2024
       due to hepatotoxicity that was hypothesized to be due to an off-target
       effect'"
    - "PMC12829784 L103: 'The phase 1 trial of LY3509754 in patients with
       psoriasis (NCT04152382) was discontinued due to safety concerns that are
       not related to IL-17 biology'"
    - "PMC10501871 L17: 'LY3509754 has been discontinued due to the apparent
       hepatotoxicity induction'"
  primary_source_not_in_corpus: >
    Datta-Mannan A, Regev A, Coutant DE, et al. 'Safety, tolerability, and
    pharmacokinetics of an oral small molecule inhibitor of IL-17A (LY3509754):
    a phase I randomized placebo-controlled study.' Clin Pharmacol Ther
    2024;115(5):1152–61. doi:10.1002/cpt.3185, PMID 38294091.
    Verified absent from the corpus (`SELECT … WHERE doi='10.1002/cpt.3185'`
    returns 0 rows across PMC + bioRxiv + arXiv). All claims above are from
    citing papers plus the registry results tables, never from the primary. F10.
  contradicting_source:
    - "PMC13149041: lists LY3509754 among agents 'currently in clinical
       development'. Both NCTs have read Terminated since 2022. Stale."
  significance: >
    The registry AE counts and the literature agree independently, and both
    sources that comment on causality call the DILI off-target. This weakens
    LY3509754 as evidence against IL-17A as a target — it is chemotype evidence.

- program: DC-806 (DICE Therapeutics, now Eli Lilly)
  outcome: NOT_TERMINATED — active/positive comparator, listed for contrast
  class: n/a
  status: >
    Oral IL-17A/IL-17RA PPI inhibitor. Reached proof of concept in psoriasis and
    is the first oral small molecule to block the IL-17/receptor interaction with
    demonstrated clinical activity and a favourable safety profile
    (PMC12829784 L103). No terminated or withdrawn record in ctgov.
  registry:
    - NCT05896527  Phase 2  Completed  n=229  2023-05-02 → 2024-03-25  (plaque psoriasis; dose-ranging)
    - NCT05994807  Phase 1  Completed  n=33
    - NCT06045000  Phase 1  Completed  n=8    (incl. [14C]-DC-806 — see F7)
    - NCT06092931  Phase 1  Completed  n=28
    - NCT06808815  Phase 1  Completed  n=104  2021-09-22 → 2022-08-23
  note: >
    ctgov labels NCT05896527 Phase 2; a preprint (bio_5d455b8f77b9) calls it a
    'Phase IIb dose-ranging trial'. Successor LY4100511 (DC-853) is in
    development (PMC13149041). Registry phase labels and literature phase labels
    disagree routinely; quote both.

excluded_false_positive:
  - compound: icotrokinra / JNJ-77242113 / JNJ-2113 / PN-235
    why_excluded: >
      NOT an IL-17A small molecule. It is an oral PEPTIDE antagonist of the
      IL-23 RECEPTOR. 'JNJ-77242113, a highly potent, selective peptide targeting
      the IL-23 receptor' (PMC11289455 title); 'an oral IL-23R antagonist
      peptide' (PMC11279831). Active and expanding — four Phase 2/3 studies
      Recruiting or Active Not Recruiting (NCT06807424, NCT06878404, NCT07196722,
      NCT07196748), zero with why_stopped. Belongs in an IL23R dossier. F4.
```

## TNF — UniProt P01375

```yaml
target: {uniprot: P01375, symbol: TNF}
terminated_programs:

- program: SPD304
  sponsor: Shire (originally); no clinical sponsor
  modality: small molecule, TNF trimer dissociator (first-in-class, He et al. 2005)
  highest_phase_reached: NONE — never entered clinical development
  outcome: NEVER_DEVELOPED
  class: PK_FORMULATION + SAFETY_PRECLINICAL (chemotype)
  reasons_found:
    - "Metabolic liability of the chemotype: 'as the 3-alkylindole moiety of
       SPD304 can be metabolized by cytochrome P450s to produce toxic
       electrophilic intermediates, its further applications in vivo is limited
       (Sun and Yost, 2008).' PMC5893771 L13"
    - "Cytotoxicity precluding in vivo use: 'SPD304 cannot be used in vivo due to
       his [sic] high toxicity.' PMC5469758"
    - "Weak potency: 'IC50 of 22 μM by ELISA' — PMC5893771 L13; described as
       having 'low potency and poor physicochemical properties'
       (bio_11a66dad6f70) and 'general' affinity (PMC10584158)"
    - "Backup series also failed: 'Analogs of SPD304 have been investigated to
       reduce toxicity, but none has produced optimal anti-TNFα activity for
       further clinical studies.' PMC7859422"
    - "Same liability in the RANKL context: 'SPD304 was first reported to promote
       the dissociation of RANKL trimer … but it is interrupted due to the high
       toxicity.' PMC10175743"
  registry: "ZERO records. `SELECT COUNT(*) FROM ctgov.interventions WHERE
    name ILIKE '%SPD304%'` = 0. Consistent with never reaching the clinic."
  interpretation: >
    All four documented reasons are properties of the molecule (3-alkylindole
    bioactivation, µM potency, cytotoxicity), not of TNF. SPD304 is evidence that
    this chemotype failed, and — via the failed analog campaign — weak evidence
    that the trimer-dissociation site is hard to drug from this scaffold. It is
    not evidence that TNF is undruggable by small molecules; balinatunfib
    disproves that.

- program: balinatunfib / SAR441566
  sponsor: Sanofi (with UCB Pharma)
  modality: oral small molecule; stabilises an asymmetric soluble-TNF trimer that
            engages only two of three TNFR1 sites, blocking TNFR1 but not TNFR2
            (PMC11413425, PMC13331996)
  highest_phase_reached: Phase 2
  outcome: CONTESTED — DO NOT RESOLVE
  class: n/a (see below)
  claims:
    - assertion: >
        Sanofi discontinued development of balinatunfib as MONOTHERAPY after a
        Phase 2 failure.
      source: "PMC12795581 L174 (ref 22) — Beaney A, 'Sanofi axes development of
        oral TNF inhibitor as monotherapy following Phase II fail',
        clinicaltrialsarena.com/news/balinatunfib-sanofi-monotherapy-fail/, 2025.
        Citing paper: Artif Intell Life Sci, pub 2025-12-27,
        doi:10.1016/j.ailsci.2025.100143"
      source_date: 2025-12-27 (citing paper); news item dated 2025
      source_type: trade press, cited secondhand — the news item itself is not in
        the corpus and could not be read directly
      registry_consistency: >
        Consistent-ish. The two monotherapy Phase 2s both read Completed, not
        Terminated, with why_stopped NULL: NCT06073119 (psoriasis, n=221,
        completed 2024-12-11) and NCT06073093 (rheumatoid arthritis, n=264,
        completed 2025-07-02). A quiet post-completion shelving is exactly the
        state ctgov cannot express (see F1 of the why_stopped limits).
    - assertion: >
        'Patients with mild-to-moderate psoriasis treated with balinatunfib
        reported no severe or serious TEAEs and showed promising clinical
        responses, suggesting that further evaluation of TNFR1 signal inhibition
        in inflammatory diseases is warranted.'
      source: "PMC13206447 — Nassr N et al., 'Phase 1 study of balinatunfib, an
        oral inhibitor of TNFR1 signal in mild-to-moderate psoriasis', JEADV,
        doi:10.1111/jdv.70262, PMID 41671079"
      source_date: 2026-02-11
      source_type: peer-reviewed primary trial report
      registry_consistency: >
        Consistent. Three balinatunfib studies are live or recently active as of
        the 2026-06 registry snapshot:
          NCT06637631  Phase 2  Recruiting               Crohn's disease,      n=260, start 2024-12-10, last update 2026-06-12
          NCT06867094  Phase 2  Recruiting               ulcerative colitis,   n=204, start 2025-03-28, last update 2026-06-23
          NCT07222189  Phase 2  Enrolling By Invitation  CD + UC long-term ext, n=325, start 2026-05-19, last update 2026-05-22
          NCT07272629  Phase 1  Completed                cardiac repolarisation (TQT), n=48, 2025-12-04 → 2026-06-08
  resolution: NOT_RESOLVED — reported as time- and indication-dependent.
  why_not_resolved: >
    The two claims are not actually in contradiction once indication is added,
    and the registry shows the shape: the PSORIASIS and RHEUMATOID ARTHRITIS
    monotherapy Phase 2s completed and (per the 2025 trade report) the
    monotherapy path was dropped, while the IBD programme (Crohn's, ulcerative
    colitis) is actively recruiting with a long-term extension enrolling as of
    2026-05, and a dedicated TQT study ran to 2026-06. Both sources are accurate
    about different indications at different dates. Any single-sentence verdict
    ('balinatunfib failed' / 'balinatunfib is advancing') is wrong. Reported with
    both dates per Rule 3.
  additional_sources:
    - "PMC12166256 — Nassr N et al., first-in-human SAD/MAD, Clin Pharmacol Ther
       2025, doi:10.1002/cpt.3655: 'demonstrated a good safety profile along with
       favorable PK/PD characteristics', complete TNFα occupancy at all
       timepoints, t½ 22–30 h"
    - "PMC12276028 / PMC12651885 — 'Balinatunfib: A Clinical Oral Small Molecule
       TNFα Inhibitor', ChemMedChem 2025;20:e202500258,
       doi:10.1002/cmdc.202500258"
    - "PMC11734553 — records a China Phase 2, CTR20241078 (Aptuit (Verona) Srl)"

target_level_read: >
  TNF is not a failed small-molecule target; it is a target where the first
  chemotype (SPD304) failed on its own metabolic liabilities and a later,
  mechanistically different chemotype reached Phase 2 in four indications. The
  live question is efficacy sufficiency as monotherapy, not tractability.
```

## KRAS — UniProt P01116 (control)

```yaml
target: {uniprot: P01116, symbol: KRAS}
approved: [sotorasib (AMG 510), adagrasib (MRTX849)]   # KRAS G12C, both approved
sweep_query: |
  SELECT DISTINCT s.nct_id, i.name, s.overall_status, s.phase, s.why_stopped,
         s.start_date, sp.name AS sponsor
  FROM ctgov.interventions i
  JOIN ctgov.studies s ON s.nct_id = i.nct_id
  LEFT JOIN ctgov.sponsors sp
         ON sp.nct_id = s.nct_id AND sp.lead_or_collaborator = 'lead'
  WHERE s.why_stopped IS NOT NULL
    AND i.name ~* '(KRAS|sotorasib|adagrasib|AMG 510|MRTX|divarasib|GDC-6036|JDQ443|opnurasib|glecirasib|olomorasib|LY3537982|D-1553|BI 1701963|BI 1823911|RMC-6236|RMC-6291|garsorasib|JAB-21822|zoldonrasib|daraxonrasib)'
  ORDER BY s.start_date
  -- regex must be one line; PostgreSQL ~* does not ignore embedded whitespace.
  -- Build the alternation from the precedent-lookup drug list, not from memory.
terminated_programs:

- program: MRTX1133
  sponsor: Mirati Therapeutics (Bristol-Myers Squibb)
  modality: non-covalent KRAS-G12D-selective inhibitor
  highest_phase_reached: Phase 1/2
  outcome: TERMINATED
  class: PK_FORMULATION
  reason: >
    Poor drug-like properties, not target failure.
  evidence:
    - "NCT05737706  Phase 1  TERMINATED  start 2023-03-06
       why_stopped: 'Formulation challenges'"
    - "PMC12352898 L47: 'clinical evaluation of MRTX1133 (NCT0537706) was
       recently terminated because the drug exhibited high pharmacokinetic
       variability and failed to meet thresholds for advancement.'"
    - "PMC12462266: 'Although MRTX1133 was a milestone non-covalent KRAS G12D
       targeting inhibitor, its clinical development has been terminated.'"
  significance: >
    KRAS G12D remains an active target (HRS-4642, zoldonrasib, daraxonrasib).
    This is the cleanest available example of F13 — a molecule died, the target
    did not.

- program: JDQ443 / opnurasib
  sponsor: Novartis
  modality: covalent KRAS G12C 'OFF-state' inhibitor
  highest_phase_reached: Phase 3 (KontRASt-02, NCT05132075)
  outcome: DISCONTINUED
  class: BUSINESS_STRATEGIC
  reason: >
    Competitive/pipeline decision in a crowded G12C field, not a safety or
    primary-efficacy failure. ORR 57.1% in KontRASt-01; Grade 3 AEs 5.9%, no
    Grade 4–5 (PMC12529572).
  evidence:
    - "PMC12660860 L106: KontRASt-02 (NCT05132075) phase III vs docetaxel was
       'halted when the manufacturer discontinued the development of opnurasib'"
    - "PMC12658495 (table): 'JDQ443 (opnurasib) | Novartis | Covalent OFF
       inhibitor | Discontinued (2024) | Pipeline termination'"
    - "PMC13308734: 'in an increasingly crowded therapeutic space. As a result,
       the further development of JDQ443 was recently announced to be
       discontinued.'"
    - "PMC12399206: attributes discontinuation partly to 'extremely poor brain
       penetration'"
    - "NCT05999357  Phase 2  WITHDRAWN
       why_stopped: 'Novartis discontinued develpment of JDQ443' [sic — sponsor
       typo preserved verbatim; note this is an investigator-sponsored study
       reporting a pharma decision]"
  contradiction_note: >
    Three sources say pipeline/competitive; one (PMC12399206) says poor brain
    penetration. Both recorded; not resolved.

- program: BI 1701963 (SOS1::KRAS inhibitor, partner agent)
  sponsor: Boehringer Ingelheim
  outcome: TERMINATED (three trials)
  class: BUSINESS_STRATEGIC / UNKNOWN
  evidence:
    - "NCT04627142  Phase 1  Terminated  2020-11-23  why_stopped: 'Sponsor decision'"
    - "NCT04835714  Phase 1  Terminated  2021-04-20  why_stopped: 'Sponsor decision'"
    - "NCT04975256  Phase 1  Terminated  2021-07-28  why_stopped: 'The decision
       was made to terminate this study to further ...' [TRUNCATED at 250 chars —
       F15; detailed_descriptions not resolved]"
  note: "SOS1, not KRAS itself. Included because it appears in any KRAS
    intervention sweep and would otherwise be miscounted against KRAS."

- program: sotorasib, NCT04933695 (Amgen)
  outcome: TERMINATED
  class: EFFICACY_FUTILITY
  evidence: "NCT04933695  Phase 2  Terminated  start 2022-01-28
    why_stopped: 'Futility'  — the single genuinely efficacy-driven stop among
    28 stopped KRAS-agent trials"

control_finding — THIS IS THE POINT OF THE CONTROL: >
  A why_stopped sweep over KRAS-agent interventions (sotorasib, adagrasib/MRTX849,
  MRTX0902, MRTX1133, JDQ443, D-1553, BI 1701963, AMG 510, KRAS TCR products)
  returns 32 intervention rows / 28 distinct NCTs (F7 in miniature).
  Classified by distinct NCT:
    BUSINESS_STRATEGIC   15  'Sponsor decision' ×3, 'Business Reasons',
                             'Business decision', 'business objectives have
                             changed' ×3, 'Funder Decision', 'Sponsor's decision
                             not related to any safety concern', 'This was a
                             strategic business decision. There were no saf...',
                             'Adjustment of drug development strategy.',
                             'Novartis discontinued develpment of JDQ443',
                             'Sponsor discontinued the trial', 'Study drug was
                             discontinued by manufacturer for business ...'
    OPERATIONAL           8  'Lack of enrollment', '0 participant accrual',
                             'No participants enrolled', 'no participants
                             enrolled across sites', 'Grantor withdrew support
                             due to slow accrual', 'Study was closed due to lack
                             of accrual and study outcome...', 'PI left NIH.'
    PK_FORMULATION        1  MRTX1133, 'Formulation challenges'
    EFFICACY_FUTILITY     1  NCT04933695, 'Futility'
    ADMINISTRATIVE        2  NCT05815186 (record merged into NCT05815173);
                             NCT06061523 withdrawn because 'FDA and EMA agreed
                             existing data are appropriate for appr...' —
                             a SUCCESS wearing a Withdrawn status (F8)
    OTHER                 1  NCT06068153, 'In view of the recent developments in
                             the field (the incr...' [truncated]
    SAFETY                0  — and two entries go out of their way to say so
                             ('not related to any safety concern', 'There were no
                             saf[ety]...')
  Zero safety terminations across an approved, heavily-trialled target. A raw
  count of 28 "terminated KRAS trials" would be the single most misleading number
  in this dossier. Always report the class breakdown, never the count. F5.
```

---

# NOT_FOUND ledger

Stated as NOT_FOUND, not estimated:

- **VTP-43742 Phase 2a registry record** — NOT_FOUND. The trial that terminated
  the program has no ClinicalTrials.gov entry. Termination evidence is
  literature-only. (F9)
- **VTP-43742 / TAK-828F backup series fate** — NOT_FOUND (patents unavailable).
  Gege C, *Expert Opin Ther Pat* 2017;27:1–8 (doi:10.1080/13543776.2017.1262350,
  PMID 27852111), cited by ≥5 corpus papers as the record of the VTP-43742
  back-up compounds, is not retrievable: the patents source is disabled and the
  article body is not in PMC.
- **TAK-828F `why_stopped`** — NOT_FOUND in ctgov (NULL on the Terminated
  NCT02817516). Recovered from literature only.
- **TAK-828F termination date** — NOT_FOUND. Trial completion 2016-08-22 is the
  latest firm date; no source dates the program decision.
- **JTE-451 oral discontinuation reason** — NOT_FOUND. Date given (February 2021,
  PMC10334362 L91); reason not stated in any retrieved source.
- **ABBV-157 / cedirogant Phase 2 termination reason** — NOT_FOUND.
  PMC10334362 Table 4 records it verbatim as "Terminated: due to unknown reason".
- **LY3509754 primary trial report** — NOT IN CORPUS (doi:10.1002/cpt.3185,
  0 rows). All detail is secondhand or from ctgov results tables. (F10)
- **Balinatunfib Phase 2 monotherapy failure, primary source** — NOT_FOUND. The
  only evidence is a trade-press URL appearing in a reference list
  (PMC12795581 L174). The news item is not in the corpus; no peer-reviewed report
  of the Phase 2 psoriasis or RA topline result was retrieved. Neither
  NCT06073119 nor NCT06073093 has posted results.
- **Whether the balinatunfib monotherapy decision extends to IBD** — NOT_FOUND
  and NOT_RESOLVED. The IBD trials were recruiting as of the 2026-06 snapshot;
  no source connects the two decisions.
- **SPD304 formal development-stop record** — NOT_FOUND. No corporate or
  registry record exists; it appears never to have had one. The reasons in the
  block above are the literature's retrospective account, not a sponsor
  statement.
- **BI 1701963 NCT04975256 full termination reason** — NOT_FOUND (truncated at
  250 chars in `why_stopped`; `detailed_descriptions` not resolved). (F15)
- **Patent-only programs, all four targets** — NOT_FOUND (patents unavailable in
  this deployment). (F12)

---

# Output contract

```yaml
terminated_programs:
- program: str                 # compound code + INN if it has one
  sponsor: str
  modality: str                # small molecule / peptide / mAb / other
  target_confirmed: bool       # F4 gate — did you verify the molecular target?
  highest_phase_reached: str
  outcome: TERMINATED | DISCONTINUED_AFTER_COMPLETION | NEVER_DEVELOPED
         | CONTESTED | NOT_TERMINATED
  class: <taxonomy>            # one of the nine
  reason: str | NOT_FOUND
  reason_date: date | NOT_FOUND
  evidence: [str]              # each entry carries PMC+line, DOI, or NCT
  registry: [str]              # nct_id / status / phase / n / dates / why_stopped
  registry_note: str | null    # where registry and literature diverge
  contradicting_source: [str] | null   # with dates; never silently dropped
  resolution: RESOLVED | NOT_RESOLVED
not_found: [str]               # explicit, never estimated
patents_checked: false         # always — patents unavailable
```

Hard rules for the writer:

- Every claim carries a source. A claim without PMC+line, DOI, or NCT does not go
  in the block.
- `NOT_FOUND` is a valid, expected, and often correct value. An estimate is not.
- Contradictions are reported, dated, and left unresolved unless the registry
  settles them. Never silently pick.
- `class` is mandatory. "Terminated for safety" without the on-target /
  off-target call is not an answer.
- Never report a bare count of terminated trials.

## Mapping into the dossier JSON

The dossier template's `target_precedent.terminated_programs[]` is narrower:

```json
{"program": "", "year": null, "stated_reason": "", "source": ""}
```

Collapse each block above into it as follows, and keep the full block as the
working record behind it:

| dossier key | from |
|---|---|
| `program` | `program` (code + INN) |
| `year` | `reason_date` year, else the year of the terminating trial's `completion_date`; `null` if `NOT_FOUND` — never the paper's publication year |
| `stated_reason` | `class` + `:` + `reason`, e.g. `"SAFETY_CLINICAL: reversible transaminase elevations in Phase 2a psoriasis"`. Where `outcome: CONTESTED`, prefix `"CONTESTED (see note): "` and give both claims with their dates. |
| `source` | the single strongest identifier — NCT if `why_stopped` carried it, else `PMC<id> L<n>`, else DOI |

Also emit, per the parent contract:

- programs excluded as false positives (F4) → `not_found[]` with the reason, so
  the reader can see icotrokinra was considered and rejected rather than missed;
- every `NOT_FOUND` from the ledger → `not_found[]` verbatim;
- `falsification.checks_run` gets `"terminated-programs sweep: ctgov why_stopped
  + literature grep"` **whether or not it found anything** — a target with no
  terminated programs is a real finding and silence is not.

**This skill never changes a number.** Per dossier rule 7, clinical failure is
not evidence against tractability: RORγt is simultaneously the most convincing
small-molecule-tractable target in this fixture set and the most comprehensively
clinically failed. Record both. Do not let a termination lower
`tractability.*`, and do not let a clean sweep raise it.
