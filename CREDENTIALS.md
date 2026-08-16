# Credentials: where they live, what breaks, how to rotate

Short answer to the question that prompted this file: **the deployed agent holds
no API key, and it cannot be given one by this deployment.** Credentials stay on
the operator's machine. What follows is why that is structurally true rather
than merely currently true, and what to do when a key changes.

> **Setting this up on a machine that has never run it? Read [§8](#8-onboarding-a-second-operator--the-whole-path-from-a-clean-machine) first**,
> then come back here. Sections 1–7 assume the pipeline already works and tell
> you how it is wired and how to rotate keys; §8 is the path from nothing.
> It also states plainly which parts **cannot** be shared, and what breaks
> without them, so you find that out at the start rather than forty minutes
> into a run.

## 1. Why the deployed agent holds no key

Every one of this agent's nine custom tools is answered by a handler that runs in
the **local process** — the laptop or server that called `runTask` — not in the
cloud sandbox. `lib/claude-managed-agent.ts` parks the session at
`status_idle` / `stop_reason: requires_action`, this process runs the matching
handler, and posts the result back as `user.custom_tool_result`.

What actually uploads is a short list, and no credential is on it:

| Artifact | Uploaded as | Contains a key? |
| --- | --- | --- |
| `CLAUDE.md` | the agent's `system` prompt | no |
| `rubric.md` | runtime rubric, outcome mode (`claude-managed-agent.ts:294`) | no |
| `.claude/skills/<dir>/**` | zipped whole, Skills API | no |
| each tool's `name`, `description`, `input_schema` | agent `tools[]` (`scripts/deploy.ts:166-174`) | no |
| `manifest.json` `name`, `description`, `model` | agent config | no |

**Handler bodies never leave the machine.** `deploy.ts` reads `tool.handler` only
to *call* it locally; it is never serialised into `agentConfig`. `.env` is loaded
by dotenvx into this process and is uploaded by nothing. `acl.ts`, `fixtures/`
and `pipeline.html` do not upload at all.

So the sandbox could not use a Paperclip key if it had one — it has no
`paperclip` binary to use it with. The credential and the binary are on the same
machine, and that machine is not the sandbox.

### No vault is needed, and none should be provisioned

`manifest.json` sets `"mcp_servers": []`. The vault mechanism in the starter
(credentials of type `static_bearer`, keyed by MCP server URL, attached via
`vault_ids`) exists **only** to let a deployed agent authenticate to a remote MCP
server. This agent talks to no MCP server.

**Do not provision a vault credential for this agent.** Doing so would place a
real key in cloud-side storage to serve a code path that does not execute, which
is strictly worse than the current position. If this agent ever gains an
`mcp_servers` entry, revisit this section — until then the correct number of
vault credentials is zero.

## 2. Credential and binary inventory

Two of these are `.env` variables. The third is not, and that catches people out:
**Modal does not authenticate from `.env`.**

| What | Where it lives | What breaks without it |
| --- | --- | --- |
| `ANTHROPIC_API_KEY` | repo-root `.env` | Nothing runs. No session is created at all. |
| `PAPERCLIP_API_KEY` | repo-root `.env` | The entire retrieved-precedent axis: `paperclip_sql`, `_search`, `_grep`, `_read`, and the Paperclip half of `neighbour_precedent`. |
| Modal token (`token_id` / `token_secret`) | `~/.modal.toml`, under `[rafwiewiora]` | The entire computed-tractability axis — `pocket_scan`, which is the only route to fpocket/mdpocket. |
| `MODAL_BIN` | env var (see §3) | Same as above; there is no local fallback. |
| `MICROMAMBA_BIN` | env var, defaults to `~/.local/bin/micromamba` | `cryptic_analysis`, `interface_analysis`, `disorder_scan`, `neighbour_precedent`. |
| `PAPERCLIP_BIN` | env var, defaults to `~/.local/bin/paperclip` | Same as `PAPERCLIP_API_KEY`. |
| `DRUGGABILITY_ENV` | env var, defaults to `druggability` | The gemmi/numpy scripts. |
| `MODAL_PROFILE` | env var, defaults to `rafwiewiora` | Wrong workspace, or none. |
| Membership of the `rafwiewiora` Modal **workspace** | granted by its owner, not by any file in this repo | A Modal token you can mint yourself does not help: it authenticates you to *your* workspace, and the profile guard rejects it by name. See §8.2. |
| `PAPERCLIP_CONFIG_DIR` | env var, defaults to a per-run temp dir the handler creates | Not a credential — an **isolation** setting. Left pointing at the shared `~/.paperclip`, every SQL query becomes a ~15 ms no-op that exits 0. See §8.4. |

**A fourth thing that is not in that table and catches people out on this
machine specifically:** the `paperclip` CLI also reads a stored login from
`~/.paperclip/credentials.json`. On the laptop this pipeline has always run on,
that file exists, so the CLI keeps working **even with `PAPERCLIP_API_KEY`
unset** — verified by running the preflight with the variable deleted: it
reports the missing variable but its liveness probe still returns rows. The
`requireEnv` guard is therefore what protects a *clean* machine, and the
liveness probe on *this* machine cannot prove a clean machine would work. Do
not read a green preflight here as evidence that a colleague's setup is
complete.

## 3. Making the Modal binary durable

**Done — this section used to describe a live defect and now describes the
fix.** The `modal` binary was previously only present at
`/private/tmp/foldarium-modal-test-venv/bin/modal`, which does not survive a
reboot and is not on PATH, so the next reboot would have turned every
`pocket_scan` into a hard failure. It now lives in the same `druggability`
micromamba env as everything else:

    /Users/bb/micromamba/envs/druggability/bin/modal      # modal 1.5.4

and the repo-root `.env` pins `MODAL_BIN` to that absolute path. `MODAL_BIN` is
honoured ahead of PATH, so setting it is always sufficient; §8.3 builds the env
that contains it from nothing.

Two things about that failure worth keeping, because they generalise: the
breakage was **loud** (`resolveBin` throws naming `MODAL_BIN`) rather than a
silent fallback, and loud is what made it a five-minute fix instead of a
mid-run mystery. And the binary living in the same env as `fpocket`, `gemmi` and
`numpy` means one `micromamba create` reproduces all of it, rather than four
independent installs that can each drift.

### The profile is enforced, not merely defaulted

`~/.modal.toml` contains three profiles: `foldariumtest`,
`molspace-production` and `rafwiewiora` (the active one). Only `rafwiewiora` may
be used by this pipeline.

Two guards now enforce that, because defaulting alone was not enough:

1. A blank `MODAL_PROFILE` is treated as unset. The previous
   `process.env.MODAL_PROFILE ?? "rafwiewiora"` would have let `MODAL_PROFILE=""`
   through — empty string is not nullish — and Modal would then have silently
   selected its own active profile, which is whatever `modal profile activate`
   last set. That is a wrong-workspace bug with no error message.
2. A profile that is not `rafwiewiora` is **rejected by name**, not merely
   checked for existence. Existence was the wrong test: the forbidden workspaces
   are in the same file, so "is it a real profile" waves
   `molspace-production` straight through. If you ever genuinely need another
   workspace, set `MODAL_PROFILE_OVERRIDE` to the same value to acknowledge it.

## 4. Rotation procedure

### Why this is mandatory, not ceremony

Both the Anthropic key and the Paperclip key were **pasted into the session
transcript**. That transcript is the primary input to `/managed-agent-deploy`,
which mines it for lessons and writes `CLAUDE.md`, `rubric.md`, skills and
`manifest.json` from it. So the keys exist in at least one place that is read by
a program whose job is to copy things out of it and upload them.

Nothing in the current artifacts contains a key — that was swept and is clean
(§6). The exposure is the transcript itself, not the output. But a key that has
been pasted into a document processed by a compiler is a key that must be
treated as disclosed, regardless of whether this particular compile happened to
copy it. Rotate both.

### Anthropic key

1. Revoke and re-issue at <https://platform.claude.com/settings/keys>.
2. Update `ANTHROPIC_API_KEY` in the repo-root `.env`. Nowhere else — the
   `.claude/get-api-key.sh` helper reads that same file, so there is exactly one
   copy.
3. Verify:

       bun run console druggability-dossier -- --once "reply with OK"

   A stale key fails immediately with `401` from the Agents API.

### Paperclip key

1. Re-issue through Paperclip.
2. Update `PAPERCLIP_API_KEY` in the repo-root `.env`.
3. Verify — and verify with a query whose answer you already know, not one that
   could legitimately be empty:

       paperclip sql -s proteins "SELECT COUNT(*) FROM chembl_v.drugs_by_accession WHERE accession = 'P23458'"

   JAK1 has 11 approved rows. A zero or an error means the key did not take.

### Modal token

Rotation only — this assumes you are **already a member** of the `rafwiewiora`
workspace. If you are not, no `modal token new` will get you in; see §8.2.

1. `modal token new --profile rafwiewiora` — this rewrites `~/.modal.toml`.
2. Confirm `active = true` still sits under `[rafwiewiora]` and that you have not
   been switched to another workspace.
3. Verify: `"$MODAL_BIN" profile current` should print `rafwiewiora`.

### After any rotation

Run the preflight. It checks all three credential sources and all four binaries
in one pass:

    bun -e 'import {preflight} from "@/managed/druggability-dossier/tools.ts"; await preflight(); console.log("preflight OK")'

## 5. Preflight: failing at second zero

`preflight()` is exported from `tools.ts` and is called by
`agent/tools/druggability-dossier.ts` **before the session is created**. It
verifies `ANTHROPIC_API_KEY`, `PAPERCLIP_API_KEY`, the `paperclip`,
`micromamba` and `modal` binaries, the Modal profile, and that the
`DRUGGABILITY_ENV` conda env actually imports gemmi and numpy.

Two design points that are the whole reason it exists:

- **It aggregates.** Every problem is reported in one throw, not just the first.
  Failing one at a time turns setup into a guess-and-recheck loop.
- **It runs before the run, not at first use.** `pocket_scan` is typically
  reached tens of minutes into a dossier, after the precedent queries. A missing
  `MODAL_BIN` discovered there costs the entire run.

`DOSSIER_SKIP_PREFLIGHT=1` exists for driving the Paperclip tools by hand on a
machine with no Modal. Never set it for a real dossier.

### The failure mode this is all defending against

A missing credential must never look like a negative result. A `paperclip` call
with no valid key returning zero rows is indistinguishable from *a target with no
precedent* — and telling those two apart is this agent's entire job. Every
credential path therefore raises with a named variable rather than returning an
empty set.

That includes a key that is **present but dead**: `requireEnv` cannot catch an
expired key, so failed `paperclip` runs are additionally checked for auth
signatures (401/403/unauthorized/invalid key) and converted into a throw that
says, in words, "this is an authentication failure, NOT an empty result". The
check is scoped to non-zero exits on purpose — a successful literature search can
legitimately return document text containing the word "unauthorized", and turning
retrieved evidence into a hard failure would be the same bug wearing a different
hat.

## 6. Pre-upload artifact scan — PROPOSED THEN APPLIED, and applied better

**This section used to open "There is currently no guard." That is no longer
true.** It is corrected rather than deleted, because the reasoning still matters
and because a document that tells you you are unprotected when you are is how a
second, weaker copy of the same guard gets written.

The guard now exists, in shared code, and it landed in a **stronger** place than
this section originally proposed:

| where | what it covers |
| --- | --- |
| `lib/credential-scan.ts` | the scanner itself — `findCredentials`, `redact`, `assertNoCredentials`. One copy, so the patterns cannot drift. |
| `scripts/deploy.ts:27, 87, 390` | deploy-time artifacts, including the skill bundles, which only deploy can see. |
| `lib/claude-managed-agent.ts` → `loadManagedAgent()` | **the real chokepoint.** Every path converges here — console, the router wrapper, and deploy. |
| `lib/claude-managed-agent.ts:516` | custom-tool *results* are redacted in flight, with a note telling you to rotate. |

**Why the chokepoint matters more than the deploy-time scan**, and this is the
part the original proposal got wrong: `rubric.md` never passes through
`deploy.ts` on its way to the API at all. It ships at *runtime*, from
`runTask()`'s `user.define_outcome` event — and so do several `manifest.json`
fields (session title, `vault_ids`, `memory.instructions`, sent at
`sessions.create`). A key added to either file *after* a successful deploy would
have reached the API on the next `bun run console` without `deploy.ts` ever
seeing it. Scanning inside `loadManagedAgent()` closes that.

Nothing further is owed here. What follows is the original proposal, kept for
the reasoning only.

Scan every artifact that leaves the machine — the skill zips, `instructions`
(CLAUDE.md), `rubric.md`, and the serialised tool descriptions — and **fail the
deploy**, not warn. A warning in a deploy log is a warning nobody reads.

Suggested patch to `scripts/deploy.ts`, before the `agents.create` /
`agents.update` calls:

```ts
const CREDENTIAL_PATTERNS = [
  /sk-ant-[A-Za-z0-9_-]{20,}/,       // Anthropic
  /gxl_[A-Za-z0-9_-]{16,}/,          // Paperclip
  /AKIA[0-9A-Z]{16}/,                // AWS
  /-----BEGIN [A-Z ]*PRIVATE KEY-----/,
];

function scanArtifact(label: string, text: string): void {
  for (const pattern of CREDENTIAL_PATTERNS) {
    if (pattern.test(text)) {
      throw new Error(
        `refusing to deploy: ${label} matches a credential pattern ` +
          `(${pattern}). The compiler mined a transcript containing a live ` +
          "key. Remove it from the artifact AND rotate the key — it is " +
          "disclosed either way."
      );
    }
  }
}

scanArtifact("CLAUDE.md (system prompt)", instructions);
scanArtifact("rubric.md", rubric ?? "");
scanArtifact("tool definitions", JSON.stringify(agentConfig.tools));
// and, per skill bundle, over the zip's text members before upload
```

Two notes for whoever applies it:

- Scan the **zip contents**, not the zip bytes — compression hides the pattern.
  The same pass should exclude `__pycache__`/`*.pyc`, which currently upload and
  carry absolute local paths (already routed in `manifest.json`; this agent's
  `.gitignore` stops them being committed but does not stop the zip).
- A match means rotate, not just edit. Removing the key from the artifact does
  not un-disclose it.

## 7. A trap for anyone re-running the credential sweep

`grep` in this shell is **not** `/usr/bin/grep`. It is a function wrapping
`ugrep` with `--ignore-files`, which honours `.gitignore`. So:

    grep -r 'sk-ant-' .        # silently skips .env, *.pyc, and everything gitignored

That will report clean on a repo whose `.env` is full of live keys, because
`.env` is gitignored. Any credential sweep must use `/usr/bin/grep` directly, or
pass `--no-ignore-files`, and should include `-a` so compiled bytecode is
searched too. A scan that cannot find a key you know is there is not evidence of
absence — verify the scanner against a known positive before trusting a negative.

## 8. Onboarding a second operator — the whole path from a clean machine

**Why this section exists.** Until now the computed axis has run on exactly one
laptop, and everything it needs was acquired incidentally there: a Modal
profile, a Paperclip login, a conda env, and an environment variable that stops
a shared config file from silently voiding every query. None of it was written
down, so the pipeline had a single point of failure and it was a person. This
section is the path a second person follows. Read §8.6 **first** if you only
have five minutes — it says which parts you cannot obtain for yourself.

Nothing here changes the profile guard. `MODAL_PROFILE` really is restricted to
`rafwiewiora` and other profiles really are rejected by name (§3), because the
founder's instruction is that the other workspaces must not be billed or read,
and because a blank value used to fall through to whatever profile Modal had
active. The guard is correct. What was missing is the answer to "so how do *I*
get in", which is §8.2.

### 8.1 What you are actually setting up

The dossier has two axes and they have completely separate dependencies:

| axis | needs | if it is missing |
| --- | --- | --- |
| retrieved precedent | the `paperclip` CLI, an authenticated Paperclip identity, and an unpoisoned config dir | `paperclip_sql` / `_search` / `_grep` / `_read` and half of `neighbour_precedent` raise. They **raise rather than return zero rows**, on purpose — see §5. |
| computed tractability | the `modal` CLI, membership of the `rafwiewiora` Modal workspace, and the `druggability` micromamba env | `pocket_scan` fails. There is no local fallback wired into the tool, so this is the entire axis. |

`ANTHROPIC_API_KEY` sits above both: without it no session is created at all.

### 8.2 Modal — the thing you need is workspace membership, not a token

**The app must live in one workspace, and it is `rafwiewiora`. You cannot use
your own.** This is not a defaulting preference that a flag relaxes:

- `tools.ts` holds `EXPECTED_MODAL_PROFILE = "rafwiewiora"` and tests
  **identity**, not existence — existence was the wrong test, because the
  forbidden workspaces (`molspace-production`, `foldariumtest`) are in the same
  `~/.modal.toml` and would have been waved straight through.
- A blank `MODAL_PROFILE` is treated as unset rather than passed on, so it
  cannot fall through to whatever `modal profile activate` last selected.

So minting yourself a Modal token does not help. A token you can create on your
own authenticates you to *your* workspace, and the guard rejects it by name.
**What you need is to be added to the existing workspace**, which only its
owner can do:

1. Ask the `rafwiewiora` workspace owner to invite you (Modal dashboard →
   workspace **Settings → Members**). This is the step with a human in it, and
   it is the one that used to be undocumented.
2. On your machine, create a **local profile literally named `rafwiewiora`**
   and point it at that workspace:

       modal token new --profile rafwiewiora

   The browser flow asks which workspace to authenticate against — choose the
   shared one. This writes a `[rafwiewiora]` block into `~/.modal.toml`.
3. Confirm you landed in the right place, because step 2 is where it goes
   wrong silently:

       modal profile current      # → rafwiewiora
       modal app list             # → should list the druggability-* apps

   An empty app list means you named a local profile `rafwiewiora` but pointed
   it at your own workspace. **The guard checks the profile name, not the
   workspace behind it** — it defends against picking the wrong profile by
   accident, not against wiring the right name to the wrong place. `modal app
   list` is the check that actually catches this.
4. Leave `MODAL_PROFILE` unset, or set it to `rafwiewiora`. Both work; unset is
   less to get wrong. `MODAL_PROFILE_OVERRIDE` exists only to acknowledge a
   deliberate switch and is not part of onboarding.

Your Modal token is personal. **Do not copy `~/.modal.toml` between machines** —
it is a credential file, and the profile guard is about which workspace is
billed, not about which human is at the keyboard.

One practical consequence of the shared workspace worth knowing: the pocket-scan
image is **built and cached there**, so your first `pocket_scan` starts a
container instead of building a ~2 GB image. In a fresh workspace it would build
from scratch, and `POCKET_SCAN_TIMEOUT_S` is 1800 s for the whole `modal run`
subprocess — a build plus a scan can exceed it. That is a second, quieter reason
the workspace is shared rather than per-person.

### 8.3 The `druggability` micromamba env, from nothing

**`micromamba` only. Never `conda`, never `mamba`** — that is a standing rule in
this project and it applies to every command here without exception.

This one env carries the whole local stack: the pocket binaries, the structure
libraries, the Modal client, and the optional Foldseek path.

    /Users/bb/.local/bin/micromamba create -y -n druggability -c conda-forge \
        python=3.14 fpocket gemmi numpy pip

    /Users/bb/.local/bin/micromamba run -n druggability python -m pip install \
        modal \
        "proto-tools[mcp] @ git+https://github.com/evo-design/proto-tools.git@1e0bd8f5a8f4525eb5e5c736cbf25c1366929e73"

What that gives you, and what each piece is for:

| package | source | why it is here |
| --- | --- | --- |
| `fpocket` 4.2.3 | conda-forge | pocket detection. **`mdpocket` ships inside this same package**, along with `tpocket` and `dpocket` — there is no separate `mdpocket` install and looking for one is a dead end. |
| `gemmi` 0.7.5 | conda-forge | mmCIF parsing. It is the *only* structure format read. |
| `numpy` 2.5.2 | conda-forge | Kabsch superposition, `.dx` grid handling, `cryptic_analysis`, `interface_analysis`. |
| `python` 3.14.6, `pip` | conda-forge | the interpreter the analysis scripts run under. |
| `modal` 1.5.4 | pip | the client `pocket_scan` shells out to. Installing it *here* rather than in a throwaway venv is what made `MODAL_BIN` durable (§3). |
| `proto-tools` 0.1.0 | pip, pinned to `1e0bd8f5…` | Foldseek search for `neighbour_precedent`. **Optional** — see below. |

Then point the repo-root `.env` at it:

    MICROMAMBA_BIN=/Users/<you>/.local/bin/micromamba
    MODAL_BIN=/Users/<you>/micromamba/envs/druggability/bin/modal
    DRUGGABILITY_ENV=druggability

Verify the env directly, the same way the preflight does:

    /Users/bb/.local/bin/micromamba run -n druggability python -c "import gemmi, numpy"
    /Users/bb/.local/bin/micromamba run -n druggability which fpocket mdpocket

**`proto_tools` is optional and the preflight treats it as a warning, not a
failure.** That asymmetry is deliberate: CLAUDE.md rule 13 says a missing
`proto_tools` **nulls** the `structural_neighbour_precedent` axis with a stated
reason, which is a legal dossier — so failing the preflight on it would refuse
runs the agent is specified to complete. But discovering it at the point of use
looks identical to "this fold has no neighbours", which is exactly the confusion
this whole pipeline exists to prevent, so it is said out loud at second zero.
The pinned SHA above is the same commit the Modal image runs; keep the two the
same or the local and remote halves of that axis diverge.

### 8.4 Paperclip — the key, and the config-directory isolation

**The key.** `PAPERCLIP_API_KEY` (a `gxl_…` value) goes in the repo-root `.env`
and nowhere else, and it authenticates the CLI non-interactively.

**Issue your own; do not copy the one in this repo's `.env`.** §4 records that
both keys in that file were pasted into a session transcript that is the primary
input to a compiler whose job is to copy things out of it, so they are treated
as **disclosed and pending rotation**. Propagating a key that is already queued
for revocation onto a second machine gives you a setup that breaks the moment
§4 is finally carried out, and widens the disclosure in the meantime.

**The isolation, and why it is not optional.** `~/.paperclip/config.json` is
sticky client state, and it carries a persistent `cli_cwd` that the CLI passes
as the **working directory of every command, `sql` included**. When that
directory is not readable the query never runs:

    cli_cwd = "/papers/"   →  vsh: cd: /papers/: Permission denied
                              printed on STDOUT, exit code 0, ~15 ms

A well-formed no-op that exits successfully. It renders like an empty table, and
an empty table is indistinguishable from a target with no precedent — the one
error this dossier cannot survive. One navigation command, by anything on the
machine, poisons every later query. **This has been re-poisoned twice today by
something outside the session that runs the pipeline**, so treat it as a
recurring condition and not a one-off.

The handler already defends itself: it runs Paperclip against **its own** config
directory (`$TMPDIR/druggability-dossier-paperclip`, mode 0700) and rewrites
`cli_cwd` to `/` before every call. So you get isolation for free — with two
ways to lose it:

- **Do not set `PAPERCLIP_CONFIG_DIR` to `~/.paperclip`.** Leave it unset unless
  you have a specific reason. If you set it, set it to a directory this pipeline
  owns. If the poisoned-`cli_cwd` symptom reappears *through the tools*, that is
  the first thing to check.
- **Anything you run by hand is on the shared config.** `paperclip sql …` typed
  into a terminal does not go through the handler and will hit the poisoned
  `cli_cwd`. Prefix it:

      PAPERCLIP_CONFIG_DIR=$(mktemp -d) paperclip sql -s proteins "SELECT 1 AS ok"

**And the rule that generalises past this one bug: with this tool, a zero exit
is not success.** Exit 0 is also what you get from a schema error (`ERR: sql:
unknown column …`, printed on stdout), from a silently row-capped result set,
and from a display cap that renders 5 of 100 rows while truthfully reporting
100. CLAUDE.md rules 14 and 15 enumerate all seven signatures. Check the output,
not the status code, and reconcile every count against an independently issued
`COUNT`.

### 8.5 The `.env` a second operator ends up with

Repo-root `.env`, on the new machine, with your own values:

    ANTHROPIC_API_KEY=…            # yours, from platform.claude.com
    PAPERCLIP_API_KEY=gxl_…        # yours, issued through Paperclip
    MODAL_BIN=/Users/<you>/micromamba/envs/druggability/bin/modal
    MICROMAMBA_BIN=/Users/<you>/.local/bin/micromamba
    DRUGGABILITY_ENV=druggability
    # MODAL_PROFILE — leave unset; it defaults to rafwiewiora
    # PAPERCLIP_CONFIG_DIR — leave unset; the handler manages its own
    # DOSSIER_SKIP_PREFLIGHT — never set this for a real run

`.env.example` at the repo root documents every one of these inline.

### 8.6 What cannot be shared, stated plainly

| thing | can a second person get it themselves? | what happens without it |
| --- | --- | --- |
| **Membership of the `rafwiewiora` Modal workspace** | **No.** Only the workspace owner can add you. A token for your own workspace is rejected by name and is not a workaround. | **The entire computed-tractability axis is unavailable.** `pocket_scan` is the only route to fpocket/mdpocket in this pipeline and there is no local fallback wired into the tool. The dossier can still be produced from retrieved precedent alone, with the computed axis nulled and the reason recorded — but half of what this agent exists to do is gone. |
| Paperclip corpus access | Probably not on your own — the key is issued through Paperclip against an account with corpus access. Ask whoever administers it. | **The entire retrieved-precedent axis.** The handlers raise instead of returning zero rows, so this fails loudly rather than becoming a false "no precedent found". |
| `ANTHROPIC_API_KEY` | **Yes** — self-service at <https://platform.claude.com/settings/keys>. | No session at all; nothing runs. |
| The `druggability` micromamba env | **Yes** — §8.3 builds it from nothing in two commands. | `cryptic_analysis`, `interface_analysis`, `disorder_scan`, `neighbour_precedent` cannot run; `MODAL_BIN` has nowhere durable to live. |
| `~/.modal.toml` tokens | Yes, once you are a workspace member (§8.2 step 2). | `pocket_scan` cannot authenticate. **Never copy this file between machines.** |

The honest summary: **one of the five is a genuine bottleneck and it is Modal
workspace membership.** Everything else a second person can either obtain
self-service or build from these instructions. If the workspace owner is
unavailable, a colleague can still run the retrieved-precedent half of the
dossier today, and that is worth knowing before spending an afternoon on setup.

### 8.7 Verify, in this order

Cheapest first, so a failure costs seconds rather than a whole run.

1. **Preflight.** It aggregates every credential and binary into one throw, so
   one run tells you everything that is wrong rather than the first thing:

       bun -e 'import {preflight} from "@/managed/druggability-dossier/tools.ts"; await preflight(); console.log("preflight OK")'

   Every message names the variable that fixes it and the whole failure points
   back at this file.

2. **A Paperclip query whose answer you already know** — not one that could
   legitimately be empty:

       paperclip sql -s proteins "SELECT COUNT(*) FROM chembl_v.drugs_by_accession WHERE accession = 'P23458'"

   JAK1 has 11 approved rows. Zero, an error, or a suspiciously instant reply
   means the key did not take or the config is poisoned (§8.4).

3. **A real pocket scan**, because a green preflight only proves the binaries
   resolve. 1TNF with `chains: {"1TNF": ["A","B"]}` is the cheapest known-answer
   case: 7 pockets, volume ~155.9–156.0 Å³, druggability 0.201.

One caution on step 1 for anyone verifying on **this** laptop rather than a
clean one: the preflight passes with `PAPERCLIP_API_KEY` deleted, because the
CLI falls back to the stored login in `~/.paperclip/credentials.json` (§2). A
green preflight here is not evidence that a colleague's machine is complete.
