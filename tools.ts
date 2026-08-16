/**
 * Custom tools for the `druggability-dossier` managed agent.
 *
 * Every handler here runs in *this* process — the laptop/server that called
 * `runTask` — not in the deployed agent's sandbox. `lib/claude-managed-agent.ts`
 * parks the session on `requires_action`, runs the matching handler locally,
 * and posts a `user.custom_tool_result` back. So these handlers get the local
 * PATH, the local micromamba envs, the local Modal credentials, and the `.env`
 * dotenvx already loaded into `process.env`.
 *
 * That is the *only* reason these tools exist. The sandbox has unrestricted
 * outbound network (verified on the `mvp-shared` environment:
 * `config.networking.type === "unrestricted"`), so nothing here bridges
 * connectivity. What it bridges is **binaries and conda packages the sandbox
 * cannot install**: `paperclip` is a private binary in none of the six
 * declarable registries (apt/cargo/gem/go/npm/pip), and `fpocket`/`mdpocket`
 * are conda-forge only with no conda member in the schema.
 *
 * Handlers shell out to the exact scripts under `.claude/skills/`. None of the
 * Python logic is reimplemented here — a reimplementation would be a second,
 * untested copy of a calibrated measurement.
 *
 * Credentials come from `process.env` and every miss is a hard failure with a
 * named variable. A silent fallback would produce a run with no data that reads
 * exactly like a run with no results, which is the worst possible output for
 * this agent.
 */
import { execFile } from "node:child_process";
import {
  existsSync,
  mkdirSync,
  mkdtempSync,
  readFileSync,
  writeFileSync,
} from "node:fs";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { promisify } from "node:util";
import type { CustomToolSpec } from "@/lib/claude-managed-agent.ts";
import { repoRoot } from "@/lib/claude-managed-agent.ts";

const execFileAsync = promisify(execFile);

const SKILLS_DIR = join(
  repoRoot,
  "managed",
  "druggability-dossier",
  ".claude",
  "skills"
);

/** Truncate handler output so one wide result cannot blow the session budget. */
const MAX_OUTPUT_CHARS = 180_000;
const MS_PER_SECOND = 1000;
const DEFAULT_TIMEOUT_S = 300;
const POCKET_SCAN_TIMEOUT_S = 1800;
const NEIGHBOUR_TIMEOUT_S = 900;
const EXEC_MAX_BUFFER = 64 * 1024 * 1024;

// ---------------------------------------------------------------------------
// Environment resolution — loud on every miss, never a silent fallback
// ---------------------------------------------------------------------------

function requireEnv(name: string, why: string): string {
  const value = process.env[name];
  if (!value) {
    throw new Error(
      `${name} is not set in this process's environment. ${why} ` +
        "Add it to the repo-root .env (dotenvx loads that file into the " +
        "process that answers custom tools) and re-run the task. " +
        "If you do not know where this value comes from, " +
        "managed/druggability-dossier/CREDENTIALS.md section 8 is the " +
        "end-to-end setup path and section 8.6 says what breaks without it."
    );
  }
  return value;
}

/** Find `name` on PATH, or undefined. Used so resolveBin can keep its promise. */
function onPath(name: string): string | undefined {
  const parts = (process.env.PATH ?? "").split(":").filter(Boolean);
  return parts
    .map((part) => join(part, name))
    .find((candidate) => existsSync(candidate));
}

/**
 * Resolve an executable: explicit env override first, then a known install
 * location, then PATH. Never returns a name that does not exist — a bare name
 * that is not actually on PATH is thrown here, naming the variable that would
 * fix it, rather than deferring to an ENOENT that names no variable.
 */
function resolveBin(args: {
  candidates: string[];
  envVar: string;
  fallbackOnPath: string;
  why: string;
}): string {
  const override = process.env[args.envVar];
  if (override) {
    if (!existsSync(override)) {
      throw new Error(
        `${args.envVar} points at ${override}, which does not exist. ${args.why} ` +
          "Fix the path, or see managed/druggability-dossier/CREDENTIALS.md " +
          "section 8.3 for building the toolchain from nothing."
      );
    }
    return override;
  }
  const found = args.candidates.find((candidate) => existsSync(candidate));
  if (found) {
    return found;
  }
  const fromPath = onPath(args.fallbackOnPath);
  if (fromPath) {
    return fromPath;
  }
  throw new Error(
    `${args.fallbackOnPath} was not found: not at any known install location ` +
      `and not on this process's PATH. ${args.why} ` +
      `Install it, then set ${args.envVar} to its absolute path and re-run the ` +
      "task. managed/druggability-dossier/CREDENTIALS.md section 8.3 builds " +
      "the whole toolchain (fpocket/mdpocket, gemmi, numpy, modal) into one " +
      "micromamba env from nothing."
  );
}

function micromamba(): string {
  return resolveBin({
    candidates: [join(homeDir(), ".local", "bin", "micromamba")],
    envVar: "MICROMAMBA_BIN",
    fallbackOnPath: "micromamba",
    why: "It runs the gemmi/numpy analysis scripts (cryptic_analysis, interface_analysis, disorder, neighbour_precedent) inside their conda env.",
  });
}

function homeDir(): string {
  return process.env.HOME ?? process.env.USERPROFILE ?? "";
}

/** The env carrying gemmi + numpy for the local structure-analysis scripts. */
function analysisEnvName(): string {
  return process.env.DRUGGABILITY_ENV ?? "druggability";
}

// ---------------------------------------------------------------------------
// Process runner
// ---------------------------------------------------------------------------

type RunResult = { code: number; stderr: string; stdout: string };

async function run(
  file: string,
  argv: string[],
  opts: { env?: NodeJS.ProcessEnv; timeoutSeconds?: number } = {}
): Promise<RunResult> {
  try {
    const { stdout, stderr } = await execFileAsync(file, argv, {
      env: opts.env ?? process.env,
      maxBuffer: EXEC_MAX_BUFFER,
      timeout: (opts.timeoutSeconds ?? DEFAULT_TIMEOUT_S) * MS_PER_SECOND,
    });
    return { code: 0, stderr, stdout };
  } catch (error) {
    const err = error as {
      code?: number | string;
      stderr?: string;
      stdout?: string;
    };
    if (err.code === "ENOENT") {
      throw new Error(
        `executable not found: ${file}. Install it, or point the matching ` +
          "*_BIN environment variable at it.",
        { cause: error }
      );
    }
    return {
      code: typeof err.code === "number" ? err.code : 1,
      stderr: err.stderr ?? String(error),
      stdout: err.stdout ?? "",
    };
  }
}

function clip(text: string): string {
  return text.length > MAX_OUTPUT_CHARS
    ? `${text.slice(0, MAX_OUTPUT_CHARS)}\n…[truncated by the local tool handler at ${MAX_OUTPUT_CHARS} characters — narrow the query]`
    : text;
}

/** Render a completed run for the agent: stdout on success, both on failure. */
function report(label: string, result: RunResult): string {
  if (result.code === 0) {
    return clip(result.stdout || result.stderr || "(no output)");
  }
  return clip(
    `${label} exited ${result.code}.\n\n--- stdout ---\n${result.stdout}\n\n--- stderr ---\n${result.stderr}`
  );
}

// ---------------------------------------------------------------------------
// Input coercion (the model sends JSON; keep every read explicit)
// ---------------------------------------------------------------------------

function str(input: Record<string, unknown>, key: string): string | undefined {
  const value = input[key];
  return typeof value === "string" && value.length > 0 ? value : undefined;
}

function requiredStr(input: Record<string, unknown>, key: string): string {
  const value = str(input, key);
  if (!value) {
    throw new Error(`"${key}" is required and must be a non-empty string`);
  }
  return value;
}

function num(input: Record<string, unknown>, key: string): number | undefined {
  const value = input[key];
  return typeof value === "number" ? value : undefined;
}

function bool(input: Record<string, unknown>, key: string): boolean {
  return input[key] === true;
}

/**
 * Tri-state read for a switch whose *default is on*. `bool` cannot express the
 * three stage switches (`run_disorder`/`run_cryptic`/`run_mdpocket`): they
 * default to true in `modal_app.py`, so "absent" and "false" must not collapse.
 * Only an explicit `false` may emit a `--no-…` flag.
 */
function optBool(
  input: Record<string, unknown>,
  key: string
): boolean | undefined {
  const value = input[key];
  return typeof value === "boolean" ? value : undefined;
}

function list(input: Record<string, unknown>, key: string): string[] {
  const value = input[key];
  if (Array.isArray(value)) {
    return value.map(String).filter((item) => item.length > 0);
  }
  return typeof value === "string" && value.length > 0
    ? value
        .split(",")
        .map((item) => item.trim())
        .filter((item) => item.length > 0)
    : [];
}

const NUM_SEPARATOR = /[\s,]+/;

/** Residue numbers, accepted as an array or as a comma/space-separated string. */
function numList(input: Record<string, unknown>, key: string): number[] {
  const value = input[key];
  let raw: unknown[] = [];
  if (Array.isArray(value)) {
    raw = value;
  } else if (typeof value === "string") {
    raw = value.split(NUM_SEPARATOR);
  }
  return raw
    .map((item) => Number(item))
    .filter((item) => Number.isFinite(item));
}

/**
 * Render the `chains` mapping into the CLI's own encoding.
 *
 * `modal_app.py`'s local entrypoint takes one string and splits it on `;` then
 * `=`: `{"1TNF": ["A","B"], "2ZJC": ["A"]}` becomes `1TNF=A,B;2ZJC=A`. The
 * model sends the object because that is the shape `pocket_scan()` itself
 * documents; a pre-formatted string is accepted too so a caller reading the
 * SKILL.md's shell example can paste it straight through.
 */
function chainSpec(
  input: Record<string, unknown>,
  key: string
): string | undefined {
  const value = input[key];
  if (typeof value === "string") {
    return value.trim().length > 0 ? value.trim() : undefined;
  }
  if (!value || typeof value !== "object" || Array.isArray(value)) {
    return;
  }
  const parts = Object.entries(value as Record<string, unknown>)
    .map(([pdbId, picked]) => {
      const ids = Array.isArray(picked)
        ? picked.map(String)
        : String(picked).split(",");
      const clean = ids
        .map((item) => item.trim())
        .filter((item) => item.length > 0);
      return clean.length > 0 ? `${pdbId.trim()}=${clean.join(",")}` : "";
    })
    .filter((part) => part.length > 0);
  return parts.length > 0 ? parts.join(";") : undefined;
}

/** Append `flag value` only when there is a value. Keeps argv builders flat. */
function pushFlag(argv: string[], flag: string, value: string | undefined) {
  if (value) {
    argv.push(flag, value);
  }
}

// ---------------------------------------------------------------------------
// paperclip
// ---------------------------------------------------------------------------

function paperclipBin(): string {
  return resolveBin({
    candidates: [join(homeDir(), ".local", "bin", "paperclip")],
    envVar: "PAPERCLIP_BIN",
    fallbackOnPath: "paperclip",
    why: "It is the only route to the Paperclip corpus, which carries the entire retrieved-precedent axis.",
  });
}

/**
 * Everything Paperclip can say that is a failure rather than an answer.
 *
 * This is the distinction the whole agent exists to protect: "no rows" and
 * "the query died" render almost identically once they reach the model, and one
 * of them is a finding about the target while the other is a finding about the
 * laptop. For this pipeline a failure that reads as an empty result is the
 * worst possible outcome — it is the difference between "no family precedent
 * exists" and "three queries died".
 *
 * The list started as auth-only and that was far too narrow. A full dry run on
 * IL-6 lost 11 of 30 `paperclip_sql` calls across FOUR signatures, three of
 * them undocumented, and the auth guard would not have caught a single one:
 *
 *   [error] Request timed out                    (server-side, ~120 s)
 *   [error] Something went wrong. Please try again.   (~76 s)
 *   vsh: cd: /papers/: Permission denied         (a SQL query, exit 0, 15 ms)
 *   sql: unknown column / relation does not exist    (recoverable)
 *
 * TWO MEASURED FACTS ABOUT HOW THEY ARRIVE, both of which the old guard got
 * wrong. First, they arrive on STDOUT WITH EXIT 0 — measured here on
 * `paperclip sql -s proteins`, which printed `vsh: cd: /papers/: Permission
 * denied` and exited 0 three times in a row. The old guard only inspected runs
 * with a non-zero exit, so it was structurally incapable of seeing them, and
 * `report()` handed the text back as a successful result. Second, the
 * `vsh: cd:` case is not random: `~/.paperclip/config.json` carries sticky
 * client state (`{"cli_cwd": "/papers/"}`) left behind by an earlier navigation
 * command, and every later command — `sql` included — is short-circuited by it.
 *
 * So the guard now reads every run, and on a clean exit it anchors on the FIRST
 * non-empty line. A grep hit over the literature begins with a document path,
 * never with `[error]` or `vsh:`, so retrieved evidence that quotes an error
 * string in its text still comes back as evidence.
 */
type FailureSignature = { label: string; pattern: RegExp; remedy: string };

const PAPERCLIP_FAILURES: FailureSignature[] = [
  {
    label: "AUTHENTICATION FAILURE",
    pattern:
      /\b401\b|\b403\b|unauthori[sz]ed|forbidden|invalid api key|invalid token|authentication failed|not authenticated|expired (?:api )?(?:key|token)|missing api key/i,
    remedy:
      "PAPERCLIP_API_KEY was present but Paperclip rejected it. Rotate or " +
      "re-issue the key, put it in the repo-root .env, and re-run.",
  },
  {
    label: "SERVER-SIDE TIMEOUT",
    pattern: /\[error\]\s*request timed out|statement timeout/i,
    remedy:
      "The query reached Paperclip and was killed there. Do NOT assume the " +
      "query was too broad: the three shapes a dry run saw time out here " +
      "(`SELECT 1 AS ok`, a GROUP BY over an assay table, a 17-id IN list " +
      "over pdb_v.entry_ligands) all answer in 7-31 ms once the CLI's sticky " +
      "cli_cwd is out of the way, so try the same query again before " +
      "narrowing it. If it repeats, narrow to one key per query and loop, or " +
      "record the gap in `not_found` — never as 'no rows'.",
  },
  {
    label: "TRANSIENT SERVER ERROR",
    pattern: /\[error\]\s*something went wrong/i,
    remedy:
      "Paperclip failed internally, typically after ~75 s. Retry once; if it " +
      "repeats, record the gap in `not_found` naming this error.",
  },
  {
    label: "STICKY CLIENT STATE (vsh)",
    pattern: /^vsh:|vsh: cd: .*permission denied/im,
    remedy:
      "The Paperclip CLI cds into a sticky working directory before running " +
      "any command, and the one it was given is not readable, so the command " +
      "never ran at all — note the ~15 ms round trip. This handler already " +
      "runs Paperclip against a config directory of its own to prevent " +
      "exactly this, so seeing it means PAPERCLIP_CONFIG_DIR is set in the " +
      "environment and points somewhere poisoned: unset it, or reset " +
      "`cli_cwd` to `/` in that directory's config.json. It is an operator " +
      "fix, not a query fix. Nothing about the query is wrong and NOTHING " +
      "about the target has been learned.",
  },
  {
    label: "QUERY ERROR",
    pattern:
      /sql:\s*unknown column|relation .* does not exist|column .* does not exist|syntax error at or near|no such table/i,
    remedy:
      "The schema does not have what the query asked for. This one is " +
      "recoverable: fix the column or table name and re-run. It is still not " +
      "an empty result.",
  },
  {
    label: "UNCLASSIFIED PAPERCLIP ERROR",
    pattern: /^\[error\]/im,
    remedy:
      "Paperclip printed an error this handler does not have a signature " +
      "for. Treat it as a failed lookup, not as an absence of data, and quote " +
      "the raw text in `not_found`.",
  },
];

function firstLine(text: string): string {
  return (
    text
      .split("\n")
      .map((line) => line.trim())
      .find((line) => line.length > 0) ?? ""
  );
}

/** How the CLI opens a line it wrote itself, as opposed to corpus text. */
const CLI_ERROR_PREFIX = /^(?:\[error\]|vsh:|sql:|err:|error:)/i;

/**
 * Run Paperclip against a config directory this pipeline owns.
 *
 * `~/.paperclip/config.json` is sticky client state — `{"cli_cwd": "/papers/"}`
 * — and `cli/app.py` passes it as the `cwd` of EVERY command, `sql` included.
 * When that directory is not readable the command never runs: the server
 * answers `vsh: cd: /papers/: Permission denied`, the CLI exits 0, and the
 * whole query is a ~15 ms no-op that renders like an empty table. One earlier
 * navigation command, by this agent or by anything else on the machine, poisons
 * every later query.
 *
 * Measured today, both directions, same key, minutes apart:
 *
 *   shared ~/.paperclip (cli_cwd=/papers/)  every sql call → vsh: cd: …, exit 0
 *   PAPERCLIP_CONFIG_DIR=<empty dir>        SELECT 1 AS ok → 7 ms
 *                                           uniprot_v.proteins by accession → 6 ms
 *                                           GROUP BY assay_type → 31 ms
 *                                           17-id IN-list over entry_ligands → rows
 *
 * — i.e. THREE of the four query shapes the dry run reported as 120-second
 * timeouts answer in milliseconds once the sticky cwd is out of the way. So the
 * handler keeps its own directory and rewrites `cli_cwd` to `/` before every
 * call. `PAPERCLIP_API_KEY` is required anyway and authenticates a fresh
 * directory on its own (verified), so no credential file is copied anywhere.
 */
const PAPERCLIP_CONFIG_DIR = join(tmpdir(), "druggability-dossier-paperclip");
const OWNER_ONLY = 0o700;

function paperclipEnv(): NodeJS.ProcessEnv {
  const dir = process.env.PAPERCLIP_CONFIG_DIR ?? PAPERCLIP_CONFIG_DIR;
  mkdirSync(dir, { mode: OWNER_ONLY, recursive: true });
  writeFileSync(join(dir, "config.json"), JSON.stringify({ cli_cwd: "/" }));
  return { ...process.env, PAPERCLIP_CONFIG_DIR: dir };
}

function paperclipFailure(result: RunResult): FailureSignature | undefined {
  const all = `${result.stderr}\n${result.stdout}`;
  if (result.code !== 0) {
    return PAPERCLIP_FAILURES.find((sig) => sig.pattern.test(all));
  }
  // Exit 0 needs the stricter test, and this is the half the old auth guard
  // was right about: a successful grep over the literature can quote an error
  // string as document text, and turning retrieved evidence into a hard
  // failure would be its own version of the bug this guard exists to prevent.
  // Verified against a grep hit whose first line contains both `[error]
  // Request timed out` and `vsh: cd: /papers/: Permission denied` inside a
  // quoted passage: it must come back as evidence, and it does. So the test is
  // the FIRST non-empty line, and only when that line OPENS the way the CLI
  // opens its own errors — a grep hit opens with a document path.
  const head = firstLine(result.stderr) || firstLine(result.stdout);
  if (!CLI_ERROR_PREFIX.test(head)) {
    return;
  }
  return PAPERCLIP_FAILURES.find((sig) => sig.pattern.test(head));
}

async function paperclip(
  argv: string[],
  timeoutSeconds = DEFAULT_TIMEOUT_S
): Promise<RunResult> {
  requireEnv(
    "PAPERCLIP_API_KEY",
    "The Paperclip CLI authenticates non-interactively with it, and without " +
      "it every query would return an auth error that reads like an empty result."
  );
  const result = await run(paperclipBin(), argv, {
    env: paperclipEnv(),
    timeoutSeconds,
  });
  const failure = paperclipFailure(result);
  if (failure) {
    throw new Error(
      `PAPERCLIP ${failure.label} (exit ${result.code}) — THIS IS NOT AN ` +
        "EMPTY RESULT. Nothing was retrieved and nothing about the target " +
        "has been learned, so it must not be recorded as 'no precedent " +
        "found', 'no rows', 0, or an empty list; if the gap survives, it goes " +
        `in \`not_found\` quoting this error. ${failure.remedy}\n\n` +
        `--- paperclip ${argv[0]} said ---\n` +
        `${clip(result.stderr || result.stdout || "(no output at all)")}`
    );
  }
  return result;
}

/**
 * `report()` renders an empty run as "(no output)", which on this tool reads
 * exactly like a query that returned no rows. Say what it actually is instead.
 */
function reportPaperclip(label: string, result: RunResult): string {
  if (result.code === 0 && !(result.stdout.trim() || result.stderr.trim())) {
    return (
      `${label} exited 0 and printed nothing at all — no table, no header, no ` +
      "error. THIS IS AMBIGUOUS AND MUST NOT BE RECORDED AS ZERO ROWS. Re-run " +
      "a query you know returns rows (a single-accession lookup against " +
      "`uniprot_v.proteins` answers in ~7 ms) before writing any count, and " +
      "if that one is also silent the corpus is not answering and the gap " +
      "belongs in `not_found`."
    );
  }
  return report(label, result);
}

const paperclipSql: CustomToolSpec = {
  description:
    "Run one read-only SELECT against a Paperclip database and return the rendered table. " +
    "Use it for every structured precedent lookup — `-s proteins` reaches chembl_v/pdb_v/uniprot_v (drugs by accession, bioactivities, structures, Pfam), `-s trials` reaches the AACT-style ctgov schema, and omitting the source hits the paper corpus (documents, content_blocks, figures). " +
    "Three caveats are measured, not guessed: results are hard-capped at 200 rows, a server-side statement timeout kills long queries, and wide cells are truncated with a literal `...` at roughly 880 characters — which silently destroys json_agg output, so aggregate into separate columns instead of one JSON blob. " +
    "THERE IS NO CURSED TABLE AND BREADTH IS NOT THE DIAGNOSED CAUSE — both characterisations are withdrawn. A dry run lost 11 of 30 calls to what looked like timeouts on `SELECT 1 AS ok`, on a `GROUP BY assay_type` and on a 17-id `IN` list over `pdb_v.entry_ligands`. Re-measured with the CLI's sticky client state neutralised, those same three shapes answer in 7 ms, 31 ms and immediately. The cause was `cli_cwd` in the Paperclip config, which the CLI cds into before EVERY command and which makes every query a ~15 ms no-op printing `vsh: cd: /papers/: Permission denied` at exit 0. The handler now runs against its own config directory, so that class is closed; if you still see it, the operator has PAPERCLIP_CONFIG_DIR pointed somewhere poisoned. " +
    "Rewriting an `IN` list as a subquery is one measurement (a Pfam cross-reference join: 85.1 s inline, 2.2 s as a subquery), not a rule. Try both forms before concluding the data is absent, and prefer looping one key per query to widening a predicate. " +
    "FOUR FAILURE SIGNATURES, AND EVERY ONE OF THEM CAN ARRIVE ON STDOUT WITH EXIT 0: `[error] Request timed out`, `[error] Something went wrong. Please try again.`, `vsh: cd: …: Permission denied`, and `ERR: sql: unknown column` / `relation … does not exist`. The handler throws on all four with the raw text attached. A thrown error here is a failed lookup and never an empty result — it must never become a zero, an empty list, or 'no precedent found'; if the gap survives, it goes in `not_found` quoting the error.",
  async handler(input) {
    const query = requiredStr(input, "query");
    const source = str(input, "source");
    const argv = source ? ["sql", "-s", source, query] : ["sql", query];
    return reportPaperclip("paperclip sql", await paperclip(argv));
  },
  input_schema: {
    properties: {
      query: {
        description:
          "A single SELECT statement. Only SELECT is accepted; no DDL or DML.",
        type: "string",
      },
      source: {
        description:
          "Database to query, e.g. `proteins` (chembl_v/pdb_v/uniprot_v) or `trials` (ctgov). Omit for the paper corpus.",
        type: "string",
      },
    },
    required: ["query"],
    type: "object",
  },
  name: "paperclip_sql",
};

const paperclipGrep: CustomToolSpec = {
  description:
    "Full-text regex search over the Paperclip document corpus — `/papers/` for the literature and `/trials/` for registry records — returning matching lines with their document IDs. " +
    "Use it to find what a database column cannot express: a compound code, a termination reason, a stated adverse event. " +
    "The flags are short-form only and this is measured: `-C NUM` and `-m NUM` work, while the long forms `--context` and `--limit` are not parsed as flags and produce the misleading error `Cannot read path: /papers` — which reads like a missing corpus rather than a bad flag. " +
    "Two more measured traps: a hyphenated compound code collides with document UUIDs (3 of 20 hits for `DC-806` were substrings inside UUIDs, not mentions — always read the surrounding line before counting), and `/trials/` result paths are frequently not readable afterwards while the tool still labels its output as 'matching papers', so treat a trials hit as a pointer to check in `paperclip_sql -s trials`, not as a retrievable document.",
  async handler(input) {
    const pattern = requiredStr(input, "pattern");
    const path = str(input, "path") ?? "/papers/";
    const argv = ["grep"];
    if (bool(input, "ignore_case")) {
      argv.push("-i");
    }
    if (bool(input, "line_numbers")) {
      argv.push("-n");
    }
    if (bool(input, "count")) {
      argv.push("-c");
    }
    if (bool(input, "list_files")) {
      argv.push("-l");
    }
    if (bool(input, "fixed_string")) {
      argv.push("-F");
    }
    if (bool(input, "whole_word")) {
      argv.push("-w");
    }
    if (bool(input, "bool_mode")) {
      argv.push("--bool");
    }
    const context = num(input, "context");
    if (context !== undefined) {
      argv.push("-C", String(context));
    }
    const maxMatches = num(input, "max_matches");
    if (maxMatches !== undefined) {
      argv.push("-m", String(maxMatches));
    }
    const section = str(input, "section");
    if (section) {
      argv.push("--section", section);
    }
    const blockType = str(input, "block_type");
    if (blockType) {
      argv.push("--block-type", blockType);
    }
    argv.push(pattern, path);
    return reportPaperclip("paperclip grep", await paperclip(argv));
  },
  input_schema: {
    properties: {
      block_type: {
        description: "Restrict matches to a block type.",
        type: "string",
      },
      bool_mode: {
        description: "Whole-document Boolean regex (`--bool`).",
        type: "boolean",
      },
      context: {
        description:
          "Lines of context around each match. Sent as `-C NUM`; the long form does not work.",
        type: "number",
      },
      count: { description: "Count matches (`-c`).", type: "boolean" },
      fixed_string: {
        description: "Literal match, no regex (`-F`).",
        type: "boolean",
      },
      ignore_case: { description: "Case-insensitive (`-i`).", type: "boolean" },
      line_numbers: {
        description: "Show line numbers (`-n`).",
        type: "boolean",
      },
      list_files: {
        description: "List only matching document paths (`-l`).",
        type: "boolean",
      },
      max_matches: {
        description:
          "Stop after N matches per corpus. Sent as `-m NUM`; the long form does not work.",
        type: "number",
      },
      path: {
        description:
          "Corpus path: `/papers/`, `/trials/`, or a specific document such as `/papers/PMC8080595/content.lines`. Defaults to `/papers/`.",
        type: "string",
      },
      pattern: {
        description: "Regex (or literal, with fixed_string) to search for.",
        type: "string",
      },
      section: {
        description: "Restrict matches to a named section, e.g. `Methods`.",
        type: "string",
      },
      whole_word: { description: "Whole words only (`-w`).", type: "boolean" },
    },
    required: ["pattern"],
    type: "object",
  },
  name: "paperclip_grep",
};

const paperclipRead: CustomToolSpec = {
  description:
    "Read one file out of the Paperclip virtual filesystem — a document body (`/papers/<id>/content.lines`), its metadata (`/papers/<id>/meta.json`, giving doi, journal and pub_date), or a trial record (`/trials/us/<NCT>/meta.json`). " +
    "Use it after `paperclip_grep` returns a document ID, to read the passage in context and pull the citation you will put in a `source` field. " +
    "Pass `numbered: true` to get line numbers so a claim can be pinned to a line, which is what the dossier's provenance rule asks for. " +
    'Caveat, measured: `/trials/` paths that appear in grep output are frequently not readable here and return `Cannot read path`, so treat a failed trials read as expected and go to `paperclip_sql` with `source: "trials"` instead of retrying.',
  async handler(input) {
    const path = requiredStr(input, "path");
    const argv = ["cat"];
    if (bool(input, "numbered")) {
      argv.push("-n");
    }
    argv.push(path);
    return reportPaperclip("paperclip cat", await paperclip(argv));
  },
  input_schema: {
    properties: {
      numbered: {
        description: "Number the output lines (`-n`).",
        type: "boolean",
      },
      path: {
        description:
          "Path in the Paperclip VFS, e.g. `/papers/PMC8080595/meta.json`.",
        type: "string",
      },
    },
    required: ["path"],
    type: "object",
  },
  name: "paperclip_read",
};

const paperclipSearch: CustomToolSpec = {
  description:
    "Semantic + BM25 search over a named Paperclip source, returning ranked documents with their IDs. " +
    "Use it when you do not know the exact string to grep for — a mechanism, an indication, a programme description — and then follow up with `paperclip_grep` or `paperclip_read` on the IDs it returns. " +
    "`source` is mandatory on every search (unlike grep) and accepts comma-separated values such as `pmc`, `biorxiv`, `medrxiv`, `fda`, `trials/us` or `proteins`. " +
    "Caveat: `patents` is advertised but not provisioned and returns `Patents sources are not available.` for search, sql and ls alike, so patent counts must be reported as `null` with the reason in `not_found` rather than guessed.",
  async handler(input) {
    const query = requiredStr(input, "query");
    const source = requiredStr(input, "source");
    const argv = ["search", "-s", source];
    const limit = num(input, "limit");
    if (limit !== undefined) {
      argv.push("-n", String(limit));
    }
    argv.push(query);
    return reportPaperclip("paperclip search", await paperclip(argv));
  },
  input_schema: {
    properties: {
      limit: { description: "Maximum results (`-n`).", type: "number" },
      query: {
        description: "Natural-language or keyword query.",
        type: "string",
      },
      source: {
        description:
          "Required. One or more comma-separated sources, e.g. `pmc,biorxiv` or `trials/us`.",
        type: "string",
      },
    },
    required: ["query", "source"],
    type: "object",
  },
  name: "paperclip_search",
};

// ---------------------------------------------------------------------------
// pocket_scan — the Modal app
// ---------------------------------------------------------------------------

/**
 * The only Modal workspace this agent may bill or read. The founder was
 * explicit that the other workspaces in ~/.modal.toml — `molspace-production`
 * and `foldariumtest` — cannot be used, so this is a hard default with no
 * fallback, not a preference.
 */
const EXPECTED_MODAL_PROFILE = "rafwiewiora";

/**
 * `modal` now lives in the `druggability` micromamba env alongside fpocket,
 * gemmi and numpy, and `.env` pins MODAL_BIN to it (CREDENTIALS.md section 3).
 * It previously lived only in a throwaway venv under /private/tmp that would
 * not have survived a reboot, which is why there is still deliberately no
 * candidate path here: MODAL_BIN or PATH, and a loud throw otherwise. A guessed
 * fallback is what turns a five-minute fix into a mid-run mystery.
 */
function modalBin(): string {
  return resolveBin({
    candidates: [],
    envVar: "MODAL_BIN",
    fallbackOnPath: "modal",
    why: "It runs the fpocket/mdpocket image, which is the entire computed-tractability axis; there is no local fallback for it.",
  });
}

/**
 * Modal does not authenticate from `.env`; its token_id/token_secret live in
 * `~/.modal.toml`, one block per profile. So the credential check here is
 * "does the named profile exist in that file", not a `requireEnv`.
 *
 * Resolving the profile with `??` would have been a silent-fallback bug:
 * `MODAL_PROFILE=""` is not nullish, so an empty value would have passed
 * through and let Modal pick its own active profile — which is a different
 * workspace whenever someone has run `modal profile activate` elsewhere.
 * Treat blank as unset, and refuse an unknown profile by name.
 */
function modalProfile(): string {
  const raw = process.env.MODAL_PROFILE?.trim();
  const profile = raw && raw.length > 0 ? raw : EXPECTED_MODAL_PROFILE;
  const configPath = join(homeDir(), ".modal.toml");
  if (!existsSync(configPath)) {
    throw new Error(
      `no Modal config at ${configPath}, so profile "${profile}" cannot be ` +
        "authenticated. Run `modal token new --profile " +
        `${EXPECTED_MODAL_PROFILE}` +
        "`, or point MODAL_BIN at a Modal install whose config has it. " +
        "NOTE: that command only works once you are a MEMBER of the " +
        `"${EXPECTED_MODAL_PROFILE}" workspace — only its owner can add you, ` +
        "and a token for your own workspace is rejected by name. Setup path: " +
        "managed/druggability-dossier/CREDENTIALS.md section 8.2."
    );
  }
  // Existence is the wrong test: the forbidden workspaces are *also* in this
  // file, so "is it a real profile" would wave `molspace-production` straight
  // through. The test is identity.
  if (
    profile !== EXPECTED_MODAL_PROFILE &&
    process.env.MODAL_PROFILE_OVERRIDE !== profile
  ) {
    throw new Error(
      `MODAL_PROFILE is "${profile}", but this agent runs only in the ` +
        `"${EXPECTED_MODAL_PROFILE}" Modal workspace. The other profiles in ` +
        `${configPath} belong to different workspaces and must not be ` +
        "billed or read by this pipeline. Unset MODAL_PROFILE to use the " +
        "correct one. If you genuinely mean to switch, set " +
        "MODAL_PROFILE_OVERRIDE to the same value to acknowledge it. " +
        "If you have never used this workspace: you cannot substitute your " +
        "own. You must be ADDED to it by its owner, then run `modal token " +
        `new --profile ${EXPECTED_MODAL_PROFILE}` +
        "` and confirm with `modal app list` that you actually landed in the " +
        "shared workspace. Full path: " +
        "managed/druggability-dossier/CREDENTIALS.md section 8.2."
    );
  }
  const config = readFileSync(configPath, "utf8");
  if (!new RegExp(`^\\[${profile}\\]`, "m").test(config)) {
    throw new Error(
      `MODAL_PROFILE is "${profile}", which is not a profile in ${configPath}. ` +
        "Run `modal token new --profile " +
        `${profile}` +
        "` to create it, choosing the shared workspace in the browser flow, " +
        "then check `modal app list` is not empty — an empty list means the " +
        "profile name is right and points at the wrong workspace, which this " +
        "check cannot see. managed/druggability-dossier/CREDENTIALS.md " +
        "section 8.2."
    );
  }
  return profile;
}

/** Stage switches, in the order `modal_app.py`'s entrypoint declares them. */
const POCKET_SCAN_STAGES = ["run_disorder", "run_cryptic", "run_mdpocket"];

// ---------------------------------------------------------------------------
// pocket_scan payload reduction
//
// THE INVARIANT: this handler never returns a truncated JSON document. Not
// once, not with a warning. A payload cut mid-string parses as nothing, and a
// model handed nothing after a paid Modal run has no way to tell a wide result
// from a broken one.
//
// The previous version broke that invariant in three ways at once, measured on
// a three-structure IL-6 (P05231) run that came back at 180,209 characters:
//
//   1. It dropped `pockets` and nothing else. `pockets` is not where the bulk
//      is — with all six `pockets` arrays gone (98 objects) the payload was
//      STILL over the cap, because two `structures.<ID>` blocks are ~118 kB on
//      their own and `pocket_vs_interface.per_structure` is another ~60 kB.
//   2. It then called `clip`, which hard-truncated the reduced document and
//      handed the model invalid JSON anyway.
//   3. `_handler_note` was appended as the LAST key, so `JSON.stringify` put
//      it at the very end and `clip` deleted the explanation first.
//
// The consequence was concrete and expensive: on IL-6 the tool could not return
// `mdpocket.sites` (rule 4b) and `pocket_vs_interface` (rule 2b) in one
// parseable payload, so the axis cost two paid Modal runs.
//
// So: a LADDER of named reductions, cheapest first, re-measured after every
// rung, with `_handler_note` written FIRST, and a hard throw naming the size if
// the ladder runs out. No rung recomputes, rounds or summarises a number —
// every rung deletes whole keys, and every deletion is named in the note.
// ---------------------------------------------------------------------------

/** Prose keys shorter than this are flags, not commentary; they stay. */
const PROSE_MIN_CHARS = 160;
/** `_why`, `_note`, `tier_note`, `_aggregation_rule`, `pockets_omitted_note`… */
const PROSE_KEY = /^_|_(?:note|why|warning|caveat|rule)$/;

type JsonObject = Record<string, unknown>;

function isProse(key: string, value: unknown): boolean {
  return (
    typeof value === "string" &&
    value.length >= PROSE_MIN_CHARS &&
    PROSE_KEY.test(key)
  );
}

/** Depth-first walk over every plain object in the tree, including in arrays. */
function walkObjects(node: unknown, visit: (obj: JsonObject) => void): void {
  if (Array.isArray(node)) {
    for (const item of node) {
      walkObjects(item, visit);
    }
    return;
  }
  if (!node || typeof node !== "object") {
    return;
  }
  const obj = node as JsonObject;
  visit(obj);
  for (const value of Object.values(obj)) {
    walkObjects(value, visit);
  }
}

/** Replace `obj[key]` with a marker naming what left. Returns 1 if it fired. */
function dropKey(obj: JsonObject, key: string, what: string): number {
  const value = obj[key];
  if (value === undefined || typeof value === "string") {
    return 0;
  }
  const n = Array.isArray(value)
    ? value.length
    : Object.keys(value ?? {}).length;
  obj[key] =
    `[${n} ${what} dropped by the local tool handler to fit the output cap — NOT absent, NOT zero]`;
  return 1;
}

/** Rung 1: method commentary. Constant overhead, carries no measurement. */
function dropProse(root: JsonObject, names: Set<string>): number {
  let n = 0;
  walkObjects(root, (obj) => {
    for (const [key, value] of Object.entries(obj)) {
      if (isProse(key, value)) {
        delete obj[key];
        names.add(key);
        n += 1;
      }
    }
  });
  return n;
}

/**
 * Rung 2: `pocket_vs_interface.per_structure`.
 *
 * ~60 kB on the IL-6 run and it is the one large block that is a DUPLICATE:
 * modal_app.py copies every field the dossier reads into
 * `structures.<ID>.pocket_vs_interface.<D>` and summarises the same entries in
 * `pocket_vs_interface.per_structure_consensus`. Rule 2b survives this rung.
 */
function dropInterfacePerStructure(root: JsonObject): number {
  const block = root.pocket_vs_interface as JsonObject | undefined;
  if (!block) {
    return 0;
  }
  return dropKey(
    block,
    "per_structure",
    "raw per-structure classification dicts (the same fields are in structures.<ID>.pocket_vs_interface.<D> and in per_structure_consensus)"
  );
}

/** Rung 3: the top-30 pocket lists. `site_pocket` and the ranks stay. */
function dropPocketLists(root: JsonObject): number {
  let n = 0;
  for (const entry of Object.values(
    (root.structures as JsonObject | undefined) ?? {}
  )) {
    const byClustering = (entry as JsonObject)?.by_clustering as
      | Record<string, JsonObject>
      | undefined;
    for (const block of Object.values(byClustering ?? {})) {
      n += dropKey(block, "pockets", "pocket objects");
    }
  }
  return n;
}

/** Rung 4: per-rank interface classification. This one DOES cost rule 2b. */
function dropPerRankClassification(root: JsonObject): number {
  let n = 0;
  for (const entry of Object.values(
    (root.structures as JsonObject | undefined) ?? {}
  )) {
    const pvi = (entry as JsonObject)?.pocket_vs_interface as
      | Record<string, JsonObject>
      | undefined;
    for (const block of Object.values(pvi ?? {})) {
      n += dropKey(block, "by_fpocket_rank", "per-rank classifications");
    }
  }
  return n;
}

/** Rung 5: residue-name lists that duplicate or annotate `residues`. */
const BULK_LISTS = ["lining_residue_names", "missing_residues"];

function dropBulkLists(root: JsonObject): number {
  let n = 0;
  walkObjects(root, (obj) => {
    for (const key of BULK_LISTS) {
      if (Array.isArray(obj[key])) {
        n += dropKey(obj, key, `${key} entries`);
      }
    }
  });
  return n;
}

/**
 * Rung for `cryptic_analysis`: `displacement.per_residue`.
 *
 * One record per compared residue. Every figure rule 5 reads off this block is
 * computed by the script and sits beside the list — `max_ca_displacement`,
 * `max_ca_displacement_at`, `ca_rmsd`, `n_residues`, `n_residues_compared` —
 * so the cryptic-mechanism classification survives this rung intact.
 */
function dropPerResidueDisplacement(root: JsonObject): number {
  const block = root.displacement as JsonObject | undefined;
  if (!block) {
    return 0;
  }
  return dropKey(
    block,
    "per_residue",
    "per-residue C-alpha displacement records (max_ca_displacement, max_ca_displacement_at, ca_rmsd and n_residues_compared are computed from these and are still here)"
  );
}

/** The two clash reports `cryptic_analysis` emits, at 2.0 A and at 2.5 A. */
const CONTACT_BLOCKS = ["contacts", "contacts_wide"];

/**
 * Rung for `cryptic_analysis`: the per-residue contact lists.
 *
 * `contact_pairs`, `n_ligand_atoms`, `n_protein_atoms`, `min_distance` and
 * `by_category` stay on each block, and the backbone / side-chain /
 * displaced-chain split rule 5 attributes clashes with is a separate top-level
 * `counts` key that this rung does not touch.
 */
function dropContactResidueLists(root: JsonObject): number {
  let n = 0;
  for (const key of CONTACT_BLOCKS) {
    const block = root[key] as JsonObject | undefined;
    if (block) {
      n += dropKey(
        block,
        "residues",
        `contacting residues in \`${key}\` (contact_pairs, by_category and min_distance are still here, and the backbone/sidechain/displaced_chain split is the top-level \`counts\` key)`
      );
    }
  }
  return n;
}

/**
 * Rung for `neighbour_precedent`: per-neighbour `rejected_ligands`.
 *
 * `ligand_filter`'s diagnostics saying why a candidate component was not
 * counted as drug-like. The four fields the axis is actually read from —
 * `has_druglike_holo`, `ligands`, `holo_determined` and `undetermined_ligands`,
 * which are what separate a real apo call from a failed lookup — are untouched.
 */
function dropRejectedLigands(root: JsonObject): number {
  const { neighbours } = root;
  if (!Array.isArray(neighbours)) {
    return 0;
  }
  let n = 0;
  for (const entry of neighbours) {
    if (entry && typeof entry === "object") {
      n += dropKey(
        entry as JsonObject,
        "rejected_ligands",
        "rejected-ligand diagnostics (holo_determined, undetermined_ligands, has_druglike_holo and ligands are untouched, so `n_holo = 0` is still readable against `n_undetermined`)"
      );
    }
  }
  return n;
}

type Rung = {
  apply: (root: JsonObject, names: Set<string>) => number;
  what: string;
};

/**
 * One ladder for every handler, ordered cheapest-first.
 *
 * Each rung looks for keys only some payloads have and returns 0 on the rest,
 * so `pocket_scan`, `cryptic_analysis` and `neighbour_precedent` share this
 * list rather than keeping one apiece. That is deliberate: a duplicated
 * reducer drifts, and the rung that deletes the note explaining the deletion
 * is the one nobody notices drifting.
 */
const REDUCTION_LADDER: Rung[] = [
  {
    apply: dropProse,
    what: "prose-only keys (string-valued `_why`/`_note`/`_warning` commentary, no numbers in any of them)",
  },
  {
    apply: dropInterfacePerStructure,
    what: "`pocket_vs_interface.per_structure` (duplicated into structures.<ID>.pocket_vs_interface and per_structure_consensus)",
  },
  {
    apply: dropPerResidueDisplacement,
    what: "`displacement.per_residue` (its max, argmax, RMSD and residue counts are computed by the script and are still here)",
  },
  { apply: dropPocketLists, what: "`by_clustering.<D>.pockets`" },
  {
    apply: dropContactResidueLists,
    what: "`contacts.residues` and `contacts_wide.residues` (contact_pairs, by_category, min_distance and the top-level backbone/sidechain/displaced_chain `counts` are still here)",
  },
  {
    apply: dropRejectedLigands,
    what: "per-neighbour `rejected_ligands` (holo_determined and undetermined_ligands are untouched)",
  },
  {
    apply: dropPerRankClassification,
    what: "`structures.<ID>.pocket_vs_interface.<D>.by_fpocket_rank` — THIS ONE COSTS RULE 2b's per-pocket classification; the selected site pocket's classification is still there",
  },
  {
    apply: dropBulkLists,
    what: "`lining_residue_names` and `missing_residues`",
  },
];

/**
 * What a reduction or a refusal needs to know about the handler it is for.
 *
 * `provenance` is named purely as a pointer for a human: every handler here
 * runs off-sandbox, so the complete payload is readable by the operator and
 * NOT by the agent. It is never a retrieval route for the model.
 */
type FitTarget = {
  /** Named in the note and in every refusal. */
  label: string;
  /** Where the complete payload is. Reads after "The complete payload is on". */
  provenance: string;
  /** How to make the next call smaller. Appended to every refusal. */
  advice: string;
};

function handlerNote(args: {
  fullChars: number;
  applied: string[];
  proseKeys: Set<string>;
  target: FitTarget;
}): string {
  const prose =
    args.proseKeys.size > 0
      ? ` The prose keys removed were: ${[...args.proseKeys].sort().join(", ")}.`
      : "";
  return (
    "READ THIS FIRST — THIS PAYLOAD IS COMPLETE JSON BUT IT IS NOT THE WHOLE " +
    `RESULT. ${args.target.label} returned ${args.fullChars} characters, ` +
    `above this handler's ${MAX_OUTPUT_CHARS}-character cap, so whole keys ` +
    "were deleted, cheapest first, until it fit. NOTHING WAS RECOMPUTED, " +
    "ROUNDED OR SUMMARISED — every number still here is the tool's own, and " +
    "no string was truncated. Removed, in order: " +
    args.applied.map((item, i) => `(${i + 1}) ${item}`).join("; ") +
    "." +
    prose +
    " A key whose value is now a `[… dropped by the local tool handler …]` " +
    "string was PRESENT AND NON-EMPTY; it is not absent and not zero. The " +
    `complete payload is on ${args.target.provenance}. ${args.target.advice}`
  );
}

/**
 * The reduction did not have to delete anything — the payload was only over
 * the cap because it pretty-prints. Say so, so nobody reads a note as a loss.
 */
function whitespaceNote(fullChars: number, target: FitTarget): string {
  return (
    "READ THIS FIRST — THIS PAYLOAD IS COMPLETE AND NOTHING WAS DELETED. " +
    `${target.label} returned ${fullChars} characters, above this handler's ` +
    `${MAX_OUTPUT_CHARS}-character cap, only because it pretty-prints its ` +
    "JSON. It has been re-serialised without the indentation and now fits. No " +
    "key was dropped, no string was truncated and no number was recomputed — " +
    "the only bytes removed are insignificant whitespace."
  );
}

function withNote(parsed: JsonObject, note: string): string {
  // `_handler_note` FIRST, never last: JSON.stringify emits string keys in
  // insertion order, and the explanation of a reduction must not be the first
  // thing any downstream size limit deletes.
  return JSON.stringify({ _handler_note: note, ...parsed });
}

function sizeCensus(parsed: JsonObject): string {
  return Object.entries(parsed)
    .map(([key, value]) => `${key}=${JSON.stringify(value)?.length ?? 0}`)
    .sort()
    .join(", ");
}

/** A plain JSON object. An array or a scalar cannot carry `_handler_note`. */
function asJsonObject(value: unknown): JsonObject | undefined {
  return value && typeof value === "object" && !Array.isArray(value)
    ? (value as JsonObject)
    : undefined;
}

/**
 * Fit an already-parsed payload under the cap by deleting whole keys, or refuse.
 *
 * `lead` is a note that must reach the model whether or not anything was
 * reduced (the chain-rename note); it is written FIRST for exactly the reason
 * the reduction note is.
 */
function fitParsed(args: {
  parsed: JsonObject;
  fullChars: number;
  target: FitTarget;
  lead?: string;
}): string {
  const { parsed, fullChars, target } = args;
  const applied: string[] = [];
  const proseKeys = new Set<string>();
  const emit = (note: string | undefined, indent: number) => {
    const joined = [args.lead, note].filter(Boolean).join(" ");
    return JSON.stringify(
      joined ? { _handler_note: joined, ...parsed } : parsed,
      null,
      indent
    );
  };

  // Rung 0 — whitespace. Deletes nothing, so it runs before any rung that
  // deletes a key. `cryptic_analysis` prints with indent=2 and roughly halves
  // on this rung alone.
  const pretty = emit(undefined, 2);
  if (pretty.length <= MAX_OUTPUT_CHARS) {
    return pretty;
  }
  const compact = emit(whitespaceNote(fullChars, target), 0);
  if (compact.length <= MAX_OUTPUT_CHARS) {
    return compact;
  }

  for (const rung of REDUCTION_LADDER) {
    const hits = rung.apply(parsed, proseKeys);
    if (hits > 0) {
      applied.push(rung.what);
    }
    const candidate = emit(
      handlerNote({ applied, fullChars, proseKeys, target }),
      0
    );
    // Re-measured after EVERY rung. The old version measured once, before
    // reducing, and returned whatever came out.
    if (candidate.length <= MAX_OUTPUT_CHARS) {
      return candidate;
    }
  }
  const remaining = withNote(parsed, "").length;
  throw new Error(
    `the ${target.label} payload is ${fullChars} characters and every ` +
      `reduction this handler has still leaves ${remaining}, above the ` +
      `${MAX_OUTPUT_CHARS}-character cap. It is NOT being truncated: a JSON ` +
      "document cut mid-string parses as nothing and would read as a failed " +
      "run rather than a wide one, so this is a refusal, not a result. " +
      `Everything already removed: ${applied.join("; ") || "nothing — no rung on this handler's ladder matched a key in this payload"}. ` +
      `Top-level key sizes in what is left: ${sizeCensus(parsed)}. The ` +
      `complete payload is on ${target.provenance}. ${target.advice}`
  );
}

/**
 * THE INVARIANT, for any handler: reduce to a document that still parses, or
 * throw with the size. Never return a payload cut mid-string.
 *
 * A payload that is not a JSON object has no key to delete, so it is refused
 * outright with its byte count. That is the correct answer and not a
 * degradation: a truncated line-oriented listing reads exactly like a complete
 * shorter one, which on this agent is the single error the dossier cannot
 * survive — a failed retrieval wearing the costume of a finding.
 */
function fitText(text: string, target: FitTarget): string {
  let parsed: JsonObject | undefined;
  try {
    parsed = asJsonObject(JSON.parse(text));
  } catch {
    parsed = undefined;
  }
  if (!parsed) {
    throw new Error(
      `${target.label} returned ${text.length} characters, above this ` +
        `handler's ${MAX_OUTPUT_CHARS}-character cap, and that text is not a ` +
        "JSON object, so it cannot be reduced key by key. THIS IS A REFUSAL, " +
        "NOT AN EMPTY RESULT AND NOT A SHORT ONE. Returning a copy cut " +
        "mid-line would hand you a listing indistinguishable from a complete " +
        "shorter listing, so the handler refuses instead. The complete " +
        `output is on ${target.provenance}. ${target.advice}\n\n` +
        `--- first 2000 characters ---\n${text.slice(0, 2000)}`
    );
  }
  return fitParsed({ fullChars: text.length, parsed, target });
}

/**
 * The never-truncate output path for a completed local run.
 *
 * A NON-ZERO exit still goes through `report`, which clips, and that asymmetry
 * is deliberate: a failure dump is a stdout/stderr transcript for a human to
 * read, not a document anything parses, and clipping it behind its own marker
 * keeps the error text reachable where refusing would delete it. `pocket_scan`
 * has always done the same — see `readScanPayload`.
 */
function fitReport(result: RunResult, target: FitTarget): string {
  if (result.code !== 0) {
    return report(target.label, result);
  }
  const text = result.stdout || result.stderr || "(no output)";
  return text.length > MAX_OUTPUT_CHARS ? fitText(text, target) : text;
}

const CRYPTIC_FIT: FitTarget = {
  advice:
    "Narrow the comparison rather than re-running it whole: pass `holo_chains` and `apo_chains` to compare one protomer pair instead of the entire assembly, or set `no_free_volume`.",
  label: "cryptic_analysis.py",
  provenance:
    "the machine that ran this handler, which did not save it to a file — the sandbox cannot see it, so ask the operator to re-run the script there",
};

const INTERFACE_FIT: FitTarget = {
  advice:
    "Ask for one accession at a time. `selftest_dir` output is a diagnostic harness dump and is never evidence about a target, so never widen the cap for it.",
  label: "interface_analysis.py",
  provenance:
    "the machine that ran this handler, which did not save it to a file — the sandbox cannot see it, so ask the operator to re-run the script there",
};

const DISORDER_FIT: FitTarget = {
  advice:
    "Send fewer accessions per call — the output is one line per accession, so splitting the list loses nothing at all.",
  label: "disorder.py",
  provenance:
    "the machine that ran this handler, which did not save it to a file — the sandbox cannot see it, so ask the operator to re-run the script there",
};

const NEIGHBOUR_FIT: FitTarget = {
  advice:
    "Lower `max_neighbours`, or raise `min_alignment_length`, and say in `not_found` that the neighbour list was capped by this handler rather than by the fold — a capped list is not a count of neighbours.",
  label: "neighbour_precedent.py",
  provenance:
    "the machine that ran this handler, which did not save it to a file — the sandbox cannot see it, so ask the operator to re-run the script there",
};

/** Fit the `pocket_scan` payload under the cap by deleting whole keys, or refuse. */
function reduceScanPayload(text: string, outFile: string): string {
  return fitText(text, {
    advice:
      "Re-run with a smaller `pdb_ids` list, or split the stages — one call with `run_mdpocket` for rule 4b's `mdpocket.sites` and one with `partner_structures` for rule 2b's `pocket_vs_interface` — and say in `tractability.caveat` that the axis was assembled from two runs.",
    label: "pocket_scan",
    provenance: `${outFile} on the machine running this handler (the sandbox cannot read that path), so ask the operator for it rather than re-running the scan, which costs Modal credits`,
  });
}

/**
 * Retrieve the payload `--out` wrote.
 *
 * The CLI no longer puts the JSON on stdout, and this is not a preference:
 * `modal run` frames stdout with its own progress banner and a trailing
 * "Stopping app…", so the old handler's `report()` of stdout returned banner
 * text with a JSON body embedded in it — or, once the entrypoint moved the
 * payload to stderr, no payload at all on a *successful* run. A clean exit with
 * nothing in the file is therefore a hard failure here rather than an empty
 * string that reads like a scan finding no pockets.
 */
function readScanPayload(outFile: string, result: RunResult): string {
  if (result.code !== 0) {
    return report("modal run modal_app.py", result);
  }
  if (existsSync(outFile)) {
    const text = readFileSync(outFile, "utf8");
    if (text.trim().length > 0) {
      return text.length > MAX_OUTPUT_CHARS
        ? reduceScanPayload(text, outFile)
        : text;
    }
  }
  throw new Error(
    "modal run exited 0 but wrote no JSON to its --out file. This is NOT an " +
      "empty scan result and must not be recorded as one — the run either " +
      "never reached the entrypoint's write, or Modal killed it after the " +
      "function returned.\n\n" +
      `--- stdout ---\n${clip(result.stdout)}\n\n--- stderr ---\n${clip(result.stderr)}`
  );
}

const pocketScan: CustomToolSpec = {
  description:
    "Run the whole computed-tractability half of the dossier in one Modal invocation: fpocket + PRANK at D = 1.6 and 2.4 (the sweep is done for you — there is no clustering_d argument; read `clustering_swept` in the method block), plus the disorder, cryptic, interface and mdpocket stages. " +
    "Use it for every pocket measurement — fpocket and mdpocket are conda-forge binaries that exist only inside this Modal image, so there is no other way to get a volume, a druggability range, or a site fixed by construction. " +
    "Pass the full ensemble at once (`pdb_ids`) rather than one structure per call: one invocation pays one cold start, and same-site tracking only works when a holo structure sits in the same run as the apo ones it anchors. " +
    "READ VOLUME, NOT DRUGGABILITY. `pocket_volume_a3` at D=1.6 is the primary number (target-level AUC 1.000 over 15 targets); fpocket's druggability score does NOT separate druggable from hard — AUC 0.720 with a 95% CI of 0.44-0.94 at D=1.6 and 0.520 (chance) at D=2.4, and on 37 holo structures with a drug-like ligand physically bound 41% score below 0.1 (EGFR with osimertinib in it scores 0.013). Report it as a range, never let it carry a verdict, and do not substitute persistence (AUC 0.500) for it. mdpocket's druggability field is null BY DESIGN, not by failure: fpocket's score is min-max normalised across the other pockets of the same structure, and a fixed grid has a population of one, so the quantity is undefined there. " +
    "CHAIN SELECTION IS THE ARGUMENT THAT CHANGES THE ANSWER, so assert it. `chains` takes `{\"1TNF\": [\"A\",\"B\"]}` and `site_residues` takes residue numbers; together they express rule 2b directly. KRAS 4OBE gives druggability 0.442 at rank 1 on chain A and 0.257 at rank 6 on chains A+B — same structure, same clustering, different verdict. They also unlock the subunit-removed control: TNF-alpha's SPD304 site measures 0.00 A^3 intact and ~280-550 A^3 with one protomer deleted, which is the experiment that separates 'the cavity is too small' from 'a protomer is standing in it'. A chain flag alone is not always enough — 3V2Y's T4-lysozyme fusion sits INSIDE chain A at 1002-1161 beside the receptor at 16-330, which needs `site_residues`. Record what you passed in `tractability.method.chains_used`. " +
    "PASS `uniprot_accession` OR THE DISORDER NUMBER IS ABOUT A DIFFERENT MOLECULE. Without it the stage can fall back to the crystallised construct, which is the ordered part of the protein by selection: IRAK4 came back 0.0 over 284 residues against a true 0.1413 over 460, and a bare 0.0 in `disorder_fraction` reads as 'no disorder', not as 'not measured'. A construct-scoped number is reported in `construct_disorder_fraction` with `is_full_length_sequence: false`, never in `disorder_fraction` — if you see that key, you omitted the accession. Carry `disorder.method` beside any fraction you quote; metapredict and the MobiDB fallback differ by 23% on the same target. " +
    "HOLO/APO IS DECIDED FROM SMILES, AND WITHOUT SMILES IT FAILS SILENTLY. Ligands are classified on their SMILES graph, so a record source with no SMILES returns `unknown` for every component, nothing is `druglike`, and the whole ensemble comes back apo-free and holo-free while every `<stage>_status` still says ok. Sources that carry SMILES: RCSB REST `data.rcsb.org/rest/v1/core/chemcomp/<ID>` (this image's source), Paperclip `pdb_v.chemcomps`, and the CCD ligand file. The entry's OWN mmCIF `_chem_comp` block does NOT — it carries id, type, name, formula and formula_weight and nothing else, and it is the obvious one to reach for because the file is already on disk. This app refuses such a source with a run-killing `LigandSourceError`, so that error is a misconfiguration, never a result. " +
    '`unknown` IS NOT `apo`, AND A LOOKUP TIMEOUT IS NOT A CHEMISTRY MISS. Read `holo_call(ids)["determined"]` before calling anything apo; components whose lookup failed carry `lookup_failed` and land in `["undetermined"]`, a third tier. A genuine CCD 404 caches as absent and is not an error. ' +
    "Remaining caveats that decide how you read the result: `ligand_codes` is an override and not a requirement (naming one code across four structures left three falling back to the weaker signature path and moved 7JRA from 0.000/306.9 A^3 to 0.926/1542.9 A^3), `mdpocket.sites` returns up to two definitions of which only `site_from_ligand` is the ligand site — read `distance_to_donor_ligand_centroid_a` on every entry before quoting a number — `site_pocket_selected_by` says whether a spread describes one site at all, and every stage after fpocket is non-fatal and reports its own `<stage>_status`, so check those before treating a missing block as a null result. A run takes minutes and costs real credits; do not re-run it to retry a formatting question.",
  async handler(input) {
    const pdbIds = list(input, "pdb_ids");
    if (pdbIds.length === 0) {
      throw new Error(
        "pdb_ids is required and must name at least one PDB entry"
      );
    }
    const modal = modalBin();
    const profile = modalProfile();
    const outFile = join(
      mkdtempSync(join(tmpdir(), "pocket-scan-")),
      "scan.json"
    );
    const argv = [
      "run",
      join(SKILLS_DIR, "pocket-scan", "modal_app.py"),
      "--pdb-ids",
      pdbIds.join(","),
      // The payload no longer reaches stdout at all; see readScanPayload.
      "--out",
      outFile,
    ];
    pushFlag(argv, "--chains", chainSpec(input, "chains"));
    const siteResidues = numList(input, "site_residues");
    pushFlag(
      argv,
      "--site-residues",
      siteResidues.length > 0 ? siteResidues.join(",") : undefined
    );
    const ligandCodes = list(input, "ligand_codes");
    pushFlag(
      argv,
      "--ligand-codes",
      ligandCodes.length > 0 ? ligandCodes.join(",") : undefined
    );
    pushFlag(argv, "--uniprot-accession", str(input, "uniprot_accession"));
    const partners = list(input, "partner_structures");
    pushFlag(
      argv,
      "--partner-structures",
      partners.length > 0 ? partners.join(",") : undefined
    );
    pushFlag(argv, "--mdpocket-site-donor", str(input, "mdpocket_site_donor"));
    for (const stage of POCKET_SCAN_STAGES) {
      if (optBool(input, stage) === false) {
        argv.push(`--no-${stage.replaceAll("_", "-")}`);
      }
    }
    const result = await run(modal, argv, {
      env: { ...process.env, MODAL_PROFILE: profile },
      timeoutSeconds: POCKET_SCAN_TIMEOUT_S,
    });
    return readScanPayload(outFile, result);
  },
  input_schema: {
    properties: {
      chains: {
        additionalProperties: { items: { type: "string" }, type: "array" },
        description:
          'Chains to keep, per PDB entry: `{"1TNF": ["A","B"], "2ZJC": ["A","B"]}`. This is rule 2b\'s chain selection and it changes the answer (KRAS 4OBE: 0.442 at rank 1 on chain A, 0.257 at rank 6 on A+B), so assert it rather than defaulting to the whole assembly. Deleting a protomer here IS the subunit-removed control. Omit only when no mechanism hypothesis was supplied, and say so in `tractability.caveat`.',
        type: "object",
      },
      ligand_codes: {
        description:
          'Optional chemical component IDs used to anchor the site, e.g. `["MOV"]`. An override, not a requirement: a structure carrying its own drug-like ligand anchors itself when no supplied code matches.',
        items: { type: "string" },
        type: "array",
      },
      mdpocket_site_donor: {
        description:
          "A holo PDB ID used ONLY to define the mdpocket site, not added to the ensemble. This is how a pure-apo ensemble gets a ligand-anchored site. It also donates the site signature to the fpocket pass and the holo half of the cryptic comparison.",
        type: "string",
      },
      partner_structures: {
        description:
          'PDB IDs of complexes containing the binding partner, e.g. `["3ALQ"]`. Turns the orthosteric/allosteric question into a measurement.',
        items: { type: "string" },
        type: "array",
      },
      pdb_ids: {
        description:
          'The ensemble, e.g. `["6OIM","4OBE"]`. Send every structure in one call.',
        items: { type: "string" },
        type: "array",
      },
      run_cryptic: {
        description:
          "Defaults to true. Set false only to skip the apo/holo superposition stage; it is non-fatal either way and reports `cryptic_status`.",
        type: "boolean",
      },
      run_disorder: {
        description:
          "Defaults to true. Set false only to skip the disorder stage — never as a substitute for supplying `uniprot_accession`.",
        type: "boolean",
      },
      run_mdpocket: {
        description:
          "Defaults to true. Set false to skip the fixed-grid stage, which is the slowest one and the only source of `mdpocket.sites`.",
        type: "boolean",
      },
      site_residues: {
        description:
          "Residue numbers defining the site signature, e.g. `[57,58,59,60,61]`. Use when a chain flag cannot express the selection — a fusion chaperone inside the same chain (3V2Y: T4 lysozyme at 1002-1161, receptor at 16-330), or an allosteric domain picked by range. Matched chain-agnostically, so it is unreliable on a homo-oligomer, where the run reports `site_signature_unreliable_homooligomer`.",
        items: { type: "integer" },
        type: "array",
      },
      uniprot_accession: {
        description:
          "REQUIRED IN PRACTICE. Drives the disorder stage. Without it the module can measure the crystallised construct instead of the protein — IRAK4 returned 0.0 over 284 residues against a true 0.1413 over 460 — and that 0.0 is a wrong number, not a null.",
        type: "string",
      },
    },
    required: ["pdb_ids"],
    type: "object",
  },
  name: "pocket_scan",
};

// ---------------------------------------------------------------------------
// Local structure-analysis scripts (gemmi + numpy)
// ---------------------------------------------------------------------------

function pythonArgv(script: string, scriptArgs: string[]): string[] {
  return ["run", "-n", analysisEnvName(), "python", script, ...scriptArgs];
}

const PDB_ID = /^[0-9][A-Za-z0-9]{3}$/;
const STRUCTURE_CACHE = join(tmpdir(), "druggability-dossier-structures");

/**
 * These handlers run locally, so a structure the *sandbox* downloaded is not
 * reachable here. Accept a 4-character PDB ID and materialise the biological
 * assembly locally; anything else is passed through as a path untouched. This
 * is not a general fetch tool — the sandbox has open egress and can curl RCSB
 * itself; it is the local half of a local-only script's input.
 */
async function resolveStructure(value: string): Promise<string> {
  if (!PDB_ID.test(value)) {
    return value;
  }
  const id = value.toUpperCase();
  mkdirSync(STRUCTURE_CACHE, { recursive: true });
  const target = join(STRUCTURE_CACHE, `${id}-assembly1.cif`);
  if (existsSync(target)) {
    return target;
  }
  const url = `https://files.rcsb.org/download/${id}-assembly1.cif`;
  const response = await fetch(url);
  if (!response.ok) {
    throw new Error(
      `could not download ${url} (HTTP ${response.status}). Pass an explicit ` +
        "file path instead, or check the entry has a deposited assembly 1."
    );
  }
  writeFileSync(target, Buffer.from(await response.arrayBuffer()));
  return target;
}

// ---------------------------------------------------------------------------
// Assembly chain names vs the PDB format's one-character column
//
// `resolveStructure` downloads `<ID>-assembly1.cif`, and an assembly CIF names
// chains in ways the PDB format cannot hold: RCSB disambiguates symmetry copies
// with a suffix (1ALU: `A`, `A-2`) and renames large assemblies to multi-letter
// ids (7NXZ: `AAA`). `neighbour_precedent.py`'s `count_chains` converts to PDB
// before anything else runs, so both entries died on arrival:
//
//   1ALU: chain name too long for the PDB format: A-2
//   7NXZ: chain name too long for the PDB format: AAA
//
// WHY NOT THE ASYMMETRIC UNIT. The standing rule is the biological assembly,
// never the ASU, and it still holds here — the ASU is the wrong input for a
// Foldseek query specifically, because the multimer search runs only when the
// file has more than one chain, and an ASU that happens to be one protomer of a
// dimeric assembly would silently route a multimer query down the single-chain
// path and return neighbours of a monomer. What has to change is the chain
// NAMES, and only the names: Foldseek's own `target_id` is a filename blob whose
// chain token the parser already strips, so no downstream number is keyed on
// them. `gemmi.Structure.shorten_chain_names()` is the operation the format
// conversion needs and it renames nothing else — coordinates, residues,
// entities and the assembly's chain count are untouched.
//
// Scoped to `neighbour_precedent` on purpose. `cryptic_analysis` takes
// `--holo-chains`/`--apo-chains` as deposited chain ids, so renaming underneath
// it would silently re-point a caller's chain selection at a different protomer.
// ---------------------------------------------------------------------------

const PDB_CHAIN_RENAME = `
import json, sys
import gemmi
src, dst = sys.argv[1], sys.argv[2]
st = gemmi.read_structure(src)
st.setup_entities()
before = [ch.name for ch in st[0]]
if all(len(n) <= 1 for n in before):
    print(json.dumps({"renamed": None}))
else:
    st.shorten_chain_names()
    after = [ch.name for ch in st[0]]
    st.make_mmcif_document().write_file(dst)
    print(json.dumps({"renamed": dict(zip(before, after))}))
`;

const FILE_EXTENSION = /\.[^.]+$/;

type ChainRename = { path: string; renamed: Record<string, string> | null };

/**
 * Return a copy of `path` whose chain names fit the PDB format's one character,
 * or `path` itself when they already do.
 */
async function pdbSafeChainNames(path: string): Promise<ChainRename> {
  mkdirSync(STRUCTURE_CACHE, { recursive: true });
  const base = path.split("/").pop() ?? "structure.cif";
  const dst = join(
    STRUCTURE_CACHE,
    `${base.replace(FILE_EXTENSION, "")}-pdbchains.cif`
  );
  const result = await run(
    micromamba(),
    [
      "run",
      "-n",
      analysisEnvName(),
      "python",
      "-c",
      PDB_CHAIN_RENAME,
      path,
      dst,
    ],
    { timeoutSeconds: DEFAULT_TIMEOUT_S }
  );
  if (result.code !== 0) {
    throw new Error(
      `could not read chain names out of ${path} with gemmi, so this handler ` +
        "cannot tell whether the file will survive the CIF→PDB conversion " +
        "`neighbour_precedent.py` does before it searches. This is a broken " +
        `analysis environment ("${analysisEnvName()}"), NOT unavailability of ` +
        "the structural-neighbour axis and NOT a target with no neighbours — " +
        "do not null the axis for it.\n\n" +
        `--- stderr ---\n${clip(result.stderr)}`
    );
  }
  const parsed = JSON.parse(result.stdout.trim()) as {
    renamed: Record<string, string> | null;
  };
  return parsed.renamed
    ? { path: dst, renamed: parsed.renamed }
    : { path, renamed: null };
}

/** A CIF→PDB chain-name death, which rule 13 must NOT be applied to. */
const CHAIN_NAME_TOO_LONG = /chain name too long for the PDB format/i;

/**
 * `query_structure` must be the ENTRY, not the file the script was handed.
 *
 * `neighbour_precedent.py` reports `"query_structure": str(structure_path)`,
 * which was defensible when the caller passed a path and is not any more: this
 * handler fetches biological assembly 1 into a temp directory and then caches a
 * chain-renamed copy beside it, so the script now echoes
 * `/var/folders/.../1ALU-assembly1-pdbchains.cif` — a path that is not the
 * entry, names a rename the deposited file never had, and is unreadable from
 * the sandbox anyway. The dossier's `structural_neighbour_precedent.query_structure`
 * wants `1ALU`.
 *
 * The handler knows which entry it was called with, so it echoes that and moves
 * the script's path to `query_structure_file`, where it survives as provenance
 * and cannot be mistaken for an ID. A caller who passed a path gets that path
 * back verbatim in both places — there is no entry to recover, and inventing
 * one would be worse than the bug.
 */
function stampQueryStructure(doc: JsonObject, requested: string): void {
  const echoed = PDB_ID.test(requested) ? requested.toUpperCase() : requested;
  const scriptValue = doc.query_structure;
  doc.query_structure = echoed;
  if (typeof scriptValue === "string" && scriptValue !== echoed) {
    doc.query_structure_file = scriptValue;
  }
}

function renameNote(
  requested: string,
  renamed: Record<string, string>
): string {
  return (
    `${requested}'s biological assembly names its chains ` +
    `${Object.keys(renamed).join(", ")}, which the PDB format this script ` +
    "converts to cannot hold, so the local handler renamed them to " +
    `${Object.values(renamed).join(", ")} with gemmi.shorten_chain_names(). ` +
    "SAME ASSEMBLY, SAME COORDINATES, SAME CHAIN COUNT — only the labels " +
    "changed, and no reported number is keyed on them. If you quote a chain " +
    "id from this result, say it is the handler's relabelling of the " +
    "deposited id."
  );
}

/**
 * Render a successful `neighbour_precedent` run: stamp the entry into
 * `query_structure`, attach the chain-rename note when there was one, and put
 * the whole thing through the same never-truncate ladder every other handler
 * uses.
 *
 * The order matters. The stamp and the note are applied to the PARSED document
 * before it is measured, so a reduced payload is measured as it will actually
 * be sent — re-serialising a reduced document afterwards to add a note is how
 * a handler ships something over the cap while believing it fits.
 */
function renderNeighbourPayload(
  text: string,
  requested: string,
  renamed: Record<string, string> | null
): string {
  const lead = renamed ? renameNote(requested, renamed) : undefined;
  let parsed: JsonObject | undefined;
  try {
    parsed = asJsonObject(JSON.parse(text));
  } catch {
    parsed = undefined;
  }
  if (!parsed) {
    // Not the JSON document. Never fabricate one around it, and never cut it:
    // if it is also oversized, `fitText` refuses with the byte count.
    if (text.length > MAX_OUTPUT_CHARS) {
      return fitText(text, NEIGHBOUR_FIT);
    }
    return lead ? `${lead}\n\n${text}` : text;
  }
  stampQueryStructure(parsed, requested);
  return fitParsed({
    fullChars: text.length,
    lead,
    parsed,
    target: NEIGHBOUR_FIT,
  });
}

const crypticAnalysis: CustomToolSpec = {
  description:
    "Classify the cryptic-pocket mechanism from one apo/holo structure pair: core C-alpha superposition excluding the mobile region, max backbone displacement at the site, clash attribution split into backbone / side-chain / displaced-subunit, a self-control against the holo structure itself, and the ligand's free-volume fraction in the apo frame. " +
    "Run it whenever both an apo and a holo structure exist, because `cryptic_pocket_risk` must be measured rather than set from structure tier. " +
    "Classify on C-alpha displacement and NOT on which atoms clash: KRAS switch-II moves 8.8 A yet zero of its 12 clashing atoms at 2.0 A are backbone, so keying on clash composition labels the canonical nanomolar target as side-chain occlusion and hands it a micromolar prognosis. " +
    "Read the self-control first — it must come back near zero, and if it does not, the superposition or the ligand placement is broken and every other number in the result is meaningless. " +
    "PASS PDB IDs, NOT PATHS. This handler runs off-sandbox and cannot see a file you downloaded there; a 4-character ID has biological assembly 1 fetched for you on the machine that actually runs the script. A path is accepted only when it is readable by that machine.",
  async handler(input) {
    const argv = [
      await resolveStructure(requiredStr(input, "apo")),
      await resolveStructure(requiredStr(input, "holo")),
      requiredStr(input, "ligand_comp_id"),
    ];
    const exclude = numList(input, "exclude");
    if (exclude.length > 0) {
      argv.push("--exclude", ...exclude.map(String));
    }
    const holoChains = list(input, "holo_chains");
    if (holoChains.length > 0) {
      argv.push("--holo-chains", ...holoChains);
    }
    const apoChains = list(input, "apo_chains");
    if (apoChains.length > 0) {
      argv.push("--apo-chains", ...apoChains);
    }
    const ligandChain = str(input, "ligand_chain");
    if (ligandChain) {
      argv.push("--ligand-chain", ligandChain);
    }
    const excludeRadius = num(input, "exclude_radius");
    if (excludeRadius !== undefined) {
      argv.push("--exclude-radius", String(excludeRadius));
    }
    if (bool(input, "no_trim")) {
      argv.push("--no-trim");
    }
    if (bool(input, "no_free_volume")) {
      argv.push("--no-free-volume");
    }
    const result = await run(
      micromamba(),
      pythonArgv(join(SKILLS_DIR, "pocket-scan", "cryptic_analysis.py"), argv),
      { timeoutSeconds: DEFAULT_TIMEOUT_S }
    );
    return fitReport(result, CRYPTIC_FIT);
  },
  input_schema: {
    properties: {
      apo: {
        description:
          "Apo structure: a 4-character PDB ID (biological assembly 1 is fetched automatically on the machine running the script) or a path readable by that machine. A sandbox path is not readable here.",
        type: "string",
      },
      apo_chains: {
        description: "Chain IDs to keep from the apo structure.",
        items: { type: "string" },
        type: "array",
      },
      exclude: {
        description:
          "Residue numbers to exclude from the superposition fit by hand, when `exclude_radius` cannot express the mobile region.",
        items: { type: "integer" },
        type: "array",
      },
      exclude_radius: {
        description:
          "Radius around the site used to exclude mobile residues from the superposition fit.",
        type: "number",
      },
      holo: {
        description:
          "Holo structure: a 4-character PDB ID (biological assembly 1 is fetched automatically on the machine running the script) or a path readable by that machine. A sandbox path is not readable here.",
        type: "string",
      },
      holo_chains: {
        description: "Chain IDs to keep from the holo structure.",
        items: { type: "string" },
        type: "array",
      },
      ligand_chain: {
        description: "Chain the reference ligand sits in.",
        type: "string",
      },
      ligand_comp_id: {
        description:
          "Chemical component ID of the holo ligand, e.g. `MOV` or `307`. Use the component ID, not the paper's name for the compound.",
        type: "string",
      },
      no_free_volume: {
        description: "Skip the free-volume calculation.",
        type: "boolean",
      },
      no_trim: {
        description:
          "Disable auto-trim. The hand-calibration protocol used this; the deployed default does not, and the two differ by 0.1-0.2 A.",
        type: "boolean",
      },
    },
    required: ["apo", "holo", "ligand_comp_id"],
    type: "object",
  },
  name: "cryptic_analysis",
};

const interfaceAnalysis: CustomToolSpec = {
  description:
    "PPI-interface support for `tractability.pocket_vs_interface`: with `partners_accession` it lists the deposited complex structures that contain a binding partner for that UniProt accession, which is the input `pocket_scan`'s `partner_structures` argument needs. " +
    "Call it before `pocket_scan` whenever the mechanism hypothesis is orthosteric or oligomer-destabilisation, so the pocket-versus-interface label is measured rather than assumed. " +
    "With `selftest_dir` it instead runs the bundled fixture harness (IL-17A, TNF-alpha, KRAS) and caches CIFs there — use that only to check the module still behaves, never as evidence about a target. " +
    "Caveat: an asymmetric unit is not a biological assembly, and this module's own docstring records the case — 2AZ5's `assembly1` is a crystallographic tetramer of two independent TNF-alpha dimers, and scoring all four chains fuses sites across a packing contact, so always state which chains you meant.",
  async handler(input) {
    const partners = str(input, "partners_accession");
    const selftest = str(input, "selftest_dir");
    if (!(partners || selftest)) {
      throw new Error(
        "supply either partners_accession (to list partner complexes) or selftest_dir (to run the fixture harness)"
      );
    }
    const argv = partners
      ? ["--partners", partners]
      : ["--selftest", selftest ?? ""];
    const result = await run(
      micromamba(),
      pythonArgv(
        join(SKILLS_DIR, "pocket-scan", "interface_analysis.py"),
        argv
      ),
      { timeoutSeconds: DEFAULT_TIMEOUT_S }
    );
    return fitReport(result, INTERFACE_FIT);
  },
  input_schema: {
    properties: {
      partners_accession: {
        description:
          "UniProt accession whose partner-containing complex structures should be listed, e.g. `P01375`.",
        type: "string",
      },
      selftest_dir: {
        description:
          "Directory to cache CIFs in while running the fixture harness. Diagnostic only.",
        type: "string",
      },
    },
    required: [],
    type: "object",
  },
  name: "interface_analysis",
};

const disorderScan: CustomToolSpec = {
  description:
    "Predict the intrinsic-disorder fraction for one or more UniProt accessions and report, per target, the fraction, the method that produced it, a confidence flag and the disordered regions. " +
    "Use it to fill `tractability.disorder_fraction` when you are not already running `pocket_scan` (whose disorder stage is the preferred source, because the Modal image carries metapredict). " +
    "A disorder fraction is not comparable to any other disorder fraction unless you carry `method` beside it: the Modal image returned 0.3419 on a target where a local environment without metapredict fell back to MobiDB and returned 0.277 — a 23% difference from the method alone. " +
    "The cardinal rule of this module is that a folded protein and a failed prediction must never look identical: 0.000 is a real answer (CDK2 and KRAS both score it) and failure is reported as FAILED, so never read a missing number as zero. " +
    "IT TAKES ACCESSIONS AND ONLY ACCESSIONS, WHICH IS THE POINT. A disorder fraction measured on a crystallised construct is an answer to a different question, not a smaller version of the same one — a deposited construct is the ordered part of the protein by selection. IRAK4 measured on its construct gives 0.0 over 284 residues against 0.1413 over 460 full-length, and that 0.0 reads as `no disorder` rather than `not measured`. So supply `uniprot_accession` to `pocket_scan` too; if its output carries `construct_disorder_fraction` and `is_full_length_sequence: false` instead of `disorder_fraction`, the accession was missing and the number must not go into `tractability.disorder_fraction`.",
  async handler(input) {
    const accessions = list(input, "accessions");
    if (accessions.length === 0) {
      throw new Error(
        "accessions is required and must name at least one UniProt accession"
      );
    }
    const envName = process.env.DISORDER_ENV ?? analysisEnvName();
    const result = await run(
      micromamba(),
      [
        "run",
        "-n",
        envName,
        "python",
        join(SKILLS_DIR, "pocket-scan", "disorder.py"),
        ...accessions,
      ],
      { timeoutSeconds: DEFAULT_TIMEOUT_S }
    );
    return fitReport(result, DISORDER_FIT);
  },
  input_schema: {
    properties: {
      accessions: {
        description: 'UniProt accessions to score, e.g. `["P01116","P01375"]`.',
        items: { type: "string" },
        type: "array",
      },
    },
    required: ["accessions"],
    type: "object",
  },
  name: "disorder_scan",
};

// ---------------------------------------------------------------------------
// structure-select: Foldseek neighbour precedent
// ---------------------------------------------------------------------------

/**
 * `neighbour_precedent.py` merges a `.env` into the child environment and
 * raises FileNotFoundError when the path is absent — and its argparse default
 * is one contributor's absolute laptop path. Derive the repo's own `.env`
 * instead, and hand it an empty file rather than that default when there is
 * none, since PAPERCLIP_API_KEY already reaches the child through process.env.
 */
function envFileArg(): string {
  const repoEnv = join(repoRoot, ".env");
  if (existsSync(repoEnv)) {
    return repoEnv;
  }
  const staging = mkdtempSync(join(tmpdir(), "dossier-env-"));
  const empty = join(staging, "empty.env");
  writeFileSync(empty, "");
  return empty;
}

const neighbourPrecedent: CustomToolSpec = {
  description:
    "Structural-neighbour precedent: Foldseek the query structure, then ask Paperclip whether any neighbour fold has ever had a drug-like small molecule put into it, filling the dossier's `structural_neighbour_precedent` axis. " +
    "Use it to answer 'what other folds look like mine and has anyone drugged one', which is a different and much stronger question than Pfam family membership — TNF-alpha and IL-17A are both cytokines and share nothing mechanically. " +
    "`structure` takes a 4-character PDB ID — biological assembly 1 is fetched on the machine that runs the script, because a file you downloaded in the sandbox is not readable there — or a path that machine can already read. `accession` is excluded from its own results. " +
    "ZERO HOLO NEIGHBOURS IS A FINDING; A FAILED LOOKUP WEARING THAT COSTUME IS THE WORST CONFUSION AVAILABLE ON THIS AXIS. Holo is decided by `ligand_filter` from each component's SMILES graph, so read `neighbour_entry_summary.n_undetermined` and `undetermined_pdb_ids` — and per neighbour, `holo_determined` and `undetermined_ligands` — before writing `no drug-like ligand among the fold neighbours`. A component whose lookup failed carries `lookup_failed` and lands in `undetermined`, which is a THIRD tier and not apo. `n_holo = 0` is only a finding when `n_undetermined = 0` beside it. " +
    "Records without SMILES classify as `unknown`, and `unknown` is not `druglike`, so a record source with no SMILES silently renders every neighbour apo. The sources that carry it are RCSB REST `data.rcsb.org/rest/v1/core/chemcomp/<ID>`, Paperclip `pdb_v.chemcomps` (this script's own source) and the CCD ligand file; an entry's own mmCIF `_chem_comp` block does not. " +
    "Three Foldseek column caveats are already handled inside the script and you should not re-correct them: remote mode mislabels columns so `evalue` is really the probability and `bit_score` is really the E-value, a TM-score only exists via tmalign mode, and `target_id` is a filename-plus-title blob rather than an ID. The load-bearing caveat is in the output: ligands are attributed at entry level, so every holo count comes back twice — an entry-level upper bound and a single-protein-entry lower bound — and you must report the gap rather than picking one. " +
    "PASS THE 4-CHARACTER PDB ID AND NOTHING ELSE. The handler fetches biological assembly 1 and, when that assembly names its chains in a way the PDB format cannot hold (1ALU deposits `A` and `A-2`; 7NXZ deposits `AAA`), renames them with gemmi before the script's CIF-to-PDB conversion sees them. Both entries used to die on arrival with `chain name too long for the PDB format`. Same assembly, same coordinates, same chain count, different labels — and a run that was renamed says so in a `_handler_note` at the top of its result. " +
    "`query_structure` COMES BACK AS THE ENTRY YOU ASKED FOR, e.g. `1ALU`, and it is what `structural_neighbour_precedent.query_structure` takes verbatim. The script itself reports the file it was handed, which since the assembly fetch and the chain-rename cache is a temp path like `/var/folders/.../1ALU-assembly1-pdbchains.cif` — not an entry, and not readable from the sandbox — so this handler overwrites it with the requested ID and keeps the script's path beside it in `query_structure_file` as provenance. Never copy `query_structure_file` into the dossier. " +
    "TWO DIFFERENT FAILURES, ONLY ONE OF THEM IS RULE 13. A `ModuleNotFoundError` for `proto_tools` is unavailability of the whole axis: null `structural_neighbour_precedent`, record it in `not_found`, and never write it up as `no structural neighbours found`. A chain-name or file-format error is NOT that case — the tool is available and one input could not be read — so it throws with that distinction spelled out, and the answer is another entry for the same target (recorded in `not_found` as a substitution), not a nulled axis.",
  async handler(input) {
    const requested = requiredStr(input, "structure");
    // Assembly CIFs carry chain names the PDB format cannot hold, and
    // `count_chains` converts to PDB before the search runs. Normalise the
    // names here — same assembly, same coordinates — or 1ALU and 7NXZ die
    // before Foldseek is ever called. See pdbSafeChainNames.
    const structure = await pdbSafeChainNames(
      await resolveStructure(requested)
    );
    const argv = [
      structure.path,
      requiredStr(input, "accession"),
      "--env-file",
      envFileArg(),
    ];
    const multimer = str(input, "multimer");
    if (multimer) {
      argv.push("--multimer", multimer);
    }
    const maxNeighbours = num(input, "max_neighbours");
    if (maxNeighbours !== undefined) {
      argv.push("--max-neighbours", String(maxNeighbours));
    }
    const minAlignmentLength = num(input, "min_alignment_length");
    if (minAlignmentLength !== undefined) {
      argv.push("--min-alignment-length", String(minAlignmentLength));
    }
    const cache = str(input, "cache");
    if (cache) {
      argv.push("--cache", cache);
    }
    if (bool(input, "no_tm")) {
      argv.push("--no-tm");
    }
    requireEnv(
      "PAPERCLIP_API_KEY",
      "The neighbour lookup shells out to the Paperclip CLI for each neighbour's holo counts."
    );
    const result = await run(
      micromamba(),
      pythonArgv(
        join(SKILLS_DIR, "structure-select", "neighbour_precedent.py"),
        argv
      ),
      // The script shells out to `paperclip` itself for every neighbour's holo
      // counts, and those calls are outside this file's guard, so hand it the
      // same isolated config directory: a sticky `cli_cwd` would otherwise
      // kill the axis from inside a process this handler cannot inspect.
      { env: paperclipEnv(), timeoutSeconds: NEIGHBOUR_TIMEOUT_S }
    );
    if (CHAIN_NAME_TOO_LONG.test(`${result.stderr}\n${result.stdout}`)) {
      throw new Error(
        `${requested} still carries a chain name the PDB format cannot hold ` +
          "after this handler normalised the assembly's chain names, so the " +
          "Foldseek query never ran. THIS IS A FILE-FORMAT FAILURE, NOT THE " +
          "rule-13 CASE: rule 13 nulls this axis on a `proto_tools` " +
          "ModuleNotFoundError, which means the tool is unavailable. This is " +
          "an available tool that could not read one input. Do not null " +
          "`structural_neighbour_precedent` and do not write 'no structural " +
          "neighbours found'. Try another entry for the same target and " +
          "record the substitution in `not_found`, or hand the operator a " +
          "path to a file whose chains are single-character.\n\n" +
          `--- stderr ---\n${clip(result.stderr)}`
      );
    }
    if (result.code !== 0) {
      return report("neighbour_precedent.py", result);
    }
    return renderNeighbourPayload(
      result.stdout || result.stderr || "(no output)",
      requested,
      structure.renamed
    );
  },
  input_schema: {
    properties: {
      accession: {
        description:
          "The query's UniProt accession. Excluded from its own results.",
        type: "string",
      },
      cache: {
        description: "Path to cache Foldseek results in across calls.",
        type: "string",
      },
      max_neighbours: {
        description: "Cap on distinct neighbour PDB entries carried forward.",
        type: "number",
      },
      min_alignment_length: {
        description:
          "Override the alignment-length floor. Leave unset for the verified default of 120 with automatic relaxation on short queries.",
        type: "number",
      },
      multimer: {
        description:
          "`auto` (default) runs a multimer search when the file has more than one chain; `yes`/`no` force it. Set `no` when a fetched assembly's extra chains are symmetry copies you do not want matched.",
        enum: ["auto", "yes", "no"],
        type: "string",
      },
      no_tm: {
        description:
          "Skip the tmalign pass. Faster, but then no TM-score is available at all.",
        type: "boolean",
      },
      structure: {
        description:
          "The query: a 4-character PDB ID (biological assembly 1 is fetched on the machine running the script, and its chain names are normalised to single characters if the assembly uses `A-2` or `AAA` style ids) or a path that machine can read. A sandbox path is not readable here.",
        type: "string",
      },
    },
    required: ["structure", "accession"],
    type: "object",
  },
  name: "neighbour_precedent",
};

// ---------------------------------------------------------------------------
// Preflight
// ---------------------------------------------------------------------------

/**
 * Verify every credential and binary this agent needs, before the run starts.
 *
 * The per-handler checks above are already loud, but they fire at *use* time —
 * and `pocket_scan` is typically reached tens of minutes into a dossier, after
 * the precedent queries. A missing MODAL_BIN discovered there costs the whole
 * run. This runs at second zero instead.
 *
 * It deliberately collects *every* problem before throwing. Failing on the
 * first one turns "fix your setup" into a guess-and-recheck loop; one run
 * should tell you everything that is wrong.
 *
 * Set `DOSSIER_SKIP_PREFLIGHT=1` only to drive the Paperclip tools by hand on a
 * machine with no Modal — never for a real dossier, because the run will then
 * lose the computed-tractability axis at the point of use instead of here.
 */
/**
 * Ask the corpus one question it must be able to answer.
 *
 * Resolving the binary and finding the key in `.env` proves nothing about
 * whether Paperclip will answer — measured today, a machine that passed both
 * checks returned `vsh: cd: /papers/: Permission denied` at exit 0 for every
 * query in the run. That reads like an empty corpus, and an empty corpus is a
 * dossier with no retrieved-precedent axis. One 7 ms query at second zero
 * separates the two.
 */
async function checkPaperclipAnswers(problems: string[]): Promise<void> {
  let bin: string;
  try {
    bin = paperclipBin();
  } catch {
    return; // already reported by the binary check
  }
  if (!process.env.PAPERCLIP_API_KEY) {
    return; // already reported by the credential check
  }
  const ask = () =>
    run(
      bin,
      [
        "sql",
        "-s",
        "proteins",
        "SELECT accession FROM uniprot_v.proteins LIMIT 1",
      ],
      { env: paperclipEnv(), timeoutSeconds: 180 }
    );
  // Twice, and the second one is not superstition: the first query against a
  // freshly created config directory also fetches feature flags, and that
  // first call has been measured timing out on a corpus that answers the very
  // next call in 7 ms. Failing preflight on it would refuse a healthy run.
  let probe = await ask();
  if (paperclipFailure(probe) || !probe.stdout.includes("accession")) {
    probe = await ask();
  }
  const failure = paperclipFailure(probe);
  if (failure || !probe.stdout.includes("accession")) {
    problems.push(
      "  - paperclip liveness: a one-row probe against uniprot_v.proteins did " +
        `not come back with rows (${failure?.label ?? "unrecognised output"}), ` +
        "so the retrieved-precedent axis is not reachable and an empty answer " +
        "would read exactly like a target with no precedent. If this was a " +
        "transient timeout, re-run; the same probe answers in ~7 ms when the " +
        `corpus is healthy. (${clip(probe.stdout || probe.stderr)
          .slice(0, 200)
          .trim()})`
    );
  }
}

export async function preflight(): Promise<void> {
  if (process.env.DOSSIER_SKIP_PREFLIGHT === "1") {
    return;
  }
  const problems: string[] = [];
  const check = (label: string, probe: () => unknown) => {
    try {
      probe();
    } catch (error) {
      problems.push(`  - ${label}: ${(error as Error).message}`);
    }
  };

  check("ANTHROPIC_API_KEY", () =>
    requireEnv(
      "ANTHROPIC_API_KEY",
      "It authenticates this process to the Agents API; without it there is no session at all."
    )
  );
  check("PAPERCLIP_API_KEY", () =>
    requireEnv(
      "PAPERCLIP_API_KEY",
      "It carries the entire retrieved-precedent axis."
    )
  );
  check("paperclip binary", paperclipBin);
  check("micromamba binary", micromamba);
  await checkPaperclipAnswers(problems);
  check("modal binary", modalBin);
  check("modal profile", modalProfile);

  // The conda env is the one thing that cannot be checked without running it:
  // micromamba resolving is not the same as the env existing with gemmi+numpy.
  let micromambaPath: string | undefined;
  try {
    micromambaPath = micromamba();
  } catch {
    micromambaPath = undefined;
  }
  if (micromambaPath) {
    const env = analysisEnvName();
    const probe = await run(
      micromambaPath,
      ["run", "-n", env, "python", "-c", "import gemmi, numpy"],
      { timeoutSeconds: 120 }
    );
    if (probe.code !== 0) {
      problems.push(
        `  - DRUGGABILITY_ENV: micromamba env "${env}" is missing or lacks gemmi/numpy, ` +
          "so cryptic_analysis, interface_analysis and neighbour_precedent cannot run. " +
          "Create it with micromamba (never conda/mamba): `micromamba create -y " +
          `-n ${env} -c conda-forge python=3.14 fpocket gemmi numpy pip\`, then ` +
          "`pip install modal` and the pinned proto-tools into it — the full " +
          "recipe is managed/druggability-dossier/CREDENTIALS.md section 8.3. " +
          "Or set DRUGGABILITY_ENV to an env that already has them. " +
          `(${clip(probe.stderr).split("\n").slice(-3).join(" ").trim()})`
      );
    }
    // `proto_tools` is a WARNING, not a problem, and the asymmetry is
    // deliberate. CLAUDE.md rule 13 says a missing `proto_tools` nulls the
    // structural-neighbour axis with a stated reason, which is a legal dossier
    // — so failing preflight on it would refuse runs the agent is specified to
    // complete. But discovering it at the point of use looks identical to
    // "this fold has no neighbours", so it is worth saying at second zero.
    const proto = await run(
      micromambaPath,
      ["run", "-n", analysisEnvName(), "python", "-c", "import proto_tools"],
      { timeoutSeconds: 120 }
    );
    if (proto.code !== 0) {
      console.error(
        `druggability-dossier preflight WARNING: "${analysisEnvName()}" cannot ` +
          "import proto_tools, so `neighbour_precedent` will return a " +
          "ModuleNotFoundError. That is UNAVAILABILITY of the " +
          "structural_neighbour_precedent axis, which rule 13 says to null " +
          "with a reason in not_found — it must never be written up as 'no " +
          "structural neighbours found'. Not fatal; the run can complete."
      );
    }
  }

  if (problems.length > 0) {
    throw new Error(
      `druggability-dossier preflight failed — ${problems.length} problem(s) ` +
        "must be fixed before a run can produce a complete dossier. Each one " +
        "below would otherwise surface as a missing axis mid-run, which is " +
        "indistinguishable from a target with no data:\n" +
        `${problems.join("\n")}\n\n` +
        "See managed/druggability-dossier/CREDENTIALS.md."
    );
  }
}

export const tools: CustomToolSpec[] = [
  paperclipSql,
  paperclipSearch,
  paperclipGrep,
  paperclipRead,
  pocketScan,
  crypticAnalysis,
  interfaceAnalysis,
  disorderScan,
  neighbourPrecedent,
];
