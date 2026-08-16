#!/usr/bin/env python3
"""Assemble the target-intake bundle from an upstream evidence graph.

Mechanical only. Graph traversal, edge classification, evidence tiering from
each link's `basis`, and — where the verb is unrecognised — an adjudication
packet carrying every deterministic signal the graph offers. Reading a mechanism
out of a quote is judgment and stays with the agent; see SKILL.md.

Nothing here decides tractability, ranks targets, or picks an accession.

Usage:
    python3 graph_read.py <graph.json>
    python3 graph_read.py <graph.json> --thing t1
    python3 graph_read.py <graph.json> --allow-fixture
    python3 graph_read.py <graph.json> --ask-context
    python3 graph_read.py <graph.json> --check-ask '<ask json>'
"""

import argparse
import json
import re
import sys

# SCHEMA.md v1.1 gives `kind` six values. Both of these name a protein target;
# papers say "IRAK4 knockdown" as readily as "IRAK4 protein", so the extractor
# can legitimately type the same target either way.
TARGET_KINDS = {"protein", "gene"}

# Tiers that must not set chain selection. `hedged_only` is every finding saying
# "may" or "suggests"; `background_only` is every finding restating someone
# else's work. Both produce a confident-looking answer from evidence that has
# not asserted anything.
NON_ACTIONABLE_BASIS = {"background_only", "hedged_only"}

# `how` has NO enum in SCHEMA.md. Every other categorical field there carries an
# explicit a|b|c comment; `how` does not. It is open vocabulary written by the
# upstream extraction model, so these sets can never be complete. An unmatched
# verb is NOT dropped -- it goes to `needs_adjudication` with the signals below.
DIRECT_ACTION = {
    "inhibits", "binds", "blocks", "antagonises", "antagonizes",
    "agonises", "agonizes", "degrades", "stabilises", "stabilizes",
    "activates", "engages", "occupies", "targets", "modulates",
    "inactivates", "disrupts",
}

DOWNSTREAM_EFFECT = {
    "reduces", "increases", "improves", "worsens", "lowers", "raises",
    "suppresses", "restores", "prevents", "attenuates", "ameliorates",
    "induces", "normalises", "normalizes",
}

# Quote-level signals. These read the EVIDENCE, not the verb, so they survive an
# extractor that invents new relation words. Matched as whole words.
DIRECT_TERMS = [
    "ic50", "ec50", " ki ", " kd ", "affinity", "binds", "binding",
    "target engagement", "occupancy", "kinase activity", "enzymatic activity",
    "catalytic", "co-crystal", "cocrystal", "biochemical", "displacement",
    "atp-competitive", "allosteric", "active site",
]
DOWNSTREAM_TERMS = [
    "secretion", "release", "levels", "expression", "production", "output",
    "score", "response rate", "acr20", "acr50", "pasi", "serum", "plasma",
    "symptom", "endpoint", "placebo",
]

# `where` values that place a measurement in a direct-binding context.
DIRECT_CONTEXTS = ["biochemical", "cell-free", "cell free", "purified", "in vitro binding"]

NS_WORD = re.compile(r"[^a-z0-9]+")

# --- Second nomination route: symbols buried in entity names -----------------
#
# The kind-based route needs a `protein` or `gene` node. A real upstream graph
# may have none: "IRAK4 inhibition" is typed `small_molecule` because it names an
# intervention, and the protein exists only as a substring of that name.
#
# This route PROPOSES ONLY. A regex cannot know that a token is a gene, so every
# symbol below is emitted for verification against uniprot_v.proteins (SKILL.md
# step 4) and nothing here enters `nominations`. The script stays stdlib-only and
# offline by construction -- it must not call paperclip.

# One token, 2-10 chars, alphanumeric with internal hyphens.
SYMBOL_SHAPE = re.compile(r"^[A-Za-z][A-Za-z0-9]*(?:-[A-Za-z0-9]+)*$")

# Compound codes are symbol-shaped (ST2825, PF-06650833, KIC-0101). No human gene
# symbol carries a run of four digits, so that run separates the two.
COMPOUND_CODE = re.compile(r"\d{4,}")

# Separators that join DISTINCT symbols. Hyphen is NOT one of them: NF-kB, IL-6
# and IRAK-4 are single symbols that happen to carry a hyphen.
SYMBOL_SEPARATORS = re.compile(r"[/,+&]")

TRIM = " \t\"'()[]{}.,;:!?"

# The action word sitting next to a symbol seeds `interaction_to_disrupt`:
# "MyD88 dimerization inhibition" implies disrupting a dimerization interface,
# "IRAK4 inhibition" implies catalytic function. Carried verbatim -- turning it
# into a mechanism is the agent's call, not this script's.
ACTION_WORDS = {
    "inhibition", "inhibitor", "inhibitors", "inhibiting", "inhibited",
    "blockade", "blocking", "blocker", "block",
    "knockdown", "knockout", "silencing", "depletion", "ablation", "deletion",
    "degradation", "degrader", "degrading",
    "agonism", "agonist", "antagonism", "antagonist",
    "dimerization", "dimerisation", "oligomerization", "oligomerisation",
    "activation", "activator", "stabilization", "stabilisation",
    "engagement", "occupancy", "disruption", "suppression", "modulation",
    "deficiency", "loss",
}

# Stop-list: symbol-SHAPED tokens that are never gene symbols -- disease and
# tissue abbreviations, cell lines, reagents, assays, clinical endpoints, plus
# the action vocabulary itself. Matched case-insensitively, so an extractor that
# writes AXIS or SIGNALLING in caps still proposes nothing.
#
# Deliberately NOT here: TNF, TLR, IL-6 and friends. They are real symbols; the
# UniProt lookup, not this set, decides whether they resolve.
NOT_SYMBOLS = {
    "ra", "oa", "sle", "ibd", "copd", "gvhd", "as", "ms", "cd",
    "acr20", "acr50", "acr70", "das28", "pasi", "sdai", "cdai", "rct",
    "fls", "sf", "sfs", "pbmc", "pbmcs", "thp1", "thp", "hek", "hek293",
    "hela", "jurkat", "u937", "k562", "cho", "mcf7", "a549", "raw264",
    "bmdm", "huvec", "ipsc",
    "lps", "pma", "cfa", "atp", "adp", "gtp", "dmso", "pbs", "fbs",
    "dna", "rna", "mrna", "sirna", "shrna", "crispr", "cas9",
    "ic50", "ec50", "kd", "ki", "elisa", "facs", "pcr", "qpcr", "nmr",
    "hplc", "lcms", "msd", "spr", "itc", "auc", "cmax",
    "wt", "ko", "usa", "uk", "eu", "fda", "ema", "nih",
    "axis", "pathway", "signalling", "signaling", "inflammation", "disease",
    "protein", "kinase", "receptor", "complex", "cells", "cell",
}
NOT_SYMBOLS |= ACTION_WORDS

# The five lists SCHEMA.md guarantees, plus `rounds`. A missing list is not the
# same as an empty one and neither is `null`: `status: empty` with `things: []`
# is a real answer, a graph with no `things` key at all is a malformed file, and
# both used to produce the same silent zero-nomination output.
GRAPH_LISTS = ("things", "papers", "findings", "links", "gaps", "rounds")


class GraphShapeError(Exception):
    """A structural defect that must stop the run rather than degrade it."""


def _reject_constant(name):
    raise GraphShapeError(
        f"graph contains the JSON-invalid constant {name}. Python accepts it on "
        "read and re-emits it on write, so it would leave here as JSON no strict "
        "parser downstream can read."
    )


def read_list(graph, name, require_id=True):
    """Fetch one top-level list, refusing the shapes that silently zero us out.

    `rounds` rows are keyed by `n`, not `id`, so they are read with
    require_id=False.
    """
    if name not in graph:
        raise GraphShapeError(
            f"`{name}` is absent. SCHEMA.md v1.1 makes all five lists (plus "
            f"`rounds`) mandatory -- an exhausted search writes `{name}: []`, it "
            f"does not omit the key. An absent list and an empty one produce the "
            f"same zero nominations, so this is refused rather than reported."
        )
    val = graph[name]
    if val is None:
        raise GraphShapeError(
            f"`{name}` is null. SCHEMA.md note 7: failure is still a graph and "
            f"the lists are empty, never null. Treating null as empty would read "
            f"a broken producer as an absence of literature."
        )
    if not isinstance(val, list):
        raise GraphShapeError(
            f"`{name}` is {type(val).__name__}, not a list."
        )
    for i, row in enumerate(val):
        if not isinstance(row, dict):
            raise GraphShapeError(
                f"`{name}[{i}]` is {type(row).__name__}, not an object."
            )
        if require_id and "id" not in row:
            raise GraphShapeError(f"`{name}[{i}]` has no `id`.")
    return val


def index(graph):
    idx = {}
    dupes = []
    for key in ("things", "papers", "findings", "links", "gaps"):
        rows = read_list(graph, key)
        seen = {}
        for row in rows:
            rid = row["id"]
            if rid in seen:
                dupes.append({
                    "list": key,
                    "id": rid,
                    "kept": _identity(key, row),
                    "discarded": _identity(key, seen[rid]),
                })
            seen[rid] = row
        idx[key] = seen
    idx["_duplicate_ids"] = dupes
    return idx


def _identity(key, row):
    """Enough of a row to see that two rows sharing an id are not the same row."""
    if key == "things":
        return {"name": row.get("name"), "kind": row.get("kind")}
    if key == "links":
        return {"from": row.get("from"), "how": row.get("how"), "to": row.get("to")}
    if key == "findings":
        return {"quote": (row.get("quote") or "")[:80], "paper": row.get("paper")}
    if key == "papers":
        return {"title": (row.get("title") or "")[:80], "year": row.get("year")}
    return {"missing": row.get("missing")}


def integrity(graph, idx):
    """Every id that does not resolve, and every finding no link references.

    SCHEMA.md's guarantees section promises `from`, `to`, `paper` and the ids
    inside `yes`/`no`/`no_effect`/`implied_by` all resolve. That promise is not
    enforced upstream -- the one real graph we have (g_1a4f) carries
    `rounds[1].target: "g3"` against a gaps list holding g1, g2, g4. So the
    promise is checked here, and a dangling id is reported rather than dropped.
    Before this, a link whose `to` did not resolve vanished from every output:
    no nomination, no rejection, no adjudication, exit 0.
    """
    dangling, orphans = [], []
    things, papers = idx["things"], idx["papers"]
    findings, links, gaps = idx["findings"], idx["links"], idx["gaps"]

    def check(kind, rid, field, value, table, tname):
        if value is None or value in table:
            return
        dangling.append({"in": kind, "id": rid, "field": field,
                         "unresolved": value, "expected_in": tname})

    for lid, l in links.items():
        check("link", lid, "from", l.get("from"), things, "things")
        check("link", lid, "to", l.get("to"), things, "things")
        for arr in ("yes", "no", "no_effect"):
            for fid in l.get(arr) or []:
                check("link", lid, arr, fid, findings, "findings")

    for fid, f in findings.items():
        check("finding", fid, "paper", f.get("paper"), papers, "papers")
        # SCHEMA.md v1.1 puts from/how/to on findings as well as on links; the
        # real graph populates them. They are not read for classification, but a
        # dangling one is still a broken row.
        for field in ("from", "to"):
            if field in f:
                check("finding", fid, field, f.get(field), things, "things")

    for gid, gp in gaps.items():
        for tid in gp.get("missing") or []:
            check("gap", gid, "missing", tid, things, "things")
        for lid in gp.get("implied_by") or []:
            check("gap", gid, "implied_by", lid, links, "links")

    reachable = {"expand_node": things, "resolve_link": links, "test_gap": gaps}
    for r in read_list(graph, "rounds", require_id=False):
        tgt, verb = r.get("target"), r.get("ask")
        if tgt and verb in reachable and tgt not in reachable[verb]:
            dangling.append({"in": "round", "id": r.get("n"), "field": "target",
                             "unresolved": tgt, "expected_in": verb + " index"})

    cited = {fid for l in links.values() for arr in ("yes", "no", "no_effect")
             for fid in (l.get(arr) or [])}
    for fid in sorted(set(findings) - cited):
        orphans.append({"finding": fid, "paper": findings[fid].get("paper"),
                        "quote": (findings[fid].get("quote") or "")[:120]})

    return {
        "dangling_ids": dangling,
        "duplicate_ids": idx["_duplicate_ids"],
        "findings_referenced_by_no_link": orphans,
        "note": (
            "SCHEMA.md guarantees every id resolves. It is not enforced upstream "
            "-- g_1a4f ships a dangling rounds[].target. A dangling id used to be "
            "dropped silently; a duplicate id used to overwrite last-wins, which "
            "relabels a whole evidence neighbourhood with another protein's name. "
            "Neither is a parse error, so neither stops the run -- but a "
            "nomination touching one of these rows is not trustworthy."
        ) if (dangling or idx["_duplicate_ids"] or orphans) else None,
    }


def classify(link, things):
    """Edge class depends on the SUBJECT's kind, not on the verb alone.

    `activates` from a small molecule is an agonist; `activates` from a receptor
    is pathway biology. Same verb, different class.
    """
    subject = things.get(link["from"], {})
    verb = (link.get("how") or "").lower().strip()

    if subject.get("kind") != "small_molecule":
        return "biological_relation"
    if verb in DIRECT_ACTION:
        return "direct_action"
    if verb in DOWNSTREAM_EFFECT:
        return "downstream_effect"
    return "unclassified"


def paper_ref(paper):
    if not paper:
        return None
    study = (paper.get("study_type") or "").replace("_", " ")
    ref = f"{paper.get('first_author')} et al., {paper.get('journal')} {paper.get('year')}"
    return f"{ref} ({study})" if study else ref


def expand(finding_id, link, idx):
    f = idx["findings"].get(finding_id)
    if not f:
        return None
    things = idx["things"]
    return {
        "finding": finding_id,
        "link": link["id"],
        "relation": (
            f"{things.get(link['from'], {}).get('name')} {link.get('how')} "
            f"{things.get(link['to'], {}).get('name')}"
        ),
        "quote": f.get("quote"),
        "where": f.get("where"),
        "section": f.get("section"),
        "says": f.get("says"),
        "paper": f.get("paper"),
        "paper_ref": paper_ref(idx["papers"].get(f.get("paper"))),
        "retracted": (idx["papers"].get(f.get("paper")) or {}).get("retracted"),
        "is_own_result": f.get("is_own_result"),
        "hedged": f.get("hedged"),
        "finding_confidence": f.get("confidence"),
        "flags": f.get("flags", []),
        "link_basis": link.get("basis"),
        "link_state": link.get("state"),
        "link_confidence": (link.get("confidence") or {}).get("overall"),
    }


def link_findings(link):
    return (
        list(link.get("yes", []))
        + list(link.get("no", []))
        + list(link.get("no_effect", []))
    )


def matched(text, terms):
    if not text:
        return []
    padded = " " + NS_WORD.sub(" ", text.lower()) + " "
    return [t.strip() for t in terms if t.strip() and (" " + t.strip() + " ") in padded
            or (" " not in t.strip() and " " + t.strip() + " " in padded)]


def signals(link, idx):
    """Deterministic evidence for a target-vs-readout call, read from the
    quotes and the graph shape rather than from the verb string.

    Returned as evidence, never as a verdict. The agent adjudicates.
    """
    things = idx["things"]
    obj_id = link.get("to")
    rows = [expand(f, link, idx) for f in link_findings(link)]
    rows = [r for r in rows if r]

    quotes = " ".join(r["quote"] or "" for r in rows)
    wheres = [r["where"] for r in rows if r.get("where")]

    # A readout usually carries the causal chain onward to a disease. A target
    # usually does not. Weak alone, useful alongside the quote terms.
    onward_to_disease = [
        l["id"] for l in idx["links"].values()
        if l.get("from") == obj_id
        and things.get(l.get("to"), {}).get("kind") == "disease"
    ]

    return {
        "object_kind": things.get(obj_id, {}).get("kind"),
        "object_has_edge_to_disease": onward_to_disease,
        "assay_contexts": wheres,
        "direct_context": [
            w for w in wheres if any(c in (w or "").lower() for c in DIRECT_CONTEXTS)
        ],
        "direct_terms_in_quotes": matched(quotes, DIRECT_TERMS),
        "downstream_terms_in_quotes": matched(quotes, DOWNSTREAM_TERMS),
    }


def neighbourhood(thing_id, idx):
    """Every link touching this thing, with its findings expanded and tiered.

    Tier comes from the LINK's `basis`, not from the finding's own confidence.
    A 0.88-confidence quote from a review is still background.
    """
    tiers = {"primary": [], "mixed": [], "background_only": []}
    for link in idx["links"].values():
        if thing_id not in (link.get("from"), link.get("to")):
            continue
        basis = link.get("basis") or "unknown"
        bucket = tiers.setdefault(basis, [])
        # Three arrays, not two. `no_effect` is a measured null result, which is
        # not the same as `no` (evidence against) -- on a direct-action edge it
        # is real tractability evidence that the compound does not engage.
        for fid in link_findings(link):
            row = expand(fid, link, idx)
            if row:
                row["actionable"] = basis not in NON_ACTIONABLE_BASIS
                bucket.append(row)
    return {k: v for k, v in tiers.items() if v}


def signal_reading(sig):
    """Which way the QUOTES point, independent of the verb. Mechanical.

    Returns `direct`, `downstream` or `silent`. This is not a verdict and it does
    not override anything -- it exists so that a verb and its own evidence
    disagreeing becomes visible instead of being decided by the verb alone.

    `silent` is the common case and it is deliberately not an answer: a blank
    means the field was blank far more often than it means the experiment was
    absent (see SKILL.md's signals table).
    """
    direct = bool(sig["direct_terms_in_quotes"]) or bool(sig["direct_context"])
    downstream = bool(sig["downstream_terms_in_quotes"])
    if direct and not downstream:
        return "direct"
    if downstream and not direct:
        return "downstream"
    return "silent"


def contested_packet(link, idx, cls, reading, obj, things):
    return {
        "link": link["id"],
        "how": link.get("how"),
        "why_contested": (
            f"verb {link.get('how')!r} classifies as {cls}, but the quotes on this "
            f"edge read as {reading}. The verb and its own evidence disagree."
        ),
        "verb_class": cls,
        "quote_reading": reading,
        "subject": {"id": link.get("from"),
                    "name": things.get(link.get("from"), {}).get("name")},
        "object": {"id": link.get("to"), "name": obj.get("name"),
                   "kind": obj.get("kind")},
        "eligible_kind": obj.get("kind") in TARGET_KINDS,
        "signals": signals(link, idx),
        "findings": [r for r in (expand(f, link, idx)
                                 for f in link_findings(link)) if r],
        "decide": (
            "Is this a direct action on a target, or a downstream effect on a "
            "readout? SKILL.md 'Adjudicating an unknown verb' applies unchanged -- "
            "read the quotes first, the signals only as corroboration. Refusing is "
            "allowed; guessing is not."
        ),
    }


def nominate(graph, idx):
    """A thing is a target candidate if its kind is protein or gene AND either

      (a) the object of a direct-action edge from a small molecule, or
      (b) named in a gap.

    (b) is what carries the undrugged candidates. Without it the intake can only
    ever return targets somebody already made a molecule against.

    Two things this deliberately does NOT do any more:

    A direct-action verb no longer nominates on its own when the quotes on that
    same edge read downstream, and a downstream verb no longer rejects on its own
    when the quotes read direct. `classify()` is a pure function of (subject
    kind, verb) and the quote never reached it, so `blocks IL-6 secretion`
    nominated a cytokine as a target and `suppresses IRAK4 catalytic activity,
    Ki 1.1 nM` rejected a kinase as a readout -- both silently. Those now go to
    `needs_adjudication` where the existing procedure decides them.

    A gap no longer erases a rejection. Upstream `gaps` are structural open
    triangles (assemble.py `find_gaps`: A-B and B-C exist, A-C does not), not a
    curated list of undrugged candidates -- run their generator over our own RA
    fixture and it emits `missing: [IRAK4, IL-6]`, because both connect to
    zimlovisertib. The old `rejected.pop(tid)` then deleted IL-6's "readout, not
    target" reason and nominated it. A thing that is both is now reported as both.
    """
    things = idx["things"]
    nominated, rejected, adjudicate = {}, {}, []

    for link in idx["links"].values():
        cls = classify(link, things)
        obj = things.get(link.get("to"), {})
        if not obj:
            # `to` does not resolve. Reported by integrity(); not silently gone.
            continue
        reading, sig = None, None
        if cls in ("direct_action", "downstream_effect"):
            sig = signals(link, idx)
            reading = signal_reading(sig)
            if (cls == "direct_action" and reading == "downstream") or \
               (cls == "downstream_effect" and reading == "direct"):
                adjudicate.append(
                    contested_packet(link, idx, cls, reading, obj, things))
                continue
        if cls == "direct_action":
            if obj.get("kind") not in TARGET_KINDS:
                rejected.setdefault(obj["id"], []).append(
                    f"direct-action edge {link['id']} ({link.get('how')}) but kind is "
                    f"'{obj.get('kind')}', not one of {sorted(TARGET_KINDS)}"
                )
                continue
            # `says` has three values and a link has three arrays to match. A
            # direct-action edge whose findings ALL say `no` is the literature
            # reporting that the compound does not act on this thing, and it used
            # to produce a nomination indistinguishable from a positive one.
            support = {
                "yes": len(link.get("yes") or []),
                "no": len(link.get("no") or []),
                "no_effect": len(link.get("no_effect") or []),
            }
            nominated.setdefault(obj["id"], []).append({
                "via": link["id"],
                "evidence_class": "direct_action_edge",
                "support": support,
                # The quote reading travels WITH the nomination, including when
                # it is `silent`. A silent reading is the residual boundary: the
                # quote carries no vocabulary either way, so nothing but the verb
                # supports this nomination. `targets IL-6 driven inflammation in
                # synovium` reads silent and nominates a cytokine -- the verb is
                # doing all the work and the agent has to know that.
                "quote_reading": reading,
                "quote_reading_note": (
                    "Quotes corroborate a direct action." if reading == "direct" else
                    "NO quote on this edge carries direct-binding vocabulary. The "
                    "verb alone put this thing forward as a target. Read the quotes "
                    "(SKILL.md step 2) before accepting it -- this is the shape "
                    "that turns a readout into a target."
                ),
                "signals": sig,
                "direction_warning": (
                    None if support["yes"] else
                    f"edge {link['id']} carries NO `yes` findings "
                    f"({support['no']} no, {support['no_effect']} no_effect). "
                    "SCHEMA.md: a negation lives in `says`, never in the verb, so "
                    f"'{link.get('how')}' here states the opposite of what it reads "
                    "like. A no_effect on a direct-action edge is real tractability "
                    "evidence; a `no` is evidence against the relation itself."
                ),
                "why": (
                    f"object of direct-action edge from small_molecule "
                    f"{link['from']} ({link.get('how')})"
                ),
            })
        elif cls == "downstream_effect":
            rejected.setdefault(obj["id"], []).append(
                f"reached only by downstream-effect edge {link['id']} "
                f"({link.get('how')}) -- readout, not target"
            )
        elif cls == "unclassified":
            # NOT dropped. `how` is open vocabulary, so an unmatched verb is a
            # target the intake could not classify -- a decision to make, not a
            # rare edge to ignore.
            adjudicate.append({
                "link": link["id"],
                "how": link.get("how"),
                "subject": {"id": link.get("from"),
                            "name": things.get(link.get("from"), {}).get("name")},
                "object": {"id": link.get("to"), "name": obj.get("name"),
                           "kind": obj.get("kind")},
                "eligible_kind": obj.get("kind") in TARGET_KINDS,
                "signals": signals(link, idx),
                "findings": [r for r in (expand(f, link, idx)
                                         for f in link_findings(link)) if r],
                "decide": (
                    "Is this a direct action on a target, or a downstream effect on a "
                    "readout? See SKILL.md 'Adjudicating an unknown verb'. Refusing is "
                    "allowed; guessing is not."
                ),
            })

    for gap in idx["gaps"].values():
        for tid in gap.get("missing", []):
            thing = things.get(tid, {})
            if not thing:
                continue
            if thing.get("kind") not in TARGET_KINDS:
                rejected.setdefault(tid, []).append(
                    f"named in gap {gap['id']} but kind is '{thing.get('kind')}', "
                    f"not one of {sorted(TARGET_KINDS)}"
                )
                continue
            nominated.setdefault(tid, []).append({
                "via": gap["id"],
                "evidence_class": "structural_gap_only",
                "why": (
                    f"named in gap {gap['id']} (missing pair "
                    f"{gap.get('missing')}). Upstream gaps are open triangles, "
                    "computed from graph shape, NOT curated undrugged candidates: "
                    "this thing shares a neighbour with the other half of the pair "
                    "and has no edge to it. It carries no evidence of its own."
                ),
            })

    # A gap nomination no longer deletes a rejection. Both are reported, and the
    # collision is named, because "this was rejected as a readout AND nominated
    # by a gap" is exactly the case the old pop() hid.
    contested_nominations = {}
    for tid in nominated:
        if tid in rejected:
            contested_nominations[tid] = rejected[tid]

    return nominated, rejected, adjudicate, contested_nominations


def coverage_notes(graph):
    """Everything in `coverage` that bears on reading an absence as a result.

    Three things used to be invisible here. An ABSENT `coverage` block produced
    no warning at all, so a graph missing the block read exactly like a
    `complete` one. `no_quote_discarded` -- the count of claims upstream threw
    away for having no verbatim sentence, and its only documented removal -- was
    never surfaced. And `depth: quick` was never surfaced, though SCHEMA.md note
    2 says quick reads page 1 and page 1 lies, so at quick an absence means
    nothing even when stop_reason is `complete`.
    """
    if "coverage" not in graph:
        return {
            "coverage_present": False,
            "note": (
                "`coverage` is ABSENT. SCHEMA.md: 'what was NOT read. Always "
                "present.' With no coverage block there is no basis on which to "
                "read any absence in this graph as a result -- treat every "
                "unstated thing as unknown, not as absent."
            ),
        }
    coverage = graph["coverage"]
    if coverage is None or not isinstance(coverage, dict):
        raise GraphShapeError(
            f"`coverage` is {type(coverage).__name__}, not an object. SCHEMA.md "
            "makes it always present and always real, including on a failed graph."
        )

    stop = coverage.get("stop_reason")
    depth = coverage.get("depth")
    discarded = coverage.get("no_quote_discarded")
    warnings = []

    if coverage.get("truncated"):
        warnings.append("`truncated` is true: this is a sample, not the literature.")
    if stop is None:
        warnings.append("`stop_reason` is absent -- cannot tell budget from exhaustion.")
    elif stop != "complete":
        warnings.append(
            f"`stop_reason` is {stop!r}. Only 'complete' means the literature was "
            "exhausted; the other four mean the run ran out of budget, so an "
            "absent mechanism statement here is a budget limit, not an "
            "established absence (SCHEMA.md note 6)."
        )
    if depth == "quick":
        warnings.append(
            "`depth` is 'quick'. SCHEMA.md note 2: quick reads page 1, and page 1 "
            "lies -- at this depth absence means unknown, whatever stop_reason says."
        )
    if discarded:
        warnings.append(
            f"`no_quote_discarded` is {discarded}: upstream dropped {discarded} "
            "extracted claim(s) for carrying no verbatim sentence. That is the "
            "pipeline's only documented removal, it is silent in the graph, and "
            "those claims are not recoverable from this file. Anything this "
            "intake reports as unstated may be among them."
        )

    return {
        "coverage_present": True,
        "truncated": coverage.get("truncated"),
        "stop_reason": stop,
        "depth": depth,
        "no_quote_discarded": discarded,
        "literature_exhausted": stop == "complete" and depth != "quick"
                                and not coverage.get("truncated"),
        "warnings": warnings,
        "note": None if warnings else
                "coverage is clean: untruncated, stop_reason 'complete', nothing "
                "discarded. An absence in this graph may be read as an absence.",
    }


def selection_note(out):
    """The dossier's Contract takes ONE `uniprot_accession`. The graph can
    nominate any number and offers no basis to rank them -- and this skill does
    not rank (SKILL.md 'What this skill does not do'). Zero and many are both
    real outcomes and both used to look like a plain list.
    """
    n = len(out)
    with_evidence = [
        o["thing"] for o in out
        if any(r.get("evidence_class") == "direct_action_edge"
               for r in o["nominated_by"])
    ]
    gap_only = [o["thing"] for o in out if o["thing"] not in with_evidence]
    if n == 0:
        note = ("No candidate. Read `rejected`, `needs_adjudication` and "
                "`integrity` before reporting this as 'no targets in the "
                "literature' -- a granularity mismatch, a dangling id or an "
                "unknown verb all land here too.")
    elif n == 1:
        note = "One candidate. No choice to make."
    else:
        note = (
            f"{n} candidates and NO basis in the graph to choose between them. "
            "The dossier Contract takes one `uniprot_accession`, and this skill "
            "does not rank. Do not take the first: the order here is by `thing` "
            "id, which is insertion order upstream and carries no meaning. Either "
            "run the dossier once per candidate, or state the ambiguity and issue "
            "an ask. Picking one silently is the failure this field exists to stop."
        )
    return {
        "n_candidates": n,
        "basis_to_choose": None,
        "with_direct_action_evidence": with_evidence,
        "gap_only_no_evidence": gap_only,
        "note": note,
    }


def symbol_shaped(token):
    """Shape test only. Says nothing about whether the token names a gene."""
    if not (2 <= len(token) <= 10):
        return False
    if not SYMBOL_SHAPE.match(token):
        return False
    if COMPOUND_CODE.search(token):
        return False
    # Two capitals is the floor. One is ordinary prose capitalisation (Rho,
    # Toll-like, Matrigel); symbols carry their case (MyD88, NF-kB, IRAK4).
    if sum(1 for c in token if c.isupper()) < 2:
        return False
    return token.lower() not in NOT_SYMBOLS


def symbol_key(symbol):
    """Dedupe key only. IRAK-4 and IRAK4 are one candidate; the spelling the
    extractor used is kept verbatim on every mention."""
    return symbol.upper().replace("-", "")


def query_forms(symbol):
    """Spellings to put in the SQL `IN` list, in order. Not a rewrite of the
    symbol -- each form is a separate lookup that may return nothing."""
    forms = []
    for form in (symbol, symbol.upper(), symbol.upper().replace("-", "")):
        if form not in forms:
            forms.append(form)
    return forms


def action_near(tokens, i):
    """The action word adjacent to tokens[i], as (text, position).

    Suffix first ("MyD88 dimerization inhibition"), then prefix, optionally
    across "of" ("knockdown of MYD88"). Contiguous runs only, so an action word
    elsewhere in the phrase is not attached to this symbol.
    """
    after = []
    j = i + 1
    while j < len(tokens) and tokens[j].strip(TRIM).lower() in ACTION_WORDS:
        after.append(tokens[j].strip(TRIM))
        j += 1
    if after:
        return " ".join(after), "suffix"

    j = i - 1
    if j >= 0 and tokens[j].strip(TRIM).lower() == "of":
        j -= 1
    before = []
    while j >= 0 and tokens[j].strip(TRIM).lower() in ACTION_WORDS:
        before.insert(0, tokens[j].strip(TRIM))
        j -= 1
    if before:
        return " ".join(before), "prefix"
    return None, None


def scan_phrase(phrase):
    """Every symbol in one verbatim string, with its adjacent action word.

    A multi-symbol phrase returns EVERY symbol. "TLR/MyD88/NF-kB signalling
    axis" is three candidates; collapsing it to one is the failure this route
    exists to avoid.
    """
    tokens = phrase.split()
    hits, seen = [], set()
    for i, token in enumerate(tokens):
        for part in SYMBOL_SEPARATORS.split(token):
            part = part.strip(TRIM)
            if not symbol_shaped(part) or symbol_key(part) in seen:
                continue
            seen.add(symbol_key(part))
            action, position = action_near(tokens, i)
            hits.append({"symbol": part, "action": action,
                         "action_position": position})

    spellings = [h["symbol"] for h in hits]
    for hit in hits:
        hit["co_occurring_symbols"] = [
            s for s in spellings if symbol_key(s) != symbol_key(hit["symbol"])
        ]
    return hits


def thing_symbols(thing):
    """Scan `name` and every alias, whatever the thing's `kind`. Kind is exactly
    what this route cannot trust -- the target may be typed small_molecule."""
    fields = [("name", thing.get("name"))]
    fields += [("aliases[%d]" % n, a)
               for n, a in enumerate(thing.get("aliases") or [])]

    found, order = {}, []
    for field, phrase in fields:
        if not phrase:
            continue
        for hit in scan_phrase(phrase):
            mention = {
                "as_written": hit["symbol"],
                "field": field,
                # Whole-field means the extractor gave the symbol; parsed-out
                # means this regex inferred it from a longer phrase.
                "whole_field": hit["symbol"] == phrase.strip(TRIM),
                "phrase": phrase,
                "action": hit["action"],
                "action_position": hit["action_position"],
                "co_occurring_symbols": hit["co_occurring_symbols"],
            }
            key = symbol_key(hit["symbol"])
            if key not in found:
                found[key] = {"symbol": hit["symbol"], "mentions": []}
                order.append(key)
            found[key]["mentions"].append(mention)
    return [found[k] for k in order]


def symbol_candidates(idx, nominated_ids, only=None):
    """Symbols proposed from entity names, for UniProt verification.

    Never a nomination and never asserted. A candidate that resolves to no row
    is the lookup answering -- "NF-kB" names a complex, not a gene, so it is
    EXPECTED to fail and that failure is the result, not an error.
    """
    candidates = []
    for thing in idx["things"].values():
        if only and thing["id"] != only:
            continue
        for entry in thing_symbols(thing):
            mentions = entry["mentions"]
            # An action word is the point of this route, so a mention carrying
            # one leads even if a bare mention came first.
            lead = next((m for m in mentions if m["action"]), mentions[0])
            co = []
            for m in mentions:
                for s in m["co_occurring_symbols"]:
                    if symbol_key(s) not in {symbol_key(x) for x in co}:
                        co.append(s)
            candidates.append({
                "symbol": entry["symbol"],
                "query_forms": query_forms(entry["symbol"]),
                "action": lead["action"],
                "action_position": lead["action_position"],
                "thing": thing["id"],
                "thing_kind": thing.get("kind"),
                "thing_name": thing.get("name"),
                "already_nominated": thing["id"] in nominated_ids,
                "field": lead["field"],
                "whole_field": lead["whole_field"],
                "phrase": lead["phrase"],
                # True whenever the symbol shared a phrase with another symbol.
                # The agent resolves which one the dossier is about; picking one
                # here would be a guess.
                "ambiguous": bool(co),
                "co_occurring_symbols": co,
                "other_mentions": [m for m in mentions if m is not lead],
                # Filled by the agent from the SQL. Left null on purpose.
                "verified": None,
                "uniprot_accession": None,
            })

    # Symbols scraped from ONE thing's name + aliases are alternative names for
    # what the extractor thought was one entity, so verification must confirm
    # they agree on an accession -- not stop at the first that resolves.
    #
    # Measured on the modality-trap fixture: t1 yields TL1A, TNFSF15, TNF and
    # VEGI. TL1A and VEGI return no row, TNFSF15 returns O95150, and TNF -- a
    # substring of the descriptive alias "TNF ligand superfamily member 15" --
    # returns P01375, a different and heavily drugged protein. A resolver that
    # takes the first symbol that verifies picks TNF and hands the dossier an
    # accession with thousands of small-molecule actives against it. That is
    # failure mode 3 with a clean-looking row behind it.
    per_thing = {}
    for c in candidates:
        per_thing.setdefault(c["thing"], []).append(c["symbol"])
    for c in candidates:
        siblings = [x for x in per_thing[c["thing"]] if x != c["symbol"]]
        c["sibling_symbols"] = siblings
        # Sharing a phrase means different proteins; sharing a thing means
        # possibly-different names for one. Both block a unilateral pick, for
        # different reasons, so keep them as separate fields.
        c["needs_agreement_check"] = bool(siblings)

    return {
        "note": (
            "PROPOSED, NOT CONFIRMED. Regex over thing `name` and `aliases`, run "
            "because a graph can carry its proteins only inside intervention "
            "names -- 'IRAK4 inhibition' is typed small_molecule. Nothing here is "
            "a nomination. Verify every symbol against uniprot_v.proteins before "
            "using it; a symbol naming a complex rather than a gene (NF-kB) is "
            "expected to return no row, and that is the answer, not a failure."
        ),
        "verify_with": (
            "SELECT accession, gene_name, protein_name, organism, sequence_length "
            "FROM uniprot_v.proteins WHERE gene_name IN (<query_forms>) "
            "AND organism = 'Homo sapiens'"
        ),
        "agreement_rule": (
            "Symbols from one thing are alternative names for one entity, so "
            "verify ALL of them and require they agree on an accession. Two "
            "symbols off the same thing resolving to DIFFERENT accessions is an "
            "unresolved conflict, never a majority vote -- one of them is a "
            "substring of a descriptive alias. Do not stop at the first that "
            "resolves."
        ),
        "ambiguous_things": sorted({c["thing"] for c in candidates
                                    if c["ambiguous"]}),
        "things_needing_agreement": sorted({c["thing"] for c in candidates
                                            if c["needs_agreement_check"]}),
        "candidates": candidates,
    }

def build(graph, only=None):
    idx = index(graph)
    things = idx["things"]
    nominated, rejected, adjudicate, contested = nominate(graph, idx)
    # Captured before the --thing filter: `already_nominated` reports the graph,
    # not the slice being printed.
    nominated_ids = set(nominated)

    if only:
        if only not in nominated:
            raise GraphShapeError(
                f"--thing {only!r} is not a nomination. Nominated: "
                f"{sorted(nominated) or '(none)'}. Filtering to a thing that was "
                "never nominated returns an empty result that reads like a finding."
            )
        nominated = {k: v for k, v in nominated.items() if k == only}

    out = []
    for tid, reasons in sorted(nominated.items()):
        t = things[tid]
        out.append({
            "thing": tid,
            "name": t.get("name"),
            "kind": t.get("kind"),
            "aliases": t.get("aliases", []),
            "mentions": t.get("mentions"),
            "nominated_by": reasons,
            # Populated when this same thing also carries a rejection. The old
            # code deleted the rejection; a readout named in a structural gap
            # then became an unqualified target.
            "also_rejected_because": contested.get(tid),
            # Filled by the agent. Left null on purpose -- see SKILL.md.
            "gene_symbol": None,
            "uniprot_accession": None,
            "ambiguity": None,
            "interaction_to_disrupt": None,
            "mechanism_hypothesis": None,
            "evidence": neighbourhood(tid, idx),
        })

    coverage = graph.get("coverage") or {}
    status = graph.get("status")
    # `complete` is the ONLY stop_reason meaning the literature was exhausted.
    # The other four mean the run ran out of budget (SCHEMA.md note 6).
    stop_reason = coverage.get("stop_reason")
    retracted = [p["id"] for p in idx["papers"].values() if p.get("retracted")]

    return {
        "integrity": integrity(graph, idx),
        "coverage": coverage_notes(graph),
        "selection": selection_note(out),
        "graph_id": graph.get("graph_id"),
        "round": graph.get("round"),
        "question": graph.get("question"),
        # `status` is never an error blob -- an `empty` or `failed` graph parses
        # fine and yields zero nominations, which reads as "no targets found".
        "status": status,
        "status_warning": (
            None if status == "ok"
            else f"graph status is '{status}' -- lists may be empty for reasons that "
                 f"are not evidence. Do not read zero nominations as a result."
        ),
        # Kept for callers that read it. `coverage` above is the fuller block and
        # is the one to read: this one is null on an ABSENT coverage key, which
        # was the bug.
        "coverage_warning": (
            {
                "truncated": coverage.get("truncated"),
                "stop_reason": stop_reason,
                "depth": coverage.get("depth"),
                "note": (
                    "Only stop_reason 'complete' means the literature was exhausted. "
                    "An absent mechanism statement here is a budget limit, not an "
                    "established absence."
                ),
            }
            if coverage.get("truncated") or (stop_reason and stop_reason != "complete")
            else None
        ),
        "retracted_papers": retracted,
        "nominations": out,
        "rejected": [
            {"thing": tid, "name": things.get(tid, {}).get("name"), "why": why,
             "also_nominated": tid in contested}
            for tid, why in sorted(rejected.items())
        ],
        "contested_nominations": [
            {"thing": tid, "name": things.get(tid, {}).get("name"),
             "nominated_by": [r["via"] for r in nominated.get(tid, [])],
             "rejected_because": why,
             "note": (
                 "This thing is BOTH nominated and rejected. Most often: reached "
                 "by a downstream-effect edge (readout) and separately named in a "
                 "structural gap. The gap is not evidence -- it is an open "
                 "triangle. Decide before handing this to the dossier."
             )}
            for tid, why in sorted(contested.items())
        ],
        "needs_adjudication": adjudicate,
        # Second route. Independent of `nominations` -- it neither adds to nor
        # subtracts from them, and a graph with entity nodes still nominates on
        # kind exactly as before.
        "symbol_candidates": symbol_candidates(idx, nominated_ids, only),
    }


# ---------------------------------------------------------------------------
# Ask-back gating. MECHANICAL ONLY.
#
# SCHEMA.md defines FOUR ask verbs and nothing here adds a fifth. What this
# section does is refuse the asks that are mechanically wrong -- malformed,
# unactionable, already issued, or pointed at a graph that has nothing left to
# give. It CANNOT check the three gates that matter most (does the claim affect
# the dossier, is its support secondary-only, did we try to resolve it
# ourselves). Those are judgment and they stay with the agent; see SKILL.md
# "Asking back after intake". A green result here is permission to consider an
# ask, never permission to issue one.
# ---------------------------------------------------------------------------

ASK_TYPES = {"expand_node", "resolve_link", "test_gap", "new_question"}

# Which index an ask's `target` must resolve into. `new_question` points at
# nothing by design -- it is the ask for a claim the graph has no row for.
ASK_TARGET_INDEX = {
    "expand_node": "things",
    "resolve_link": "links",
    "test_gap": "gaps",
    "new_question": None,
}

# An ask has to be answerable by somebody who cannot see our run. "is obefazimod
# a TL1A agent?" is not; naming the sources on both sides is. We cannot judge
# prose, but we can insist the question carries identifiers -- two of them, so
# both sides of the disagreement are reachable.
SOURCE_TOKEN = re.compile(r"(PMC\d{5,}|NCT\d{8}|10\.\d{4,}/\S+|CHEMBL\d+|\b[0-9][A-Z0-9]{3}\b)")

MIN_SOURCE_TOKENS = 2

# An ask that carries OUR answer and contradicts a row the graph asserts. It is
# a correction, not a request for work, and the judgment gates were written for
# requests -- see SKILL.md "The one ask that skips gates 2 and 3". Detected from
# the `question` text because SKILL.md already requires such an ask to say so in
# prose ("Mark it plainly as post-resolution"), and because adding a field to the
# ask object would be a schema change upstream never agreed to.
POST_RESOLUTION = re.compile(r"post[\s_-]?resolution", re.I)


def is_secondary(finding, papers):
    """A finding that restates someone else's work rather than reporting one.

    Three independent signals, any of which is enough. `flags: [background]` is
    the extractor's own call, `is_own_result: false` is the schema's, and a
    review paper is the journal's. They disagree often enough that reading only
    one of them misses cases -- f1/f2 in the ask-back fixture carry all three,
    but a class table inside a primary trial report carries only the first.
    """
    if "background" in (finding.get("flags") or []):
        return True
    if finding.get("is_own_result") is False:
        return True
    study = (papers.get(finding.get("paper"), {}) or {}).get("study_type") or ""
    return study in {"review", "meta_analysis"}


def question_identity(question):
    """What makes two TARGETLESS asks the same ask.

    The identity is the SET OF SOURCE IDENTIFIERS the question names, not the
    question text and not a hash of it. Two reasons, and the second is the one
    that decides:

    - Text is not stable across the round trip. The upstream team rewords a
      question when it services it, and `rounds` is the only record of what was
      asked, so a hash of the wording calls a rephrasing a different ask and lets
      the same question through twice.
    - The identifiers are what make an ask routable at all -- QUESTION_IS_ACTIONABLE
      already refuses a question that names fewer than two of them -- and they
      survive rewording, translation and reordering.

    Falls back to normalised text only when a question names no identifier, which
    is a question that cannot pass QUESTION_IS_ACTIONABLE anyway.

    Inherits SOURCE_TOKEN's looseness deliberately, so that "source identifier"
    means one thing in this file. That regex's PDB-code branch also matches any
    bare four-digit number, so a quoted measurement or a year joins the identity.
    The error runs toward calling two asks DIFFERENT, i.e. toward letting an ask
    through -- the opposite direction from the bug this replaced, and the safer
    one, since a duplicate ask costs a round and a false match costs the verb.
    """
    tokens = tuple(sorted(set(SOURCE_TOKEN.findall(question or ""))))
    if tokens:
        return ("sources", tokens)
    return ("text", NS_WORD.sub(" ", (question or "").lower()).strip())


def already_asked(graph, ask_type, target, question=None):
    """Rounds already issued against this graph. Returns (hits, unmatchable).

    Re-asking is not harmless. It costs a round, returns the same evidence, and
    -- because `rounds` is the only record of what was tried -- makes the graph
    look like it was interrogated twice as hard as it was.

    Matched on (verb, target) for the three verbs that HAVE a target.
    `new_question` does not: SCHEMA.md gives it `target: null` by design, so
    (verb, target) is `("new_question", None)` for every new_question ever
    issued and the first one asked retires the verb for the life of the graph.
    Failure mode 18. A targetless ask is matched on `question_identity` instead.

    `unmatchable` carries the targetless rounds that record no `question` text,
    because those cannot be compared to anything: they are prior asks this gate
    is structurally unable to see, and reporting them as "no match" would be a
    false all-clear. The caller says so rather than deciding it.
    """
    hits, unmatchable = [], []
    targetless = target is None and ASK_TARGET_INDEX.get(ask_type, False) is None
    for r in read_list(graph, "rounds", require_id=False):
        if r.get("ask") != ask_type:
            continue
        if not targetless:
            if r.get("target") == target:
                hits.append(r)
            continue
        prior_question = r.get("question")
        if not prior_question:
            unmatchable.append(r)
        elif question_identity(prior_question) == question_identity(question):
            hits.append(r)
    return hits, unmatchable


def check_ask(graph, ask):
    """Mechanical gates on one proposed ask. Returns (gates, unchecked, exempt).

    `gates` is a list of {gate, ok, detail}. `unchecked` names the judgment
    gates this function deliberately does not evaluate, so a caller cannot read
    an all-green result as approval. `exempt` names the judgment gates that do
    not apply to THIS ask and says why -- empty for every ordinary ask.
    """
    idx = index(graph)
    gates = []

    def gate(name, ok, detail):
        gates.append({"gate": name, "ok": bool(ok), "detail": detail})

    ask_type = ask.get("ask")
    target = ask.get("target")
    question = (ask.get("question") or "").strip()

    gate(
        "ASK_TYPE_IS_ONE_OF_FOUR",
        ask_type in ASK_TYPES,
        f"{ask_type!r}; SCHEMA.md defines exactly {sorted(ASK_TYPES)}. "
        "Inventing a fifth verb makes the ask unconsumable upstream.",
    )

    if ask_type in ASK_TARGET_INDEX:
        want = ASK_TARGET_INDEX[ask_type]
        if want is None:
            gate(
                "TARGET_SHAPE",
                target is None,
                "new_question takes target: null -- it is the ask for a claim "
                f"the graph has no row for. Got {target!r}.",
            )
        else:
            gate(
                "TARGET_RESOLVES_TO_A_ROW",
                target in idx[want],
                f"{ask_type} must point at an id in `{want}`. "
                f"{target!r} {'resolves' if target in idx[want] else 'does NOT resolve'}. "
                "Pointing in prose instead of by id is how an ask becomes unroutable.",
            )
    else:
        gate("TARGET_RESOLVES_TO_A_ROW", False, "unknown ask type; target not checked")

    tokens = sorted(set(SOURCE_TOKEN.findall(question)))
    gate(
        "QUESTION_IS_ACTIONABLE",
        len(tokens) >= MIN_SOURCE_TOKENS,
        f"found {len(tokens)} source identifier(s) {tokens} in `question`; need "
        f"at least {MIN_SOURCE_TOKENS}. An ask must name the sources on BOTH "
        "sides and what would settle it, or the answering team reruns our search.",
    )

    prior, unmatchable = already_asked(graph, ask_type, target, question)
    matched_on = "target" if target is not None else "the source ids in `question`"
    gate(
        "NOT_ALREADY_ASKED",
        not prior,
        f"rounds already carries {len(prior)} ask(s) of {ask_type} at {target!r}, "
        f"matched on {matched_on}"
        + (f" (round {prior[0].get('n')}, outcome {prior[0].get('outcome')!r})" if prior else "")
        + (
            f". NOTE: {len(unmatchable)} prior {ask_type} round(s) "
            f"{[r.get('n') for r in unmatchable]} record no `question` text, so this "
            "gate cannot compare against them -- a targetless ask has no target to "
            "match on and upstream does not have to write the question into `rounds`. "
            "Read those rounds before issuing this."
            if unmatchable else ""
        ),
    )

    coverage = graph.get("coverage", {}) or {}
    stop = coverage.get("stop_reason")
    gate(
        "LITERATURE_NOT_EXHAUSTED",
        stop != "complete",
        f"coverage.stop_reason is {stop!r}. Only 'complete' means the literature "
        "was exhausted; against a complete graph another round returns the same "
        "evidence, so the ask is noise.",
    )

    # Gates 2 and 3 are judgment and are not evaluated here. But both are written
    # for an ask that REQUESTS work, and one ask type does not: the post-resolution
    # contradiction ask, which carries our own settled answer against a row the
    # graph asserts. Naming it as unchecked-and-therefore-required routed the one
    # ask type that has proven valuable nowhere -- gate 2's own rationale ("either
    # usable or contradicted by something we can measure") names the contradiction
    # case and then had no branch for it.
    post_resolution = bool(POST_RESOLUTION.search(question))
    unchecked = [
        "AFFECTS_THE_DOSSIER -- does this claim change a value in the template? "
        "An efficacy or clinical-outcome argument does not touch a tractability "
        "number (dossier rule 7) and fails this gate."
    ]
    exempt = []
    if post_resolution:
        exempt += [
            "SUPPORT_IS_SECONDARY_ONLY (gate 2) -- EXEMPT. `question` declares this "
            "ask post-resolution, so it carries our answer rather than a request for "
            "one. Gate 2 exists because a primary-supported claim is 'either usable "
            "or contradicted by something we can measure' -- this IS the second "
            "branch, and a wrong `primary` row is more damaging than a wrong "
            "`background_only` one, not less. The basis of the target link is "
            "irrelevant to an ask that is a correction.",
            "WE_TRIED_AND_FAILED (gate 3) -- EXEMPT, and for the same reason it "
            "always was: gate 3 stops us outsourcing work we did not do, and here "
            "the work is done. What replaces both gates is stricter, not looser -- "
            "the ask MUST state our answer, its source and its date, and must not "
            "block any dossier field. An ask claiming post-resolution and carrying "
            "only a doubt is the abuse this exemption creates; see failure mode 19.",
        ]
    else:
        unchecked.append(
            "SUPPORT_IS_SECONDARY_ONLY -- is every source a review, class table or "
            "citing paraphrase? A primary or mixed basis is not an ask, UNLESS the "
            "ask is a post-resolution contradiction carrying our own answer, which "
            "is exempt and must say so in `question`."
        )
        unchecked.append(
            "WE_TRIED_AND_FAILED -- were ChEMBL, the registry, the structure and a "
            "corpus grep on the exact identifiers all run and all silent, and is "
            "each null recorded in not_found BEFORE the ask was written?"
        )
    return gates, unchecked, exempt


def ask_context(graph):
    """Per-link evidence for the ask-back decision. Evidence, not a verdict.

    Same discipline as `signals()`: every field here is a deterministic fact
    about the graph, and none of them decides. In particular
    `support_is_secondary_only` is necessary and nowhere near sufficient -- most
    secondary-only links are resolvable from ChEMBL or the registry in one
    query, and those must never become asks.
    """
    idx = index(graph)
    coverage = graph.get("coverage", {}) or {}
    stop = coverage.get("stop_reason")

    rows = []
    for link in idx["links"].values():
        fids = link_findings(link)
        findings = [idx["findings"][f] for f in fids if f in idx["findings"]]
        secondary = [f["id"] for f in findings if is_secondary(f, idx["papers"])]
        prior, _ = already_asked(graph, "resolve_link", link["id"])
        secondary_only = bool(findings) and len(secondary) == len(findings)
        # Everything except gate 2. Split out because gate 2 is the one gate a
        # post-resolution contradiction ask is exempt from, and folding it into a
        # single boolean made L4 (basis primary, the ask that WORKED) and L2
        # (basis primary, the ask that must never fire) indistinguishable in this
        # output -- both just `mechanical_gates_clear: false`.
        clear_except_gate2 = bool(findings) and not prior and stop != "complete"
        rows.append({
            "link": link["id"],
            "relation": (
                f"{idx['things'].get(link.get('from'), {}).get('name')} "
                f"{link.get('how')} "
                f"{idx['things'].get(link.get('to'), {}).get('name')}"
            ),
            "basis": link.get("basis"),
            "n_findings": len(findings),
            "secondary_findings": secondary,
            "support_is_secondary_only": secondary_only,
            "distinct_papers": sorted({f.get("paper") for f in findings if f.get("paper")}),
            "already_asked_rounds": [r.get("n") for r in prior],
            # The basis test comes FIRST and is the one that decides, exactly as
            # in step 3: the tier comes from the link's `basis`, never from the
            # findings' own properties. The two can disagree -- in the RA graph
            # L6 (IL-6 drives RA) is `basis: primary` while its only finding is
            # a meta-analysis, so `support_is_secondary_only` is true and the
            # basis is not. `basis` wins, L6 is not an ask, and the divergence
            # is worth reading rather than smoothing over.
            "mechanical_gates_clear": (
                link.get("basis") in NON_ACTIONABLE_BASIS
                and secondary_only
                and clear_except_gate2
            ),
            # The post-resolution branch. True where every mechanical gate EXCEPT
            # gate 2 is clear -- i.e. this link is only blocked by carrying primary
            # or mixed support. That blocks a request for work, and it must not
            # block a correction: if we have settled this row ourselves and our
            # answer contradicts it, the ask is legitimate here and a `primary`
            # basis makes the wrong row MORE worth correcting, not less.
            # Still not permission: it requires that we actually did the work and
            # that the ask carries our answer and its source.
            "clear_if_post_resolution_contradiction": (
                clear_except_gate2 and link.get("basis") not in NON_ACTIONABLE_BASIS
            ),
        })

    return {
        "graph_id": graph.get("graph_id"),
        "round": graph.get("round"),
        "stop_reason": stop,
        "literature_not_exhausted": stop != "complete",
        "asks_already_issued": [
            {"n": r.get("n"), "ask": r.get("ask"), "target": r.get("target"),
             "outcome": r.get("outcome")}
            for r in read_list(graph, "rounds", require_id=False)
        ],
        "links": sorted(rows, key=lambda r: r["link"]),
        "_warning": (
            "mechanical_gates_clear is NOT permission to ask. It says only that "
            "the link is secondary-only, unasked, and against a non-exhausted "
            "graph. The three gates that matter -- does it affect the dossier, "
            "and did WE try ChEMBL / the registry / the structure first, and is "
            "each failed attempt recorded -- are not evaluated here and cannot "
            "be. See SKILL.md 'What must never become an ask'."
        ),
        "_post_resolution_note": (
            "clear_if_post_resolution_contradiction is the OTHER direction, and it "
            "is the one that has fired in practice: we settled a claim ourselves "
            "and our answer contradicts this row. Such an ask is a CORRECTION, not "
            "a request, so gates 2 and 3 do not apply -- it must instead state our "
            "answer, its source and its date, say 'post-resolution' in `question`, "
            "and block no dossier field. Everything else still applies: an id to "
            "point at, no prior round, a non-exhausted graph, and two source "
            "identifiers. A true value here on a link we have NOT actually measured "
            "against is the abuse this field creates."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("graph", help="path to the upstream evidence graph JSON")
    ap.add_argument("--thing", help="restrict nominations to one thing id")
    ap.add_argument(
        "--allow-fixture",
        action="store_true",
        help="permit a graph carrying _fixture: true (test runs only)",
    )
    ap.add_argument(
        "--ask-context",
        action="store_true",
        help="per-link mechanical facts bearing on whether to ask back",
    )
    ap.add_argument(
        "--check-ask",
        metavar="JSON",
        help="gate one proposed ask; exits 1 if any mechanical gate fails",
    )
    args = ap.parse_args()

    try:
        with open(args.graph) as fh:
            # parse_constant fires on NaN/Infinity, which json.load otherwise
            # accepts and json.dump re-emits as invalid JSON for a strict consumer.
            graph = json.load(fh, parse_constant=_reject_constant)
    except GraphShapeError as exc:
        sys.exit(f"refusing: {exc}")

    if not isinstance(graph, dict):
        sys.exit(f"refusing: top level is {type(graph).__name__}, not an object.")

    if graph.get("_fixture") and not args.allow_fixture:
        sys.exit(
            "refusing: graph carries _fixture: true, so its papers and quotes are "
            "synthetic. Re-run with --allow-fixture for a test."
        )

    try:
        if args.check_ask:
            gates, unchecked, exempt = check_ask(graph, json.loads(args.check_ask))
            json.dump({"gates": gates, "not_checked_here": unchecked,
                       "exempt_for_this_ask": exempt},
                      sys.stdout, indent=2)
            sys.stdout.write("\n")
            sys.exit(0 if all(g["ok"] for g in gates) else 1)

        if args.ask_context:
            payload = ask_context(graph)
        else:
            payload = build(graph, args.thing)
    except GraphShapeError as exc:
        # Loud and specific. The alternative -- a bare TypeError, or worse a
        # clean exit with zero nominations -- reads as "no targets in this
        # literature", which is the one wrong answer this file must never give.
        sys.exit(f"refusing: {exc}")

    # allow_nan=False so a NaN that got this far cannot leave as invalid JSON.
    json.dump(payload, sys.stdout, indent=2, allow_nan=False)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
