"""Build the per-run interpretability object carried at ``output.interpretability``.

This is a PURE function of a druggability dossier: every value it emits is pulled
from a dossier field, nothing is fabricated, and it is robust to nulls and missing
keys. Its output validates against ``schemas/interpretability.schema.json`` for any
schema-valid dossier.

The mapping is ported from the display heuristic in ``INTERPRETABILITY_PANEL.md``
(the spec that governs the ``interpretability`` block shipped in
``examples/output.json``): dossier fields -> headline / figure / two axes / trace /
not_found / provenance, following the section-by-section rules there. The object is
a WORKFLOW reasoning trace (which stages ran, their intermediate scores, provenance
and caveats) -- never LLM chain-of-thought.

No numbers are invented here. A field that is null in the dossier is rendered as
null (or omitted), never as a zero.
"""

from __future__ import annotations

from typing import Any

SCHEMA_VERSION = "1.0"

_VERDICT_LABELS = {
    "small_molecule_tractable": "Small-molecule tractable",
    "not_tractable": "Not tractable",
    "insufficient_evidence": "Insufficient evidence",
}


def _get(obj: Any, *path: str, default: Any = None) -> Any:
    """Safe nested lookup: returns ``default`` on any missing key or non-dict."""
    cur = obj
    for key in path:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(key, default)
    return cur


def _fmt_num(value: Any) -> str | None:
    """Format a number for display without inventing precision. Non-numbers -> None."""
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, (int, float)):
        # %g drops trailing zeros while preserving the value's own precision.
        return f"{value:g}"
    return None


def _headline(dossier: dict) -> dict:
    verdict = dossier.get("verdict")
    verdict_basis = dossier.get("verdict_basis")

    tp = dossier.get("target_precedent") or {}
    distinct_actives = tp.get("distinct_actives")
    best_potency = tp.get("best_potency_nm")
    best_modality = tp.get("best_potency_modality")
    approved_count = tp.get("approved_small_molecules_count")
    clinical = tp.get("clinical_stage_small_molecules") or []

    precedent_bits: list[str] = []
    if isinstance(distinct_actives, int):
        precedent_bits.append(f"{distinct_actives} small-molecule actives")
    bp = _fmt_num(best_potency)
    if bp is not None:
        suffix = ""
        if best_modality and best_modality != "small_molecule":
            suffix = f" ({best_modality})"
        precedent_bits.append(f"best {bp} nM{suffix}")
    if isinstance(approved_count, int) and approved_count > 0:
        precedent_bits.append(f"{approved_count} approved small molecules")
    if clinical:
        precedent_bits.append(f"{len(clinical)} clinical-stage small molecules")

    site_phrase = None
    pdb = _get(dossier, "tractability", "site_pocket_rank", "structure_pdb_id")
    if _get(dossier, "tractability", "mdpocket_site_definition_used") == "site_from_ligand" and pdb:
        site_phrase = f"ligand-anchored pocket at {pdb}"

    segments: list[str] = []
    if precedent_bits:
        segments.append(", ".join(precedent_bits))
    if site_phrase:
        segments.append(site_phrase)

    if segments:
        one_line = "; ".join(segments) + "."
    else:
        label = _VERDICT_LABELS.get(verdict, verdict or "Verdict")
        one_line = (
            f"{label} — see the two axes and the 'Not retrieved' list for what "
            "was and was not established."
        )

    headline: dict[str, Any] = {
        "verdict": verdict,
        "verdict_label": _VERDICT_LABELS.get(verdict),
        "one_line": one_line,
    }
    if verdict_basis in ("retrieved_precedent", "computed_tractability", "both", "none"):
        headline["basis"] = verdict_basis
    return headline


def _figure(dossier: dict) -> dict:
    pdb = _get(dossier, "tractability", "site_pocket_rank", "structure_pdb_id")
    mdpdef = _get(dossier, "tractability", "mdpocket_site_definition_used")
    if mdpdef == "site_from_ligand" and pdb:
        rank = _get(dossier, "tractability", "site_pocket_rank", "fpocket")
        n = _get(dossier, "tractability", "site_pocket_rank", "n_pockets")
        caption = f"{pdb} pocket"
        if isinstance(rank, int) and isinstance(n, int):
            caption = f"{pdb} pocket, rank {rank} of {n}"
        return {
            "kind": "pocket_render",
            "path": f"figures/{pdb}_pocket.png",
            "alt": f"{pdb} ligand-anchored pocket",
            "caption": caption,
            "placement": "right_of_headline",
        }
    # A precedent-carried or unestablished-site verdict often has no figure.
    return {"kind": "none", "path": None, "alt": None, "caption": None, "placement": "none"}


def _banner(dossier: dict) -> dict | None:
    axis_conflict = dossier.get("axis_conflict")
    tier = _get(dossier, "structure", "tier")
    messages: list[str] = []
    level = "info"
    if isinstance(axis_conflict, str) and axis_conflict.strip():
        messages.append(axis_conflict)
        level = "warning"
    if tier == "predicted":
        plddt = _get(dossier, "structure", "predicted_plddt")
        note = "Structure is a predicted model, not experimental"
        if plddt is not None:
            note += f" (pLDDT {_fmt_num(plddt)})"
        note += "; treat the computed axis as flagged, not measured."
        messages.append(note)
        level = "warning"
    if not messages:
        return None
    return {"level": level, "messages": messages}


def _retrieved_axis(dossier: dict) -> dict:
    tp = dossier.get("target_precedent") or {}
    verdict_basis = dossier.get("verdict_basis")

    distinct_actives = tp.get("distinct_actives")
    best_potency = tp.get("best_potency_nm")
    best_modality = tp.get("best_potency_modality")
    best_characterised = tp.get("best_potency_characterised")
    best_assay = tp.get("best_potency_assay")
    approved_count = tp.get("approved_small_molecules_count")
    clinical = tp.get("clinical_stage_small_molecules") or []

    has_precedent = bool(
        (isinstance(distinct_actives, int) and distinct_actives > 0)
        or best_potency is not None
        or (isinstance(approved_count, int) and approved_count > 0)
        or clinical
    )
    if verdict_basis in ("retrieved_precedent", "both") or has_precedent:
        status = "supports"
    else:
        status = "insufficient"

    bp = _fmt_num(best_potency)
    potency_caveat = None
    if best_characterised is False:
        potency_caveat = "uncharacterised assay — unusable for a potency claim (rule 6)"
    elif best_modality and best_modality != "small_molecule":
        potency_caveat = f"modality is {best_modality}, not a small-molecule claim (rule 1d)"

    headline_number = None
    if bp is not None:
        headline_number = {"label": "best potency", "value": f"{bp} nM", "caveat": potency_caveat}
    elif isinstance(approved_count, int):
        headline_number = {
            "label": "approved small molecules",
            "value": str(approved_count),
            "caveat": None,
        }

    rows: list[dict[str, Any]] = []
    if bp is not None:
        rows.append(
            {
                "label": "best potency",
                "value": f"{bp} nM",
                "caveat": potency_caveat,
                "provenance": best_assay if isinstance(best_assay, str) else None,
            }
        )
    if isinstance(distinct_actives, int):
        rows.append(
            {
                "label": "distinct small-molecule actives",
                "value": distinct_actives,
                "caveat": "small molecules only; a raw count measures assays, not the target (rules 1b, 6, 14)",
                "provenance": None,
            }
        )
    if isinstance(approved_count, int):
        rows.append(
            {
                "label": "approved small molecules",
                "value": approved_count,
                "caveat": None,
                "provenance": None,
            }
        )
    rows.append(
        {
            "label": "clinical-stage small molecules",
            "value": len(clinical),
            "caveat": None,
            "provenance": None,
        }
    )

    axis: dict[str, Any] = {
        "id": "retrieved_precedent",
        "title": "Retrieved precedent",
        "status": status,
        "rows": rows,
        "display": "column",
    }
    if headline_number is not None:
        axis["headline_number"] = headline_number
    return axis


def _computed_axis(dossier: dict) -> dict:
    tract = dossier.get("tractability") or {}
    site_basis = tract.get("site_hypothesis_basis")
    rank = tract.get("site_pocket_rank") or {}
    fpocket = rank.get("fpocket")
    prank = rank.get("prank")
    n_pockets = rank.get("n_pockets")
    pdb = rank.get("structure_pdb_id")
    vol = _get(tract, "pocket_volume_a3", "primary_d1_6_a3")
    jaccard = tract.get("ligand_site_jaccard")
    cryptic = tract.get("cryptic_pocket_risk")

    site_not_established = isinstance(site_basis, str) and "not_established" in site_basis.lower()
    has_rank = isinstance(fpocket, int) and isinstance(n_pockets, int) and bool(pdb)

    if not tract or site_not_established or not has_rank:
        # Present-but-not-reportable, or no site at all: the computed axis carries no verdict.
        status = "not_run"
        rows = [
            {
                "label": "site hypothesis",
                "value": "not established" if site_not_established else "no reportable site",
                "caveat": (
                    "no ligand-anchored site established; the computed axis carries no "
                    "verdict (rules 4b, 4.0)"
                ),
                "provenance": None,
            }
        ]
        caveat = tract.get("caveat")
        if isinstance(caveat, str) and caveat.strip():
            rows.append({"label": "caveat", "value": caveat, "caveat": None, "provenance": None})
        return {
            "id": "computed_tractability",
            "title": "Computed tractability",
            "status": status,
            "rows": rows,
            "display": "column",
        }

    # A reportable within-structure rank exists. Druggability is never load-bearing,
    # so the computed axis is 'insufficient' (present, informative, but not verdict-bearing).
    rank_value = f"rank {fpocket} of {n_pockets} in {pdb}"
    rows = [
        {
            "label": "druggability rank",
            "value": rank_value,
            "caveat": "within-structure rank, never compared across structures (rule 4)",
            "provenance": pdb,
        }
    ]
    if isinstance(prank, int):
        rows.append(
            {
                "label": "PRANK rank",
                "value": f"rank {prank}",
                "caveat": "site-finding aid, not a quality score (rule 4d)",
                "provenance": None,
            }
        )
    vol_fmt = _fmt_num(vol)
    if vol_fmt is not None:
        rows.append(
            {
                "label": "pocket volume (D=1.6)",
                "value": f"{vol_fmt} A^3",
                "caveat": "absolute quantity; the 210/240 A^3 guide is withdrawn (rule 4a)",
                "provenance": None,
            }
        )
    if isinstance(jaccard, (int, float)) and not isinstance(jaccard, bool):
        rows.append(
            {"label": "ligand-site Jaccard", "value": jaccard, "caveat": None, "provenance": None}
        )
    if cryptic is not None:
        rows.append(
            {"label": "cryptic risk", "value": cryptic, "caveat": None, "provenance": None}
        )

    return {
        "id": "computed_tractability",
        "title": "Computed tractability",
        "status": "insufficient",
        "headline_number": {
            "label": "druggability rank",
            "value": rank_value,
            "caveat": "within-structure only",
        },
        "rows": rows,
        "display": "column",
    }


def _trace(dossier: dict) -> list[dict]:
    verdict = dossier.get("verdict")
    verdict_basis = dossier.get("verdict_basis")

    resolved = _get(dossier, "target", "uniprot_accession") or _get(
        dossier, "input", "uniprot_accession"
    )
    gene = _get(dossier, "target", "gene_symbol")
    resolved_summary = f"Resolved {resolved}" if resolved else "Intake"
    if gene:
        resolved_summary += f" ({gene})"

    tier = _get(dossier, "structure", "tier")
    pdb = _get(dossier, "structure", "pdb_id")
    holo = _get(dossier, "structure", "holo_count")
    apo = _get(dossier, "structure", "apo_count")
    struct_summary = f"{tier} — {pdb}" if tier and pdb else (tier or "no structure")

    rank = _get(dossier, "tractability", "site_pocket_rank") or {}
    fpocket = rank.get("fpocket")
    n_pockets = rank.get("n_pockets")
    rank_pdb = rank.get("structure_pdb_id")
    if isinstance(fpocket, int) and isinstance(n_pockets, int) and rank_pdb:
        scan_summary = f"rank {fpocket} of {n_pockets} in {rank_pdb}"
    else:
        scan_summary = "no reportable site"
    vol = _fmt_num(_get(dossier, "tractability", "pocket_volume_a3", "primary_d1_6_a3"))
    cryptic = _get(dossier, "tractability", "cryptic_pocket_risk")

    tp = dossier.get("target_precedent") or {}
    actives = tp.get("distinct_actives")
    bp = _fmt_num(tp.get("best_potency_nm"))
    clinical = tp.get("clinical_stage_small_molecules") or []
    prec_summary = f"{actives if actives is not None else 'None'} actives"
    if bp is not None:
        prec_summary += f", best {bp} nM"

    fals = dossier.get("falsification") or {}
    checks = fals.get("checks_run") or []
    survived = fals.get("survived")
    fals_summary = "survived" if survived else ("did not survive" if survived is False else "run")

    open_stages: set[str] = set()
    if verdict_basis == "retrieved_precedent":
        open_stages = {"precedent-lookup"}
    elif verdict_basis == "computed_tractability":
        open_stages = {"pocket-scan"}
    elif verdict_basis == "both":
        open_stages = {"pocket-scan", "precedent-lookup"}
    else:
        open_stages = {"assemble-dossier"}

    def disp(stage: str) -> str:
        return "open" if stage in open_stages else "collapsed"

    stages: list[dict[str, Any]] = [
        {
            "stage": "graph-intake",
            "order": 1,
            "summary": resolved_summary,
            "status": "ok",
            "display": disp("graph-intake"),
            "fields": [{"label": "accession", "value": resolved, "caveat": None}],
        },
        {
            "stage": "structure-select",
            "order": 2,
            "summary": struct_summary,
            "status": "ok" if tier and tier != "none" else "not_run",
            "display": disp("structure-select"),
            "fields": [
                {"label": "tier", "value": tier, "caveat": None},
                {
                    "label": "holo / apo",
                    "value": f"{holo} / {apo}" if holo is not None and apo is not None else None,
                    "caveat": None,
                },
            ],
        },
        {
            "stage": "pocket-scan",
            "order": 3,
            "summary": scan_summary,
            "status": "ok",
            "display": disp("pocket-scan"),
            "fields": [
                {"label": "druggability rank", "value": scan_summary, "caveat": "within-structure"},
                {"label": "volume D=1.6", "value": f"{vol} A^3" if vol else None, "caveat": None},
                {"label": "cryptic", "value": cryptic, "caveat": None},
            ],
        },
        {
            "stage": "precedent-lookup",
            "order": 4,
            "summary": prec_summary,
            "status": "ok",
            "display": disp("precedent-lookup"),
            "fields": [
                {"label": "distinct actives", "value": actives, "caveat": None},
                {"label": "best potency", "value": f"{bp} nM" if bp else None, "caveat": None},
                {"label": "clinical-stage SMs", "value": len(clinical), "caveat": None},
            ],
        },
        {
            "stage": "falsification",
            "order": 5,
            "summary": fals_summary,
            "status": "ok",
            "display": disp("falsification"),
            "fields": [{"label": "checks run", "value": len(checks), "caveat": None}],
        },
        {
            "stage": "assemble-dossier",
            "order": 6,
            "summary": f"{verdict} / {verdict_basis}",
            "status": "ok",
            "display": disp("assemble-dossier"),
            "fields": [
                {"label": "verdict", "value": verdict, "caveat": None},
                {"label": "basis", "value": verdict_basis, "caveat": None},
            ],
        },
    ]

    # Attach the pocket render to the pocket-scan stage when there is one.
    fig = _figure(dossier)
    if fig.get("kind") == "pocket_render":
        stages[2]["figure"] = {
            "kind": "pocket_render",
            "path": fig["path"],
            "caption": fig["caption"],
        }
    return stages


def _not_found(dossier: dict) -> list[dict]:
    out: list[dict[str, Any]] = []
    for item in dossier.get("not_found") or []:
        if isinstance(item, str):
            out.append({"field": None, "reason": item})
        elif isinstance(item, dict):
            out.append(
                {
                    "field": item.get("field"),
                    "reason": item.get("reason") or item.get("signature"),
                }
            )
    return out


def _provenance(dossier: dict) -> dict:
    sources: list[str] = []
    seen: set[str] = set()
    for block in (
        "target",
        "target_precedent",
        "family_precedent",
        "structural_neighbour_precedent",
        "pocket_neighbour_precedent",
        "structure",
        "tractability",
        "affinity",
    ):
        for src in _get(dossier, block, "sources", default=[]) or []:
            if isinstance(src, str) and src and src not in seen:
                seen.add(src)
                sources.append(src)
    return {"as_of_date": dossier.get("as_of_date"), "sources": sources}


def build_interpretability(dossier: dict) -> dict:
    """Map a druggability dossier to its interpretability side-panel object.

    Pure function. Reads only fields present on ``dossier``; robust to nulls and
    missing keys. The result validates against
    ``schemas/interpretability.schema.json`` for any schema-valid dossier.
    """
    if not isinstance(dossier, dict):
        raise TypeError("build_interpretability expects a dossier dict")

    return {
        "module": "simulation",
        "schema_version": SCHEMA_VERSION,
        "headline": _headline(dossier),
        "figure": _figure(dossier),
        "banner": _banner(dossier),
        "axes": [_retrieved_axis(dossier), _computed_axis(dossier)],
        "trace": _trace(dossier),
        "not_found": _not_found(dossier),
        "provenance": _provenance(dossier),
    }
