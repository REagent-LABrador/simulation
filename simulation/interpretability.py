"""Build the per-run ``output.interpretability`` object — the LABrador shared
interpretability contract (message.txt, section E: Tractability/Simulation).

This is a PURE, DETERMINISTIC function of a druggability dossier. It MAPS existing
authoritative dossier fields into the common contract; it never recomputes a
scientific value and never fabricates evidence, formulas, provenance, confidence
or precision. A field that is null/absent in the dossier stays null and earns a
limitation — it is never silently turned into 0, false, or an empty finding.

Common-contract shape (schemas/interpretability.schema.json is authoritative):
    schema_version, headline{title,result,plain_language,status,basis[]},
    metrics[], steps[], evidence[], assumptions[],
    uncertainty{method,intervals[],seed,draws,limitations[]},
    limitations[], counterfactuals[], lineage[], extensions{}.

Section-E specifics that live here:
  * the two evidence axes are preserved SEPARATELY under extensions.axes — never
    collapsed into one scalar;
  * the ordered trace, next_experiment, structured PDB/ChEMBL/tool identifiers and
    the figure reference live under extensions;
  * every not_found item becomes a limitation with a machine-readable field_path;
  * "actives not retrieved" is stated explicitly, never as "None actives" or 0;
  * extensions.cache_key is a hash of the COMPLETE input, not only the accession.

The previous axes/trace/figure-shaped object is SUPERSEDED by this contract; its
content is retained under extensions so no information is lost.
"""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any

SCHEMA_VERSION = "1.0.0"

_VERDICT_LABELS = {
    "small_molecule_tractable": "Small-molecule tractable",
    "not_tractable": "Not tractable",
    "insufficient_evidence": "Insufficient evidence",
}

# Honest verdict -> headline.status mapping.
#   small_molecule_tractable -> SUPPORTED   (a positive claim the axes support)
#   insufficient_evidence    -> INCONCLUSIVE
#   not_tractable            -> QUALIFIED    (a negative conclusion reported
#       conservatively: this station's design forbids druggability alone from
#       carrying a negative (rule 4.0/4a), so a not_tractable verdict is
#       evidence-qualified rather than asserted as strongly SUPPORTED).
_STATUS = {
    "small_molecule_tractable": "SUPPORTED",
    "not_tractable": "QUALIFIED",
    "insufficient_evidence": "INCONCLUSIVE",
}


def _get(obj: Any, *path: str, default: Any = None) -> Any:
    cur = obj
    for key in path:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(key, default)
    return cur


def _isnum(v: Any) -> bool:
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def _fmt(v: Any) -> str | None:
    if not _isnum(v):
        return None
    return f"{v:g}"


# --------------------------------------------------------------------------- #
# Evidence and assumptions registries — built first so metrics/steps can only
# reference ids that resolve. IDs are stable/semantic, never array-position.
# --------------------------------------------------------------------------- #
def _evidence(dossier: dict) -> list[dict]:
    tp = dossier.get("target_precedent") or {}
    out: list[dict] = []

    chembl = tp.get("chembl_target_id")
    if isinstance(chembl, str) and chembl:
        out.append({
            "id": "evidence.chembl_target",
            "claim": "Measured bioactivity records for the target in ChEMBL.",
            "source_type": "database",
            "source_id": f"ChEMBL:{chembl}",
            "source_url": f"https://www.ebi.ac.uk/chembl/target_report_card/{chembl}/",
            "locator": None,
            "quote": None,          # fail closed: no source text verified
            "grade": "HIGH",        # measured assay data
            "synthetic": False,
        })

    assay = tp.get("best_potency_assay")
    if isinstance(assay, str) and assay:
        characterised = tp.get("best_potency_characterised")
        out.append({
            "id": "evidence.best_potency_assay",
            "claim": "Assay underlying the best reported potency.",
            "source_type": "assay",
            "source_id": assay,
            "source_url": None,
            "locator": None,
            "quote": None,
            "grade": "HIGH" if characterised else "LOW",
            "synthetic": False,
        })

    # Structures actually scored / selected.
    seen_pdb: set[str] = set()
    for pdb in (
        _get(dossier, "tractability", "site_pocket_rank", "structure_pdb_id"),
        _get(dossier, "structure", "pdb_id"),
    ):
        if isinstance(pdb, str) and pdb and pdb not in seen_pdb:
            seen_pdb.add(pdb)
            out.append({
                "id": f"evidence.pdb.{pdb}",
                "claim": f"Experimental structure {pdb} used for the computed axis.",
                "source_type": "structure",
                "source_id": f"PDB:{pdb}",
                "source_url": f"https://www.rcsb.org/structure/{pdb}",
                "locator": None,
                "quote": None,
                "grade": "HIGH",
                "synthetic": False,
            })

    if _get(dossier, "structure", "tier") == "predicted":
        out.append({
            "id": "evidence.predicted_model",
            "claim": "ESMFold predicted model used because no usable experimental structure exists.",
            "source_type": "prediction",
            "source_id": "ESMFold",
            "source_url": None,
            "locator": None,
            "quote": None,
            "grade": "LOW",         # a model, not an observation
            "synthetic": False,     # a real model run, not fabricated data
        })
    return out


def _assumptions(dossier: dict) -> list[dict]:
    inp = dossier.get("input") or {}
    out: list[dict] = []
    mech = inp.get("mechanism_hypothesis")
    out.append({
        "id": "assumption.mechanism_hypothesis",
        "path": "input.mechanism_hypothesis",
        "value": mech,
        "unit": "category",
        "basis": "Caller-supplied structural question; not inferred by this module.",
        "synthetic": False,
    })
    out.append({
        "id": "assumption.as_of_date",
        "path": "input.as_of_date",
        "value": inp.get("as_of_date"),
        "unit": "date",
        "basis": "Retrieval cutoff; null means as-of-today with no cutoff enforced.",
        "synthetic": False,
    })
    site_basis = _get(dossier, "tractability", "site_hypothesis_basis")
    if site_basis is not None:
        out.append({
            "id": "assumption.site_basis",
            "path": "tractability.site_hypothesis_basis",
            "value": site_basis,
            "unit": None,
            "basis": "How the scored pocket's site was defined (ligand-anchored, signature, or not established).",
            "synthetic": False,
        })
    return out


# --------------------------------------------------------------------------- #
# Metrics — every numeric metric carries a unit; each references evidence and/or
# assumptions that resolve, else the caller adds an UNTAGGED_VALUE limitation.
# --------------------------------------------------------------------------- #
def _metrics(dossier: dict, evidence_ids: set[str]) -> list[dict]:
    tp = dossier.get("target_precedent") or {}
    tract = dossier.get("tractability") or {}
    out: list[dict] = []

    def ev(*ids: str) -> list[str]:
        return [i for i in ids if i in evidence_ids]

    bp = tp.get("best_potency_nm")
    if _isnum(bp):
        mod = tp.get("best_potency_modality")
        characterised = tp.get("best_potency_characterised")
        meaning = "Best measured potency of a small molecule against the target."
        if characterised is False:
            meaning = "Best raw potency; assay uncharacterised (rule 6) — not a usable potency claim."
        elif mod and mod != "small_molecule":
            meaning = f"Best measured potency, but modality is {mod}, not a small molecule (rule 1d)."
        out.append({
            "id": "metric.best_potency_nm",
            "label": "Best potency",
            "value": bp,
            "unit": "nM",
            "display": f"{_fmt(bp)} nM",
            "meaning": meaning,
            "direction": "positive" if bp is not None else "unknown",
            "evidence_ids": ev("evidence.best_potency_assay", "evidence.chembl_target"),
            "assumption_ids": [],
        })

    da = tp.get("distinct_actives")
    if isinstance(da, int):
        out.append({
            "id": "metric.distinct_actives",
            "label": "Distinct small-molecule actives",
            "value": da,
            "unit": "compounds",
            "display": f"{da} compounds",
            "meaning": "Distinct small molecules with measured activity; a count of compounds, not of assays (rules 1b, 14).",
            "direction": "positive",
            "evidence_ids": ev("evidence.chembl_target"),
            "assumption_ids": [],
        })

    ac = tp.get("approved_small_molecules_count")
    if isinstance(ac, int):
        out.append({
            "id": "metric.approved_small_molecules",
            "label": "Approved small molecules",
            "value": ac,
            "unit": "molecules",
            "display": f"{ac} molecules",
            "meaning": "Approved small-molecule drugs acting directly on the target.",
            "direction": "positive" if ac > 0 else "neutral",
            "evidence_ids": ev("evidence.chembl_target"),
            "assumption_ids": [],
        })

    clinical = tp.get("clinical_stage_small_molecules")
    if isinstance(clinical, list):
        out.append({
            "id": "metric.clinical_stage_small_molecules",
            "label": "Clinical-stage small molecules",
            "value": len(clinical),
            "unit": "molecules",
            "display": f"{len(clinical)} molecules",
            "meaning": "Small molecules that reached clinical development against the target.",
            "direction": "positive" if clinical else "neutral",
            "evidence_ids": ev("evidence.chembl_target"),
            "assumption_ids": [],
        })

    rank = tract.get("site_pocket_rank") or {}
    fp, n, rpdb = rank.get("fpocket"), rank.get("n_pockets"), rank.get("structure_pdb_id")
    if isinstance(fp, int) and isinstance(n, int) and rpdb:
        out.append({
            "id": "metric.druggability_rank",
            "label": "Druggability rank",
            "value": fp,
            "unit": "rank",
            "display": f"rank {fp} of {n} in {rpdb}",
            "meaning": "fpocket rank of the site pocket AMONG that structure's pockets — a within-structure ordering, never compared across structures (rule 4). Not load-bearing.",
            "direction": "neutral",
            "evidence_ids": ev(f"evidence.pdb.{rpdb}"),
            "assumption_ids": [i for i in ("assumption.site_basis",) ],
        })

    vol = _get(tract, "pocket_volume_a3", "primary_d1_6_a3")
    if _isnum(vol):
        pdb = rank.get("structure_pdb_id")
        out.append({
            "id": "metric.pocket_volume_d16",
            "label": "Pocket volume (D=1.6)",
            "value": vol,
            "unit": "A^3",
            "display": f"{_fmt(vol)} A^3",
            "meaning": "Absolute site volume at clustering D=1.6. Comparable across structures, but its clustering sensitivity travels with it; the 210/240 A^3 tractability guide is withdrawn (rule 4a).",
            "direction": "neutral",
            "evidence_ids": ev(f"evidence.pdb.{pdb}") if isinstance(pdb, str) else [],
            "assumption_ids": [],
        })

    jac = tract.get("ligand_site_jaccard")
    if _isnum(jac):
        pdb = rank.get("structure_pdb_id")
        out.append({
            "id": "metric.ligand_site_jaccard",
            "label": "Ligand-site Jaccard",
            "value": jac,
            "unit": "fraction",
            "display": _fmt(jac),
            "meaning": "Overlap between the scored pocket and the co-crystal ligand site.",
            "direction": "positive",
            "evidence_ids": ev(f"evidence.pdb.{pdb}") if isinstance(pdb, str) else [],
            "assumption_ids": [],
        })

    dis = tract.get("disorder_fraction")
    if _isnum(dis):
        out.append({
            "id": "metric.disorder_fraction",
            "label": "Disorder fraction",
            "value": dis,
            "unit": "fraction",
            "display": _fmt(dis),
            "meaning": "Fraction of the sequence predicted disordered; high disorder disfavours a stable small-molecule pocket.",
            "direction": "negative",
            "evidence_ids": [],
            "assumption_ids": [],
        })

    plddt = _get(dossier, "structure", "predicted_plddt")
    if _isnum(plddt):
        out.append({
            "id": "metric.predicted_plddt",
            "label": "Predicted-model pLDDT",
            "value": plddt,
            "unit": "score",
            "display": _fmt(plddt),
            "meaning": "ESMFold per-residue confidence (0-100). A model-confidence SCORE, not a probability of tractability.",
            "direction": "neutral",
            "evidence_ids": ev("evidence.predicted_model"),
            "assumption_ids": [],
        })
    return out


def _steps(dossier: dict, evidence_ids: set[str]) -> list[dict]:
    tract = dossier.get("tractability") or {}
    rank = tract.get("site_pocket_rank") or {}
    out: list[dict] = []

    def ev(*ids: str) -> list[str]:
        return [i for i in ids if i in evidence_ids]

    # The verdict decision (a rule, not arithmetic -> formula null).
    out.append({
        "id": "step.verdict_decision",
        "label": "Two-axis verdict",
        "method": "Independent axes (retrieved precedent, computed tractability) combined by rule; never averaged into a scalar.",
        "formula": None,
        "inputs": [
            {"path": "verdict_basis", "value": dossier.get("verdict_basis"), "unit": None},
            {"path": "target_precedent.best_potency_nm", "value": _get(dossier, "target_precedent", "best_potency_nm"), "unit": "nM"},
            {"path": "tractability.site_pocket_rank.fpocket", "value": rank.get("fpocket"), "unit": "rank"},
        ],
        "result": {"value": dossier.get("verdict"), "unit": None},
        "evidence_ids": [],
        "assumption_ids": ["assumption.mechanism_hypothesis"],
    })

    fp, n, rpdb = rank.get("fpocket"), rank.get("n_pockets"), rank.get("structure_pdb_id")
    if isinstance(fp, int) and isinstance(n, int) and rpdb:
        out.append({
            "id": "step.druggability_rank",
            "label": "Within-structure druggability rank",
            "method": "fpocket pocket scoring at clustering D=1.6; report the site pocket's rank among the structure's pockets.",
            "formula": "rank of site pocket among n_pockets by fpocket score",
            "inputs": [
                {"path": "tractability.site_pocket_rank.n_pockets", "value": n, "unit": "pockets"},
                {"path": "tractability.site_hypothesis_basis", "value": tract.get("site_hypothesis_basis"), "unit": None},
            ],
            "result": {"value": fp, "unit": "rank"},
            "evidence_ids": ev(f"evidence.pdb.{rpdb}"),
            "assumption_ids": ["assumption.site_basis"] if _get(dossier, "tractability", "site_hypothesis_basis") is not None else [],
        })
    return out


def _uncertainty(dossier: dict, metric_ids: set[str]) -> dict:
    tract = dossier.get("tractability") or {}
    vol = tract.get("pocket_volume_a3") or {}
    intervals: list[dict] = []
    limitations: list[str] = []

    lo, ce, hi = vol.get("min"), vol.get("primary_d1_6_a3"), vol.get("max")
    if _isnum(lo) and _isnum(hi) and "metric.pocket_volume_d16" in metric_ids:
        intervals.append({
            "metric_id": "metric.pocket_volume_d16",
            "low": lo,
            "central": ce if _isnum(ce) else None,
            "high": hi,
            "unit": "A^3",
            "confidence_level": None,   # not a CI
        })
        limitations.append(
            "The pocket-volume interval is a clustering-parameter SCENARIO range "
            "(fpocket D sweep), not a confidence interval or a percentile."
        )

    cons = tract.get("ensemble_consensus_fraction") or {}
    frac = cons.get("fraction_with_strong_pocket")
    if _isnum(frac):
        limitations.append(
            f"Ensemble consensus: fraction_with_strong_pocket={_fmt(frac)} is a HEURISTIC "
            "anti-cherry-picking control, not a probability of tractability."
        )

    if not intervals:
        limitations.append(
            "No probabilistic uncertainty model was run; ranges, where present, are "
            "clustering scenarios or heuristic scores, never confidence intervals."
        )

    return {
        "method": "fpocket clustering-scenario sweep and ensemble-consensus heuristic; no Monte Carlo / probabilistic model.",
        "intervals": intervals,
        "seed": None,
        "draws": None,
        "limitations": limitations,
    }


def _limitations(dossier: dict, metrics: list[dict]) -> list[dict]:
    out: list[dict] = []

    # not_found -> a limitation with a machine-readable field_path.
    for item in dossier.get("not_found") or []:
        if isinstance(item, str):
            field, reason = None, item
            if ":" in item:
                head, _, tail = item.partition(":")
                if " " not in head.strip():
                    field, reason = head.strip(), tail.strip()
            out.append({
                "code": "NOT_RETRIEVED",
                "severity": "WARNING",
                "message": reason,
                "field_path": f"output.{field}" if field else None,
            })
        elif isinstance(item, dict):
            field = item.get("field")
            msg = item.get("reason") or item.get("signature") or "not retrieved"
            sev = "ERROR" if isinstance(msg, str) and "error" in msg.lower() else "WARNING"
            out.append({
                "code": item.get("code") or "NOT_RETRIEVED",
                "severity": sev,
                "message": msg,
                "field_path": f"output.{field}" if isinstance(field, str) and field else None,
            })

    # Actives not retrieved — stated explicitly, never "None actives"/0.
    if _get(dossier, "target_precedent", "distinct_actives") is None:
        out.append({
            "code": "ACTIVES_NOT_RETRIEVED",
            "severity": "WARNING",
            "message": "Distinct small-molecule actives were not retrieved (null, not zero); the count is unknown, not absent.",
            "field_path": "output.target_precedent.distinct_actives",
        })

    # Drug-bound-pocket false negative: a low within-structure druggability on a
    # holo structure is the documented false negative, not evidence against a site.
    drug = _get(dossier, "tractability", "pocket_druggability") or {}
    dmax = drug.get("max")
    if _isnum(dmax) and dmax < 0.1 and _get(dossier, "structure", "tier") == "holo_experimental":
        out.append({
            "code": "DRUGGABILITY_DRUG_BOUND_FALSE_NEGATIVE",
            "severity": "INFO",
            "message": "Low fpocket druggability on a ligand-bound pocket is the documented drug-bound false negative (rule 4.0); it is not evidence against the site.",
            "field_path": "output.tractability.pocket_druggability",
        })

    # Predicted-structure flag.
    if _get(dossier, "structure", "tier") == "predicted":
        out.append({
            "code": "STRUCTURE_PREDICTED",
            "severity": "WARNING",
            "message": "The computed axis ran on a predicted model, not an experimental structure; treat it as flagged, not measured (rule 4e).",
            "field_path": "output.structure.tier",
        })

    # A general tractability caveat, if present.
    caveat = _get(dossier, "tractability", "caveat")
    if isinstance(caveat, str) and caveat.strip():
        out.append({
            "code": "TRACTABILITY_CAVEAT",
            "severity": "INFO",
            "message": caveat.strip(),
            "field_path": "output.tractability.caveat",
        })

    # UNTAGGED_VALUE for any metric with neither evidence nor assumption links.
    for m in metrics:
        if not m.get("evidence_ids") and not m.get("assumption_ids"):
            out.append({
                "code": "UNTAGGED_VALUE",
                "severity": "INFO",
                "message": f"Metric {m['id']} carries no evidence or assumption link.",
                "field_path": None,
            })
    return out


def _counterfactuals(dossier: dict) -> list[dict]:
    out: list[dict] = []
    ne = dossier.get("next_experiment") or {}
    if isinstance(ne.get("description"), str) and ne["description"].strip():
        out.append({
            "change": ne.get("description"),
            "result": ne.get("resolves") or "Would resolve the open question the verdict rests on.",
            "meaning": ne.get("rationale") or "The experiment most likely to move the verdict.",
        })
    for chk in _get(dossier, "falsification", "checks_run") or []:
        if isinstance(chk, str) and chk.strip():
            out.append({
                "change": f"If the falsification check '{chk}' had found a problem",
                "result": "the precedent/tractability claim would be weakened or withdrawn",
                "meaning": "Falsification checks are the conditions under which the conclusion changes.",
            })
    if not out:
        out.append({
            "change": "No counterfactual recorded",
            "result": "unchanged",
            "meaning": "No next experiment or falsification check was available to vary.",
        })
    return out


def _lineage(dossier: dict) -> list[dict]:
    out = [{
        "output_path": "output.verdict",
        "input_paths": ["output.target_precedent", "output.tractability"],
        "transformation": "Two independent axes combined by rule (rules 4/4a/4b); the axes are reported separately and never averaged.",
    }]
    rank = _get(dossier, "tractability", "site_pocket_rank") or {}
    if rank.get("structure_pdb_id"):
        out.append({
            "output_path": "output.tractability.site_pocket_rank.fpocket",
            "input_paths": ["output.structure.pdb_id", "output.tractability.site_hypothesis_basis"],
            "transformation": "fpocket within-structure ranking at D=1.6 over the selected structure's pockets.",
        })
    if _isnum(_get(dossier, "tractability", "pocket_volume_a3", "primary_d1_6_a3")):
        out.append({
            "output_path": "output.tractability.pocket_volume_a3.primary_d1_6_a3",
            "input_paths": ["output.structure.pdb_id"],
            "transformation": "mdpocket site volume at clustering D=1.6.",
        })
    return out


# ---- extensions: preserve the two axes separately + trace + identifiers ---- #
def _axis_retrieved(dossier: dict) -> dict:
    tp = dossier.get("target_precedent") or {}
    vb = dossier.get("verdict_basis")
    has = bool(
        (isinstance(tp.get("distinct_actives"), int) and tp["distinct_actives"] > 0)
        or _isnum(tp.get("best_potency_nm"))
        or (isinstance(tp.get("approved_small_molecules_count"), int) and tp["approved_small_molecules_count"] > 0)
        or (tp.get("clinical_stage_small_molecules") or [])
    )
    status = "supports" if (vb in ("retrieved_precedent", "both") or has) else "insufficient"
    return {
        "id": "retrieved_precedent",
        "title": "Retrieved precedent",
        "status": status,
        "chembl_target_id": tp.get("chembl_target_id"),
        "best_potency_nm": tp.get("best_potency_nm"),
        "best_potency_characterised": tp.get("best_potency_characterised"),
        "distinct_actives": tp.get("distinct_actives"),
        "approved_small_molecules_count": tp.get("approved_small_molecules_count"),
        "clinical_stage_small_molecules": len(tp.get("clinical_stage_small_molecules") or []),
        "terminated_programs": len(tp.get("terminated_programs") or []),
    }


def _axis_computed(dossier: dict) -> dict:
    tract = dossier.get("tractability") or {}
    rank = tract.get("site_pocket_rank") or {}
    site_basis = tract.get("site_hypothesis_basis")
    not_est = isinstance(site_basis, str) and "not_established" in site_basis.lower()
    has_rank = isinstance(rank.get("fpocket"), int) and isinstance(rank.get("n_pockets"), int) and bool(rank.get("structure_pdb_id"))
    return {
        "id": "computed_tractability",
        "title": "Computed tractability",
        "status": "not_run" if (not tract or not_est or not has_rank) else "insufficient",
        "site_hypothesis_basis": site_basis,
        "structure_pdb_id": rank.get("structure_pdb_id"),
        "fpocket_rank": rank.get("fpocket"),
        "prank_rank": rank.get("prank"),
        "n_pockets": rank.get("n_pockets"),
        "pocket_volume_d16_a3": _get(tract, "pocket_volume_a3", "primary_d1_6_a3"),
        "ligand_site_jaccard": tract.get("ligand_site_jaccard"),
        "cryptic_pocket_risk": tract.get("cryptic_pocket_risk"),
        "load_bearing": _get(tract, "pocket_druggability", "load_bearing"),
    }


def _figure_ext(dossier: dict) -> dict:
    pdb = _get(dossier, "tractability", "site_pocket_rank", "structure_pdb_id")
    if _get(dossier, "tractability", "mdpocket_site_definition_used") == "site_from_ligand" and pdb:
        return {"kind": "pocket_render", "path": f"figures/{pdb}_pocket.png",
                "caption": f"{pdb} ligand-anchored pocket", "present": False,
                "note": "Reference only; the image file is not packaged with this dossier."}
    return {"kind": "none", "path": None, "caption": None, "present": False,
            "note": "No figure applies to this run."}


def _trace(dossier: dict) -> list[dict]:
    tier = _get(dossier, "structure", "tier")
    rank = _get(dossier, "tractability", "site_pocket_rank") or {}
    tp = dossier.get("target_precedent") or {}
    fp, n, rpdb = rank.get("fpocket"), rank.get("n_pockets"), rank.get("structure_pdb_id")
    scan = f"rank {fp} of {n} in {rpdb}" if isinstance(fp, int) and isinstance(n, int) and rpdb else "no reportable site"
    return [
        {"stage": "graph-intake", "order": 1, "summary": f"Resolved {_get(dossier, 'target', 'uniprot_accession')}", "status": "ok"},
        {"stage": "structure-select", "order": 2, "summary": f"{tier} — {_get(dossier, 'structure', 'pdb_id')}", "status": "ok" if tier and tier != "none" else "not_run"},
        {"stage": "pocket-scan", "order": 3, "summary": scan, "status": "ok"},
        {"stage": "precedent-lookup", "order": 4,
         "summary": (f"{tp.get('distinct_actives')} actives" if tp.get("distinct_actives") is not None else "Actives not retrieved"), "status": "ok"},
        {"stage": "falsification", "order": 5,
         "summary": ("survived" if _get(dossier, "falsification", "survived") else "checks recorded"), "status": "ok"},
        {"stage": "assemble-dossier", "order": 6, "summary": f"{dossier.get('verdict')} / {dossier.get('verdict_basis')}", "status": "ok"},
    ]


def _cache_key(dossier: dict) -> str:
    """Hash of the COMPLETE input, not just the accession (rule: cache binding)."""
    inp = dossier.get("input") or {}
    chains = _get(dossier, "tractability", "method", "chains_used")
    canonical = {
        "uniprot_accession": inp.get("uniprot_accession"),
        "as_of_date": inp.get("as_of_date"),
        "disease_context": inp.get("disease_context"),
        "interaction_to_disrupt": inp.get("interaction_to_disrupt"),
        "mechanism_hypothesis": inp.get("mechanism_hypothesis"),
        "chains_used": chains,
    }
    blob = json.dumps(canonical, sort_keys=True, separators=(",", ":"))
    return "sha256:" + hashlib.sha256(blob.encode("utf-8")).hexdigest()


def _identifiers(dossier: dict) -> dict:
    pdbs: list[str] = []
    for p in (
        _get(dossier, "structure", "pdb_id"),
        _get(dossier, "tractability", "site_pocket_rank", "structure_pdb_id"),
        *(_get(dossier, "tractability", "method", "ensemble_pdb_ids") or []),
    ):
        if isinstance(p, str) and p and p not in pdbs:
            pdbs.append(p)
    return {
        "pdb_ids": pdbs,
        "chembl_target_id": _get(dossier, "target_precedent", "chembl_target_id"),
        "uniprot_accession": _get(dossier, "target", "uniprot_accession"),
        "tools": [t for t in [_get(dossier, "tractability", "method", "tool"), "fpocket", "foldseek"] if isinstance(t, str)],
    }


def _extensions(dossier: dict) -> dict:
    return {
        "axes": {
            "retrieved_precedent": _axis_retrieved(dossier),
            "computed_tractability": _axis_computed(dossier),
            "atomistic_simulation": {"status": "not_run", "note": "No molecular-dynamics / atomistic free-energy simulation is run by this module."},
        },
        "trace": _trace(dossier),
        "next_experiment": dossier.get("next_experiment"),
        "figure": _figure_ext(dossier),
        "identifiers": _identifiers(dossier),
        "cache_key": _cache_key(dossier),
    }


def _headline(dossier: dict) -> dict:
    verdict = dossier.get("verdict")
    tp = dossier.get("target_precedent") or {}
    bits: list[str] = []
    if isinstance(tp.get("distinct_actives"), int):
        bits.append(f"{tp['distinct_actives']} actives")
    elif tp.get("distinct_actives") is None:
        bits.append("actives not retrieved")
    if _isnum(tp.get("best_potency_nm")):
        bits.append(f"best {_fmt(tp['best_potency_nm'])} nM")
    clin = tp.get("clinical_stage_small_molecules") or []
    if clin:
        bits.append(f"{len(clin)} clinical-stage small molecules")
    acc = _get(dossier, "target", "uniprot_accession") or _get(dossier, "input", "uniprot_accession")
    label = _VERDICT_LABELS.get(verdict, verdict or "Verdict")

    basis: list[str] = []
    if bits and any(k in (tp or {}) for k in ("distinct_actives", "best_potency_nm", "approved_small_molecules_count")):
        basis.append("OBSERVED")
    if _get(dossier, "tractability", "site_pocket_rank", "fpocket") is not None or _get(dossier, "structure", "tier") not in (None, "none"):
        basis.append("MODELED")
    if _get(dossier, "structure", "tier") == "predicted" and "MODELED" not in basis:
        basis.append("MODELED")
    if not basis:
        basis = ["OBSERVED"]  # we observed the databases, even if empty

    plain = f"{label} for {acc}" + (": " + "; ".join(bits) + "." if bits else ".")
    return {
        "title": f"Small-molecule tractability: {acc}" if acc else "Small-molecule tractability",
        "result": verdict,
        "plain_language": plain,
        "status": _STATUS.get(verdict, "INCONCLUSIVE"),
        "basis": basis,
    }


def build_interpretability(dossier: dict) -> dict:
    """Map a druggability dossier to the LABrador common interpretability contract.

    Pure and deterministic. Reads only fields present on ``dossier``; unknowns stay
    null and earn a limitation. Validates against schemas/interpretability.schema.json
    for any schema-valid dossier, including abstaining/degraded runs.
    """
    if not isinstance(dossier, dict):
        raise TypeError("build_interpretability expects a dossier dict")

    evidence = _evidence(dossier)
    assumptions = _assumptions(dossier)
    evidence_ids = {e["id"] for e in evidence}
    metrics = _metrics(dossier, evidence_ids)
    metric_ids = {m["id"] for m in metrics}
    steps = _steps(dossier, evidence_ids)
    limitations = _limitations(dossier, metrics)

    # An empty metrics/evidence set is allowed only with a limitation explaining it.
    if not metrics:
        limitations.append({
            "code": "NO_METRICS",
            "severity": "INFO",
            "message": "No reportable metrics for this run (e.g. an insufficient-evidence abstention); see the not-retrieved limitations.",
            "field_path": None,
        })
    if not evidence:
        limitations.append({
            "code": "NO_EVIDENCE",
            "severity": "WARNING",
            "message": "No structured evidence sources were retrieved for this run.",
            "field_path": None,
        })

    return {
        "schema_version": SCHEMA_VERSION,
        "headline": _headline(dossier),
        "metrics": metrics,
        "steps": steps,
        "evidence": evidence,
        "assumptions": assumptions,
        "uncertainty": _uncertainty(dossier, metric_ids),
        "limitations": limitations,
        "counterfactuals": _counterfactuals(dossier),
        "lineage": _lineage(dossier),
        "extensions": _extensions(dossier),
    }
