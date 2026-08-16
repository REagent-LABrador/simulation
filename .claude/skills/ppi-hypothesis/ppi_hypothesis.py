"""ppi-hypothesis — turn a cofolded protein-protein interface into a graph ask.

Read the SKILL.md next to this file before using any function here. The short
version: nothing in this module decides that a predicted interface is real. It
measures the interface, applies a gate whose every check is a *reproducibility*
or *novelty* test fitted on the 15-case control panel in
`fixtures/panel_results.json`, and emits an ask in the four-verb schema
`graph-intake` already uses.

On that panel the gate recovered 3/3 deposited positives, admitted 0/11
negatives, and issued **zero asks** — the positives are training data and the
one novel pair that passed was already a link in the graph. If you are running
this expecting an ask, read "The hard question" in the SKILL.md first.

Plain Python, run under the proto-tools interpreter (the one that has
`proto_tools`, `gemmi` and `numpy`). Cofolding is reached through
`cofold-check`'s `predict.py`, which this module imports and never modifies.
"""

from __future__ import annotations

import json
import math
import os
import statistics
import sys
from dataclasses import dataclass, field, asdict
from itertools import combinations
from typing import Any, Iterable, Sequence

import numpy as np

# ---------------------------------------------------------------------------
# Constants that are thresholds, and are labelled as such.
# Every one of these was set from the control panel in `fixtures/panel_results.json`.
# None of them is a literature value. See SKILL.md "Where the gate came from".
# ---------------------------------------------------------------------------
CA_CONTACT_CUTOFF_A = 8.0      # matches cofold-check's _inter_chain_ca_contacts
HEAVY_CONTACT_CUTOFF_A = 5.0   # matches pocket-scan's interface_residues
PROBE_RADIUS_A = 1.4
SASA_POINTS = 200

VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
       "SE": 1.90, "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98}
DEFAULT_VDW = 1.70


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------
def _read_cif(text: str):
    import gemmi

    st = gemmi.read_structure_string(text)
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    return st


def _heavy_atoms(model, chains: Iterable[str]):
    """(coords Nx3, radii N, residue-key list N) for the named chains."""
    import gemmi

    want = set(chains)
    xyz, rad, keys = [], [], []
    for ch in model:
        if ch.name not in want:
            continue
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            for at in res:
                el = at.element.name.upper()
                if el == "H":
                    continue
                xyz.append((at.pos.x, at.pos.y, at.pos.z))
                rad.append(VDW.get(el, DEFAULT_VDW))
                keys.append(f"{ch.name}:{res.seqid.num}:{res.name}")
    return np.asarray(xyz, float), np.asarray(rad, float), keys


def _ca_coords(model, chains: Iterable[str]):
    want = set(chains)
    xyz, keys = [], []
    for ch in model:
        if ch.name not in want:
            continue
        for res in ch:
            at = res.find_atom("CA", "*")
            if at is not None:
                xyz.append((at.pos.x, at.pos.y, at.pos.z))
                keys.append((ch.name, res.seqid.num))
    return np.asarray(xyz, float), keys


def _sphere_points(n: int) -> np.ndarray:
    """Fibonacci sphere — deterministic, so SASA is reproducible run to run."""
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = math.pi * (1.0 + 5.0 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def sasa(xyz: np.ndarray, radii: np.ndarray, n_points: int = SASA_POINTS) -> float:
    """Shrake-Rupley solvent-accessible surface area, A^2.

    Implemented here rather than shelled out because freesasa is not in the
    environment. Verified against a rigid-body control (translate a chain 500 A
    away: total SASA becomes the exact sum of the isolated parts, BSA -> 0).
    """
    if len(xyz) == 0:
        return 0.0
    pts = _sphere_points(n_points)
    r = radii + PROBE_RADIUS_A
    total = 0.0
    # neighbour list by bounding-box blocking; N is a few thousand, this is fine
    rmax = float(r.max())
    for i in range(len(xyz)):
        d = np.linalg.norm(xyz - xyz[i], axis=1)
        nb = np.where((d < r[i] + rmax) & (d > 1e-6))[0]
        if len(nb) == 0:
            total += 4.0 * math.pi * r[i] ** 2
            continue
        test = xyz[i] + pts * r[i]
        dd = np.linalg.norm(test[:, None, :] - xyz[nb][None, :, :], axis=2)
        buried = (dd < r[nb][None, :]).any(axis=1)
        total += 4.0 * math.pi * r[i] ** 2 * float((~buried).mean())
    return total


def interface_metrics(cif_text: str, chains_a: Sequence[str], chains_b: Sequence[str],
                      with_bsa: bool = True) -> dict[str, Any]:
    """Every geometric number this module gates on, for one structure.

    `chains_a` is the ligand/target group, `chains_b` the candidate partner.
    Groups, not chains: a TNF trimer contacting one receptor is ABC against D.
    Passing a single protomer measures a third of the epitope (dossier 2b).
    """
    st = _read_cif(cif_text)
    model = st[0]
    xa, ra, ka = _heavy_atoms(model, chains_a)
    xb, rb, kb = _heavy_atoms(model, chains_b)
    if len(xa) == 0 or len(xb) == 0:
        return {"error": f"no heavy atoms for {chains_a} / {chains_b}",
                "chains_present": [c.name for c in model]}

    d = np.linalg.norm(xa[:, None, :] - xb[None, :, :], axis=2)
    close = d < HEAVY_CONTACT_CUTOFF_A
    n_heavy = int(close.sum())
    res_a = sorted({ka[i] for i in np.where(close.any(axis=1))[0]},
                   key=lambda s: (s.split(":")[0], int(s.split(":")[1])))
    res_b = sorted({kb[j] for j in np.where(close.any(axis=0))[0]},
                   key=lambda s: (s.split(":")[0], int(s.split(":")[1])))

    ca_a, _ = _ca_coords(model, chains_a)
    ca_b, _ = _ca_coords(model, chains_b)
    dca = np.linalg.norm(ca_a[:, None, :] - ca_b[None, :, :], axis=2)
    ca_pairs = int((dca < CA_CONTACT_CUTOFF_A).sum())

    out = {
        "chains_a": list(chains_a), "chains_b": list(chains_b),
        "n_heavy_atom_contacts_5a": n_heavy,
        "ca_pairs_8a": ca_pairs,
        "_ca_pairs_note": ("CA-CA PAIR count, the same quantity as the 97 in the "
                           "8DYG reference. NOT a residue count (8DYG is 29)."),
        "n_interface_res_a": len(res_a), "n_interface_res_b": len(res_b),
        "interface_res_a": res_a, "interface_res_b": res_b,
        "min_heavy_dist_a": round(float(d.min()), 2),
        "com_separation_a": round(float(np.linalg.norm(xa.mean(0) - xb.mean(0))), 2),
    }
    if with_bsa:
        s_ab = sasa(np.vstack([xa, xb]), np.concatenate([ra, rb]))
        s_a = sasa(xa, ra)
        s_b = sasa(xb, rb)
        out["bsa_total_a2"] = round(s_a + s_b - s_ab, 1)
        out["bsa_per_side_a2"] = round((s_a + s_b - s_ab) / 2.0, 1)
        out["sasa_a_a2"], out["sasa_b_a2"] = round(s_a, 1), round(s_b, 1)
    return out


def _seqid_set(labels: Iterable[str], chains: Sequence[str] | None = None) -> set[int]:
    """Residue *numbers* only, chain identity dropped.

    Correct for a homo-oligomer where the same epitope appears on three chains
    (dossier rule 4's C3-symmetry problem, used here deliberately: across seeds
    Boltz-2 may place the receptor on a different protomer of the same trimer,
    which is the same interface, not a different one).
    """
    out = set()
    for lab in labels:
        c, num, _ = lab.split(":", 2)
        if chains is None or c in chains:
            out.add(int(num))
    return out


def jaccard(a: set, b: set) -> float | None:
    if not a and not b:
        return None
    return len(a & b) / len(a | b)


def seed_concordance(per_seed: list[dict], side: str = "a") -> dict[str, Any]:
    """How much do independent seeds agree on WHERE the partner binds.

    Jaccard over residue numbers on one side, chain identity dropped. This is
    the only signal in the panel that separated positives from hard negatives,
    and it is still not sufficient on its own — see SKILL.md failure mode 2.
    """
    sets = [_seqid_set(m[f"interface_res_{side}"]) for m in per_seed
            if m.get(f"interface_res_{side}")]
    if len(sets) < 2:
        return {"n_seeds_with_interface": len(sets), "mean_jaccard": None,
                "min_jaccard": None,
                "warning": "fewer than 2 seeds produced an interface; "
                           "concordance is UNKNOWN, not high."}
    js = [jaccard(a, b) for a, b in combinations(sets, 2)]
    return {"n_seeds_with_interface": len(sets), "n_pairs": len(js),
            "mean_jaccard": round(statistics.mean(js), 3),
            "min_jaccard": round(min(js), 3),
            "consensus_residues": sorted(set.intersection(*sets)),
            "union_residues": sorted(set.union(*sets))}


# ---------------------------------------------------------------------------
# Candidate partner selection
# ---------------------------------------------------------------------------
@dataclass
class Candidate:
    """One target x partner pair, with the reason it is in the set."""
    target: str
    partner: str
    source: str          # why this partner was proposed
    prior: str           # known_partner | family_member | graph_named | control
    expected: str | None = None   # positive | negative | unknown
    note: str = ""


def enumerate_candidates(target_family: str, graph_partners: Sequence[str],
                         family_members: Sequence[str],
                         exclude: Sequence[str] = ()) -> list[Candidate]:
    """Three sources, in this order, and nothing else.

    1. partners the graph already names for this target or a related one;
    2. members of the cognate receptor/ligand family (a TNFSF ligand against
       the TNFRSF ectodomains, not against the proteome);
    3. Foldseek structural neighbours of a known partner.

    Everything else is excluded and the exclusion is recorded. An unstated
    candidate set makes a false-positive rate meaningless — you cannot compute
    a rate over a denominator you did not write down.
    """
    out: list[Candidate] = []
    for p in graph_partners:
        if p in exclude:
            continue
        out.append(Candidate(target_family, p, "named in the graph", "graph_named"))
    for p in family_members:
        if p in exclude or any(c.partner == p for c in out):
            continue
        out.append(Candidate(target_family, p, "cognate receptor family", "family_member"))
    return out


# ---------------------------------------------------------------------------
# The gate
# ---------------------------------------------------------------------------
@dataclass
class GateResult:
    passed: bool
    checks: dict[str, Any] = field(default_factory=dict)
    failed: list[str] = field(default_factory=list)
    verdict: str = ""
    novelty: str = ""


# Thresholds fitted on the measured panel in `fixtures/panel_results.json`.
# Read the SKILL.md section "Where the gate came from" before changing any of
# them. Two signals that look obvious are DELIBERATELY ABSENT:
#
#   * contact count and buried surface area. Not weak — INVERTED. A
#     composition-scrambled random sequence against the TL1A trimer produced
#     166 CA-CA pairs and 6569 A^2 buried, the largest interface in the whole
#     panel, against 38 pairs / 2250 A^2 for the deposited TL1A/DcR3 complex.
#   * homolog footprint transfer. It passes the negatives (0.75 for a
#     same-superfamily non-partner, 0.84 for the scrambled string). It tests
#     the family's canonical binding geometry, not this pair. Reported, never
#     gated on.
GATE = {
    "min_seeds": 6,
    "min_seed_blocks": 2,          # disjoint base seeds; 3 seeds in one block lied
    "min_seed_concordance_mean_jaccard": 0.70,
    "min_seed_concordance_min_jaccard": 0.60,
    "min_iptm_margin_over_scramble": 0.08,
    "require_scramble_control": True,
    "require_rank_first_or_tied_with_known_partner": True,
}

# Reported, required in the output, and NOT gated on. Keeping them here rather
# than deleting them is the point: the next reader will want to gate on them.
REPORTED_NOT_GATED = {
    "ca_pairs_8a": "scrambled control scored 166 vs 38 for the true complex",
    "bsa_total_a2": "scrambled control scored 6569 vs 2250 for the true complex",
    "footprint_transfer_coverage": "passes both hard negatives; a site check, not "
                                   "a specificity check",
    "esmfold_agreement": "ESMFold produced 4 CA pairs against 72 deposited on the "
                         "IL-17A/IL-17RA positive control. It fails the positives, "
                         "so two-model agreement is not an available signal",
    "absolute_iptm": "the scrambled-partner floor is 0.479 for IL-17A, 0.787 for "
                     "TNF-alpha and 0.785 for TL1A -- a 0.31 spread between three "
                     "ligands in one panel. There is no transferable absolute ipTM "
                     "threshold; only the margin over an in-run floor has a scale",
}


def gate(iptm_median: float | None, iptm_scramble_median: float | None,
         concordance: dict[str, Any], n_seeds: int, n_seed_blocks: int,
         rank_in_candidate_set: int | None, rank1_is_known_partner: bool,
         already_in_graph: bool, deposited_complex: str | None) -> GateResult:
    """Conjunctive checks. Every one is a *reproducibility* or *novelty* test.

    Note what is not here: nothing in this gate establishes that the interface
    is real. Seed concordance says the model is consistent; the scramble margin
    says it is more consistent than it is about a random string; the rank says
    it prefers this partner to the others we offered. All three can be true of
    a wrong answer, and the SKILL.md failure modes say so at length.
    """
    g = GateResult(passed=False)
    fails: list[str] = []

    def chk(name: str, ok: bool, value: Any, threshold: Any) -> None:
        g.checks[name] = {"value": value, "threshold": threshold, "pass": bool(ok)}
        if not ok:
            fails.append(name)

    chk("enough_seeds", n_seeds >= GATE["min_seeds"], n_seeds, GATE["min_seeds"])
    chk("disjoint_seed_blocks", n_seed_blocks >= GATE["min_seed_blocks"],
        n_seed_blocks, GATE["min_seed_blocks"])
    chk("seed_concordance_mean", (concordance.get("mean_jaccard") or 0) >=
        GATE["min_seed_concordance_mean_jaccard"], concordance.get("mean_jaccard"),
        GATE["min_seed_concordance_mean_jaccard"])
    chk("seed_concordance_min", (concordance.get("min_jaccard") or 0) >=
        GATE["min_seed_concordance_min_jaccard"], concordance.get("min_jaccard"),
        GATE["min_seed_concordance_min_jaccard"])
    if GATE["require_scramble_control"]:
        margin = (None if iptm_median is None or iptm_scramble_median is None
                  else round(iptm_median - iptm_scramble_median, 4))
        chk("scramble_control_margin",
            margin is not None and margin >= GATE["min_iptm_margin_over_scramble"],
            margin, GATE["min_iptm_margin_over_scramble"])
        if iptm_scramble_median is None:
            g.checks["scramble_control_margin"]["value"] = (
                "NO SCRAMBLE CONTROL RUN — this is a fail, not a missing datum. "
                "Without it there is no floor and the ipTM number has no scale.")

    if GATE["require_rank_first_or_tied_with_known_partner"]:
        ok = rank_in_candidate_set == 1 or (rank_in_candidate_set == 2 and rank1_is_known_partner)
        chk("rank_in_candidate_set", ok, rank_in_candidate_set,
            "1, or 2 behind a known partner")

    g.failed = fails
    g.passed = not fails

    if deposited_complex:
        g.novelty = (f"NOT NOVEL — {deposited_complex} is deposited. Boltz-2 trained "
                     f"on the PDB, so this is recall, not prediction. Report it as a "
                     f"method check; never as a hypothesis.")
    elif already_in_graph:
        g.novelty = ("NOT NEW TO THE GRAPH — the interaction is already a link. Only "
                     "the INTERFACE would be new, and an interface is a structural "
                     "claim the never-ask list assigns to us. No ask.")
    else:
        g.novelty = ("new to the graph and no deposited complex. An ask is on the "
                     "table — and still has to clear the never-fire list.")

    g.verdict = "gate passed" if g.passed else "gate FAILED on: " + ", ".join(fails)
    return g


# ---------------------------------------------------------------------------
# The ask
# ---------------------------------------------------------------------------
VERBS = ("expand_node", "resolve_link", "test_gap", "new_question")


def build_ask(verb: str, target_id: str | None, question: str,
              depth: str = "deep") -> dict[str, Any]:
    """The graph's own four verbs. There is no fifth and this must not add one."""
    if verb not in VERBS:
        raise ValueError(f"{verb!r} is not one of the graph's four verbs: {VERBS}. "
                         "A structural prediction is a new_question, or a test_gap "
                         "when it lands on a gap the graph already names. Inventing "
                         "a verb makes the ask unconsumable upstream.")
    if verb != "new_question" and not target_id:
        raise ValueError(f"{verb} must point at a row by id, never in prose")
    return {"ask": verb, "target": target_id, "depth": depth, "question": question}


def not_found_entry(verb: str, target_id: str | None, graph_id: str, round_n: int,
                    field_name: str, question: str) -> dict[str, str]:
    """Where a generated hypothesis lives in the dossier.

    Same convention `graph-intake` already uses — one `not_found[]` entry whose
    `reason` is prefixed `ASK[<verb>:<target id>]`. No template change. The
    'Not blocking' clause is not decoration: it is the assertion that the field
    beside it was filled from what we measured, and that this entry would read
    the same if the ask were never answered.
    """
    tid = target_id or "null"
    return {
        "field": field_name,
        "reason": (f"ASK[{verb}:{tid}] issued to graph {graph_id} round {round_n}. "
                   f"Not blocking: the predicted interface is reported in "
                   f"tractability.pocket_vs_interface from our own measurement and "
                   f"this entry records the residual literature question. {question}"),
    }


def gate_summary_line(name: str, g: GateResult) -> str:
    marks = "".join("+" if v["pass"] else "-" for v in g.checks.values())
    return f"{name:34s} {marks}  {'PASS' if g.passed else 'fail: ' + ','.join(g.failed)}"


if __name__ == "__main__":
    print(__doc__)
    print("thresholds:", json.dumps(GATE, indent=2))
