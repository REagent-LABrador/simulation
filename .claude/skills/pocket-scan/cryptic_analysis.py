"""Cryptic-pocket mechanism classification from an apo/holo structure pair.

A site that is invisible in an apo structure is invisible for one of several
physically distinct reasons, and those reasons carry different priors on
potency.  This module measures which one applies.

The public entry point is :func:`analyze_cryptic_mechanism`.

What it measures
----------------
1. A core C-alpha superposition of holo onto apo that *excludes the mobile
   region*, so the moving part cannot drag the fit.
2. Max backbone C-alpha displacement at the site.
3. Clash attribution -- the holo ligand is placed in the apo frame and every
   heavy-atom contact under a cutoff is attributed to backbone atoms of a
   retained chain, side-chain atoms of a retained chain, or atoms of a chain
   the ligand would have to displace.
4. A self-control: the same test run against the holo structure itself.  This
   must come back near zero.  If it does not, the superposition or the ligand
   placement is broken and every other number here is meaningless.
5. Ligand free-volume fraction in the apo frame, recomputed with the displaced
   chain removed and with clashing side chains trimmed to CB, so the
   contribution of each obstruction is separable.

Mechanism classes
-----------------
``loop_or_backbone_motion``
    The backbone itself moved.  The site is genuinely absent in apo.
``sidechain_occlusion``
    The backbone is where it should be; side chains sit in the site.
``subunit_occlusion``
    Another chain of the oligomer sits in the site.
``none``
    The site is open in apo.

Why the distinction is a prior on potency, not a label: across the CryptoSite
set, 25 of 27 loop-motion sites reached nanomolar, while every side-chain-motion
site with affinity data bound low-micromolar at best.  Side chains reorient on
~10^-11 s and simply compete with the ligand; loops move slowly enough that a
ligand can wedge them open.

Crypticity
----------
Mechanism and crypticity are different questions and are reported separately.
Vajda's definition: a site is cryptic only if the pocket is absent in all or
nearly all unbound structures.  CryptoBench operationalises this as
pocket-residue RMSD > 2 A.  A site that is fully formed in apo but has another
subunit or a rotamer sitting in it is *occluded*, not cryptic, and the
escalation it needs is different.

Dependencies: gemmi and numpy only.
"""

from __future__ import annotations

import itertools
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

try:  # pragma: no cover - import guard
    import gemmi
except ImportError as exc:  # pragma: no cover
    raise ImportError(
        "cryptic_analysis requires gemmi (pip install gemmi). gemmi is used "
        "rather than a legacy PDB parser because it handles mmCIF, 5-character "
        "chemical component IDs and multi-character chain IDs."
    ) from exc


__all__ = [
    "analyze_cryptic_mechanism",
    "BACKBONE_ATOMS",
    "VDW_RADII",
    "MECHANISMS",
]


# --------------------------------------------------------------------------
# constants
# --------------------------------------------------------------------------

BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

#: Bondi van der Waals radii, angstrom.
VDW_RADII: Dict[str, float] = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80,
    "S": 1.80, "CL": 1.75, "BR": 1.85, "I": 1.98, "SE": 1.90, "B": 1.92,
    "SI": 2.10, "NA": 2.27, "MG": 1.73, "K": 2.75, "CA": 2.31, "ZN": 1.39,
    "FE": 2.05, "MN": 2.05, "CU": 1.40, "NI": 1.63, "CO": 2.00,
}
_DEFAULT_VDW = 1.70

MECHANISMS = (
    "loop_or_backbone_motion",
    "sidechain_occlusion",
    "subunit_occlusion",
    "none",
)

# Fraction of the vdW radius used as the steric-obstruction criterion in the
# free-volume calculation. 0.85 permits the ~15% vdW interpenetration that real
# contacts show; at a hard 1.0 every genuine van der Waals contact reads as an
# obstruction and free volume is systematically understated.
_DEFAULT_VDW_SCALE = 0.85


# --------------------------------------------------------------------------
# structure I/O
# --------------------------------------------------------------------------

def _load_structure(path: str) -> "gemmi.Structure":
    """Read a PDB or mmCIF file and normalise it for geometric comparison."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    return st


def _is_amino_acid(resname: str) -> bool:
    info = gemmi.find_tabulated_residue(resname)
    return bool(info and info.is_amino_acid())


class _Atoms:
    """Flat, vectorised view of one model of a structure."""

    __slots__ = ("chain", "resi", "resn", "name", "elem", "xyz", "is_aa", "n")

    def __init__(self, st: "gemmi.Structure", model: int = 0) -> None:
        chain: List[str] = []
        resi: List[int] = []
        resn: List[str] = []
        name: List[str] = []
        elem: List[str] = []
        xyz: List[Tuple[float, float, float]] = []
        is_aa: List[bool] = []
        aa_cache: Dict[str, bool] = {}
        for ch in st[model]:
            for res in ch:
                aa = aa_cache.get(res.name)
                if aa is None:
                    aa = _is_amino_acid(res.name)
                    aa_cache[res.name] = aa
                for at in res:
                    chain.append(ch.name)
                    resi.append(res.seqid.num)
                    resn.append(res.name)
                    name.append(at.name)
                    elem.append(at.element.name.upper())
                    xyz.append((at.pos.x, at.pos.y, at.pos.z))
                    is_aa.append(aa)
        self.chain = np.array(chain, dtype=object)
        self.resi = np.array(resi, dtype=int)
        self.resn = np.array(resn, dtype=object)
        self.name = np.array(name, dtype=object)
        self.elem = np.array(elem, dtype=object)
        self.xyz = np.array(xyz, dtype=float).reshape(-1, 3)
        self.is_aa = np.array(is_aa, dtype=bool)
        self.n = len(self.chain)

    def mask(self, m: np.ndarray) -> "_AtomSubset":
        return _AtomSubset(self, np.asarray(m, dtype=bool))

    def protein(self, chains: Optional[Sequence[str]] = None) -> "_AtomSubset":
        m = self.is_aa.copy()
        if chains is not None:
            m &= np.isin(self.chain.astype(str), list(chains))
        return self.mask(m)

    def component(self, comp_id: str,
                  chains: Optional[Sequence[str]] = None) -> "_AtomSubset":
        m = self.resn.astype(str) == str(comp_id)
        if chains is not None:
            m &= np.isin(self.chain.astype(str), list(chains))
        return self.mask(m)

    def chain_names(self) -> List[str]:
        seen: List[str] = []
        for c in self.chain.astype(str):
            if c not in seen:
                seen.append(c)
        return seen


class _AtomSubset:
    __slots__ = ("parent", "idx")

    def __init__(self, parent: _Atoms, m: np.ndarray) -> None:
        self.parent = parent
        self.idx = np.where(m)[0]

    def __len__(self) -> int:
        return len(self.idx)

    @property
    def xyz(self) -> np.ndarray:
        return self.parent.xyz[self.idx]

    @property
    def chain(self) -> np.ndarray:
        return self.parent.chain[self.idx]

    @property
    def resi(self) -> np.ndarray:
        return self.parent.resi[self.idx]

    @property
    def resn(self) -> np.ndarray:
        return self.parent.resn[self.idx]

    @property
    def name(self) -> np.ndarray:
        return self.parent.name[self.idx]

    @property
    def elem(self) -> np.ndarray:
        return self.parent.elem[self.idx]

    def radii(self, scale: float = 1.0) -> np.ndarray:
        return np.array([VDW_RADII.get(str(e), _DEFAULT_VDW) for e in self.elem]) * scale

    def subset(self, m: np.ndarray) -> "_AtomSubset":
        out = _AtomSubset.__new__(_AtomSubset)
        out.parent = self.parent
        out.idx = self.idx[np.asarray(m, dtype=bool)]
        return out


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------

def _kabsch(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Rigid transform mapping P onto Q.  Returns (R, t) with ``P @ R.T + t ~= Q``."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def _apply(R: np.ndarray, t: np.ndarray, X: np.ndarray) -> np.ndarray:
    return np.asarray(X, dtype=float).reshape(-1, 3) @ R.T + t


def _rmsd(A: np.ndarray, B: np.ndarray) -> float:
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def _pairwise_min(A: np.ndarray, B: np.ndarray, chunk: int = 512) -> np.ndarray:
    """Per-row minimum distance from A to B."""
    out = np.full(len(A), np.inf)
    for i in range(0, len(B), chunk):
        d = np.linalg.norm(A[:, None, :] - B[None, i:i + chunk, :], axis=2)
        out = np.minimum(out, d.min(1))
    return out


def _ca_map(atoms: _Atoms, chain: str) -> Dict[int, np.ndarray]:
    m = (atoms.name.astype(str) == "CA") & (atoms.chain.astype(str) == str(chain)) & atoms.is_aa
    idx = np.where(m)[0]
    return {int(atoms.resi[i]): atoms.xyz[i] for i in idx}


def _res_map(atoms: _Atoms, chain: str) -> Dict[int, str]:
    m = (atoms.name.astype(str) == "CA") & (atoms.chain.astype(str) == str(chain)) & atoms.is_aa
    idx = np.where(m)[0]
    return {int(atoms.resi[i]): str(atoms.resn[i]) for i in idx}


# --------------------------------------------------------------------------
# chain correspondence
# --------------------------------------------------------------------------

def _pair_coords(holo: _Atoms, apo: _Atoms,
                 mapping: Sequence[Tuple[str, str]],
                 match_names: bool = True,
                 ) -> Tuple[np.ndarray, np.ndarray, List[Tuple[str, str, int]],
                            List[Dict[str, Any]]]:
    """Equivalent CA coordinate pairs for a holo-chain -> apo-chain mapping.

    Residues are matched on author sequence number, and (by default) rejected
    when the residue names disagree, so point mutants and engineered
    constructs between the two entries drop out of the fit instead of quietly
    biasing it.  The mismatches are always reported, whether or not they are
    excluded, because they tell you the two entries are not the same construct.
    """
    P: List[np.ndarray] = []
    Q: List[np.ndarray] = []
    keys: List[Tuple[str, str, int]] = []
    mismatches: List[Dict[str, Any]] = []
    for hc, ac in mapping:
        hm, am = _ca_map(holo, hc), _ca_map(apo, ac)
        hn, an = _res_map(holo, hc), _res_map(apo, ac)
        for r in sorted(set(hm) & set(am)):
            if hn[r] != an[r]:
                mismatches.append({"holo_chain": hc, "apo_chain": ac, "resi": r,
                                   "holo_resn": hn[r], "apo_resn": an[r]})
                if match_names:
                    continue
            P.append(hm[r])
            Q.append(am[r])
            keys.append((hc, ac, r))
    if not P:
        return np.zeros((0, 3)), np.zeros((0, 3)), keys, mismatches
    return np.array(P), np.array(Q), keys, mismatches


def _best_chain_mapping(holo: _Atoms, apo: _Atoms,
                        holo_chains: Sequence[str],
                        apo_chains: Sequence[str],
                        max_permutations: int = 5040,
                        match_names: bool = True,
                        ) -> Tuple[List[Tuple[str, str]], float, List[Dict[str, Any]]]:
    """Choose the apo chain for each holo chain by lowest CA RMSD.

    Exhaustive over permutations when the assembly is small; greedy per-chain
    otherwise, so this does not blow up on large oligomers.
    """
    k = len(holo_chains)
    if k > len(apo_chains):
        raise ValueError(
            f"apo has {len(apo_chains)} candidate chains but the holo assembly "
            f"needs {k}; pass apo_chains= explicitly."
        )

    n_perm = 1
    for i in range(k):
        n_perm *= (len(apo_chains) - i)

    tried: List[Dict[str, Any]] = []
    if n_perm <= max_permutations:
        best: Optional[Tuple[float, List[Tuple[str, str]]]] = None
        for perm in itertools.permutations(apo_chains, k):
            mapping = list(zip(holo_chains, perm))
            P, Q, _, _ = _pair_coords(holo, apo, mapping, match_names)
            if len(P) < 3:
                continue
            R, t = _kabsch(P, Q)
            r = _rmsd(_apply(R, t, P), Q)
            tried.append({"mapping": {h: a for h, a in mapping},
                          "n_ca": int(len(P)), "rmsd": round(r, 3)})
            if best is None or r < best[0]:
                best = (r, mapping)
        if best is None:
            raise ValueError("no usable chain correspondence between apo and holo")
        return best[1], best[0], tried

    # greedy fallback for large assemblies
    mapping: List[Tuple[str, str]] = []
    remaining = list(apo_chains)
    for hc in holo_chains:
        scored = []
        for ac in remaining:
            P, Q, _, _ = _pair_coords(holo, apo, [(hc, ac)], match_names)
            if len(P) < 3:
                continue
            R, t = _kabsch(P, Q)
            scored.append((_rmsd(_apply(R, t, P), Q), ac))
        if not scored:
            raise ValueError(f"no apo chain matches holo chain {hc}")
        scored.sort()
        mapping.append((hc, scored[0][1]))
        remaining.remove(scored[0][1])
    P, Q, _, _ = _pair_coords(holo, apo, mapping, match_names)
    R, t = _kabsch(P, Q)
    r = _rmsd(_apply(R, t, P), Q)
    tried.append({"mapping": {h: a for h, a in mapping},
                  "n_ca": int(len(P)), "rmsd": round(r, 3), "greedy": True})
    return mapping, r, tried


# --------------------------------------------------------------------------
# superposition with mobile-region exclusion
# --------------------------------------------------------------------------

def _normalise_exclusions(spec: Any, holo_chains: Sequence[str]
                          ) -> Dict[str, set]:
    """Accept ``{chain: iterable[int]}``, ``iterable[int]``, or None."""
    out: Dict[str, set] = {c: set() for c in holo_chains}
    if spec is None:
        return out
    if isinstance(spec, Mapping):
        for c, rs in spec.items():
            out.setdefault(str(c), set()).update(int(r) for r in rs)
    else:
        rs = {int(r) for r in spec}
        for c in holo_chains:
            out[c] = set(rs)
    return out


def _core_superposition(holo: _Atoms, apo: _Atoms,
                        mapping: Sequence[Tuple[str, str]],
                        ligand_xyz: np.ndarray,
                        exclude_residues: Any,
                        exclude_radius: Optional[float],
                        auto_trim: bool,
                        trim_k: float,
                        trim_floor: float,
                        fit_residue_range: Optional[Tuple[int, int]],
                        min_fit_fraction: float,
                        match_names: bool,
                        ) -> Dict[str, Any]:
    """Fit holo onto apo on core C-alpha, excluding the mobile region.

    The mobile region is the union of three optional sources:

    * ``exclude_residues`` -- residues the caller knows to be mobile
      (e.g. the switch regions of a GTPase), in holo numbering;
    * ``exclude_radius``   -- residues whose CA lies within this distance of
      the ligand, i.e. the site itself;
    * ``auto_trim``        -- iterative rejection of displacement outliers,
      which finds mobile regions the caller did not name.
    """
    holo_chains = [hc for hc, _ in mapping]
    P, Q, keys, mismatches = _pair_coords(holo, apo, mapping, match_names)
    if len(P) < 3:
        raise ValueError("fewer than 3 equivalent C-alpha atoms; cannot superpose")

    excl = _normalise_exclusions(exclude_residues, holo_chains)
    keep = np.ones(len(P), dtype=bool)
    reasons: Dict[str, int] = {"explicit": 0, "near_site": 0, "outlier": 0}

    for i, (hc, _ac, r) in enumerate(keys):
        if r in excl.get(hc, ()):
            keep[i] = False
            reasons["explicit"] += 1
        elif fit_residue_range is not None and not (
                fit_residue_range[0] <= r <= fit_residue_range[1]):
            keep[i] = False
            reasons["explicit"] += 1

    if exclude_radius is not None and len(ligand_xyz):
        d = _pairwise_min(P, ligand_xyz)
        near = d <= float(exclude_radius)
        reasons["near_site"] += int((near & keep).sum())
        keep &= ~near

    if keep.sum() < 3:
        raise ValueError(
            "mobile-region exclusion removed nearly every C-alpha; loosen "
            "exclude_radius / exclude_residues"
        )

    n_iter = 0
    if auto_trim:
        sel = np.where(keep)[0]
        for n_iter in range(1, 21):
            R, t = _kabsch(P[sel], Q[sel])
            dev_all = np.linalg.norm(_apply(R, t, P) - Q, axis=1)
            rms = float(np.sqrt((dev_all[sel] ** 2).mean()))
            cutoff = max(trim_floor, trim_k * rms)
            new = np.array([i for i in sel if dev_all[i] <= cutoff])
            if len(new) < max(3, int(min_fit_fraction * len(P))) or len(new) == len(sel):
                break
            sel = new
        reasons["outlier"] += int(keep.sum() - len(sel))
        keep = np.zeros(len(P), dtype=bool)
        keep[sel] = True

    fit_idx = np.where(keep)[0]
    R, t = _kabsch(P[fit_idx], Q[fit_idx])
    core_rmsd = _rmsd(_apply(R, t, P[fit_idx]), Q[fit_idx])
    all_rmsd = _rmsd(_apply(R, t, P), Q)

    # diagnostic: the naive fit on every equivalent CA, no exclusion at all
    R0, t0 = _kabsch(P, Q)
    naive_rmsd = _rmsd(_apply(R0, t0, P), Q)

    excluded = [{"holo_chain": keys[i][0], "apo_chain": keys[i][1], "resi": keys[i][2]}
                for i in range(len(P)) if not keep[i]]

    return {
        "R": R,
        "t": t,
        "n_equivalent_ca": int(len(P)),
        "n_fitted_ca": int(len(fit_idx)),
        "core_ca_rmsd": round(core_rmsd, 3),
        "all_ca_rmsd_after_core_fit": round(all_rmsd, 3),
        "all_ca_rmsd_naive_fit": round(naive_rmsd, 3),
        "n_excluded_ca": int(len(P) - len(fit_idx)),
        "excluded_by": reasons,
        "excluded_residues": excluded,
        "trim_iterations": n_iter,
        "chain_mapping": {h: a for h, a in mapping},
        "match_residue_names": bool(match_names),
        "n_residue_name_mismatches": len(mismatches),
        "residue_name_mismatches": mismatches,
    }


# --------------------------------------------------------------------------
# contacts
# --------------------------------------------------------------------------

def _classify_atom(chain: str, atom_name: str, displaced: Sequence[str]) -> str:
    if chain in displaced:
        return "displaced_chain"
    if atom_name in BACKBONE_ATOMS:
        return "backbone"
    return "sidechain"


def _contact_report(prot: _AtomSubset, prot_xyz: np.ndarray,
                    lig_xyz: np.ndarray, cutoff: float,
                    displaced: Sequence[str]) -> Dict[str, Any]:
    """Attribute every protein/ligand heavy-atom contact under ``cutoff``."""
    if len(prot) == 0 or len(lig_xyz) == 0:
        return {
            "cutoff": cutoff, "contact_pairs": 0, "n_ligand_atoms": 0,
            "n_protein_atoms": 0, "min_distance": None,
            "by_category": {"backbone": 0, "sidechain": 0, "displaced_chain": 0},
            "residues": [],
        }
    D = np.linalg.norm(prot_xyz[:, None, :] - lig_xyz[None, :, :], axis=2)
    hits = np.argwhere(D < cutoff)
    cats = {"backbone": 0, "sidechain": 0, "displaced_chain": 0}
    per_res: Dict[Tuple[str, str, int, str], Dict[str, Any]] = {}
    chains = prot.chain.astype(str)
    names = prot.name.astype(str)
    resns = prot.resn.astype(str)
    resis = prot.resi
    for i, j in hits:
        cat = _classify_atom(chains[i], names[i], displaced)
        cats[cat] += 1
        key = (cat, chains[i], int(resis[i]), resns[i])
        rec = per_res.setdefault(
            key, {"category": cat, "chain": chains[i], "resi": int(resis[i]),
                  "resn": resns[i], "atoms": set(), "n_contacts": 0,
                  "min_distance": np.inf})
        rec["atoms"].add(names[i])
        rec["n_contacts"] += 1
        rec["min_distance"] = min(rec["min_distance"], float(D[i, j]))
    residues = []
    for key in sorted(per_res, key=lambda k: (k[0], k[1], k[2])):
        rec = per_res[key]
        rec["atoms"] = sorted(rec["atoms"])
        rec["min_distance"] = round(rec["min_distance"], 2)
        residues.append(rec)
    return {
        "cutoff": cutoff,
        "contact_pairs": int(len(hits)),
        "n_ligand_atoms": int((D.min(0) < cutoff).sum()),
        "n_protein_atoms": int((D.min(1) < cutoff).sum()),
        "min_distance": round(float(D.min()), 2),
        "by_category": cats,
        "residues": residues,
    }


# --------------------------------------------------------------------------
# free volume
# --------------------------------------------------------------------------

def _ligand_volume_points(lig_xyz: np.ndarray, lig_elem: Sequence[str],
                          spacing: float) -> Tuple[np.ndarray, float]:
    """Grid points inside the ligand's own van der Waals envelope."""
    radii = np.array([VDW_RADII.get(str(e), _DEFAULT_VDW) for e in lig_elem])
    lo = lig_xyz.min(0) - (radii.max() + spacing)
    hi = lig_xyz.max(0) + (radii.max() + spacing)
    axes = [np.arange(lo[i], hi[i] + 1e-9, spacing) for i in range(3)]
    G = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    inside = np.zeros(len(G), dtype=bool)
    for i in range(0, len(lig_xyz), 256):
        d = np.linalg.norm(G[:, None, :] - lig_xyz[None, i:i + 256, :], axis=2)
        inside |= (d < radii[None, i:i + 256]).any(1)
    return G[inside], spacing ** 3


def _free_fraction(points: np.ndarray, prot_xyz: np.ndarray,
                   prot_radii: np.ndarray) -> float:
    """Fraction of ligand-volume points not inside any protein vdW sphere."""
    if len(points) == 0:
        return float("nan")
    if len(prot_xyz) == 0:
        return 1.0
    pad = float(prot_radii.max()) + 1.0
    lo, hi = points.min(0) - pad, points.max(0) + pad
    m = (prot_xyz > lo).all(1) & (prot_xyz < hi).all(1)
    X, Rr = prot_xyz[m], prot_radii[m]
    if len(X) == 0:
        return 1.0
    blocked = np.zeros(len(points), dtype=bool)
    for i in range(0, len(X), 256):
        d = np.linalg.norm(points[:, None, :] - X[None, i:i + 256, :], axis=2)
        blocked |= (d < Rr[None, i:i + 256]).any(1)
    return float(1.0 - blocked.mean())


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------

def _classify(contacts: Dict[str, Any],
              max_site_ca_disp: float,
              free_fraction_apo: Optional[float],
              backbone_motion_threshold: float,
              subunit_fraction_threshold: float,
              open_free_fraction: float,
              ) -> Dict[str, Any]:
    """Assign the mechanism.

    Order matters, and the order encodes a physical claim:

    1. No contacts at all (or a site already essentially free) -> ``none``.
    2. A chain the ligand would displace dominates -> ``subunit_occlusion``.
       This is checked before backbone motion because a displaced subunit
       contributes backbone atoms of its own, which would otherwise be
       misread as local backbone motion.
    3. C-alpha displacement at the site above threshold ->
       ``loop_or_backbone_motion``.  **Displacement, not the presence of
       backbone atoms in the contact list, is the discriminator.**  A loop
       that swings 9 A carries its side chains with it, so the atoms actually
       occupying the site can be entirely side-chain even though the cause is
       backbone motion.  Counting backbone contact atoms instead would
       misclassify exactly the cases this module exists to catch.
    4. Otherwise -> ``sidechain_occlusion``.
    """
    cats = contacts["by_category"]
    total = contacts["contact_pairs"]
    n_bb, n_sc, n_dp = cats["backbone"], cats["sidechain"], cats["displaced_chain"]

    site_open = total == 0 or (
        free_fraction_apo is not None and free_fraction_apo >= open_free_fraction
    )

    frac_displaced = (n_dp / total) if total else 0.0
    ranked = sorted(
        [("subunit_occlusion", n_dp), ("sidechain_occlusion", n_sc),
         ("loop_or_backbone_motion", n_bb)],
        key=lambda kv: -kv[1])

    if site_open:
        mechanism = "none"
        rationale = (
            f"{total} heavy-atom contacts under {contacts['cutoff']} A"
            + ("" if free_fraction_apo is None
               else f"; ligand free volume in apo {free_fraction_apo * 100:.1f}%")
            + " - the site is open in apo."
        )
    elif frac_displaced >= subunit_fraction_threshold:
        mechanism = "subunit_occlusion"
        rationale = (
            f"{n_dp} of {total} contacts ({frac_displaced * 100:.0f}%) come from a "
            f"chain the ligand would displace; max site C-alpha displacement "
            f"{max_site_ca_disp:.2f} A."
        )
    elif max_site_ca_disp > backbone_motion_threshold:
        mechanism = "loop_or_backbone_motion"
        rationale = (
            f"max site C-alpha displacement {max_site_ca_disp:.2f} A exceeds "
            f"{backbone_motion_threshold} A; the backbone itself moved "
            f"({n_bb} backbone / {n_sc} side-chain contact atoms)."
        )
    else:
        mechanism = "sidechain_occlusion"
        rationale = (
            f"max site C-alpha displacement {max_site_ca_disp:.2f} A is below "
            f"{backbone_motion_threshold} A and {n_sc} of {total} contacts are "
            f"side-chain atoms; the backbone is in place."
        )

    secondary = None
    if mechanism != "none":
        for name, count in ranked:
            if name != mechanism and count > 0:
                secondary = name
                break

    return {
        "mechanism": mechanism,
        "secondary_mechanism": secondary,
        "rationale": rationale,
        "contact_fraction_displaced_chain": round(frac_displaced, 3),
        "counts": {"backbone": n_bb, "sidechain": n_sc,
                   "displaced_chain": n_dp, "total": total},
    }


# --------------------------------------------------------------------------
# public API
# --------------------------------------------------------------------------

def analyze_cryptic_mechanism(
    apo_path: str,
    holo_path: str,
    ligand_comp_id: str,
    *,
    holo_chains: Optional[Sequence[str]] = None,
    apo_chains: Optional[Sequence[str]] = None,
    ligand_chain: Optional[str] = None,
    exclude_residues: Any = None,
    exclude_radius: Optional[float] = None,
    auto_trim: bool = True,
    trim_k: float = 2.5,
    trim_floor: float = 1.0,
    fit_residue_range: Optional[Tuple[int, int]] = None,
    min_fit_fraction: float = 0.5,
    match_residue_names: bool = True,
    min_chain_contact_fraction: float = 0.15,
    site_radius: float = 5.0,
    clash_cutoff: float = 2.0,
    wide_cutoff: float = 3.0,
    backbone_motion_threshold: float = 2.0,
    subunit_fraction_threshold: float = 0.5,
    cryptic_rmsd_threshold: float = 2.0,
    open_free_fraction: float = 0.95,
    compute_free_volume: bool = True,
    grid_spacing: float = 0.4,
    vdw_scale: float = _DEFAULT_VDW_SCALE,
    self_control_max_contacts: int = 3,
    pocket_present_in_apo_fraction: Optional[float] = None,
) -> Dict[str, Any]:
    """Classify why a ligand-binding site is not visible in an apo structure.

    Parameters
    ----------
    apo_path, holo_path
        Structure files (PDB or mmCIF; gemmi reads both).
    ligand_comp_id
        Chemical component ID of the holo ligand, e.g. ``"MOV"`` or ``"307"``.
        Note that PDB entry titles often name a compound differently from its
        component ID -- use the component ID.
    holo_chains
        Chains forming the holo biological assembly.  Default: the chains that
        contact the chosen ligand copy within ``site_radius``.  An asymmetric
        unit is not a biological assembly; state the assembly deliberately.
    apo_chains
        Chains of the apo assembly to consider.  Default: all protein chains.
        Any apo chain not mapped to a holo chain is treated as a chain the
        ligand would have to displace -- so restrict this if the apo file
        contains unrelated crystallographic copies.
    min_chain_contact_fraction
        Only used when ``holo_chains`` is inferred.  A contacting chain is kept
        only if it contributes at least this fraction of the atom count of the
        largest contributor, which drops chains that merely brush the ligand
        across a crystal contact.  State ``holo_chains`` explicitly instead
        wherever you can.
    ligand_chain
        Which copy of the ligand to use when the file has several.  Default:
        the copy with the most protein contacts.
    exclude_residues
        Mobile-region residues to keep out of the superposition, in holo
        numbering.  Either an iterable of residue numbers (applied to every
        chain) or ``{chain: iterable}``.
    exclude_radius
        Also exclude residues whose CA is within this distance of the ligand.
        ``None`` (default) disables.
    auto_trim
        Iteratively reject displacement outliers from the fit.  This finds
        mobile regions the caller did not name.
    match_residue_names
        Keep out of the superposition any position where the apo and holo
        residue names disagree, i.e. mutations and engineered constructs.
        Default True.  The mismatching positions are reported either way as
        ``superposition.residue_name_mismatches``; a non-empty list means the
        two entries are not the same construct, which is worth knowing.  Set
        False to reproduce a hand fit that matched on residue number alone.
    cryptic_rmsd_threshold
        CryptoBench's operational definition of cryptic: pocket-residue RMSD
        above this (A) is cryptic.
    pocket_present_in_apo_fraction
        Optional.  Fraction of unbound structures in which the pocket is
        already present.  Vajda's definition requires the pocket to be absent
        in all or nearly all unbound structures, so a value above 0.5 forces
        ``is_cryptic=False`` regardless of the pairwise RMSD.

    Returns
    -------
    dict
        JSON-serialisable.  Key fields: ``mechanism``, ``is_cryptic``,
        ``self_control.passed``.  **If ``self_control.passed`` is False the
        result is not interpretable** -- the superposition or the ligand
        placement is broken, and every downstream number is noise.
    """
    holo_st, apo_st = _load_structure(holo_path), _load_structure(apo_path)
    holo, apo = _Atoms(holo_st), _Atoms(apo_st)

    # ---- ligand -------------------------------------------------------
    lig_all = holo.component(ligand_comp_id)
    if len(lig_all) == 0:
        present = sorted({str(r) for r in holo.resn[~holo.is_aa]})
        raise ValueError(
            f"component {ligand_comp_id!r} not found in {holo_path}. "
            f"Non-polymer components present: {present}"
        )
    copies: Dict[str, np.ndarray] = {}
    for c in sorted({str(x) for x in lig_all.chain}):
        copies[c] = np.where(lig_all.chain.astype(str) == c)[0]

    holo_prot_all = holo.protein()
    if ligand_chain is not None:
        if str(ligand_chain) not in copies:
            raise ValueError(
                f"ligand {ligand_comp_id} has no copy in chain {ligand_chain}; "
                f"copies are in chains {sorted(copies)}"
            )
        chosen = str(ligand_chain)
    elif len(copies) == 1:
        chosen = next(iter(copies))
    else:
        best_c, best_n = None, -1
        for c, sel in copies.items():
            L = lig_all.xyz[sel]
            d = _pairwise_min(holo_prot_all.xyz, L)
            n = int((d < site_radius).sum())
            if n > best_n:
                best_c, best_n = c, n
        chosen = best_c

    lig = lig_all.subset(lig_all.chain.astype(str) == chosen)
    lig_xyz = lig.xyz
    lig_elem = [str(e) for e in lig.elem]

    # ---- assemblies ---------------------------------------------------
    # An asymmetric unit is not a biological assembly. Auto-derivation below is
    # a convenience, not a substitute for stating the assembly: 2AZ5's ASU is
    # two independent TNF-alpha dimers, and the ligand in chain A brushes three
    # atoms of chain D across a crystal contact. Taking every contacting chain
    # would build an A/B/D "assembly" that has no biological meaning and, when
    # mapped onto the apo trimer, leaves no chain free to be displaced -- which
    # silently converts subunit_occlusion into loop_or_backbone_motion.
    assembly_inferred = holo_chains is None
    dropped_minor: List[Dict[str, Any]] = []
    if holo_chains is None:
        d = _pairwise_min(holo_prot_all.xyz, lig_xyz)
        near = np.where(d < site_radius)[0]
        counts: Dict[str, int] = {}
        for i in near:
            c = str(holo_prot_all.chain[i])
            counts[c] = counts.get(c, 0) + 1
        if not counts:
            raise ValueError("no holo protein chain contacts the ligand")
        top = max(counts.values())
        holo_chains = sorted(c for c, n in counts.items()
                             if n >= min_chain_contact_fraction * top)
        dropped_minor = [
            {"chain": c, "n_atoms_near_ligand": n,
             "fraction_of_top": round(n / top, 3)}
            for c, n in sorted(counts.items()) if c not in holo_chains]
    holo_chains = [str(c) for c in holo_chains]
    if not holo_chains:
        raise ValueError("no holo protein chain contacts the ligand")

    apo_all_chains = [c for c in apo.chain_names()
                      if bool(apo.is_aa[apo.chain.astype(str) == c].any())]
    if apo_chains is None:
        apo_chains = apo_all_chains
    apo_chains = [str(c) for c in apo_chains]

    mapping, map_rmsd, tried = _best_chain_mapping(
        holo, apo, holo_chains, apo_chains, match_names=match_residue_names)
    mapped_apo = [a for _, a in mapping]
    unmapped_chains = [c for c in apo_chains if c not in mapped_apo]

    # ---- superposition ------------------------------------------------
    sup = _core_superposition(
        holo, apo, mapping, lig_xyz, exclude_residues, exclude_radius,
        auto_trim, trim_k, trim_floor, fit_residue_range, min_fit_fraction,
        match_residue_names)
    R, t = sup.pop("R"), sup.pop("t")
    sup["best_mapping_rmsd_all_ca"] = round(map_rmsd, 3)
    sup["mappings_tried"] = tried

    lig_in_apo = _apply(R, t, lig_xyz)

    # A chain only counts as one the ligand *displaces* if the ligand actually
    # runs into it. An apo file often carries extra crystallographic copies of
    # the same protomer (4OBE has two); calling those "displaced subunits"
    # would be a mislabel, and on a file where such a copy happened to sit near
    # the site it would fabricate a subunit_occlusion call.
    displaced_chains: List[str] = []
    bystander_chains: List[Dict[str, Any]] = []
    for c in unmapped_chains:
        sub = apo.protein([c])
        dmin = (float(_pairwise_min(lig_in_apo, sub.xyz).min())
                if len(sub) else float("inf"))
        if dmin < wide_cutoff:
            displaced_chains.append(c)
        else:
            bystander_chains.append(
                {"chain": c, "min_distance_to_ligand": (
                    None if dmin == float("inf") else round(dmin, 2))})

    # ---- site definition and C-alpha displacement ---------------------
    holo_prot = holo.protein(holo_chains)
    d_site = _pairwise_min(holo_prot.xyz, lig_xyz)
    site_keys = {(str(holo_prot.chain[i]), int(holo_prot.resi[i]))
                 for i in np.where(d_site < site_radius)[0]}
    h2a = {h: a for h, a in mapping}
    per_res_disp = []
    for hc, r in sorted(site_keys):
        hm, am = _ca_map(holo, hc), _ca_map(apo, h2a[hc])
        if r in hm and r in am:
            dd = float(np.linalg.norm(_apply(R, t, hm[r][None])[0] - am[r]))
            per_res_disp.append({"holo_chain": hc, "apo_chain": h2a[hc],
                                 "resi": r, "ca_displacement": round(dd, 2)})
    devs = np.array([p["ca_displacement"] for p in per_res_disp]) if per_res_disp \
        else np.zeros(0)
    max_disp = float(devs.max()) if len(devs) else 0.0
    site_rmsd = float(np.sqrt((devs ** 2).mean())) if len(devs) else 0.0
    arg_max = (max(per_res_disp, key=lambda p: p["ca_displacement"])
               if per_res_disp else None)

    # ---- clash attribution --------------------------------------------
    apo_prot = apo.protein(apo_chains)
    contacts = _contact_report(apo_prot, apo_prot.xyz, lig_in_apo,
                               clash_cutoff, displaced_chains)
    contacts_wide = _contact_report(apo_prot, apo_prot.xyz, lig_in_apo,
                                    wide_cutoff, displaced_chains)

    # ---- self-control --------------------------------------------------
    self_contacts = _contact_report(holo_prot, holo_prot.xyz, lig_xyz,
                                    clash_cutoff, [])
    self_wide = _contact_report(holo_prot, holo_prot.xyz, lig_xyz,
                                wide_cutoff, [])
    passed = self_contacts["contact_pairs"] <= self_control_max_contacts
    self_control = {
        "description": ("the same clash test run against the holo structure "
                        "itself; must be near zero or the superposition / "
                        "ligand placement is broken"),
        "contact_pairs": self_contacts["contact_pairs"],
        "n_protein_atoms": self_contacts["n_protein_atoms"],
        "n_ligand_atoms": self_contacts["n_ligand_atoms"],
        "min_distance": self_contacts["min_distance"],
        "contact_pairs_wide": self_wide["contact_pairs"],
        "n_protein_atoms_wide": self_wide["n_protein_atoms"],
        "residues": self_contacts["residues"],
        "threshold": self_control_max_contacts,
        "passed": bool(passed),
    }

    # ---- free volume ---------------------------------------------------
    free_volume: Dict[str, Any] = {"computed": False}
    free_apo = None
    if compute_free_volume:
        pts_holo, cell = _ligand_volume_points(lig_xyz, lig_elem, grid_spacing)
        pts_apo, _ = _ligand_volume_points(lig_in_apo, lig_elem, grid_spacing)
        lig_volume = len(pts_holo) * cell

        f_holo = _free_fraction(pts_holo, holo_prot.xyz, holo_prot.radii(vdw_scale))
        f_apo = _free_fraction(pts_apo, apo_prot.xyz, apo_prot.radii(vdw_scale))

        keep_no_sub = ~np.isin(apo_prot.chain.astype(str), displaced_chains)
        sub_removed = apo_prot.subset(keep_no_sub)
        f_no_sub = _free_fraction(pts_apo, sub_removed.xyz,
                                  sub_removed.radii(vdw_scale))

        clash_res = {(r["chain"], r["resi"]) for r in contacts["residues"]
                     if r["category"] == "sidechain"}
        keep_sc = np.array([
            (c, int(i)) not in clash_res or n in BACKBONE_ATOMS or n == "CB"
            for c, i, n in zip(apo_prot.chain.astype(str), apo_prot.resi,
                               apo_prot.name.astype(str))])
        trimmed_only = apo_prot.subset(keep_sc)
        f_trim_only = _free_fraction(pts_apo, trimmed_only.xyz,
                                     trimmed_only.radii(vdw_scale))
        both = apo_prot.subset(keep_sc & keep_no_sub)
        f_both = _free_fraction(pts_apo, both.xyz, both.radii(vdw_scale))

        free_volume = {
            "computed": True,
            "ligand_vdw_volume_A3": round(lig_volume, 1),
            "grid_spacing": grid_spacing,
            "vdw_scale": vdw_scale,
            "holo": round(f_holo, 4),
            "apo_intact": round(f_apo, 4),
            "apo_minus_displaced_chains": round(f_no_sub, 4),
            "apo_minus_clashing_sidechains": round(f_trim_only, 4),
            "apo_minus_displaced_chains_and_sidechains": round(f_both, 4),
            "trimmed_sidechain_residues": sorted(
                f"{c}/{i}" for c, i in clash_res),
            "note": ("fraction of the ligand's own van der Waals volume that is "
                     "sterically unobstructed; protein radii scaled by "
                     f"{vdw_scale}"),
        }
        free_apo = f_apo

    # ---- classification ------------------------------------------------
    cls = _classify(contacts, max_disp, free_apo, backbone_motion_threshold,
                    subunit_fraction_threshold, open_free_fraction)

    # ---- crypticity -----------------------------------------------------
    is_cryptic = site_rmsd > cryptic_rmsd_threshold
    reason = (
        f"pocket-residue C-alpha RMSD {site_rmsd:.2f} A "
        f"{'>' if is_cryptic else '<='} {cryptic_rmsd_threshold} A "
        f"(CryptoBench criterion); max displacement {max_disp:.2f} A"
    )
    if pocket_present_in_apo_fraction is not None:
        if pocket_present_in_apo_fraction > 0.5:
            is_cryptic = False
            reason += (
                f"; overridden to not-cryptic because the pocket is already "
                f"present in {pocket_present_in_apo_fraction * 100:.0f}% of "
                f"unbound structures (Vajda: cryptic requires absence in all "
                f"or nearly all unbound structures)"
            )
    if not is_cryptic and cls["mechanism"] != "none":
        reason += (
            f" - the site is not cryptic but is OCCLUDED "
            f"({cls['mechanism']}): present in apo, with something standing in it"
        )

    result: Dict[str, Any] = {
        "mechanism": cls["mechanism"],
        "secondary_mechanism": cls["secondary_mechanism"],
        "rationale": cls["rationale"],
        "is_cryptic": bool(is_cryptic),
        "crypticity": {
            "is_cryptic": bool(is_cryptic),
            "site_ca_rmsd": round(site_rmsd, 2),
            "max_site_ca_displacement": round(max_disp, 2),
            "threshold": cryptic_rmsd_threshold,
            "criterion": ("CryptoBench: pocket-residue RMSD > 2 A. Vajda: "
                          "absent in all or nearly all unbound structures."),
            "reason": reason,
        },
        "inputs": {
            "apo": str(apo_path), "holo": str(holo_path),
            "ligand_comp_id": str(ligand_comp_id),
            "ligand_chain": chosen,
            "ligand_heavy_atoms": int(len(lig_xyz)),
            "ligand_copies_in_holo": sorted(copies),
            "holo_chains": holo_chains,
            "holo_assembly_inferred": bool(assembly_inferred),
            "chains_dropped_as_crystal_contacts": dropped_minor,
            "apo_chains": apo_chains,
            "displaced_apo_chains": displaced_chains,
            "unmapped_apo_chains": unmapped_chains,
            "bystander_apo_chains": bystander_chains,
        },
        "superposition": sup,
        "site": {
            "radius": site_radius,
            "n_residues": len(site_keys),
            "n_residues_compared": len(per_res_disp),
            "max_ca_displacement": round(max_disp, 2),
            "max_ca_displacement_at": arg_max,
            "ca_rmsd": round(site_rmsd, 2),
            "per_residue": per_res_disp,
        },
        "contacts": contacts,
        "contacts_wide": contacts_wide,
        "self_control": self_control,
        "free_volume": free_volume,
        "counts": cls["counts"],
    }
    warnings: List[str] = []
    if not passed:
        warnings.append(
            f"SELF-CONTROL FAILED: {self_contacts['contact_pairs']} contacts "
            f"under {clash_cutoff} A of the ligand with its own holo structure "
            f"(threshold {self_control_max_contacts}). The superposition or the "
            f"ligand placement is broken and this result is not interpretable."
        )
    if assembly_inferred:
        warnings.append(
            f"holo_chains was inferred from ligand contacts as {holo_chains}, "
            f"not stated. An asymmetric unit is not a biological assembly; "
            f"confirm the assembly and pass holo_chains explicitly."
            + (f" Dropped as crystal contacts: {dropped_minor}."
               if dropped_minor else "")
        )
    if len(unmapped_chains) == 0 and len(apo_chains) == len(holo_chains):
        warnings.append(
            "every apo chain was mapped to a holo chain, so no chain is "
            "available to be displaced; subunit_occlusion cannot be detected "
            "with this chain selection."
        )
    if bystander_chains:
        warnings.append(
            f"apo chains {[b['chain'] for b in bystander_chains]} are unmapped "
            f"but do not reach the ligand, so they are reported as bystanders "
            f"(likely extra crystallographic copies) rather than displaced "
            f"subunits."
        )
    if warnings:
        result["warnings"] = warnings
    return result


# --------------------------------------------------------------------------

def _main(argv: Optional[Sequence[str]] = None) -> int:  # pragma: no cover
    import argparse
    import json

    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("apo")
    p.add_argument("holo")
    p.add_argument("ligand_comp_id")
    p.add_argument("--holo-chains", nargs="*", default=None)
    p.add_argument("--apo-chains", nargs="*", default=None)
    p.add_argument("--ligand-chain", default=None)
    p.add_argument("--exclude", nargs="*", type=int, default=None)
    p.add_argument("--exclude-radius", type=float, default=None)
    p.add_argument("--no-trim", action="store_true")
    p.add_argument("--no-free-volume", action="store_true")
    a = p.parse_args(argv)

    res = analyze_cryptic_mechanism(
        a.apo, a.holo, a.ligand_comp_id,
        holo_chains=a.holo_chains, apo_chains=a.apo_chains,
        ligand_chain=a.ligand_chain, exclude_residues=a.exclude,
        exclude_radius=a.exclude_radius, auto_trim=not a.no_trim,
        compute_free_volume=not a.no_free_volume)
    print(json.dumps(res, indent=2, default=str))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(_main())
