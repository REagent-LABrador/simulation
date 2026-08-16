"""PPI-interface analysis and pocket classification.

Answers the question the dossier's `tractability.pocket_vs_interface` block asks:
*given a pocket, what is its relationship to the interaction we were asked to
block?*  The site you block is frequently not the site the partner binds, so the
relationship is measured here rather than assumed from the mechanism hypothesis.

Four labels, matching CLAUDE.md rule 2b:

    orthosteric_candidate    pocket overlaps the partner interface
    allosteric_candidate     pocket is distal from it
    destabiliser_candidate   pocket is buried *within* the oligomer -- lined by
                             two or more subunits and enclosed, not a solvent
                             facing groove at the interface rim
    no_partner_structure     no complex was supplied, and the pocket is not an
                             enclosed inter-subunit cavity either, so nothing
                             can be said

The label is never returned on its own.  Every result carries the overlap
fraction, the enclosure, the per-chain lining composition and the distance to
the interface, because the label is a thresholded view of those numbers and the
thresholds are judgement calls.

Requires gemmi (mmCIF, five-character ligand codes such as ``A1JPS``, and
multi-character chain names -- legacy PDB column parsing handles none of these)
and numpy.  Nothing else; no GPU.

Always feed this module a **biological assembly**
(``https://files.rcsb.org/download/<ID>-assembly1.cif``), and check that the
deposited assembly is actually the biological unit: 2AZ5's ``assembly1`` is a
crystallographic tetramer of two independent TNF-alpha dimers, and scoring all
four chains fuses sites across a packing contact.  Every function here takes an
explicit ``chains`` argument for exactly that reason.

Typical use
-----------
    st = load_structure("2AZ5-assembly1.cif")
    sites = ligand_site_residues(st, "307", chains=("A", "B"))
    iface = interface_residues(st_complex, ["A", "B", "C"], ["R", "S", "T"])
    res = classify_pocket(
        sites[0].residues, iface.side_a, st,
        target_chains=("A", "B"), probe_points=sites[0].positions,
        match_by="seqid",
    )
    print(res.classification, res.overlap_fraction, res.enclosure)

Run ``python interface_analysis.py --selftest <cif_dir>`` for the fixture
harness (IL-17A, TNF-alpha, KRAS).

What was measured while building this
-------------------------------------
Enclosure, per site, from the module's own ray casting:

    SPD304 in the 1TNF apo trimer      1.000   gain from other subunits 0.393
    A1JPS in the 9SQX IL-17A dimer     0.937   gain 0.365
    sotorasib in 6VJJ KRAS (monomer)   0.918   gain 0.000
    sotorasib in 6OIM KRAS (monomer)   0.830   gain 0.000
    U5Q in the 8DYG IL-17A dimer       0.702   gain 0.202-0.213
    SPD304 in the 2AZ5 holo dimer      0.689   gain 0.293

The first surprise: **absolute enclosure does not separate an inter-subunit
cavity from an ordinary surface pocket.** The KRAS switch-II pocket -- one
chain, on the surface -- is *more* enclosed (0.83) than the TNF-alpha axial
channel in 2AZ5 (0.69), because in 2AZ5 the ligand has already pushed the third
subunit out and the channel is half open.  So the destabiliser test is built on
the *gain* in enclosure contributed by subunits other than the dominant one,
which is 0.00 for any monomeric pocket by construction.

The second surprise, found by a control that failed: **PDB entries for one
protein do not share a residue numbering.**  IL-17A entries 7UWM/8DYG/9SQX
number from the UniProt precursor, 4HSA numbers the mature chain, offset 23.
Comparing their epitopes by residue number without correcting gave Jaccard 0.22
where the truth is 0.85.  `detect_numbering_offset` exists because of that, and
should be run before any cross-entry ``match_by="seqid"`` comparison.

What this module cannot do
--------------------------
* It measures **geometry, not mechanism.**  ``destabiliser_candidate`` means
  "buried between subunits", which is a necessary but not sufficient condition
  for subunit displacement.  A ligand in such a pocket may equally distort the
  oligomer without disassembling it.
* A pocket can be both buried between subunits and under the partner epitope
  (IL-17A's dimer-groove site is), and no amount of geometry separates the two
  mechanisms.  ``also_overlaps_interface`` marks that case; it is a prompt to
  go and look, not an answer.
* The thresholds are calibrated on four sites.  The destabiliser boundary in
  particular has true positives at 0.29 and true negatives at 0.21 -- a gap of
  0.08.  Treat anything in the `borderline` band as undecided.
* Interfaces are taken from crystal or cryo-EM complexes, which are one
  conformational state each.  A partner-bound epitope in one entry is not
  necessarily the epitope in solution.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

try:
    import gemmi
except ImportError as exc:  # pragma: no cover - environment problem, not logic
    raise ImportError(
        "interface_analysis requires gemmi (conda-forge gemmi). "
        "Legacy PDB parsing cannot handle 5-character ligand codes."
    ) from exc


# ---------------------------------------------------------------------------
# Tunable constants.  Every threshold in this module is named here; none are
# buried in the code.  The rationale for each is stated because these are
# judgement calls, not measured constants, and a caller may reasonably disagree.
# ---------------------------------------------------------------------------

#: Heavy-atom distance defining the interface contact shell.  5.0 A is the
#: common convention for a PPI contact shell (roughly a first solvation shell
#: beyond van der Waals contact at ~4.0 A); it is deliberately more generous
#: than a 4.0 A contact definition because we are asking "would a ligand here
#: sit under the partner", not "is this a hydrogen bond".
CONTACT_CUTOFF_A = 5.0

#: Heavy-atom distance defining the lining of a pocket derived from a bound
#: ligand.  Tighter than the interface cutoff: these are the residues actually
#: touching the ligand.  4.5 A is the value used elsewhere in pocket-scan for
#: ligand-site definition, kept identical so Jaccard numbers are comparable.
LIGAND_SITE_CUTOFF_A = 4.5

#: Fraction of pocket-lining residues that must also be interface residues
#: before the pocket is called orthosteric.  0.25 means a quarter of the pocket
#: wall is under the partner.  Below this a ligand can plausibly bind without
#: contesting the epitope.  This is the least defensible number in the module
#: -- it is a convention, chosen so that a pocket merely *adjacent* to an
#: interface (KRAS switch-II beside the effector site) is not called
#: orthosteric, while a pocket cut into the epitope (a BH3-groove pocket) is.
ORTHOSTERIC_OVERLAP_MIN = 0.25

#: Below this overlap the pocket is confidently distal.  Between
#: ALLOSTERIC_OVERLAP_MAX and ORTHOSTERIC_OVERLAP_MIN the label is still
#: allosteric_candidate but `borderline` is set, because one or two shared
#: residues at a pocket rim is not evidence of epitope engagement and is also
#: not evidence against it.
ALLOSTERIC_OVERLAP_MAX = 0.10

#: Absolute enclosure floor (see `enclosure`).  This is a *floor*, not the
#: discriminator -- absolute enclosure turned out not to separate the cases.
#: Measured: sotorasib in the KRAS switch-II pocket, a single-chain surface
#: site, scores 0.83, while SPD304 in the 2AZ5 TNF-alpha dimer -- the textbook
#: buried inter-subunit site -- scores only 0.69, because the ligand has
#: already displaced the third subunit that would close the channel.  Any
#: absolute threshold that admits the second admits most surface pockets.  The
#: floor only excludes shallow patches.
BURIAL_ENCLOSURE_MIN = 0.60

#: The actual discriminator for "buried within the oligomer": how much of the
#: pocket's enclosure is contributed by subunits *other than* the one that
#: contributes most of it.  Computed as
#: ``enclosure(all target chains) - max_c enclosure(chain c alone)``.
#: A pocket that a single protomer would enclose on its own gains nothing; a
#: pocket cut into the oligomer core gains a lot.
#:
#: Measured on the fixtures (see module docstring / selftest):
#:     SPD304, 1TNF apo trimer        0.393
#:     A1JPS, 9SQX IL-17A dimer site  0.365
#:     SPD304, 2AZ5 holo dimer        0.293   <- lowest true positive
#:     ----------------- 0.25 -----------------
#:     U5Q, 8DYG IL-17A C-term site   0.202-0.213
#:     sotorasib, 6OIM KRAS monomer   0.000
#: The boundary therefore sits between 0.21 and 0.29, which is a gap of 0.08 on
#: a four-site calibration.  It is not comfortably separated and callers should
#: treat gains between 0.20 and 0.30 as undecided; `borderline` is set there.
SUBUNIT_ENCLOSURE_GAIN_MIN = 0.25

#: Half-width of the undecided band around `SUBUNIT_ENCLOSURE_GAIN_MIN`.
BORDERLINE_GAIN_BAND = 0.05

#: A pocket that shares no residues with the interface but sits within this
#: distance of it is *adjacent*, not distal, and gets `adjacent_to_interface`
#: set on top of the allosteric label.  The KRAS switch-II pocket is the case:
#: overlap 0.00 against the RAF1 RBD epitope, yet 1.3 A from it -- a ligand
#: there locks the inactive state next door to the effector site rather than
#: acting from across the protein.  Calling that "allosteric" without the flag
#: would be true and misleading.  5.0 A is one contact shell.
ADJACENT_DISTANCE_MAX_A = 5.0

#: Minimum share of the pocket lining a subunit must contribute before it
#: counts as genuinely lining the pocket.  Prevents a single glancing residue
#: from a neighbouring chain (or from a crystal contact that survived assembly
#: selection) from promoting a one-protomer pocket to "inter-subunit".
MULTICHAIN_MIN_SHARE = 0.20

#: Ray casting parameters for `enclosure`.  512 directions on a Fibonacci
#: sphere gives ~0.02 resolution on the fraction, well under the threshold
#: spacing.  12 A is roughly the distance at which a ligand-sized pocket's
#: walls end; a clearance of 2.6 A treats a ray as blocked when it passes
#: within about a carbon van der Waals radius plus a small probe of an atom
#: centre.
RAY_COUNT = 512
RAY_LENGTH_A = 12.0
RAY_CLEARANCE_A = 2.6
RAY_MIN_A = 1.0  # ignore the atoms the probe point is sitting on top of

#: Elements never treated as heavy atoms.
_LIGHT_ELEMENTS = {"H", "D"}

#: Residues excluded from every residue set (solvent and common cryo additives
#: that are not part of the protein surface for this purpose).
_SOLVENT = {"HOH", "DOD", "WAT"}


# ---------------------------------------------------------------------------
# Small value types
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ResidueKey:
    """A residue, labelled with its chain.  Hashable, so sets work."""

    chain: str
    seqid: int
    icode: str
    name: str

    @property
    def label(self) -> str:
        ic = self.icode.strip()
        return f"{self.chain}:{self.name}{self.seqid}{ic}"

    def __str__(self) -> str:  # pragma: no cover - cosmetic
        return self.label

    def seq_key(self) -> Tuple[int, str]:
        """Chain-agnostic key, for cross-structure / symmetric-oligomer use."""
        return (self.seqid, self.icode.strip())


@dataclass
class InterfaceResult:
    """Contact shell between two groups of chains."""

    side_a: List[ResidueKey]
    side_b: List[ResidueKey]
    chains_a: Tuple[str, ...]
    chains_b: Tuple[str, ...]
    cutoff_a: float
    n_atom_contacts: int

    @property
    def all_residues(self) -> List[ResidueKey]:
        return list(self.side_a) + list(self.side_b)

    def as_dict(self) -> dict:
        return {
            "chains_a": list(self.chains_a),
            "chains_b": list(self.chains_b),
            "cutoff_a": self.cutoff_a,
            "n_atom_contacts": self.n_atom_contacts,
            "side_a": [r.label for r in self.side_a],
            "side_b": [r.label for r in self.side_b],
            "n_side_a": len(self.side_a),
            "n_side_b": len(self.side_b),
        }


@dataclass
class LigandSite:
    """A pocket defined by a bound ligand copy."""

    comp_id: str
    ligand_key: ResidueKey
    residues: List[ResidueKey]
    positions: np.ndarray  # (n_atoms, 3) heavy-atom coordinates of the ligand

    def chain_counts(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for r in self.residues:
            out[r.chain] = out.get(r.chain, 0) + 1
        return out


@dataclass
class PocketClassification:
    """Result of `classify_pocket`.  The label is never meaningful alone."""

    classification: str
    overlap_fraction: Optional[float]
    interface_coverage: Optional[float]
    n_pocket_residues: int
    n_interface_residues: int
    n_shared_residues: int
    shared_residues: List[str]
    enclosure: Optional[float]
    enclosure_single_subunit_max: Optional[float]
    subunit_enclosure_gain: Optional[float]
    enclosure_basis: str
    lining_chains: List[str]
    lining_chain_counts: Dict[str, int]
    n_lining_subunits: int
    min_distance_to_interface_a: Optional[float]
    borderline: bool
    also_overlaps_interface: bool
    adjacent_to_interface: bool
    match_by: str
    thresholds: Dict[str, float]
    notes: List[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        d = dict(self.__dict__)
        return d

    def summary(self) -> str:
        ov = "n/a" if self.overlap_fraction is None else f"{self.overlap_fraction:.2f}"
        en = "n/a" if self.enclosure is None else f"{self.enclosure:.2f}"
        dist = (
            "n/a"
            if self.min_distance_to_interface_a is None
            else f"{self.min_distance_to_interface_a:.1f}A"
        )
        gain = (
            "n/a"
            if self.subunit_enclosure_gain is None
            else f"{self.subunit_enclosure_gain:.2f}"
        )
        flags = []
        if self.adjacent_to_interface:
            flags.append("adjacent")
        if self.also_overlaps_interface:
            flags.append("also-overlaps-interface")
        if self.borderline:
            flags.append("borderline")
        tail = ("  [" + ", ".join(flags) + "]") if flags else ""
        return (
            f"{self.classification}  overlap={ov}  enclosure={en} (gain={gain})  "
            f"subunits={self.n_lining_subunits}{sorted(self.lining_chains)}  "
            f"d(interface)={dist}{tail}"
        )


# ---------------------------------------------------------------------------
# Structure loading and atom selection
# ---------------------------------------------------------------------------


def load_structure(path: str, model_index: int = 0) -> gemmi.Structure:
    """Read an mmCIF/PDB file, strip waters and hydrogens, set up entities.

    Hydrogens are stripped because riding hydrogens in high-resolution entries
    change every distance-based cutoff.  Waters are stripped because they are
    not part of either the interface or the pocket wall for this purpose.
    """
    st = gemmi.read_structure(path)
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_waters()
    if model_index != 0:
        while len(st) > 1 and model_index > 0:
            del st[0]
            model_index -= 1
    return st


def _is_amino_acid(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_amino_acid())


def _is_nucleic(name: str) -> bool:
    info = gemmi.find_tabulated_residue(name)
    return bool(info and info.is_nucleic_acid())


def _residue_key(chain: gemmi.Chain, res: gemmi.Residue) -> ResidueKey:
    return ResidueKey(
        chain=chain.name,
        seqid=int(res.seqid.num),
        icode=(res.seqid.icode or " ").strip(),
        name=res.name,
    )


def _chain_set(chains: Optional[Iterable[str] | str]) -> Optional[set]:
    if chains is None:
        return None
    if isinstance(chains, str):
        return {chains}
    return set(chains)


def _collect_atoms(
    st: gemmi.Structure,
    chains: Optional[Iterable[str] | str] = None,
    polymer_only: bool = True,
    exclude_comp_ids: Optional[Iterable[str]] = None,
    model_index: int = 0,
) -> Tuple[np.ndarray, List[ResidueKey]]:
    """Return (coords Nx3, per-atom residue keys) for heavy atoms.

    `polymer_only` keeps amino-acid and nucleotide residues only, which is what
    you want for interface and burial calculations -- a bound ligand or a
    glycan should not be counted as a pocket wall.
    """
    want = _chain_set(chains)
    skip = {c.upper() for c in (exclude_comp_ids or ())}
    coords: List[Tuple[float, float, float]] = []
    keys: List[ResidueKey] = []
    model = st[model_index]
    for chain in model:
        if want is not None and chain.name not in want:
            continue
        for res in chain:
            if res.name in _SOLVENT or res.name.upper() in skip:
                continue
            if polymer_only and not (_is_amino_acid(res.name) or _is_nucleic(res.name)):
                continue
            key = _residue_key(chain, res)
            for atom in res:
                if atom.element.name in _LIGHT_ELEMENTS:
                    continue
                if atom.altloc not in ("\0", "", "A", " "):
                    continue
                coords.append((atom.pos.x, atom.pos.y, atom.pos.z))
                keys.append(key)
    if not coords:
        return np.zeros((0, 3), dtype=float), []
    return np.asarray(coords, dtype=float), keys


def _pairs_within(a: np.ndarray, b: np.ndarray, cutoff: float, chunk: int = 4096):
    """Yield (i_index_array, j_index_array) of pairs within `cutoff`.

    Brute force in numpy chunks.  Interfaces are a few thousand atoms a side,
    so this is milliseconds and avoids a spatial-index dependency.
    """
    if len(a) == 0 or len(b) == 0:
        return
    c2 = cutoff * cutoff
    for start in range(0, len(a), chunk):
        block = a[start : start + chunk]
        d2 = ((block[:, None, :] - b[None, :, :]) ** 2).sum(axis=-1)
        ii, jj = np.nonzero(d2 <= c2)
        if len(ii):
            yield ii + start, jj


# ---------------------------------------------------------------------------
# 1. Interface residues
# ---------------------------------------------------------------------------


def interface_residues(
    structure: gemmi.Structure,
    chain_a: Iterable[str] | str,
    chain_b: Iterable[str] | str,
    cutoff: float = CONTACT_CUTOFF_A,
    model_index: int = 0,
) -> InterfaceResult:
    """Heavy-atom contact shell between two chain groups.

    Parameters
    ----------
    structure : gemmi.Structure
        Preferably a biological assembly, already loaded with `load_structure`.
    chain_a, chain_b : chain name or iterable of chain names
        Groups, not single chains: a cytokine homodimer contacting one receptor
        is ``["A", "B"]`` against ``["C"]``.  Passing a group is the difference
        between the real epitope and half of it.
    cutoff : float
        Heavy-atom distance in Angstrom.  Default `CONTACT_CUTOFF_A` (5.0).

    Returns
    -------
    InterfaceResult with residues on each side, each labelled with its chain.

    Notes
    -----
    Only polymer residues are considered on both sides, so a glycan or a bound
    ligand at the interface does not contribute residues.  Hydrogens and
    alternate conformations B+ are excluded by `load_structure` / selection.
    """
    ca, ka = _collect_atoms(structure, chain_a, model_index=model_index)
    cb, kb = _collect_atoms(structure, chain_b, model_index=model_index)
    if len(ca) == 0 or len(cb) == 0:
        raise ValueError(
            f"no polymer heavy atoms for chains {chain_a!r} / {chain_b!r}; "
            f"check the chain names against the assembly"
        )
    side_a, side_b = set(), set()
    n_contacts = 0
    for ii, jj in _pairs_within(ca, cb, cutoff):
        n_contacts += len(ii)
        for i in np.unique(ii):
            side_a.add(ka[int(i)])
        for j in np.unique(jj):
            side_b.add(kb[int(j)])
    order = lambda r: (r.chain, r.seqid, r.icode)
    return InterfaceResult(
        side_a=sorted(side_a, key=order),
        side_b=sorted(side_b, key=order),
        chains_a=tuple(sorted(_chain_set(chain_a) or ())),
        chains_b=tuple(sorted(_chain_set(chain_b) or ())),
        cutoff_a=cutoff,
        n_atom_contacts=n_contacts,
    )


# ---------------------------------------------------------------------------
# Pocket definition helpers
# ---------------------------------------------------------------------------


def ligand_site_residues(
    structure: gemmi.Structure,
    comp_id: str,
    chains: Optional[Iterable[str] | str] = None,
    cutoff: float = LIGAND_SITE_CUTOFF_A,
    lining_chains: Optional[Iterable[str] | str] = None,
    model_index: int = 0,
) -> List[LigandSite]:
    """Every copy of ligand `comp_id`, with the residues lining it.

    `comp_id` may be five characters (``A1JPS``); gemmi handles the extended
    CCD codes that legacy PDB column parsing truncates.

    `chains` restricts which *ligand copies* are taken; `lining_chains`
    restricts which chains may contribute lining residues.  Use the latter to
    keep a crystal-packing neighbour out of the pocket wall -- in 2AZ5's
    deposited tetramer, ligand ``307`` in chain A picks up two residues from
    chain D, which belongs to the other biological dimer.
    """
    want = _chain_set(chains)
    comp = comp_id.upper()
    prot_coords, prot_keys = _collect_atoms(
        structure, lining_chains, polymer_only=True, model_index=model_index
    )
    out: List[LigandSite] = []
    model = structure[model_index]
    for chain in model:
        if want is not None and chain.name not in want:
            continue
        for res in chain:
            if res.name.upper() != comp:
                continue
            pos = np.asarray(
                [
                    (a.pos.x, a.pos.y, a.pos.z)
                    for a in res
                    if a.element.name not in _LIGHT_ELEMENTS
                    and a.altloc in ("\0", "", "A", " ")
                ],
                dtype=float,
            )
            if len(pos) == 0:
                continue
            lining = set()
            for _, jj in _pairs_within(pos, prot_coords, cutoff):
                for j in np.unique(jj):
                    lining.add(prot_keys[int(j)])
            out.append(
                LigandSite(
                    comp_id=comp,
                    ligand_key=_residue_key(chain, res),
                    residues=sorted(lining, key=lambda r: (r.chain, r.seqid, r.icode)),
                    positions=pos,
                )
            )
    return out


def residues_within(
    structure: gemmi.Structure,
    points: np.ndarray,
    cutoff: float = LIGAND_SITE_CUTOFF_A,
    chains: Optional[Iterable[str] | str] = None,
    model_index: int = 0,
) -> List[ResidueKey]:
    """Polymer residues within `cutoff` of an arbitrary point cloud.

    Use this to turn fpocket alpha-sphere centres (the ``*_vert.pqr`` file) or
    a transferred ligand into a pocket-residue list.
    """
    coords, keys = _collect_atoms(structure, chains, model_index=model_index)
    hits = set()
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    for _, jj in _pairs_within(pts, coords, cutoff):
        for j in np.unique(jj):
            hits.add(keys[int(j)])
    return sorted(hits, key=lambda r: (r.chain, r.seqid, r.icode))


# ---------------------------------------------------------------------------
# Burial / enclosure
# ---------------------------------------------------------------------------


def _fibonacci_directions(n: int) -> np.ndarray:
    """`n` roughly uniform unit vectors on the sphere."""
    i = np.arange(n, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    golden = math.pi * (1.0 + 5.0**0.5)
    theta = golden * i
    return np.stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1
    )


def enclosure(
    structure: gemmi.Structure,
    probe_points: np.ndarray,
    chains: Optional[Iterable[str] | str] = None,
    n_rays: int = RAY_COUNT,
    ray_length: float = RAY_LENGTH_A,
    clearance: float = RAY_CLEARANCE_A,
    model_index: int = 0,
) -> float:
    """Fraction of directions from the pocket that are blocked by protein.

    For each probe point (ligand heavy atoms, or alpha-sphere centres) cast
    `n_rays` rays.  A ray is blocked if some protein heavy-atom centre lies
    within `clearance` of the ray, between `RAY_MIN_A` and `ray_length` along
    it.  The returned value is the mean blocked fraction over probe points.

    1.0 is a sealed cavity, ~0.5 a flat surface patch, and a deep open groove
    lands in between.  This is a geometric enclosure, not a solvent-accessible
    surface calculation -- it is cheap, has no radii table to get wrong, and it
    is the quantity the destabiliser test actually needs: *can this pocket see
    the solvent, or is it inside the oligomer?*
    """
    coords, _ = _collect_atoms(structure, chains, polymer_only=True, model_index=model_index)
    if len(coords) == 0:
        raise ValueError("no protein atoms for enclosure calculation")
    pts = np.atleast_2d(np.asarray(probe_points, dtype=float))
    dirs = _fibonacci_directions(n_rays)
    fractions = []
    c2 = clearance * clearance
    for p in pts:
        rel = coords - p  # (N, 3)
        d2 = (rel * rel).sum(axis=1)
        near = rel[d2 <= (ray_length + clearance) ** 2]
        if len(near) == 0:
            fractions.append(0.0)
            continue
        # projection of each atom onto each ray direction
        t = near @ dirs.T  # (N_atoms, n_rays)
        perp2 = (near * near).sum(axis=1)[:, None] - t * t
        blocked = ((t >= RAY_MIN_A) & (t <= ray_length) & (perp2 <= c2)).any(axis=0)
        fractions.append(float(blocked.mean()))
    return float(np.mean(fractions))


# ---------------------------------------------------------------------------
# 2. Pocket classification
# ---------------------------------------------------------------------------


def _match_keys(keys: Sequence[ResidueKey], match_by: str):
    if match_by == "chain_seqid":
        return {(k.chain, k.seqid, k.icode) for k in keys}
    if match_by == "seqid":
        return {k.seq_key() for k in keys}
    raise ValueError(f"match_by must be 'chain_seqid' or 'seqid', got {match_by!r}")


def classify_pocket(
    pocket_residues: Sequence[ResidueKey],
    interface_residues_: Optional[Sequence[ResidueKey]],
    structure: Optional[gemmi.Structure] = None,
    target_chains: Optional[Iterable[str] | str] = None,
    probe_points: Optional[np.ndarray] = None,
    match_by: str = "chain_seqid",
    orthosteric_overlap_min: float = ORTHOSTERIC_OVERLAP_MIN,
    allosteric_overlap_max: float = ALLOSTERIC_OVERLAP_MAX,
    burial_enclosure_min: float = BURIAL_ENCLOSURE_MIN,
    subunit_enclosure_gain_min: float = SUBUNIT_ENCLOSURE_GAIN_MIN,
    multichain_min_share: float = MULTICHAIN_MIN_SHARE,
    model_index: int = 0,
) -> PocketClassification:
    """Classify a pocket by its relationship to a partner interface.

    Parameters
    ----------
    pocket_residues : sequence of ResidueKey
        Lining residues of the pocket, chain-labelled.  From
        `ligand_site_residues` (holo), `residues_within` (fpocket alpha-sphere
        centres), or any other route.
    interface_residues_ : sequence of ResidueKey or None
        The **target side** of the interface with the binding partner, i.e.
        ``InterfaceResult.side_a`` when side A was the target.  ``None`` or an
        empty sequence means no complex structure was available.
    structure : gemmi.Structure, optional
        Needed for enclosure and for the distance-to-interface measurement.
        Without it neither is computed and the destabiliser test cannot fire.
    target_chains : chains forming the target oligomer
        Enclosure is measured against these chains only, so a bound partner or
        a packing neighbour does not make a surface pocket look buried.
    probe_points : (n, 3) array, optional
        Points inside the pocket -- ligand heavy atoms or alpha-sphere centres.
        Strongly preferred.  Without them the centroid of the lining residues'
        side-chain atoms is used, which for a shallow surface groove falls
        *inside* the protein and overstates burial; the result records
        ``enclosure_basis == "lining_centroid"`` so this is visible.
    match_by : {"chain_seqid", "seqid"}
        How pocket residues are matched against interface residues.  Use
        ``"seqid"`` when the two sets come from different PDB entries, or when
        the target is a symmetric homo-oligomer and chain labelling is
        arbitrary; the price is that the interface set becomes the union over
        protomers, which inflates overlap.  Use ``"chain_seqid"`` when both
        sets were computed on the same assembly.

    Decision order (documented, and deliberately not a score)
    --------------------------------------------------------
    1. **destabiliser_candidate** -- two or more subunits each contribute at
       least `multichain_min_share` of the lining, enclosure clears the
       `burial_enclosure_min` floor, *and* the subunit enclosure gain is at
       least `subunit_enclosure_gain_min`, i.e. the other subunits are what
       make the pocket enclosed.  Checked first, and it does not need a partner
       structure, because burial inside an oligomer is measurable on the
       oligomer alone.  TNF-alpha's SPD304 channel is the case; a ligand there
       displaces a subunit rather than blocking TNF/TNFR.
    2. **no_partner_structure** -- no interface supplied.
    3. **orthosteric_candidate** -- overlap >= `orthosteric_overlap_min`.
    4. **allosteric_candidate** -- otherwise; `borderline` is set when overlap
       is between `allosteric_overlap_max` and `orthosteric_overlap_min`, and
       `adjacent_to_interface` is set when the pocket shares no residues with
       the interface but lies within `ADJACENT_DISTANCE_MAX_A` of it.

    When rule 1 fires *and* the pocket also overlaps the partner interface,
    ``also_overlaps_interface`` is set.  Both facts are true and the caller
    must see both: an inter-subunit pocket that is also under the epitope can
    work by either mechanism, and geometry cannot separate them.
    """
    notes: List[str] = []
    pocket = list(pocket_residues)
    if not pocket:
        raise ValueError("pocket_residues is empty")

    # --- lining composition, restricted to the target oligomer -------------
    want = _chain_set(target_chains)
    lining = [r for r in pocket if want is None or r.chain in want]
    if want is not None and len(lining) != len(pocket):
        dropped = sorted({r.chain for r in pocket} - want)
        notes.append(
            f"{len(pocket) - len(lining)} lining residues from chains {dropped} "
            f"outside target_chains were ignored"
        )
    if not lining:
        raise ValueError("no pocket residues fall inside target_chains")

    counts: Dict[str, int] = {}
    for r in lining:
        counts[r.chain] = counts.get(r.chain, 0) + 1
    total = float(len(lining))
    major = sorted(c for c, n in counts.items() if n / total >= multichain_min_share)
    n_subunits = len(major)

    # --- enclosure, and how much of it the other subunits provide ----------
    enc: Optional[float] = None
    enc_single: Optional[float] = None
    gain: Optional[float] = None
    enc_basis = "not_computed"
    if structure is not None:
        pts = probe_points
        if pts is None:
            pts = _lining_centroid(structure, lining, model_index=model_index)
            enc_basis = "lining_centroid"
            notes.append(
                "enclosure computed from the lining centroid, not from probe "
                "points; for a shallow surface pocket this overstates burial"
            )
        else:
            enc_basis = "probe_points"
        chains_for_enc = sorted(want) if want else sorted(counts)
        if want is None:
            notes.append(
                f"target_chains not given; enclosure measured against the "
                f"lining chains only {chains_for_enc}, which underestimates "
                f"burial if other subunits of the oligomer close the pocket"
            )
        enc = enclosure(structure, pts, chains=chains_for_enc, model_index=model_index)
        if len(chains_for_enc) > 1:
            per_chain = [
                enclosure(structure, pts, chains=[c], model_index=model_index)
                for c in chains_for_enc
            ]
            enc_single = max(per_chain)
            gain = enc - enc_single
        else:
            enc_single = enc
            gain = 0.0
    else:
        notes.append(
            "no structure supplied: enclosure and interface distance were not "
            "computed, so the destabiliser test could not run and this result "
            "cannot rule out a buried inter-subunit pocket"
        )

    # --- overlap with the interface ---------------------------------------
    overlap: Optional[float] = None
    coverage: Optional[float] = None
    shared: List[str] = []
    n_iface = 0
    iface = list(interface_residues_) if interface_residues_ else []
    if iface:
        n_iface = len(iface)
        iface_set = _match_keys(iface, match_by)
        shared_keys = [r for r in lining if _match_keys([r], match_by) & iface_set]
        shared = [r.label for r in shared_keys]
        overlap = len(shared_keys) / total
        coverage = len(_match_keys(shared_keys, match_by)) / max(len(iface_set), 1)

    # --- distance to interface --------------------------------------------
    dist: Optional[float] = None
    if structure is not None and iface:
        dist = _min_distance(structure, lining, iface, match_by, model_index=model_index)

    # --- decision ----------------------------------------------------------
    borderline = False
    also_overlaps = False
    buried = (
        enc is not None
        and gain is not None
        and enc >= burial_enclosure_min
        and gain >= subunit_enclosure_gain_min
    )
    if gain is not None and abs(gain - subunit_enclosure_gain_min) <= BORDERLINE_GAIN_BAND:
        borderline = True
        notes.append(
            f"subunit enclosure gain {gain:.2f} is within "
            f"{BORDERLINE_GAIN_BAND} of the {subunit_enclosure_gain_min} "
            f"boundary; the buried-within-oligomer call is undecided here"
        )
    if n_subunits >= 2 and buried:
        label = "destabiliser_candidate"
        if overlap is not None and overlap >= orthosteric_overlap_min:
            also_overlaps = True
            notes.append(
                "pocket is buried between subunits AND overlaps the partner "
                "interface; geometry cannot separate destabilisation from "
                "epitope blockade here"
            )
    elif not iface:
        label = "no_partner_structure"
        if n_subunits >= 2 and not buried:
            notes.append(
                f"pocket is lined by {n_subunits} subunits but does not clear "
                f"the burial test (enclosure "
                f"{'n/a' if enc is None else round(enc, 2)} >= "
                f"{burial_enclosure_min}, gain "
                f"{'n/a' if gain is None else round(gain, 2)} >= "
                f"{subunit_enclosure_gain_min}): an inter-subunit groove at the "
                f"rim, not an internal cavity"
            )
    elif overlap is not None and overlap >= orthosteric_overlap_min:
        label = "orthosteric_candidate"
    else:
        label = "allosteric_candidate"
        if overlap is not None and overlap >= allosteric_overlap_max:
            borderline = True
            notes.append(
                f"overlap {overlap:.2f} sits between {allosteric_overlap_max} "
                f"and {orthosteric_overlap_min}: the pocket touches the "
                f"interface rim, which is neither engagement nor distance"
            )

    adjacent = False
    if (
        label == "allosteric_candidate"
        and not shared  # no residues in common, so "distance" is meaningful
        and dist is not None
        and dist <= ADJACENT_DISTANCE_MAX_A
    ):
        adjacent = True
        notes.append(
            f"pocket shares no residues with the interface but lies {dist:.1f} A "
            f"from it: adjacent, not distal -- a ligand here may act by locking "
            f"a state next to the epitope rather than from a remote site"
        )

    return PocketClassification(
        classification=label,
        overlap_fraction=None if overlap is None else round(overlap, 3),
        interface_coverage=None if coverage is None else round(coverage, 3),
        n_pocket_residues=len(lining),
        n_interface_residues=n_iface,
        n_shared_residues=len(shared),
        shared_residues=shared,
        enclosure=None if enc is None else round(enc, 3),
        enclosure_single_subunit_max=None if enc_single is None else round(enc_single, 3),
        subunit_enclosure_gain=None if gain is None else round(gain, 3),
        enclosure_basis=enc_basis,
        lining_chains=major,
        lining_chain_counts=counts,
        n_lining_subunits=n_subunits,
        min_distance_to_interface_a=None if dist is None else round(dist, 2),
        borderline=borderline,
        also_overlaps_interface=also_overlaps,
        adjacent_to_interface=adjacent,
        match_by=match_by,
        thresholds={
            "orthosteric_overlap_min": orthosteric_overlap_min,
            "allosteric_overlap_max": allosteric_overlap_max,
            "burial_enclosure_min": burial_enclosure_min,
            "subunit_enclosure_gain_min": subunit_enclosure_gain_min,
            "multichain_min_share": multichain_min_share,
            "adjacent_distance_max_a": ADJACENT_DISTANCE_MAX_A,
            "contact_cutoff_a": CONTACT_CUTOFF_A,
        },
        notes=notes,
    )


def _residue_coords(
    structure: gemmi.Structure,
    keys: Sequence[ResidueKey],
    match_by: str = "chain_seqid",
    model_index: int = 0,
) -> np.ndarray:
    wanted = _match_keys(keys, match_by)
    out = []
    for chain in structure[model_index]:
        for res in chain:
            if res.name in _SOLVENT:
                continue
            k = _residue_key(chain, res)
            if next(iter(_match_keys([k], match_by))) in wanted:
                for a in res:
                    if a.element.name not in _LIGHT_ELEMENTS:
                        out.append((a.pos.x, a.pos.y, a.pos.z))
    return np.asarray(out, dtype=float) if out else np.zeros((0, 3))


def _lining_centroid(
    structure: gemmi.Structure, keys: Sequence[ResidueKey], model_index: int = 0
) -> np.ndarray:
    c = _residue_coords(structure, keys, "chain_seqid", model_index=model_index)
    if len(c) == 0:
        raise ValueError("lining residues not found in the supplied structure")
    return c.mean(axis=0).reshape(1, 3)


def _min_distance(
    structure: gemmi.Structure,
    pocket: Sequence[ResidueKey],
    iface: Sequence[ResidueKey],
    match_by: str,
    model_index: int = 0,
) -> Optional[float]:
    a = _residue_coords(structure, pocket, "chain_seqid", model_index=model_index)
    b = _residue_coords(structure, iface, match_by, model_index=model_index)
    if len(a) == 0 or len(b) == 0:
        return None
    best = np.inf
    for start in range(0, len(a), 2048):
        blk = a[start : start + 2048]
        d2 = ((blk[:, None, :] - b[None, :, :]) ** 2).sum(-1)
        best = min(best, float(d2.min()))
    return math.sqrt(best)


# ---------------------------------------------------------------------------
# Cross-structure transfer
# ---------------------------------------------------------------------------


def superpose_chains(
    mobile: gemmi.Structure,
    reference: gemmi.Structure,
    chain_pairs: Sequence[Tuple[str, str]],
    atom_name: str = "CA",
    model_index: int = 0,
) -> Tuple[gemmi.Transform, float, int]:
    """Least-squares superposition of `mobile` onto `reference`.

    Matches residues by sequence number within each chain pair, on `atom_name`.
    Returns (transform, rmsd, n_atoms_used).  Use it to move a ligand from a
    holo entry into the frame of a complex, so that pocket and interface can be
    compared with real chain identities instead of by residue number.
    """
    p_mob, p_ref = [], []
    for mob_chain, ref_chain in chain_pairs:
        mob_map = _ca_map(mobile, mob_chain, atom_name, model_index)
        ref_map = _ca_map(reference, ref_chain, atom_name, model_index)
        for key in sorted(set(mob_map) & set(ref_map)):
            p_mob.append(mob_map[key])
            p_ref.append(ref_map[key])
    if len(p_mob) < 3:
        raise ValueError("fewer than 3 matched atoms; cannot superpose")
    sup = gemmi.superpose_positions(p_ref, p_mob)
    return sup.transform, float(sup.rmsd), len(p_mob)


def _ca_map(st, chain_name, atom_name, model_index):
    out = {}
    for chain in st[model_index]:
        if chain.name != chain_name:
            continue
        for res in chain:
            if not _is_amino_acid(res.name):
                continue
            at = res.find_atom(atom_name, "*")
            if at is not None:
                out[(res.seqid.num, (res.seqid.icode or " ").strip())] = at.pos
    return out


def detect_numbering_offset(
    ref: gemmi.Structure,
    ref_chain: str,
    mobile: gemmi.Structure,
    mobile_chain: str,
    max_offset: int = 250,
    model_index: int = 0,
) -> Tuple[int, float, int]:
    """Constant residue-numbering offset between two entries of one protein.

    Returns ``(offset, identity, n_matched)`` such that
    ``mobile_seqid + offset == ref_seqid``.

    **This is not optional book-keeping.**  PDB entries for the same protein
    routinely use different conventions: IL-17A entries 7UWM, 8DYG and 9SQX
    number from the UniProt precursor (34-155) while 4HSA numbers the mature
    chain (17-131), an offset of 23.  Comparing residue *numbers* across those
    two without correcting gives a silently wrong answer -- a Jaccard of 0.22
    where the real agreement is 0.9.  Correct before any ``match_by="seqid"``
    comparison across entries, or superpose and work in one frame instead.

    Assumes a single constant offset, i.e. no insertions or deletions between
    the two constructs.  `identity` is the fraction of overlapping positions
    whose residue names agree; treat anything below ~0.9 as a failed detection
    and fall back to superposition.
    """
    a = {
        r.seqid.num: r.name
        for c in ref[model_index]
        if c.name == ref_chain
        for r in c
        if _is_amino_acid(r.name)
    }
    b = {
        r.seqid.num: r.name
        for c in mobile[model_index]
        if c.name == mobile_chain
        for r in c
        if _is_amino_acid(r.name)
    }
    if not a or not b:
        raise ValueError("empty chain in numbering-offset detection")
    best = (0, 0.0, 0)
    for off in range(-max_offset, max_offset + 1):
        shared = [n for n in b if (n + off) in a]
        if len(shared) < 10:
            continue
        same = sum(1 for n in shared if b[n] == a[n + off])
        ident = same / len(shared)
        # prefer more matched positions at equal identity
        if (ident, len(shared)) > (best[1], best[2]):
            best = (off, ident, len(shared))
    return best


def renumber_residues(
    keys: Sequence[ResidueKey], offset: int, chain_map: Optional[Dict[str, str]] = None
) -> List[ResidueKey]:
    """Shift residue numbers by `offset` and optionally rename chains."""
    out = []
    for k in keys:
        out.append(
            ResidueKey(
                chain=(chain_map or {}).get(k.chain, k.chain),
                seqid=k.seqid + offset,
                icode=k.icode,
                name=k.name,
            )
        )
    return out


def apply_transform(points: np.ndarray, tr: gemmi.Transform) -> np.ndarray:
    """Apply a gemmi Transform to an (n, 3) coordinate array."""
    out = []
    for p in np.atleast_2d(np.asarray(points, dtype=float)):
        v = tr.apply(gemmi.Position(*p))
        out.append((v.x, v.y, v.z))
    return np.asarray(out, dtype=float)


# ---------------------------------------------------------------------------
# 3. Finding complexes for an accession
# ---------------------------------------------------------------------------

#: The query is verified working against the Paperclip protein database.  It is
#: kept here as the documented route by which a caller finds a partner-bearing
#: structure; the analysis functions above take structures, not accessions.
PARTNER_STRUCTURE_SQL = """
SELECT s.entry_id, s.resolution,
       STRING_AGG(DISTINCT pe.uniprot_accession,' + ') AS partners,
       COUNT(DISTINCT pe.uniprot_accession) AS n_prot
FROM pdb_v.structures_by_accession s
JOIN pdb_v.polymer_entities pe ON pe.entry_id = s.entry_id
WHERE s.accession = '{accession}' AND pe.uniprot_accession IS NOT NULL
GROUP BY s.entry_id, s.resolution
HAVING COUNT(DISTINCT pe.uniprot_accession) > 1
ORDER BY s.resolution NULLS LAST
LIMIT {limit}
""".strip()


def find_partner_structures(
    accession: str,
    limit: int = 8,
    paperclip: str = "paperclip",
    timeout: int = 60,
) -> List[dict]:
    """Deposited structures of `accession` that also contain another protein.

    Shells out to the paperclip CLI (`paperclip sql -s proteins`), which needs
    the repo `.env` sourced for credentials:

        set -a; . <repo>/.env; set +a

    Returns a list of dicts with keys entry_id, resolution, partners, n_prot,
    best resolution first.  Rows include the target's own accession in
    `partners`; the partner is the other accession.  A homo-oligomer will NOT
    appear here -- it is one accession -- which is exactly why the destabiliser
    test in `classify_pocket` does not depend on this query.
    """
    sql = PARTNER_STRUCTURE_SQL.format(accession=accession, limit=int(limit))
    proc = subprocess.run(
        [paperclip, "sql", "-s", "proteins", sql],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"paperclip failed: {proc.stderr.strip()[:400]}")
    return _parse_paperclip_table(proc.stdout)


def _parse_paperclip_table(text: str) -> List[dict]:
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    header = None
    rows = []
    for ln in lines:
        if re.match(r"^-+\+", ln) or re.match(r"^\(\d+ rows?", ln) or ln.startswith("["):
            continue
        parts = [p.strip() for p in ln.split("|")]
        if header is None:
            if len(parts) > 1:
                header = parts
            continue
        if len(parts) != len(header):
            continue
        row = dict(zip(header, parts))
        if row.get("resolution"):
            try:
                row["resolution"] = float(row["resolution"])
            except ValueError:
                row["resolution"] = None
        if row.get("n_prot"):
            try:
                row["n_prot"] = int(row["n_prot"])
            except ValueError:
                pass
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Self-test harness (IL-17A and TNF-alpha fixtures)
# ---------------------------------------------------------------------------


def _fetch(entry_id: str, directory: str) -> str:
    """Biological assembly 1, never the asymmetric unit."""
    path = os.path.join(directory, f"{entry_id.upper()}-assembly1.cif")
    if not os.path.exists(path):
        import urllib.request

        url = f"https://files.rcsb.org/download/{entry_id.upper()}-assembly1.cif"
        os.makedirs(directory, exist_ok=True)
        urllib.request.urlretrieve(url, path)
    return path


def selftest(directory: str) -> int:
    """Run the two documented fixtures and print the numbers."""
    failures = 0

    def check(name, ok):
        nonlocal failures
        print(f"    [{'PASS' if ok else 'FAIL'}] {name}")
        if not ok:
            failures += 1

    print("=" * 72)
    print("IL-17A (Q16552) -- interface vs IL-17RA (Q96F46), two ligand sites")
    print("=" * 72)
    st_7uwm = load_structure(_fetch("7UWM", directory))
    # 7UWM assembly 1 holds two copies of the 2:1 complex (A,B + C and D,E + F).
    # Take one: the IL-17A homodimer A+B against IL-17RA chain C.
    iface = interface_residues(st_7uwm, ["A", "B"], ["C"])
    print(f"  7UWM  IL-17A(A,B) vs IL-17RA(C) @ {iface.cutoff_a} A")
    print(f"    IL-17A side : {len(iface.side_a)} residues")
    print("      " + ", ".join(r.label for r in iface.side_a))
    print(f"    IL-17RA side: {len(iface.side_b)} residues")

    # Independent complex, and a numbering trap: 4HSA numbers the mature chain
    # (17-131) while 7UWM numbers the UniProt precursor (34-155).
    st_4hsa = load_structure(_fetch("4HSA", directory))
    iface_4hsa = interface_residues(st_4hsa, ["A", "B"], ["C"])
    off, ident, n_match = detect_numbering_offset(st_7uwm, "A", st_4hsa, "A")
    a_nums = {r.seqid for r in iface.side_a}
    raw = {r.seqid for r in iface_4hsa.side_a}
    fixed = {r.seqid for r in renumber_residues(iface_4hsa.side_a, off)}
    jac_raw = len(a_nums & raw) / len(a_nums | raw)
    jac = len(a_nums & fixed) / len(a_nums | fixed)
    print(f"  4HSA  IL-17A side: {len(iface_4hsa.side_a)} residues")
    print(f"    numbering offset vs 7UWM: +{off} "
          f"(identity {ident:.2f} over {n_match} positions)")
    print(f"    seqid Jaccard vs 7UWM: {jac_raw:.2f} uncorrected -> {jac:.2f} corrected")
    check("numbering offset 4HSA -> 7UWM detected as +23", off == 23 and ident > 0.9)
    check("7UWM and 4HSA epitopes agree once renumbered (Jaccard > 0.5)", jac > 0.5)

    for entry, comp, label in (
        ("9SQX", "A1JPS", "dimer-interface site (Trp90/Leu120 = Trp67/Leu97 mature)"),
        ("8DYG", "U5Q", "C-terminal site (His109/Met110)"),
    ):
        st = load_structure(_fetch(entry, directory))
        sites = ligand_site_residues(st, comp, lining_chains=["A", "B"])
        print(f"\n  {entry} ligand {comp} -- {label}: {len(sites)} cop(ies)")
        for site in sites:
            res = classify_pocket(
                site.residues,
                iface.side_a,
                st,
                target_chains=["A", "B"],
                probe_points=site.positions,
                match_by="seqid",  # different entry from the complex
            )
            print(f"    copy {site.ligand_key.label}: {res.summary()}")
            print(f"      lining: {', '.join(r.label for r in site.residues)}")
            if res.shared_residues:
                print(f"      shared with epitope: {', '.join(res.shared_residues)}")
            for n in res.notes:
                print(f"      note: {n}")

    print()
    print("=" * 72)
    print("TNF-alpha (P01375) -- SPD304 site must classify as destabiliser")
    print("=" * 72)
    # Receptor epitope from a real complex: 3ALQ is TNF-alpha(A,B,C) + TNFR2
    # (R,S,T) in mature 1-157 numbering, the same numbering as 1TNF and 2AZ5.
    st_3alq = load_structure(_fetch("3ALQ", directory))
    tnfr = interface_residues(st_3alq, ["A", "B", "C"], ["R", "S", "T"])
    print(f"  3ALQ  TNF(A,B,C) vs TNFR2(R,S,T): {len(tnfr.side_a)} TNF residues")
    print("    " + ", ".join(r.label for r in tnfr.side_a))

    # 2AZ5 assembly 1 is a crystallographic tetramer -- two independent TNF
    # dimers.  Restrict to one biological dimer (A,B) for the copy in chain A.
    st_2az5 = load_structure(_fetch("2AZ5", directory))
    sites = ligand_site_residues(
        st_2az5, "307", chains=["A"], lining_chains=["A", "B"]
    )
    site = sites[0]
    res = classify_pocket(
        site.residues,
        tnfr.side_a,
        st_2az5,
        target_chains=["A", "B"],
        probe_points=site.positions,
        match_by="seqid",
    )
    print(f"\n  2AZ5 ligand 307 (holo dimer A+B): {res.summary()}")
    print(f"    lining: {', '.join(r.label for r in site.residues)}")
    check("2AZ5 SPD304 -> destabiliser_candidate",
          res.classification == "destabiliser_candidate")

    # And in the apo trimer, with the ligand transferred by superposition --
    # the state the mechanism is actually about.
    st_1tnf = load_structure(_fetch("1TNF", directory))
    tr, rmsd, n = superpose_chains(st_2az5, st_1tnf, [("A", "A"), ("B", "B")])
    moved = apply_transform(site.positions, tr)
    lining = residues_within(st_1tnf, moved, chains=["A", "B", "C"])
    res_apo = classify_pocket(
        lining,
        tnfr.side_a,
        st_1tnf,
        target_chains=["A", "B", "C"],
        probe_points=moved,
        match_by="seqid",
    )
    print(f"\n  1TNF apo trimer, 307 transferred (rmsd {rmsd:.2f} A over {n} CA):")
    print(f"    {res_apo.summary()}")
    print(f"    lining: {', '.join(r.label for r in lining)}")
    check("1TNF transferred SPD304 site -> destabiliser_candidate",
          res_apo.classification == "destabiliser_candidate")
    check("SPD304 site does not overlap the TNFR epitope (< 0.25)",
          (res_apo.overlap_fraction or 0.0) < ORTHOSTERIC_OVERLAP_MIN)

    print()
    print("=" * 72)
    print("KRAS (P01116) -- controls for the other two labels")
    print("=" * 72)
    st_6oim = load_structure(_fetch("6OIM", directory))
    mov = ligand_site_residues(st_6oim, "MOV", lining_chains=["A"])[0]

    # (a) no complex supplied -> no_partner_structure, and a single-chain
    #     surface pocket must not trip the burial test even though its
    #     absolute enclosure (0.83) exceeds that of the TNF-alpha channel.
    res_np = classify_pocket(
        mov.residues, None, st_6oim,
        target_chains=["A"], probe_points=mov.positions,
    )
    print(f"  6OIM sotorasib site, no partner supplied: {res_np.summary()}")
    check("no interface -> no_partner_structure",
          res_np.classification == "no_partner_structure")
    check("monomeric surface pocket is not called a destabiliser despite "
          f"enclosure {res_np.enclosure}",
          res_np.classification != "destabiliser_candidate")

    # (b) with the effector complex, the switch-II pocket is adjacent to the
    #     RAF1 interface, not in it -- the fourth fixture mechanism.
    st_6vjj = load_structure(_fetch("6VJJ", directory))
    raf = interface_residues(st_6vjj, ["A"], ["B"])
    tr2, rmsd2, n2 = superpose_chains(st_6oim, st_6vjj, [("A", "A")])
    mov_moved = apply_transform(mov.positions, tr2)
    mov_lining = residues_within(st_6vjj, mov_moved, chains=["A"])
    res_kras = classify_pocket(
        mov_lining, raf.side_a, st_6vjj,
        target_chains=["A"], probe_points=mov_moved, match_by="chain_seqid",
    )
    print(f"\n  6VJJ KRAS(A) vs RAF1 RBD(B): {len(raf.side_a)} KRAS residues")
    print(f"    " + ", ".join(r.label for r in raf.side_a))
    print(f"  sotorasib site transferred (rmsd {rmsd2:.2f} A over {n2} CA):")
    print(f"    {res_kras.summary()}")
    check("KRAS switch-II is adjacent to the effector interface, not in it "
          "-> not orthosteric",
          res_kras.classification == "allosteric_candidate")
    check("KRAS switch-II flagged adjacent_to_interface (overlap 0 but "
          f"{res_kras.min_distance_to_interface_a} A away)",
          res_kras.adjacent_to_interface)

    print()
    print(f"{'ALL CHECKS PASSED' if failures == 0 else str(failures) + ' CHECK(S) FAILED'}")
    return failures


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--selftest", metavar="DIR", help="run fixture harness, caching CIFs in DIR")
    ap.add_argument("--partners", metavar="ACCESSION", help="list complex structures for an accession")
    args = ap.parse_args(argv)
    if args.partners:
        for row in find_partner_structures(args.partners):
            print(row)
        return 0
    if args.selftest:
        return selftest(args.selftest)
    ap.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
