"""Transferred-homolog site anchoring, with three machine-checked guards.

WHY THIS FILE EXISTS
--------------------
An anchor-agreement test put sixteen **ligand-free** site anchors against the
known ligand site on four targets (TNF-alpha, IL-17A, NLRP3, S1PR1). Four could
not be built at all. Of the twelve that could, **four found the site**:

    transferred homolog      2 of 3 constructible   <- the only positive record
    interface                1 of 6
    symmetry axis            1 of 3
    annotated function       0 of 4

Four targets is not a rate. What it is enough to say is that transferred
homolog is the one worth building on, and that **none of the four is safe
unaided**. This module is the "building on it" half. The other three anchors
are handled by reporting, not selection, and that lives in `pocket-scan`
(`_annotate_pocket_labels`) — see `SKILL.md`, section "The other three anchors".

`pocket-scan/modal_app.py::_annotate_pocket_labels` documents
`transferred_homolog_site` as "NOT AVAILABLE HERE — it needs Foldseek, which
lives in structure-select/neighbour_precedent". This is that missing piece. It
emits the label plus the provenance block that must travel with it.

THE THREE GUARDS, AND THE MEASURED FAILURE EACH ONE CLOSES
----------------------------------------------------------
Every one of these comes from a transfer that ran, produced a confident answer,
and was wrong. None is a hypothetical.

**Guard 1 — the donor's ligand must pass the drug-like classifier.**
NLRP3's only constructible donor was NOD2 (5IRN) with **ADP**, a cofactor. The
transfer itself was excellent — validated to 0.68 A against NLRP3's own ADP —
and it still selected the nucleotide lobe rather than the drug site. *A perfect
transfer of the wrong ligand is a wrong answer carrying high confidence*, which
is worse than a failed transfer, because a failed transfer announces itself.
The check is `ligand_filter.classify_ligand`, which this skill owns, and the
rejected verdicts are every one that is not `druglike`.

**Guard 2 — domain attribution: the donor ligand's contact shell must overlap
the ALIGNED REGION, not merely sit on the same chain.**
Measured on 7KRZ: bortezomib (`BO2`) is on the correct LONP1 chain, at auth
768-898, while the NACHT-aligned region is auth 506-721. Right chain, wrong
domain. And the proteasome hits never carry their entry-title drug on the AAA+
chain that produced the fold match at all. **Chain-level attribution — which we
already do, and which `SKILL.md` spends two failure modes on — does not catch
either case.** A ligand 250 residues outside the aligned span is transferred by
a rotation fitted somewhere else entirely; the transfer is an extrapolation and
its coordinates mean nothing.

**Guard 3 — a declared TM-score floor, an RMSD ceiling, and a steric check.**
IL-2 forced onto IL-17A at **TM 0.254-0.274** put the transferred ligand
**21.59 A away and inside the protein** — 29 atoms within 2.0 A of protein
heavy atoms, closest contact **0.51 A**. That is physically impossible and it
was produced silently. The floor is discussed at `MIN_TM_SCORE` below; the
short version is that **the steric check is the stronger guard than any RMSD
number**, because it tests the thing we actually care about (did the ligand
land somewhere a ligand can be) rather than a proxy for it.

But the steric check has to be counted against **backbone**, and that was
forced by the reproductions rather than assumed — see
`MAX_BACKBONE_CLASH_ATOMS`. On a total-heavy-atom count the S1PR1 success is
*worse* than the IL-2 failure (11 clashing atoms of 38 at 0.16 A, against 12
of 33 at 0.65 A) while landing 1.79 A from the reference ligand rather than
21.41 A away. Transferring any ligand between two crystal structures produces
side-chain overlaps, because the two structures disagree about rotamers. On
backbone the classes separate completely: 0, 0, 0 for the three good transfers
and 7 for the failure.

NOTHING MAY FIRE SILENTLY
-------------------------
Every result carries `donor_pdb_id`, `donor_ligand`, `tm_score`, `rmsd_a`,
`aligned_length` and `clash_count`, on rejection as well as acceptance. A
transfer without its provenance is not usable evidence and
`TransferResult.provenance` is not optional. `assert_reportable()` raises if any
of the six is missing.

WHAT THIS MODULE DOES NOT DO
----------------------------
It does not score pockets, it does not choose among pockets, and it does not
decide tractability. It produces one anchor — a set of coordinates and a set of
target residues — with the evidence needed to discount it.

PURE STDLIB, ON PURPOSE
-----------------------
The deployed sandbox has no `gemmi`, no `numpy`, no `fpocket` and no
`paperclip` binary (see `managed/druggability-dossier/CLAUDE.md`, "Your tools,
and what the sandbox can and cannot do"). It does have outbound network. So the
mmCIF parser, the sequence alignment, the superposition and the TM-score are
all implemented here in stdlib, and chem-comp records are fetched from RCSB
rather than Paperclip. That is not a preference: on 2026-08-15 `paperclip sql`
against `pdb_v.chemcomps` returned `[error] Request timed out` after 122 s
**with exit code 0**, from a config carrying `{"cli_cwd": "/"}`.
"""

from __future__ import annotations

import json
import math
import os
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

__all__ = [
    "MIN_TM_SCORE",
    "MAX_ALIGNMENT_RMSD_A",
    "MARGINAL_ALIGNMENT_RMSD_A",
    "CLASH_DISTANCE_A",
    "MAX_CLASH_ATOMS",
    "MIN_SHELL_IN_ALIGNED_FRACTION",
    "REJECTED_LIGAND_VERDICTS",
    "Structure",
    "TransferResult",
    "GuardResult",
    "transfer_homolog_site",
    "load_structure",
    "fetch_mmcif",
    "RcsbChemComps",
]


# --------------------------------------------------------------------------
# Guard thresholds. Every one of these is DECLARED HERE and nowhere else.
# --------------------------------------------------------------------------

#: Verdicts from `ligand_filter` that disqualify a donor ligand (guard 1).
#: This is every verdict the classifier can return except `druglike`, plus
#: `unknown`. `unknown` is rejected because a donor we cannot classify is a
#: donor we cannot vouch for, and the whole point of guard 1 is that the
#: transfer machinery is *more* confident than the evidence supports.
REJECTED_LIGAND_VERDICTS = frozenset({
    "cofactor",
    "ion_or_solvent",
    "lipid_or_detergent",
    "sugar_or_glycan",
    "crystallisation_additive",
    "peptide_or_polymer",
    "polymer_conjugate",
    "unknown",
})

#: Guard 2. Fraction of the donor ligand's contact-shell residues that must lie
#: inside the aligned region.
#:
#: DECLARED, NOT FITTED. The rule being enforced is "the ligand sits in the part
#: of the donor that we actually superposed", and a *majority* is the weakest
#: reading of that sentence that still means anything. The measured failures are
#: not near this line — 7KRZ's bortezomib shell (auth 768-898) and the NACHT
#: aligned region (auth 506-721) are disjoint, so it scores 0.0, and the
#: proteasome hits score 0.0 by carrying no ligand on the matched chain at all.
#: A threshold anywhere in (0.0, 1.0) would reject both. 0.5 is chosen because
#: it is the one value in that range that can be justified without reference to
#: the four cases.
MIN_SHELL_IN_ALIGNED_FRACTION = 0.5

#: Guard 3, part 1. TM-score floor.
#:
#: **This is the literature fold-identity threshold, deliberately NOT a number
#: read off our four cases.** TM-score >= 0.5 is the standard cut at which two
#: structures are held to share a fold (Xu & Zhang 2010: below 0.5 the pair is
#: statistically indistinguishable from a random pair). We adopt it because our
#: own data cannot support a tuned value and would flatter itself if we tried:
#: n=2 successes is not a distribution. What our data does do is confirm the
#: standard cut discriminates the cases we have — IL-2 onto IL-17A sits at
#: 0.254-0.274, well under it, and both successes sit well over it.
MIN_TM_SCORE = 0.5

#: Guard 3, part 2. Hard RMSD ceiling over the aligned CA core, in angstrom.
#:
#: **Deliberately loose, and the reason is in the data.** TNF's donor fit at
#: 1.35 A over 111 CA pairs; S1PR1's fit at a mediocre 3.04 A over 115 and
#: still landed. But S1PR1's donor (LPA1) is a ~45%-identical paralogue, not a
#: distant fold neighbour, so **3.04 A must not be read as a safe general
#: floor** — it is one paralogue's tolerance, and a distant neighbour at 3.0 A
#: is a different and much weaker claim. We therefore refuse to set the ceiling
#: anywhere near the observed values. 5.0 A rejects only alignments that are
#: incoherent by any standard, and the real work is done by the steric check.
MAX_ALIGNMENT_RMSD_A = 5.0

#: Above this the transfer is accepted but flagged `alignment_marginal` in the
#: provenance. S1PR1's 3.04 A is a marginal fit that happened to land; a reader
#: is entitled to know that without us pretending it was a good one.
MARGINAL_ALIGNMENT_RMSD_A = 2.0

#: Guard 3, part 3. Heavy-atom separation below which a transferred ligand atom
#: is counted as clashing with a target protein heavy atom.
#:
#: 2.0 A is far inside any heavy-atom van der Waals contact (a C--C contact
#: bottoms out near 3.4 A, and even a strong hydrogen bond holds its heavy
#: atoms near 2.7 A). An atom pair at 2.0 A is not a tight fit; it is the same
#: piece of space claimed twice. The IL-2/IL-17A failure produced 29 such atoms
#: with a closest approach of 0.51 A.
CLASH_DISTANCE_A = 2.0

#: Backbone atom names. A clash against one of these is what voids a transfer.
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O", "OXT"})

#: How many transferred ligand atoms may clash with protein BACKBONE. Zero.
#:
#: **This is the guard that is stronger than any RMSD number, and it is the one
#: allowed to void a transfer that every other check passed.** But it must be
#: counted against backbone, not against all heavy atoms, and that distinction
#: was forced by the measurements — a total-atom count does not separate the
#: cases at all:
#:
#:     case                    backbone      side-chain     wanted
#:     TNF   2AZ5 <- 3LKJ  LKJ    0 (2.40 A)    4 (1.15 A)   accept
#:     S1PR1 3V2Y <- 4Z34  ON7    0 (2.50 A)   11 (0.16 A)   accept
#:     NLRP3 7ALV <- 5IRN  ADP    0 (2.59 A)    0 (2.99 A)   guard 1
#:     IL17A 9SQX <- 1M48  FRG    7 (0.65 A)    7 (0.80 A)   REJECT
#:
#: On totals, S1PR1 (11 of 38 atoms, closest 0.16 A) is *worse* than the IL-2
#: failure (12 of 33, closest 0.65 A) — and S1PR1 lands 1.79 A from the
#: reference ligand while IL-2 lands 21.41 A away. A total-atom threshold
#: would therefore have to reject the second success to reject the failure.
#: On backbone the two classes separate completely, 0/0/0 against 7.
#:
#: The reason is physical and is already in `CLAUDE.md` rule 5: side-chain
#: occlusion and backbone occlusion are different mechanisms with different
#: prognoses. A ligand overlapping a side chain is a **rotamer** problem —
#: two crystal structures of the same pocket disagree on rotamers, so a
#: transferred ligand almost always overlaps some — and it is resolvable by
#: repacking. A ligand overlapping N, CA, C or O is not resolvable by
#: anything: no rotamer moves a backbone. So the transfer is void on backbone
#: contact and merely *flagged* on side-chain contact.
MAX_BACKBONE_CLASH_ATOMS = 0

#: Retained name for the total-atom count, which is REPORTED and never decides.
#: The measured IL-2/IL-17A failure is on record as "29 atoms within 2.0 A of
#: protein heavy atoms, closest contact 0.51 A", which is this quantity, so it
#: stays reportable and comparable.
MAX_CLASH_ATOMS = MAX_BACKBONE_CLASH_ATOMS

#: Donor protein residues with any heavy atom within this distance of a donor
#: ligand heavy atom form that ligand's contact shell (guard 2). 4.5 A is the
#: conventional first-shell contact radius and is what `pocket-scan` uses for
#: lining residues, so the two are comparable.
CONTACT_SHELL_A = 4.5

#: A CA pair closer than this after the TM-optimal superposition counts toward
#: the reported `aligned_length` and `rmsd_a`. This is the "core" convention:
#: RMSD over every seeded pair including the ones that did not superpose is not
#: a meaningful number.
CORE_PAIR_A = 5.0

#: Polymer chains shorter than this are ignored when building the alignment.
MIN_CHAIN_RESIDUES = 25


# --------------------------------------------------------------------------
# mmCIF: a fast dedicated `_atom_site` reader
# --------------------------------------------------------------------------

_AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    # The modified residues common enough to matter for an alignment seed.
    "MSE": "M", "SEP": "S", "TPO": "T", "PTR": "Y", "CSO": "C", "HYP": "P",
    "MLY": "K", "M3L": "K", "KCX": "K", "LLP": "K", "CME": "C", "PCA": "E",
    "SEC": "C", "PYL": "K", "FME": "M", "ABA": "A", "AIB": "A", "NLE": "L",
    "ORN": "K", "DAL": "A", "DLE": "L", "DVA": "V", "DPR": "P", "SAR": "G",
}

_WATERS = frozenset({"HOH", "DOD", "WAT", "H2O"})


@dataclass(frozen=True)
class Residue:
    """One residue or one het group. `atoms` is `{atom_id: (x, y, z)}`."""

    chain: str
    seq: str          # auth_seq_id, kept as a string: insertion codes exist
    comp_id: str
    is_polymer: bool
    atoms: dict[str, tuple[float, float, float]]

    @property
    def key(self) -> str:
        return f"{self.chain}/{self.seq}"

    @property
    def one_letter(self) -> str:
        return _AA3.get(self.comp_id.upper(), "X")

    def ca(self) -> tuple[float, float, float] | None:
        return self.atoms.get("CA")

    def heavy(self) -> list[tuple[float, float, float]]:
        return list(self.atoms.values())

    def seq_int(self) -> int | None:
        s = "".join(c for c in self.seq if c.isdigit() or c == "-")
        try:
            return int(s)
        except ValueError:
            return None


@dataclass
class Structure:
    """The parts of an mmCIF this module needs, and nothing else."""

    pdb_id: str
    residues: list[Residue]

    # ---- polymer views -------------------------------------------------
    def polymer_chains(self) -> dict[str, list[Residue]]:
        out: dict[str, list[Residue]] = {}
        for r in self.residues:
            if r.is_polymer and r.ca() is not None and r.one_letter != "X":
                out.setdefault(r.chain, []).append(r)
        return {
            c: rs for c, rs in out.items() if len(rs) >= MIN_CHAIN_RESIDUES
        }

    def protein_heavy_atoms(
        self, chains: Iterable[str] | None = None
    ) -> list[tuple[float, float, float]]:
        keep = set(chains) if chains is not None else None
        out: list[tuple[float, float, float]] = []
        for r in self.residues:
            if not r.is_polymer:
                continue
            if keep is not None and r.chain not in keep:
                continue
            out.extend(r.heavy())
        return out

    # ---- het views -----------------------------------------------------
    def het_instances(self, comp_id: str) -> list[Residue]:
        cid = comp_id.upper()
        return [
            r for r in self.residues
            if not r.is_polymer and r.comp_id.upper() == cid
        ]

    def het_comp_ids(self) -> list[str]:
        seen: dict[str, int] = {}
        for r in self.residues:
            if r.is_polymer or r.comp_id.upper() in _WATERS:
                continue
            seen[r.comp_id.upper()] = seen.get(r.comp_id.upper(), 0) + 1
        return sorted(seen, key=lambda c: -seen[c])


_ATOM_FIELDS = (
    "group_PDB", "label_atom_id", "label_alt_id", "label_comp_id",
    "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
    "Cartn_x", "Cartn_y", "Cartn_z", "type_symbol", "pdbx_PDB_model_num",
    "label_asym_id", "label_seq_id",
)


def parse_atom_site(text: str, pdb_id: str = "") -> Structure:
    """Read `_atom_site` out of an mmCIF document.

    A dedicated reader rather than `ligand_filter.read_mmcif_categories`
    because `_atom_site` is the one category that is always a plain `loop_`
    with no multi-line values, and it is also the one with 10^4-10^6 rows. The
    general tokenizer is correct on it and roughly an order of magnitude
    slower.

    Model 1 only. Altloc: blank/`.` or the first letter seen for that atom.
    Hydrogens and deuteriums are dropped — every distance in this module is a
    heavy-atom distance, and a file that happens to carry riding hydrogens must
    not produce different clash counts from one that does not.
    """
    lines = text.splitlines()
    n = len(lines)
    i = 0
    header: list[str] = []
    start = -1
    while i < n:
        if lines[i].startswith("loop_"):
            j = i + 1
            tags: list[str] = []
            while j < n and lines[j].lstrip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            if tags and tags[0].startswith("_atom_site."):
                header = [t.split(".", 1)[1] for t in tags]
                start = j
                break
            i = j
            continue
        i += 1
    if start < 0:
        return Structure(pdb_id=pdb_id, residues=[])

    idx = {f: header.index(f) if f in header else -1 for f in _ATOM_FIELDS}
    need = ("group_PDB", "label_comp_id", "auth_asym_id", "auth_seq_id",
            "Cartn_x", "Cartn_y", "Cartn_z", "label_atom_id")
    for f in need:
        if idx[f] < 0:
            raise ValueError(f"_atom_site is missing {f}")

    i_model = idx["pdbx_PDB_model_num"]
    i_alt = idx["label_alt_id"]
    i_elem = idx["type_symbol"]
    i_ins = idx["pdbx_PDB_ins_code"]

    # residue key -> Residue-under-construction
    order: list[tuple[str, str, str, bool]] = []
    build: dict[tuple[str, str, str, bool], dict[str, tuple[float, float, float]]] = {}
    alt_seen: dict[tuple[str, str, str, bool], str] = {}

    k = start
    while k < n:
        raw = lines[k]
        k += 1
        if not raw or raw[0] == "#":
            break
        if raw.startswith("loop_") or raw.lstrip().startswith("_"):
            break
        parts = raw.split()
        if len(parts) < len(header):
            # Quoted atom names such as "O5'" contain no space, so a short row
            # here is a genuinely malformed line; skip rather than guess.
            continue
        if i_model >= 0 and parts[i_model] not in ("1", ".", "?"):
            continue
        if i_elem >= 0 and parts[i_elem].upper() in ("H", "D"):
            continue
        grp = parts[idx["group_PDB"]]
        comp = parts[idx["label_comp_id"]].strip("'\"")
        if comp.upper() in _WATERS:
            continue
        chain = parts[idx["auth_asym_id"]].strip("'\"")
        seq = parts[idx["auth_seq_id"]]
        if i_ins >= 0 and parts[i_ins] not in (".", "?"):
            seq = seq + parts[i_ins]
        is_poly = grp == "ATOM"
        rkey = (chain, seq, comp, is_poly)
        alt = parts[i_alt] if i_alt >= 0 else "."
        if alt not in (".", "?", ""):
            first = alt_seen.setdefault(rkey, alt)
            if alt != first:
                continue
        try:
            xyz = (
                float(parts[idx["Cartn_x"]]),
                float(parts[idx["Cartn_y"]]),
                float(parts[idx["Cartn_z"]]),
            )
        except ValueError:
            continue
        name = parts[idx["label_atom_id"]].strip("'\"")
        d = build.get(rkey)
        if d is None:
            d = build[rkey] = {}
            order.append(rkey)
        d.setdefault(name, xyz)

    residues = [
        Residue(chain=c, seq=s, comp_id=comp, is_polymer=p, atoms=build[key])
        for key in order
        for (c, s, comp, p) in (key,)
    ]
    return Structure(pdb_id=pdb_id, residues=residues)


# --------------------------------------------------------------------------
# Fetching
# --------------------------------------------------------------------------

def _cache_dir() -> Path:
    """Where downloaded mmCIFs live.

    NOT inside the skill directory, deliberately. `scripts/deploy.ts` zips each
    `.claude/skills/<dir>/` whole with no exclusions, so a structure cache
    parked beside this file would ship to the Skills API on every deploy — the
    same trap the ligand_filter test fixtures had to be trimmed to escape (see
    `tests/README.md`). Entry files here run to megabytes.
    """
    p = os.environ.get("STRUCTURE_SELECT_CACHE")
    if p:
        return Path(p)
    return Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) \
        / "structure-select" / "mmcif"


def fetch_mmcif(pdb_id: str, *, timeout: float = 120.0) -> str:
    """mmCIF text for a PDB ID, cached on disk."""
    pid = pdb_id.strip().upper()
    cache = _cache_dir()
    f = cache / f"{pid}.cif"
    if f.is_file():
        return f.read_text()
    url = f"https://files.rcsb.org/download/{pid}.cif"
    req = urllib.request.Request(  # noqa: S310
        url, headers={"User-Agent": "structure-select/homolog_transfer"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as fh:  # noqa: S310
        text = fh.read().decode("utf-8", "replace")
    try:
        cache.mkdir(parents=True, exist_ok=True)
        f.write_text(text)
    except OSError:
        pass
    return text


def load_structure(pdb_id: str, *, path: str | os.PathLike[str] | None = None
                   ) -> Structure:
    """Parse a local mmCIF if given one, else fetch by ID."""
    if path is not None:
        return parse_atom_site(Path(path).read_text(), pdb_id.upper())
    return parse_atom_site(fetch_mmcif(pdb_id), pdb_id.upper())


class RcsbChemComps:
    """Chem-comp records from RCSB's REST API, shaped for `ligand_filter`.

    `ligand_filter.ChemCompSource` reads Paperclip and Paperclip only, which is
    correct for the pipeline and unusable here for two reasons measured on
    2026-08-15: the deployed sandbox has no `paperclip` binary at all (calls go
    through the `paperclip_sql` tool), and on the operator's machine the shared
    config was carrying `{"cli_cwd": "/"}`, under which
    `paperclip sql -s proteins "SELECT comp_id, type FROM pdb_v.chemcomps
    WHERE comp_id = 'ADP'"` returned `[error] Request timed out` after 122 s
    **with exit code 0**.

    Guard 1 must not be silently skippable, so it gets a source that has no
    credentials and no shared mutable state. The verdicts still come from
    `ligand_filter.classify_record` — only the row does not.
    """

    _URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"

    def __init__(self, *, timeout: float = 30.0,
                 cache_path: str | os.PathLike[str] | None = None) -> None:
        self._timeout = timeout
        self._cache: dict[str, dict[str, Any] | None] = {}
        self._cache_path = Path(cache_path) if cache_path else None
        if self._cache_path and self._cache_path.is_file():
            try:
                self._cache.update(json.loads(self._cache_path.read_text()))
            except (OSError, ValueError):
                pass
        #: comp_id -> why the lookup failed. A LOOKUP FAILURE IS NOT A MISS —
        #: same contract as `ligand_filter.ChemCompSource.fetch_errors`, and
        #: for the same reason: an unreachable CCD row must not render as
        #: "not drug-like", because that would make guard 1 fail *open* on
        #: exactly the runs where the network is the problem.
        self.fetch_errors: dict[str, str] = {}

    def preload(self, records: dict[str, dict[str, Any]]) -> None:
        for k, v in records.items():
            self._cache[k.upper()] = dict(v)

    def get(self, comp_id: str) -> dict[str, Any] | None:
        cid = comp_id.upper()
        if cid in self._cache:
            return self._cache[cid]
        try:
            req = urllib.request.Request(  # noqa: S310
                self._URL.format(cid),
                headers={"User-Agent": "structure-select/homolog_transfer"},
            )
            with urllib.request.urlopen(req, timeout=self._timeout) as fh:  # noqa: S310
                d = json.loads(fh.read().decode())
        except (urllib.error.URLError, OSError, ValueError) as exc:
            self.fetch_errors[cid] = f"{type(exc).__name__}: {exc}"
            return None
        cc = d.get("chem_comp") or {}
        desc = d.get("rcsb_chem_comp_descriptor") or {}
        rel = d.get("rcsb_chem_comp_related") or []
        drugbank = next(
            (r.get("resource_accession_code") for r in rel
             if (r.get("resource_name") or "").lower() == "drugbank"), None,
        )
        rec = {
            "comp_id": cc.get("id") or cid,
            "type": cc.get("type"),
            "formula": cc.get("formula"),
            "formula_weight": cc.get("formula_weight"),
            "drugbank_id": drugbank,
            "inchikey": desc.get("InChIKey"),
            "smiles": desc.get("SMILES_stereo") or desc.get("SMILES"),
            "name": cc.get("name"),
        }
        self._cache[cid] = rec
        return rec

    def save_cache(self) -> None:
        if self._cache_path:
            self._cache_path.parent.mkdir(parents=True, exist_ok=True)
            self._cache_path.write_text(json.dumps(self._cache, sort_keys=True))


# --------------------------------------------------------------------------
# Sequence alignment (Gotoh, affine gaps, BLOSUM62)
# --------------------------------------------------------------------------

_B62_ORDER = "ARNDCQEGHILKMFPSTWYVBZX*"
_B62_ROWS = """
  4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0  -2  -1   0  -4
 -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3  -1   0  -1  -4
 -2   0   6   1  -3   0   0   0   1  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3   3   0  -1  -4
 -2  -2   1   6  -3   0   2  -1  -1  -3  -4  -1  -3  -3  -1   0  -1  -4  -3  -3   4   1  -1  -4
  0  -3  -3  -3   9  -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1  -3  -3  -2  -4
 -1   1   0   0  -3   5   2  -2   0  -3  -2   1   0  -3  -1   0  -1  -2  -1  -2   0   3  -1  -4
 -1   0   0   2  -4   2   5  -2   0  -3  -3   1  -2  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
  0  -2   0  -1  -3  -2  -2   6  -2  -4  -4  -2  -3  -3  -2   0  -2  -2  -3  -3  -1  -2  -1  -4
 -2   0   1  -1  -3   0   0  -2   8  -3  -3  -1  -2  -1  -2  -1  -2  -2   2  -3   0   0  -1  -4
 -1  -3  -3  -3  -1  -3  -3  -4  -3   4   2  -3   1   0  -3  -2  -1  -3  -1   3  -3  -3  -1  -4
 -1  -2  -3  -4  -1  -2  -3  -4  -3   2   4  -2   2   0  -3  -2  -1  -2  -1   1  -4  -3  -1  -4
 -1   2   0  -1  -3   1   1  -2  -1  -3  -2   5  -1  -3  -1   0  -1  -3  -2  -2   0   1  -1  -4
 -1  -1  -2  -3  -1   0  -2  -3  -2   1   2  -1   5   0  -2  -1  -1  -1  -1   1  -3  -1  -1  -4
 -2  -3  -3  -3  -2  -3  -3  -3  -1   0   0  -3   0   6  -4  -2  -2   1   3  -1  -3  -3  -1  -4
 -1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4   7  -1  -1  -4  -3  -2  -2  -1  -2  -4
  1  -1   1   0  -1   0   0   0  -1  -2  -2   0  -1  -2  -1   4   1  -3  -2  -2   0   0   0  -4
  0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1   1   5  -2  -2   0  -1  -1   0  -4
 -3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1   1  -4  -3  -2  11   2  -3  -4  -3  -2  -4
 -2  -2  -2  -3  -2  -1  -2  -3   2  -1  -1  -2  -1   3  -3  -2  -2   2   7  -1  -3  -2  -1  -4
  0  -3  -3  -3  -1  -2  -2  -3  -3   3   1  -2   1  -1  -2  -2   0  -3  -1   4  -3  -2  -1  -4
 -2  -1   3   4  -3   0   1  -1   0  -3  -4   0  -3  -3  -2   0  -1  -4  -3  -3   4   1  -1  -4
 -1   0   0   1  -3   3   4  -2   0  -3  -3   1  -1  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
  0  -1  -1  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2   0   0  -2  -1  -1  -1  -1  -1  -4
 -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4   1
"""


def _build_blosum62() -> dict[tuple[str, str], int]:
    m: dict[tuple[str, str], int] = {}
    rows = [r for r in _B62_ROWS.strip("\n").split("\n") if r.strip()]
    for a, row in zip(_B62_ORDER, rows):
        vals = [int(x) for x in row.split()]
        if len(vals) != len(_B62_ORDER):
            raise ValueError(f"BLOSUM62 row {a} has {len(vals)} columns")
        for b, v in zip(_B62_ORDER, vals):
            m[(a, b)] = v
    return m


_B62 = _build_blosum62()
_GAP_OPEN = -11
_GAP_EXTEND = -1


def _sub(a: str, b: str) -> int:
    return _B62.get((a, b), _B62.get(("X", "X"), -1))


def align_sequences(s1: str, s2: str) -> list[tuple[int, int]]:
    """Global affine-gap alignment. Returns aligned index pairs `(i, j)`.

    Gotoh, three matrices, traceback via direction codes. Used only to *seed*
    the structural superposition — the TM-score refinement below is free to
    throw any of these pairs away, which is what makes a bad seed survivable.
    """
    n, m = len(s1), len(s2)
    if n == 0 or m == 0:
        return []
    neg = float("-inf")
    # M: ends aligned; X: gap in s2; Y: gap in s1
    prev_m = [neg] * (m + 1)
    prev_x = [neg] * (m + 1)
    prev_y = [neg] * (m + 1)
    prev_m[0] = 0.0
    for j in range(1, m + 1):
        prev_y[j] = _GAP_OPEN + _GAP_EXTEND * (j - 1)
    # traceback: 0 = M, 1 = X, 2 = Y
    tb_m = [[0] * (m + 1) for _ in range(n + 1)]
    tb_x = [[0] * (m + 1) for _ in range(n + 1)]
    tb_y = [[0] * (m + 1) for _ in range(n + 1)]
    for j in range(1, m + 1):
        tb_y[0][j] = 2

    for i in range(1, n + 1):
        cur_m = [neg] * (m + 1)
        cur_x = [neg] * (m + 1)
        cur_y = [neg] * (m + 1)
        cur_x[0] = _GAP_OPEN + _GAP_EXTEND * (i - 1)
        tb_x[i][0] = 1
        a = s1[i - 1]
        for j in range(1, m + 1):
            sc = _sub(a, s2[j - 1])
            best, src = prev_m[j - 1], 0
            if prev_x[j - 1] > best:
                best, src = prev_x[j - 1], 1
            if prev_y[j - 1] > best:
                best, src = prev_y[j - 1], 2
            cur_m[j] = best + sc
            tb_m[i][j] = src

            o, e = prev_m[j] + _GAP_OPEN, prev_x[j] + _GAP_EXTEND
            if o >= e:
                cur_x[j], tb_x[i][j] = o, 0
            else:
                cur_x[j], tb_x[i][j] = e, 1

            o, e = cur_m[j - 1] + _GAP_OPEN, cur_y[j - 1] + _GAP_EXTEND
            if o >= e:
                cur_y[j], tb_y[i][j] = o, 0
            else:
                cur_y[j], tb_y[i][j] = e, 2
        prev_m, prev_x, prev_y = cur_m, cur_x, cur_y

    i, j = n, m
    state = max(((prev_m[m], 0), (prev_x[m], 1), (prev_y[m], 2)))[1]
    pairs: list[tuple[int, int]] = []
    while i > 0 or j > 0:
        if state == 0:
            if i == 0 or j == 0:
                state = 1 if j == 0 else 2
                continue
            pairs.append((i - 1, j - 1))
            state = tb_m[i][j]
            i, j = i - 1, j - 1
        elif state == 1:
            if i == 0:
                state = 2
                continue
            state = tb_x[i][j]
            i -= 1
        else:
            if j == 0:
                state = 1
                continue
            state = tb_y[i][j]
            j -= 1
    pairs.reverse()
    return pairs


# --------------------------------------------------------------------------
# Superposition: Horn/Kearsley quaternion, Jacobi eigensolver
# --------------------------------------------------------------------------

def _jacobi(a: list[list[float]], iters: int = 100
            ) -> tuple[list[float], list[list[float]]]:
    """Eigenvalues and eigenvectors of a real symmetric matrix.

    Cyclic Jacobi. Chosen over an SVD because it is 40 lines of stdlib, it is
    unconditionally stable on the 4x4 symmetric key matrix below, and the
    quaternion route it enables **cannot return an improper rotation** — the
    reflection trap that a hand-rolled SVD Kabsch falls into and that would
    silently mirror a ligand into the wrong enantiomeric site.
    """
    n = len(a)
    m = [row[:] for row in a]
    v = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(iters):
        off = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                off += m[i][j] * m[i][j]
        if off < 1e-24:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(m[p][q]) < 1e-18:
                    continue
                theta = (m[q][q] - m[p][p]) / (2.0 * m[p][q])
                t = (1.0 if theta >= 0 else -1.0) / (
                    abs(theta) + math.sqrt(theta * theta + 1.0)
                )
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    mkp, mkq = m[k][p], m[k][q]
                    m[k][p] = c * mkp - s * mkq
                    m[k][q] = s * mkp + c * mkq
                for k in range(n):
                    mpk, mqk = m[p][k], m[q][k]
                    m[p][k] = c * mpk - s * mqk
                    m[q][k] = s * mpk + c * mqk
                for k in range(n):
                    vkp, vkq = v[k][p], v[k][q]
                    v[k][p] = c * vkp - s * vkq
                    v[k][q] = s * vkp + c * vkq
    return [m[i][i] for i in range(n)], v


Vec3 = tuple[float, float, float]
Rot = tuple[tuple[float, float, float], tuple[float, float, float],
            tuple[float, float, float]]


def superpose(mob: Sequence[Vec3], ref: Sequence[Vec3]
              ) -> tuple[Rot, Vec3, float]:
    """Rotation+translation carrying `mob` onto `ref`, and the RMSD.

    Apply as `x' = R @ (x - centroid_mob) + centroid_ref`, which is what
    `apply_transform` below does. Returns `(R, t, rmsd)` with `t` already
    folded so that `x' = R @ x + t`.
    """
    n = len(mob)
    if n < 3 or n != len(ref):
        raise ValueError(f"need >= 3 matched points, got {n} vs {len(ref)}")
    cm = [sum(p[k] for p in mob) / n for k in range(3)]
    cr = [sum(p[k] for p in ref) / n for k in range(3)]
    sxx = sxy = sxz = syx = syy = syz = szx = szy = szz = 0.0
    for p, q in zip(mob, ref):
        px, py, pz = p[0] - cm[0], p[1] - cm[1], p[2] - cm[2]
        qx, qy, qz = q[0] - cr[0], q[1] - cr[1], q[2] - cr[2]
        sxx += px * qx; sxy += px * qy; sxz += px * qz
        syx += py * qx; syy += py * qy; syz += py * qz
        szx += pz * qx; szy += pz * qy; szz += pz * qz
    k = [
        [sxx + syy + szz, syz - szy,       szx - sxz,       sxy - syx],
        [syz - szy,       sxx - syy - szz, sxy + syx,       szx + sxz],
        [szx - sxz,       sxy + syx,      -sxx + syy - szz, syz + szy],
        [sxy - syx,       szx + sxz,       syz + szy,      -sxx - syy + szz],
    ]
    vals, vecs = _jacobi(k)
    best = max(range(4), key=lambda i: vals[i])
    q0, q1, q2, q3 = (vecs[r][best] for r in range(4))
    nrm = math.sqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3) or 1.0
    q0, q1, q2, q3 = q0 / nrm, q1 / nrm, q2 / nrm, q3 / nrm
    r: Rot = (
        (q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3,
         2.0 * (q1 * q2 - q0 * q3),
         2.0 * (q1 * q3 + q0 * q2)),
        (2.0 * (q1 * q2 + q0 * q3),
         q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3,
         2.0 * (q2 * q3 - q0 * q1)),
        (2.0 * (q1 * q3 - q0 * q2),
         2.0 * (q2 * q3 + q0 * q1),
         q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3),
    )
    t: Vec3 = tuple(  # type: ignore[assignment]
        cr[i] - sum(r[i][j] * cm[j] for j in range(3)) for i in range(3)
    )
    ss = 0.0
    for p, q in zip(mob, ref):
        x = tuple(sum(r[i][j] * p[j] for j in range(3)) + t[i] for i in range(3))
        ss += sum((x[i] - q[i]) ** 2 for i in range(3))
    return r, t, math.sqrt(ss / n)


def apply_transform(r: Rot, t: Vec3, p: Vec3) -> Vec3:
    return (
        r[0][0] * p[0] + r[0][1] * p[1] + r[0][2] * p[2] + t[0],
        r[1][0] * p[0] + r[1][1] * p[1] + r[1][2] * p[2] + t[1],
        r[2][0] * p[0] + r[2][1] * p[1] + r[2][2] * p[2] + t[2],
    )


def _dist(a: Vec3, b: Vec3) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                     + (a[2] - b[2]) ** 2)


def _centroid(pts: Sequence[Vec3]) -> Vec3:
    n = len(pts) or 1
    return (sum(p[0] for p in pts) / n,
            sum(p[1] for p in pts) / n,
            sum(p[2] for p in pts) / n)


# --------------------------------------------------------------------------
# TM-score
# --------------------------------------------------------------------------

def _d0(length: int) -> float:
    if length <= 15:
        return 0.5
    return max(0.5, 1.24 * (length - 15) ** (1.0 / 3.0) - 1.8)


@dataclass
class Alignment:
    """The structural superposition and everything read off it."""

    rot: Rot
    trans: Vec3
    tm_score: float
    rmsd_a: float
    aligned_length: int
    #: donor residue index -> target residue index, for the core pairs only
    core_pairs: list[tuple[int, int]]
    norm_length: int
    #: The TM-score NUMERATOR, sum 1/(1+(d/d0)^2) over seeded pairs — the
    #: "number of structurally equivalent residues". Normalisation-free, so it
    #: is what chain maps are ranked on. See `align_structures`.
    mass: float = 0.0
    chain_map: dict[str, str] = field(default_factory=dict)
    #: Set when several chain maps tied on mass and something else chose.
    tie_broken_on: str | None = None


def tm_superpose(
    donor_ca: Sequence[Vec3],
    target_ca: Sequence[Vec3],
    seed_pairs: Sequence[tuple[int, int]],
    norm_length: int,
    d0: float | None = None,
) -> Alignment:
    """TM-score-optimal superposition of donor onto target over `seed_pairs`.

    The published TM-score search: seed on a contiguous fragment, superpose,
    keep the pairs inside a distance cut, re-superpose, iterate to a fixed
    point, and sweep fragment length and starting offset. The correspondence is
    *fixed* by `seed_pairs` (it comes from the sequence alignment); what is
    searched is which subset of it defines the rotation. That is the part that
    matters here — it is what lets a fusion construct, a disordered terminus or
    a hinged domain fall out of the fit instead of dragging it.

    `d0` may be supplied so that several candidate chain maps are scored on one
    scale; otherwise it is derived from `norm_length` as usual.
    """
    if len(seed_pairs) < 3:
        raise ValueError("fewer than 3 seed pairs")
    if d0 is None:
        d0 = _d0(norm_length)
    n = len(seed_pairs)
    best_tm = -1.0
    best: tuple[Rot, Vec3] | None = None

    lengths = sorted({n, max(4, n // 2), max(4, n // 4), max(4, n // 8), 4},
                     reverse=True)
    attempts = 0
    for li in lengths:
        if li > n:
            continue
        stride = max(1, li // 2)
        for start in range(0, n - li + 1, stride):
            if attempts > 80:
                break
            attempts += 1
            sub = list(seed_pairs[start:start + li])
            prev_key: tuple[int, ...] | None = None
            for _ in range(30):
                if len(sub) < 3:
                    break
                try:
                    r, t, _ = superpose(
                        [donor_ca[i] for i, _ in sub],
                        [target_ca[j] for _, j in sub],
                    )
                except ValueError:
                    break
                dists = [
                    _dist(apply_transform(r, t, donor_ca[i]), target_ca[j])
                    for i, j in seed_pairs
                ]
                tm = sum(1.0 / (1.0 + (d / d0) ** 2) for d in dists)
                if tm > best_tm:
                    best_tm, best = tm, (r, t)
                cut = d0 + 1.0
                nxt: list[tuple[int, int]] = []
                while cut < 20.0:
                    nxt = [p for p, d in zip(seed_pairs, dists) if d < cut]
                    if len(nxt) >= 4:
                        break
                    cut += 0.5
                key = tuple(i for i, _ in nxt)
                if key == prev_key or len(nxt) < 3:
                    break
                prev_key, sub = key, nxt

    if best is None:
        raise ValueError("superposition search found no solution")
    r, t = best
    core: list[tuple[int, int]] = []
    ss = 0.0
    for i, j in seed_pairs:
        d = _dist(apply_transform(r, t, donor_ca[i]), target_ca[j])
        if d < CORE_PAIR_A:
            core.append((i, j))
            ss += d * d
    rmsd = math.sqrt(ss / len(core)) if core else float("nan")
    return Alignment(
        rot=r, trans=t, tm_score=round(best_tm / norm_length, 4),
        rmsd_a=round(rmsd, 3) if core else float("nan"),
        aligned_length=len(core), core_pairs=core, norm_length=norm_length,
        mass=best_tm,
    )


# --------------------------------------------------------------------------
# Chain correspondence
# --------------------------------------------------------------------------

def _seq_of(residues: Sequence[Residue]) -> str:
    return "".join(r.one_letter for r in residues)


def _pair_identity(a: Sequence[Residue], b: Sequence[Residue]) -> float:
    pairs = align_sequences(_seq_of(a), _seq_of(b))
    if not pairs:
        return 0.0
    sa, sb = _seq_of(a), _seq_of(b)
    same = sum(1 for i, j in pairs if sa[i] == sb[j])
    return same / max(1, min(len(sa), len(sb)))


#: Hard cap on how many chain assignments are superposed. 60 runs to a few
#: seconds; the cap exists so a 6-chain-versus-6-chain input cannot silently
#: turn into 40 minutes of work.
MAX_CHAIN_MAPS = 60

#: Two chains are held to be in the same biological unit when at least this
#: many residues of one contact the other within `CONTACT_SHELL_A`.
#:
#: PROPOSED, NOT CALIBRATED, and stated with the measurement behind it. An
#: asymmetric unit is a crystallographic object, not a biological one, and the
#: distinction is load-bearing here: 2AZ5 is **two TNF dimers**, not one
#: tetramer, and a chain map allowed to span them put the transferred CD40L
#: ligand 8.79 A from SPD304 with a contact-shell Jaccard of 0.171 while
#: passing all three guards. Measured interface sizes, in residues:
#:
#:     2AZ5  A-B 24   C-D 23   |  B-D  9   A-C 0   A-D 0   B-C 0
#:     3LKJ  A-B 22   B-C 18   |  A-C  3
#:     9SQX  A-B 52   C-D 41   |  A-D 12   A-C 2   B-D 1   B-C 0
#:
#: Every biological interface in that set is >= 18 and every packing contact
#: is <= 12. 15 sits in the gap. It rests on three entries, so it is a
#: proposal — but note the failure mode is asymmetric and mild: too high
#: splits a real oligomer into single chains, which loses an interface anchor
#: and says so; too low re-admits the packing artifact above.
MIN_INTERFACE_RESIDUES = 15

#: Chain maps whose aligned mass is within this fraction of the best are held
#: to be TIED, and the tie is broken on donor-ligand contact-shell coverage.
#:
#: This is the homo-oligomer ambiguity one level down, and it is not cosmetic.
#: Superposing a CD40L trimer onto a TNF dimer has three distinct donor-chain
#: pairs and two orientations each, and on a homotrimer they are near-identical
#: by every sequence and backbone measure — but they place the ligand in
#: different protomer gaps. Ranking on mass alone picked `A->B, B->A`, which
#: leaves the ligand's principal contact chain (donor C, 238 contacting atoms
#: against 85 and 63) out of the fit entirely and lands the ligand 10.3 A from
#: SPD304 with guard 2 sitting exactly on its threshold at 0.50.
#:
#: So among tied maps, prefer the one that actually superposes the chains the
#: donor ligand sits against. That is guard 2's own criterion used to *choose*
#: rather than only to *judge*, which is legitimate here because the donor
#: ligand is the donor's, not the target's — the anchor stays ligand-free with
#: respect to the target, which is the property that matters.
MASS_TIE_FRACTION = 0.90

#: Radius around the donor ligand, in angstrom, whose aligned CA pairs are used
#: to re-fit the transform that actually MOVES the ligand.
#:
#: A global superposition minimises error everywhere and therefore nowhere in
#: particular. The number that matters for a transferred site is the fit *where
#: the ligand is*, and on a two-domain or hinged pair the two can differ by
#: several angstrom. The fold match stays global — `tm_score` and `rmsd_a` are
#: reported from it, and guard 3's TM floor is applied to it — but the
#: coordinates handed downstream come from the local re-fit, and its RMSD is
#: reported separately as `local_rmsd_a`. Reporting only the global number
#: while transferring on the local one would be the same kind of mismatch
#: `CLAUDE.md` rule 5 records for the two C-alpha displacement protocols.
LOCAL_REFINE_A = 15.0

#: Below this many pairs the local re-fit is not attempted — a transform fitted
#: to a handful of atoms is free to say anything. The global one is used and
#: `local_refit` records that it was skipped.
MIN_LOCAL_REFINE_PAIRS = 20


def _chain_maps(donor_chains: list[str], target_chains: list[str],
                limit: int = MAX_CHAIN_MAPS) -> list[dict[str, str]]:
    """Every injective donor-chain -> target-chain assignment, PARTIAL MAPS
    INCLUDED, largest first.

    Two things this has to get right, both measured on 2AZ5:

    **Sequence identity cannot rank chain maps for a homo-oligomer.** Every
    pairing of a homotrimer's chains scores identically by construction, so
    "take the top 3 by identity" is an arbitrary tie-break wearing a score.
    This is the same defect class `CLAUDE.md` records as
    `site_signature_unreliable_homooligomer`. The maps are therefore ranked by
    the *structural* fit, not the sequence one — see `align_structures`.

    **The map must be allowed to be partial.** 2AZ5 is not a TNF trimer: it is
    **two TNF dimers**, A+B and C+D (CA-centroid separations 22.1 A within a
    pair against 31-39 A between), because SPD304 occupies the third
    protomer's place — that *is* the mechanism. Forcing a 3-chain CD40L donor
    onto 3 of the 4 target chains produced the map A->C, B->B, C->D, which
    **spans both dimers**: TM 0.373, RMSD 3.09 A, 36 clashes, and a rejected
    transfer on the one case that is supposed to pass. A donor whose oligomer
    is more complete than the target's must be allowed to leave a chain
    unmatched.
    """
    from itertools import combinations, permutations

    if not donor_chains or not target_chains:
        return []
    out: list[dict[str, str]] = []
    kmax = min(len(donor_chains), len(target_chains))
    for k in range(kmax, 0, -1):
        for dsub in combinations(donor_chains, k):
            for perm in permutations(target_chains, k):
                out.append(dict(zip(dsub, perm)))
                if len(out) >= limit:
                    return out
    return out


def _atom_grid(pts: Sequence[tuple[Any, Vec3]], cell: float
               ) -> dict[tuple[int, int, int], list[tuple[Any, Vec3]]]:
    g: dict[tuple[int, int, int], list[tuple[Any, Vec3]]] = {}
    for tag, p in pts:
        g.setdefault(
            (int(p[0] // cell), int(p[1] // cell), int(p[2] // cell)), []
        ).append((tag, p))
    return g


def chain_components(chains: dict[str, list[Residue]]) -> list[list[str]]:
    """Group chains into biological units by interface size.

    Connected components over the graph whose edges are chain pairs sharing at
    least `MIN_INTERFACE_RESIDUES` contacting residues. This is what stops a
    chain map spanning two copies of the same oligomer in one asymmetric unit —
    see `MIN_INTERFACE_RESIDUES` for the measurement and for what it cost when
    it was missing.
    """
    names = sorted(chains)
    if len(names) < 2:
        return [names] if names else []
    cell = CONTACT_SHELL_A
    pts: list[tuple[tuple[str, str], Vec3]] = []
    for c in names:
        for r in chains[c]:
            for a in r.atoms.values():
                pts.append(((c, r.key), a))
    grid = _atom_grid(pts, cell)
    r2 = cell * cell
    touch: dict[tuple[str, str], set[str]] = {}
    for (c, rkey), a in pts:
        gi = (int(a[0] // cell), int(a[1] // cell), int(a[2] // cell))
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for (c2, _), b in grid.get(
                        (gi[0] + dx, gi[1] + dy, gi[2] + dz), ()
                    ):
                        if c2 == c:
                            continue
                        if ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                                + (a[2] - b[2]) ** 2) <= r2:
                            touch.setdefault((c, c2), set()).add(rkey)
    parent = {c: c for c in names}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for (a_, b_), rs in touch.items():
        if len(rs) >= MIN_INTERFACE_RESIDUES:
            ra, rb = find(a_), find(b_)
            if ra != rb:
                parent[ra] = rb
    groups: dict[str, list[str]] = {}
    for c in names:
        groups.setdefault(find(c), []).append(c)
    return sorted(groups.values(), key=lambda g: (-len(g), g))


def align_structures(donor: Structure, target: Structure, *,
                     max_chains: int = 4,
                     donor_chains: Sequence[str] | None = None,
                     target_chains: Sequence[str] | None = None,
                     prefer_ligand: Residue | None = None) -> Alignment:
    """Best superposition of `donor` onto `target`.

    Multi-chain by construction. The single-chain shortcut is wrong for every
    interface and every axial site — the two shapes this anchor is most often
    used on — because a rotation fitted to one protomer places an interface
    ligand relative to that protomer alone.

    **Chain maps are ranked on the TM-score NUMERATOR, not on TM-score.** The
    numerator, `sum 1/(1+(d/d0)^2)` over seeded pairs, is the number of
    structurally equivalent residues and carries no normalisation, so maps
    covering different numbers of chains are comparable on it — and because it
    rewards coverage as well as accuracy, it will not collapse to a
    one-chain map the way maximising a per-map-normalised TM-score would.
    `d0` is fixed from the whole target so that every map is scored on one
    scale. The reported `tm_score` is then the numerator over the residues of
    the *mapped* target chains, which is the honest denominator: it answers
    "how well does the donor explain the part of the target it matched",
    and the part it did not match is visible in `chain_map`.
    """
    dch = donor.polymer_chains()
    tch = target.polymer_chains()
    if donor_chains is not None:
        dch = {c: r for c, r in dch.items() if c in set(donor_chains)}
    if target_chains is not None:
        tch = {c: r for c, r in tch.items() if c in set(target_chains)}
    if not dch or not tch:
        raise ValueError("no polymer chain of usable length in donor or target")
    dnames = sorted(dch, key=lambda c: -len(dch[c]))[:max_chains]
    tnames = sorted(tch, key=lambda c: -len(tch[c]))[:max_chains]

    seq_pairs: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for d in dnames:
        for t in tnames:
            seq_pairs[(d, t)] = align_sequences(_seq_of(dch[d]), _seq_of(tch[t]))

    # Flat index over donor / target residues, so pairs from several chains
    # can be pooled into one superposition.
    d_index: list[Residue] = []
    d_off: dict[str, int] = {}
    for c in dnames:
        d_off[c] = len(d_index)
        d_index.extend(dch[c])
    t_index: list[Residue] = []
    t_off: dict[str, int] = {}
    for c in tnames:
        t_off[c] = len(t_index)
        t_index.extend(tch[c])
    d_ca = [r.ca() for r in d_index]
    t_ca = [r.ca() for r in t_index]

    # A chain map may not span two biological units. See
    # `MIN_INTERFACE_RESIDUES`: 2AZ5 is two TNF dimers and a map allowed to
    # cross them silently produced a passing-but-wrong anchor.
    d_units = chain_components({c: dch[c] for c in dnames})
    t_units = chain_components({c: tch[c] for c in tnames})
    candidate_maps: list[dict[str, str]] = []
    seen: set[tuple[tuple[str, str], ...]] = set()
    for du in d_units:
        for tu in t_units:
            for m in _chain_maps(
                sorted(du, key=lambda c: -len(dch[c])),
                sorted(tu, key=lambda c: -len(tch[c])),
            ):
                key = tuple(sorted(m.items()))
                if key not in seen:
                    seen.add(key)
                    candidate_maps.append(m)

    d0 = _d0(len(t_index))
    scored: list[Alignment] = []
    for cmap in candidate_maps:
        pooled: list[tuple[int, int]] = []
        for d, t in cmap.items():
            for i, j in seq_pairs[(d, t)]:
                pooled.append((d_off[d] + i, t_off[t] + j))
        if len(pooled) < 3:
            continue
        pooled.sort()
        norm = sum(len(tch[t]) for t in cmap.values())
        try:
            al = tm_superpose(d_ca, t_ca, pooled, norm, d0)  # type: ignore[arg-type]
        except ValueError:
            continue
        al.chain_map = cmap
        al._donor_index = d_index    # type: ignore[attr-defined]
        al._target_index = t_index   # type: ignore[attr-defined]
        scored.append(al)
    if not scored:
        raise ValueError("no chain assignment produced a superposition")

    top = max(a.mass for a in scored)
    tied = [a for a in scored if a.mass >= top * MASS_TIE_FRACTION]
    if prefer_ligand is not None and len(tied) > 1:
        shell = {r.key for r in contact_shell(donor, prefer_ligand)}
        if shell:
            def _cover(a: Alignment) -> tuple[float, float]:
                aligned = {d_index[i].key for i, _ in a.core_pairs}
                return (len(shell & aligned) / len(shell), a.mass)
            tied.sort(key=_cover, reverse=True)
            best = tied[0]
            best.tie_broken_on = (   # type: ignore[attr-defined]
                f"donor-ligand shell coverage among {len(tied)} chain maps "
                f"tied within {MASS_TIE_FRACTION:.0%} of the best aligned mass"
            )
            return best
    return max(tied, key=lambda a: a.mass)


# --------------------------------------------------------------------------
# Guards
# --------------------------------------------------------------------------

@dataclass
class GuardResult:
    name: str
    passed: bool
    detail: str
    measured: dict[str, Any] = field(default_factory=dict)

    def as_dict(self) -> dict[str, Any]:
        return {"guard": self.name, "passed": self.passed,
                "detail": self.detail, "measured": self.measured}


def guard_donor_ligand_druglike(
    comp_id: str, *, chemcomps: Any = None
) -> GuardResult:
    """GUARD 1. The donor's ligand must pass the drug-like classifier.

    Measured failure: NLRP3's only constructible donor was NOD2 (5IRN) carrying
    **ADP**. The transfer validated to 0.68 A against NLRP3's own ADP — it was
    a *good* transfer — and it selected the nucleotide lobe, not the drug site.
    A cofactor site is a real site; it is not the site a small-molecule dossier
    is asking about, and nothing downstream of a transfer can tell the two
    apart once the coordinates exist.
    """
    import ligand_filter as lf

    src = chemcomps if chemcomps is not None else RcsbChemComps()
    rec = src.get(comp_id)
    if rec is None:
        why = getattr(src, "fetch_errors", {}).get(comp_id.upper())
        return GuardResult(
            name="donor_ligand_druglike", passed=False,
            detail=(
                f"{comp_id}: chem-comp record could not be retrieved"
                + (f" ({why})" if why else "")
                + " — NOT CHECKED, not a pass. A lookup failure must not"
                  " render as drug-like."
            ),
            measured={"comp_id": comp_id, "verdict": None,
                      "lookup_failed": True},
        )
    v = lf.classify_record(rec, comp_id)
    ok = v.verdict not in REJECTED_LIGAND_VERDICTS
    return GuardResult(
        name="donor_ligand_druglike", passed=ok,
        detail=(
            f"{comp_id} classified {v.verdict} ({v.confidence}): {v.reason}"
            + ("" if ok else "  -> REJECTED: a perfect transfer of a"
                             " non-drug-like ligand is a confident wrong site")
        ),
        measured={
            "comp_id": comp_id, "verdict": v.verdict,
            "confidence": v.confidence, "flags": list(v.flags or ()),
            "rejected_verdicts": sorted(REJECTED_LIGAND_VERDICTS),
        },
    )


def contact_shell(donor: Structure, lig: Residue,
                  *, radius: float = CONTACT_SHELL_A) -> list[Residue]:
    """Donor polymer residues with a heavy atom within `radius` of `lig`."""
    lig_atoms = lig.heavy()
    if not lig_atoms:
        return []
    xs = [p[0] for p in lig_atoms]
    ys = [p[1] for p in lig_atoms]
    zs = [p[2] for p in lig_atoms]
    lo = (min(xs) - radius, min(ys) - radius, min(zs) - radius)
    hi = (max(xs) + radius, max(ys) + radius, max(zs) + radius)
    r2 = radius * radius
    out: list[Residue] = []
    for res in donor.residues:
        if not res.is_polymer:
            continue
        hit = False
        for a in res.atoms.values():
            if not (lo[0] <= a[0] <= hi[0] and lo[1] <= a[1] <= hi[1]
                    and lo[2] <= a[2] <= hi[2]):
                continue
            for b in lig_atoms:
                if ((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                        + (a[2] - b[2]) ** 2) <= r2:
                    hit = True
                    break
            if hit:
                break
        if hit:
            out.append(res)
    return out


def guard_domain_attribution(
    donor: Structure, lig: Residue, al: Alignment,
) -> GuardResult:
    """GUARD 2. The donor ligand's contact shell must overlap the ALIGNED REGION.

    Not the same chain — the aligned region. Two measured failures that
    chain-level attribution passes:

      * **7KRZ.** Bortezomib (`BO2`) is on the correct LONP1 chain, at auth
        768-898. The NACHT-aligned region is auth 506-721. Right chain, wrong
        domain, and the rotation that would carry the ligand was fitted 250
        residues away from it.
      * **The proteasome hits.** They never carry their entry-title drug on the
        AAA+ chain that produced the fold match at all.

    `SKILL.md` already documents that `entry_ligands` cannot attribute a ligand
    to a chain, and that the same bug exists one level up on the protein side.
    This is the third level: a chain is not a domain.
    """
    d_index: list[Residue] = getattr(al, "_donor_index", [])
    aligned_keys = {d_index[i].key for i, _ in al.core_pairs}
    shell = contact_shell(donor, lig)
    if not shell:
        return GuardResult(
            name="domain_attribution", passed=False,
            detail=(f"{lig.comp_id} {lig.key}: no donor polymer residue within "
                    f"{CONTACT_SHELL_A} A — the ligand contacts nothing, so "
                    "there is no site to transfer"),
            measured={"shell_size": 0, "fraction_in_aligned_region": None},
        )
    inside = [r for r in shell if r.key in aligned_keys]
    frac = len(inside) / len(shell)
    ok = frac >= MIN_SHELL_IN_ALIGNED_FRACTION

    def _rng(rs: Sequence[Residue]) -> str:
        nums = sorted(x for x in (r.seq_int() for r in rs) if x is not None)
        return f"{nums[0]}-{nums[-1]}" if nums else "n/a"

    aligned_res = [d_index[i] for i, _ in al.core_pairs]
    return GuardResult(
        name="domain_attribution", passed=ok,
        detail=(
            f"{lig.comp_id} {lig.key}: {len(inside)}/{len(shell)} "
            f"({frac:.2f}) contact-shell residues lie in the aligned region; "
            f"shell spans auth {_rng(shell)} on chain(s) "
            f"{sorted({r.chain for r in shell})}, aligned region spans auth "
            f"{_rng(aligned_res)} on chain(s) "
            f"{sorted({r.chain for r in aligned_res})}"
            + ("" if ok else
               f"  -> REJECTED: below {MIN_SHELL_IN_ALIGNED_FRACTION:.2f}."
               " Right chain is not right domain; the transform was fitted"
               " somewhere the ligand is not.")
        ),
        measured={
            "shell_size": len(shell),
            "shell_in_aligned_region": len(inside),
            "fraction_in_aligned_region": round(frac, 3),
            "threshold": MIN_SHELL_IN_ALIGNED_FRACTION,
            "shell_auth_range": _rng(shell),
            "aligned_auth_range": _rng(aligned_res),
            "shell_chains": sorted({r.chain for r in shell}),
            "aligned_chains": sorted({r.chain for r in aligned_res}),
        },
    )


def guard_alignment_and_sterics(
    al: Alignment, placed: Sequence[Vec3], target: Structure,
    *, chains: Iterable[str] | None = None,
) -> GuardResult:
    """GUARD 3. TM floor, RMSD ceiling, and — the strong one — a steric check.

    Measured failure: IL-2 forced onto IL-17A at TM 0.254-0.274 put the
    transferred ligand 21.59 A from the site and **inside the protein** — 29
    atoms within 2.0 A of protein heavy atoms, closest contact 0.51 A. Nothing
    announced it.

    The three tests are not redundant and are reported separately even when one
    of them already voids the transfer, because *which* one fired is the
    finding. An alignment that clears the TM floor and still clashes is a
    different fault from one that never had a fold match.

    **The clash count is attributed by target chain, and this matters.** On an
    oligomer-destabilisation mechanism, clashing with a protomer is not an
    error — it *is* the mechanism. `CLAUDE.md` rule 5 records TNF-alpha's site
    as 40 of 66 clashes coming from the subunit the ligand displaces. So the
    guard still fails, because a set of coordinates overlapping protein is not
    a site as it stands, but `clash_by_chain` and
    `clashes_confined_to_one_chain` are reported so a reader can tell "the
    ligand is inside the backbone" from "the ligand is where a subunit this
    hypothesis proposes to displace currently sits". Those need different
    escalations and only one of them is a bug.
    """
    keep = set(chains) if chains is not None else None
    prot: list[tuple[str, bool, Vec3]] = []
    for res in target.residues:
        if not res.is_polymer:
            continue
        if keep is not None and res.chain not in keep:
            continue
        for nm, a in res.atoms.items():
            prot.append((res.chain, nm in BACKBONE_ATOMS, a))

    clashing = bb_clash = sc_clash = 0
    closest = closest_bb = float("inf")
    by_chain: dict[str, int] = {}
    if prot and placed:
        cell = CLASH_DISTANCE_A
        grid: dict[tuple[int, int, int], list[tuple[str, bool, Vec3]]] = {}
        for ch, isbb, p in prot:
            grid.setdefault(
                (int(p[0] // cell), int(p[1] // cell), int(p[2] // cell)), []
            ).append((ch, isbb, p))
        for a in placed:
            gi = (int(a[0] // cell), int(a[1] // cell), int(a[2] // cell))
            hit: set[str] = set()
            hit_bb = False
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for ch, isbb, p in grid.get(
                            (gi[0] + dx, gi[1] + dy, gi[2] + dz), ()
                        ):
                            d = _dist(a, p)
                            if d < closest:
                                closest = d
                            if isbb and d < closest_bb:
                                closest_bb = d
                            if d < CLASH_DISTANCE_A:
                                hit.add(ch)
                                hit_bb = hit_bb or isbb
            if hit:
                clashing += 1
                if hit_bb:
                    bb_clash += 1
                else:
                    sc_clash += 1
                for ch in hit:
                    by_chain[ch] = by_chain.get(ch, 0) + 1

    tm_ok = al.tm_score >= MIN_TM_SCORE
    rmsd_ok = (al.rmsd_a == al.rmsd_a) and al.rmsd_a <= MAX_ALIGNMENT_RMSD_A
    steric_ok = bb_clash <= MAX_BACKBONE_CLASH_ATOMS
    ok = tm_ok and rmsd_ok and steric_ok
    fails = []
    if not tm_ok:
        fails.append(f"TM {al.tm_score:.3f} < {MIN_TM_SCORE}")
    if not rmsd_ok:
        fails.append(f"RMSD {al.rmsd_a} A > {MAX_ALIGNMENT_RMSD_A}")
    if not steric_ok:
        fails.append(
            f"{bb_clash} transferred ligand atom(s) within "
            f"{CLASH_DISTANCE_A} A of protein BACKBONE atoms "
            f"(closest backbone contact {closest_bb:.2f} A) — no rotamer moves "
            "a backbone, so this is physically impossible and the transfer is "
            "void regardless of how the alignment scored"
        )
    return GuardResult(
        name="alignment_and_sterics", passed=ok,
        detail=(
            f"TM {al.tm_score:.3f} (floor {MIN_TM_SCORE}), RMSD {al.rmsd_a} A "
            f"over {al.aligned_length} CA pairs (ceiling "
            f"{MAX_ALIGNMENT_RMSD_A} A), backbone clashes {bb_clash} "
            f"(max {MAX_BACKBONE_CLASH_ATOMS}), side-chain clashes {sc_clash} "
            "(reported, never rejects)"
            + ("" if ok else "  -> REJECTED: " + "; ".join(fails))
            + (""
               if al.rmsd_a != al.rmsd_a or al.rmsd_a <= MARGINAL_ALIGNMENT_RMSD_A
               else f"  [alignment_marginal: RMSD exceeds "
                    f"{MARGINAL_ALIGNMENT_RMSD_A} A]")
            + ("" if not sc_clash else
               f"  [sidechain_occlusion: {sc_clash} atom(s) overlap side "
               "chains, resolvable by repacking — but CLAUDE.md rule 5 puts "
               "the prognosis for a side-chain-occluded site at "
               "micromolar-at-best, so this is a finding, not noise]")
        ),
        measured={
            "tm_score": al.tm_score, "tm_floor": MIN_TM_SCORE,
            "rmsd_a": al.rmsd_a, "rmsd_ceiling_a": MAX_ALIGNMENT_RMSD_A,
            "aligned_length": al.aligned_length,
            "clash_count": clashing,
            "backbone_clash_count": bb_clash,
            "sidechain_clash_count": sc_clash,
            "clash_distance_a": CLASH_DISTANCE_A,
            "max_backbone_clash_atoms": MAX_BACKBONE_CLASH_ATOMS,
            "closest_contact_a": (round(closest, 3) if closest < float("inf")
                                  else None),
            "closest_backbone_contact_a": (
                round(closest_bb, 3) if closest_bb < float("inf") else None
            ),
            "clash_by_chain": dict(sorted(by_chain.items())),
            "clashes_confined_to_one_chain": (
                len(by_chain) == 1 if by_chain else None
            ),
            "tm_passed": tm_ok, "rmsd_passed": rmsd_ok,
            "steric_passed": steric_ok,
            "alignment_marginal": bool(
                al.rmsd_a == al.rmsd_a and al.rmsd_a > MARGINAL_ALIGNMENT_RMSD_A
            ),
            "sidechain_occlusion": sc_clash > 0,
        },
    )


# --------------------------------------------------------------------------
# The transfer
# --------------------------------------------------------------------------

@dataclass
class TransferResult:
    accepted: bool
    guards: list[GuardResult]
    provenance: dict[str, Any]
    #: Target residue keys ("A/57") lining the transferred ligand. The anchor.
    site_residues: list[str] = field(default_factory=list)
    transferred_centroid: Vec3 | None = None
    transferred_atoms: list[Vec3] = field(default_factory=list)
    validation: dict[str, Any] = field(default_factory=dict)

    #: The six fields that must travel with every transferred-homolog anchor.
    REQUIRED_PROVENANCE = (
        "donor_pdb_id", "donor_ligand", "tm_score", "rmsd_a",
        "aligned_length", "clash_count", "backbone_clash_count",
    )

    @property
    def rejected_by(self) -> list[str]:
        return [g.name for g in self.guards if not g.passed]

    def assert_reportable(self) -> None:
        """Raise unless the six mandatory provenance fields are present.

        **Nothing may fire silently.** A transferred-homolog anchor without its
        donor, its ligand, its fit quality and its clash count is not weak
        evidence — it is unusable, because there is no way for a reader to
        discount it. Call this before writing an anchor into a dossier.
        """
        missing = [k for k in self.REQUIRED_PROVENANCE
                   if k not in self.provenance]
        if missing:
            raise ValueError(
                "transferred_homolog_site is not reportable — missing "
                f"provenance: {missing}"
            )

    def as_dict(self) -> dict[str, Any]:
        self.assert_reportable()
        return {
            "anchor": "transferred_homolog_site",
            "accepted": self.accepted,
            "rejected_by": self.rejected_by,
            "provenance": self.provenance,
            "guards": [g.as_dict() for g in self.guards],
            "site_residues": self.site_residues,
            "transferred_centroid": (
                [round(c, 3) for c in self.transferred_centroid]
                if self.transferred_centroid else None
            ),
            "validation": self.validation or None,
        }


def transfer_homolog_site(
    target_pdb_id: str,
    donor_pdb_id: str,
    *,
    donor_ligand: str | None = None,
    target_structure: Structure | None = None,
    donor_structure: Structure | None = None,
    chemcomps: Any = None,
    reference_ligand: str | None = None,
    donor_chains: Sequence[str] | None = None,
    target_chains: Sequence[str] | None = None,
) -> TransferResult:
    """Transfer a donor's ligand site onto the target, under all three guards.

    `donor_ligand` may be omitted, in which case the most-copied non-water het
    component is used and guard 1 judges it. That default is deliberate: the
    common real failure is not "we picked the wrong ligand", it is "the entry
    had exactly one ligand and it was ADP".

    `reference_ligand` is a ligand *in the target*, used only to VALIDATE the
    result after the guards have run. It is never an input to a guard — a
    ligand-free anchor that consults the target's own ligand is not a
    ligand-free anchor.
    """
    target = target_structure or load_structure(target_pdb_id)
    donor = donor_structure or load_structure(donor_pdb_id)

    prov: dict[str, Any] = {
        "donor_pdb_id": donor_pdb_id.upper(),
        "donor_ligand": None,
        "tm_score": None,
        "rmsd_a": None,
        "aligned_length": None,
        "clash_count": None,
        "backbone_clash_count": None,
        "target_pdb_id": target_pdb_id.upper(),
        "method": ("stdlib TM-score superposition seeded by BLOSUM62 affine "
                   "global alignment, multi-chain, TM normalised by target"),
        "thresholds": {
            "min_tm_score": MIN_TM_SCORE,
            "max_alignment_rmsd_a": MAX_ALIGNMENT_RMSD_A,
            "clash_distance_a": CLASH_DISTANCE_A,
            "max_clash_atoms": MAX_CLASH_ATOMS,
            "min_shell_in_aligned_fraction": MIN_SHELL_IN_ALIGNED_FRACTION,
        },
    }
    guards: list[GuardResult] = []

    # ---- pick the donor ligand ----------------------------------------
    if donor_ligand is None:
        cands = donor.het_comp_ids()
        if not cands:
            prov["donor_ligand"] = None
            guards.append(GuardResult(
                name="donor_ligand_druglike", passed=False,
                detail=f"{donor_pdb_id}: no non-water het component in the entry",
                measured={"comp_id": None},
            ))
            prov["clash_count"] = None
            return TransferResult(False, guards, prov)
        donor_ligand = cands[0]
    prov["donor_ligand"] = donor_ligand.upper()

    # ---- GUARD 1 -------------------------------------------------------
    g1 = guard_donor_ligand_druglike(donor_ligand, chemcomps=chemcomps)
    guards.append(g1)

    # ---- alignment (needed for provenance even on a guard-1 rejection) --
    # The chain map is chosen with the donor ligand in view, so that a tie
    # between symmetry-equivalent maps is broken toward the one that actually
    # superposes the chains the ligand sits against. See `MASS_TIE_FRACTION`.
    _insts = donor.het_instances(donor_ligand)
    _prefer = (max(_insts, key=lambda r: len(contact_shell(donor, r)))
               if _insts else None)
    try:
        al = align_structures(donor, target, donor_chains=donor_chains,
                              target_chains=target_chains,
                              prefer_ligand=_prefer)
    except ValueError as exc:
        prov["alignment_error"] = str(exc)
        return TransferResult(False, guards, prov)
    prov["chain_map_tie_broken_on"] = al.tie_broken_on
    prov["tm_score"] = al.tm_score
    prov["rmsd_a"] = al.rmsd_a
    prov["aligned_length"] = al.aligned_length
    prov["chain_map_donor_to_target"] = al.chain_map
    prov["tm_normalised_by"] = (
        f"{al.norm_length} residues in mapped target chain(s) "
        f"{sorted(set(al.chain_map.values()))}"
    )

    instances = donor.het_instances(donor_ligand)
    if not instances:
        prov["clash_count"] = None
        guards.append(GuardResult(
            name="domain_attribution", passed=False,
            detail=f"{donor_ligand} does not occur in {donor_pdb_id}",
            measured={"shell_size": 0},
        ))
        return TransferResult(False, guards, prov)

    # ---- GUARD 2, over every copy; keep the best-attributed one ---------
    scored = []
    for lig in instances:
        g = guard_domain_attribution(donor, lig, al)
        scored.append((g.measured.get("fraction_in_aligned_region") or 0.0,
                       lig, g))
    scored.sort(key=lambda x: -x[0])
    _, lig, g2 = scored[0]
    prov["donor_ligand_instance"] = lig.key
    prov["donor_ligand_copies"] = len(instances)
    guards.append(g2)

    # ---- place it, on a transform re-fitted AROUND THE LIGAND -----------
    d_index: list[Residue] = getattr(al, "_donor_index", [])
    t_index: list[Residue] = getattr(al, "_target_index", [])
    rot, trans = al.rot, al.trans
    lig_c = _centroid(lig.heavy())
    local = [
        (i, j) for i, j in al.core_pairs
        if _dist(d_index[i].ca(), lig_c) <= LOCAL_REFINE_A  # type: ignore[arg-type]
    ]
    if len(local) >= MIN_LOCAL_REFINE_PAIRS:
        rot, trans, lrms = superpose(
            [d_index[i].ca() for i, _ in local],   # type: ignore[misc]
            [t_index[j].ca() for _, j in local],   # type: ignore[misc]
        )
        prov["local_refit"] = {
            "n_pairs": len(local), "radius_a": LOCAL_REFINE_A,
            "local_rmsd_a": round(lrms, 3),
            "_note": ("the transform that moved the ligand was re-fitted on "
                      "the aligned pairs within this radius of the donor "
                      "ligand; tm_score and rmsd_a above are the GLOBAL fold "
                      "match and are what guard 3's TM floor is applied to"),
        }
    else:
        prov["local_refit"] = {
            "n_pairs": len(local), "radius_a": LOCAL_REFINE_A,
            "local_rmsd_a": None,
            "_note": (f"skipped: fewer than {MIN_LOCAL_REFINE_PAIRS} aligned "
                      "pairs near the ligand, global transform used"),
        }
    placed = [apply_transform(rot, trans, a) for a in lig.heavy()]
    centroid = _centroid(placed) if placed else None

    # ---- GUARD 3 -------------------------------------------------------
    # Clashes are counted against the MAPPED target chains. A chain the donor
    # never matched is, as far as this anchor is concerned, a bystander — and
    # `CLAUDE.md` rule 5 records the measured cost of treating one as part of
    # the assembly (TNF chain D touches the chain-A ligand with 3 atoms against
    # 44 and 39 for the real partners, and counting it consumed all three apo
    # chains). The count over every chain is reported beside it.
    mapped = sorted(set(al.chain_map.values()))
    g3 = guard_alignment_and_sterics(al, placed, target, chains=mapped)
    guards.append(g3)
    g3_all = guard_alignment_and_sterics(al, placed, target, chains=None)
    g3.measured["clash_count_all_target_chains"] = (
        g3_all.measured.get("clash_count")
    )
    g3.measured["clash_by_chain_all_target_chains"] = (
        g3_all.measured.get("clash_by_chain")
    )
    prov["clash_count"] = g3.measured.get("clash_count")
    prov["backbone_clash_count"] = g3.measured.get("backbone_clash_count")
    prov["sidechain_clash_count"] = g3.measured.get("sidechain_clash_count")
    prov["closest_backbone_contact_a"] = (
        g3.measured.get("closest_backbone_contact_a")
    )
    prov["clash_by_chain"] = g3.measured.get("clash_by_chain")
    prov["clash_count_all_target_chains"] = (
        g3.measured.get("clash_count_all_target_chains")
    )
    prov["closest_contact_a"] = g3.measured.get("closest_contact_a")
    prov["alignment_marginal"] = g3.measured.get("alignment_marginal")

    accepted = all(g.passed for g in guards)

    # ---- the anchor it produces ---------------------------------------
    site: list[str] = []
    if placed:
        r2 = CONTACT_SHELL_A * CONTACT_SHELL_A
        for res in target.residues:
            if not res.is_polymer:
                continue
            for a in res.atoms.values():
                if any(((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2
                        + (a[2] - b[2]) ** 2) <= r2 for b in placed):
                    site.append(res.key)
                    break

    # ---- validation, AFTER the guards ---------------------------------
    validation: dict[str, Any] = {}
    if reference_ligand:
        refs = target.het_instances(reference_ligand)
        if refs and centroid:
            best_d = min(_dist(centroid, _centroid(r.heavy())) for r in refs)
            ref = min(refs, key=lambda r: _dist(centroid, _centroid(r.heavy())))
            ref_shell = {r.key for r in contact_shell(target, ref)}
            inter = len(set(site) & ref_shell)
            union = len(set(site) | ref_shell) or 1
            validation = {
                "reference_ligand": reference_ligand.upper(),
                "reference_instance": ref.key,
                "centroid_distance_a": round(best_d, 2),
                "contact_shell_jaccard": round(inter / union, 3),
                "_note": ("VALIDATION ONLY — computed after the guards and "
                          "never fed into one. A ligand-free anchor that "
                          "consulted the target's own ligand would not be "
                          "ligand-free."),
            }
        else:
            validation = {"reference_ligand": reference_ligand.upper(),
                          "centroid_distance_a": None,
                          "_note": "reference ligand not present in target"}

    return TransferResult(
        accepted=accepted, guards=guards, provenance=prov,
        site_residues=sorted(site), transferred_centroid=centroid,
        transferred_atoms=placed, validation=validation,
    )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def _main(argv: list[str]) -> int:
    import argparse

    ap = argparse.ArgumentParser(
        prog="homolog_transfer",
        description="Transfer a homolog's ligand site onto a target, under "
                    "three machine-checked guards.",
    )
    ap.add_argument("target", help="target PDB ID")
    ap.add_argument("donor", help="donor PDB ID")
    ap.add_argument("--donor-ligand", default=None)
    ap.add_argument("--reference-ligand", default=None,
                    help="ligand in the TARGET, for post-hoc validation only")
    a = ap.parse_args(argv)
    res = transfer_homolog_site(
        a.target, a.donor, donor_ligand=a.donor_ligand,
        reference_ligand=a.reference_ligand,
    )
    print(json.dumps(res.as_dict(), indent=2))
    return 0 if res.accepted else 1


if __name__ == "__main__":
    import sys

    raise SystemExit(_main(sys.argv[1:]))
