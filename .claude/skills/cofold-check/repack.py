"""repack.py — turn a backbone(+CB) frame into a full-atom structure fpocket can score.

Why this exists
---------------
`bioemu_ensemble` returns **backbone + C-beta only** (~4.9 atoms/residue, residues
zero-indexed, all B-factors 0.00). fpocket and mdpocket define a pocket from
side-chain atoms, so a BioEmu ensemble cannot enter pocket detection at all until
side chains are rebuilt. This module is that step.

Engine: **FASPR** (Huang, Pearce & Zhang, Bioinformatics 2020) — a
backbone-dependent (Dunbrack 2010) rotamer packer solved by a
deterministic/annealed combinatorial search. Single self-contained C++ program,
no runtime dependencies, ~0.1-1 s per 170-residue chain, MIT LICENSE file in the
upstream repo.

NO force-field minimisation is performed and no RDKit conformer work is done.
FASPR is a discrete rotamer search, not a continuous optimiser.

Contract
--------
    repack_pdb(pdb_text, sequence=None) -> dict
    repack_frames(frames, sequence=None) -> list[dict]

See ``SKILL.md`` for the full contract.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile
import urllib.request
import zipfile
from pathlib import Path
from typing import Any

__all__ = [
    "RepackError",
    "faspr_binary",
    "ensure_faspr",
    "repack_pdb",
    "repack_frames",
    "sidechain_rmsd",
]

FASPR_URL = "https://github.com/tommyhuangthu/FASPR/archive/refs/heads/master.zip"
FASPR_LICENCE = "MIT (LICENSE file, upstream repo). README additionally says 'free to academic users' — the two are in tension; the LICENSE file is MIT."

# Backbone atoms FASPR requires present for every residue.
BACKBONE = ("N", "CA", "C", "O")

_THREE = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
_ONE = {v: k for k, v in _THREE.items()}
# Common non-standard names BioEmu / PDB may emit for the 20 canonicals.
_ONE.update({"MSE": "M", "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "SEC": "C"})

# Heavy side-chain atoms per residue, used by sidechain_rmsd. CB included.
SIDECHAIN_ATOMS = {
    "ALA": ["CB"],
    "ARG": ["CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"],
    "ASN": ["CB", "CG", "OD1", "ND2"],
    "ASP": ["CB", "CG", "OD1", "OD2"],
    "CYS": ["CB", "SG"],
    "GLN": ["CB", "CG", "CD", "OE1", "NE2"],
    "GLU": ["CB", "CG", "CD", "OE1", "OE2"],
    "GLY": [],
    "HIS": ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    "ILE": ["CB", "CG1", "CG2", "CD1"],
    "LEU": ["CB", "CG", "CD1", "CD2"],
    "LYS": ["CB", "CG", "CD", "CE", "NZ"],
    "MET": ["CB", "CG", "SD", "CE"],
    "PHE": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"],
    "PRO": ["CB", "CG", "CD"],
    "SER": ["CB", "OG"],
    "THR": ["CB", "OG1", "CG2"],
    "TRP": ["CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"],
    "TYR": ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"],
    "VAL": ["CB", "CG1", "CG2"],
}
# Side chains with a symmetry-equivalent atom swap; RMSD takes the better of the
# two assignments. Without this, PHE/TYR/ASP/GLU ring or carboxylate flips
# inflate the reported error by an amount that has nothing to do with packing.
SYMMETRIC_SWAPS = {
    "PHE": [("CD1", "CD2"), ("CE1", "CE2")],
    "TYR": [("CD1", "CD2"), ("CE1", "CE2")],
    "ASP": [("OD1", "OD2")],
    "GLU": [("OE1", "OE2")],
    "ARG": [("NH1", "NH2")],
    "VAL": [("CG1", "CG2")],
    "LEU": [("CD1", "CD2")],
}


class RepackError(RuntimeError):
    """Raised when repacking cannot be performed or its input is unusable."""


# --------------------------------------------------------------------------
# PDB parsing — deliberately stdlib-only, so this module imports anywhere.
# --------------------------------------------------------------------------
def _parse_atoms(pdb_text: str) -> list[dict[str, Any]]:
    atoms = []
    for line in pdb_text.splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        atoms.append(
            {
                "record": line[:6],
                "serial": line[6:11],
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21],
                "resseq": line[22:26].strip(),
                "icode": line[26],
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
                "element": line[76:78].strip() or line[12:16].strip()[:1],
                "line": line,
            }
        )
    return atoms


def _residues(atoms: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Group atoms into residues, preserving file order."""
    out: list[dict[str, Any]] = []
    index: dict[tuple[str, str, str], dict[str, Any]] = {}
    for a in atoms:
        key = (a["chain"], a["resseq"], a["icode"])
        res = index.get(key)
        if res is None:
            res = {"chain": a["chain"], "resseq": a["resseq"], "icode": a["icode"],
                   "resname": a["resname"], "atoms": {}}
            index[key] = res
            out.append(res)
        res["atoms"][a["name"]] = a
    return out


def _fmt_atom(serial: int, name: str, resname: str, chain: str, resseq: int,
              x: float, y: float, z: float, element: str, bfac: float = 0.0) -> str:
    an = f" {name:<3s}" if len(name) < 4 else name
    return (
        f"ATOM  {serial:5d} {an:4s} {resname:>3s} {chain}{resseq:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{bfac:6.2f}          {element:>2s}"
    )


# --------------------------------------------------------------------------
# FASPR bootstrap
# --------------------------------------------------------------------------
def _default_home() -> Path:
    """Where the FASPR executable and its 13.7 MB rotamer library are cached.

    Deliberately OUTSIDE the source tree: this is a build artifact, and dropping
    a 14 MB binary blob next to a SKILL.md that gets uploaded to the Skills API
    is not a thing to do by default. ``FASPR_HOME`` overrides.
    """
    env = os.environ.get("FASPR_HOME")
    if env:
        return Path(env)
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg) if xdg else Path.home() / ".cache"
    return base / "faspr"


def faspr_binary(home: Path | str | None = None) -> Path | None:
    """Return the FASPR executable if it is already present, else None."""
    h = Path(home) if home else _default_home()
    exe = h / "FASPR"
    lib = h / "dun2010bbdep.bin"
    if exe.exists() and lib.exists() and os.access(exe, os.X_OK):
        return exe
    which = shutil.which("FASPR")
    if which:
        p = Path(which)
        if (p.parent / "dun2010bbdep.bin").exists():
            return p
    return None


def _compile_commands(src: list[str], out: Path) -> list[list[str]]:
    """Compiler invocations to try, in order.

    The plain one first. The explicit-libc++ one is the fallback for macOS boxes
    where CommandLineTools is installed but its libc++ headers live only under
    the SDK — on the machine this was written on, `#include <cmath>` failed
    under a bare `clang++` and succeeded with `-nostdinc++ -isystem <sdk>/usr/include/c++/v1`.
    """
    cmds = [["c++", "-O3", "-o", str(out), *src]]
    if sys.platform == "darwin":
        try:
            sdk = subprocess.run(
                ["xcrun", "--sdk", "macosx", "--show-sdk-path"],
                capture_output=True, text=True, timeout=30, check=False,
            ).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            sdk = ""
        if sdk:
            cmds.append([
                "c++", "-O3", "-nostdinc++",
                "-isystem", f"{sdk}/usr/include/c++/v1", "-isysroot", sdk,
                "-o", str(out), *src,
            ])
    return cmds


def ensure_faspr(home: Path | str | None = None, *, allow_download: bool = True) -> Path:
    """Return a usable FASPR executable, downloading and building it if needed."""
    existing = faspr_binary(home)
    if existing:
        return existing
    h = Path(home) if home else _default_home()
    h.mkdir(parents=True, exist_ok=True)
    if not allow_download:
        raise RepackError(f"FASPR not found under {h} and allow_download=False")

    with tempfile.TemporaryDirectory() as td:
        zp = Path(td) / "faspr.zip"
        try:
            urllib.request.urlopen(FASPR_URL, timeout=120)  # noqa: S310
        except Exception as exc:  # pragma: no cover
            raise RepackError(f"cannot reach FASPR source at {FASPR_URL}: {exc}") from exc
        urllib.request.urlretrieve(FASPR_URL, zp)  # noqa: S310
        with zipfile.ZipFile(zp) as z:
            z.extractall(td)
        root = next(Path(td).glob("FASPR-*"))
        shutil.copy(root / "dun2010bbdep.bin", h / "dun2010bbdep.bin")
        src = sorted(str(p) for p in (root / "src").glob("*.cpp"))
        errs = []
        for cmd in _compile_commands(src, h / "FASPR"):
            r = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if r.returncode == 0 and (h / "FASPR").exists():
                (h / "FASPR").chmod(0o755)
                return h / "FASPR"
            errs.append(f"$ {' '.join(cmd[:4])} ...\n{r.stderr[-600:]}")
        raise RepackError("FASPR build failed:\n" + "\n".join(errs))


# --------------------------------------------------------------------------
# The repack step
# --------------------------------------------------------------------------
def _prepare_input(pdb_text: str, sequence: str | None) -> tuple[str, list[dict], dict]:
    """Rewrite a frame into a PDB FASPR will accept, and report what was changed.

    BioEmu frames are 0-indexed and B-factor-free. FASPR does not care about
    numbering, but everything downstream does, so the mapping back is returned.
    """
    residues = _residues(_parse_atoms(pdb_text))
    if not residues:
        raise RepackError("no ATOM/HETATM records parsed from the frame")

    if sequence is not None:
        seq = "".join(sequence.split()).upper()
        if len(seq) != len(residues):
            raise RepackError(
                f"sequence length {len(seq)} does not match {len(residues)} residues "
                "in the frame. Refusing to pack — a silent off-by-one here renames "
                "every residue downstream of the gap."
            )
        for res, one in zip(residues, seq, strict=True):
            if one not in _THREE:
                raise RepackError(f"non-canonical residue '{one}' in sequence")
            res["resname"] = _THREE[one]
    else:
        for res in residues:
            if res["resname"] not in SIDECHAIN_ATOMS:
                mapped = _ONE.get(res["resname"])
                if mapped is None:
                    raise RepackError(
                        f"residue name '{res['resname']}' is not a canonical amino acid "
                        "and no sequence was supplied to override it"
                    )
                res["resname"] = _THREE[mapped]

    missing = []
    lines = []
    serial = 1
    for i, res in enumerate(residues, start=1):
        absent = [b for b in BACKBONE if b not in res["atoms"]]
        if absent:
            missing.append({"index": i, "resname": res["resname"], "missing": absent})
            continue
        for name in BACKBONE:
            a = res["atoms"][name]
            lines.append(_fmt_atom(serial, name, res["resname"], "A", i,
                                   a["x"], a["y"], a["z"], a["element"] or name[0]))
            serial += 1
        # CB is passed through when present: FASPR ignores input side chains, but
        # keeping it costs nothing and makes the input file self-describing.
        cb = res["atoms"].get("CB")
        if cb is not None and res["resname"] != "GLY":
            lines.append(_fmt_atom(serial, "CB", res["resname"], "A", i,
                                   cb["x"], cb["y"], cb["z"], "C"))
            serial += 1
    if missing:
        raise RepackError(
            f"{len(missing)} residue(s) lack a complete N/CA/C/O backbone, which FASPR "
            f"requires: {missing[:5]}"
        )
    lines.append("TER")
    lines.append("END")
    prep_meta = {
        "n_residues": len(residues),
        "original_numbering": [
            {"index": i, "chain": r["chain"], "resseq": r["resseq"], "icode": r["icode"].strip()}
            for i, r in enumerate(residues, start=1)
        ],
    }
    return "\n".join(lines) + "\n", residues, prep_meta


def repack_pdb(
    pdb_text: str,
    sequence: str | None = None,
    *,
    faspr: Path | str | None = None,
    renumber: str = "one_indexed",
    timeout: int = 600,
) -> dict[str, Any]:
    """Rebuild side chains on a backbone(+CB) frame.

    Args:
        pdb_text: the frame, as PDB text. Needs N, CA, C, O per residue.
        sequence: one-letter target sequence. Supply it whenever the frame's
            residue names may not be trustworthy. Length must equal the residue
            count; a mismatch raises rather than packing a shifted sequence.
        faspr: path to the FASPR executable. Default bootstraps one.
        renumber: ``one_indexed`` (default, what fpocket and every PDB tool
            expect) or ``preserve`` to write the input's own numbering back.
        timeout: seconds.

    Returns:
        dict with ``pdb`` (repacked full-atom PDB text) and a ``repack`` block
        carrying the engine, its licence, atom counts before and after, and the
        per-residue completeness check.
    """
    exe = Path(faspr) if faspr else ensure_faspr()
    prepared, residues, prep_meta = _prepare_input(pdb_text, sequence)
    n_in = len(_parse_atoms(prepared))

    with tempfile.TemporaryDirectory() as td:
        inp = Path(td) / "in.pdb"
        out = Path(td) / "out.pdb"
        inp.write_text(prepared)
        cmd = [str(exe), "-i", str(inp), "-o", str(out)]
        if sequence is not None:
            sf = Path(td) / "seq.txt"
            sf.write_text("".join(sequence.split()).upper() + "\n")
            cmd += ["-s", str(sf)]
        r = subprocess.run(
            cmd, capture_output=True, text=True, timeout=timeout,
            cwd=str(exe.parent), check=False,
        )
        if r.returncode != 0 or not out.exists():
            raise RepackError(
                f"FASPR failed (rc={r.returncode}).\nstdout: {r.stdout[-800:]}\n"
                f"stderr: {r.stderr[-800:]}"
            )
        packed = out.read_text()

    packed_res = _residues(_parse_atoms(packed))
    if len(packed_res) != len(residues):
        raise RepackError(
            f"FASPR returned {len(packed_res)} residues for {len(residues)} in — refusing "
            "to hand a silently truncated structure downstream"
        )

    incomplete = []
    for i, res in enumerate(packed_res, start=1):
        want = set(BACKBONE) | set(SIDECHAIN_ATOMS.get(res["resname"], []))
        got = set(res["atoms"])
        gap = sorted(want - got)
        if gap:
            incomplete.append({"index": i, "resname": res["resname"], "missing": gap})

    if renumber == "preserve":
        packed = _restore_numbering(packed_res, prep_meta["original_numbering"])
    elif renumber != "one_indexed":
        raise RepackError(f"renumber must be 'one_indexed' or 'preserve' (got {renumber!r})")

    n_out = len(_parse_atoms(packed))
    return {
        "pdb": packed,
        "repack": {
            "engine": "FASPR",
            "engine_version": "20200309",
            "method": "backbone-dependent (Dunbrack 2010) rotamer packing, discrete search",
            "licence": FASPR_LICENCE,
            "minimisation_performed": False,
            "force_field_optimisation": None,
            "n_residues": len(packed_res),
            "atoms_in": n_in,
            "atoms_out": n_out,
            "atoms_per_residue_in": round(n_in / len(residues), 2),
            "atoms_per_residue_out": round(n_out / len(packed_res), 2),
            "sequence_supplied": sequence is not None,
            "residues_with_missing_sidechain_atoms": incomplete,
            "complete": not incomplete,
            "numbering": renumber,
            "caveat": (
                "Rebuilt side chains are a rotamer-library prediction conditioned on this "
                "backbone, not a measurement. They carry the packer's own error "
                "(see SKILL.md for the measured chi1 / RMSD numbers on this pipeline). "
                "Buried side chains are far better predicted than surface ones, and a "
                "pocket lining is mostly buried — but any pocket volume computed off a "
                "repacked frame inherits this error and must be reported as such."
            ),
        },
    }


def _restore_numbering(packed_res: list[dict], numbering: list[dict]) -> str:
    lines = []
    serial = 1
    for res, num in zip(packed_res, numbering, strict=True):
        for name, a in res["atoms"].items():
            lines.append(
                _fmt_atom(serial, name, res["resname"], num["chain"] or "A",
                          int(num["resseq"]), a["x"], a["y"], a["z"],
                          a["element"] or name[0])
            )
            serial += 1
    lines += ["TER", "END"]
    return "\n".join(lines) + "\n"


def repack_frames(
    frames: list[str] | list[Path],
    sequence: str | None = None,
    **kwargs: Any,
) -> list[dict[str, Any]]:
    """Repack a list of frames (PDB text or paths). One dict per frame, in order.

    A frame that fails is returned as ``{"error": ...}`` rather than aborting the
    ensemble — one unphysical frame should not cost you the other thirty-one.
    """
    exe = kwargs.pop("faspr", None) or ensure_faspr()
    out = []
    for i, f in enumerate(frames):
        text = Path(f).read_text() if (isinstance(f, Path) or (isinstance(f, str) and "\n" not in f)) else f
        try:
            res = repack_pdb(text, sequence, faspr=exe, **kwargs)
        except (RepackError, subprocess.SubprocessError) as exc:
            res = {"pdb": None, "error": str(exc), "repack": {"engine": "FASPR", "complete": False}}
        res["frame_index"] = i
        out.append(res)
    return out


# --------------------------------------------------------------------------
# The step's own error, measured — not asserted
# --------------------------------------------------------------------------
def sidechain_rmsd(
    predicted_pdb: str,
    reference_pdb: str,
    *,
    subset_residues: list[int] | None = None,
    buried_only: bool = False,
    burial_cutoff: int = 16,
) -> dict[str, Any]:
    """Heavy-atom side-chain RMSD of a repacked structure against ground truth.

    Residues are paired **by order**, and the pairing is refused if the residue
    names disagree — a mis-paired RMSD is worse than no RMSD.

    Symmetry-equivalent atom pairs (PHE/TYR ring, ASP/GLU carboxylate, ARG
    guanidinium, VAL/LEU methyls) are resolved to the better assignment.

    Args:
        subset_residues: 1-based indices to restrict to (e.g. a pocket lining).
        buried_only: restrict to residues with >= ``burial_cutoff`` CB neighbours
            within 10 A. Surface side chains are genuinely mobile and their RMSD
            is not a packing error; pocket linings are buried.
    """
    import math

    pred = _residues(_parse_atoms(predicted_pdb))
    ref = _residues(_parse_atoms(reference_pdb))
    if len(pred) != len(ref):
        return {"error": f"residue count mismatch: predicted {len(pred)} vs reference {len(ref)}",
                "rmsd_a": None}

    mismatch = [i for i, (p, q) in enumerate(zip(pred, ref, strict=True), start=1)
                if p["resname"] != q["resname"]]
    if mismatch:
        return {"error": f"residue name mismatch at 1-based indices {mismatch[:10]} "
                         f"({len(mismatch)} total) — refusing to report an RMSD",
                "rmsd_a": None}

    # Burial from the reference CB/CA cloud.
    centres = []
    for r in ref:
        a = r["atoms"].get("CB") or r["atoms"].get("CA")
        centres.append((a["x"], a["y"], a["z"]) if a else None)
    burial = []
    for i, ci in enumerate(centres):
        if ci is None:
            burial.append(0)
            continue
        n = sum(
            1 for j, cj in enumerate(centres)
            if j != i and cj is not None
            and (ci[0] - cj[0]) ** 2 + (ci[1] - cj[1]) ** 2 + (ci[2] - cj[2]) ** 2 < 100.0
        )
        burial.append(n)

    sq_total, n_total = 0.0, 0
    per_residue = []
    chi1_ok, chi1_n = 0, 0
    keep = set(subset_residues) if subset_residues else None

    for idx, (p, q) in enumerate(zip(pred, ref, strict=True), start=1):
        if keep is not None and idx not in keep:
            continue
        if buried_only and burial[idx - 1] < burial_cutoff:
            continue
        rn = q["resname"]
        names = SIDECHAIN_ATOMS.get(rn, [])
        # Beyond-CB only: CB is determined by the backbone, so including it
        # flatters the packer. ALA and GLY therefore contribute nothing.
        names = [n for n in names if n != "CB"]
        if not names:
            continue
        if not all(n in p["atoms"] and n in q["atoms"] for n in names):
            continue

        def sqd(mapping: dict[str, str]) -> float:
            s = 0.0
            for n in names:
                a, b = p["atoms"][mapping.get(n, n)], q["atoms"][n]
                s += (a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2
            return s

        best = sqd({})
        for a, b in SYMMETRIC_SWAPS.get(rn, []):
            if a in p["atoms"] and b in p["atoms"]:
                cand = sqd({a: b, b: a})
                best = min(best, cand)
        sq_total += best
        n_total += len(names)
        per_residue.append({"index": idx, "resname": rn,
                            "rmsd_a": round(math.sqrt(best / len(names)), 3),
                            "cb_neighbours_10a": burial[idx - 1]})

        # chi1: correct if the first side-chain dihedral is within 40 deg.
        d = _chi1_delta(p, q, rn)
        if d is not None:
            chi1_n += 1
            chi1_ok += int(d <= 40.0)

    return {
        "rmsd_a": round(math.sqrt(sq_total / n_total), 3) if n_total else None,
        "n_residues_scored": len(per_residue),
        "n_atoms_scored": n_total,
        "chi1_within_40deg_fraction": round(chi1_ok / chi1_n, 3) if chi1_n else None,
        "n_chi1_scored": chi1_n,
        "buried_only": buried_only,
        "burial_cutoff_cb_neighbours_10a": burial_cutoff if buried_only else None,
        "excludes_cb": True,
        "symmetry_corrected": True,
        "per_residue": per_residue,
    }


_CHI1_ATOM = {
    "ARG": "CG", "ASN": "CG", "ASP": "CG", "CYS": "SG", "GLN": "CG", "GLU": "CG",
    "HIS": "CG", "ILE": "CG1", "LEU": "CG", "LYS": "CG", "MET": "CG", "PHE": "CG",
    "PRO": "CG", "SER": "OG", "THR": "OG1", "TRP": "CG", "TYR": "CG", "VAL": "CG1",
}


def _chi1_delta(p: dict, q: dict, resname: str) -> float | None:
    import math

    g = _CHI1_ATOM.get(resname)
    if g is None:
        return None
    need = ("N", "CA", "CB", g)
    if not all(n in p["atoms"] and n in q["atoms"] for n in need):
        return None

    def dihedral(res: dict) -> float:
        pts = [(res["atoms"][n]["x"], res["atoms"][n]["y"], res["atoms"][n]["z"]) for n in need]
        b0 = [pts[0][i] - pts[1][i] for i in range(3)]
        b1 = [pts[2][i] - pts[1][i] for i in range(3)]
        b2 = [pts[3][i] - pts[2][i] for i in range(3)]
        n1 = _cross(b0, b1)
        n2 = _cross(b1, b2)
        m1 = _cross(n1, [c / (_norm(b1) or 1e-9) for c in b1])
        x = sum(n1[i] * n2[i] for i in range(3))
        y = sum(m1[i] * n2[i] for i in range(3))
        return math.degrees(math.atan2(y, x))

    d = abs(dihedral(p) - dihedral(q)) % 360.0
    return min(d, 360.0 - d)


def _cross(a: list[float], b: list[float]) -> list[float]:
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]


def _norm(a: list[float]) -> float:
    return sum(c * c for c in a) ** 0.5


def strip_to_backbone_cb(pdb_text: str, *, zero_index: bool = True) -> str:
    """Reduce a full-atom structure to N/CA/C/O/CB, 0-indexed, B-factors 0.00.

    This is the *inverse* of the problem — it manufactures a BioEmu-shaped frame
    from a structure whose real side chains are known, which is how the repacker's
    own error is measured without needing a held-out set.
    """
    res_list = _residues(_parse_atoms(pdb_text))
    lines, serial = [], 1
    for i, res in enumerate(res_list):
        rn = res["resname"]
        for name in (*BACKBONE, "CB"):
            a = res["atoms"].get(name)
            if a is None:
                continue
            if name == "CB" and rn == "GLY":
                continue
            lines.append(_fmt_atom(serial, name, rn, "A", i if zero_index else i + 1,
                                   a["x"], a["y"], a["z"], name[0], 0.0))
            serial += 1
    lines += ["TER", "END"]
    return "\n".join(lines) + "\n"
