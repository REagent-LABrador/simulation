"""Footprint transfer — the one gate signal that is not the model's own opinion.

A predicted interface is a set of residues on the target. So is the interface in
a deposited complex of the *same target* with a *homologous* partner. Comparing
them asks a question the confidence metrics cannot: does the model put the new
partner where the family's receptors actually bind?

Both sides are converted to UniProt numbering first, because the two structures
do not share a numbering scheme and a Jaccard over mismatched numbering is a
number with no meaning.
"""
from __future__ import annotations

import os
from typing import Sequence

import gemmi

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, "data")

_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U", "PYL": "O",
}


def uniprot_seq(acc: str) -> str:
    with open(os.path.join(_DATA, f"{acc}.fasta")) as fh:
        return "".join(l.strip() for l in fh if not l.startswith(">"))


def detect_offset(structure_path: str, chain_name: str, acc: str,
                  window: int = 25) -> int:
    """offset such that uniprot_number = auth_seqid + offset.

    Found by string-matching a run of `window` consecutive observed residues
    against the UniProt sequence, not assumed. Raises if the match is not
    unique — a silently wrong offset would fabricate a footprint overlap.
    """
    st = gemmi.read_structure(structure_path)
    st.setup_entities()
    ch = next(c for c in st[0] if c.name == chain_name)
    obs = [(r.seqid.num, _AA3TO1.get(r.name)) for r in ch if _AA3TO1.get(r.name)]
    seq = uniprot_seq(acc)
    for i in range(len(obs) - window):
        chunk = obs[i:i + window]
        nums = [n for n, _ in chunk]
        if nums != list(range(nums[0], nums[0] + window)):
            continue                      # gap; try the next window
        frag = "".join(a for _, a in chunk)
        hits = [j for j in range(len(seq) - window + 1) if seq[j:j + window] == frag]
        if len(hits) == 1:
            return (hits[0] + 1) - nums[0]
    raise ValueError(f"no unique {window}-mer match for {structure_path}:{chain_name} "
                     f"against {acc}; refusing to guess an offset")


def to_uniprot(labels: Sequence[str], chains: Sequence[str] | None,
               offset: int) -> set[int]:
    """Residue labels 'CHAIN:NUM:NAME' -> a set of UniProt numbers.

    Chain identity is dropped on purpose: for a homotrimeric ligand the same
    epitope exists on three chains and seeds may pick a different protomer.
    """
    out = set()
    for lab in labels:
        c, num, _ = lab.split(":", 2)
        if chains is None or c in chains:
            out.add(int(num) + offset)
    return out


def jaccard(a: set, b: set):
    if not a and not b:
        return None
    return len(a & b) / len(a | b)


def coverage(pred: set, ref: set):
    """Fraction of the reference footprint the prediction recovers."""
    if not ref:
        return None
    return len(pred & ref) / len(ref)


def compare(pred_labels: Sequence[str], pred_construct_start: int,
            ref_labels_uniprot: set[int], label: str) -> dict:
    """Predicted footprint (construct numbering) against a UniProt-numbered ref."""
    pred = {int(l.split(":")[1]) + pred_construct_start - 1 for l in pred_labels}
    j = jaccard(pred, ref_labels_uniprot)
    return {
        "reference": label,
        "n_pred": len(pred), "n_ref": len(ref_labels_uniprot),
        "n_shared": len(pred & ref_labels_uniprot),
        "jaccard": None if j is None else round(j, 3),
        "deviation": None if j is None else round(1.0 - j, 3),
        "ref_coverage": None if not ref_labels_uniprot else
        round(coverage(pred, ref_labels_uniprot), 3),
        "pred_uniprot_residues": sorted(pred),
        "ref_uniprot_residues": sorted(ref_labels_uniprot),
    }
