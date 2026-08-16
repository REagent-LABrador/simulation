"""Needleman-Wunsch with BLOSUM62, stdlib+numpy, used only to map one
receptor's contact residues onto its paralogue's numbering.

Written out rather than shelled to a package because the mapping is
load-bearing: it is what lets a footprint measured on a deposited complex be
compared against a prediction for a receptor with no structure at all. A wrong
alignment silently fabricates agreement, so the caller is given the identity
and the per-residue map to inspect.
"""
from __future__ import annotations

import numpy as np

_B62_RAW = """
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62 = {}
for _line in _B62_RAW.strip().splitlines():
    _p = _line.split()
    for _j, _v in enumerate(_p[1:]):
        _B62[(_p[0], _ORDER[_j])] = int(_v)


def _score(a: str, b: str) -> int:
    return _B62.get((a, b), -4)


def needleman_wunsch(s1: str, s2: str, gap_open: int = -11, gap_extend: int = -1):
    """Affine-gap global alignment. Returns (aligned1, aligned2, identity)."""
    n, m = len(s1), len(s2)
    NEG = -10 ** 9
    M = np.full((n + 1, m + 1), NEG, float)
    Ix = np.full((n + 1, m + 1), NEG, float)   # gap in s2
    Iy = np.full((n + 1, m + 1), NEG, float)   # gap in s1
    ptr = np.zeros((n + 1, m + 1, 3), np.int8)
    M[0, 0] = 0
    for i in range(1, n + 1):
        Ix[i, 0] = gap_open + gap_extend * (i - 1)
    for j in range(1, m + 1):
        Iy[0, j] = gap_open + gap_extend * (j - 1)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sc = _score(s1[i - 1], s2[j - 1])
            opts = (M[i - 1, j - 1], Ix[i - 1, j - 1], Iy[i - 1, j - 1])
            k = int(np.argmax(opts)); M[i, j] = opts[k] + sc; ptr[i, j, 0] = k
            opts = (M[i - 1, j] + gap_open, Ix[i - 1, j] + gap_extend, Iy[i - 1, j] + gap_open)
            k = int(np.argmax(opts)); Ix[i, j] = opts[k]; ptr[i, j, 1] = k
            opts = (M[i, j - 1] + gap_open, Ix[i, j - 1] + gap_open, Iy[i, j - 1] + gap_extend)
            k = int(np.argmax(opts)); Iy[i, j] = opts[k]; ptr[i, j, 2] = k
    i, j = n, m
    state = int(np.argmax((M[n, m], Ix[n, m], Iy[n, m])))
    a1, a2 = [], []
    while i > 0 or j > 0:
        if state == 0 and i > 0 and j > 0:
            nxt = ptr[i, j, 0]; a1.append(s1[i - 1]); a2.append(s2[j - 1]); i -= 1; j -= 1
        elif state == 1 and i > 0:
            nxt = ptr[i, j, 1]; a1.append(s1[i - 1]); a2.append("-"); i -= 1
        elif state == 2 and j > 0:
            nxt = ptr[i, j, 2]; a1.append("-"); a2.append(s2[j - 1]); j -= 1
        elif i > 0:
            a1.append(s1[i - 1]); a2.append("-"); i -= 1; nxt = 1
        else:
            a1.append("-"); a2.append(s2[j - 1]); j -= 1; nxt = 2
        state = int(nxt)
    a1, a2 = "".join(reversed(a1)), "".join(reversed(a2))
    pairs = [(x, y) for x, y in zip(a1, a2) if x != "-" and y != "-"]
    ident = sum(1 for x, y in pairs if x == y) / len(pairs) if pairs else 0.0
    return a1, a2, ident


def residue_map(s1: str, start1: int, s2: str, start2: int):
    """{numbering-in-1: numbering-in-2} for aligned, non-gap positions."""
    a1, a2, ident = needleman_wunsch(s1, s2)
    i, j, out = start1 - 1, start2 - 1, {}
    for x, y in zip(a1, a2):
        if x != "-":
            i += 1
        if y != "-":
            j += 1
        if x != "-" and y != "-":
            out[i] = j
    return out, ident, (a1, a2)
