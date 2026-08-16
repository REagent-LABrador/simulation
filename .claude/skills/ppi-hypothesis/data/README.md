UniProt canonical sequences for the panel, plus the four deposited assemblies
every reference number is measured on. Construct ranges live in
`../constructs.py`, one line of provenance each.

`make_references.py` re-downloads any missing assembly, so these are a cache,
not a dependency.

**Biological assembly, never the ASU.** 3K51's asymmetric unit is one TL1A
protomer with one DcR3; the trimer-groove interface this skill measures does not
exist in it. `-assembly1.cif` is the only correct input, and 3K51's assembly
carries chains `A / A-2 / A-3` for the trimer — chain naming that a script
expecting single-letter chains will silently drop.

`8DYG-assembly1.cif` is the convention control: chains A/B must return exactly
**97** CA-CA pairs under 8 Å. `make_references.py` asserts it. If that assertion
ever fails, no other number in this skill is comparable to the rest of the
station's.
