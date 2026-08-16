"""
Batched pocket scan on Modal — the computed-tractability half of the dossier.

Runs the whole ensemble in ONE invocation: every PDB entry, every clustering
value, and holo ligand-site derivation. One call, one cold start. Four separate
calls would pay four.

WHAT THIS RETURNS, in five independently-reported stages. Each stage after
fpocket is NON-FATAL and carries its own `<stage>_status` in {ok, failed,
not_run} plus a reason, following the `prank_status` pattern — a stage that
dies must cost that stage and nothing else:

  1. fpocket + PRANK per structure per clustering value (`by_clustering`).
  2. `disorder`      — metapredict disorder fraction for the target sequence.
  3. `cryptic`       — apo/holo superposition, C-alpha displacement, clash
                       attribution, free volume (`cryptic_analysis.py`).
  4. `interface`     — pocket vs partner epitope, orthosteric / allosteric /
                       destabiliser (`interface_analysis.py`).
  5. `mdpocket`      — THE SITE FIXED BY CONSTRUCTION. One grid definition
                       applied to every superposed structure, replacing
                       post-hoc residue-number matching.

STILL NOT IMPLEMENTED — do not assume the output contains it:
  * site transfer from a structural neighbour. Requires a residue-numbering
    equivalence policy; 6OIM and 4OBE happen to share numbering, which is not
    general.

Pooled ensemble volume/druggability from the fpocket path may still span
DIFFERENT sites; see `_pooling_caveat` in the returned `ensemble` block. The
`mdpocket` block does not have that problem, because the site there is a fixed
set of grid points rather than a per-structure match — on the five apo TNF-alpha
structures that cut the across-ensemble volume CV from ~28% to ~10% (measured
28.1% at D=1.6 against 9.9%). Both figures carry about 1 percentage point of
fpocket Monte-Carlo volume noise, so quote them to two significant figures at
most; the improvement is real, the third digit is not.

BUT `mdpocket` RETURNS TWO SITE DEFINITIONS AND ONLY ONE IS THE LIGAND SITE.
`sites.site_from_ligand` is the site the dossier asks about.
`sites.site_from_density` is the most PERSISTENT cavity, which on apo TNF-alpha
is 7.73 A away — the on-axis cavity, i.e. exactly the pocket the retracted
residue-number matcher reported as "the SPD304 site". Read
`ligand_anchored` and `distance_to_donor_ligand_centroid_a` on every site entry
before calling it the pocket. See `_mdpocket_ensemble`.

Deploy (workspace MUST be rafwiewiora):
    MODAL_PROFILE=rafwiewiora modal deploy modal_app.py

Call from the eve tool handler:
    import modal
    fn = modal.Function.lookup("druggability-pocket-scan", "pocket_scan")
    result = fn.remote(pdb_ids=["6OIM", "4OBE"], ligand_codes=["MOV"])

Nothing here runs on the Anthropic sandbox or on a laptop; the agent only
decides to call it.
"""

import json
import os
import shutil
import subprocess
from collections.abc import Sequence
import sys
import urllib.request
from pathlib import Path

# `ligand_filter` lives in the sibling structure-select skill, not here. It has
# to be IMPORTABLE AT DEPLOY TIME for `add_local_python_source` below to find
# and ship it, so its directory goes on sys.path before the image is built.
# Nothing at runtime depends on this line: in the container the module is
# mounted at the top level.
sys.path.insert(
    0, str(Path(__file__).resolve().parent.parent / "structure-select")
)

import modal

# fpocket from conda-forge; proto-tools for the in-process CPU tools
# (vina-docking, foldseek-*, pyrosetta-*) that have no Proto Modal app.
# gemmi (0.7.5) arrives as a proto-tools dependency and is REQUIRED, not
# optional — the mmCIF parse is the whole input path now. If proto-tools ever
# stops depending on it this image needs an explicit `gemmi` pin.
P2RANK_VERSION = "2.5.1"
P2RANK_HOME = f"/opt/p2rank_{P2RANK_VERSION}"

image = (
    modal.Image.micromamba(python_version="3.12")
    # fpocket the conda-forge package ships mdpocket, tpocket and dpocket
    # alongside fpocket, all on PATH. mdpocket is what makes the site
    # fixed-by-construction measurement possible; there is no separate install.
    .micromamba_install("fpocket", channels=["conda-forge"])
    # JDK 17 is a hard floor for P2Rank 2.5.1 — Java 11 dies with
    # UnsupportedClassVersionError (class file v61).
    #
    # build-essential is for metapredict, NOT for anything Java. metapredict
    # publishes no Linux wheel, so pip builds it from the sdist and the build
    # needs a C compiler; without gcc the pip_install below fails and, because
    # it is its own layer, takes the whole deploy down.
    .apt_install("git", "curl", "openjdk-17-jre-headless", "build-essential")
    .run_commands(
        f"curl -sL -o /tmp/p2rank.tar.gz https://github.com/rdk/p2rank/releases/"
        f"download/{P2RANK_VERSION}/p2rank_{P2RANK_VERSION}.tar.gz",
        "tar xzf /tmp/p2rank.tar.gz -C /opt",
        f"{P2RANK_HOME}/prank --help > /dev/null 2>&1 || true",
    )
    .pip_install(
        # gemmi is a HARD requirement — mmCIF is the only structure format read.
        # It also arrives transitively via proto-tools, but pinning it directly
        # means a proto-tools dependency change cannot fail every structure at
        # stage "prepare".
        "gemmi>=0.7",
        # numpy is likewise transitive but load-bearing: the Kabsch
        # superposition, the .dx grid handling, cryptic_analysis and
        # interface_analysis are all numpy and nothing else.
        "numpy>=1.26",
        # PINNED TO A COMMIT, deliberately. Unpinned, this line resolved to
        # whatever the default branch happened to be at build time, so a
        # rebuild could silently swap the foldseek/vina/pyrosetta code that
        # every computed-tractability number was validated against — and an
        # image rebuild is not an event anyone watches.
        #
        # This SHA was NOT taken from HEAD. It is the commit the already-built
        # image is actually running, read out of the CACHED image with
        #   importlib.metadata.distribution("proto-tools")
        #       .read_text("direct_url.json")
        # and independently corroborated by the local `druggability`
        # micromamba env, which reports the same commit. (It happens to
        # coincide with the default branch tip on 2026-08-15, so this pin
        # changes nothing about what installs — only about what a FUTURE
        # rebuild installs.)
        #
        # Bumping it is a deliberate, re-validated act: change the SHA only
        # together with a re-run of the pocket-scan fixtures.
        "proto-tools[mcp] @ git+https://github.com/evo-design/proto-tools.git"
        "@1e0bd8f5a8f4525eb5e5c736cbf25c1366929e73",
    )
    # CPU TORCH, ITS OWN LAYER, BEFORE metapredict. A bare `pip install torch`
    # resolves to the CUDA build and drags in >500 MB of nvidia-* wheels that
    # this function can never use — it has no GPU. The pytorch CPU index serves
    # a torch that is ~10x smaller. metapredict declares `torch` as a plain
    # dependency, so it MUST be satisfied first from this index or pip will
    # pull the CUDA one while resolving metapredict.
    .pip_install("torch", index_url="https://download.pytorch.org/whl/cpu")
    # metapredict is the right disorder tool on merit — MIT (IUPred3's licence
    # forbids commercial use), CPU-only, and it separates MYC P01106 (~0.83
    # disordered) from CDK2 P24941 (~0.00). An earlier version of this file
    # dropped it because a bare `metapredict` pin failed to build a wheel and
    # broke the deploy. The cause was the missing compiler above, not the
    # package: it has no Linux wheel and builds from sdist. Pinned to 3.0.2 so a
    # new sdist upstream cannot silently change the numbers or the build needs.
    .pip_install("metapredict==3.0.2")
    # Modal 1.x does not automount sibling modules. These three are imported
    # INSIDE the function body, but they still have to exist in the container.
    # `ligand_filter` replaces the comp_id denylist and the heavy-atom floor
    # that used to decide holo vs apo here. It is stdlib-only by design — this
    # image has no RDKit and a toolkit dependency would make the verdict vary
    # by environment.
    .add_local_python_source(
        "cryptic_analysis", "disorder", "interface_analysis", "ligand_filter"
    )
)

app = modal.App("druggability-pocket-scan", image=image)

# Sweep, never pin — see the skill's failure modes. -D 1.6 gives a FALSE
# NEGATIVE (0.002) on TNF-alpha's co-crystallised site; -D 2.4 fuses KRAS's
# nucleotide and switch-II sites into one meaningless 1540 A^3 mega-pocket.
D_VALUES = (1.6, 2.4)

# ---------------------------------------------------------------------------
# HOLO vs APO IS A CHEMISTRY QUESTION AND IT IS ANSWERED BY CHEMISTRY.
#
# What used to be here: DRUGLIKE_MIN_HEAVY_ATOMS = 18, a NON_LIGANDS denylist of
# buffers and ions, and a COFACTORS denylist of nucleotides, sugars and lipids.
# All three are DELETED. A size threshold cannot work and no list can be
# complete, and both halves of that are measured:
#
#   * ADP has 27 heavy atoms. So does `A1IPJ`, the genuine inhibitor in 9GU4.
#     No floor separates them, ever.
#   * Identity filtering gave 16 holo / 8 apo on NLRP3 where a naive size window
#     gave 19 / 5 — three false holo entries, a 19% overstatement. And `CPS`
#     (CHAPS, 615 Da) was simply missing from the list and sailed through.
#   * The same shape produced wrong answers on CD20 (sterol tails), KRAS
#     neighbours (2UK: purine + ribose + phosphate) and IL-17A neighbours
#     (L44's 21-carbon chain).
#
# `ligand_filter.classify_record` reads the actual structure of the component:
# 259/262 on ground truth, 61/70 on a blind held-out set with ZERO false
# positives, and it reproduces the deleted COFACTORS set without having been
# shown it. Every remaining error is conservative — it calls a drug a cofactor
# rather than the reverse.
#
# Two of its behaviours are load-bearing here and must not be collapsed:
#   1. `unknown` IS NOT `apo`. An unclassifiable component leaves the entry's
#      state UNDETERMINED, which is a third tier, not a quiet apo.
#   2. A lookup failure is not a CCD miss. A component whose record could not be
#      retrieved carries `lookup_failed` and lands in `undetermined`, so a flaky
#      network cannot render a holo structure apo.

# Waters are the one thing still dropped by name, and only because a structure
# has hundreds of them: one payload row per water copy is a payload problem, not
# a classification problem. Everything else goes through the classifier.
WATER_COMP_IDS = frozenset({"HOH", "DOD", "WAT"})

# RCSB's public Chemical Component REST endpoint. The classifier's own default
# source shells out to the `paperclip` binary, which is NOT in this image; this
# is the same data (type, name, formula, formula_weight and — the field that
# matters — SMILES) over the network path this function already uses for
# structures. Verified: MOV druglike, GDP cofactor, ADP cofactor, CPS
# lipid_or_detergent, 307 druglike + promiscuity_advisory, A1JPS druglike,
# N5S (the 5QQE fragment) druglike, Y01/CLR/PC1/L44 lipid_or_detergent.
RCSB_CHEMCOMP_URL = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"

# How many pockets come back per structure per clustering value. fpocket routinely
# detects >100 (134 on an IRAK4 assembly at D=1.6), the great majority of which
# are sub-100 A^3 surface dimples, and returning all of them for every structure
# at every D makes the payload unusable. The cut is by fpocket rank, the selected
# site pocket is always included whatever its rank, and `pockets_omitted_summary`
# bounds the volume, druggability and site overlap of everything left out — so
# the truncation is checkable rather than merely declared. PROPOSED, NOT
# CALIBRATED: it is a payload-size choice, not a statement about pockets.
MAX_POCKETS_RETURNED = 30

# How many pockets get the full interface classification. Lower than
# MAX_POCKETS_RETURNED because enclosure casts 512 rays per probe point per
# chain and is the expensive step, whereas returning a pocket is nearly free.
# The selected site pocket is always classified whatever its rank.
MAX_POCKETS_CLASSIFIED = 10

# Legal single-character PDB chain identifiers, in the order they get handed
# out when an mmCIF chain name will not fit column 22.
_PDB_CHAIN_POOL = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"

# Polymer-chain count above which whole-assembly scoring stops being a
# measurement of a binding site. PROPOSED, NOT CALIBRATED: set above every
# assembly in the calibration set (NLRP3's octamer is the largest at 8) and far
# below the 60-mers that produced the failure. It FLAGS and never filters.
LARGE_ASSEMBLY_CHAINS = 12


def _fetch(pdb_id: str, dest: Path) -> Path:
    """Fetch the mmCIF. Always the mmCIF, never the legacy PDB.

    Legacy PDB is not a subset of the truth, it is a lossy encoding of it, and
    three of its losses bite this module directly:

      * the chemical component ID has three columns. The PDB ran out of 3-char
        codes, so 2024+ depositions carry five-character comp_ids — 9SQX's
        ligand is `A1JPS`. Parsed out of columns 18-20 that reads as `A1J`,
        the ligand is never found, the ligand site comes back empty, and the
        run silently degrades to "most druggable pocket anywhere", which on
        9SQX picks a 3606 A^3 merge artifact. That was a real wrong answer on
        IL-17A, not a hypothetical.
      * chain IDs have one column, and >99999 atoms cannot be numbered.
      * RCSB no longer issues it at all for newer entries — verified, 9SQX.pdb
        is HTTP 404 while 9SQX.cif is 200. Recent structures are exactly where
        new chemistry lives, so a pdb-first fetcher fails on the most
        interesting targets.

    So mmCIF is the single source of truth for the whole per-structure pass,
    and the only PDB file in play is the one written for fpocket, which accepts
    nothing else.

    THE BIOLOGICAL ASSEMBLY, NOT THE ASYMMETRIC UNIT.

    The asymmetric unit is a crystallographic artifact. It may hold several
    copies of the biological unit, or only part of one, and both errors are
    silent and severe for pocket detection. Measured on our own runs:

      * 9SQX's preferred assembly is a DIMER, but its ASU holds two of them.
        Scoring all four chains fused them and produced a 3606 A^3 "pocket"
        that no molecule occupies.
      * 5HI3, the IL-17A macrocycle structure, has a HEPTAMER as its preferred
        assembly while the small-molecule site lies in the dimer groove.
      * 8USS's preferred assembly is a MONOMER, so a site that spans the
        IL-17A dimer interface is only half present — which is the likely
        reason it recovered 15 site residues at Jaccard 0.29 while 8DYG
        managed 0.69.

    So fetch `<ID>-assembly1.cif` first and record what was used. Falling back
    to the ASU is allowed but must be visible in the output, never silent.
    """
    cif = dest / f"{pdb_id}.cif"
    if cif.exists():
        return cif

    # Preferred biological assembly first.
    try:
        with urllib.request.urlopen(  # noqa: S310
            f"https://files.rcsb.org/download/{pdb_id}-assembly1.cif", timeout=60
        ) as r:
            cif.write_bytes(r.read())
        (dest / f"{pdb_id}.source").write_text("assembly1")
        return cif
    except urllib.error.HTTPError:
        pass  # fall through to the asymmetric unit

    try:
        with urllib.request.urlopen(  # noqa: S310
            f"https://files.rcsb.org/download/{pdb_id}.cif", timeout=60
        ) as r:
            cif.write_bytes(r.read())
        (dest / f"{pdb_id}.source").write_text("asymmetric_unit")
    except urllib.error.HTTPError as exc:
        # Never let the raw HTTPError escape: it holds a BufferedReader, which
        # cannot pickle, so Modal reports an opaque SerializationError instead
        # of the 404 that actually happened.
        raise RuntimeError(f"{pdb_id}: no CIF at RCSB (HTTP {exc.code})") from None
    return cif


def _load(cif: Path) -> tuple[object, list[str], dict[str, str]]:
    """Parse the mmCIF ONCE into the object everything else is derived from.

    Returns (structure, missing_residues, chain_renaming).

    Cleaning happens here, not in each consumer, so that the ligand inventory,
    the ligand contact site and the file handed to fpocket can never disagree
    about which atoms exist:
      * hydrogens and deuteriums dropped (element symbol, not a name guess —
        mmCIF carries `type_symbol`, so HG the mercury is never a hydrogen);
      * altloc kept at blank or A, then blanked, so nothing downstream sees a
        partly occupied "A" as a distinct conformer.

    Chain names are forced to be single-character here as well, BEFORE anything
    reads them. This is the one place the mmCIF -> PDB round trip can silently
    corrupt a run, and it is worse than "gemmi renames the chain": measured on a
    2-character chain name, `make_pdb_string()` writes BOTH characters, into
    columns 21-22, eating the space after resName —

        ATOM      1  N   THRAA  44      -1.396  21.115   8.728

    so column 22 (the chain ID everything downstream slices) ends up holding the
    SECOND character. Two chains AA and BA would both come back as "A". fpocket
    reports pocket residues in the chain IDs of the file it was given, so those
    IDs would not match the site residues derived from the CIF, every Jaccard
    would be 0.0, and the module would report "no pocket overlapped the ligand
    site" on a structure where one plainly does. Renaming up front, once, with
    the map returned, makes the two sides consistent by construction; `_prep`
    re-reads the written file and asserts it.
    """
    import gemmi

    doc = gemmi.cif.read(str(cif))
    block = doc.sole_block()
    st = gemmi.make_structure_from_block(block)
    # `make_structure_from_block` emits one Chain per SUBCHAIN, so an entry with
    # polymer + ligand + waters under auth chain A comes back as three separate
    # Chain objects all called "A" — verified on 8DYG, six chains for a dimer.
    # `read_structure()` merges them and this must too, or the rename map below
    # would key several chains on one name.
    st.merge_chain_parts()
    st.setup_entities()
    st.remove_hydrogens()

    for chain in st[0]:
        for res in chain:
            for i in range(len(res) - 1, -1, -1):
                if res[i].altloc not in ("\x00", "A"):
                    del res[i]
                else:
                    res[i].altloc = "\x00"

    # Single-character chain IDs, keeping every name that already fits.
    used = {c.name for c in st[0] if len(c.name) == 1 and c.name != " "}
    pool = [c for c in _PDB_CHAIN_POOL if c not in used]
    renamed: dict[str, str] = {}
    for chain in st[0]:
        if (len(chain.name) == 1 and chain.name != " ") or chain.name in renamed:
            continue
        if not pool:
            raise RuntimeError(
                f"{cif.stem}: more chains than PDB chain IDs "
                f"({len(_PDB_CHAIN_POOL)}); cannot write an fpocket input "
                "without losing chain identity. NOTE THIS IS AN IMPLEMENTATION "
                "LIMIT, NOT A JUDGEMENT ABOUT THE STRUCTURE: an assembly just "
                "under it is not thereby suitable for whole-assembly scoring. "
                "1JH5, a 60-mer, fits and returned 378 pockets with the "
                "selected one 60.28 A from the protein centre. See "
                "LARGE_ASSEMBLY_CHAINS, which flags both sides of this "
                "boundary the same way."
            )
        renamed[chain.name] = pool.pop(0)
    for old, new in renamed.items():
        st.rename_chain(old, new)

    # `_pdbx_unobs_or_zero_occ_residues` is the mmCIF category RCSB generates
    # REMARK 465 from — verified row-for-row against 6OIM's 16 REMARK 465 lines.
    missing: list[str] = []
    tab = block.find(
        "_pdbx_unobs_or_zero_occ_residues.",
        ["auth_asym_id", "auth_seq_id", "?PDB_model_num"],
    )
    for row in tab:
        if row.has(2) and row.str(2) not in ("", ".", "?", "1"):
            continue
        ch = row.str(0)
        missing.append(f"{renamed.get(ch, ch)}/{row.str(1)}")
    return st, missing, renamed


def _prep(
    st, dest: Path, stem: str, chains: list[str] | None,
    drop_chains: Sequence[str] = (),
) -> tuple[Path, list[str], list[str]]:
    """The fpocket input: polymer only, from the same object as everything else.

    ... AND NOT THE POLYMER THAT IS THE LIGAND. `het_flag == 'A'` keeps EVERY
    polymer, so rule 4's strip-every-ligand-before-scoring requirement has been
    silently violated for any entry whose ligand is a peptide, a nanobody or a
    designed mini-binder. Measured on 8QFZ (TSLP + a 12-residue Bicycle
    peptide): with the peptide kept, the ligand-anchored pocket is 283.6 A^3 and
    SIX of its TEN lining residues are the peptide itself. With the peptide
    chain dropped the site does not exist at all and the whole target has ONE
    pocket, 147.8 A^3, 15.4 A away. That is a 136 A^3 swing straight across the
    discriminating band of the (retracted) rule 4a volume guide, on a decision
    this file was making by accident.

    `drop_chains` is that decision made on purpose. It is NOT a size heuristic —
    see `_classify_polymer_chains` for the four cases and the order they are
    decided in — and the third return value names what was dropped, so nothing
    is silent.

    Chain selection is per-target and deliberate: KRAS is a monomer, TNF-alpha's
    site sits on the trimer axis and disappears if you keep one chain.

    Polymer means het_flag == 'A', which is exactly what used to be selected by
    `line.startswith("ATOM")`. Modified residues (MSE and friends) are HETATM in
    both encodings and are dropped here as they always were.

    The written file's chain IDs are re-read and checked against the chains we
    think we wrote. fpocket's residue lists are only comparable to the ligand
    site because those two agree, so this is asserted rather than assumed.
    """
    import gemmi

    sel = gemmi.Structure()
    sel.name = st.name
    sel.cell = st.cell
    sel.spacegroup_hm = st.spacegroup_hm
    model = gemmi.Model("1")
    seen_chains = set()
    dropped: list[str] = []
    drop = set(drop_chains or ())
    for chain in st[0]:
        if chains and chain.name not in chains:
            continue
        if chain.name in drop:
            dropped.append(chain.name)
            continue
        keep = gemmi.Chain(chain.name)
        for res in chain:
            if res.het_flag != "A" or not len(res):
                continue
            keep.add_residue(res)
        if len(keep):
            model.add_chain(keep)
            seen_chains.add(chain.name)
    sel.add_model(model)
    sel.setup_entities()

    # Only the coordinate records, as before — no CRYST1, no headers, nothing
    # for fpocket to have an opinion about.
    lines = [
        ln
        for ln in sel.make_pdb_string().splitlines()
        if ln.startswith(("ATOM", "TER", "END"))
    ]
    out = dest / f"{stem}_prep.pdb"
    out.write_text("\n".join(lines) + "\n")

    written = {ln[21] for ln in lines if ln.startswith("ATOM") and len(ln) >= 22}
    if written != seen_chains:
        raise RuntimeError(
            f"{stem}: chain IDs changed on PDB write ({sorted(seen_chains)} -> "
            f"{sorted(written)}); fpocket residues would not map to the "
            "ligand site"
        )
    return out, sorted(seen_chains), sorted(dropped)


# A non-target polymer chain at or below this many monomers is a LIGAND, not a
# partner. PROPOSED, NOT CALIBRATED. The Bicycle peptide in 8QFZ is 12 monomers;
# the smallest thing anyone would call a binding partner in the fixture set is
# 9Q8N's LptE at 170+. Nothing in between has been measured, and the test is
# third of four for that reason — chemistry (`polymer_conjugate`) and the
# caller's own assertion both decide before size does.
POLYMER_LIGAND_MAX_MONOMERS = 50

# A non-target polymer chain contributing at least this share of the ANCHORED
# pocket's lining is that pocket's ligand whatever its size. PROPOSED, NOT
# CALIBRATED. On 8QFZ the peptide contributes 6 of 10 (0.60); an antibody Fab
# lining a cavity reported as the target's site contributes 1.00.
POLYMER_LIGAND_MIN_LINING_FRACTION = 0.25


def _chain_monomer_counts(st) -> dict[str, int]:
    """{chain -> polymer residue count}, on the same objects `_prep` selects."""
    out: dict[str, int] = {}
    for chain in st[0]:
        n = sum(1 for r in chain if r.het_flag == "A" and len(r))
        if n:
            out[chain.name] = n
    return out


def _polymer_conjugate_host_chains(holo_call: dict, context, renamed: dict) -> set[str]:
    """Chains of the polymer entities a `polymer_conjugate` is a part of.

    `ligand_filter` returns `polymer_ligand_precedent` naming the HOST entity of
    every covalently conjugated component — for 8QFZ, the entity of the bicyclic
    peptide `LFI` is bonded three times into. `_entity_poly.pdbx_strand_id` gives
    that entity's chains, which is the thing `_prep` has to drop.

    Empty when there is no context, which is the honest answer: with
    `_struct_conn` absent nothing is known to be bonded to anything, and this
    must not be read as "nothing is".
    """
    if context is None:
        return set()
    ents = {
        str(p.get("entity_id"))
        for p in (holo_call.get("polymer_ligand_precedent") or [])
        if p.get("entity_id")
    }
    out: set[str] = set()
    for eid in ents:
        ent = (getattr(context, "polymer_entities", {}) or {}).get(eid)
        for strand in (getattr(ent, "strand_ids", ()) or ()) if ent else ():
            out.add(renamed.get(strand, strand))
    return out


def _classify_polymer_chains(
    st, tgt_chains: Sequence[str], verified: bool,
    caller_chains: Sequence[str] | None, conjugate_chains: set[str],
    lining_share: dict[str, float] | None = None,
) -> dict[str, dict]:
    """Which polymer chains are the target, a partner, or the LIGAND.

    Rule 4 says strip every ligand before scoring, and a Bicycle peptide is a
    ligand. Four cases, DECIDED IN THIS ORDER, and the order is the whole point
    — chemistry and the caller's assertion both outrank size:

      target           its `_struct_ref` accession is the run's target      KEEP
      caller_asserted  named in the `chains` argument (rule 2b)             KEEP
      polymer_ligand   non-target AND (hosts a `polymer_conjugate`
                       component, OR <= 50 monomers, OR lines >= 25% of
                       the anchored pocket)                                 STRIP
      partner          anything else non-target                             KEEP

    The homo-oligomer cases are safe BY CONSTRUCTION and this is why the test is
    accession-based rather than sequence- or size-based: TNF-alpha's three
    subunits are the SAME entity with the SAME accession, so all three are
    `target` and all three are kept, and the trimer-axis site survives. KRAS
    4OBE A vs A+B is unaffected. An obligate hetero-oligomer whose site genuinely
    spans two proteins (9Q8N's LptD/LptE barrel) lands in `partner` and is kept.

    FAILS SAFE ON AN UNVERIFIED CHAIN SET. `verified` False means the entry's
    UniProt mapping could not be read, so no chain can be shown to be non-target
    and NOTHING is stripped. Stripping on an unreadable header would delete the
    target itself.
    """
    counts = _chain_monomer_counts(st)
    tgt = set(tgt_chains or ())
    caller = set(caller_chains or ())
    out: dict[str, dict] = {}
    for ch, n in sorted(counts.items()):
        if ch in tgt:
            cls, why = "target", "accession matches the run's target"
        elif ch in caller:
            cls, why = "caller_asserted", "named in the `chains` argument (rule 2b)"
        elif not verified:
            cls, why = "partner", (
                "the entry's UniProt mapping could not be read, so this chain "
                "cannot be shown to be non-target; kept"
            )
        elif ch in conjugate_chains:
            cls, why = "polymer_ligand", (
                "hosts a component ligand_filter classified `polymer_conjugate`"
            )
        elif n <= POLYMER_LIGAND_MAX_MONOMERS:
            cls, why = "polymer_ligand", (
                f"non-target polymer of {n} monomers "
                f"(<= {POLYMER_LIGAND_MAX_MONOMERS})"
            )
        elif (lining_share or {}).get(ch, 0.0) >= POLYMER_LIGAND_MIN_LINING_FRACTION:
            cls, why = "polymer_ligand", (
                f"non-target chain contributing {lining_share[ch]:.0%} of the "
                "anchored pocket's lining residues"
            )
        else:
            cls, why = "partner", "non-target, kept and its lining share reported"
        out[ch] = {
            "class": cls,
            "n_monomers": n,
            "basis": why,
            "lining_fraction_of_anchored_pocket": (lining_share or {}).get(ch),
        }
    return out


def _fpocket_once(prep: Path, work: Path, tag: str, d: float) -> list[dict]:
    """fpocket on one prepared file at one clustering value. Never raises."""
    run = work / f"{tag}_D{d}"
    run.mkdir(parents=True, exist_ok=True)
    tgt = run / prep.name
    tgt.write_text(prep.read_text())
    out_dir = run / f"{tgt.stem}_out"
    # A warm Modal container reuses /tmp; a stale _out would be parsed as this
    # run's result. Same reason as the main pass.
    shutil.rmtree(out_dir, ignore_errors=True)
    try:
        subprocess.run(  # noqa: S603
            ["fpocket", "-f", str(tgt), "-D", str(d)],  # noqa: S607
            check=False, capture_output=True,
        )
    except Exception:  # noqa: BLE001
        return []
    return _parse_pockets(out_dir)


def _anchored_pocket(pockets: list[dict], anchor_residues: list[str]):
    """The pocket that IS the anchor ligand's site, or None. Jaccard, > 0."""
    if not anchor_residues:
        return None
    best = max(
        pockets,
        key=lambda p: _jaccard(p.get("residues", []), anchor_residues) or 0.0,
        default=None,
    )
    if best is None or not (_jaccard(best.get("residues", []), anchor_residues) or 0.0):
        return None
    return best


class LigandSourceError(RuntimeError):
    """The chemical-component source cannot support a holo/apo call at all.

    RAISED, NOT RETURNED, AND DELIBERATELY NOT CAUGHT BY THE PER-STRUCTURE
    HANDLER. This is the most dangerous failure the ligand stage has, because
    its symptom is a clean-looking run:

    `ligand_filter` classifies on the component's SMILES graph. Handed records
    with no SMILES it correctly returns `unknown` for EVERY component — and
    `unknown` is not `druglike`, so every structure comes back apo (or, under
    the lookup-failure rule, `undetermined`). The payload is well-formed, every
    status says ok, and the entire ensemble is silently holo-free. Verified
    directly: MOV, GDP, ADP, CPS and `307` all return `unknown` with
    "the CCD row has no SMILES, so no chemistry test can run" when the record
    carries only type/name/formula/weight.

    WHICH SOURCES CARRY SMILES:
      * RCSB `data.rcsb.org/rest/v1/core/chemcomp/<ID>` — YES, in
        `pdbx_chem_comp_descriptor` (type SMILES_CANONICAL or SMILES). This is
        what this module uses.
      * Paperclip `pdb_v.chemcomps` — YES, in the `smiles` column. The
        classifier's own default source.
      * The RCSB CCD ligand file `files.rcsb.org/ligands/download/<ID>.cif` —
        YES.
      * THE ENTRY'S OWN mmCIF `_chem_comp` BLOCK — **NO**. It carries id, type,
        name, formula and formula_weight and nothing else. It is the obvious
        source to reach for, because the file is already parsed and on disk, and
        it is the one that does not work.

    So the check is not "did classification succeed" but "does this source
    return the field classification needs", and it is answered by looking at the
    records rather than at the verdicts.
    """


def _assert_records_carry_smiles(src, distinct: list[str]) -> None:
    """Refuse a record source that returns rows without SMILES.

    Only fires when records WERE retrieved and none of them has a SMILES string.
    A 404 (component genuinely absent from the CCD) caches as `None` and is not
    a record, so a structure whose components are all unknown to the CCD does
    not trip this; nor does a network failure, which has its own `lookup_failed`
    path and its own `undetermined` tier.
    """
    have = [
        r for r in (src.get_many(distinct) or {}).values() if isinstance(r, dict)
    ]
    if not have:
        return
    if any((r.get("smiles") or "").strip() for r in have):
        return
    raise LigandSourceError(
        f"the chemical-component source returned {len(have)} record(s) for "
        f"{', '.join(distinct[:8])} and NOT ONE carries a SMILES string. "
        "ligand_filter classifies on the SMILES graph, so every component "
        "would come back `unknown`, nothing would be `druglike`, and this run "
        "would report an entirely holo-free ensemble while looking healthy. "
        "That is a misconfigured record source, not a result. The entry's own "
        "mmCIF `_chem_comp` block is the usual cause: it has type, name, "
        "formula and weight but no SMILES. Use RCSB "
        "data.rcsb.org/rest/v1/core/chemcomp/<ID>, Paperclip "
        "pdb_v.chemcomps, or the CCD ligand file — all three carry it."
    )


def _lf_verdicts():
    """`ligand_filter` itself, for provenance in the method block."""
    import ligand_filter as LF

    return LF


def _chemcomp_source():
    """`ligand_filter.ChemCompSource` backed by RCSB instead of Paperclip.

    The classifier's default source shells out to the `paperclip` binary, which
    this image does not carry. Overriding the fetch keeps every other behaviour
    of the class — the process cache, the batching, and above all the
    `fetch_errors` bookkeeping that separates A LOOKUP THAT FAILED from A
    COMPONENT THAT IS NOT IN THE CCD. Collapsing those two is the same
    fail-open shape as reporting a credential error as "no data", and it would
    turn a flaky network into a run full of apo structures.
    """
    import ligand_filter as LF

    class _RcsbChemComps(LF.ChemCompSource):
        def _fetch_batch(self, batch: list[str], *, attempts: int = 3) -> None:
            for cid in batch:
                err = None
                doc = None
                for _ in range(max(1, attempts)):
                    try:
                        with urllib.request.urlopen(  # noqa: S310
                            RCSB_CHEMCOMP_URL.format(cid), timeout=30
                        ) as r:
                            doc = json.load(r)
                        err = None
                        break
                    except urllib.error.HTTPError as exc:
                        if exc.code == 404:
                            # A GENUINE CCD MISS. Cached as absent, no error.
                            err, doc = None, None
                            break
                        err = f"HTTP {exc.code}"
                    except Exception as exc:  # noqa: BLE001
                        err = f"{type(exc).__name__}: {exc}"
                if err:
                    self.last_error = err
                    self.fetch_errors[cid] = err
                    self._cache.setdefault(cid, None)
                    continue
                if doc is None:
                    self._cache[cid] = None
                    continue
                cc = doc.get("chem_comp") or {}
                desc = doc.get("pdbx_chem_comp_descriptor") or []
                def _pick(*types):
                    for t in types:
                        for d in desc:
                            if d.get("type") == t and d.get("descriptor"):
                                return d["descriptor"]
                    return None
                self._cache[cid] = {
                    "comp_id": cc.get("id") or cid,
                    "type": cc.get("type"),
                    "formula": cc.get("formula"),
                    "formula_weight": cc.get("formula_weight"),
                    "drugbank_id": None,
                    "inchikey": _pick("InChIKey"),
                    "smiles": _pick("SMILES_CANONICAL", "SMILES"),
                    "name": cc.get("name"),
                }

    return _RcsbChemComps()


def _ligands(st, src, context=None) -> tuple[list[dict], dict]:
    """Nonpolymer components, CLASSIFIED BY CHEMISTRY rather than by list.

    `context` is a `ligand_filter.StructureContext` for THIS entry, and it is
    what turns the covalent rules on. Without it a component that is a covalent
    constituent of a polymer ligand — the crosslinker of a bicyclic peptide, a
    warhead on a nanobody — classifies on its own SMILES graph and comes back
    `druglike`, which is the wrong modality and, worse, makes it the thing the
    site is anchored on. Measured on 8QFZ: `LFI` is `druglike` with no context
    and `polymer_conjugate` with one. BUILD IT FROM THE HEADER, never from the
    assembly file — see `_structure_context`.

    comp_id comes from the mmCIF, so it is the FULL component ID: `A1JPS`, not
    the first three characters of it — the legacy PDB truncation is a
    documented wrong answer on IL-17A.

    Returns (per-copy ligand list, entry-level holo call). The `druglike` and
    `cofactor` keys are kept so nothing downstream changes shape, but they are
    now derived from `ligand_filter`'s verdict:

        druglike  <- verdict == "druglike"
        cofactor  <- verdict == "cofactor"

    The other six verdicts (lipid_or_detergent, crystallisation_additive,
    sugar_or_glycan, ion_or_solvent, peptide_or_polymer, unknown) are neither,
    and each is reported with the reason it was given.
    """
    import ligand_filter as LF

    counts: dict[tuple[str, str, str], int] = {}
    for chain in st[0]:
        for res in chain:
            if res.het_flag != "H" or not len(res):
                continue
            if res.name.upper() in WATER_COMP_IDS:
                continue
            key = (res.name, chain.name, str(res.seqid.num))
            counts[key] = counts.get(key, 0) + len(res)
    distinct = sorted({c for c, _ch, _rs in counts})
    if distinct:
        # BEFORE any verdict is read. A source with no SMILES produces a
        # perfectly well-formed, entirely holo-free run; see LigandSourceError.
        _assert_records_carry_smiles(src, distinct)
    verdicts = (LF.classify_ligands(distinct, chemcomps=src, context=context)
                if distinct else {})
    holo = LF.holo_call(distinct, chemcomps=src, context=context) if distinct else {
        "is_holo": False, "druglike_ligands": [], "by_verdict": {},
        "unknown_ligands": [], "undetermined": [], "determined": True,
        "verdicts": {}, "flags": [],
        # Shape-stable with the classified branch: a consumer reading
        # `polymer_ligand_precedent` must not have to distinguish "no components
        # at all" from "no polymer ligand".
        "polymer_conjugates": [], "polymer_ligand_precedent": [],
        "context_applied": False,
    }
    ligs = []
    for (c, ch, rs), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        v = verdicts.get(c.upper())
        ligs.append({
            "comp_id": c,
            "chain": ch,
            "resseq": rs,
            "heavy_atoms": n,
            # Reported, not dropped: the dossier has to be able to say "apo,
            # but carrying GDP" rather than "apo" full stop.
            "cofactor": bool(v and v.verdict == "cofactor"),
            "druglike": bool(v and v.verdict == "druglike"),
            "verdict": v.verdict if v else "unknown",
            "verdict_reason": v.reason if v else "no classification attempted",
            "verdict_confidence": v.confidence if v else "none",
            "verdict_flags": list(v.flags) if v else [],
            "verdict_source": v.source if v else None,
        })
    # Ions and ordered solvent collapse to a distinct-comp_id summary: a
    # structure can carry fifty sulfates and one row each is payload, not
    # information. Nothing is silently dropped — the count is stated.
    ions = [lig for lig in ligs if lig["verdict"] == "ion_or_solvent"]
    ligs = [lig for lig in ligs if lig["verdict"] != "ion_or_solvent"]
    holo = dict(holo)
    holo["ion_or_solvent_copies"] = len(ions)
    holo["ion_or_solvent_comp_ids"] = sorted({lig["comp_id"] for lig in ions})
    holo["waters_excluded_by_name"] = sorted(WATER_COMP_IDS)
    return ligs, holo


def _ligand_site(
    st,
    comp_id: str,
    chains: list[str] | None = None,
    cutoff: float = 5.0,
) -> tuple[list[str], str | None]:
    """Residues within `cutoff` of ONE copy of the ligand — the only ground
    truth for whether the detected pocket is the pocket that matters.

    One copy, and only the chains that were actually scored. Both restrictions
    are load-bearing: 2AZ5 has two copies of ligand 307 across an A/B/C/D
    two-dimer asymmetric unit, and pooling them returns 43 residues spanning
    four chains instead of the 19-residue A/B site. Jaccard against that union
    is meaningless — a pocket found in the A/B dimer can never exceed ~0.44
    against it, so the wrong pocket wins.

    Same structure object as `_prep`, so the chain/resseq labels here are the
    ones fpocket will hand back.

    Returns (residues, "chain/resseq of the copy used").
    """
    copies: dict[tuple[str, str], list] = {}
    grid: dict[tuple[int, int, int], list] = {}
    for chain in st[0]:
        for res in chain:
            if res.het_flag == "H" and res.name == comp_id:
                copies.setdefault((chain.name, str(res.seqid.num)), []).extend(
                    (a.pos.x, a.pos.y, a.pos.z) for a in res
                )
            elif res.het_flag == "A":
                if chains and chain.name not in chains:
                    continue
                tag = f"{chain.name}/{res.seqid.num}"
                for a in res:
                    cell = (
                        int(a.pos.x // cutoff),
                        int(a.pos.y // cutoff),
                        int(a.pos.z // cutoff),
                    )
                    grid.setdefault(cell, []).append((a.pos.x, a.pos.y, a.pos.z, tag))
    c2 = cutoff * cutoff
    offsets = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]

    def contacts(lig: list) -> set:
        # Cell-hashed at the cutoff, so the 27 neighbouring cells hold every
        # atom that can possibly be within it. Same answer as all-pairs.
        hits = set()
        for lx, ly, lz in lig:
            base = (int(lx // cutoff), int(ly // cutoff), int(lz // cutoff))
            for di, dj, dk in offsets:
                for px, py, pz, tag in grid.get(
                    (base[0] + di, base[1] + dj, base[2] + dk), ()
                ):
                    if tag in hits:
                        continue
                    if (px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2 <= c2:
                        hits.add(tag)
        return hits

    # The copy best engaged by the chains we kept. For a single-copy structure
    # this is a no-op; for a multi-copy one it picks the site we can score.
    best: set = set()
    best_key: tuple[str, str] | None = None
    for key, lig in copies.items():
        hits = contacts(lig)
        if len(hits) > len(best):
            best, best_key = hits, key
    return (
        sorted(best, key=lambda s: (s.split("/")[0], int(s.split("/")[1]))),
        f"{best_key[0]}/{best_key[1]}" if best_key else None,
    )


def _ligand_contact_chains(
    st,
    comp_id: str,
    copy_key: str | None,
    cutoff: float = 5.0,
    min_fraction: float = 0.15,
) -> tuple[list[str], dict[str, int]]:
    """Chains that really line one ligand copy, counted in ATOMS not residues.

    THE UNIT MATTERS AND GETTING IT WRONG IS A REAL BUG, not a nicety. Counting
    RESIDUES within the shell and keeping any chain above 15% of the top
    contributor kept 2AZ5's chain B — a crystal contact worth 4 atoms out of 46
    — because 2 residues against 11 clears 15% while 4 atoms against 46 (8.7%)
    does not. The donor then had THREE chains, no 2-of-3 mapping onto the
    TNF-alpha trimer existed, the best superposition was 17.3 A, the ligand was
    never transferred and the whole ligand-anchored mdpocket site vanished with
    a fit-failed message. Atom counts are what cryptic_analysis uses for the
    same decision, and they are what this uses.

    Returns (kept chains, per-chain atom counts).
    """
    lig: list[tuple[float, float, float]] = []
    for chain in st[0]:
        for res in chain:
            if res.het_flag != "H" or res.name != comp_id:
                continue
            if copy_key and f"{chain.name}/{res.seqid.num}" != copy_key:
                continue
            lig.extend((a.pos.x, a.pos.y, a.pos.z) for a in res)
    if not lig:
        return [], {}
    c2 = cutoff * cutoff
    counts: dict[str, int] = {}
    for chain in st[0]:
        n = 0
        for res in chain:
            if res.het_flag != "A":
                continue
            for a in res:
                px, py, pz = a.pos.x, a.pos.y, a.pos.z
                if any((px - lx) ** 2 + (py - ly) ** 2 + (pz - lz) ** 2 <= c2
                       for lx, ly, lz in lig):
                    n += 1
        if n:
            counts[chain.name] = n
    if not counts:
        return [], {}
    top = max(counts.values())
    return sorted(c for c, n in counts.items() if n >= min_fraction * top), counts


def _chain_sequences(st, chains: list[str] | None = None) -> dict[str, dict[int, str]]:
    """Per-chain polymer sequence as {seqid -> residue name}.

    Keyed on author residue number rather than on position, because that is the
    numbering everything else in this module compares on, and because the
    protomers of a biological assembly of a homo-oligomer carry identical
    numbering by construction.
    """
    seqs: dict[str, dict[int, str]] = {}
    for chain in st[0]:
        if chains and chain.name not in chains:
            continue
        seq = {r.seqid.num: r.name for r in chain if r.het_flag == "A" and len(r)}
        if len(seq) >= 20:
            seqs[chain.name] = seq
    return seqs


def _homo_oligomer(
    st, chains: list[str] | None = None, min_identity: float = 0.9
) -> dict:
    """Detect several chains with identical or near-identical sequence.

    This is not a curiosity, it invalidates a specific measurement. The site
    signature used to track the SAME site across an ensemble is a set of residue
    NUMBERS with chain identity discarded. On a homotrimer the three protomers
    triplicate every number, so a 19-residue reference site collapses to 11
    distinct numbers and a C3-symmetric site cannot be resolved even in
    principle — any of the three symmetry copies, and pockets that touch none of
    them, score identically. Measured on apo TNF-alpha: 4 of 5 structures matched
    a pocket 7.7 A off-site sharing only residues 61 and 119, and the fifth
    matched a pocket 12.2 A away from those four. All five were reported as "the
    same site".

    So when this returns True the caller must NOT present
    `site_signature_overlap` as a same-site basis.
    """
    seqs = _chain_sequences(st, chains)
    names = sorted(seqs)
    groups: list[list[str]] = []
    for name in names:
        cur = seqs[name]
        for g in groups:
            ref = seqs[g[0]]
            shared = ref.keys() & cur.keys()
            if not shared or len(shared) < 0.5 * min(len(ref), len(cur)):
                continue
            if sum(ref[k] == cur[k] for k in shared) / len(shared) >= min_identity:
                g.append(name)
                break
        else:
            groups.append([name])
    biggest = max(groups, key=len) if groups else []
    return {
        "is_homo_oligomer": len(biggest) > 1,
        "n_identical_chains": len(biggest),
        "identical_chains": sorted(biggest) if len(biggest) > 1 else [],
        "n_polymer_chains": len(names),
        "sequence_identity_threshold": min_identity,
    }


def _target_polymer_chains(
    st, cif: Path | None, renamed: dict, accession: str | None
) -> set[str] | None:
    """Chain names of the polymer entities that ARE the target.

    `_chain_accessions` already maps chains to accessions from `_struct_ref`;
    this is the same read, expressed as a set and matched through UniProt's own
    merge history so a deposition naming a since-merged accession still counts
    (TL1A's Q8NFE9 IS O95150 — see `_accession_matches`).

    RETURNS None WHEN THE ACCESSION IS UNKNOWN, and None must be handled by the
    caller as "cannot filter", NEVER as "filter to nothing". Filtering to
    nothing would empty every site signature on an entry whose header would not
    parse, which converts an unreadable mapping into a confident wrong answer —
    the same fail-open/fail-closed distinction `_target_chains` makes.
    """
    if not accession or cif is None:
        return None
    names = [c.name for c in st[0]]
    chain_acc, status = _chain_accessions(cif, renamed, names)
    if status != "ok" or not chain_acc:
        return None
    out: set[str] = set()
    for ch, accs in chain_acc.items():
        for a in accs:
            if _accession_matches(a, accession)[0]:
                out.add(ch)
                break
    return out or None


def _polymer_ligand_control(
    st, work: Path, pid: str, want: list[str] | None,
    holo_call: dict, context, tgt_chains: Sequence[str], verified: bool,
    caller_chains: Sequence[str] | None, renamed: dict,
    anchor_comp: str | None, d: float = 1.6,
) -> dict | None:
    """THE PAIRED MEASUREMENT. Score the site with and without its polymer ligand.

    Returns None when there is nothing to decide — every polymer chain is the
    target or was named by the caller, which is the TNF-alpha and KRAS case.

    DO NOT REPLACE ONE NUMBER WITH THE OTHER. `volume_a3_stripped` is what gets
    reported, but a site that COLLAPSES when its polymer ligand is removed is an
    INDUCED-FIT / occluded site, not an absent one, and this project already has
    the calibration for exactly that: KRAS switch-II is druggability 0.708 on
    holo 6OIM and 0.000 on apo 4OBE AT THE SAME SITE, and reading the apo number
    as a verdict is the thirty-year KRAS error. On 8QFZ the pair is 283.6 A^3
    with the peptide and NO SITE AT ALL without it, and the pair is far more
    informative than either number. Neither number alone is interpretable; the
    difference between them is the finding.

    `induced_fit_signal` true must force `cryptic_pocket_risk: high` and a
    `tractability.caveat` saying the geometric number cannot be read in either
    direction without a holo SMALL-MOLECULE structure.
    """
    conj = _polymer_conjugate_host_chains(holo_call, context, renamed)
    prelim = _classify_polymer_chains(
        st, tgt_chains, verified, caller_chains, conj, None
    )
    undecided = [
        ch for ch, v in prelim.items()
        if v["class"] in ("polymer_ligand", "partner")
    ]
    if not undecided:
        return None

    anchor_res, _copy = (
        _ligand_site(st, anchor_comp, None) if anchor_comp else ([], None)
    )
    try:
        prep_wp, _kept_wp, _drop_wp = _prep(st, work, f"{pid}_withpolymer", want)
    except Exception as exc:  # noqa: BLE001
        return {"error": f"{type(exc).__name__}: {exc}", "stage": "prep_withpolymer"}
    pockets_wp = _fpocket_once(prep_wp, work, f"{pid}_withpolymer", d)
    anchored = _anchored_pocket(pockets_wp, anchor_res)

    share: dict[str, float] = {}
    if anchored and anchored.get("residues"):
        n = len(anchored["residues"])
        for ch, k in (anchored.get("lining_by_chain") or {}).items():
            share[ch] = round(k / n, 3)

    final = _classify_polymer_chains(
        st, tgt_chains, verified, caller_chains, conj, share
    )
    lig_chains = sorted(
        ch for ch, v in final.items() if v["class"] == "polymer_ligand"
    )

    if lig_chains:
        try:
            prep_s, _kept_s, _drop_s = _prep(
                st, work, f"{pid}_stripcheck", want, drop_chains=lig_chains
            )
        except Exception as exc:  # noqa: BLE001
            return {"error": f"{type(exc).__name__}: {exc}", "stage": "prep_stripped"}
        pockets_s = _fpocket_once(prep_s, work, f"{pid}_stripcheck", d)
    else:
        # Nothing to strip: the two inputs are the same file, so say so rather
        # than paying for an identical second fpocket run.
        pockets_s = pockets_wp
    anchored_s = _anchored_pocket(pockets_s, anchor_res)

    v_wp = anchored.get("volume") if anchored else None
    v_s = anchored_s.get("volume") if anchored_s else None
    nearest_d = nearest_v = None
    if anchored and not anchored_s and pockets_s:
        # SAME COORDINATE FRAME — one structure object, two atom selections — so
        # this distance is a real distance, unlike a cross-entry centroid
        # difference. It is how "no site" is distinguished from "the site moved".
        cands = [
            (_distance(anchored.get("centroid"), p.get("centroid")), p)
            for p in pockets_s if p.get("centroid")
        ]
        cands = [c for c in cands if c[0] is not None]
        if cands:
            nearest_d, near = min(cands, key=lambda c: c[0])
            nearest_v = near.get("volume")
    frac = (
        round(sum(share.get(c, 0.0) for c in lig_chains), 3)
        if (anchored and lig_chains) else None
    )
    return {
        "clustering_d": d,
        "anchor_comp_id": anchor_comp,
        "n_anchor_site_residues": len(anchor_res),
        "polymer_ligand_chains": lig_chains,
        "classification": final,
        # THE REPORTED NUMBER.
        "volume_a3_stripped": v_s,
        "volume_a3_with_polymer_ligand": v_wp,
        "site_present_when_stripped": v_s is not None,
        "n_pockets_with_polymer_ligand": len(pockets_wp),
        "n_pockets_stripped": len(pockets_s),
        "anchored_pocket_lining_by_chain": (
            anchored.get("lining_by_chain") if anchored else None
        ),
        "lining_fraction_from_polymer_ligand": frac,
        "nearest_pocket_when_site_absent_a": nearest_d,
        "nearest_pocket_when_site_absent_volume_a3": nearest_v,
        "induced_fit_signal": bool(v_wp is not None and v_s is None),
        "forces_cryptic_pocket_risk_high": bool(v_wp is not None and v_s is None),
        "_why": (
            "A site that exists only while its polymer ligand is present is an "
            "INDUCED-FIT / occluded site, not an absent one. This project has "
            "the calibration for that already: KRAS switch-II is 0.708 on holo "
            "6OIM and 0.000 on apo 4OBE at the same site, and reading the apo "
            "number as a verdict is the thirty-year KRAS error. Report the "
            "pair. Neither number alone is interpretable and the difference "
            "between them is the finding. volume_a3_stripped is what goes in "
            "tractability.pocket_volume_a3.primary_d1_6_a3; "
            "induced_fit_signal true forces cryptic_pocket_risk high and a "
            "caveat that the geometric number cannot be read in either "
            "direction without a holo SMALL-MOLECULE structure."
        ),
    }


def _centroid(coords: list[tuple[float, float, float]]) -> list[float] | None:
    if not coords:
        return None
    n = len(coords)
    return [round(sum(c[i] for c in coords) / n, 2) for i in range(3)]


def _pairs(items: list):
    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            yield items[i], items[j]


def _distance(a: list[float] | None, b: list[float] | None) -> float | None:
    if not a or not b:
        return None
    return round(sum((x - y) ** 2 for x, y in zip(a, b, strict=False)) ** 0.5, 2)


def _pdb_coords(path: Path, records: tuple[str, ...] = ("ATOM",)) -> list[tuple]:
    """Coordinates out of a PDB file, by fixed column, no parser dependency."""
    coords: list[tuple[float, float, float]] = []
    if not path.exists():
        return coords
    for line in path.read_text().splitlines():
        if not line.startswith(records) or len(line) < 54:
            continue
        try:
            coords.append(
                (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            )
        except ValueError:
            continue
    return coords


def _prank_rescore(
    out_dir: Path, work: Path, protein: Path
) -> tuple[dict[int, int], dict]:
    """Re-rank fpocket's pockets with PRANK.

    Returns ({fpocket_rank: prank_rank}, status) where status carries
    `prank_status` in {ok, failed, not_run} plus the reason and captured stderr
    on failure.

    The status field exists because the failure was INVISIBLE. `{}` on every
    error meant the caller saw `prank_rank: null` for every pocket, which is
    also exactly what a successful run that simply did not rank this pocket
    looks like. Observed: two runs on the same 4 PDB IDs, one with all
    `prank_rank: null` and no mention of p2rank anywhere in stdout, the next with
    real ranks. An entire rescoring stage disappeared and nothing in the output
    said so. A silent optional stage is worse than a missing one.

    `protein` is the SAME prepared PDB that was handed to fpocket. It is not
    optional and it is not cosmetic: `prank rescore` reads the dataset as two
    whitespace-separated columns, and a bare list of paths is rejected with
    "Dataset must contain 'protein' and 'prediction' columns!" — non-fatally,
    which is how this silently returned {} for every structure. The format is
    taken from p2rank's own shipped test_data/fpocket3.ds:

        PARAM.PREDICTION_METHOD=fpocket
        HEADER: prediction protein
        <stem>_out/<stem>_out.pdb   <stem>.pdb

    Why this exists: fpocket's DETECTION geometry is sound but its own ranking
    is the known weak link, and the best-recall configuration in the LIGYSIS
    benchmark of 13 predictors is fpocket detection + PRANK rescoring (60% top-
    N+2 recall, ahead of DeepPocket 58% and P2Rank standalone 52%).

    Measured on our own structures:
        6OIM switch-II   fpocket rank 9  ->  PRANK rank 2
        2AZ5 SPD304      fpocket rank 2  ->  PRANK rank 1

    Two gotchas, both confirmed by direct test, both load-bearing:

      * In RESCORE mode the `probability` column is NOT calibrated — the true
        SPD304 site scored 0.011 and a large decoy scored 0.783. Only the
        RANKING is usable here. `predict` mode's probability is well-calibrated
        (0.735 on the same site), so cross-structure probability needs a
        separate `predict` run.
      * `rescore` emits no `_residues.csv`; P2Rank only lays SAS points over the
        surface in `predict` mode.

    The `<protein>_rescored.csv` header, read off a real run rather than assumed:

        name,score,rank,old_rank,change,

    so `old_rank` is fpocket's rank and `rank` is PRANK's. An earlier parser
    looked for `pocket` and `new_rank`, neither of which P2Rank has ever
    emitted; it would have returned {} even once the dataset was accepted.

    A failure here is non-fatal by design: fpocket's own ranking survives and
    the caller sees prank_rank missing rather than a dead run.
    """
    def fail(reason: str, stderr: bytes | str = b"") -> tuple[dict, dict]:
        if isinstance(stderr, bytes):
            stderr = stderr.decode(errors="replace")
        return {}, {
            "prank_status": "failed",
            "prank_reason": reason,
            "prank_stderr": stderr.strip()[-1500:] or None,
        }

    pdbs = list(out_dir.glob("*_out.pdb"))
    if not pdbs or not protein.exists():
        # Nothing was ever handed to PRANK — fpocket produced no output to
        # rescore. Distinct from PRANK itself dying, and reported as such.
        return {}, {
            "prank_status": "not_run",
            "prank_reason": (
                "no fpocket *_out.pdb to rescore"
                if not pdbs
                else f"prepared protein missing: {protein.name}"
            ),
            "prank_stderr": None,
        }
    ds = work / f"{out_dir.name}_rescore.ds"
    ds.write_text(
        "PARAM.PREDICTION_METHOD=fpocket\n"
        "HEADER: prediction protein\n"
        f"{pdbs[0]}  {protein}\n"
    )
    outdir = work / f"{out_dir.name}_prank"
    shutil.rmtree(outdir, ignore_errors=True)
    try:
        proc = subprocess.run(  # noqa: S603
            [
                f"{P2RANK_HOME}/prank", "rescore", str(ds),
                "-o", str(outdir), "-c", "rescore_2024",
            ],
            check=False, capture_output=True, timeout=600,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        # A missing prank binary, a missing JRE or a timeout must cost the
        # rescore and nothing else — fpocket's own ranking is the deliverable
        # and prank_rank is the optional extra. Still non-fatal, now audible.
        return fail(f"{type(exc).__name__}: {exc}")
    if proc.returncode != 0:
        return fail(f"prank exited {proc.returncode}", proc.stderr)
    csvs = list(outdir.rglob("*_rescored.csv"))
    if not csvs:
        return fail("prank wrote no *_rescored.csv", proc.stderr)
    mapping: dict[int, int] = {}
    lines = csvs[0].read_text().splitlines()
    if not lines:
        return fail(f"{csvs[0].name} is empty", proc.stderr)
    hdr = [h.strip().lower() for h in lines[0].split(",")]
    if "rank" not in hdr:
        return fail(f"no 'rank' column in {csvs[0].name}: {hdr}", proc.stderr)
    i_new = hdr.index("rank")
    # `old_rank` is the authoritative fpocket rank. `name` ("pocket.9") carries
    # the same number and is the fallback, so a column rename upstream costs the
    # mapping rather than silently mis-keying it.
    i_old = hdr.index("old_rank") if "old_rank" in hdr else hdr.index("name")
    for line in lines[1:]:
        parts = [c.strip() for c in line.split(",")]
        if len(parts) <= max(i_old, i_new):
            continue
        digits = "".join(ch for ch in parts[i_old] if ch.isdigit())
        if not digits:
            continue
        try:
            mapping[int(digits)] = int(float(parts[i_new]))
        except ValueError:
            continue
    return mapping, {
        "prank_status": "ok",
        "prank_reason": None,
        "prank_stderr": None,
        "prank_pockets_ranked": len(mapping),
    }


def _parse_pockets(out_dir: Path) -> list[dict]:
    # Verified against real output: fpocket writes <input_stem>_out/ and inside
    # it <input_stem>_info.txt. removesuffix, not replace — replace() would eat
    # an "_out" occurring anywhere in the stem.
    info = out_dir / f"{out_dir.name.removesuffix('_out')}_info.txt"
    if not info.exists():
        return []
    pockets, cur = [], None
    for line in info.read_text().splitlines():
        s = line.strip()
        if s.startswith("Pocket") and s.endswith(":"):
            if cur:
                pockets.append(cur)
            cur = {"rank": int(s.split()[1])}
        elif cur is not None and ":" in s:
            k, _, v = s.partition(":")
            k = k.strip().lower().replace(" ", "_")
            try:
                cur[k] = float(v.strip())
            except ValueError:
                pass
    if cur:
        pockets.append(cur)

    for p in pockets:
        # fpocket numbers the per-pocket files from 1, matching "Pocket N :" in
        # info.txt exactly — there is no pocket0_atm.pdb. Checked on every run
        # in the calibration set. The old rank-1 shifted every pocket's residues
        # onto its neighbour and silently gave rank 1 no residues at all.
        atm = out_dir / "pockets" / f"pocket{p['rank']}_atm.pdb"
        res = set()
        names: dict[str, str] = {}
        coords: list[tuple[float, float, float]] = []
        if not atm.exists():
            # Never expected. Say so rather than reporting an empty pocket.
            p["residues_unavailable"] = atm.name
        else:
            for line in atm.read_text().splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    tag = f"{line[21]}/{line[22:26].strip()}"
                    res.add(tag)
                    names.setdefault(tag, line[17:20].strip())
            coords = _pdb_coords(atm, ("ATOM", "HETATM"))
        p["residues"] = sorted(
            res, key=lambda s: (s.split("/")[0], int(s.split("/")[1]))
        )
        # WHAT THE POCKET IS LINED WITH, not only where it is. fpocket's
        # druggability regression rewards a hydrophobic, sealed shape, and the
        # hydrophobic CORE of a folded domain is exactly that shape — measured
        # on IRAK4's death domain, a buried core scored 0.890 at rank 1 of 134
        # with nine Leu/Ile/Val/Phe, one Arg and one Tyr lining it. Composition
        # is half of what tells a core apart from a site; enclosure is the
        # other half (see `_buried_core_flag`).
        p["lining_residue_names"] = [names.get(t) for t in p["residues"]]
        # WHICH CHAINS LINE IT, ON EVERY POCKET, ALWAYS. The cheapest field in
        # this file and the one that would have caught the entire wrong-protein
        # failure class on the face of the output, with no accession lookup, no
        # classification and no judgement: MYC's headline pocket is 100% MAX,
        # IL-11's is 100% Q14626 (the receptor), RORgt's 6C1P pocket is 100% an
        # ion channel, IL-13's is inside an antibody Fab, and 8QFZ's is 6 of 10
        # residues on a 12-mer bicyclic peptide. Every one of those is a
        # `lining_by_chain` a reader can see. It is deliberately NOT conditional
        # on a chain resolver having succeeded, because the resolver failing is
        # one of the ways this goes wrong (`chain_accessions` was `{}` on every
        # entry of the retracted calibration set while the adjacent `_why`
        # asserted it resolved).
        p["lining_by_chain"] = {
            ch: sum(1 for r in p["residues"] if r.split("/")[0] == ch)
            for ch in sorted({r.split("/")[0] for r in p["residues"]})
        }
        p["n_apolar_lining_residues"] = sum(
            1 for n in names.values() if n in APOLAR_RESIDUES
        )
        p["apolar_lining_fraction"] = (
            round(p["n_apolar_lining_residues"] / len(names), 3) if names else None
        )
        # WHERE the pocket is, not just which residue numbers line it. An
        # overlap fraction cannot distinguish "the same site" from "a pocket
        # 12 A away that happens to share residue numbers", and on a
        # homo-oligomer sharing numbers is close to guaranteed. The centroid is
        # what makes that check possible downstream.
        p["centroid"] = _centroid(coords)
        p["n_lining_atoms"] = len(coords)
    return pockets


def _jaccard(a: list[str], b: list[str]) -> float | None:
    sa, sb = set(a), set(b)
    if not sa or not sb:
        return None
    return round(len(sa & sb) / len(sa | sb), 3)


# Conventional apolar side chains. GLY and PRO are deliberately excluded: they
# are conformational rather than hydrophobic, and including them would let a
# flexible loop look like a hydrophobic core.
APOLAR_RESIDUES = frozenset(
    "ALA VAL LEU ILE PHE MET TRP CYS".split()
)


def _unit_interval(value, field: str, context: str = ""):
    """Refuse to emit a [0,1] score that is not in [0,1].

    THIS EXISTS BECAUSE A REAL BUG SHIPPED THROUGH THIS GAP.
    `mdpocket.sites.*.druggability_by_structure` reported fpocket's
    `volume_score` descriptor — observed 3.35 to 4.00 — under a field name that
    invites a reader to quote it as a druggability probability. Nothing checked
    the range, so a 4.00 "druggability" left the function and was written into a
    dossier. A number carrying a [0,1] name and a value of 4.00 is not a noisy
    measurement; it is a different quantity wearing the wrong label, and the
    cheapest place to catch that is at the point of emission.

    `None` passes: a missing measurement is not an out-of-range one.
    """
    if value is None:
        return None
    v = float(value)
    if not (0.0 <= v <= 1.0):
        raise RuntimeError(
            f"{field}={v!r} is outside [0,1]"
            + (f" ({context})" if context else "")
            + ". A field named as a [0,1] score must never carry a value "
            "outside it: that is the signature of a wrong descriptor column, "
            "not of a noisy estimate. Refusing to emit it."
        )
    return v


# ===========================================================================
# disorder — metapredict, with the cardinal rule preserved
# ===========================================================================


# {accession -> (aliases, resolved_ok)} for the life of one container.
_ACC_ALIAS_CACHE: dict[str, dict] = {}


def _accession_aliases(acc: str) -> dict:
    """Every accession UniProt considers the same entry as `acc`.

    THIS IS WHAT MAKES FAILING CLOSED SAFE. A PDB entry declares the accession
    that was current when it was deposited, and UniProt merges accessions
    afterwards — so an older entry for the same protein legitimately names a
    different string. Caught before shipping, on the TL1A ensemble: 2O0O, 2QE3
    and 2RE9 declare **Q8NFE9**, whose UniProt record is
    `inactiveReason: {inactiveReasonType: "MERGED", mergeDemergeTo: ["O95150"]}`
    — it IS O95150 — while 3K51 and the newer entries declare O95150 directly.
    Refusing on a literal string comparison would have thrown away three of six
    entries of the target's own ensemble and called them "not this protein".

    Returns `{aliases, genes, taxid, ok}`. `ok` is False when the lookup could
    not be made, and a caller must NOT refuse an entry on an unresolved record —
    an unanswered question is not a negative answer.
    """
    acc = (acc or "").strip()
    if not acc:
        return {"aliases": set(), "genes": set(), "taxid": None, "ok": False}
    if acc in _ACC_ALIAS_CACHE:
        return _ACC_ALIAS_CACHE[acc]
    import json as _json
    import urllib.error
    import urllib.request

    rec = {
        "aliases": {acc, acc.split("-")[0]},
        "genes": set(),
        "taxid": None,
        "ok": False,
    }
    try:
        req = urllib.request.Request(  # noqa: S310
            f"https://rest.uniprot.org/uniprotkb/{acc}.json"
            "?fields=accession,gene_names,organism_id",
            headers={"User-Agent": "pocket-scan"},
        )
        with urllib.request.urlopen(req, timeout=20) as fh:  # noqa: S310
            d = _json.loads(fh.read().decode())
        rec["ok"] = True
        if d.get("primaryAccession"):
            rec["aliases"].add(d["primaryAccession"])
        for s in d.get("secondaryAccessions") or ():
            rec["aliases"].add(s)
        # An accession that has been merged away answers with entryType
        # "Inactive" and names its successor here.
        for s in (d.get("inactiveReason") or {}).get("mergeDemergeTo") or ():
            rec["aliases"].add(s)
        for g in d.get("genes") or ():
            name = (g.get("geneName") or {}).get("value")
            if name:
                rec["genes"].add(name.strip().upper())
        rec["taxid"] = (d.get("organism") or {}).get("taxonId")
    except (urllib.error.URLError, OSError, ValueError, KeyError):
        rec["ok"] = False
    _ACC_ALIAS_CACHE[acc] = rec
    return rec


def _accession_matches(declared: str, target: str) -> tuple[bool, bool]:
    """Is `declared` the same protein as `target`? Returns (matches, checkable).

    THREE TESTS, AND THE THIRD IS THE ONE THAT KEEPS FAILING CLOSED HONEST.
    String equality alone is far too brittle to refuse an entry on:

      1. the accession, or its isoform base, is literally the same;
      2. UniProt has MERGED one into the other — TL1A's 2O0O/2QE3/2RE9 declare
         Q8NFE9, whose record says `mergeDemergeTo: ["O95150"]`;
      3. they are the SAME GENE in the SAME ORGANISM under two different
         UniProt entries. IL-13's 3BPO declares **Q4VB50**, which is not merged
         into P35225 and never will be — it is an unreviewed TrEMBL entry whose
         recommended name is "Interleukin-13" and whose gene is IL13, human.
         Refusing 3BPO on the string would have been wrong in the damaging
         direction: the entry DOES contain IL-13, and what it actually needs is
         for chain A to be recognised as the target so that chains B and C —
         IL-4R-alpha and IL-13R-alpha-1 — are recognised as partners.

    `checkable` is False when neither side could be resolved against UniProt.
    An unanswered question is not a negative answer and must not refuse.
    """
    if declared == target or declared.split("-")[0] == target.split("-")[0]:
        return True, True
    d, t = _accession_aliases(declared), _accession_aliases(target)
    if d["aliases"] & t["aliases"]:
        return True, True
    if (
        d["genes"] and t["genes"] and d["genes"] & t["genes"]
        and d["taxid"] is not None and d["taxid"] == t["taxid"]
    ):
        return True, True
    return False, (d["ok"] and t["ok"])


def _assembly_base_chain(name: str) -> str:
    """`"A-3"` -> `"A"`. The source chain an assembly-expansion copy came from.

    gemmi appends `-<n>` when it expands a biological assembly, and
    `_struct_ref_seq` only ever names the un-suffixed strand. Everything that
    joins a chain to an accession has to strip this or the copies come back
    unmapped — see `_chain_accessions`.
    """
    import re

    m = re.match(r"^(.+)-\d+$", name)
    return m.group(1) if m else name


def _chain_accessions(
    cif: Path | None,
    renamed: dict[str, str],
    chain_names: list[str] | None = None,
) -> tuple[dict[str, list[str]], str]:
    """{chain id -> UniProt accessions}, from `_struct_ref` + `_struct_ref_seq`.

    WHICH CHAIN IS THE TARGET IS A LOOKUP, NOT A GUESS. Everything that used to
    pick the target by chain LENGTH is wrong the moment a partner is bigger than
    the target, and on a GPCR-G-protein complex the partner always is. Measured
    on S1PR1: G-beta-1 is 331-338 residues against the receptor's 278-290 in all
    four entries, so `_one_letter`'s "longest chain" picked G-beta every time.
    The interface stage then split 7TD4 into target ["B"] and partner
    ["A","G","R"] — chain B is G-beta-1, chain R is S1PR1 — computed the
    G-beta/G-alpha-G-gamma interface, and reported it as the target's epitope
    with `interface_status: ok`, 93 interface residues and no warning at all.

    Chain IDs are returned in THIS module's namespace, i.e. after `_load`'s
    single-character renaming, so they can be compared against `used_chains`.

    Returns `(mapping, status)`. THE STATUS IS NOT DECORATION — it is the
    difference between "this entry says it does not contain your protein" and
    "we could not read what this entry says", and those must not be collapsed.
    Collapsing them is what made the resolver fail OPEN: an empty mapping was
    treated as "no mapping declared", the caller fell back to every chain, and
    the pocket selector was then free to pick a cavity inside an antibody.

        ok            the entry declares UniProt refs and they were read
        no_header     the header file could not be fetched
        unparsable    it was fetched and gemmi would not read it
        no_unp_refs   it parsed and declares no UNP reference at all

    ASSEMBLY-EXPANSION COPIES INHERIT THEIR SOURCE CHAIN'S ACCESSION. gemmi
    names the copies it creates when expanding a biological assembly
    `<orig>-<n>`, and `_struct_ref_seq` only ever names `<orig>`. Without this
    every copy came back unmapped: 1JH5 is a 60-mer of ONE protein in which ten
    chains resolved to Q9Y275 and the other FIFTY — the same protein, renamed by
    expansion — carried `chain_accessions: null` and were reported as
    `non_target_chains_scored`. A pocket lined by any of those fifty is on the
    target and was being described as if it were not.
    """
    import gemmi

    if cif is None:
        return {}, "no_header"
    try:
        block = gemmi.cif.read(str(cif)).sole_block()
    except Exception:  # noqa: BLE001
        return {}, "unparsable"
    ref_acc: dict[str, str] = {}
    for row in block.find("_struct_ref.", ["id", "db_name", "pdbx_db_accession"]):
        if row.str(1).strip().upper() != "UNP":
            continue
        acc = row.str(2).strip()
        if acc and acc not in ("?", "."):
            ref_acc[row.str(0).strip()] = acc
    by_strand: dict[str, list[str]] = {}
    for row in block.find("_struct_ref_seq.", ["ref_id", "pdbx_strand_id"]):
        acc = ref_acc.get(row.str(0).strip())
        if not acc:
            continue
        for raw in row.str(1).replace(",", " ").split():
            if acc not in by_strand.setdefault(raw, []):
                by_strand[raw].append(acc)
    if not by_strand:
        return {}, "no_unp_refs"
    inv = {new: old for old, new in renamed.items()}
    names = (
        list(chain_names) if chain_names is not None
        else [renamed.get(r, r) for r in by_strand]
    )
    out: dict[str, list[str]] = {}
    for ch in names:
        orig = inv.get(ch, ch)
        for key in (orig, _assembly_base_chain(orig)):
            if key in by_strand:
                out[ch] = list(by_strand[key])
                break
    return out, "ok"


def _target_chains(
    chain_acc: dict[str, list[str]],
    accession: str | None,
    fallback: list[str],
    acc_status: str = "ok",
) -> dict:
    """The chains that ARE the target, by accession. FAILS CLOSED.

    Returns {chains, basis, verified, refuse, reason, declared_accessions}.

    THIS USED TO FAIL OPEN AND THAT IS WHAT LET THE POCKET SELECTOR LOOSE.
    When the requested accession was absent from the entry it returned EVERY
    chain with the note "using every chain scored, which may include partners" —
    a warning in a string, downstream of nothing. Measured on IL-13 3BPO, whose
    `_struct_ref` declares Q4VB50, P24394 and P78552 and does NOT declare
    P35225: all three chains became "target", the longest of them (IL-13R-alpha-1
    at 314 aa) became the target SEQUENCE, and the interface stage then put the
    receptor on the target's side of its own interface and returned
    `interface_status: ok`.

    THE THREE OUTCOMES ARE DIFFERENT AND ARE NOW DISTINGUISHED:

      * the entry declares UniProt refs and one matches   -> verified, use them
      * the entry declares UniProt refs and NONE match    -> REFUSE the entry.
        It says what proteins it contains and yours is not among them. Scoring
        it anyway measures a different molecule, which is precisely how a pocket
        inside tralokinumab's Fab became IL-13's headline volume.
      * the mapping could not be READ (no header, unparsable, no UNP refs)
        -> fall back to every chain, but `verified` is False, and everything
        downstream that would otherwise trust the chain set must degrade rather
        than assume. Unreadable is not the same claim as absent and must not
        refuse an entry on our own parser's behalf.

    With no accession supplied there is nothing to verify against; that is the
    caller's choice, not a failure, and it is reported as unverified rather than
    refused.
    """
    declared = sorted({a for accs in chain_acc.values() for a in accs})
    if not accession:
        return {
            "chains": list(fallback),
            "basis": (
                "no uniprot_accession supplied; every scored chain is treated "
                "as target and NOTHING below is verified against an accession"
            ),
            "verified": False,
            "refuse": False,
            "reason": None,
            "declared_accessions": declared,
        }
    if acc_status != "ok":
        return {
            "chains": list(fallback),
            "basis": (
                f"the entry's UniProt mapping could not be read ({acc_status}); "
                "every scored chain is treated as target and the chain set is "
                "UNVERIFIED"
            ),
            "verified": False,
            "refuse": False,
            "reason": (
                "unreadable is not the same finding as absent, so this entry is "
                "not refused — but no pocket from it has been checked against "
                "the target's chains, and a pocket lined by a partner cannot be "
                "distinguished from one lined by the target here."
            ),
            "declared_accessions": declared,
        }
    # Matched through UniProt's own merge history, not by string equality — see
    # `_accession_aliases`. `checkable` False anywhere means we could not
    # establish that a declared accession is NOT the target, and an unanswered
    # question must not refuse an entry.
    hits, checkable = [], True
    for c in fallback:
        for a in chain_acc.get(c, ()):
            m, ck = _accession_matches(a, accession)
            checkable = checkable and ck
            if m:
                hits.append(c)
                break
    if not hits and not checkable:
        return {
            "chains": list(fallback),
            "basis": (
                f"no chain of this entry maps to {accession} by string, and "
                "UniProt could not be reached to check whether any of "
                f"{declared} is a merged form of it; every scored chain is "
                "treated as target and the chain set is UNVERIFIED"
            ),
            "verified": False,
            "refuse": False,
            "reason": (
                "the entry declares accessions this run could not resolve. "
                "Older depositions legitimately name accessions UniProt has "
                "since merged — TL1A's 2O0O, 2QE3 and 2RE9 declare Q8NFE9, "
                "which IS O95150 — so refusing on the string alone would "
                "discard half an ensemble. Not refused, not verified."
            ),
            "declared_accessions": declared,
        }
    if not hits:
        return {
            "chains": [],
            "basis": f"no chain of this entry maps to {accession}",
            "verified": True,
            "refuse": True,
            "reason": (
                f"this entry declares UniProt accessions {declared or '[]'} in "
                f"_struct_ref, and none of them is {accession} or a UniProt "
                "merge of it, so it does not "
                "contain the target. REFUSED rather than scored: the previous "
                "behaviour was to use every chain 'which may include partners', "
                "and what that produced was pocket volumes measured inside "
                "antibody Fabs and on receptor chains, reported as the target's."
            ),
            "declared_accessions": declared,
        }
    return {
        "chains": hits,
        "basis": f"chains mapping to {accession} in _struct_ref_seq",
        "verified": True,
        "refuse": False,
        "reason": None,
        "declared_accessions": declared,
    }


# Fraction of a pocket's lining residues that must sit on the target's own
# chains before the pocket may be SELECTED as that target's site.
#
# PROPOSED, NOT CALIBRATED. A majority, not unanimity, deliberately: a genuine
# orthosteric pocket at a target/partner interface is legitimately lined by
# both, and requiring 1.0 would refuse exactly the pockets rule 2b exists to
# find. The failures this catches are not marginal — they are pockets lined
# ENTIRELY by a partner, at on-target fraction 0.00: cavities inside the Fabs of
# tralokinumab and lebrikizumab reported as IL-13's site, a cavity inside
# belimumab reported as BAFF's, cavities inside rituximab's Fab reported as
# CD20's. Anything near this boundary should be read off
# `on_target_residue_fraction` directly rather than off the flag.
POCKET_MIN_ON_TARGET_FRACTION = 0.5


def _annotate_on_target(
    pockets: list[dict], target_chains: list[str] | None, verified: bool
) -> None:
    """Mark each pocket with how much of it is actually on the target.

    THE RESOLVER WORKED AND THE SELECTOR IGNORED IT. `target_chains` was
    resolved by UniProt accession, announced in `target_chains_basis` with a
    `_why` naming the case it was built for — and then the site pocket was
    chosen as the most druggable pocket ANYWHERE in the file, with no check that
    a single lining residue was on those chains. Measured, selected pockets that
    were actually on the target:

        IL-13   1 of 9    the rest inside the Fabs of tralokinumab and
                          lebrikizumab (3L5X, 5L6Y, 3L5W, 4PS4, 3G6D) and on the
                          receptor chain (3LB6)
        BAFF    2 of 5    5Y9J lined by belimumab; no fully on-target pocket
                          exists among its 22
        CD20    4 of 7    6Y90 and 6Y97 lined by rituximab's Fab

    Filtering on this inverts a verdict-relevant number: IL-13's median volume
    moves 312.3 -> 106.8 A^3, from above the (now suspended) druggable bound to
    below the hard bound. BAFF 258.3 -> 177.4. CD20 281.0 -> 242.3. Nothing in
    the payload flagged any of it — the figure was uniformly PRESENT and quietly
    measuring a different molecule, which is the twin of a field that is
    uniformly null and reads as "not measured".

    `on_target` is None when the chain set is unverified, and None never
    excludes: an unreadable accession mapping must not silently drop pockets. It
    is a False that excludes, and a False requires a verified chain set.
    """
    tgt = set(target_chains or ())
    for p in pockets:
        res = p.get("residues") or []
        chains = [r.split("/")[0] for r in res]
        on = [c for c in chains if c in tgt]
        p["n_on_target_lining_residues"] = len(on)
        p["on_target_residue_fraction"] = (
            round(len(on) / len(res), 3) if res else None
        )
        p["off_target_lining_chains"] = sorted({c for c in chains if c not in tgt})
        # The same number as `on_target_residue_fraction` read from the other
        # end, and it is the one a reader scans for. Reported whether or not the
        # chain set is verified: an unverified `1.0` here says "every lining
        # residue is on a chain this run could not attribute", which is a
        # finding, not a blank.
        p["lining_fraction_non_target"] = (
            round(1.0 - (len(on) / len(res)), 3) if res else None
        )
        p["on_target"] = (
            None if (not verified or not res)
            else (len(on) / len(res)) >= POCKET_MIN_ON_TARGET_FRACTION
        )


# A pocket is called symmetry-axis-anchored when at least this many residue
# NUMBERS are contributed to it by two or more sequence-identical chains.
# PROPOSED, NOT CALIBRATED. BAFF's axial site carries four (Gln144, Phe194,
# Leu282, Leu284, each from three protomers) with no ligand anywhere; one shared
# number is a chain contact, not an axis.
SYMMETRY_AXIS_MIN_SHARED_RESIDUES = 2

_UNP_SITE_CACHE: dict[str, tuple[dict[int, list[str]], bool]] = {}


def _uniprot_functional_sites(acc: str) -> tuple[dict[int, list[str]], bool]:
    """{UniProt sequence position -> feature names} for binding/active/site.

    THE ANNOTATION THAT DOES NOT NEED CHEMISTRY. Every other external label a
    pocket can carry needs something bound, or a partner, or a homolog with
    something bound. This one is a curated statement about the protein itself,
    which is exactly what a target with no ligand and no complex still has.
    """
    acc = (acc or "").strip()
    if not acc:
        return {}, False
    if acc in _UNP_SITE_CACHE:
        return _UNP_SITE_CACHE[acc]
    import json as _json
    import urllib.error
    import urllib.request

    out: dict[int, list[str]] = {}
    ok = False
    try:
        req = urllib.request.Request(  # noqa: S310
            f"https://rest.uniprot.org/uniprotkb/{acc}.json"
            "?fields=ft_binding,ft_act_site,ft_site",
            headers={"User-Agent": "pocket-scan"},
        )
        with urllib.request.urlopen(req, timeout=20) as fh:  # noqa: S310
            d = _json.loads(fh.read().decode())
        ok = True
        for ft in d.get("features") or ():
            kind = ft.get("type")
            if kind not in ("Binding site", "Active site", "Site"):
                continue
            loc = ft.get("location") or {}
            beg = (loc.get("start") or {}).get("value")
            end = (loc.get("end") or {}).get("value")
            if beg is None:
                continue
            label = ft.get("description") or kind
            for pos in range(int(beg), int(end or beg) + 1):
                if label not in out.setdefault(pos, []):
                    out[pos].append(f"{kind}: {label}"[:60])
    except (urllib.error.URLError, OSError, ValueError, KeyError, TypeError):
        ok = False
    _UNP_SITE_CACHE[acc] = (out, ok)
    return out, ok


def _chain_unp_offsets(
    cif: Path | None, renamed: dict[str, str], chain_names: list[str],
) -> dict[str, int]:
    """{chain -> offset} such that `auth_residue_number = unp_position + offset`.

    Read straight out of `_struct_ref_seq`'s own alignment columns, so a
    construct numbered from its own start and a construct numbered on the
    precursor both land on the UniProt sequence correctly. Without this every
    UniProt feature position would be compared against a residue number that
    means something else — the same class of error as `match_by="seqid"` across
    two entries, and it has already bitten this file once on TL1A.
    """
    import gemmi

    if cif is None:
        return {}
    try:
        block = gemmi.cif.read(str(cif)).sole_block()
    except Exception:  # noqa: BLE001
        return {}
    inv = {new: old for old, new in renamed.items()}
    by_strand: dict[str, int] = {}
    for row in block.find(
        "_struct_ref_seq.",
        ["pdbx_strand_id", "db_align_beg", "pdbx_auth_seq_align_beg"],
    ):
        try:
            off = int(row.str(2)) - int(row.str(1))
        except (ValueError, TypeError):
            continue
        for raw in row.str(0).replace(",", " ").split():
            by_strand.setdefault(raw, off)
    out: dict[str, int] = {}
    for ch in chain_names:
        orig = inv.get(ch, ch)
        for key in (orig, _assembly_base_chain(orig)):
            if key in by_strand:
                out[ch] = by_strand[key]
                break
    return out


def _annotate_pocket_labels(
    pockets: list[dict],
    target_chains: list[str] | None,
    chain_acc: dict[str, list[str]],
    identical_chains: list[str] | None,
    unp_sites: dict[int, list[str]],
    unp_offsets: dict[str, int],
) -> None:
    """Tag every pocket with EVERY external label that applies to it.

    ANCHORING IS AN ANNOTATION, NOT AN ELECTION. The old design chose one pocket
    as "the" site and, when nothing external applied, fell back to whichever
    scored highest — and that fallback is where all four bad calibration anchors
    were born: MYC's pocket was on MAX, IL-11's on IL-11 receptor alpha, IL-13's
    inside tralokinumab, CD20's on a cholesterol-hemisuccinate site. A pocket
    that carries no external label is not "the site by default"; it is a pocket
    with no external label, and saying so is a true and useful statement about a
    protein.

    A pocket may carry several labels or none. The labels are:

        ligand_site               overlaps a drug-like co-crystallised ligand
        interface                 overlaps a partner epitope (added later, by
                                  the interface stage, which is where the
                                  epitope exists)
        symmetry_axis             built from equivalent residues of two or more
                                  identical chains, as BAFF's axial site is with
                                  no ligand anywhere
        annotated_functional_site overlaps a UniProt binding/active/site feature
        buried_core               the existing geometry flag (added later, with
                                  enclosure, by the interface stage)
        transferred_homolog_site  NOT AVAILABLE HERE — it needs Foldseek, which
                                  lives in `structure-select`/`neighbour_
                                  precedent`. Its absence is reported rather
                                  than left to look like a negative.
    """
    tgt = set(target_chains or ())
    ident = set(identical_chains or ())
    for p in pockets:
        res = p.get("residues") or []
        labels: list[str] = []
        detail: dict = {}

        jac = p.get("jaccard_vs_ligand_site")
        if jac:
            labels.append("ligand_site")
            detail["ligand_site_jaccard"] = jac

        # --- symmetry axis ---------------------------------------------
        if len(ident) > 1:
            by_num: dict[str, set[str]] = {}
            for r in res:
                ch, _, num = r.partition("/")
                if ch in ident:
                    by_num.setdefault(num, set()).add(ch)
            shared = {n: cs for n, cs in by_num.items() if len(cs) > 1}
            if len(shared) >= SYMMETRY_AXIS_MIN_SHARED_RESIDUES:
                labels.append("symmetry_axis")
                detail["symmetry_axis"] = {
                    "n_shared_residue_numbers": len(shared),
                    "n_chains": len(set().union(*shared.values())),
                    "residues": sorted(shared, key=lambda x: int(x))[:12],
                }

        # --- UniProt functional features -------------------------------
        if unp_sites:
            hits: list[str] = []
            for r in res:
                ch, _, num = r.partition("/")
                if ch not in tgt or ch not in unp_offsets:
                    continue
                try:
                    pos = int(num) - unp_offsets[ch]
                except ValueError:
                    continue
                for lab in unp_sites.get(pos, ()):
                    tag = f"{num}:{lab}"
                    if tag not in hits:
                        hits.append(tag)
            if hits:
                labels.append("annotated_functional_site")
                detail["annotated_functional_site"] = hits[:8]

        p["anchor_labels"] = labels
        p["anchor_detail"] = detail
        p["lining_chains"] = sorted({r.partition("/")[0] for r in res})
        p["lining_chain_accessions"] = {
            c: chain_acc.get(c) for c in p["lining_chains"]
        }


def _pocket_table(pockets: list[dict]) -> list[dict]:
    """The compact per-pocket record. THIS IS THE PRIMARY OUTPUT NOW.

    One row per returned pocket, small enough that thirty of them fit inside a
    payload cap: rank, size, score, where it is, what it is made of, and which
    external labels it carries. NO PROSE PER POCKET — the explanations live once
    per clustering value in `on_target_selection` and `anchor_summary`, because
    a payload that truncates deletes its own trailing explanation first.

    A distribution replaces a selection. Reporting one elected pocket is a
    maximum over N draws AND it is where every bad anchor came from; a table of
    thirty with nine of them on chain B cannot hide either.
    """
    return [
        {
            "rank": p.get("rank"),
            "prank_rank": p.get("prank_rank"),
            "volume_a3": p.get("volume"),
            "druggability": p.get("druggability_score"),
            "score": p.get("score"),
            "n_lining_residues": len(p.get("residues") or []),
            "chains": p.get("lining_chains"),
            "chain_accessions": p.get("lining_chain_accessions"),
            "on_target_fraction": p.get("on_target_residue_fraction"),
            "on_target": p.get("on_target"),
            "anchors": p.get("anchor_labels"),
            "anchor_detail": p.get("anchor_detail") or None,
            "centroid": p.get("centroid"),
        }
        for p in pockets
    ]


def _one_letter(st, chains: list[str] | None = None) -> tuple[str | None, str | None]:
    """Longest polymer chain as a one-letter sequence, and which chain it was.

    `chains` MUST be the target's chains, not the whole assembly — see
    `_chain_accessions`. Longest-wins is only correct inside one accession.
    """
    import gemmi

    best_seq, best_chain = "", None
    for chain in st[0]:
        if chains and chain.name not in chains:
            continue
        names = [r.name for r in chain if r.het_flag == "A" and len(r)]
        if len(names) < 20:
            continue
        seq = gemmi.one_letter_code(names).upper()
        seq = "".join(c for c in seq if c.isalpha())
        if len(seq) > len(best_seq):
            best_seq, best_chain = seq, chain.name
    return (best_seq or None), best_chain


def _fetch_header(pdb_id: str, dest: Path) -> Path | None:
    """The HEADER-ONLY mmCIF, which is the only copy carrying `_struct_ref`.

    RCSB STRIPS `_struct_ref` FROM ASSEMBLY FILES. Verified on 6OIM: the
    assembly1 CIF this module fetches for coordinates has no `_struct_ref`
    category at all, so the UniProt accession — the thing that decides which
    chain is the target and which sequence disorder is measured on — is simply
    not in the file the rest of this module reads. `files.rcsb.org/header/<ID>.cif`
    is ~100 kB, carries `_struct_ref` and `_struct_ref_seq` in full, and is
    fetched once per entry.

    Returns None rather than raising: an entry with no retrievable header costs
    the accession and nothing else, and the loss is reported by its callers.
    """
    hdr = dest / f"{pdb_id}_header.cif"
    if hdr.exists():
        return hdr
    try:
        with urllib.request.urlopen(  # noqa: S310
            f"https://files.rcsb.org/header/{pdb_id}.cif", timeout=60
        ) as r:
            text = r.read().decode("utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        return None
    # RCSB'S HEADER FILE IS NOT VALID mmCIF AS SERVED. It is the full entry with
    # the coordinate loops deleted, and the deletion leaves the `_atom_site`
    # loop's `loop_` keyword behind with no tags and no rows. gemmi rejects the
    # whole document on it — "parse error" at the last line — so every
    # accession lookup came back empty and every target chain fell back to
    # longest-wins, which is the bug this file exists to close. Trimming the
    # dangling keyword costs nothing and there is no other way to read the
    # category we need.
    # They occur mid-file as well as at the end — 4OBE has three consecutive
    # `loop_` / `#` pairs at line 1031 where the coordinate, anisotropic-B and
    # one other loop were removed — so this drops EVERY `loop_` not followed by
    # a tag, not just a trailing one.
    src_lines = text.splitlines()
    keep: list[str] = []
    for i, line in enumerate(src_lines):
        if line.strip() == "loop_":
            nxt = next(
                (s for s in (x.strip() for x in src_lines[i + 1:])
                 if s and s != "#"),
                "",
            )
            if not nxt.startswith("_"):
                continue  # a loop with no tags: not mmCIF, drop the keyword
        keep.append(line)
    hdr.write_text("\n".join(keep) + "\n")
    return hdr


def _structure_context(pdb_id: str, dest: Path, accession: str | None):
    """`ligand_filter.StructureContext` for one entry, or None.

    FROM THE HEADER, NEVER FROM THE ASSEMBLY. RCSB strips `_struct_conn` from
    assembly files exactly as it strips `_struct_ref` — verified on 8QFZ, whose
    `-assembly1.cif` carries 23 categories and neither of those two. A context
    built from the coordinate file this module already holds comes back with an
    empty link table that is INDISTINGUISHABLE FROM "nothing is covalently
    bonded", and `LFI` goes straight back to `druglike`. `ligand_filter` detects
    this itself (`StructureContext.has_struct_conn_category` is False, every
    verdict is flagged `struct_conn_absent_from_context` and confidence drops),
    but the right fix is to hand it the header.

    `_fetch_header` is already called once per entry for the accession, and it
    caches on disk, so this costs NO additional network call.
    """
    import ligand_filter as LF

    hdr = _fetch_header(pdb_id, dest)
    if hdr is None:
        return None
    try:
        return LF.StructureContext.from_mmcif_path(
            hdr, entry_id=pdb_id, target_accession=accession)
    except Exception:  # noqa: BLE001
        return None


def _uniprot_from_cif(cif: Path | None) -> list[str]:
    """UniProt accessions the entry itself declares, from `_struct_ref`.

    THE ACCESSION IS IN THE FILE. Not resolving it is what made disorder
    silently measure the wrong molecule: with no `uniprot_accession` argument
    the stage fell back to the crystallised construct, and IRAK4 came back
    0.0 over 284 residues — the kinase domain — where the full 460-residue
    protein is 0.1413 with a disordered region at 101-162. A deposited entry
    is the ORDERED part of a protein by selection, so that fallback does not
    understate the answer by a little, it answers a different question, and a
    0.0 reads as "no disorder" rather than "not measured".

    Returns [] when the entry declares none (a synthetic construct, a peptide),
    which is a real state and is reported as such rather than guessed at.
    """
    import gemmi

    if cif is None:
        return []
    try:
        block = gemmi.cif.read(str(cif)).sole_block()
    except Exception:  # noqa: BLE001
        return []
    found: list[str] = []
    for row in block.find("_struct_ref.", ["db_name", "pdbx_db_accession"]):
        if row.str(0).strip().upper() != "UNP":
            continue
        acc = row.str(1).strip()
        if acc and acc not in ("?", ".") and acc not in found:
            found.append(acc)
    return found


def _disorder_block(
    accession: str | None, sequence: str | None, sequence_source: str | None,
    accession_source: str | None = None,
) -> dict:
    """Disorder fraction for the target, never fatal, never silently zero.

    THE CARDINAL RULE, which is disorder.py's and is preserved here: 0.0 means
    genuinely folded. It must never mean "something went wrong". predict_disorder
    returns None on failure, so a failure becomes disorder_status=failed with
    disorder_fraction=None — never 0.0.

    Validation targets: MYC P01106 ~0.83 disordered, CDK2 P24941 ~0.00. A run
    that returns 0.0 for MYC has a broken metapredict, not a folded MYC.

    An accession is strongly preferred over a structure-derived sequence. A
    crystallised construct is, by selection, the ordered part of the protein —
    MYC's deposited fragments are short helices — so a sequence lifted from a
    PDB entry systematically UNDERSTATES disorder. When that path is used it is
    labelled, and the caveat travels with the number.
    """
    out: dict = {
        "disorder_status": "not_run",
        "disorder_reason": None,
        "disorder_fraction": None,
        "accession": accession,
        "accession_source": accession_source,
        "sequence_source": sequence_source,
        # WHAT WAS MEASURED, always, before any number is read. A fraction with
        # no scope is not interpretable: 0.0 over a 284-residue crystallised
        # kinase domain and 0.0 over a 460-residue protein are different claims
        # and only one of them is about the target.
        "scope": None,
        "is_full_length_sequence": None,
        "n_residues_measured": None,
        "construct_disorder_fraction": None,
    }
    if not accession and not sequence:
        out["disorder_reason"] = (
            "no uniprot_accession supplied and no polymer sequence could be "
            "read from any structure"
        )
        return out
    try:
        import disorder as _disorder

        res = _disorder.predict_disorder(
            sequence=sequence if not accession else None, accession=accession
        )
    except Exception as exc:  # noqa: BLE001
        out.update(disorder_status="failed",
                   disorder_reason=f"{type(exc).__name__}: {exc}")
        return out
    if res is None:
        out.update(
            disorder_status="failed",
            disorder_reason=(
                "predict_disorder returned None — every method failed. NOT 0.0: "
                "a failed prediction is not a folded protein."
            ),
        )
        return out
    fraction = res.get("disorder_fraction")
    full_length = bool(accession)
    out.update(
        disorder_status="ok",
        method=res.get("method"),
        confidence=res.get("confidence"),
        length=res.get("length"),
        n_residues_measured=res.get("length"),
        disordered_regions=res.get("disordered_regions"),
        n_disordered_regions=res.get("n_disordered_regions"),
        longest_disordered_region=res.get("longest_disordered_region"),
        fallback_from=res.get("fallback_from"),
        source=res.get("source"),
        is_full_length_sequence=full_length,
        scope="full_length_uniprot" if full_length else "crystallised_construct",
    )
    if full_length:
        out["disorder_fraction"] = _unit_interval(
            fraction, "disorder_fraction", f"accession {accession}"
        )
        return out

    # ---- construct-only path ---------------------------------------------
    # `disorder_fraction` STAYS NULL. This is not squeamishness: the field is
    # read straight into the dossier's `tractability.disorder_fraction`, which
    # is a statement about the protein, and a construct measurement placed there
    # is a statement about a different molecule. Measured on IRAK4: the
    # crystallised kinase domain gives 0.0 over 284 residues while the full
    # 460-residue protein is 0.1413 with a disordered region at 101-162. A bare
    # 0.0 in that slot reads as "no disorder", not as "not measured", and it is
    # the second reading that is true.
    out["construct_disorder_fraction"] = _unit_interval(
        fraction, "construct_disorder_fraction", str(sequence_source)
    )
    out["disorder_fraction"] = None
    out["disorder_reason"] = (
        f"NOT MEASURED ON THE FULL PROTEIN. No UniProt accession was supplied "
        f"and none could be read from any entry's _struct_ref, so the only "
        f"sequence available was the crystallised construct "
        f"({sequence_source}), {res.get('length')} residues. That is the "
        f"ordered part of the protein BY SELECTION, so its disorder fraction "
        f"({fraction}) is a lower bound on the construct and says nothing "
        f"about the rest of the chain. It is reported as "
        f"construct_disorder_fraction; disorder_fraction is null because a "
        f"number in that field would be read as the protein's."
    )
    out["_caveat"] = (
        "Construct-only measurement. Quote it as "
        "'disorder <value> over <n> residues of the crystallised construct "
        "(<source>), not the full-length protein', or supply "
        "uniprot_accession and re-run. Never quote it as the target's disorder "
        "fraction, and never read a construct 0.0 as 'this protein is ordered'."
    )
    return out


# ===========================================================================
# cryptic mechanism — apo vs holo, the rule-5 measurement
# ===========================================================================


# Core C-alpha RMSD above which two entries are not superposed and every number
# derived from the fit is measured in the wrong frame. THE SAME VALUE mdpocket
# uses (`MDPOCKET_MAX_ACCEPTABLE_RMSD_A`), deliberately: the two stages were able
# to return opposite verdicts on one pair in one payload because only one of them
# had a gate.
CRYPTIC_MAX_CORE_CA_RMSD_A = 5.0

# A fit on a handful of atoms has a low RMSD because it has nothing to
# disagree with. S1PR1's receptor was mapped onto a 25-residue peptide and
# fitted on FIVE equivalent C-alpha; RMSD alone would never have caught it.
#
# THIS IS NOW A CEILING ON THE FLOOR, NOT THE FLOOR. See
# `_min_fitted_ca_floor`: an absolute 20 was never exercised against a small
# target — every pair in the regression carried 162-476 equivalent C-alpha and
# 135-422 fitted — and moving the gated count from `n_equivalent_ca` to
# `n_fitted_ca` took it closer to biting, because auto_trim's `min_fit_fraction`
# is 0.5 and a 30-residue pair can present 15 fitted. A 30-residue peptide and a
# single small domain are real cases in this set (TL1A's entries are 111-270
# residues; interface partners run 25-63), so the floor scales with the smaller
# of the two MAPPED CHAINS below 40 residues and is 20 at or above it.
CRYPTIC_MIN_EQUIVALENT_CA = 20

# Below `CRYPTIC_MIN_EQUIVALENT_CA` the floor becomes this fraction of the
# SMALLER MAPPED CHAIN — not of `n_equivalent_ca`, which is the exploitable one:
# S1PR1's bad mapping produced 5 equivalent positions onto a 25-residue peptide,
# and scaling by that would have made 5 of 5 pass. Chain sizes are read from the
# two coordinate files and cannot be narrowed by the fit.
CRYPTIC_MIN_FITTED_CA_CHAIN_FRACTION = 0.5

# And a hard bottom, whatever the chains are. PROPOSED, NOT CALIBRATED, on a
# geometric argument rather than a measurement: three points determine a rigid
# body exactly, so a fit on a handful of C-alpha reports an RMSD near zero by
# construction and measures nothing. Eight positions is 24 coordinates against
# six degrees of freedom.
CRYPTIC_ABS_MIN_FITTED_CA = 8


def _ca_counts_by_chain(path) -> dict[str, int]:
    """{chain -> amino-acid C-alpha count} for a PDB **or mmCIF** file.

    THROUGH gemmi, NOT BY COLUMN. The two paths the cryptic gate hands this are
    the raw mmCIF entries (`cif_by_pid`, `donor["cif"]`), not the prepared PDB —
    a fixed-column reader returns `{}` on them, the floor silently falls back to
    the absolute 20, and the scaling this exists for never happens. Read with
    the same library `cryptic_analysis._load_structure` uses so the chain names
    match the `chain_mapping` they are looked up by.
    """
    import gemmi

    out: dict[str, int] = {}
    p = Path(path)
    if not p.exists():
        return out
    try:
        st = gemmi.read_structure(str(p))
        st.setup_entities()
        st.remove_waters()
    except Exception:  # noqa: BLE001
        return out
    for ch in st[0]:
        n = sum(
            1 for res in ch
            if res.find_atom("CA", "*") is not None
            and (gemmi.find_tabulated_residue(res.name) or None)
            and gemmi.find_tabulated_residue(res.name).is_amino_acid()
        )
        if n:
            out[ch.name] = n
    return out


def _smaller_mapped_chain_ca(holo_path, apo_path, mapping: dict | None) -> int | None:
    """C-alpha count of the SMALLER of the two chains that were superposed.

    `chain_mapping` is {holo_chain: apo_chain}. Sums each side over the mapped
    chains and returns the smaller total, so a multi-chain mapping is measured
    as the two assemblies actually fitted. None when it cannot be established,
    and None must fall back to the absolute floor rather than to no floor.
    """
    if not mapping:
        return None
    h_counts, a_counts = _ca_counts_by_chain(holo_path), _ca_counts_by_chain(apo_path)
    if not h_counts or not a_counts:
        return None
    h = sum(h_counts.get(c, 0) for c in mapping)
    a = sum(a_counts.get(c, 0) for c in mapping.values())
    if not h or not a:
        return None
    return min(h, a)


def _min_fitted_ca_floor(n_smaller_chain_ca: int | None) -> tuple[int, str]:
    """How many fitted C-alpha this pair must carry. Returns (floor, basis).

    A DEBT THAT WAS RECORDED AND IS NOW PAID. The count gate reads
    `n_fitted_ca`, which is right — scoring the floor on `n_equivalent_ca` let a
    narrowed fit clear every check while the pair sat 25.6 A out of frame — but
    an ABSOLUTE floor of 20 on the fitted count refuses a legitimately small
    target: 30 equivalent C-alpha through auto_trim's `min_fit_fraction` of 0.5
    is 15 fitted, and 15 < 20. That is a valid comparison on a real case (a
    30-residue peptide, a single small domain), refused by a threshold that has
    never been exercised anywhere near itself.

    So the floor is `min(20, max(8, 0.5 x smaller mapped chain))`:

        smaller chain   floor   effect
        >= 40           20      IDENTICAL to today. The whole regression set
                                (162-476 equivalent, 135-422 fitted) is here, so
                                this change cannot move a validated result.
        30              15      a 30-residue pair at the trim limit now PASSES
        25              12      S1PR1's CD69 mis-mapping still REFUSES at 5
        <= 16           8       the hard bottom takes over; a fit this small is
                                not a superposition whatever the chains are

    Scaled on the CHAIN, never on `n_equivalent_ca`: the equivalent count is
    itself a product of the mapping, and the S1PR1 failure was 5 equivalent onto
    a 25-residue peptide, which a self-referential floor would have waved
    through at 5 of 5.
    """
    if not n_smaller_chain_ca:
        return CRYPTIC_MIN_EQUIVALENT_CA, (
            f"absolute floor {CRYPTIC_MIN_EQUIVALENT_CA}: the mapped chains' "
            "sizes could not be read, so the floor is not scaled"
        )
    scaled = round(CRYPTIC_MIN_FITTED_CA_CHAIN_FRACTION * n_smaller_chain_ca)
    floor = min(CRYPTIC_MIN_EQUIVALENT_CA, max(CRYPTIC_ABS_MIN_FITTED_CA, scaled))
    if floor == CRYPTIC_MIN_EQUIVALENT_CA:
        basis = (
            f"absolute floor {CRYPTIC_MIN_EQUIVALENT_CA}; the smaller mapped "
            f"chain has {n_smaller_chain_ca} C-alpha, so the scaled floor "
            f"({scaled}) does not bind"
        )
    elif floor == CRYPTIC_ABS_MIN_FITTED_CA:
        basis = (
            f"hard bottom {CRYPTIC_ABS_MIN_FITTED_CA}: the smaller mapped chain "
            f"has only {n_smaller_chain_ca} C-alpha and "
            f"{CRYPTIC_MIN_FITTED_CA_CHAIN_FRACTION:.0%} of it ({scaled}) is "
            "below the point where a rigid-body fit is over-determined at all"
        )
    else:
        basis = (
            f"{CRYPTIC_MIN_FITTED_CA_CHAIN_FRACTION:.0%} of the smaller mapped "
            f"chain ({n_smaller_chain_ca} C-alpha) = {floor}; the absolute "
            f"{CRYPTIC_MIN_EQUIVALENT_CA} would refuse a legitimately small "
            "target that auto_trim has fitted at its own 0.5 limit"
        )
    return floor, basis

# Fraction of fitted positions allowed to name a different residue in the two
# entries. The same S1PR1 fit carried 15 name mismatches out of 5 positions'
# worth of signal. Construct differences (KRAS G12C/C51S/C80L/C118S, TNF L143D)
# are a handful; a tenth of the fit is a different sequence.
CRYPTIC_MAX_NAME_MISMATCH_FRACTION = 0.1

# `cryptic_analysis.analyze_cryptic_mechanism`'s own default for the
# displacement that separates loop_or_backbone_motion from sidechain_occlusion.
# Mirrored here ONLY to report the margin; this file never passes it and never
# overrides it, so if the default there moves this must move with it.
CRYPTIC_BACKBONE_MOTION_THRESHOLD_A = 2.0

# A mechanism label decided by less than this much displacement is a threshold
# crossing, not a measurement of a mechanism. PROPOSED, NOT CALIBRATED: it is
# set from the one case that produced the complaint (S1PR1 3V2Y -> 7TD4 clears
# 2.00 A by 0.16 A on one residue and is thereby labelled loop_or_backbone_motion
# with a nanomolar ceiling). It flags and never filters.
CRYPTIC_MECHANISM_MARGIN_A = 0.5

# A protein-wide C-alpha displacement at or above this is a rearrangement worth
# naming even when the site itself is still. PROPOSED, NOT CALIBRATED: set below
# S1PR1's 14.6 A TM6 swing and above the largest site-local motion in the
# calibration set (KRAS 8.65 A is a SITE motion and is reported as one).
CRYPTIC_GLOBAL_MOTION_NOTABLE_A = 5.0


def _global_ca_displacement(
    apo_path: Path, holo_path: Path, sup: dict
) -> dict:
    """Protein-wide maximum C-alpha displacement, beside the site-local one.

    THE PAYLOAD HAD NO PROTEIN-WIDE MOTION FIELD AT ALL, and one control proved
    that is a hole rather than a simplification. S1PR1 inactive 3V2Y against
    active 7TD4 passes the gate, returns is_cryptic false with a site C-alpha
    RMSD of 1.04 A — which is right, the site is pre-formed — and TM6's 14.6 A
    activation swing appears NOWHERE. `result["global"]` is null; the only trace
    is all_ca_rmsd_after_core_fit 2.035, which reads as "fine". A dossier built
    from that payload would state the site is pre-formed and never mention that
    the two conformers differ by an activation-state rearrangement.

    A large global motion beside a still site is a real, reportable state. It is
    not crypticity and it is not a failure; it is the thing a reader most needs
    to know about the pair, and it was invisible.

    Recomputed here rather than taken from `cryptic_analysis`, which reports
    per-residue displacement only over the SITE. The exact fit is reconstructed
    from the superposition block — same chain mapping, same name-matching
    setting, same excluded positions — and `reconstructed_core_ca_rmsd_a` is
    returned beside the reported `core_ca_rmsd` as a self-check: if those two
    disagree the reconstruction is wrong and the number must not be used.
    """
    out: dict = {
        "max_ca_displacement_a": None,
        "max_ca_displacement_at": None,
        "n_ca_compared": None,
        "reconstructed_core_ca_rmsd_a": None,
        "reconstruction_agrees": None,
        "error": None,
    }
    try:
        import numpy as np

        from cryptic_analysis import (
            _Atoms,
            _apply,
            _kabsch,
            _load_structure,
            _pair_coords,
        )

        mapping = list((sup.get("chain_mapping") or {}).items())
        if not mapping:
            out["error"] = "no chain_mapping in the superposition block"
            return out
        holo = _Atoms(_load_structure(str(holo_path)))
        apo = _Atoms(_load_structure(str(apo_path)))
        P, Q, keys, _mm = _pair_coords(
            holo, apo, mapping, bool(sup.get("match_residue_names", True))
        )
        if len(P) < 3:
            out["error"] = "fewer than 3 equivalent C-alpha to compare"
            return out
        excl = {
            (str(e.get("holo_chain")), int(e.get("resi")))
            for e in (sup.get("excluded_residues") or [])
            if e.get("resi") is not None
        }
        keep = np.array([(str(k[0]), int(k[2])) not in excl for k in keys])
        if keep.sum() < 3:
            keep = np.ones(len(P), dtype=bool)
            out["error"] = (
                "the fitted subset could not be reconstructed from "
                "excluded_residues; fitted on every equivalent position instead"
            )
        R, t = _kabsch(P[keep], Q[keep])
        dev = np.linalg.norm(_apply(R, t, P) - Q, axis=1)
        i = int(dev.argmax())
        recon = round(float(np.sqrt((dev[keep] ** 2).mean())), 3)
        reported = sup.get("core_ca_rmsd")
        out.update(
            max_ca_displacement_a=round(float(dev.max()), 2),
            max_ca_displacement_at=f"{keys[i][0]}:{keys[i][2]}",
            n_ca_compared=int(len(P)),
            n_ca_in_fit=int(keep.sum()),
            reconstructed_core_ca_rmsd_a=recon,
            reconstruction_agrees=(
                None if reported is None else abs(recon - float(reported)) <= 0.02
            ),
        )
    except Exception as exc:  # noqa: BLE001
        out["error"] = f"{type(exc).__name__}: {exc}"
    return out


def _mechanism_margin(r: dict) -> dict:
    """How close was the mechanism label to being the other label?

    `loop_or_backbone_motion` and `sidechain_occlusion` are separated by ONE
    comparison — max site C-alpha displacement against 2.0 A — and the two sides
    of it carry opposite potency priors under rule 5 (nanomolar against
    micromolar-at-best). A label produced 0.16 A from the boundary is not a
    measurement of a mechanism, and nothing in the payload said how far from the
    boundary any label was.

    `subunit_occlusion` is decided earlier and by a different quantity, and
    `none` by a third, so this block states which mechanisms it applies to
    rather than pretending to score all four.
    """
    site = r.get("site") or {}
    disp = site.get("max_ca_displacement")
    mech = r.get("mechanism")
    applies = mech in ("loop_or_backbone_motion", "sidechain_occlusion")
    margin = (
        round(float(disp) - CRYPTIC_BACKBONE_MOTION_THRESHOLD_A, 2)
        if disp is not None else None
    )
    narrow = (
        margin is not None and applies
        and abs(margin) < CRYPTIC_MECHANISM_MARGIN_A
    )
    return {
        "mechanism": mech,
        "applies_to_this_mechanism": applies,
        "deciding_quantity": "max site C-alpha displacement, one residue",
        "value_a": disp,
        "at": site.get("max_ca_displacement_at"),
        "threshold_a": CRYPTIC_BACKBONE_MOTION_THRESHOLD_A,
        "margin_a": margin if applies else None,
        "margin_warn_a": CRYPTIC_MECHANISM_MARGIN_A,
        "decided_by_a_narrow_margin": narrow if applies else None,
        "note": (
            None if not narrow else
            f"THE LABEL IS {mech} BY {abs(margin)} A ON ONE RESIDUE "
            f"({site.get('max_ca_displacement_at')}). Rule 5 turns that label "
            "into a potency prior, so a margin this small is the whole basis of "
            "the prior. Measured on S1PR1 3V2Y -> 7TD4: 2.16 A against 2.00 A, "
            "site C-alpha RMSD 1.04 A, is_cryptic false — the right answer "
            "reached by the wrong reasoning. Quote the displacement, not the "
            "label, when this is set."
        ),
        "_why": (
            "PROPOSED, NOT CALIBRATED. The warn band is set from the one case "
            "that produced the complaint and flags rather than filters: no "
            "label is changed, suppressed or re-derived here."
        ),
    }


def _motion_scope(r: dict, sup: dict, glob: dict) -> dict:
    """Site-local motion beside protein-wide motion, in one place.

    The pair (site still, protein moved) is a real and common state — an
    activation-state rearrangement with a pre-formed site — and until this block
    existed the payload could not express it. See `_global_ca_displacement`.
    """
    site = r.get("site") or {}
    site_max = site.get("max_ca_displacement")
    site_rmsd = site.get("ca_rmsd")
    gmax = glob.get("max_ca_displacement_a")
    # "STILL" IS TESTED ON THE SITE RMSD, NOT ON THE SITE MAXIMUM, and
    # deliberately so: the maximum is one residue and is exactly the knife-edge
    # quantity `mechanism_margin` exists to warn about. S1PR1 3V2Y -> 7TD4 has a
    # site C-alpha RMSD of 1.04 A and a maximum of 2.16 A, and keying on the
    # maximum would have declared the site NOT still by 0.16 A — hiding the very
    # case this block was added for. 2.0 A is also CryptoBench's own
    # pocket-residue RMSD criterion, so the same number means the same thing on
    # both sides.
    still_site = (
        site_rmsd is not None
        and float(site_rmsd) <= CRYPTIC_BACKBONE_MOTION_THRESHOLD_A
    )
    big_global = (
        gmax is not None and float(gmax) >= CRYPTIC_GLOBAL_MOTION_NOTABLE_A
    )
    return {
        "site_max_ca_displacement_a": site_max,
        "still_site_tested_on": "site_ca_rmsd_a",
        "still_site_threshold_a": CRYPTIC_BACKBONE_MOTION_THRESHOLD_A,
        "site_max_ca_displacement_at": site.get("max_ca_displacement_at"),
        "site_ca_rmsd_a": site.get("ca_rmsd"),
        "global_max_ca_displacement_a": gmax,
        "global_max_ca_displacement_at": glob.get("max_ca_displacement_at"),
        "global_n_ca_compared": glob.get("n_ca_compared"),
        "all_ca_rmsd_after_core_fit_a": sup.get("all_ca_rmsd_after_core_fit"),
        "global_notable_threshold_a": CRYPTIC_GLOBAL_MOTION_NOTABLE_A,
        "global_motion_with_still_site": (
            bool(big_global and still_site)
            if (gmax is not None and site_rmsd is not None) else None
        ),
        "note": (
            None if not (big_global and still_site) else
            f"THE SITE IS STILL (C-alpha RMSD {site_rmsd} A, maximum "
            f"{site_max} A) AND THE PROTEIN IS NOT ({gmax} A at "
            f"{glob.get('max_ca_displacement_at')}). This is a reportable "
            "state, not a defect and not crypticity: the two conformers differ "
            "by a rearrangement that does not pass through the site. S1PR1 "
            "inactive 3V2Y -> active 7TD4 is the case — site C-alpha RMSD "
            "1.04 A, is_cryptic false, and a 14.4 A TM6 swing at R:248 that had "
            "no field to appear in: result['global'] was null and the only "
            "protein-wide number was all_ca_rmsd_after_core_fit 2.035, which "
            "reads as agreement. Say BOTH in the dossier; 'the site is "
            "pre-formed' alone omits the activation-state change."
        ),
        "reconstruction": {
            k: glob.get(k) for k in (
                "reconstructed_core_ca_rmsd_a", "reconstruction_agrees",
                "n_ca_in_fit", "error",
            )
        },
        "_why": (
            "An RMSD is not a maximum. all_ca_rmsd_after_core_fit was the only "
            "protein-wide number in the payload and on S1PR1 it reads 2.035, "
            "which looks like agreement; the largest single displacement behind "
            "it is 14.6 A. Both are reported because they answer different "
            "questions."
        ),
    }


def _cryptic_block(
    apo_path: Path, holo_path: Path, comp_id: str,
    apo_pid: str, holo_pid: str,
    apo_chains: list[str] | None, ligand_chain: str | None,
) -> dict:
    """Classify why the site is not visible in this apo structure.

    `holo_chains` is deliberately NOT passed. The module infers it from ligand
    contacts with a 15%-of-top-contributor floor, and that floor is exactly what
    is needed here: 2AZ5's assembly is two independent TNF-alpha dimers and the
    chain-A ligand brushes 3 atoms of chain D against 44 and 39 for the real
    partners. Passing the four chains we scored would build an A/B/D
    "assembly", leave no apo chain free to be displaced, and silently turn
    subunit_occlusion into a confident loop_or_backbone_motion + cryptic:true on
    a target that is neither.

    `apo_chains` IS passed, and is the assembly this module actually scored:
    any apo chain not mapped to a holo chain is treated as a chain the ligand
    must displace, so leaving crystallographic extras in would invent subunits.

    Regression, AS THIS FUNCTION RUNS IT (the module's zero-knowledge default:
    auto-trim on, residue-name matching on):
        4OBE vs 6OIM (MOV) -> loop_or_backbone_motion, cryptic, ~8.65 A
        1TNF vs 2AZ5 (307) -> subunit_occlusion, NOT cryptic, ~1.55-1.58 A

    THE HAND-CALIBRATION FIGURES ARE 8.83 A AND 1.62 A AND THIS PROTOCOL DOES
    NOT REPRODUCE THEM. They came from a protocol that disabled both switches
    and named the mobile regions by hand; the default lands 0.1-0.2 A below
    them. Mechanism and is_cryptic are IDENTICAL under both, so nothing
    downstream of the label changes — but the displacement figures are not
    interchangeable and must not be quoted as each other. The calibration
    protocol is re-run below and reported in `calibration_protocol`.
    """
    out: dict = {
        "cryptic_status": "not_run",
        "cryptic_reason": None,
        "apo_pdb_id": apo_pid,
        "holo_pdb_id": holo_pid,
        "ligand_comp_id": comp_id,
    }
    try:
        from cryptic_analysis import analyze_cryptic_mechanism

        r = analyze_cryptic_mechanism(
            str(apo_path), str(holo_path), comp_id,
            apo_chains=apo_chains or None,
            ligand_chain=ligand_chain,
        )
    except Exception as exc:  # noqa: BLE001
        out.update(cryptic_status="failed",
                   cryptic_reason=f"{type(exc).__name__}: {exc}")
        return out

    sc = r.get("self_control", {})
    # ---- SUPERPOSITION-QUALITY GATE ---------------------------------------
    # The self-control alone is not enough, and a real run proved it. On NLRP3
    # the module returned cryptic_status "ok", is_cryptic true, 21.6 A
    # displacement, mechanism loop_or_backbone_motion and a NANOMOLAR potency
    # prior — on top of its own reported `core_ca_rmsd: 16.627` over 487 CA
    # with `n_excluded_ca: 0`, and all four chain mappings at 16.629. A 16.6 A
    # core RMSD is not a superposition. mdpocket REFUSED the identical pair in
    # the same run ("8SWF: best chain mapping RMSD 16.22 A exceeds 5.0 A"), so
    # one module rejected the alignment while this one built a confident
    # mechanistic call on top of it.
    #
    # The same pair against a different apo entry (7ZGU) superposed at 1.248 A
    # and gave the OPPOSITE answer: 0.95 A displacement, mechanism none,
    # is_cryptic false. Without a manual control the dossier would have called a
    # validated, clinically drugged, open ATP site cryptic with a 21.6 A
    # conformational change.
    #
    # Same threshold as mdpocket's, for the same reason and so that the two
    # stages cannot disagree about whether a pair is superposable.
    # FOUR GATES, NOT ONE — and the fourth is below, at `sup_all`, because three
    # were not enough either. RMSD alone would not have caught S1PR1, which is
    # the worst instance measured: the module mapped the S1PR1 receptor onto
    # 8G94 chain F — CD69, a 25-residue peptide — fitted it on FIVE equivalent
    # C-alpha with FIFTEEN residue-name mismatches, and reported the resulting
    # `max_backbone_ca_displacement_a: 0.00` as a measurement with
    # `mechanism: subunit_occlusion` and status ok. Rule 5 then maps
    # subunit_occlusion to a micromolar-at-best ceiling — on a target with 600
    # sub-nanomolar compounds and five approved drugs. Four log units wrong, in
    # the damaging direction, from a block that carried every diagnostic of its
    # own failure and gated on none of them. Hand re-measurement over 257 core
    # C-alpha at 1.03 A gives 1.33 A displacement, zero clashes, mechanism none.
    sup = r.get("superposition") or {}
    sup_rmsd = sup.get("core_ca_rmsd")
    # THE FITTED SUBSET IS NOT THE PAIR, AND THE GATE READ ONLY THE SUBSET.
    # `core_ca_rmsd` is scored over the positions the fit actually stood on;
    # `all_ca_rmsd_after_core_fit` applies the SAME rotation to every equivalent
    # C-alpha. Narrowing the fit drives the first down and leaves the second
    # where it was, so a gate on the first alone is defeated by any exclusion.
    # Demonstrated on 8SWF vs 9HG4 with fit_residue_range=(130,370):
    #
    #     core_ca_rmsd 1.472 over n_fitted_ca 202   ->  the old gate PASSED
    #     n_excluded_ca 274, n_equivalent_ca 476, 0 name mismatches
    #     all_ca_rmsd_after_core_fit 25.619        <- in the same block
    #     emitted: 41.7 A, is_cryptic true, loop_or_backbone_motion, nanomolar
    #
    # That is a WORSE confident answer than the 21.6 A this gate was built to
    # stop, produced from a fit the payload itself scores at 25.6 A. The field
    # that catches it was already computed and already in the output; the gate
    # simply did not read it. Both numbers are read now, and the pair-wide one
    # is the one that cannot be narrowed away.
    sup_all = sup.get("all_ca_rmsd_after_core_fit")
    n_fitted = sup.get("n_fitted_ca")
    n_equivalent = sup.get("n_equivalent_ca")
    # `n_fitted_ca`, NOT `n_equivalent_ca`: the count that has to clear the floor
    # is the number of positions the fit stood on, not the number that were
    # equivalent BEFORE the exclusions. Tested against None rather than written
    # as `n_fitted_ca or n_equivalent_ca`, because a legitimate zero is falsy and
    # would fall straight through to the pre-fit count — silently gating the
    # wrong field in exactly the case (an empty fit) the floor exists for.
    n_equiv = n_fitted if n_fitted is not None else n_equivalent
    n_mismatch = sup.get("n_residue_name_mismatches") or 0
    # A FRACTION NEEDS ITS OWN DENOMINATOR. `n_residue_name_mismatches` is
    # counted over the pre-filter overlap while `n_equivalent_ca` counts the
    # SURVIVORS of that filter, so mismatches/equivalent is not a fraction at
    # all — on the S1PR1 CD69 fit it printed "15 of 5 fitted positions", a ratio
    # of 3.0. The denominator is the overlap the mismatches were counted over,
    # which is survivors + mismatches while `match_residue_names` is on. With it
    # off nothing is dropped and the survivors ARE the overlap.
    n_overlap = (
        (n_equivalent or 0) + n_mismatch
        if sup.get("match_residue_names", True)
        else (n_equivalent or 0)
    )
    gate_fails: list[str] = []
    core_failed = (
        sup_rmsd is not None and sup_rmsd > CRYPTIC_MAX_CORE_CA_RMSD_A
    )
    if core_failed:
        gate_fails.append(
            f"core C-alpha RMSD {sup_rmsd} A over the {n_fitted} fitted "
            f"positions exceeds {CRYPTIC_MAX_CORE_CA_RMSD_A} A"
        )
    # Only when the core check did NOT already fire. With nothing excluded the
    # two numbers are the same number, and printing "the fit describes a
    # fragment" over an unnarrowed fit would misdiagnose the refusal exactly the
    # way this patch exists to stop. The value is in the gate block either way.
    if (
        not core_failed
        and sup_all is not None
        and sup_all > CRYPTIC_MAX_CORE_CA_RMSD_A
    ):
        gate_fails.append(
            f"after the core fit, RMSD over ALL {n_equivalent} equivalent "
            f"C-alpha is {sup_all} A (the fitted subset was {sup_rmsd} A over "
            f"{n_fitted}, with {sup.get('n_excluded_ca')} excluded); the fit "
            "describes a fragment, not the pair"
        )
    # THE FLOOR SCALES WITH THE SMALLER MAPPED CHAIN. An absolute 20 on the
    # FITTED count refuses a legitimately small target — 30 equivalent C-alpha
    # at auto_trim's own 0.5 limit is 15 fitted — and it has never been
    # exercised against one. See `_min_fitted_ca_floor`; at 40 residues and
    # above, which is the entire regression set, nothing changes.
    n_smaller_chain = _smaller_mapped_chain_ca(
        holo_path, apo_path, sup.get("chain_mapping")
    )
    ca_floor, ca_floor_basis = _min_fitted_ca_floor(n_smaller_chain)
    if n_equiv is not None and n_equiv < ca_floor:
        gate_fails.append(
            f"only {n_equiv} C-alpha were fitted ({ca_floor} needed — "
            f"{ca_floor_basis}); this is a fit onto the wrong chain, not a "
            "superposition"
        )
    if n_overlap and n_mismatch / n_overlap > CRYPTIC_MAX_NAME_MISMATCH_FRACTION:
        gate_fails.append(
            f"{n_mismatch} of {n_overlap} equivalent positions name a DIFFERENT "
            "residue in the two entries; they are not the same sequence"
        )
    fit_ok = not gate_fails
    out.update(
        # `self_control.passed` False means the superposition or the ligand
        # placement is broken and EVERY number below is noise. Reported as a
        # failure rather than as a result, because a mechanism label from a
        # broken fit is worse than no label.
        cryptic_status="ok" if (sc.get("passed") and fit_ok) else "failed",
        cryptic_reason=(
            None if (sc.get("passed") and fit_ok)
            else "cryptic_analysis self-control failed; result not interpretable"
            if not sc.get("passed")
            else (
                "SUPERPOSITION REFUSED: " + "; ".join(gate_fails)
                + ". Every number below the fit — the displacement, the clash "
                "attribution, the free volume and the mechanism label — was "
                "measured in this frame, so none of them is interpretable. "
                "Refusing. WHAT THIS DOES NOT SAY IS WHY, AND THE TWO REASONS "
                "ARE NOT THE SAME FINDING. A large core RMSD has (a) mismatch "
                "causes — a mis-mapped chain, a different protein, a peptide "
                "fragment, a numbering offset — and (b) HINGE causes, where the "
                "two entries are the same protein in two rigid-body states and "
                "each domain superposes perfectly on its own. THIS GATE CANNOT "
                "TELL THEM APART and must not be read as saying (a). Test (b) "
                "before discarding the pair: fit each domain separately and see "
                "whether the RMSD collapses. Measured on NLRP3 8SWF vs 9HG4 — "
                "16.503 A over 476 equivalent C-alpha with ZERO positions "
                "trimmed, and 16.507 A restricted to a single apo chain, so "
                "neither a bad chain mapping nor a trimmed outlier explains it "
                "— yet the same pair restricted to the NBD (130-370) fits at "
                "1.472 A over 202 C-alpha, to 220-370 at 1.377 A over 151, and "
                "to HD1 (371-430) at 0.922 A over 58. An order of magnitude, "
                "per domain. That is a genuine NACHT hinge rotation between an "
                "open octamer and a closed NACHT, not a broken alignment; "
                "calling it 'not superposed' misdiagnoses it exactly the way "
                "the TL1A numbering-offset message used to. A large core RMSD "
                "WITH well-superposing subdomains is a hinge and is reportable "
                "as one. THERE IS NO FLAG FOR THE FIX: the domain-restricted fit "
                "that would recover a hinged pair is cryptic_analysis's "
                "fit_residue_range / exclude_residues, and pocket_scan does not "
                "expose either — see superposition_gate._why for why exposing "
                "them is not free. So a hinge is currently refused rather than "
                "measured, and that is a known gap, not a verdict. What IS "
                "reachable: try a different apo entry, or pass apo_chains. On "
                "NLRP3 the rejected 8SWF pair gave 21.6 A / cryptic and the "
                "superposable 7ZGU pair gave 0.95 A / not cryptic."
            )
        ),
        superposition_gate={
            "core_ca_rmsd_a": sup_rmsd,
            "all_ca_rmsd_after_core_fit_a": sup_all,
            "max_acceptable_a": CRYPTIC_MAX_CORE_CA_RMSD_A,
            "n_fitted_ca": n_fitted,
            "n_excluded_ca": sup.get("n_excluded_ca"),
            "n_equivalent_ca": n_equivalent,
            "n_ca_gated": n_equiv,
            "n_smaller_mapped_chain_ca": n_smaller_chain,
            "min_fitted_ca": ca_floor,
            "min_fitted_ca_basis": ca_floor_basis,
            "min_fitted_ca_ceiling": CRYPTIC_MIN_EQUIVALENT_CA,
            "min_fitted_ca_chain_fraction": CRYPTIC_MIN_FITTED_CA_CHAIN_FRACTION,
            "min_fitted_ca_hard_bottom": CRYPTIC_ABS_MIN_FITTED_CA,
            # Back-compatible alias. The key used to name the constant; it now
            # names the floor that was actually applied, which is the number a
            # reader of a refusal needs.
            "min_equivalent_ca": ca_floor,
            "min_equivalent_ca_status": (
                "SCALED, NOT ABSOLUTE, AS OF THIS VERSION. The floor is "
                "min(20, max(8, 0.5 x the smaller MAPPED CHAIN's C-alpha)). It "
                "was an absolute 20 and it was never exercised against a small "
                "target: every pair in the regression carried 162-476 "
                "equivalent C-alpha and 135-422 fitted, so the threshold had "
                "only ever been tested far away from itself, and moving the "
                "gated count to n_fitted_ca (correctly — a narrowed fit could "
                "otherwise turn every check green while the pair sat 25.6 A out "
                "of frame) took it CLOSER to biting, since auto_trim's "
                "min_fit_fraction of 0.5 lets a 30-residue pair present 15 "
                "fitted. A 30-residue peptide and a single small domain are "
                "real cases here — TL1A's entries run 111-270 residues and "
                "interface partners 25-63 — so refusing them was wrong. At 40 "
                "residues and above the floor is still exactly 20 and nothing "
                "in the regression moves. Scaling is on the CHAIN and never on "
                "n_equivalent_ca, because the equivalent count is a product of "
                "the mapping: S1PR1's receptor mapped onto a 25-residue peptide "
                "gave FIVE equivalent positions, and a self-referential floor "
                "would have passed 5 of 5. Scaled on the chain it needs 13 and "
                "still refuses. The 0.5 fraction and the hard bottom of 8 are "
                "PROPOSED, NOT CALIBRATED; the bottom rests on a geometric "
                "argument (three points determine a rigid body exactly, so a "
                "handful of C-alpha reports a near-zero RMSD by construction), "
                "not on a measurement."
            ),
            "n_residue_name_mismatches": n_mismatch,
            "name_mismatch_denominator": n_overlap,
            "name_mismatch_fraction": (
                round(n_mismatch / n_overlap, 3) if n_overlap else None
            ),
            "max_name_mismatch_fraction": CRYPTIC_MAX_NAME_MISMATCH_FRACTION,
            "chain_mapping": sup.get("chain_mapping"),
            "passed": fit_ok,
            "failures": gate_fails,
            "self_control_passed": sc.get("passed"),
            "_why": (
                "mdpocket applies a superposition gate and this stage did not, "
                "so the two stages could return opposite verdicts on the same "
                "pair in one payload. Three targets found it: NLRP3 (16.6 A "
                "core RMSD reported as ok), TL1A (numbering offsets), S1PR1 "
                "(5 equivalent C-alpha onto a 25-residue peptide, 15 name "
                "mismatches). Count and identity are gated, not only RMSD — "
                "the S1PR1 fit had a low RMSD precisely because it was fitted "
                "on five atoms."
            ),
            "_why_all_ca": (
                "FOUR CHECKS, NOT THREE, AND THE FOURTH IS THE ONE THAT CANNOT "
                "BE GAMED. core_ca_rmsd is scored over the FITTED SUBSET, so "
                "narrowing the fit lowers it without moving the structures: "
                "8SWF vs 9HG4 restricted to residues 130-370 gives core 1.472 A "
                "over 202 fitted C-alpha with 274 excluded and 0 name "
                "mismatches — three green lights — beside "
                "all_ca_rmsd_after_core_fit 25.619 A in the same block, and the "
                "old gate emitted 41.7 A / is_cryptic true / "
                "loop_or_backbone_motion / nanomolar prior on top of it. "
                "all_ca_rmsd_after_core_fit is the same rotation scored over "
                "every equivalent C-alpha, so it cannot be narrowed away and it "
                "is gated at the same threshold. THE ROUTE IS CURRENTLY LATENT: "
                "reaching it needs fit_residue_range or exclude_residues and "
                "pocket_scan exposes neither. It is not hypothetical, though — "
                "auto_trim writes the same n_excluded_ca by itself, so on a "
                "hinged protein this is one convergence away from happening "
                "unasked. The measured core-to-all gap over six real pairs is "
                "small (0.58->1.42, 0.96->1.10, 0.71->1.02, 1.30->2.72, "
                "1.27->2.04), which is why it had not fired yet, and is also "
                "why the check costs the controls nothing."
            ),
        },
        mechanism=r.get("mechanism"),
        secondary_mechanism=r.get("secondary_mechanism"),
        is_cryptic=r.get("is_cryptic"),
        rationale=r.get("rationale"),
        max_backbone_ca_displacement_a=(r.get("site") or {}).get(
            "max_ca_displacement"),
        max_ca_displacement_at=(r.get("site") or {}).get("max_ca_displacement_at"),
        site_ca_rmsd_a=(r.get("site") or {}).get("ca_rmsd"),
        # THE LABEL IS A THRESHOLD CROSSING AND A THRESHOLD CROSSING HAS A
        # MARGIN. S1PR1 inactive 3V2Y -> active 7TD4 comes back
        # `loop_or_backbone_motion` because ONE site residue moves 2.16 A
        # against a 2.00 A threshold — 0.16 A — and rule 5 maps that label to a
        # nanomolar ceiling. The verdict happens to be defensible and the
        # REASONING is not: it was generated by a 0.16 A margin on a single
        # residue, not by the motion anyone would describe. Reported so the
        # margin can never be read off the label alone.
        mechanism_margin=_mechanism_margin(r),
        # A LARGE GLOBAL MOTION WITH A STILL SITE IS A REAL STATE AND IT WAS
        # INVISIBLE. See `_global_ca_displacement`.
        motion_scope=_motion_scope(
            r, sup, _global_ca_displacement(apo_path, holo_path, sup)
        ),
        # CLASSIFY ON DISPLACEMENT, NOT ON CLASH COMPOSITION. Reported because
        # it is informative; it must not drive the label. KRAS's switch-II loop
        # moves 8.8 A and yet zero of the clashing atoms at 2.0 A are backbone —
        # keying on composition would hand the canonical nanomolar target a
        # micromolar prognosis.
        clash_attribution=r.get("contacts"),
        clash_attribution_wide=r.get("contacts_wide"),
        free_volume=r.get("free_volume"),
        crypticity=r.get("crypticity"),
        superposition={
            k: v for k, v in (r.get("superposition") or {}).items()
            if k != "per_residue"
        },
        displaced_apo_chains=(r.get("inputs") or {}).get("displaced_apo_chains"),
        bystander_apo_chains=(r.get("inputs") or {}).get("bystander_apo_chains"),
        holo_chains_used=(r.get("inputs") or {}).get("holo_chains"),
        self_control=sc,
        warnings=r.get("warnings"),
        cryptic_potency_prior=_potency_prior(r.get("mechanism")),
    )

    # SECOND PROTOCOL, reported alongside. The run above uses the module's
    # zero-knowledge default — auto-trim finds mobile regions nobody named, and
    # residue-name matching drops construct differences (KRAS G12C/C51S/C80L/
    # C118S, TNF L143D) out of the fit. The hand calibrations did neither, and
    # the two protocols disagree slightly on the DISPLACEMENT while agreeing
    # exactly on the mechanism and on is_cryptic:
    #
    #     KRAS   default 8.65 A   calibration 8.79 A   (hand figure 8.83 A,
    #                                                   which also excluded
    #                                                   switch I 25-40 and
    #                                                   switch II 57-75 and
    #                                                   fitted 1-166)
    #     TNF    default 1.55-1.58 A                   (hand figure 1.62 A)
    #
    # Both are shown because the default is the right production protocol and
    # the calibration number is the one the dossier's rule 5 quotes. Neither is
    # allowed to be presented as the other. The default runs 0.1-0.2 A BELOW the
    # hand figures on both targets, so a dossier that quotes 8.83 or 1.62 as
    # "what the pipeline measured" is quoting a protocol this function does not
    # run. Quote `calibration_protocol.max_backbone_ca_displacement_a` if the
    # calibration number is what is wanted, and say which protocol produced it.
    try:
        from cryptic_analysis import analyze_cryptic_mechanism as _acm

        cal = _acm(
            str(apo_path), str(holo_path), comp_id,
            apo_chains=apo_chains or None, ligand_chain=ligand_chain,
            auto_trim=False, match_residue_names=False,
            compute_free_volume=False,
        )
        out["calibration_protocol"] = {
            "_protocol": "auto_trim=False, match_residue_names=False",
            "mechanism": cal.get("mechanism"),
            "is_cryptic": cal.get("is_cryptic"),
            "max_backbone_ca_displacement_a": (cal.get("site") or {}).get(
                "max_ca_displacement"),
            "site_ca_rmsd_a": (cal.get("site") or {}).get("ca_rmsd"),
            "agrees_with_default": (
                cal.get("mechanism") == out.get("mechanism")
                and cal.get("is_cryptic") == out.get("is_cryptic")
            ),
        }
    except Exception as exc:  # noqa: BLE001
        out["calibration_protocol"] = {"error": f"{type(exc).__name__}: {exc}"}
    return out


def _potency_prior(mechanism: str | None) -> dict:
    """Mechanism is a prior on achievable potency, not a taxonomy label.

    CryptoSite set (Lazou, Kozakov, Joseph-McCarthy & Vajda, Drug Discov Today
    2024): of 27 loop-motion sites all but two reached nanomolar; of 18
    side-chain-motion sites only 10 had affinity data at all and every one of
    those bound weakly, low micromolar at best. The timescale argument is why —
    side chains reorient in 1e-11 to 1e-10 s and compete with the ligand, loops
    move in 1e-9 to 1e-6 s and can be wedged open.
    """
    table = {
        "loop_or_backbone_motion": (
            "nanomolar",
            "25 of 27 loop-motion sites in the CryptoSite set reached nanomolar",
        ),
        "sidechain_occlusion": (
            "micromolar_at_best",
            "every side-chain-occluded site in the CryptoSite set with affinity "
            "data bound in the low micromolar at best (10 of 18 had any data)",
        ),
        "subunit_occlusion": (
            "micromolar_at_best",
            "occlusion by a displaced subunit is not a loop opening; the ligand "
            "competes with a folded protein surface. TNF-alpha's SPD304 is the "
            "case and it is micromolar",
        ),
    }
    ceiling, basis = table.get(mechanism or "", ("unknown", None))
    return {"expected_ceiling": ceiling, "basis": basis}


# ===========================================================================
# interface classification — orthosteric / allosteric / destabiliser
# ===========================================================================


def _kmers(seq: str, k: int = 5) -> set:
    return {seq[i:i + k] for i in range(len(seq) - k + 1)} if len(seq) >= k else set()


def _split_partner_chains(partner_st, target_seqs: list[str]) -> tuple[list, list]:
    """Which chains of a complex are OUR protein and which are the partner.

    By 5-mer set overlap against the target's own sequences, not by chain
    letter and not by an annotation. 3ALQ's TNF-alpha chains are A,B,C and its
    TNFR2 chains are R,S,T, but nothing in the file says which is which in a
    way that generalises, and getting it backwards computes the epitope on the
    wrong side of the interface.
    """
    ref = set()
    for s in target_seqs:
        ref |= _kmers(s)
    ours, theirs = [], []
    for chain in partner_st[0]:
        names = [r.name for r in chain]
        if len(names) < 20:
            continue
        import gemmi

        seq = "".join(c for c in gemmi.one_letter_code(names).upper() if c.isalpha())
        km = _kmers(seq)
        share = (len(km & ref) / len(km)) if km else 0.0
        (ours if share >= 0.5 else theirs).append(chain.name)
    return ours, theirs


def _interface_block(
    partner_cifs: dict[str, Path], target_seqs: list[str],
    target_accession: str | None = None, work: Path | None = None,
) -> dict:
    """Partner epitope on the target side, from a real complex structure.

    WHICH CHAINS ARE OURS IS RESOLVED BY ACCESSION FIRST, SEQUENCE SECOND.
    Splitting by 5-mer overlap against `target_seqs` is only as good as
    `target_seqs`, and `target_seqs` came from `_one_letter` over the target
    chains — so when the chain resolver failed open, the "target sequence" was a
    partner's and the split inherited the error silently.

    Measured on IL-13 3BPO: `target_chains ["A","C"]` against `partner ["B"]`,
    where A is IL-13 and **C is IL-13R-alpha-1 at 314 aa**. `side_a` came back
    containing `C:ASN240, C:PHE259, C:TYR276` — receptor residues reported as
    the target's own epitope — with `interface_status: "ok"` and no warning.
    5E4E splits identically. It did not fire on BAFF only because BAFF's
    receptor fragments are 31-63 aa, so BAFF is the longest chain anyway: the
    length heuristic was still deciding, it was just guessing right.

    So the partner entry's own `_struct_ref_seq` is read and the split is made
    on it whenever it resolves. The sequence split stays as the fallback for
    entries that declare nothing, and where BOTH are available and DISAGREE the
    disagreement is reported rather than resolved.
    """
    out: dict = {
        "interface_status": "not_run",
        "interface_reason": None,
        "partner_structures": sorted(partner_cifs),
    }
    if not partner_cifs:
        out["interface_reason"] = (
            "no partner_structures supplied; pocket classification needs a "
            "complex containing the binding partner"
        )
        return out
    try:
        import interface_analysis as IA
    except Exception as exc:  # noqa: BLE001
        out.update(interface_status="failed",
                   interface_reason=f"interface_analysis import: {exc}")
        return out

    chosen = None
    per_partner: dict[str, dict] = {}
    for pid, cif in sorted(partner_cifs.items()):
        try:
            st = IA.load_structure(str(cif))
            all_ch = [c.name for c in st[0]]
            by_seq_ours, by_seq_theirs = _split_partner_chains(st, target_seqs)
            # ---- accession first ---------------------------------------
            acc_map, acc_status = (
                _chain_accessions(_fetch_header(pid, work), {}, all_ch)
                if (target_accession and work is not None) else ({}, "not_attempted")
            )
            # Matched through UniProt's merge/gene history, not literally — the
            # same rule `_target_chains` uses. Matching on the string here made
            # IL-13's 3BPO fall back to the sequence heuristic and report
            # `verified: false`, because 3BPO declares Q4VB50 (an unreviewed
            # entry for the same IL13 gene) rather than P35225.
            acc_ours = []
            for c in all_ch:
                for a in (
                    acc_map.get(c) or acc_map.get(_assembly_base_chain(c)) or ()
                ):
                    if _accession_matches(a, target_accession)[0]:
                        acc_ours.append(c)
                        break
            if acc_ours:
                ours = acc_ours
                theirs = [c for c in all_ch if c not in acc_ours]
                split_basis = (
                    f"chains mapping to {target_accession} in {pid}'s own "
                    "_struct_ref_seq"
                )
            else:
                ours, theirs = by_seq_ours, by_seq_theirs
                split_basis = (
                    "5-mer sequence overlap against the target's sequence "
                    f"(accession split unavailable: {acc_status}"
                    + ("" if acc_status != "ok" else
                       f"; {pid} declares no chain for {target_accession}")
                    + "). UNVERIFIED — this is the heuristic that put "
                    "IL-13R-alpha-1 on IL-13's side of its own interface."
                )
            disagreement = (
                None if not acc_ours or sorted(by_seq_ours) == sorted(acc_ours)
                else (
                    f"the accession split gives target chains {sorted(acc_ours)} "
                    f"and the 5-mer sequence split gives {sorted(by_seq_ours)}. "
                    "The accession is used. A chain the sequence split claimed "
                    "as target and the accession does not is a partner subunit "
                    "similar enough to fool a 5-mer overlap — which is exactly "
                    "the IL-13/IL-13R-alpha-1 failure."
                )
            )
            if not ours or not theirs:
                per_partner[pid] = {
                    "error": (
                        f"could not split {pid} into target and partner chains "
                        f"(target-like {ours}, other {theirs}); it may not be a "
                        "complex, or it may be a homo-oligomer of the target only"
                    ),
                    "target_partner_split_basis": split_basis,
                }
                continue
            iface = IA.interface_residues(st, ours, theirs)
            rec = {
                "target_chains": ours,
                "partner_chains": theirs,
                "target_partner_split_basis": split_basis,
                "target_partner_split_verified": bool(acc_ours),
                "target_partner_split_disagreement": disagreement,
                "chain_accessions": acc_map,
                **iface.as_dict(),
            }
            per_partner[pid] = rec
            if chosen is None or len(iface.side_a) > len(chosen[1].side_a):
                chosen = (pid, iface, st, ours)
        except Exception as exc:  # noqa: BLE001
            per_partner[pid] = {"error": f"{type(exc).__name__}: {exc}"}
    out["per_partner"] = per_partner
    if chosen is None:
        out.update(interface_status="failed",
                   interface_reason="no partner structure yielded an interface")
        return out
    pid, iface, _st, _ours = chosen
    # THE PARTNER ENTRY'S NUMBERING, so a seqid match can be checked instead of
    # assumed. `classify_pocket(match_by="seqid")` matches pocket residues to
    # interface residues on number alone; whether that is legal depends on the
    # two entries numbering the protein the same way, and two of five TL1A
    # structures did not. Pooled over the target-side chains only.
    resnames: dict[int, str] = {}
    for chain in _st[0]:
        if chain.name not in _ours:
            continue
        for res in chain:
            if res.het_flag == "A" and len(res):
                resnames.setdefault(res.seqid.num, res.name)
    out.update(
        interface_status="ok",
        partner_pdb_id=pid,
        n_interface_residues=len(iface.side_a),
        interface_residues=[r.label for r in iface.side_a],
        partner_target_chains=list(_ours),
        _epitope=iface.side_a,  # popped before return; not JSON-serialisable
        _partner_resnames=resnames,  # popped before return
    )
    return out


# ---------------------------------------------------------------------------
# Buried core detection. EVERY THRESHOLD BELOW IS A PROPOSAL, NOT CALIBRATED.
# ---------------------------------------------------------------------------
#
# fpocket's shipped druggability score is a logistic regression on mean local
# hydrophobic density, max alpha-sphere distance and polar VDW surface. That
# combination is maximised by a large, sealed, hydrophobic void — which is a
# precise description of the hydrophobic CORE of a folded domain, and cores are
# not binding sites. They have no solvent mouth for a ligand to enter through,
# and opening one costs the fold.
#
# Measured on IRAK4's death domain in a full two-mode run: the top-ranked pocket
# of 134 scored druggability 0.890 with enclosure 0.998 (sealed — no mouth),
# subunit_enclosure_gain 0.020 (partner chains contribute nothing to the burial,
# so it is buried within ONE chain) and interface_coverage 0.026. Its lining was
# nine Leu/Ile/Val/Phe, one Arg, one Tyr. Every supporting field said "core";
# only the headline number said "site", and the headline number is the one that
# gets quoted.
#
# THESE NUMBERS ARE A PROPOSAL FROM ONE CASE. They are not calibrated, there is
# no held-out set behind them, and they are placed just inside the one
# observation we have. They gate a FLAG and never a filter: no pocket is dropped,
# reordered or rescored because of them. A buried core that is real is still
# returned, still ranked where fpocket put it, and still carries its score.
BURIED_CORE_ENCLOSURE_MIN = 0.98
BURIED_CORE_SUBUNIT_GAIN_MAX = 0.05
BURIED_CORE_APOLAR_FRACTION_MIN = 0.7


def _buried_core_flag(
    enclosure: float | None,
    subunit_gain: float | None,
    apolar_fraction: float | None,
) -> dict:
    """Is this "pocket" the hydrophobic core of a domain rather than a site?

    Keyed on GEOMETRY (near-total enclosure with no solvent mouth, and burial
    contributed by a single subunit) plus lining composition, never on the
    druggability score itself — the score is the thing being questioned.

    Returns a dict that is always present, so "not flagged" and "could not be
    measured" are distinguishable.
    """
    have = [v for v in (enclosure, subunit_gain) if v is not None]
    if not have:
        return {
            "buried_core_suspected": None,
            "buried_core_reason": (
                "enclosure was not measured for this pocket, so the core test "
                "could not be applied. It needs the interface stage, which "
                "needs a partner structure."
            ),
            "buried_core_criteria": None,
        }
    sealed = enclosure is not None and enclosure >= BURIED_CORE_ENCLOSURE_MIN
    single = subunit_gain is not None and subunit_gain <= BURIED_CORE_SUBUNIT_GAIN_MAX
    greasy = (
        apolar_fraction is not None
        and apolar_fraction >= BURIED_CORE_APOLAR_FRACTION_MIN
    )
    hit = bool(sealed and single and greasy)
    return {
        "buried_core_suspected": hit,
        "buried_core_reason": (
            (
                f"enclosure {enclosure} >= {BURIED_CORE_ENCLOSURE_MIN} (sealed, "
                f"no solvent mouth), subunit_enclosure_gain {subunit_gain} <= "
                f"{BURIED_CORE_SUBUNIT_GAIN_MAX} (buried within one subunit, not "
                f"by the assembly), apolar lining fraction {apolar_fraction} >= "
                f"{BURIED_CORE_APOLAR_FRACTION_MIN}. This is the shape of a "
                "hydrophobic core, and fpocket's druggability regression rewards "
                "precisely that shape. TREAT THE DRUGGABILITY SCORE ON THIS "
                "POCKET AS UNINTERPRETABLE, not as a high one."
            )
            if hit
            else None
        ),
        # THE MEASURED VALUES ONLY. The cut points and their PROPOSED/NOT
        # CALIBRATED status are constants and used to be re-emitted inside every
        # classified pocket — ~500 characters x every pocket x every clustering
        # value x every structure, which is a payload cap spent on repeating
        # three numbers. They are stated once, at
        # `pocket_vs_interface._buried_core_thresholds`.
        "buried_core_criteria": {
            "enclosure": enclosure,
            "subunit_enclosure_gain": subunit_gain,
            "apolar_lining_fraction": apolar_fraction,
            "sealed": sealed,
            "buried_within_one_subunit": single,
            "hydrophobic_lining": greasy,
        },
    }


BURIED_CORE_THRESHOLDS_NOTE = {
    "_status": "PROPOSED, NOT CALIBRATED — from a single observed case (IRAK4 "
               "death domain, druggability 0.890 at rank 1 of 134, enclosure "
               "0.998, gain 0.020). Do not present these cut points as "
               "measured. They gate a flag, never a filter.",
    "enclosure_min": BURIED_CORE_ENCLOSURE_MIN,
    "subunit_enclosure_gain_max": BURIED_CORE_SUBUNIT_GAIN_MAX,
    "apolar_lining_fraction_min": BURIED_CORE_APOLAR_FRACTION_MIN,
}


# Below this fraction of residue-name agreement, two entries are not on a common
# residue numbering and a seqid-keyed comparison between them is meaningless.
# PROPOSED, NOT CALIBRATED.
NUMBERING_IDENTITY_MIN = 0.9

# A recovered numbering offset is only applied when it produces agreement over
# at least this many shared positions. Below it, "agreement" is a handful of
# residues lining up by chance and the shift would buy a legal-looking
# comparison rather than a correct one. PROPOSED, NOT CALIBRATED; the same
# number as the superposition floor, for the same reason.
NUMBERING_MIN_COMPARED_FOR_OFFSET = 20


def _numbering_agreement(
    a: dict[int, str], b: dict[int, str], max_examples: int = 8
) -> dict:
    """Do two entries agree on what residue each number refers to?

    THE CHECK THAT `match_by="seqid"` NEVER HAD. Matching pocket residues to
    interface residues on residue NUMBER alone is only valid if both entries
    number the same protein the same way, and PDB entries for one protein
    routinely do not. Measured on TL1A: 2O0O numbers from the construct start
    while 3K51 runs at +67, so 2O0O's "shared A:HIS118" was matched against
    3K51's THR118 — non-homologous — and 2RE9's THR34/PRO35/THR36 were matched
    against VAL34/VAL35/ARG36, producing a spurious overlap_fraction of 0.227
    that was reported with a `borderline` flag and no numbering warning at all.

    Two of five structures were silently wrong and nothing in the payload said
    which two. Residue names cost one comparison and catch all of it.
    """
    shared = sorted(set(a) & set(b))
    if not shared:
        return {
            "n_compared": 0,
            "n_identical": 0,
            "identity_fraction": None,
            "numbering_agrees": None,
            "mismatch_examples": [],
            "reason": (
                "the two entries share no residue number at all, so their "
                "numbering conventions could not be compared"
            ),
        }
    same = [r for r in shared if a[r] == b[r]]
    frac = round(len(same) / len(shared), 3)
    return {
        "n_compared": len(shared),
        "n_identical": len(same),
        "identity_fraction": frac,
        "numbering_agrees": frac >= NUMBERING_IDENTITY_MIN,
        "mismatch_examples": [
            f"{r}: {a[r]} vs {b[r]}" for r in shared if a[r] != b[r]
        ][:max_examples],
        "identity_threshold": NUMBERING_IDENTITY_MIN,
        "_threshold_status": "PROPOSED, NOT CALIBRATED",
        "reason": (
            None
            if frac >= NUMBERING_IDENTITY_MIN
            else (
                f"only {len(same)} of {len(shared)} shared residue numbers name "
                "the same residue in both entries. The two are not on a common "
                "numbering, so every seqid-keyed comparison between them "
                "(overlap_fraction, interface_coverage, shared_residues) is "
                "matching non-homologous positions. The geometric fields "
                "(min_distance_to_interface_a, enclosure) are unaffected."
            )
        ),
    }


def _classify_site_pocket(
    prepped: Path, used_chains: list[str], vert: Path, epitope,
    partner_resnames: dict[int, str] | None = None,
    apolar_fraction: float | None = None,
) -> dict:
    """One pocket, against the partner epitope, in the pocket's own structure.

    The prepared PDB is used rather than a fresh read of the CIF, because it is
    the exact file fpocket saw: same atoms, same chain IDs, same frame. The
    alpha-sphere centres in `pocket<N>_vert.pqr` are the probe points, which is
    what the module asks for — without them enclosure falls back to the lining
    centroid and overstates burial on a shallow groove.

    `match_by="seqid"` because the epitope comes from a DIFFERENT PDB entry.
    The price is that on a homo-oligomer the interface set becomes the union
    over protomers, which inflates overlap; an overlap of 0.00 under that
    inflation, as TNF-alpha's SPD304 site gives against the TNFR2 epitope, is
    therefore a strong statement rather than a weak one.

    THE OTHER PRICE, WHICH USED TO BE PAID SILENTLY: matching on seqid assumes
    the two entries number the protein the same way, and PDB entries for one
    protein routinely do not. `partner_resnames` is the partner entry's
    {seqid -> residue name} over its target chains, and it turns that assumption
    into a reported measurement — see `_numbering_agreement`. Without it, two of
    five TL1A structures produced overlap fractions computed against
    non-homologous residues with nothing in the payload to say so.
    """
    import numpy as np

    import interface_analysis as IA

    st = IA.load_structure(str(prepped))
    pts = np.asarray(_pdb_coords(vert, ("ATOM", "HETATM")), dtype=float)
    if len(pts) == 0:
        return {"error": f"no alpha-sphere centres in {vert.name}"}
    # Enclosure casts 512 rays per probe point per chain; a 400-sphere pocket
    # would dominate the runtime for no gain in the answer.
    if len(pts) > 60:
        pts = pts[np.linspace(0, len(pts) - 1, 60).astype(int)]
    lining = IA.residues_within(st, pts, chains=used_chains or None)
    if not lining:
        return {"error": "no residue within 4.5 A of the alpha-sphere centres"}

    # ---- PUT THE EPITOPE ON OUR NUMBERING *BEFORE* THE SEQID MATCH ---------
    #
    # THE CHECK FLAGGED AND DID NOT FIX, AND THE FIX WAS SITTING IN THE PAYLOAD.
    # `_numbering_agreement` correctly fired on exactly the two corrupt TL1A
    # structures (2O0O at 2/69 = 0.029, 2RE9 at 7/139 = 0.050, against
    # 0.993-1.000 for the three valid ones) — and the corrupt comparisons then
    # propagated upward unmarked, so `per_structure_consensus["2RE9"]` came back
    # `allosteric_candidate` derived ENTIRELY from a 0.227 overlap on
    # A:THR34/PRO35/THR36, residues that are VAL/VAL/ARG in the partner, and
    # `per_structure_consensus["2O0O"]` came back `mixed` off the
    # A:HIS118-vs-THR118 artifact. SKILL.md tells callers to quote exactly that
    # field.
    #
    # The offsets that fix it outright (+67, +71) are recoverable by a single
    # vote over residue names, which the mdpocket stage of the SAME payload
    # already runs. Recovering them here too is the better fix than flagging:
    # `interface_analysis`'s own module docstring says
    # `detect_numbering_offset` "should be run before any cross-entry
    # match_by='seqid' comparison", and this call site never ran one.
    #
    # The offset is only APPLIED when it turns an illegal comparison into a
    # legal one on a non-trivial overlap, so an entry that already agrees is
    # never gratuitously shifted and a spurious shift on a handful of positions
    # cannot buy agreement. Whatever survives unfixed is still flagged, and is
    # now also excluded from the consensus rather than merely annotated.
    ours: dict[int, str] = {}
    epitope_used = epitope
    offset_applied = 0
    numbering_note: str | None = None
    agree: dict | None = None
    raw_agree: dict | None = None
    off = 0
    if partner_resnames:
        ours = _pdb_resnames_by_seqid(prepped, used_chains or None)
        raw_agree = _numbering_agreement(ours, partner_resnames)
        # `ours[r]` is taken to be the same residue as `partner[r + off]`.
        off, _matched, _overlap = _best_numbering_offset(ours, partner_resnames)
        corrected = (
            _numbering_agreement(
                {r + off: nm for r, nm in ours.items()}, partner_resnames
            )
            if off else raw_agree
        )
        apply_it = bool(
            off
            and corrected.get("numbering_agrees")
            and (corrected.get("n_compared") or 0)
            >= NUMBERING_MIN_COMPARED_FOR_OFFSET
            and (corrected.get("n_identical") or 0)
            > (raw_agree.get("n_identical") or 0)
        )
        if apply_it:
            # The epitope is in the PARTNER's numbering; shift it into ours.
            epitope_used = IA.renumber_residues(epitope, -off)
            offset_applied = -off
            agree = corrected
            numbering_note = (
                f"NUMBERING OFFSET RECOVERED AND APPLIED: this entry's residue "
                f"r is the partner entry's r{off:+d}, so the partner epitope "
                f"was renumbered by {offset_applied:+d} before the seqid match. "
                f"Raw agreement was "
                f"{raw_agree.get('n_identical')}/{raw_agree.get('n_compared')} "
                f"= {raw_agree.get('identity_fraction')}; after the shift it is "
                f"{corrected.get('n_identical')}/{corrected.get('n_compared')} "
                f"= {corrected.get('identity_fraction')}. Every seqid-keyed "
                "field below (overlap_fraction, interface_coverage, "
                "shared_residues) is computed on the corrected numbering."
            )
        else:
            agree = raw_agree

    res = IA.classify_pocket(
        lining, epitope_used, st,
        target_chains=used_chains or None,
        probe_points=pts,
        match_by="seqid",
    )
    d = res.as_dict()
    d["summary"] = res.summary()

    # ---- is the seqid match even legal between these two entries? ----------
    if partner_resnames and agree is not None:
        d["numbering_check"] = agree
        d["numbering_offset_to_partner"] = {
            "offset_ours_to_partner": off,
            "epitope_renumbered_by": offset_applied,
            "applied": bool(offset_applied),
            "before_offset": raw_agree,
            "min_compared_to_apply": NUMBERING_MIN_COMPARED_FOR_OFFSET,
            "note": numbering_note,
            # `_why` deliberately NOT here. This dict is emitted once per
            # classified pocket per clustering value per structure, and a ~450
            # character explanation repeated that many times is how a payload
            # reaches a size cap and then truncates its own trailing
            # explanation. It is stated once, at
            # `pocket_vs_interface._numbering_offset_rule`.
        }
        if numbering_note:
            d["notes"] = list(d.get("notes") or []) + [numbering_note]
        if agree.get("numbering_agrees") is False:
            d["overlap_unreliable_numbering_mismatch"] = True
            d["notes"] = list(d.get("notes") or []) + [
                "OVERLAP IS NOT INTERPRETABLE: " + (agree.get("reason") or "")
            ]
        elif agree.get("numbering_agrees") is None:
            d["overlap_unreliable_numbering_mismatch"] = None
        else:
            d["overlap_unreliable_numbering_mismatch"] = False
        # Per-shared-residue names, so a reader can see the actual pairing
        # rather than trusting a fraction. `shared_residues` labels come from
        # the pocket side only, which is exactly how "shared A:HIS118" hid the
        # fact that the partner's 118 is a THR.
        # The shift is carried into this pairing too, or the check would report
        # the pre-correction mismatch on a comparison that was corrected.
        pshift = off if offset_applied else 0
        pairs = []
        for lbl in d.get("shared_residues") or []:
            digits = "".join(ch for ch in lbl if ch.isdigit())
            if not digits:
                continue
            n = int(digits)
            pn = n + pshift
            pairs.append({
                "seqid": n,
                "partner_seqid": pn,
                "pocket_residue": lbl,
                "partner_residue_name": partner_resnames.get(pn),
                "name_agrees": (
                    None if pn not in partner_resnames or n not in ours
                    else ours[n] == partner_resnames[pn]
                ),
            })
        d["shared_residue_name_check"] = pairs
        d["n_shared_residues_name_mismatched"] = sum(
            1 for p in pairs if p["name_agrees"] is False
        )

    # ---- buried core, keyed on geometry, never a filter --------------------
    d.update(
        _buried_core_flag(
            d.get("enclosure"), d.get("subunit_enclosure_gain"), apolar_fraction
        )
    )
    return d


# ===========================================================================
# mdpocket — the site fixed BY CONSTRUCTION
# ===========================================================================
#
# Everything above matches pockets ACROSS structures after the fact, by shared
# residue numbers. That is a heuristic and it has been measured failing: on the
# five apo TNF-alpha structures it tracked a pocket 7.7 A from the site it
# claimed, with 12.2 A of internal inconsistency between structures, and
# reported a large druggability spread that was never a measurement of one
# site. mdpocket removes the matching step entirely: ONE set of grid points is
# defined once, in a common frame, and the same points are characterised in
# every structure. Same five structures, volume CV ~28% -> ~10% (measured 28.1%
# at D=1.6 against 9.9%; ~1 percentage point of each is fpocket's Monte-Carlo
# volume noise, so do not quote a third significant figure).
#
# Four things about mdpocket that are not in its documentation and each of
# which silently corrupts the answer:
#
#   1. IT SILENTLY DROPS MISSING FRAMES AND EXITS 0. Verified directly: a
#      5-entry --pdb_list with one unreadable path prints "Identified 5
#      snapshots to analyze", then "! The pdb file ... doesn't exists.", then
#      processes "1/4 ... 4/4", writes a 4-line time.txt, and returns 0. The
#      frequency grid is normalised over however many frames actually ran, so a
#      dropped structure inflates EVERY frequency in the grid — 3/4 = 0.75
#      where the truth is 3/5 = 0.60. Nothing in the return code says so.
#      Hence `_assert_frame_count`, which is not optional.
#   2. IT DOES NOT SUPERPOSE. It treats the list as MD frames already in a
#      common frame. Deposited PDB entries are not. Everything is superposed
#      here first, Kabsch on core C-alpha.
#   3. ON A HOMO-OLIGOMER THE CHAIN MAPPING MUST BE SEARCHED. For a C3 trimer
#      the three cyclic chain mappings agree to within 0.03 A and the three
#      anticyclic ones land ~22 A out. Taking chains in file order picks an
#      anticyclic mapping roughly half the time and produces a superposition
#      that is geometrically absurd but numerically silent.
#   4. FREQUENCY IS QUANTISED AT 1/N. With N=5 the only attainable values are
#      {0, 0.2, 0.4, 0.6, 0.8, 1.0}. Quoting "this pocket is open 60% of the
#      time" off five structures is quoting the grid resolution. Frequencies
#      are therefore REFUSED below N=10 and the refusal is reported.
#
# time.txt is written in density mode only; characterisation mode does not
# write one, so there the equivalent invariant is the descriptor ROW count.

# mdpocket is on PATH in the Modal image. Locally it needs its conda env
# activated (the bare binary segfaults with rc 133 on a bare exec), hence the
# override — `MDPOCKET_CMD="micromamba run -n druggability mdpocket"`.
MDPOCKET_CMD = os.environ.get("MDPOCKET_CMD", "mdpocket").split()

# Below this many structures a frequency VALUE is grid resolution, not a
# measurement. See note 4 above.
MDPOCKET_MIN_N_FOR_FREQUENCY = 10

# Grid points within this distance of the transferred ligand define the
# ligand-anchored site. 3.0 A reproduces the calibration numbers; 2.0 A gives
# 148.4 A^3 where 3.0 A gives 154.1 A^3 on 1TNF, i.e. the answer is not
# sensitive to it.
MDPOCKET_SITE_CUTOFF_A = 3.0

# An anticyclic chain mapping on a C3 trimer lands ~22 A out against ~0.03 A
# for a cyclic one. Anything above this is not a superposition.
MDPOCKET_MAX_ACCEPTABLE_RMSD_A = 5.0

# A structure that will not superpose is DROPPED, not fatal — but a spread over
# the survivors of a partial ensemble is only a measurement if enough of them
# survived. Below this many, refuse rather than report a CV over two frames.
# An ensemble that was always going to be small (the KRAS 6OIM/4OBE pair) is
# unaffected: this floor applies only once a drop has happened.
MDPOCKET_MIN_N_AFTER_DROPS = 3


def _filter_chains(by_chain: dict, keep: list[str] | None) -> dict:
    """Restrict a {chain: ...} map to `keep`, or pass it through when unknown.

    Passing through on None is deliberate: a run with no accession has no
    verified target chains, and silently dropping chains on a guess would be the
    same class of error as scoring them on a guess.
    """
    if not keep:
        return by_chain
    kept = {c: v for c, v in by_chain.items() if c in keep}
    # Never empty the structure out. If the target chain names do not appear in
    # this prepared file at all, something upstream disagreed about naming and
    # the honest response is to fall back audibly rather than return nothing.
    return kept or by_chain


def _ca_by_chain(pdb: Path) -> dict[str, dict[int, list[float]]]:
    """{chain: {resseq: CA xyz}} straight out of a PDB, by fixed column."""
    out: dict[str, dict[int, list[float]]] = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
            continue
        if line[12:16].strip() != "CA":
            continue
        if line[16] not in (" ", "A"):
            continue
        try:
            out.setdefault(line[21], {})[int(line[22:26])] = [
                float(line[30:38]), float(line[38:46]), float(line[46:54])
            ]
        except ValueError:
            continue
    return out


def _res_names_by_chain(pdb: Path) -> dict[str, dict[int, str]]:
    """{chain: {resseq: residue name}} from the CA records of a PDB."""
    out: dict[str, dict[int, str]] = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
            continue
        if line[12:16].strip() != "CA" or line[16] not in (" ", "A"):
            continue
        try:
            out.setdefault(line[21], {})[int(line[22:26])] = line[17:20].strip()
        except ValueError:
            continue
    return out


def _pdb_resnames_by_seqid(
    pdb: Path, chains: list[str] | None = None
) -> dict[int, str]:
    """{resseq: residue name}, pooled over chains. For numbering comparisons.

    Pooling is correct here because the consumer is a chain-agnostic seqid
    match; on a homo-oligomer the protomers carry identical numbering by
    construction, so pooling loses nothing. A genuine disagreement between two
    chains of one entry would show up as an arbitrary pick, which is why the
    FIRST name seen for a number wins and later ones do not overwrite it.
    """
    out: dict[int, str] = {}
    for ch, m in _res_names_by_chain(pdb).items():
        if chains and ch not in chains:
            continue
        for num, name in m.items():
            out.setdefault(num, name)
    return out


# A numbering offset search wider than this is not a construct offset, it is a
# different protein. TL1A's measured offsets are +67 and +71.
MAX_NUMBERING_OFFSET = 2000


def _best_numbering_offset(
    mob: dict[int, str], ref: dict[int, str]
) -> tuple[int, int, int]:
    """Integer offset that best maps `mob`'s numbering onto `ref`'s.

    Returns (offset, n_matched, n_overlapping). `mob[r]` is taken to be the same
    residue as `ref[r + offset]`.

    WHY THIS EXISTS, AND IT IS THE MOST SERIOUS BUG THIS FILE HAS CARRIED.
    The superposition indexed C-alpha by RAW AUTHOR RESIDUE NUMBER. PDB entries
    for one protein routinely use different numbering conventions, and TL1A's
    ensemble carries three at once: 2O0O at offset 0, then 2QE3/2RJK/2RJL/3K51/
    3MI8 at +67, and 2RE9 at +71. Fitting on raw numbers therefore fitted
    non-homologous residues onto each other:

        numbering             core CA   best 3-chain RMSD, 2O0O vs the rest
        raw author (old)           67   18.70 - 20.06 A
        aligned (this function)   138   0.51 - 1.45 A, clean C3 split

    The ensemble superposes essentially perfectly. The tool reported
    "2QE3: best chain mapping RMSD 14.84 A exceeds 5.0 A; not a superposition",
    which reads as a conformational problem and is not one — the error message
    misdiagnosed its own failure and the whole mdpocket stage was lost.

    The old `len(core) < 20` guard did not catch it and could not: its message
    already said "the entries do not share a numbering", but it tests a COUNT.
    67 residues aligned by accident at a constant offset clear a count of 20
    comfortably. What catches it is asserting residue IDENTITY at the matched
    positions, which is what `_common_core` below does.

    Offsets, not a gapped alignment, deliberately: a constant offset is the
    failure that occurs in deposited entries, it is exactly recoverable, and the
    identity assertion downstream shrinks the core wherever the offset model is
    wrong rather than letting a bad alignment through. An entry needing a truly
    gapped alignment simply yields a smaller core and, if it gets too small, an
    explicit refusal.
    """
    from collections import Counter

    if not mob or not ref:
        return 0, 0, 0
    # Every offset that could match ANY residue by name, voted once per pair.
    # votes[off] is exactly the number of name-matched positions at that
    # offset, so this is the same answer as scanning a window and costs
    # O(n_mob * n_ref / 20) instead of O(window * n_mob).
    by_name: dict[str, list[int]] = {}
    for num, name in ref.items():
        by_name.setdefault(name, []).append(num)
    votes: Counter = Counter()
    for num, name in mob.items():
        for rnum in by_name.get(name, ()):
            off = rnum - num
            if -MAX_NUMBERING_OFFSET <= off <= MAX_NUMBERING_OFFSET:
                votes[off] += 1
    if not votes:
        return 0, 0, len(set(mob) & set(ref))
    # Most matches wins; ties go to the smallest shift, so an entry that already
    # shares the reference numbering is never gratuitously renumbered.
    off = max(votes, key=lambda o: (votes[o], -abs(o)))
    overlap = sum(1 for num in mob if (num + off) in ref)
    return off, votes[off], overlap


def _align_numbering(
    ca: dict[str, dict[int, list[float]]],
    names: dict[str, dict[int, str]],
    ref_names: dict[int, str],
) -> tuple[dict[str, dict[int, list[float]]], dict[str, dict[int, str]], dict]:
    """Renumber every chain of one structure onto the reference's numbering.

    Only the internal fit uses the renumbered dictionaries; the PDB written for
    mdpocket keeps its own deposited numbering, because mdpocket reads
    coordinates and never residue numbers.
    """
    new_ca: dict[str, dict[int, list[float]]] = {}
    new_names: dict[str, dict[int, str]] = {}
    report: dict[str, dict] = {}
    for ch, m in names.items():
        off, matched, overlap = _best_numbering_offset(m, ref_names)
        new_ca[ch] = {r + off: v for r, v in ca.get(ch, {}).items()}
        new_names[ch] = {r + off: v for r, v in m.items()}
        report[ch] = {
            "offset_applied": off,
            "n_positions_matched_by_name": matched,
            "n_positions_overlapping": overlap,
            "identity_fraction": (
                round(matched / overlap, 3) if overlap else None
            ),
        }
    offsets = sorted({v["offset_applied"] for v in report.values()})
    return new_ca, new_names, {
        "per_chain": report,
        "offsets_applied": offsets,
        "renumbered": any(o != 0 for o in offsets),
    }


def _common_core(
    names_by_pid: dict[str, dict[str, dict[int, str]]],
) -> tuple[list[int], dict]:
    """Residue numbers present in EVERY chain of EVERY structure *and naming
    the same residue in all of them*.

    THE IDENTITY HALF IS THE GUARD THAT WAS MISSING. A count of shared numbers
    cannot tell a homologous core from an accidental overlap at a constant
    offset — 67 positions of TL1A's ensemble passed `len(core) >= 20` while
    being non-homologous throughout. Requiring the residue NAME to agree at
    every core position is what makes the core a core.
    """
    shared: set[int] | None = None
    for by_chain in names_by_pid.values():
        for m in by_chain.values():
            shared = set(m) if shared is None else shared & set(m)
    shared = shared or set()
    core, dropped = [], []
    for num in sorted(shared):
        seen = {
            m[num]
            for by_chain in names_by_pid.values()
            for m in by_chain.values()
        }
        (core if len(seen) == 1 else dropped).append(num)
    return core, {
        "n_shared_numbers": len(shared),
        "n_core_positions": len(core),
        "n_dropped_on_residue_name_disagreement": len(dropped),
        "dropped_positions_sample": dropped[:20],
        "_why": (
            "A core position must be the SAME residue in every structure, not "
            "merely a number they all happen to carry. Counting numbers alone "
            "let a 67-position non-homologous overlap at a constant numbering "
            "offset pass as a valid core on TL1A and turned a 0.5 A "
            "superposition into a reported 14.8-20.1 A refusal."
        ),
    }


def _kabsch(P, Q):
    """R, t such that P @ R + t is the least-squares fit of P onto Q.

    numpy only — scipy is not in this image and is not needed for this.
    """
    import numpy as np

    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = (Vt.T @ np.diag([1.0, 1.0, d]) @ U.T).T
    return R, qc - pc @ R


def _fit_to_reference(
    mob_ca: dict[str, dict[int, list[float]]],
    ref_ca: dict[str, dict[int, list[float]]],
    core: list[int],
    mob_chains: list[str] | None = None,
    max_perms: int = 2000,
) -> dict:
    """Best chain mapping of one structure onto the reference, by RMSD.

    THE PERMUTATION SEARCH IS THE POINT. Chain letters are not a biological
    correspondence — 1TNF's A and 2ZJC's A need not be the same protomer, and on
    a C3 trimer picking wrongly costs ~22 A rather than ~0.03 A. Every mapping
    is scored and all of them are returned, so the caller can see the cyclic /
    anticyclic split rather than trusting a single number.

    BOTH SIDES are searched — ordered subsets of the mobile chains against
    unordered subsets of the reference chains. Permuting only the mobile side
    is not enough and fails on a real case: 2AZ5's asymmetric unit is two
    TNF-alpha DIMERS (4 chains) and the reference is a TRIMER (3), so every
    3-chain mapping is geometrically impossible and the best of them is 17.3 A.
    The correct fit maps the ligand's own dimer onto TWO of the trimer's three
    protomers, which requires choosing a subset of the REFERENCE.

    `mob_chains` restricts the mobile side, and for a site donor it must be
    given: it is what keeps a crystallographic second copy out of the fit.
    """
    import itertools

    import numpy as np

    ref_chains = sorted(ref_ca)
    mob_list = sorted(mob_chains) if mob_chains else sorted(mob_ca)
    mob_list = [c for c in mob_list if c in mob_ca]
    if not ref_chains or not mob_list or not core:
        return {"ok": False, "reason": "no reference chains or no core residues"}

    # Map as many chains as the smaller side has. A subset mapping is correct
    # behaviour, not a fallback: 6OIM's assembly is a monomer and 4OBE's need
    # not be.
    k = min(len(ref_chains), len(mob_list))
    if k == 0:
        return {"ok": False, "reason": "no chains to map"}

    # LARGEST MAPPING THAT IS ACTUALLY A SUPERPOSITION, then fall back. Two
    # assemblies need not be the same oligomer: 2AZ5's asymmetric unit is two
    # TNF-alpha DIMERS and 1TNF is a TRIMER, so no 3-chain mapping between them
    # exists at all — the best is 17.3 A — while the 2-chain mapping of the
    # ligand's own dimer onto two protomers of the trimer is 1.13 A. Refusing
    # at k=3 threw away a measurement that was available at k=2.
    #
    # Reducing k always lowers RMSD (fewer constraints), so the search goes
    # DOWNWARD from the largest k and stops at the first acceptable one. It
    # never trades chains for a prettier number.
    scored: list[tuple[float, tuple[str, ...], tuple[str, ...]]] = []
    truncated = False
    rms = perm = tgt = None
    for kk in range(k, 0, -1):
        pairs = [
            (p, t) for p in itertools.permutations(mob_list, kk)
            for t in itertools.combinations(ref_chains, kk)
        ]
        truncated = truncated or len(pairs) > max_perms
        level: list[tuple[float, tuple[str, ...], tuple[str, ...]]] = []
        for p, t in pairs[:max_perms]:
            try:
                Pp = np.array([mob_ca[c][r] for c in p for r in core])
                Qq = np.array([ref_ca[c][r] for c in t for r in core])
            except KeyError:
                continue
            R, t_ = _kabsch(Pp, Qq)
            level.append((
                float(np.sqrt((((Pp @ R + t_) - Qq) ** 2).sum(1).mean())), p, t
            ))
        if not level:
            continue
        level.sort(key=lambda x: x[0])
        if kk == k:
            scored = level
        if level[0][0] <= MDPOCKET_MAX_ACCEPTABLE_RMSD_A or kk == 1:
            if kk != k:
                scored = level
            rms, perm, tgt = level[0]
            k = kk
            break
    if rms is None:
        return {"ok": False, "reason": "no chain permutation shared the core"}
    Pp = np.array([mob_ca[c][r] for c in perm for r in core])
    Qq = np.array([ref_ca[c][r] for c in tgt for r in core])
    R, t = _kabsch(Pp, Qq)
    return {
        "ok": rms <= MDPOCKET_MAX_ACCEPTABLE_RMSD_A,
        "reason": (
            None
            if rms <= MDPOCKET_MAX_ACCEPTABLE_RMSD_A
            else f"best chain mapping RMSD {rms:.2f} A exceeds "
                 f"{MDPOCKET_MAX_ACCEPTABLE_RMSD_A} A; not a superposition"
        ),
        "rmsd_a": round(rms, 3),
        "chain_map": {m: r for m, r in zip(perm, tgt, strict=False)},
        "n_chains_mapped": k,
        "n_mobile_chains": len(mob_list),
        "n_reference_chains": len(ref_chains),
        "unmapped_mobile_chains": [c for c in mob_list if c not in perm],
        "n_core_ca": len(Pp),
        # Every mapping, so the C3 cyclic/anticyclic split is visible rather
        # than assumed: the three good ones agree to ~0.03 A, the three bad
        # ones sit ~22 A away.
        "all_mapping_rmsd_a": sorted(round(s[0], 3) for s in scored)[:24],
        "n_mappings_tried": len(scored),
        "permutations_truncated": truncated,
        "_R": R,
        "_t": t,
    }


# The reference search is O(n^2) fits. Above this many structures it stops
# being worth its cost and the first entry is used with the fallback recorded.
MDPOCKET_MAX_IDS_FOR_REFERENCE_SEARCH = 12


def _select_mdpocket_reference(
    ids: list[str],
    cas_raw: dict[str, dict],
    names_raw: dict[str, dict],
) -> dict:
    """Which structure the ensemble is superposed ONTO. It was `ids[0]`.

    THE REFERENCE IS A MEASUREMENT DECISION AND IT WAS BEING MADE BY DICT ORDER.
    Every other structure is fitted onto the reference and DROPPED if it will
    not go, so picking an outlier as the reference does not fail the outlier —
    it fails everything else, and the payload then reports the ensemble as
    unusable rather than reporting the outlier.

    Measured on NLRP3 (8SWF, 7ZGU, 9HG4). 8SWF came first, so 7ZGU and 9HG4 were
    dropped at 16.43 and 16.55 A, 1 of 3 survived, and the whole mdpocket stage
    was refused — while 7ZGU and 9HG4 superpose onto EACH OTHER at 2.69 A by
    this stage's own whole-assembly fit, and at 1.301 A in the cryptic stage of
    the same payload, which had already measured and printed it. One structure
    cost an entire ensemble by being listed first.

    Verified after the change, on the prepared NLRP3 trio:

        candidate   median RMSD to the rest   would keep
        8SWF        16.49                     0 of 2      <- was the reference
        7ZGU         9.559                    1 of 2      <- is now
        9HG4         9.621                    1 of 2

    Lowest MEDIAN RMSD to the rest. Not the first, not the largest, and the
    median specifically: an outlier's own median is large by construction so it
    cannot elect itself, and one bad entry cannot unseat a good reference the
    way a mean would let it. `n_would_superpose` is reported beside it and
    breaks ties, because the quantity actually being protected is how many
    structures survive.

    Note what this does NOT fix on NLRP3: the best reference is 7ZGU, 8SWF is
    then correctly dropped as the outlier, and 2 of 3 survive — still below
    MDPOCKET_MIN_N_AFTER_DROPS, so that run still refuses. It refuses for the
    right reason and names the right outlier, which is the whole difference.
    """
    import statistics

    report: dict[str, dict] = {}
    if len(ids) > MDPOCKET_MAX_IDS_FOR_REFERENCE_SEARCH:
        return {
            "reference": ids[0],
            "selected_by": (
                f"first entry — {len(ids)} structures exceeds the "
                f"{MDPOCKET_MAX_IDS_FOR_REFERENCE_SEARCH}-structure cap on the "
                "O(n^2) median-RMSD search"
            ),
            "candidates": report,
        }
    best: tuple[tuple, str] | None = None
    for cand in ids:
        if not cas_raw.get(cand):
            report[cand] = {"error": "no C-alpha", "median_rmsd_a": None}
            continue
        # Same pooling rule as `_pdb_resnames_by_seqid`, off the names we have
        # already read rather than off the file again.
        ref_names: dict[int, str] = {}
        for m in names_raw[cand].values():
            for num, nm in m.items():
                ref_names.setdefault(num, nm)
        cas_c: dict[str, dict] = {}
        names_c: dict[str, dict] = {}
        for pid in ids:
            cas_c[pid], names_c[pid], _rep = _align_numbering(
                cas_raw[pid], names_raw[pid], ref_names
            )
        core_c, _cr = _common_core(names_c)
        if len(core_c) < 20:
            report[cand] = {
                "n_core_ca": len(core_c), "median_rmsd_a": None,
                "error": (
                    f"only {len(core_c)} common core positions against this "
                    "candidate (20 needed)"
                ),
            }
            continue
        rmsds: dict[str, float | None] = {}
        for pid in ids:
            if pid == cand:
                continue
            fit = _fit_to_reference(cas_c[pid], cas_c[cand], core_c)
            rmsds[pid] = fit.get("rmsd_a")
        vals = [v for v in rmsds.values() if v is not None]
        med = round(statistics.median(vals), 3) if vals else None
        n_ok = sum(1 for v in vals if v <= MDPOCKET_MAX_ACCEPTABLE_RMSD_A)
        report[cand] = {
            "n_core_ca": len(core_c),
            "rmsd_to_others_a": rmsds,
            "median_rmsd_a": med,
            "n_would_superpose": n_ok,
            "n_others": len(rmsds),
        }
        if med is None:
            continue
        key = (med, -n_ok, cand)
        if best is None or key < best[0]:
            best = (key, cand)
    if best is None:
        return {
            "reference": ids[0],
            "selected_by": (
                "first entry — no candidate produced a usable common core, so "
                "there was nothing to rank"
            ),
            "candidates": report,
        }
    return {
        "reference": best[1],
        "selected_by": (
            "lowest MEDIAN core C-alpha RMSD to the other structures "
            "(ties to the one that keeps the most structures). The reference "
            "was previously the first entry, and on NLRP3 that was the single "
            "outlier: it dropped the two structures that superpose onto each "
            "other at 1.301 A and cost the whole ensemble."
        ),
        "candidates": report,
    }


def _write_superposed(src: Path, dst: Path, R, t, chain_map: dict[str, str]) -> int:
    """Rewrite a prepared PDB into the reference frame, relabelling chains.

    EVERY CHAIN IS KEPT, not only the mapped ones. When the chain mapping falls
    back to a smaller k — 2AZ5's two dimers against 1TNF's trimer map 2-on-2,
    not 3-on-3 — dropping the unmapped chains would delete the third protomer
    of the trimer, and TNF-alpha's site IS the trimer. The pocket would simply
    not be there and mdpocket would return a confident wrong volume.

    Unmapped chains keep their own label unless it collides with a mapped
    target label, in which case they take the next free one.
    """
    import numpy as np

    taken = set(chain_map.values())
    full = dict(chain_map)
    pool = [c for c in _PDB_CHAIN_POOL if c not in taken]
    for line in src.read_text().splitlines():
        if not line.startswith("ATOM") or len(line) < 22:
            continue
        ch = line[21]
        if ch in full:
            continue
        if ch not in taken:
            full[ch] = ch
            taken.add(ch)
        elif pool:
            full[ch] = pool.pop(0)
            taken.add(full[ch])
        else:
            full[ch] = ch  # out of letters; collision is visible, not silent

    lines = []
    for line in src.read_text().splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
            continue
        xyz = np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        )
        q = xyz @ R + t
        lines.append(
            line[:21] + full.get(line[21], line[21]) + line[22:30]
            + f"{q[0]:8.3f}{q[1]:8.3f}{q[2]:8.3f}" + line[54:]
        )
    dst.write_text("\n".join(lines) + "\nEND\n")
    return len(lines)


def _run_mdpocket(args: list[str], cwd: Path, timeout: int) -> subprocess.CompletedProcess:
    return subprocess.run(  # noqa: S603
        MDPOCKET_CMD + args,
        cwd=str(cwd), check=False, capture_output=True, timeout=timeout,
    )


def _assert_frame_count(observed: int, expected: int, what: str) -> None:
    """THE NON-NEGOTIABLE CHECK. See note 1 in the block comment above."""
    if observed != expected:
        raise RuntimeError(
            f"mdpocket processed {observed} of {expected} structures ({what}). "
            "It drops unreadable frames silently and still exits 0, and the "
            "frequency grid is then normalised over the frames that ran — every "
            "frequency in the grid would be inflated. Refusing the result."
        )


def _read_dx(path: Path):
    """OpenDX scalar grid -> (values 1-D in C order, counts, origin, spacing)."""
    import numpy as np

    counts = origin = None
    deltas: list[list[float]] = []
    vals: list[float] = []
    indata = False
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("object 1"):
                counts = tuple(int(x) for x in s.split()[-3:])
                continue
            if s.startswith("origin"):
                origin = np.array([float(x) for x in s.split()[1:4]])
                continue
            if s.startswith("delta"):
                deltas.append([float(x) for x in s.split()[1:4]])
                continue
            if s.startswith("object 3"):
                indata = True
                continue
            if s.startswith(("object", "attribute", "component")):
                continue
            if indata:
                vals.extend(float(x) for x in s.split())
    if counts is None or origin is None or len(deltas) != 3:
        raise RuntimeError(f"{path.name}: not a readable OpenDX grid")
    n = int(np.prod(counts))
    return np.asarray(vals[:n], dtype=float), counts, origin, np.diag(np.array(deltas))


def _dx_points(counts, origin, spacing):
    import numpy as np

    ix, iy, iz = np.meshgrid(
        np.arange(counts[0]), np.arange(counts[1]), np.arange(counts[2]),
        indexing="ij",
    )
    return np.stack(
        [
            origin[0] + ix.ravel() * spacing[0],
            origin[1] + iy.ravel() * spacing[1],
            origin[2] + iz.ravel() * spacing[2],
        ],
        axis=1,
    )


def _largest_grid_cluster(points, spacing):
    """Largest 26-connected cluster of grid points, BFS on grid adjacency.

    A frequency mask is not a pocket — it is a mask, and on a trimer it
    contains a large disconnected surface film as well as the real internal
    cavity. Taking the biggest connected component is what turns it into one.
    """
    from collections import deque

    import numpy as np

    if len(points) == 0:
        return np.zeros(0, dtype=int)
    # Index on integer grid coordinates so adjacency is exact rather than a
    # float distance test.
    # np.diag returns a read-only VIEW, so this must be a copy before assignment.
    step = np.array(spacing, dtype=float, copy=True)
    step[step == 0] = 1.0
    key = {}
    ijk = np.rint((points - points.min(0)) / step).astype(int)
    for n, cell in enumerate(map(tuple, ijk)):
        key[cell] = n
    offs = [
        (i, j, k)
        for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)
        if (i, j, k) != (0, 0, 0)
    ]
    seen = np.zeros(len(points), dtype=bool)
    best: list[int] = []
    for s in range(len(points)):
        if seen[s]:
            continue
        q = deque([s])
        seen[s] = True
        cur = [s]
        while q:
            u = q.popleft()
            cu = ijk[u]
            for o in offs:
                v = key.get((cu[0] + o[0], cu[1] + o[1], cu[2] + o[2]))
                if v is not None and not seen[v]:
                    seen[v] = True
                    q.append(v)
                    cur.append(v)
        if len(cur) > len(best):
            best = cur
    return np.asarray(best, dtype=int)


def _write_probe_pdb(points, path: Path) -> int:
    with path.open("w") as fh:
        for i, p in enumerate(points):
            fh.write(
                f"ATOM  {i + 1:5d}  C   PTH     1    "
                f"{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}  0.00  0.00\n"
            )
        fh.write("END\n")
    return len(points)


def _read_descriptors(path: Path) -> list[dict]:
    """mdpocket's characterisation table: one row per snapshot, in list order."""
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    if len(lines) < 2:
        return []
    hdr = lines[0].split()
    rows = []
    for ln in lines[1:]:
        parts = ln.split()
        if len(parts) < len(hdr):
            continue
        row = {}
        for k, v in zip(hdr, parts, strict=False):
            try:
                row[k] = float(v)
            except ValueError:
                row[k] = v
        rows.append(row)
    return rows


def _spread(values: list[float]) -> dict:
    """Mean / SD / CV. CV is the number rule 4 of the dossier is about."""
    vals = [float(v) for v in values]
    if not vals:
        return {"n": 0, "mean": None, "sd": None, "cv_pct": None,
                "min": None, "max": None}
    n = len(vals)
    mean = sum(vals) / n
    sd = (
        (sum((v - mean) ** 2 for v in vals) / (n - 1)) ** 0.5 if n > 1 else 0.0
    )
    return {
        "n": n,
        "mean": round(mean, 2),
        "sd": round(sd, 2),
        "cv_pct": round(100 * sd / mean, 1) if mean else None,
        "min": round(min(vals), 2),
        "max": round(max(vals), 2),
        "values": [round(v, 2) for v in vals],
    }


def _mdpocket_characterise(
    run: Path, listfile: Path, sel: Path, prefix: str, ids: list[str], timeout: int
) -> dict:
    """One selected-pocket run: same grid points, every structure."""
    proc = _run_mdpocket(
        ["--pdb_list", listfile.name, "--selected_pocket", sel.name, "-o", prefix],
        run, timeout,
    )
    desc = run / f"{prefix}_descriptors.txt"
    if not desc.exists():
        raise RuntimeError(
            f"mdpocket wrote no {desc.name} (exit {proc.returncode}): "
            f"{proc.stdout.decode(errors='replace')[-400:]}"
        )
    rows = _read_descriptors(desc)
    # Characterisation mode writes no time.txt, so the descriptor row count is
    # the frame-count invariant here. Same hazard, same refusal.
    _assert_frame_count(len(rows), len(ids), f"{desc.name} rows")
    vols = [float(r.get("pock_volume", 0.0)) for r in rows]
    present = [pid for pid, v in zip(ids, vols, strict=False) if v > 0.0]
    out = {
        "n_probe_points": sum(1 for ln in sel.read_text().splitlines()
                              if ln.startswith("ATOM")),
        "volume_a3_by_structure": dict(zip(ids, [round(v, 2) for v in vols],
                                           strict=False)),
        "volume_a3": _spread(vols),
        "volume_a3_where_present": _spread([v for v in vols if v > 0.0]),
        "n_structures_with_volume": len(present),
        "structures_with_volume": present,
        "structures_without_volume": [p for p in ids if p not in present],
        # ---------------------------------------------------------------
        # THE DESCRIPTORS UNDER THEIR OWN NAMES, and NO druggability.
        #
        # This block used to carry `druggability_by_structure`, populated from
        # the `volume_score` column. Observed values were 3.35 to 4.00, which is
        # impossible for a [0,1] score and matches `volume_score` exactly. It is
        # the worst class of bug we get: a plausible number under a field name
        # that invites it to be quoted as something else, and it was.
        #
        # IT IS NOT A COLUMN-INDEX SLIP AND THERE IS NO RIGHT COLUMN.
        # mdpocket's characterisation table is fixed at 22 descriptors plus 20
        # amino-acid counts (`M_MDP_OUTP_HEADER` in fpocket's mdpocket.h) and
        # NONE of them is a druggability score:
        #
        #   snapshot pock_volume pock_asa pock_pol_asa pock_apol_asa pock_asa22
        #   pock_pol_asa22 pock_apol_asa22 nb_AS mean_as_ray mean_as_solv_acc
        #   apol_as_prop mean_loc_hyd_dens hydrophobicity_score volume_score
        #   polarity_score charge_score prop_polar_atm as_density as_max_dst
        #   convex_hull_volume nb_abpa  + ALA..VAL
        #
        # Nor can it be reconstructed. fpocket's shipped score (pscoring.c,
        # drug_score_pocket) is
        #     sigmoid(-9.5699 + 7.4798*mean_loc_hyd_dens_norm
        #             + 0.3696*as_max_dst - 0.04672*surf_pol_vdw22)
        # and `mean_loc_hyd_dens_norm` is MIN-MAX NORMALISED ACROSS THE OTHER
        # POCKETS OF THE SAME STRUCTURE (pocket.c, set_normalized_descriptors).
        # A druggability score is therefore not a property of a pocket at all;
        # it is a property of a pocket relative to the pocket population it was
        # detected with. mdpocket characterisation has a population of exactly
        # one fixed grid, so the normalisation has no referent. Applying
        # fpocket's single-pocket fallback constants ((mlhd-8.23)/15.97) to the
        # 6OIM switch-II row gives a saturated 1.000, which is not a measurement
        # either.
        #
        # So the honest answer is null with a reason. Druggability from the
        # fpocket path (`structures.<ID>.by_clustering.<D>`) is the one that
        # exists; it is a per-structure number and must not be presented as an
        # mdpocket fixed-site number.
        "druggability_by_structure": None,
        "druggability_status": "not_available",
        "druggability_reason": (
            "mdpocket characterisation mode emits no druggability score, and "
            "fpocket's score cannot be reconstructed from what it does emit: "
            "the score's leading term is min-max normalised across the other "
            "pockets of the same structure, and a fixed grid has no such "
            "population. This field previously reported the `volume_score` "
            "descriptor (observed 3.35-4.00) under this name. Read "
            "volume_score_by_structure for that descriptor, and take "
            "druggability from the fpocket path only."
        ),
        "volume_score_by_structure": dict(
            zip(ids, [round(float(r.get("volume_score", 0.0)), 3) for r in rows],
                strict=False)
        ),
        "_volume_score_note": (
            "fpocket's `volume_score` descriptor, NOT a druggability score and "
            "NOT on [0,1]. It is one of the raw scoring descriptors and it is "
            "the value this block used to mislabel as druggability."
        ),
        "mean_local_hydrophobic_density_by_structure": dict(
            zip(ids,
                [round(float(r.get("mean_loc_hyd_dens", 0.0)), 2) for r in rows],
                strict=False)
        ),
        "hydrophobicity_score_by_structure": dict(
            zip(ids,
                [round(float(r.get("hydrophobicity_score", 0.0)), 2) for r in rows],
                strict=False)
        ),
        "descriptors": rows,
    }
    # The range assertion. `volume_a3` and the descriptors are unbounded by
    # nature; the only [0,1] quantity this function can emit is druggability,
    # and it is now explicitly null. Asserting it anyway means that if anyone
    # ever repopulates it from a column, a 4.00 cannot leave the function.
    _unit_interval(
        out["druggability_by_structure"], "druggability_by_structure",
        f"{prefix} over {len(ids)} structures",
    )
    for pid_, v in (out["volume_a3_by_structure"] or {}).items():
        if v < 0.0:
            raise RuntimeError(
                f"{prefix}: negative pocket volume {v} A^3 for {pid_}; "
                "mdpocket does not produce one and this is a parse error"
            )
    return out


# A site centroid this far from the transferred ligand centroid is a DIFFERENT
# POCKET, not a noisy estimate of the same one. NOT A CALIBRATED THRESHOLD — it
# is a proposal, chosen as roughly half the 7.73 A miss measured on apo
# TNF-alpha, i.e. comfortably below the one error we know about and above the
# ~1 A grid spacing. It gates a WARNING, never a refusal, and no number in this
# module is dropped because of it. Calibrating it needs more than one case.
MDPOCKET_SITE_OFFSITE_WARN_A = 4.0


def _ligand_distance_fields(
    centroid, lig_common, donor_pid: str | None, donor_fit: dict | None
) -> dict:
    """`distance_to_donor_ligand_centroid_a`, ALWAYS, plus why when it is null.

    The distance field only prevents a misread if it is there to be read. It
    used to be emitted only when a ligand had been transferred, which is
    backwards: the case that most needs the warning is the pure-apo run where
    nothing was transferable and `site_from_density` comes back as the only
    site, with `mdpocket_status: "ok"`, and nothing in the payload says it may
    not be the pocket the dossier asked about. So a null here is itself the
    finding, and it carries its own reason.
    """
    import numpy as np

    if lig_common is not None:
        d = round(float(np.linalg.norm(np.asarray(centroid) - lig_common.mean(0))), 2)
        return {
            "distance_to_donor_ligand_centroid_a": d,
            "distance_reason": None,
            "donor_pdb_id": donor_pid,
            "off_site_warning": (
                None if d <= MDPOCKET_SITE_OFFSITE_WARN_A else
                f"this site's centroid is {d} A from the {donor_pid} ligand "
                f"centroid, past the {MDPOCKET_SITE_OFFSITE_WARN_A} A proposed "
                "(NOT calibrated) warning distance. At 7.73 A on apo TNF-alpha "
                "this is a different pocket — the on-axis cavity, the one the "
                "retracted matcher reported as the SPD304 site. Do not report "
                "this volume or druggability as the ligand site."
            ),
        }
    if not donor_pid:
        reason = (
            "no site donor: no holo structure in the ensemble carried a "
            "drug-like ligand and no mdpocket_site_donor was supplied, so "
            "there is no ligand to measure a distance to. THIS SITE IS "
            "THEREFORE UNVERIFIED — it is the most persistent cavity, and "
            "nothing here establishes that it is the site of interest."
        )
    else:
        reason = (
            f"site donor {donor_pid} was resolved but its ligand could not be "
            f"transferred into the common frame "
            f"(fit: {(donor_fit or {}).get('reason') or 'no fit attempted'}), "
            "so no distance is available. THIS SITE IS THEREFORE UNVERIFIED."
        )
    return {
        "distance_to_donor_ligand_centroid_a": None,
        "distance_reason": reason,
        "donor_pdb_id": donor_pid,
        "off_site_warning": (
            "cannot be checked — no ligand-anchored reference was available. "
            "Treat this site as a cavity of unknown identity, not as the "
            "dossier's site, and say so in tractability.site_hypothesis_basis."
        ),
    }


def _mdpocket_ensemble(
    work: Path,
    prepped: dict[str, Path],
    donor_pid: str | None,
    donor_prepped: Path | None,
    donor_ligand_xyz: list[list[float]] | None,
    donor_chains: list[str] | None = None,
    timeout: int = 900,
    target_chains: dict[str, list[str]] | None = None,
) -> dict:
    """Superpose the ensemble, then measure ONE site definition in all of it.

    Returns a dict always carrying `mdpocket_status` in {ok, failed, not_run}.
    Two site definitions are reported and they answer different questions:

      * `site_from_ligand` — grid points around the holo ligand, transferred by
        superposition. This is the site the dossier is asking about. On the
        five apo TNF-alpha structures it returns 0.00 A^3 in four of them, and
        THAT IS THE CORRECT ANSWER: the SPD304 channel is not open in those
        crystal forms with all three subunits present. A refusal, not a
        failure. Reading it as "mdpocket failed" is the mistake this comment
        exists to prevent.
      * `site_from_density` — the largest connected cluster of grid points at
        which a pocket is present in EVERY structure. Independent of any
        ligand, so it is the one that exists for a pure apo ensemble, and it is
        the measurement behind the CV ~28% -> ~10% improvement.

    THE SECOND ONE IS NOT THE LIGAND SITE AND ON OUR BEST-CHARACTERISED TEST
    CASE IT IS THE WRONG POCKET. On the apo TNF-alpha ensemble
    `site_from_density`'s centroid sits 7.73 A from the transferred SPD304
    ligand: it is the on-axis cavity, i.e. precisely the pocket the retracted
    residue-number matcher reported as "the SPD304 site". It is a real cavity
    and a real measurement — it is just the answer to a different question.

    So EVERY site entry carries `ligand_anchored` (bool) and
    `distance_to_donor_ligand_centroid_a`, unconditionally, with
    `distance_reason` when the distance could not be computed. They are emitted
    even when there is no donor at all, because the dangerous case is exactly
    the one where a pure-apo run with no transferable ligand returns
    `site_from_density` as the ONLY site, with `mdpocket_status: "ok"` — a
    confident single answer about the wrong pocket. A consumer that reads a
    volume off a site entry without reading these two fields reproduces the
    retracted bug, and it will look like a result rather than a bug.
    """
    import numpy as np

    ids = list(prepped)
    n = len(ids)
    out: dict = {
        "mdpocket_status": "not_run",
        "mdpocket_reason": None,
        "n_structures": n,
        "_why": (
            "The site is fixed BY CONSTRUCTION: one set of grid points defined "
            "once in a common frame and characterised in every structure. This "
            "replaces post-hoc residue-number matching, which on the five apo "
            "TNF-alpha structures tracked a pocket 7.7 A from the site it "
            "claimed and inflated the volume CV from ~10% to ~28%. "
            "FIXED BY CONSTRUCTION MEANS REPRODUCIBLE, NOT CORRECT: it "
            "guarantees every structure was measured at the SAME grid points, "
            "not that those points are the site of interest. site_from_density "
            "is a by-construction grid that is 7.73 A off-site on apo "
            "TNF-alpha. Check each site's ligand_anchored and "
            "distance_to_donor_ligand_centroid_a."
        ),
    }
    if n < 2:
        out["mdpocket_reason"] = (
            f"{n} structure(s); an ensemble measurement needs at least 2"
        )
        return out

    run = work / "mdpocket"
    shutil.rmtree(run, ignore_errors=True)
    run.mkdir(parents=True, exist_ok=True)

    # ---- superposition, Kabsch on core C-alpha ---------------------------
    #
    # NUMBERING IS ALIGNED BEFORE ANYTHING IS FITTED. Indexing C-alpha by raw
    # author residue number was this file's most serious bug: PDB entries for
    # one protein routinely use different numbering conventions, and TL1A's
    # ensemble carries three (2O0O at 0, five entries at +67, 2RE9 at +71). The
    # fit then matched non-homologous residues and reported
    # "best chain mapping RMSD 14.84 A exceeds 5.0 A; not a superposition" —
    # a message that misdiagnosed its own failure as a conformational one. The
    # same ensemble superposes at 0.51-1.45 A once numbering is aligned. See
    # `_best_numbering_offset` and `_common_core`.
    # ---- THE CORE IS BUILT OVER THE TARGET'S CHAINS ONLY -------------------
    # `_common_core` requires a residue number to be present in EVERY chain of
    # EVERY structure and to name the same residue in all of them. Applied
    # across every chain in the file, an entry carrying an antibody Fab, a
    # receptor ectodomain or a fusion partner poisons the intersection with
    # chains that have nothing to do with the target and nothing to do with each
    # other. Measured on BAFF, which returned:
    #
    #     "only 0 C-alpha positions are shared by every chain of every structure
    #      AND name the same residue in all of them (131 numbers were shared,
    #      131 of them dropped because the entries disagree about which residue
    #      that number is)"
    #
    # 131 shared numbers, all 131 rejected — and the whole by-construction site
    # definition lost, on all three targets in that batch. The identity gate is
    # right and its SCOPE was wrong: two chains that are different proteins are
    # supposed to disagree about residue 131, and asking them to agree is not a
    # test the ensemble can pass. Restricted to the target's own chains it is
    # the test it was meant to be.
    #
    # Falls back to every chain per structure where target chains are unknown,
    # so a run without an accession behaves exactly as before.
    tchains = target_chains or {}
    cas_raw = {
        pid: _filter_chains(_ca_by_chain(p), tchains.get(pid))
        for pid, p in prepped.items()
    }
    names_raw = {
        pid: _filter_chains(_res_names_by_chain(p), tchains.get(pid))
        for pid, p in prepped.items()
    }
    out["core_restricted_to_target_chains"] = {
        "applied": bool(tchains),
        "per_structure": {pid: tchains.get(pid) for pid in ids},
        "_why": (
            "the common core is an intersection over chains, so a chain that is "
            "a different protein empties it. BAFF returned 131 shared numbers "
            "and dropped all 131 on residue-name disagreement, losing the "
            "mdpocket stage on every target in that batch."
        ),
    }
    # THE REFERENCE IS CHOSEN, NOT INHERITED FROM DICT ORDER. It used to be
    # `ids[0]`, and on NLRP3 that was the outlier — see
    # `_select_mdpocket_reference` for the measurement and for what it does and
    # does not fix.
    ref_choice = _select_mdpocket_reference(ids, cas_raw, names_raw)
    ref_pid = ref_choice["reference"]
    out["reference_selection"] = ref_choice
    # `_fit_all` fits onto `pids[0]`, so the chosen reference leads the list.
    ids = [ref_pid] + [p for p in ids if p != ref_pid]
    if not cas_raw[ref_pid]:
        out.update(mdpocket_status="failed",
                   mdpocket_reason=f"no C-alpha in reference {ref_pid}")
        return out
    # The NUMBERING reference, fixed for the whole run. It is allowed to diverge
    # from the FRAME reference below: pass 2 may re-elect the frame reference
    # among the survivors, but every `cas`/`names` dict has already been
    # renumbered onto this one entry and renumbering them again would invalidate
    # the core, the donor alignment and the fits computed against them.
    numbering_ref_pid = ref_pid
    ref_names_pooled = _pdb_resnames_by_seqid(prepped[numbering_ref_pid])

    cas: dict[str, dict] = {}
    names: dict[str, dict] = {}
    numbering: dict[str, dict] = {}
    for pid in ids:
        cas[pid], names[pid], numbering[pid] = _align_numbering(
            cas_raw[pid], names_raw[pid], ref_names_pooled
        )
    out["numbering_alignment"] = {
        "reference": numbering_ref_pid,
        "per_structure": numbering,
        "_why": (
            "C-alpha are matched on residue number, so the entries must first "
            "be put on ONE numbering. Offsets are recovered by voting on "
            "residue-name agreement against the reference, then every core "
            "position must name the same residue in every structure. Fitting "
            "on raw author numbers turned a 0.5 A TL1A superposition into a "
            "reported 14.8-20.1 A refusal and lost the whole mdpocket stage."
        ),
    }

    def _fit_all(pids: list[str]) -> tuple[dict, list[int], dict, list[dict]]:
        """Core over `pids`, then fit each of them. Failures are returned, not
        raised — a structure that will not superpose costs itself only."""
        core_, report_ = _common_core({p: names[p] for p in pids})
        fits_: dict[str, dict] = {}
        dropped_: list[dict] = []
        if len(core_) < 20:
            return fits_, core_, report_, dropped_
        ref_ca_ = cas[pids[0]]
        for pid_ in pids:
            fit = _fit_to_reference(cas[pid_], ref_ca_, core_)
            if fit.get("ok"):
                fits_[pid_] = fit
            else:
                dropped_.append({
                    "pdb_id": pid_,
                    "reason": fit.get("reason"),
                    "rmsd_a": fit.get("rmsd_a"),
                    "n_chains_mapped": fit.get("n_chains_mapped"),
                    "n_mobile_chains": fit.get("n_mobile_chains"),
                    "n_reference_chains": fit.get("n_reference_chains"),
                    "all_mapping_rmsd_a": fit.get("all_mapping_rmsd_a"),
                    "numbering": numbering.get(pid_, {}).get("per_chain"),
                })
        return fits_, core_, report_, dropped_

    fits, core, core_report, dropped = _fit_all(ids)
    passes = [{"pass": 1, "pids": list(ids), "n_core": len(core),
               "core_report": core_report,
               "dropped": [d["pdb_id"] for d in dropped]}]
    if len(core) < 20:
        out.update(
            mdpocket_status="failed",
            mdpocket_reason=(
                f"only {len(core)} C-alpha positions are shared by every chain "
                f"of every structure AND name the same residue in all of them "
                f"({core_report['n_shared_numbers']} numbers were shared, "
                f"{core_report['n_dropped_on_residue_name_disagreement']} of "
                "them dropped because the entries disagree about which residue "
                "that number is). Numbering was aligned first, so this is a "
                "genuine lack of common sequence rather than an offset."
            ),
            numbering_alignment=out["numbering_alignment"],
            core_selection=core_report,
        )
        return out

    # A structure that failed to fit also constrained the core it was fitted
    # against. Recompute over the survivors so a dropped entry does not shrink
    # the measurement it is no longer part of.
    if dropped:
        survivors = [p for p in ids if p not in {d["pdb_id"] for d in dropped}]
        if len(survivors) >= 2:
            # AND THE REFERENCE IS RE-ELECTED, NOT INHERITED. A drop changes
            # which structure is central to what is left, and the pass-1
            # reference is by definition the one thing every drop was measured
            # against. Re-running the selection over the survivors is what turns
            # "the reference dropped everything" into a recoverable ensemble.
            resel = _select_mdpocket_reference(survivors, cas_raw, names_raw)
            new_ref = resel.get("reference") or survivors[0]
            survivors = [new_ref] + [p for p in survivors if p != new_ref]
            fits2, core2, report2, dropped2 = _fit_all(survivors)
            passes.append({
                "pass": 2, "pids": survivors, "n_core": len(core2),
                "core_report": report2,
                "reference": new_ref,
                "reference_selection": resel,
                "dropped": [d["pdb_id"] for d in dropped2],
                "_why": (
                    "pass 1 dropped a structure, and that structure had also "
                    "been constraining the common core AND been the frame every "
                    "other structure was fitted onto. Pass 2 recomputes the core "
                    "over the survivors only and re-elects the reference among "
                    "them. Numbering is NOT realigned — every dict here is "
                    "already on the pass-1 numbering reference and re-aligning "
                    "would invalidate the core it is measured against."
                ),
            })
            if len(core2) >= 20 and fits2:
                fits, core, core_report = fits2, core2, report2
                dropped = dropped + dropped2
                ref_pid = new_ref

    kept_ids = [p for p in ids if p in fits]
    n_kept = len(kept_ids)
    out["frames_dropped"] = {
        "n_input": n,
        "n_kept": n_kept,
        "n_dropped": len(dropped),
        "dropped": dropped,
        "kept": kept_ids,
        "_why": (
            "A structure that will not superpose is DROPPED and recorded here, "
            "not fatal. One 4-chain assembly matched against a 2-chain "
            "reference (IRAK4 6UYA, 23.87 A) used to abort the entire ensemble "
            "and cost a full re-run. THE REFUSAL IS CORRECT; ABORTING IS NOT."
        ),
        "_not_the_same_as_mdpockets_silent_drop": (
            "mdpocket ALSO drops frames, silently, and still exits 0 — a "
            "5-entry list with one bad path prints 'Identified 5 snapshots', "
            "then an error, then 1/4..4/4, and returns 0, renormalising every "
            "frequency in the grid over 4. That failure is still refused, by "
            "_assert_frame_count, and this deliberate drop does not weaken it: "
            "the assertion is made against n_kept (the frames actually "
            "submitted to mdpocket), so a frame we dropped on purpose and a "
            "frame mdpocket lost on its own remain distinguishable. "
            "frame_count_check carries n_input, n_submitted and n_processed "
            "separately for exactly this reason."
        ),
    }
    if n_kept < 2 or (dropped and n_kept < MDPOCKET_MIN_N_AFTER_DROPS):
        out.update(
            mdpocket_status="failed",
            mdpocket_reason=(
                f"{len(dropped)} of {n} structures would not superpose onto "
                f"{ref_pid} and were dropped, leaving {n_kept}. THE REFERENCE "
                "WAS CHOSEN, NOT INHERITED — lowest median RMSD to the rest, "
                "re-elected after the drop — so this is a statement about the "
                "dropped structures and not about which entry happened to come "
                "first; see reference_selection.candidates for every "
                "candidate's median and for what each one would have kept. An "
                "ensemble "
                f"measurement needs at least 2 structures, and at least "
                f"{MDPOCKET_MIN_N_AFTER_DROPS} once a drop has occurred — a CV "
                "over the two survivors of a partial ensemble is a statistic, "
                "not a measurement of reproducibility. Refusing rather than "
                "reporting it. Dropped: "
                + "; ".join(f"{d['pdb_id']}: {d['reason']}" for d in dropped)
            ),
        )
        return out

    sup_paths: dict[str, Path] = {}
    for pid in kept_ids:
        fit = fits[pid]
        dst = run / f"sup_{pid}.pdb"
        _write_superposed(prepped[pid], dst, fit.pop("_R"), fit.pop("_t"),
                          fit["chain_map"])
        sup_paths[pid] = dst

    ids, n = kept_ids, n_kept
    out["n_structures"] = n
    out["n_structures_submitted"] = n
    out["superposition"] = {
        "reference": ref_pid,
        "reference_selection": ref_choice,
        "numbering_reference": numbering_ref_pid,
        "_reference_note": (
            "`reference` is the FRAME every structure was fitted onto and is "
            "chosen by lowest median RMSD to the rest; `numbering_reference` is "
            "the entry whose residue numbering the whole ensemble was put on, "
            "and it is fixed at the pass-1 choice. They differ only when a drop "
            "caused the frame reference to be re-elected."
        ),
        "n_core_ca_positions": len(core),
        "core_residue_range": [core[0], core[-1]],
        "core_selection": core_report,
        "numbering_alignment_passes": passes,
        "per_structure": fits,
        "_note": (
            "mdpocket does NOT superpose; it assumes MD frames already in a "
            "common frame. Deposited entries are not. Chain mappings are "
            "searched, not assumed: all_mapping_rmsd_a shows the split — on a "
            "C3 trimer the three cyclic mappings agree to ~0.03 A and the "
            "three anticyclic ones land ~22 A out. Residue NUMBERING is "
            "aligned before any of that; see numbering_alignment."
        ),
    }

    # ---- the site donor's ligand, moved into the common frame ------------
    lig_common = None
    donor_fit = None
    donor_numbering = None
    if donor_ligand_xyz and donor_prepped and donor_prepped.exists():
        # The donor is a separate deposition and needs the SAME numbering
        # alignment as the ensemble; a site donor at a different offset would
        # fail to fit for the same reason 2QE3 did, and take the whole
        # ligand-anchored site with it.
        donor_ca, donor_names, donor_numbering = _align_numbering(
            _ca_by_chain(donor_prepped),
            _res_names_by_chain(donor_prepped),
            ref_names_pooled,
        )
        # Only the chains that line the ligand copy being transferred. Without
        # this the second crystallographic dimer of 2AZ5 joins the fit and the
        # best mapping is 17.3 A — no transfer, no ligand-anchored site.
        dchains = [c for c in (donor_chains or sorted(donor_ca)) if c in donor_ca]
        dcore = sorted(
            set.intersection(*[set(donor_ca[c]) for c in dchains]) & set(core)
        ) if dchains else []
        # Same identity guard as the ensemble core: a shared number is not a
        # shared residue.
        dcore = [
            r for r in dcore
            if all(donor_names[c].get(r) == names[ref_pid][cc].get(r)
                   for c in dchains for cc in names[ref_pid])
        ]
        if len(dcore) >= 20:
            donor_fit = _fit_to_reference(donor_ca, cas[ref_pid], dcore, dchains)
            if donor_fit.get("ok"):
                R, t = donor_fit.pop("_R"), donor_fit.pop("_t")
                lig_common = np.asarray(donor_ligand_xyz, dtype=float) @ R + t
        out["site_donor"] = {
            "pdb_id": donor_pid,
            "chains_used": dchains,
            "n_ligand_atoms": len(donor_ligand_xyz),
            "n_core_ca_positions": len(dcore),
            "numbering_alignment": donor_numbering,
            "fit": ({k: v for k, v in (donor_fit or {}).items()
                     if not k.startswith("_")} or None),
            "transferred": lig_common is not None,
            "transfer_failed_reason": (
                None if lig_common is not None else
                (f"only {len(dcore)} core C-alpha positions are shared with the "
                 f"ensemble reference {ref_pid} and name the same residue "
                 "(20 needed)")
                if len(dcore) < 20
                else (donor_fit or {}).get("reason")
            ),
        }

    # ---- density mode ----------------------------------------------------
    listfile = run / "list.txt"
    listfile.write_text("\n".join(p.name for p in sup_paths.values()) + "\n")
    try:
        proc = _run_mdpocket(["--pdb_list", listfile.name, "-o", "dens"],
                             run, timeout)
    except (OSError, subprocess.SubprocessError) as exc:
        out.update(mdpocket_status="failed",
                   mdpocket_reason=f"{type(exc).__name__}: {exc}")
        return out
    time_txt = run / "time.txt"
    try:
        if not time_txt.exists():
            raise RuntimeError(
                f"mdpocket wrote no time.txt (exit {proc.returncode}): "
                f"{proc.stdout.decode(errors='replace')[-400:]}"
            )
        n_ran = len([ln for ln in time_txt.read_text().splitlines() if ln.strip()])
        _assert_frame_count(n_ran, n, "time.txt lines")
    except RuntimeError as exc:
        out.update(mdpocket_status="failed", mdpocket_reason=str(exc))
        return out

    try:
        freq, counts, origin, spacing = _read_dx(run / "dens_freq.dx")
    except (OSError, RuntimeError) as exc:
        out.update(mdpocket_status="failed",
                   mdpocket_reason=f"frequency grid unreadable: {exc}")
        return out
    pts = _dx_points(counts, origin, spacing)

    out["frame_count_check"] = {
        # THREE NUMBERS, NOT TWO. n_input is what the caller asked for,
        # n_submitted is what survived OUR superposition refusal, n_processed is
        # what mdpocket actually ran. The assertion is n_processed == n_submitted
        # and it is undiminished: a frame we dropped on purpose is recorded in
        # `frames_dropped`, a frame mdpocket lost on its own still fails here.
        "n_input_structures": len(prepped),
        "n_submitted_to_mdpocket": n,
        "n_processed": n_ran,
        "n_dropped_by_us_before_submission": len(prepped) - n,
        "passed": True,
        "_why": (
            "mdpocket drops unreadable frames silently and still exits 0, and "
            "the frequency grid is renormalised over the frames that ran, so a "
            "dropped structure inflates every frequency in the grid. Our own "
            "deliberate drops happen BEFORE submission and are listed in "
            "frames_dropped; they never enter this comparison."
        ),
    }
    out["frequency"] = {
        "reported": n >= MDPOCKET_MIN_N_FOR_FREQUENCY,
        "quantum": round(1.0 / n, 4),
        "min_n_to_report": MDPOCKET_MIN_N_FOR_FREQUENCY,
        "refusal_reason": (
            None if n >= MDPOCKET_MIN_N_FOR_FREQUENCY else
            f"frequency is quantised at 1/N; with N={n} the only attainable "
            f"values are {{0, {round(1.0 / n, 2)}, ..., 1.0}}, so any "
            f"fractional occupancy quoted from this ensemble would be reporting "
            f"the grid resolution and not a measurement. Refusing. "
            f"n_of_n presence below is exact at any N and is reported."
        ),
    }

    # ---- the two site definitions ----------------------------------------
    sites: dict[str, dict] = {}
    all_present = freq >= 1.0 - 1e-9
    if all_present.any():
        sel_pts = pts[all_present]
        cluster = _largest_grid_cluster(sel_pts, spacing)
        if len(cluster):
            cav = sel_pts[cluster]
            sel = run / "sel_density.pdb"
            _write_probe_pdb(cav, sel)
            entry: dict = {
                "definition": (
                    f"largest connected cluster of grid points at which a "
                    f"pocket is present in ALL {n} structures"
                ),
                "n_of_n": n,
                "centroid": [round(float(x), 2) for x in cav.mean(0)],
                # NOT the ligand site, by construction. This is the most
                # persistent cavity, which is a different question. Stated as a
                # field rather than only in prose so a consumer cannot read the
                # volume without it.
                "ligand_anchored": False,
                "_is_this_the_pocket": (
                    "NO, not necessarily. This is the most PERSISTENT cavity in "
                    "the ensemble, which is not automatically the LIGAND's site. "
                    "On the apo TNF-alpha ensemble its centroid sits 7.73 A from "
                    "the transferred SPD304 ligand — it is the on-axis cavity, "
                    "and it is precisely the pocket the retracted residue-number "
                    "matcher reported as 'the SPD304 site'. Check "
                    "distance_to_donor_ligand_centroid_a before treating this as "
                    "the site the dossier is asking about; when a "
                    "site_from_ligand entry is also present, that one is the "
                    "ligand site and this one is not."
                ),
            }
            entry.update(
                _ligand_distance_fields(
                    cav.mean(0), lig_common, donor_pid, donor_fit
                )
            )
            try:
                entry.update(
                    _mdpocket_characterise(run, listfile, sel, "site_density",
                                           ids, timeout)
                )
                sites["site_from_density"] = entry
            except (RuntimeError, OSError, subprocess.SubprocessError) as exc:
                sites["site_from_density"] = {**entry, "error": str(exc)}

    if lig_common is not None:
        d = np.linalg.norm(pts[:, None, :] - lig_common[None, :, :], axis=-1).min(1)
        near = pts[d <= MDPOCKET_SITE_CUTOFF_A]
        if len(near):
            sel = run / "sel_ligand.pdb"
            _write_probe_pdb(near, sel)
            entry = {
                "definition": (
                    f"grid points within {MDPOCKET_SITE_CUTOFF_A} A of the "
                    f"{donor_pid} ligand, transferred by superposition"
                ),
                "centroid": [round(float(x), 2) for x in near.mean(0)],
                # THIS is the site the dossier is asking about. Stated as a
                # field so the two entries can be told apart programmatically
                # and not only by their key name.
                "ligand_anchored": True,
                "_is_this_the_pocket": (
                    "YES. These grid points are defined BY the transferred "
                    "ligand, so this is the ligand site by construction. When "
                    "a site_from_density entry is also present it is a "
                    "different pocket unless its "
                    "distance_to_donor_ligand_centroid_a is small."
                ),
                "_zero_is_an_answer": (
                    "0.00 A^3 in a structure means the site is NOT OPEN in that "
                    "structure at these exact grid points. That is a refusal, "
                    "not a failure: on the five apo TNF-alpha entries the "
                    "SPD304 site returns 0.00 in four of five, which is the "
                    "correct result — the channel is occluded by the third "
                    "subunit. Do not read it as a broken run."
                ),
            }
            try:
                entry.update(
                    _mdpocket_characterise(run, listfile, sel, "site_ligand",
                                           ids, timeout)
                )
                sites["site_from_ligand"] = entry
            except (RuntimeError, OSError, subprocess.SubprocessError) as exc:
                sites["site_from_ligand"] = {**entry, "error": str(exc)}

    out["sites"] = sites
    out["_reproducibility"] = (
        "fpocket estimates pocket volume by Monte Carlo, and mdpocket inherits "
        "it: three identical reruns of this exact ensemble gave site_from_density "
        "volumes of 280.6/276.1/274.6 A^3 on the same structure and CVs of "
        "12.1/11.3/10.8%. About 1 percentage point of the reported CV is the "
        "method's own noise. Do not read a CV difference smaller than that as a "
        "difference between sites."
    )
    out["_site_from_density_caveat"] = (
        "site_from_density is the most PERSISTENT cavity, which is not "
        "automatically the LIGAND's site. On the apo TNF-alpha ensemble its "
        "centroid sits ~7.7 A from the transferred SPD304 ligand — it is the "
        "on-axis cavity, and it is precisely the pocket the residue-number "
        "matcher used to report as 'the SPD304 site'. Read "
        "distance_to_donor_ligand_centroid_a before treating the two as one."
    )
    if not sites:
        out.update(
            mdpocket_status="failed",
            mdpocket_reason=(
                "the density run produced no site definition: no grid point is "
                "occupied in all structures and no ligand was transferable"
            ),
        )
        return out
    failed = [k for k, v in sites.items() if v.get("error")]
    out["mdpocket_status"] = "failed" if len(failed) == len(sites) else "ok"
    if failed:
        out["mdpocket_reason"] = "; ".join(f"{k}: {sites[k]['error']}" for k in failed)
    return out


@app.function(cpu=4.0, timeout=1800)
def pocket_scan(
    pdb_ids: list[str],
    chains: dict[str, list[str]] | None = None,
    ligand_codes: list[str] | None = None,
    site_residues: list[int] | None = None,
    uniprot_accession: str | None = None,
    partner_structures: list[str] | None = None,
    mdpocket_site_donor: str | None = None,
    run_disorder: bool = True,
    run_cryptic: bool = True,
    run_mdpocket: bool = True,
    predicted_structures: dict[str, str] | None = None,
) -> dict:
    """Scan an ensemble at every clustering value and report the spread.

    ADDED KEYWORD PARAMETERS (the four original positional ones are unchanged):

    uniprot_accession
        Drives the disorder stage. STRONGLY PREFERRED over the structure-derived
        fallback: a deposited construct is the ordered part of the protein by
        selection, so a sequence lifted out of a PDB entry understates disorder.
        Validation: P01106 (MYC) ~0.83, P24941 (CDK2) ~0.00.
    partner_structures
        PDB IDs of complexes containing the BINDING PARTNER, e.g. `["3ALQ"]`
        for TNF-alpha + TNFR2. Turns "is this pocket orthosteric?" from an
        assumption into a measurement. Target and partner chains are separated
        by sequence, not by chain letter.
    mdpocket_site_donor
        A holo PDB ID used ONLY to define the site for the mdpocket stage, not
        added to the ensemble. This is how a pure-apo ensemble gets a
        ligand-anchored site: `pdb_ids=[5 apo TNF entries],
        mdpocket_site_donor="2AZ5", ligand_codes=["307"]`. When omitted, a holo
        structure already inside `pdb_ids` donates the site.
    run_disorder, run_cryptic, run_mdpocket
        Switches for the optional stages. Every one of them is non-fatal
        regardless: each emits `<stage>_status` in {ok, failed, not_run} with a
        reason, and none of them can kill the fpocket result.
    predicted_structures
        `{label: mmCIF text}` for a target with NO experimental structure and no
        usable homolog, folded off the default path by
        `structure-select/predicted_structure_fallback.py` (ESMFold). Each is
        pre-seeded to disk so the ordinary fpocket/mdpocket scan runs on the
        model UNCHANGED, added to a LOCAL working id list (the caller's `pdb_ids`
        is never mutated), and stamped `structure_origin: "esmfold_predicted"`.
        `uniprot_accession` is REQUIRED when this is non-empty — a predicted file
        has no `_struct_ref`, so the accession cannot be recovered from it and a
        clear `ValueError` is raised if it is missing. A predicted entry NEVER
        donates a site signature and NEVER acts as a holo or mdpocket site donor
        (it carries no ligand); any pocket found on one is a MODEL pocket, listed
        in `predicted_structures_used`. See CLAUDE.md rule 4c.

    Returns volume with its across-structure spread as the PRIMARY number, and
    druggability only as a range.

    SAME-SITE TRACKING. A spread is only a measurement if every value describes
    the SAME site. On an apo structure there is no ligand to anchor to, so
    without a site signature this falls back to "most druggable pocket
    anywhere" — and pooling that across an ensemble compares different pockets
    in different places, which measures nothing.

    The signature is a set of residue numbers, matched chain-agnostically. It
    comes from one of two places:

      * `site_residues`, supplied by the caller; or
      * automatically, from the ligand site of any HOLO structure in the same
        run. Put one holo entry in an otherwise apo ensemble and every apo
        structure is then scored at the site the holo one points at. That is
        the KRAS 6OIM/4OBE case and the TNF-alpha 2AZ5/apo case.

    THE SIGNATURE PATH IS NOT VALID ON A HOMO-OLIGOMER. Chain-agnostic residue
    numbers are triplicated by a homotrimer's protomers, so a C3-symmetric site
    is unresolvable in principle and the match lands wherever the numbers land.
    When `_homo_oligomer` fires, the basis is reported as
    `site_signature_unreliable_homooligomer` and those values must not be pooled
    as one site. See `_homo_oligomer` for the measurement.

    `ligand_codes` is an OVERRIDE, NOT A REQUIREMENT. When a structure carries a
    drug-like ligand and the caller named no code that matches it, that
    structure's OWN drug-like ligand anchors its site. Not doing this made the
    answer a function of the CLI: passing one code across four TNF structures
    left three of them falling back to the weaker signature path, and passing
    all four moved 7JRA from druggability 0.000 / 306.9 A^3 to 0.926 /
    1542.9 A^3 — same structures, same clustering, near-total reversal driven by
    an argument. `site_anchor_ligand` and `site_anchor_ligand_source` record what
    was actually used and where it came from.

    Read `site_pocket_selected_by` before quoting any spread: only
    `ligand_site_jaccard` is same-site without qualification.
    """
    work = Path("/tmp/pockets")
    work.mkdir(parents=True, exist_ok=True)
    chains = chains or {}
    # ---- PREDICTED STRUCTURES, pre-seeded to disk (the least-invasive path) --
    # A structure-less target folded by predicted_structure_fallback.py arrives
    # as {label: mmCIF text}. Writing each to work/{label}.cif (and a matching
    # header copy) BEFORE anything runs makes `_fetch`/`_fetch_header` return it
    # from disk, so the whole ordinary scan runs on the model with no code path
    # aware it is synthetic — and no doomed RCSB call for an id that is not a PDB
    # entry. A predicted file has no `_struct_ref`, so its accession cannot be
    # recovered from the file and MUST be supplied.
    predicted_structures = predicted_structures or {}
    predicted_ids: set[str] = set(predicted_structures)
    if predicted_structures and not (uniprot_accession or "").strip():
        raise ValueError(
            "predicted_structures requires uniprot_accession: a predicted "
            "structure carries no _struct_ref, so the target accession cannot "
            "be recovered from the file and must be passed explicitly."
        )
    for _label, _cif_text in predicted_structures.items():
        (work / f"{_label}.cif").write_text(_cif_text)
        # A header copy too, so _fetch_header also short-circuits to disk rather
        # than trying (and failing) to GET files.rcsb.org/header/<label>.cif. It
        # carries no _struct_ref, so accession resolution off it returns
        # no_unp_refs and _target_chains falls back to scoring the modelled
        # chain(s) — which for a single-chain fold is exactly the target.
        (work / f"{_label}_header.cif").write_text(_cif_text)
        # structure_source is the file-provenance field (assembly1 / asu / ...);
        # "esmfold_predicted" is the honest value for a model and it flows into
        # method.source_used without special-casing.
        (work / f"{_label}.source").write_text("esmfold_predicted")
    # The id list the per-structure scan iterates: real entries plus predicted
    # labels. A LOCAL list — the caller's pdb_ids is never mutated.
    scan_ids = list(pdb_ids) + list(predicted_structures)
    # ONE chemical-component source for the whole run, so a comp_id is looked up
    # once however many entries carry it. A hard import: if `ligand_filter` is
    # not in the image, the holo/apo call has no basis at all, and failing
    # loudly here is the only honest option — a fallback to a size threshold
    # would reinstate the exact bug it replaces, silently.
    chemcomp_src = _chemcomp_source()
    results: dict[str, dict] = {}
    # Carried across the per-structure loop so the four stages after fpocket can
    # reuse exactly what was scored — same files, same chains, same frame.
    prepped_by_pid: dict[str, Path] = {}
    chains_by_pid: dict[str, list[str]] = {}
    tgt_chains_by_pid: dict[str, list[str]] = {}
    cif_by_pid: dict[str, Path] = {}
    seqs_by_pid: dict[str, str] = {}
    vert_by_pid: dict[str, dict[str, Path]] = {}
    holo_by_pid: dict[str, dict] = {}
    site_signature: set[str] = {str(r) for r in (site_residues or [])}
    signature_source = "caller" if site_signature else None
    signature_donor_homo: dict | None = None
    signature_n_residues_in = len(site_residues or [])

    # ---- WHICH PROTEIN IS THE TARGET, resolved before anything reads a chain
    # Three separate wrong answers all came from not knowing this: the target
    # sequence picked by chain length (S1PR1's G-beta-1), the homo-oligomer
    # guard firing on a partner's homodimer (8G94's CD69 pair, 25 and 27
    # residues, which disqualified a site that by hand matches the holo pockets
    # at Jaccard 0.79-0.94), and disorder measured on whichever chain happened
    # to be longest. The accession is already an input, and every entry declares
    # its own in `_struct_ref`.
    #
    # MOVED ABOVE THE TWO SITE-DONOR BLOCKS, and that move is load-bearing. Both
    # donors build a `ligand_filter.StructureContext`, the context needs the
    # target accession to tell "bonded to the target" from "bonded to another
    # polymer", and this block only reads headers so it depends on nothing
    # below it. Resolved late, the SIGNATURE DONOR LOOP — the one that picks the
    # component whose contact shell becomes the site signature — ran with no
    # context at all, and on 8QFZ that made `LFI` `druglike` and let the
    # crosslinker of a bicyclic peptide define the site.
    target_accession = (uniprot_accession or "").strip() or None
    accession_basis = "caller" if target_accession else None
    accession_counts: dict[str, int] = {}
    accession_by_pid: dict[str, list[str]] = {}
    for pid in pdb_ids:
        try:
            found = _uniprot_from_cif(_fetch_header(pid, work))
        except Exception:  # noqa: BLE001, S112
            continue
        accession_by_pid[pid] = found
        for a in found:
            accession_counts[a] = accession_counts.get(a, 0) + 1
    # Predicted labels were skipped by the header loop above — they are not in
    # pdb_ids and carry no _struct_ref anyway — so seed their accession directly
    # from the caller. target_accession is guaranteed present by the ValueError
    # guard whenever predicted_ids is non-empty.
    for _label in predicted_ids:
        accession_by_pid[_label] = [target_accession]
    if not target_accession and accession_counts:
        # THE ACCESSION PRESENT IN THE MOST ENTRIES. The target is the protein
        # the ensemble was assembled around, so it appears in every entry;
        # partners, fusion chaperones and crystallisation scaffolds vary. Taking
        # the first accession of the first entry would pick whatever the
        # depositor listed first.
        top = max(accession_counts.values())
        winners = [a for a in accession_counts if accession_counts[a] == top]
        target_accession = winners[0]
        accession_basis = (
            f"mmcif:_struct_ref — present in {top} of {len(pdb_ids)} entries"
            + ("" if len(winners) == 1 else
               f"; AMBIGUOUS, tied with {', '.join(winners[1:])}, first taken. "
               "Pass uniprot_accession to disambiguate.")
        )

    # ---- an out-of-ensemble holo reference, resolved SECOND ----------------
    # `mdpocket_site_donor` names a holo structure deliberately NOT measured —
    # the five apo TNF-alpha entries with 2AZ5 donating the site. It is
    # resolved before anything else because it serves three stages: it donates
    # the site SIGNATURE for the fpocket pass below, it is the holo half of the
    # cryptic comparison, and it defines the mdpocket site. Resolving it late
    # meant a pure-apo run fell back to "most druggable pocket anywhere" in the
    # fpocket pass while a perfectly good holo reference sat one argument away.
    donor: dict | None = None
    donor_error: str | None = None
    if mdpocket_site_donor and mdpocket_site_donor not in pdb_ids:
        try:
            dcif = _fetch(mdpocket_site_donor, work)
            dst, _dmiss, _dren = _load(dcif)
            dprep, dchains, _ddrop = _prep(dst, work, mdpocket_site_donor, None)
            dligs, _dholo = _ligands(
                dst, chemcomp_src,
                _structure_context(mdpocket_site_donor, work, target_accession),
            )
            dcomp = next(
                (lig["comp_id"] for lig in dligs
                 if ligand_codes and lig["comp_id"] in ligand_codes),
                next((lig["comp_id"] for lig in dligs if lig["druglike"]), None),
            )
            if not dcomp:
                raise RuntimeError(
                    f"{mdpocket_site_donor} carries no drug-like ligand, so it "
                    "cannot donate a site"
                )
            dsite, dcopy = _ligand_site(dst, dcomp, dchains)
            dkeep, dcounts = _ligand_contact_chains(dst, dcomp, dcopy)
            donor = {
                "pdb_id": mdpocket_site_donor,
                "cif": dcif,
                "prepped": dprep,
                "in_ensemble": False,
                "comp_id": dcomp,
                "ligand_chain": dcopy.split("/")[0] if dcopy else None,
                "ligand_copy": dcopy,
                "site_chains": dkeep,
                "site_chain_atom_counts": dcounts,
                "site_residues": dsite,
                "homo_oligomer": _homo_oligomer(dst),
                # The donor's OWN target chains, for the signature filter below.
                # None means "cannot filter" and must never be read as "filter
                # to nothing" — see `_target_polymer_chains`.
                "target_polymer_chains": _target_polymer_chains(
                    dst, _fetch_header(mdpocket_site_donor, work),
                    _dren, target_accession,
                ),
                "ligand_xyz": [
                    [a.pos.x, a.pos.y, a.pos.z]
                    for chain in dst[0] for res in chain
                    if res.het_flag == "H" and res.name == dcomp
                    and (dcopy is None or
                         f"{chain.name}/{res.seqid.num}" == dcopy)
                    for a in res
                ],
            }
        except Exception as exc:  # noqa: BLE001
            donor_error = (
                f"site donor {mdpocket_site_donor}: {type(exc).__name__}: {exc}"
            )

    # RESIDUE NUMBERS FROM A DIFFERENT POLYMER ARE NOT THIS PROTEIN'S RESIDUE
    # NUMBERS, and the signature discards chain identity by design. The bicyclic
    # peptide in 8QFZ numbers from 1, so 9 of the 13 residues in `LFI`'s contact
    # shell carry numbers 11-22 that ALSO exist on TSLP and mean something else
    # entirely; matched chain-agnostically they land on a different part of the
    # protein. Chain-agnostic matching is right for an inter-subunit site on ONE
    # protein (TNF-alpha's axial channel) and wrong across two DIFFERENT
    # proteins, and the homo-oligomer guard cannot see the difference because
    # nothing collapses — `collapsed_by` was 0 on exactly this failure.
    signature_foreign_dropped = 0
    signature_foreign_residues: list[str] = []
    if not site_signature and donor and donor["site_residues"]:
        _tgt_ch = donor.get("target_polymer_chains")
        _res = donor["site_residues"]
        signature_foreign_residues = (
            [r for r in _res if r.split("/")[0] not in _tgt_ch] if _tgt_ch else []
        )
        signature_foreign_dropped = len(signature_foreign_residues)
        _res_t = [r for r in _res if r.split("/")[0] in _tgt_ch] if _tgt_ch else _res
        site_signature = {r.split("/")[-1] for r in _res_t}
        signature_source = f"{donor['pdb_id']}:{donor['comp_id']} (site donor)"
        signature_n_residues_in = len(_res)
        signature_donor_homo = donor["homo_oligomer"]

    # A holo structure anywhere in the ensemble donates its ligand site as the
    # signature for the apo ones, so order the pass holo-first.
    if not site_signature:
        for pid in pdb_ids:
            if pid in predicted_ids:
                # DEFENSIVE: a model carries no ligand, so it can never be a
                # signature donor — and this loop iterates pdb_ids, which never
                # holds a predicted label — but guard explicitly so a future
                # refactor that widens the loop cannot let a model define a site.
                continue
            try:
                raw = _fetch(pid, work)
                st, _chains_avail, _renames = _load(raw)
                # WITH THE STRUCTURAL CONTEXT. This loop is where
                # `dl[0]` picks the component whose contact shell becomes the
                # site signature, so a crosslinker misread as `druglike` defines
                # the site for the whole ensemble. On 8QFZ that component was
                # `LFI`, and with the context it is `polymer_conjugate` and does
                # not reach `dl` at all.
                ligs, _sholo = _ligands(
                    st, chemcomp_src,
                    _structure_context(pid, work, target_accession),
                )
                dl = [lig for lig in ligs if lig["druglike"]]
                if not dl:
                    continue
                comp = next(
                    (lig["comp_id"] for lig in ligs
                     if ligand_codes and lig["comp_id"] in ligand_codes),
                    dl[0]["comp_id"],
                )
                res, _copy = _ligand_site(st, comp, None)
                tgt_ch = _target_polymer_chains(
                    st, _fetch_header(pid, work), _renames, target_accession
                )
                foreign = (
                    [r for r in res if r.split("/")[0] not in tgt_ch]
                    if tgt_ch else []
                )
                res_t = (
                    [r for r in res if r.split("/")[0] in tgt_ch] if tgt_ch else res
                )
                if res:
                    site_signature = {r.split("/")[-1] for r in res_t}
                    signature_source = f"{pid}:{comp}"
                    signature_n_residues_in = len(res)
                    signature_foreign_dropped = len(foreign)
                    signature_foreign_residues = foreign
                    # How badly the donor site collapses is itself the finding:
                    # 19 residues -> 11 numbers on a homotrimer.
                    signature_donor_homo = _homo_oligomer(st)
                    break
            except Exception:  # noqa: BLE001, S112
                continue

    # THE MEASURED COLLAPSE, computed once. This is what decides whether a
    # residue-number signature can identify a site at all — not how many chains
    # the donor has. 19 residues -> 11 numbers on a homotrimer is the failure;
    # 23 -> 23 is not a failure at all. See the basis selection below.
    signature_collapsed_by = max(
        0, signature_n_residues_in - len(site_signature)
    )

    for pid in scan_ids:
        # `scan_ids` is `pdb_ids` plus any predicted labels — the per-structure
        # scan is the one place a model participates, so the model's pocket
        # geometry is measured exactly like a deposited entry's.
        # One unfetchable structure must not lose the whole ensemble, and the
        # reason must survive the trip back: exceptions holding open file
        # handles cannot pickle, so Modal replaces them with an opaque
        # SerializationError. Record the failure as data instead.
        stage = "fetch"
        try:
            cif = _fetch(pid, work)
            st, missing_res, renamed = _load(cif)
            # Everything below reads the one structure object loaded above, so
            # chain IDs and residue numbers are the same in the fpocket input,
            # the ligand list, the ligand site and the missing-residue list.
            stage = "prepare"
            want = (
                sorted({renamed.get(c, c) for c in chains[pid]})
                if chains.get(pid)
                else None
            )
            # ---- RESOLVED BEFORE PREP, because prep now has a decision to make
            # The chains fpocket is given are no longer "every polymer chain":
            # a polymer that is the LIGAND has to come out (rule 4 strips every
            # ligand before scoring, and a Bicycle peptide is a ligand). Telling
            # a ligand chain from a partner chain needs the accession mapping
            # and the ligand verdicts, so both are resolved first and `_prep`
            # runs last. `avail_chains` is exactly what `_prep` would have kept.
            avail_chains = sorted(
                c for c in _chain_monomer_counts(st) if not want or c in want
            )
            ctx = _structure_context(pid, work, target_accession)
            ligs, holo_call = _ligands(st, chemcomp_src, ctx)
            # Target chains by ACCESSION. The homo-oligomer guard asks whether
            # the SITE SIGNATURE is ambiguous, so a partner's homodimer is not
            # its business: 8G94's CD69 pair (25 and 27 residues) tripped it and
            # disqualified an apo structure whose rank-1 pocket matches the holo
            # pockets at Jaccard 0.79/0.94/0.94.
            chain_acc, acc_status = _chain_accessions(
                _fetch_header(pid, work), renamed,
                [c.name for c in st[0]],
            )
            tgt_info = _target_chains(
                chain_acc, target_accession, avail_chains, acc_status
            )
            tgt_chains, tgt_basis = tgt_info["chains"], tgt_info["basis"]
            if tgt_info["refuse"]:
                # FAIL CLOSED. The entry says which proteins it contains and the
                # target is not one of them; scoring it measures a different
                # molecule. See `_target_chains`.
                results[pid] = {
                    "error": f"entry does not contain {target_accession}",
                    "stage": "target_chain_resolution",
                    "tier": "none",
                    "by_clustering": {},
                    "target_chains": tgt_info,
                    "chain_accessions": chain_acc,
                    "_why": tgt_info["reason"],
                }
                continue
            # ---- THE POLYMER LIGAND, AND THE PAIRED MEASUREMENT ------------
            # What the site is anchored on, for the control's pocket matching:
            # the caller's code, else this entry's own drug-like ligand, else —
            # and this is the 8QFZ case — the covalent constituent of the
            # polymer ligand, which is no longer `druglike` and therefore no
            # longer anchors anything by itself.
            _dl = [lig for lig in ligs if lig["druglike"]]
            _anchor = (
                next((lig["comp_id"] for lig in ligs
                      if ligand_codes and lig["comp_id"] in ligand_codes), None)
                or (_dl[0]["comp_id"] if _dl else None)
                or next(iter(holo_call.get("polymer_conjugates") or []), None)
            )
            polymer_control = _polymer_ligand_control(
                st, work, pid, want, holo_call, ctx, tgt_chains,
                tgt_info["verified"], want, renamed, _anchor,
            )
            lig_chains = (polymer_control or {}).get("polymer_ligand_chains") or []
            prepped, used_chains, dropped_chains = _prep(
                st, work, pid, want, drop_chains=lig_chains
            )
            homo = _homo_oligomer(st, tgt_chains)
            # THE FILE'S OWN CHAIN COUNT, beside the target's. The guard runs
            # over target chains, which is right, but reporting only that count
            # made a 60-chain file announce `n_polymer_chains: 10`. Two numbers,
            # both true, neither able to be mistaken for the other.
            homo["n_polymer_chains_in_file"] = sum(
                1 for c in st[0]
                if any(r.het_flag == "A" and len(r) for r in c)
            )
            homo["n_polymer_chains_scored"] = len(used_chains)
            # UniProt's own binding/active-site features, and the alignment
            # needed to compare them against author residue numbers. Fetched
            # once per structure and cached per accession for the container.
            unp_sites, unp_sites_ok = (
                _uniprot_functional_sites(target_accession)
                if target_accession else ({}, False)
            )
            unp_offsets = _chain_unp_offsets(
                _fetch_header(pid, work), renamed, [c.name for c in st[0]]
            )
            # A LARGE ASSEMBLY IS FLAGGED WHETHER OR NOT IT HAPPENED TO FIT.
            # The only thing that used to notice one was `_load` running out of
            # single-character chain IDs, which is an implementation limit at 62
            # and not a scientific boundary — so 1OQE, 1OQD and 4V46 refused
            # while 1JH5, ALSO a 60-mer, sailed through and produced 378 pockets
            # and a selected pocket 60.28 A from the protein centre. Same shape
            # of file, opposite handling, and the one that did not refuse gave
            # plausible numbers with no error anywhere.
            # AND IT REFUSES, the way 1OQE, 1OQD and 4V46 already did. Those
            # three refused only because `_load` ran out of single-character
            # chain IDs — an implementation limit at 62, not a judgement — so
            # 1JH5, ALSO a 60-mer, sailed through and produced 378 pockets with
            # the selected one 60.28 A from the protein centre, plus a
            # `n_polymer_chains: 10` that was counting a tenth of the file.
            # Same shape, opposite handling, and the one that did not refuse
            # returned plausible numbers with no error anywhere.
            #
            # `chains` is the escape hatch and it is the right one: a caller who
            # names the protomers has asserted what is being measured, which is
            # exactly rule 2b. Refusing without one is refusing to guess.
            if (
                homo["n_polymer_chains_in_file"] >= LARGE_ASSEMBLY_CHAINS
                and not want
            ):
                stage = "large_assembly"
                raise RuntimeError(
                    f"{pid}: {homo['n_polymer_chains_in_file']} polymer chains "
                    f"(>= {LARGE_ASSEMBLY_CHAINS}) and no `chains` restriction. "
                    "Whole-assembly pocket detection on an assembly this size "
                    "returns hundreds of pockets, nearly all of them shallow "
                    "surface features and inter-protomer crevices: 1JH5, a "
                    "60-mer of ONE protein, returned 378 pockets and a selected "
                    "pocket 60.28 A from the protein centre, with no error "
                    "anywhere. REFUSING rather than returning that. 1OQE, 1OQD "
                    "and 4V46 already refused at exactly this shape, but only "
                    "because they exceeded the 62 single-character chain IDs — "
                    "an implementation limit, not a boundary that means "
                    "anything. Pass `chains` naming the protomers that carry "
                    "the site, e.g. {\"1JH5\": [\"A\",\"B\",\"C\"]}."
                )
            homo["large_assembly_warning"] = (
                None if homo["n_polymer_chains_in_file"] < LARGE_ASSEMBLY_CHAINS
                else (
                    f"this file contains {homo['n_polymer_chains_in_file']} "
                    f"polymer chains ({len(used_chains)} scored). Whole-assembly "
                    "pocket detection on an assembly this size returns hundreds "
                    "of pockets, most of them shallow surface features and "
                    "inter-protomer crevices far from anything a dossier is "
                    "asking about — 1JH5, a 60-mer, returned 378 pockets and a "
                    "selected pocket 60.28 A from the protein centre. This run "
                    "was allowed through ONLY because `chains` was passed, so "
                    "the caller has asserted which protomers carry the site; "
                    "without that restriction the entry is refused."
                )
            )
        except LigandSourceError:
            # DELIBERATELY NOT RECORDED AS A PER-STRUCTURE ERROR. Every
            # structure would carry the same message and the run would return a
            # full, well-formed, entirely holo-free payload — which is exactly
            # the silent failure this exception exists to prevent. A
            # misconfigured chemical-component source is a run-level fault and
            # must kill the run.
            raise
        except Exception as exc:  # noqa: BLE001
            results[pid] = {
                "error": f"{type(exc).__name__}: {exc}",
                "stage": stage,
                "tier": "none",
                "by_clustering": {},
            }
            continue

        if want:
            missing_res = [r for r in missing_res if r.split("/")[0] in want]
        druglike = [lig for lig in ligs if lig["druglike"]]
        cofactors = sorted({lig["comp_id"] for lig in ligs if lig["cofactor"]})

        # Ground-truth site, when a drug-like ligand is present.
        #
        # `ligand_codes` is an OVERRIDE, not a gate. The old `if ligand_codes:
        # ... elif druglike:` meant that supplying ANY code disabled
        # auto-derivation for every structure whose ligand was not in the list,
        # so a structure this module had already identified as holo
        # ("drug-like ligand VGY") was scored as though it were apo. That is not
        # a missing feature, it is a wrong answer that moves with the command
        # line: 7JRA went 0.000 -> 0.926 druggability and 306.9 -> 1542.9 A^3
        # purely on whether its code was passed.
        target_comp = None
        anchor_source = None
        if ligand_codes:
            target_comp = next(
                (
                    lig["comp_id"]
                    for lig in ligs
                    if lig["comp_id"] in ligand_codes
                ),
                None,
            )
            if target_comp:
                anchor_source = "caller"
        if target_comp is None and druglike:
            # This structure's own drug-like ligand. Already computed above and
            # already reported in tier_note; there is no reason it should not
            # anchor the site.
            target_comp = druglike[0]["comp_id"]
            anchor_source = "auto_derived"
        true_site, site_copy = (
            _ligand_site(st, target_comp, used_chains)
            if target_comp
            else ([], None)
        )
        protein_centroid = _centroid(_pdb_coords(prepped))

        # Hand the downstream stages exactly what fpocket was given.
        prepped_by_pid[pid] = prepped
        chains_by_pid[pid] = used_chains
        # The TARGET's chains, carried to the mdpocket stage so the common core
        # is an intersection over one protein rather than over whatever else the
        # entry contains. See `_mdpocket_ensemble`.
        tgt_chains_by_pid[pid] = list(tgt_chains)
        cif_by_pid[pid] = cif
        # The TARGET's sequence, not the assembly's longest chain. This feeds
        # the interface stage's target/partner split and the disorder fallback,
        # and both were wrong on every S1PR1 entry because G-beta-1 is longer
        # than the receptor.
        seq, seq_chain = _one_letter(st, tgt_chains)
        if seq:
            seqs_by_pid[pid] = seq
        target_chain_info = {
            "target_accession": target_accession,
            "target_accession_basis": accession_basis,
            "target_chains": tgt_chains,
            "target_chains_basis": tgt_basis,
            "target_chains_verified": tgt_info["verified"],
            "target_chains_note": tgt_info["reason"],
            "entry_declares_accessions": tgt_info["declared_accessions"],
            "chain_accession_status": acc_status,
            "chain_accessions": chain_acc,
            "non_target_chains_scored": [
                c for c in used_chains if c not in tgt_chains
            ],
            "sequence_chain_used": seq_chain,
            "_why": (
                "Everything that used to identify the target by chain LENGTH "
                "is wrong whenever a partner chain is longer, which on a "
                "GPCR-G-protein complex is always. Chains here are resolved by "
                "UniProt accession from the entry's own _struct_ref_seq."
            ),
        }
        if target_comp and true_site:
            # The ligand's own heavy atoms, for the mdpocket site transfer, plus
            # the chains that line it — a crystallographic second copy in the
            # same file must not join the superposition.
            lig_xyz = [
                (a.pos.x, a.pos.y, a.pos.z)
                for chain in st[0] for res in chain
                if res.het_flag == "H" and res.name == target_comp
                and (site_copy is None or
                     f"{chain.name}/{res.seqid.num}" == site_copy)
                for a in res
            ]
            keep_chains, chain_atom_counts = _ligand_contact_chains(
                st, target_comp, site_copy
            )
            holo_by_pid[pid] = {
                "comp_id": target_comp,
                "ligand_xyz": [list(x) for x in lig_xyz],
                # Atom counts, 15% of the top contributor — the same rule
                # cryptic_analysis uses to reject a chain that merely brushes
                # the ligand across a crystal contact. See
                # `_ligand_contact_chains` for why residue counts are wrong.
                "site_chains": keep_chains,
                "site_chain_atom_counts": chain_atom_counts,
                "ligand_chain": site_copy.split("/")[0] if site_copy else None,
                "ligand_copy": site_copy,
                "site_residues": true_site,
            }

        per_d = {}
        for d in D_VALUES:
            run = work / f"{pid}_D{d}"
            run.mkdir(parents=True, exist_ok=True)
            tgt = run / prepped.name
            tgt.write_text(prepped.read_text())
            out_dir = run / f"{tgt.stem}_out"
            # A warm Modal container reuses /tmp. fpocket overwrites what it
            # rewrites but leaves everything else, so a failed rerun would be
            # parsed as a successful one off the previous run's files.
            shutil.rmtree(out_dir, ignore_errors=True)
            proc = subprocess.run(  # noqa: S603
                ["fpocket", "-f", str(tgt), "-D", str(d)],  # noqa: S607
                check=False,
                capture_output=True,
            )
            pockets = _parse_pockets(out_dir)
            # `tgt`, not `prepped`: PRANK must be given the identical file
            # fpocket read, or its surface and fpocket's alpha spheres are
            # computed over different atoms.
            prank_ranks, prank_info = _prank_rescore(out_dir, run, tgt)
            for p in pockets:
                p["jaccard_vs_ligand_site"] = (
                    _jaccard(p.get("residues", []), true_site) if true_site else None
                )
                # fpocket's geometry with PRANK's ranking. Reported alongside
                # fpocket's own rank, never replacing it — a large gap between
                # the two is itself a finding about how much the ranking is
                # carrying.
                p["prank_rank"] = prank_ranks.get(p["rank"])
            # ---- ON-TARGET FILTER, BEFORE ANY SELECTION -------------------
            # The chain resolver was already right and selection ignored it.
            # See `_annotate_on_target` for what that cost. Every pocket is
            # annotated (so the fraction is readable for all of them) and only
            # the on-target ones are eligible to BE the site.
            _annotate_on_target(pockets, tgt_chains, tgt_info["verified"])
            _annotate_pocket_labels(
                pockets, tgt_chains, chain_acc,
                homo.get("identical_chains"), unp_sites, unp_offsets,
            )
            candidates = [p for p in pockets if p.get("on_target") is not False]
            off_target = [p for p in pockets if p.get("on_target") is False]
            if pockets and not candidates:
                # NOT "the best of a bad set". An entry in which no pocket sits
                # on the target contributes NOTHING rather than contributing a
                # partner's value — genuinely the case for BAFF 5Y9J, where none
                # of 22 pockets is on-target and the selected one was lined by
                # belimumab.
                best, basis = None, "no_on_target_pocket"
            elif true_site:
                best = max(
                    candidates,
                    key=lambda p: p.get("jaccard_vs_ligand_site") or 0.0,
                    default=None,
                )
                basis = "ligand_site_jaccard"
                if best and not (best.get("jaccard_vs_ligand_site") or 0.0):
                    # Nothing touched the real site at this D. Returning the
                    # arbitrary first pocket as "the site pocket" is worse than
                    # returning nothing — that is the false negative in rule 4.
                    best, basis = None, "no_pocket_overlapped_ligand_site"
            elif site_signature:
                # SAME-SITE TRACKING. Without this, an apo structure falls back
                # to "most druggable pocket anywhere", and pooling those across
                # an ensemble compares different pockets on different proteins'
                # surfaces — which is not a measurement of anything.
                #
                # Match by residue NUMBER, deliberately chain-agnostic: the
                # site of interest is often inter-subunit (TNF-alpha's axial
                # channel is lined by the same residues from all three chains),
                # and chain letters are not stable across depositions.
                for p in pockets:
                    nums = {r.split("/")[-1] for r in p.get("residues", [])}
                    p["signature_overlap"] = (
                        round(len(nums & site_signature) / len(site_signature), 3)
                        if site_signature
                        else None
                    )
                best = max(
                    candidates, key=lambda p: p.get("signature_overlap") or 0.0,
                    default=None,
                )
                basis = "site_signature_overlap"
                if best and not (best.get("signature_overlap") or 0.0):
                    best, basis = None, "no_pocket_matched_site_signature"
                elif signature_collapsed_by > 0:
                    # KEYED ON THE MEASURED COLLAPSE, NOT ON CHAIN COUNT.
                    # Chain-agnostic numbers cannot resolve a symmetric site —
                    # but only when the numbers actually collapse. The old test
                    # fired on donor chain count, so apo structures were flagged
                    # `site_signature_unreliable_homooligomer` while reporting
                    # `collapsed_by: 0`: 23 residues in, 23 distinct numbers
                    # out, nothing ambiguous about the match. An
                    # over-conservative flag gets ignored, and an ignored flag
                    # is how a real one gets missed.
                    #
                    # The real hazard is unchanged and still caught: 2AZ5's 22
                    # site residues collapse to 14 numbers across four identical
                    # chains (collapsed_by 8), and pooling those 10 of 12
                    # measurements anyway regenerates a fold_range of exactly
                    # 651.0 — the withdrawn claim reproducing itself from the
                    # identical defect. The guard is load-bearing where it fires.
                    basis = "site_signature_unreliable_homooligomer"
                elif signature_foreign_dropped and (
                    signature_foreign_dropped >= 0.33 * signature_n_residues_in
                    or len(site_signature) < 6
                ):
                    # THE SIGNATURE WAS NEVER THIS PROTEIN'S. 8QFZ: 13 residues
                    # in, 9 of them the bicyclic peptide's, 4 left. A 4-number
                    # signature matched chain-agnostically will hit something in
                    # every structure and that hit means nothing. Reported per
                    # structure, never pooled, exactly like the homo-oligomer
                    # case — and it counts a DIFFERENT failure: `collapsed_by`
                    # counts numbers lost to IDENTICAL protomers, this counts
                    # numbers imported from a DIFFERENT polymer, and only the
                    # first was ever guarded.
                    basis = "site_signature_unreliable_foreign_polymer"
            elif signature_foreign_dropped:
                # Every residue of the donor's contact shell belonged to another
                # polymer, so there is no signature left at all. Contributing
                # nothing is the honest outcome; the alternative branch below
                # would contribute "the most druggable pocket anywhere".
                best, basis = None, "site_signature_unreliable_foreign_polymer"
            else:
                # "The most druggable pocket ANYWHERE" — now at least anywhere
                # ON THE TARGET. This branch is the one that produced every
                # measured failure, because it is the branch a pure-apo entry
                # with no site signature falls into, and it ranked an antibody's
                # interior against the target's surface on equal terms.
                best = max(
                    candidates,
                    key=lambda p: p.get("druggability_score") or 0.0,
                    default=None,
                )
                basis = "max_druggability_no_ligand_site"
            site_centroid = best.get("centroid") if best else None
            ranked = sorted(pockets, key=lambda p: p["rank"])
            # Alpha-sphere centres, for the interface stage. Kept as PATHS, not
            # as coordinates: a few hundred points per pocket per D per
            # structure would bloat the payload for a consumer that is a
            # classification label.
            #
            # MORE THAN THE SELECTED POCKET. Rule 2b asks for every detected
            # pocket to be classified against the interface, and classifying
            # exactly one made that unsatisfiable. Every pocket is not
            # affordable — enclosure casts 512 rays per probe point per chain —
            # so the top ranks plus the selected pocket are classified and the
            # count that was not is reported.
            cls_ranks = [p["rank"] for p in ranked[:MAX_POCKETS_CLASSIFIED]]
            if best and best["rank"] not in cls_ranks:
                cls_ranks.append(best["rank"])
            if cls_ranks:
                vert_by_pid.setdefault(pid, {})[str(d)] = {
                    "site_rank": best["rank"] if best else None,
                    "n_pockets": len(pockets),
                    "verts": {
                        r: out_dir / "pockets" / f"pocket{r}_vert.pqr"
                        for r in sorted(cls_ranks)
                    },
                    "apolar_by_rank": {
                        p["rank"]: p.get("apolar_lining_fraction")
                        for p in ranked
                    },
                }
            # ---- EVERY POCKET, not only the selected one -------------------
            # The module used to return one pocket per structure, which makes
            # the dossier's rule 2b ("classify every detected pocket")
            # unsatisfiable from this output: the IRAK4 run had to re-run
            # fpocket locally to see the other 133, and reproduced these counts
            # exactly, so the data existed here and was being discarded. Worse
            # than lost data — on TL1A, 2RE9 reported `n_pockets: 31` while
            # carrying only rank 28, so the agent could not tell whether the
            # axial cavity was ABSENT in that structure or merely UNSELECTED,
            # and could honestly report neither a persistence nor a zero.
            #
            # Payload is bounded by rank, and WHAT WAS DROPPED IS STATED. Silent
            # truncation reads as completeness, which is the same failure in a
            # new place.
            keep_ranks = {p["rank"] for p in ranked[:MAX_POCKETS_RETURNED]}
            if best:
                keep_ranks.add(best["rank"])
            returned = [dict(p) for p in ranked if p["rank"] in keep_ranks]
            for p in returned:
                p["is_site_pocket"] = bool(best and p["rank"] == best["rank"])
            omitted = [p for p in ranked if p["rank"] not in keep_ranks]
            anchored = [p for p in returned if p.get("anchor_labels")]
            per_d[str(d)] = {
                "n_pockets": len(pockets),
                # THE REPORTABLE FORM OF DRUGGABILITY, in one object, at the
                # level a consumer reads. All four parts already existed and
                # none of them were together: the fpocket rank was nested inside
                # `site_pocket`, the PRANK rank beside it, the count here, and
                # the PDB ID was only the enclosing dict's KEY. The dossier
                # template's `tractability.site_pocket_rank` therefore had no
                # single source to read, which is the same shape of gap that let
                # `ligand_site_jaccard` be computed, used and thrown away. The
                # value may travel INSIDE this object; it may not travel out of
                # it into a cross-structure comparison.
                "site_pocket_rank": {
                    "fpocket": best["rank"] if best else None,
                    "prank": best.get("prank_rank") if best else None,
                    "n_pockets": len(pockets),
                    "structure_pdb_id": pid,
                    "clustering_d": d,
                    "druggability_score": (
                        best.get("druggability_score") if best else None
                    ),
                    "_why": (
                        "'rank 1 of 30 in 6OIM' is the claim. fpocket rank and "
                        "PRANK rank are two WITHIN-STRUCTURE orderings on the "
                        "same footing (PRANK at n=70 ligand-anchored promotes "
                        "the true site in 79% and demotes in 1%, the one "
                        "demotion being 6OIM at D=1.6); report both, replace "
                        "neither. The druggability VALUE beside them is "
                        "normalised over this file's own pocket list and must "
                        "not be compared to another structure's."
                    ),
                },
                # ---- THE DISTRIBUTION IS THE PRIMARY OUTPUT ----------------
                # One compact row per returned pocket. Read this before
                # `site_pocket`: a single elected number can silently be a
                # cavity on MAX rather than MYC, on IL-11 receptor alpha rather
                # than IL-11, or inside tralokinumab rather than on IL-13 — and
                # all four of those really were calibration anchors. A table of
                # thirty rows carrying their own chains and accessions cannot
                # hide it, and reporting a distribution removes the
                # maximum-over-N selection bias by construction.
                "pocket_table": _pocket_table(returned),
                "anchor_summary": {
                    "labels_available": [
                        "ligand_site", "interface", "symmetry_axis",
                        "annotated_functional_site", "buried_core",
                    ],
                    "n_pockets_with_any_anchor": len(anchored),
                    "anchors_seen": sorted({
                        lab for p in anchored for lab in p["anchor_labels"]
                    }),
                    "site_hypothesis_basis": (
                        "not_established" if not anchored
                        else "external_anchor_labels_present"
                    ),
                    "uniprot_features_resolved": unp_sites_ok,
                    "n_uniprot_feature_positions": len(unp_sites),
                    "transferred_homolog_site": (
                        "NOT AVAILABLE IN THIS MODULE — it needs Foldseek, "
                        "which lives in structure-select / neighbour_precedent. "
                        "Its absence from `anchors_seen` is a statement about "
                        "this module, not about the protein."
                    ),
                    "interface_and_buried_core": (
                        "added by the interface stage, which is where the "
                        "partner epitope and the enclosure calculation exist. "
                        "Absent here when no partner_structures were supplied."
                    ),
                    "_why": (
                        "ANCHORING IS AN ANNOTATION AND SEVERAL CAN COEXIST. "
                        "A pocket may carry ligand_site AND interface AND "
                        "symmetry_axis, or none of them. `not_established` now "
                        "means 'no pocket carries an external label', which is "
                        "a true and useful statement about a protein — it used "
                        "to mean 'we fell back to whatever scored highest', "
                        "which is how every bad calibration anchor was born. "
                        "The open question worth answering next: where a target "
                        "has BOTH a ligand site and a receptor epitope (TNF and "
                        "IL-17A both do), do the interface- and symmetry-"
                        "anchored labels land on the same pocket as the "
                        "ligand-anchored one? If they agree, this axis works on "
                        "targets with no chemistry at all."
                    ),
                },
                "pockets": returned,
                "pockets_returned": len(returned),
                "pockets_omitted": len(omitted),
                "pockets_omitted_note": (
                    None if not omitted else
                    f"{len(omitted)} of {len(pockets)} detected pockets are not "
                    f"in `pockets` above. They are fpocket ranks "
                    f"{min(p['rank'] for p in omitted)}-"
                    f"{max(p['rank'] for p in omitted)}, i.e. everything below "
                    f"rank {MAX_POCKETS_RETURNED} that is not the selected site "
                    "pocket. THEY WERE DETECTED, NOT ABSENT — a pocket missing "
                    "from this list is not evidence that it does not exist. "
                    "The summary below bounds what is in them."
                ),
                "pockets_omitted_summary": (
                    None if not omitted else {
                        "n": len(omitted),
                        "rank_range": [min(p["rank"] for p in omitted),
                                       max(p["rank"] for p in omitted)],
                        "max_volume_a3": round(
                            max((p.get("volume") or 0.0) for p in omitted), 2),
                        "max_druggability_score": round(
                            max((p.get("druggability_score") or 0.0)
                                for p in omitted), 3),
                        "max_jaccard_vs_ligand_site": max(
                            (p.get("jaccard_vs_ligand_site") or 0.0)
                            for p in omitted),
                        "max_signature_overlap": max(
                            (p.get("signature_overlap") or 0.0)
                            for p in omitted),
                        "_why": (
                            "so a reader can see that nothing large, nothing "
                            "druggable and nothing overlapping the site was "
                            "hidden by the truncation, without having to trust "
                            "that claim."
                        ),
                    }
                ),
                "max_pockets_returned": MAX_POCKETS_RETURNED,
                # KEPT AS AN ANNOTATION, NOT A GATE. `site_pocket` is one row
                # of `pocket_table` singled out, and it is retained so existing
                # consumers keep reading the field they always did. It is no
                # longer the answer: read the table, and read
                # `anchor_summary.site_hypothesis_basis` before treating this
                # pocket as the site.
                "site_pocket": best,
                "site_pocket_selected_by": basis,
                "_site_pocket_note": (
                    "ONE ROW OF pocket_table, NOT THE ANSWER. Electing a single "
                    "pocket is a maximum over N draws and it is where four bad "
                    "calibration anchors came from (MYC's on MAX, IL-11's on "
                    "the receptor, IL-13's inside tralokinumab, CD20's on a "
                    "cholesterol site). Quote the distribution; quote this only "
                    "with its on_target_fraction and its anchors beside it."
                ),
                # ---- IS THE SELECTED POCKET EVEN ON THE TARGET? ------------
                # The question that had no field. `on_target_selection` is the
                # answer for the pocket that was chosen; the census beside it is
                # the answer for the entry.
                "on_target_selection": {
                    "target_chains": list(tgt_chains),
                    "target_chains_verified": tgt_info["verified"],
                    "min_on_target_fraction": POCKET_MIN_ON_TARGET_FRACTION,
                    "_threshold_status": "PROPOSED, NOT CALIBRATED",
                    "site_pocket_on_target_fraction": (
                        best.get("on_target_residue_fraction") if best else None
                    ),
                    "site_pocket_off_target_chains": (
                        best.get("off_target_lining_chains") if best else None
                    ),
                    "n_pockets_on_target": len(candidates),
                    "n_pockets_off_target": len(off_target),
                    # The stricter count, reported so a consumer can apply a
                    # stricter policy than this module's without re-running it.
                    # BAFF 5Y9J has ZERO of 22 fully on-target and exactly one
                    # above the 0.5 floor (at 0.667, druggability 0.000): under
                    # `>= 0.5` that entry contributes 118.8 A^3, under `== 1.0`
                    # it contributes nothing. Both are defensible and this
                    # module does not have the evidence to pick, so it reports
                    # both counts and applies the looser one.
                    "n_pockets_fully_on_target": sum(
                        1 for p in pockets
                        if (p.get("on_target_residue_fraction") or 0.0) >= 1.0
                    ),
                    "off_target_ranks": [p["rank"] for p in off_target][:40],
                    "off_target_max_volume_a3": (
                        round(max((p.get("volume") or 0.0) for p in off_target), 2)
                        if off_target else None
                    ),
                    "_why": (
                        "The chain resolver worked and selection ignored it. "
                        "Selected pockets that were actually on the target: "
                        "IL-13 1 of 9 (the rest inside the Fabs of tralokinumab "
                        "and lebrikizumab, and on the receptor chain of 3LB6), "
                        "BAFF 2 of 5 (5Y9J lined by belimumab), CD20 4 of 7 "
                        "(6Y90/6Y97 lined by rituximab's Fab). Filtering on it "
                        "moves IL-13's median volume 312.3 -> 106.8 A^3, BAFF "
                        "258.3 -> 177.4, CD20 281.0 -> 242.3 — a verdict-"
                        "relevant inversion, previously unflagged anywhere in "
                        "the payload."
                    ),
                    "_unverified_note": (
                        None if tgt_info["verified"] else
                        "THE CHAIN SET IS UNVERIFIED, so no pocket was excluded "
                        "and on_target is null throughout. Either no "
                        "uniprot_accession was supplied or this entry's "
                        "_struct_ref could not be read. The selected pocket has "
                        "NOT been shown to be on the target; read "
                        "target_chains_basis before quoting its volume."
                    ),
                    "_off_target_note": (
                        None if not off_target else
                        f"{len(off_target)} of {len(pockets)} detected pockets "
                        f"are lined <{POCKET_MIN_ON_TARGET_FRACTION:.0%} by the "
                        "target's own chains and were EXCLUDED from selection. "
                        "They are still returned in `pockets` with their "
                        "on_target_residue_fraction, because a cavity inside a "
                        "partner is a real cavity — it is just not this "
                        "target's site and must never be quoted as its volume."
                    ),
                    "_no_candidate_note": (
                        None if basis != "no_on_target_pocket" else
                        f"NO POCKET IN THIS ENTRY IS ON THE TARGET. All "
                        f"{len(pockets)} detected pockets are lined mostly by "
                        "chains that are not the target's, so this entry "
                        "contributes NO volume and NO druggability rather than "
                        "contributing a partner's. This is a real outcome, not "
                        "a failure: BAFF 5Y9J has no fully on-target pocket "
                        "among its 22."
                    ),
                },
                # The number behind the basis in defect 7's sense: whether the
                # residue-number signature could identify a site AT ALL, stated
                # per measurement rather than only once at ensemble level.
                "site_signature_collapsed_by": signature_collapsed_by,
                # THE VALUE BEHIND THE BASIS, promoted out of `site_pocket`.
                # It was being computed, used to make the selection, and then
                # left buried one level down inside the pocket dict, where the
                # dossier's `tractability.ligand_site_jaccard` never found it —
                # the measurement was made and thrown away. A basis without its
                # number is not checkable: "selected by ligand_site_jaccard" is
                # equally true at 0.74 and at 0.02, and those are not the same
                # claim about whether the pocket is the site.
                "site_pocket_ligand_site_jaccard": (
                    best.get("jaccard_vs_ligand_site") if best else None
                ),
                "site_pocket_signature_overlap": (
                    best.get("signature_overlap") if best else None
                ),
                "merge_suspected": bool(
                    best and best.get("volume", 0) > 1000
                ),
                # The centroid control. An overlap fraction cannot tell you that
                # two pockets sharing residue numbers are 12 A apart; this can.
                "site_pocket_centroid": site_centroid,
                # Frame-independent companion to the raw centroid: the distance
                # from this structure's own protein centre. Two deposited
                # entries are not in a common frame, so a raw centroid-to-
                # centroid distance across structures also contains their rigid
                # -body offset; this radius does not.
                "site_pocket_radius_from_protein_center_a": _distance(
                    site_centroid, protein_centroid
                ),
                **prank_info,
            }
            # ---- DE-DUPLICATE: the annotation lives in `pocket_table` -------
            # A dry run hit the consumer's 180,000-char cap, dropped 98 pocket
            # objects, was STILL over and truncated mid-string — producing
            # invalid JSON and deleting the trailing explanation first. Carrying
            # the same per-pocket annotation in both the compact table and the
            # verbose list costs ~160 chars x 30 pockets x every clustering
            # value x every structure for nothing. The selected pocket keeps its
            # copy, because `site_pocket` is read on its own.
            for _p in returned:
                if best is not None and _p is best:
                    continue
                for _k in (
                    "anchor_labels", "anchor_detail", "lining_chains",
                    "lining_chain_accessions", "n_on_target_lining_residues",
                    "off_target_lining_chains",
                ):
                    _p.pop(_k, None)
            if not pockets:
                per_d[str(d)]["fpocket_failed"] = {
                    "returncode": proc.returncode,
                    "stderr": proc.stderr.decode(errors="replace")[-500:],
                }

        src_marker = work / f"{pid}.source"
        results[pid] = {
            "structure_source": (
                src_marker.read_text() if src_marker.exists() else "unknown"
            ),
            # experimental vs predicted, and NEVER absent. `structure_source`
            # above is file provenance (assembly1 / asu / esmfold_predicted);
            # THIS is the origin axis a reader keys on to discount a model
            # pocket. See CLAUDE.md rule 4c.
            "structure_origin": (
                "esmfold_predicted" if pid in predicted_ids else "experimental"
            ),
            "chains_used": used_chains,
            # WHAT CAME OUT OF THE fpocket INPUT AND WHY. Rule 4 says strip
            # every ligand before scoring; `het_flag == 'A'` kept every polymer,
            # so a peptide, nanobody or designed mini-binder was never stripped
            # and lined the pocket it was being scored in.
            "chains_dropped_as_polymer_ligand": dropped_chains,
            "polymer_ligand_control": polymer_control,
            "missing_residues": missing_res,
            "ligands": ligs,
            "cofactors_present": cofactors,
            # THREE TIERS, NOT TWO. `undetermined` is not `apo`. When a
            # component could not be classified — because its record could not
            # be retrieved, not because the chemistry was unrecognised — this
            # entry's state is unknown, and calling it apo would be the same
            # class of error as reporting a credential failure as "no data".
            "tier": (
                "holo" if druglike
                else "apo" if holo_call.get("determined", True)
                else "undetermined"
            ),
            "tier_note": (
                f"drug-like ligand {target_comp or druglike[0]['comp_id']}"
                if druglike else
                "HOLO/APO UNDETERMINED: the chemical-component lookup failed "
                f"for {', '.join(holo_call.get('undetermined') or [])}. This is "
                "a lookup failure, not an absence of ligand, and this entry "
                "must not be counted as apo."
                if not holo_call.get("determined", True) else
                "no drug-like ligand (classified by chemistry, "
                "ligand_filter.classify_record — not by a size floor or a "
                "comp_id list)"
                + (f"; cofactors present: {', '.join(cofactors)}"
                   if cofactors else "")
                + (f"; other components: "
                   f"{', '.join(k for k in (holo_call.get('by_verdict') or {}) if k not in ('druglike', 'cofactor'))}"
                   if holo_call.get("by_verdict") else "")
            ),
            # The full entry-level call with every rejected component and the
            # reason it was rejected, so "apo" is checkable rather than asserted.
            "holo_call": holo_call,
            "ligand_site_residues": true_site,
            "ligand_site_copy": site_copy,
            # What actually anchored the site, and whether the caller chose it.
            # Without these two the same structures at the same clustering can
            # return different answers and the output gives no way to tell why.
            "site_anchor_ligand": target_comp,
            "site_anchor_ligand_source": anchor_source,
            "site_anchor_available_druglike": [
                lig["comp_id"] for lig in druglike
            ],
            # Measured over the TARGET's chains only. A partner's homodimer is
            # not evidence that the target's site signature is ambiguous.
            "homo_oligomer": homo,
            "target_chains": target_chain_info,
            "protein_centroid": protein_centroid,
            "site_pocket_centroids": {
                k: v["site_pocket_centroid"] for k, v in per_d.items()
            },
            # Within ONE structure every clustering value shares a coordinate
            # frame, so this distance is exact and needs no superposition — the
            # frame-free half of the centroid control. Measured on apo TNF-alpha
            # 1TNF: the pocket called "the site" at D 1.6 and the one called
            # "the site" at D 2.4 are ~12 A apart, in the same structure. An
            # overlap fraction reported both as the same site.
            "site_pocket_centroid_spread_across_clustering_a": max(
                (
                    _distance(a, b) or 0.0
                    for a, b in _pairs(
                        [
                            v["site_pocket_centroid"]
                            for v in per_d.values()
                            if v["site_pocket_centroid"]
                        ]
                    )
                ),
                default=None,
            ),
            "by_clustering": per_d,
        }
        if renamed:
            # Only when an mmCIF chain name would not fit a PDB column. Present
            # so a caller comparing against the deposited entry can see that
            # the chain IDs in every residue list above are ours, not RCSB's.
            results[pid]["chain_renamed_from_cif"] = renamed

    # ---- holo reference from INSIDE the ensemble, if none came from outside
    if donor is None:
        pid0 = (
            mdpocket_site_donor
            if mdpocket_site_donor in holo_by_pid
            else next((p for p in pdb_ids if p in holo_by_pid), None)
        )
        if pid0:
            donor = {
                "pdb_id": pid0, "cif": cif_by_pid[pid0],
                "prepped": prepped_by_pid[pid0], "in_ensemble": True,
                **holo_by_pid[pid0],
            }

    # =======================================================================
    # STAGE 2 — DISORDER. Non-fatal; never 0.0 on failure.
    # =======================================================================
    disorder_out: dict = {
        "disorder_status": "not_run",
        "disorder_reason": "run_disorder=False",
        "disorder_fraction": None,
    }
    if run_disorder:
        # THE FULL-LENGTH PATH IS THE DEFAULT WHEREVER AN ACCESSION EXISTS, and
        # one usually does exist without the caller supplying it: every
        # deposited entry declares its UniProt accession in `_struct_ref`. The
        # old code only looked at the argument, so an omitted argument silently
        # switched the measurement onto the crystallised construct — a
        # different molecule — and IRAK4 came back 0.0 over 284 residues
        # against a true 0.1413 over 460.
        acc, acc_src = target_accession, accession_basis
        acc_candidates = accession_by_pid
        seq = seq_src = None
        if acc:
            seq_src = f"uniprot:{acc}"
        else:
            for pid in pdb_ids:
                if pid in seqs_by_pid:
                    seq, seq_src = seqs_by_pid[pid], f"structure:{pid}"
                    break
        disorder_out = _disorder_block(acc, seq, seq_src, acc_src)
        disorder_out["accession_candidates_by_structure"] = acc_candidates
        if len({a for v in acc_candidates.values() for a in v}) > 1:
            # A complex, a chimera or a fusion construct. Report it rather than
            # letting the first entry's accession stand in for the target.
            disorder_out["_accession_ambiguity"] = (
                "the ensemble's entries declare more than one UniProt "
                "accession; the first was used. If the target is not that one, "
                "pass uniprot_accession explicitly."
            )

    # =======================================================================
    # STAGE 3 — CRYPTIC MECHANISM. Needs BOTH an apo and a holo in the run.
    # =======================================================================
    cryptic_out: dict = {
        "cryptic_status": "not_run",
        "cryptic_reason": "run_cryptic=False",
    }
    if run_cryptic:
        apo_ids = [
            p for p in pdb_ids
            if p in prepped_by_pid and results.get(p, {}).get("tier") == "apo"
        ]
        if donor is None or not apo_ids:
            cryptic_out["cryptic_reason"] = (
                f"needs both an apo and a holo structure; got "
                f"{0 if donor is None else 1} holo reference and "
                f"{len(apo_ids)} apo. This comparison is a PAIRWISE measurement "
                "— there is nothing to superpose against with only one state. "
                "A holo entry outside the ensemble can be supplied as "
                "mdpocket_site_donor."
                + (f" {donor_error}" if donor_error else "")
            )
        else:
            holo_pid = donor["pdb_id"]
            info = donor
            per_apo: dict[str, dict] = {}
            for apo_pid in apo_ids:
                per_apo[apo_pid] = _cryptic_block(
                    cif_by_pid[apo_pid], donor["cif"], info["comp_id"],
                    apo_pid, holo_pid,
                    chains_by_pid.get(apo_pid), info.get("ligand_chain"),
                )
                # Per-target as well as per-pair: a reader looking at one apo
                # structure must not have to find the pairwise block.
                #
                # WHICH IS EXACTLY WHY THE QUARANTINE HAS TO REACH HERE TOO.
                # Being the block a reader lands on first makes this the block
                # most likely to be QUOTED, and it was handing out a rejected
                # pair's numbers with a potency prior attached. Measured on
                # NLRP3: `cryptic_status: "failed"` sitting beside
                # `mechanism: loop_or_backbone_motion`, `is_cryptic: true`,
                # `max_backbone_ca_displacement_a: 21.13` and
                # `cryptic_potency_prior: {expected_ceiling: "nanomolar"}` —
                # every one of them measured in a 16.5 A misfit frame. The
                # aggregate below had already been fixed against precisely this;
                # this is the same bug one nesting level down. Whatever drops a
                # structure from the call drops it from every number derived
                # from the fit, at every level that reprints them.
                _blk = per_apo[apo_pid]
                _fit_ok = _blk.get("cryptic_status") == "ok"
                _derived = (
                    "mechanism", "is_cryptic",
                    "max_backbone_ca_displacement_a", "clash_attribution",
                    "cryptic_potency_prior",
                )
                results[apo_pid]["cryptic"] = {
                    k: (_blk.get(k) if (_fit_ok or k not in _derived) else None)
                    for k in (
                        "cryptic_status", "cryptic_reason", "mechanism",
                        "is_cryptic", "max_backbone_ca_displacement_a",
                        "clash_attribution", "cryptic_potency_prior",
                        "holo_pdb_id", "ligand_comp_id",
                    )
                }
                results[apo_pid]["cryptic"]["core_ca_rmsd_a"] = (
                    (_blk.get("superposition_gate") or {}).get("core_ca_rmsd_a")
                )
                results[apo_pid]["cryptic"]["_quarantined_keys"] = (
                    [] if _fit_ok else list(_derived)
                )
                results[apo_pid]["cryptic"]["_quarantine_note"] = (
                    None if _fit_ok else
                    "cryptic_status is not ok, so every key derived from the "
                    f"fit is null here: {', '.join(_derived)}. They are not "
                    "MISSING — they were computed and refused, and a null is "
                    "the only honest form for a number measured in a frame "
                    "this module has just rejected. The refused values are in "
                    f"cryptic.per_apo_structure.{apo_pid} for diagnosis; they "
                    "are not measurements of this structure and must not be "
                    "quoted as any. Read cryptic_reason for why."
                )
            ok_items = [
                (k, v) for k, v in per_apo.items()
                if v["cryptic_status"] == "ok"
            ]
            rejected = [
                {"pdb_id": k, "reason": v.get("cryptic_reason"),
                 "core_ca_rmsd_a": (v.get("superposition_gate") or {}).get(
                     "core_ca_rmsd_a")}
                for k, v in per_apo.items() if v["cryptic_status"] != "ok"
            ]
            # EVERY DERIVED NUMBER COMES FROM ONE STRUCTURE, and it is named.
            # The old block took `mechanism` and `is_cryptic` from ok[0] but
            # `max_backbone_ca_displacement_a` from a max over all of ok, so on a
            # disagreeing ensemble the label and the displacement described
            # different structures. Measured on NLRP3: the block reported
            # is_cryptic false / mechanism none (from 7ZGU) beside
            # max_backbone_ca_displacement_a 21.6 (from the rejected 8SWF).
            # Those cannot both be true, and whatever drops a structure from the
            # call must drop it from every statistic derived from it.
            #
            # The representative is the apo entry with the BEST superposition,
            # because that is the comparison least likely to be measuring the
            # frame rather than the site.
            rep_pid, rep = (
                min(
                    ok_items,
                    key=lambda kv: (
                        (kv[1].get("superposition_gate") or {}).get(
                            "core_ca_rmsd_a") or 0.0
                    ),
                )
                if ok_items else (None, {})
            )
            mechs = sorted({v.get("mechanism") for _k, v in ok_items
                            if v.get("mechanism")})
            cryptic_out = {
                "cryptic_status": "ok" if ok_items else "failed",
                "cryptic_reason": (
                    None if ok_items else
                    "; ".join(f"{k}: {v['cryptic_reason']}"
                              for k, v in per_apo.items())
                ),
                "holo_pdb_id": holo_pid,
                "holo_in_ensemble": donor["in_ensemble"],
                "ligand_comp_id": info["comp_id"],
                "per_apo_structure": per_apo,
                # The ensemble-level call: the dossier asks one question about
                # the target, and every field below describes ONE structure.
                "representative_apo_pdb_id": rep_pid,
                "representative_selected_by": (
                    "lowest core C-alpha RMSD among apo entries that passed the "
                    "superposition gate; all headline numbers below come from "
                    "this one structure so they cannot contradict each other"
                ),
                "representative_core_ca_rmsd_a": (
                    (rep.get("superposition_gate") or {}).get("core_ca_rmsd_a")
                ),
                "mechanism": rep.get("mechanism") if ok_items else None,
                "is_cryptic": rep.get("is_cryptic") if ok_items else None,
                "max_backbone_ca_displacement_a": (
                    rep.get("max_backbone_ca_displacement_a")
                    if ok_items else None
                ),
                "cryptic_potency_prior": (
                    rep.get("cryptic_potency_prior") if ok_items else None
                ),
                "displacement_by_apo_structure": {
                    k: v.get("max_backbone_ca_displacement_a")
                    for k, v in ok_items
                },
                "apo_structures_rejected": rejected,
                "_rejected_note": (
                    None if not rejected else
                    "These apo entries did not superpose onto the holo "
                    "reference and contribute NOTHING to any field above — not "
                    "the mechanism, not the displacement, not the census. A "
                    "rejected structure's displacement used to survive into the "
                    "aggregate while its mechanism did not."
                ),
                "mechanisms_across_apo": mechs,
                "mechanisms_agree": (len(mechs) <= 1) if ok_items else None,
                "_disagreement_note": (
                    None if len(mechs) <= 1 else
                    f"the apo entries disagree ({', '.join(mechs)}). The "
                    "headline mechanism is the best-superposed one and the rest "
                    "are in per_apo_structure. Do not average them, and do not "
                    "quote the headline without saying which structure it is."
                ),
                # DERIVED FROM THE INPUT, NOT HARDCODED. This note used to say
                # "With one apo entry this cannot be applied" unconditionally.
                # On a run with TWO apo entries (IRAK4: 2OIB and 2O8Y) it said
                # so anyway, contradicting the very census printed beside it.
                # A caveat that does not track its own data is worse than no
                # caveat, because it looks like it was checked.
                "n_apo_examined": len(apo_ids),
                "apo_examined": list(apo_ids),
                "n_apo_ok": len(ok_items),
                "_vajda_note": (
                    "Vajda's stringent definition requires the pocket to be "
                    "absent in ALL or nearly all unbound structures. "
                    + (
                        f"This run examined {len(apo_ids)} apo entr"
                        f"{'y' if len(apo_ids) == 1 else 'ies'} "
                        f"({', '.join(apo_ids)})"
                    )
                    + (
                        ", which is one structure — the definition cannot be "
                        "applied and the label below is a PAIRWISE RMSD result "
                        "only, not a cryptic call in Vajda's sense."
                        if len(apo_ids) < 2 else
                        f", of which {len(ok_items)} produced an interpretable "
                        "comparison. The definition can be applied to the "
                        "extent that this ensemble represents the unbound "
                        "states: read n_apo_examined as the denominator and "
                        "mechanisms_across_apo for whether they agree. A site "
                        "absent from some but not nearly all of them is "
                        "low-scoring, not cryptic."
                    )
                    + " Note the TNF-alpha case: ~1.6 A displacement (this "
                    "default protocol measures ~1.55-1.58 A; 1.62 A is the "
                    "hand-calibration figure and is NOT what this run "
                    "produces), site recovered in all five apo structures once "
                    "the third subunit is removed — occluded, NOT cryptic."
                ),
            }

    # =======================================================================
    # STAGE 4 — INTERFACE. Needs a structure containing the binding partner.
    # =======================================================================
    interface_out: dict = {
        "interface_status": "not_run",
        "interface_reason": (
            "no partner_structures supplied; without a complex containing the "
            "partner, orthosteric/allosteric is an assumption and this module "
            "refuses to assert it"
        ),
        "classification": "no_partner_structure",
    }
    if partner_structures:
        partner_cifs: dict[str, Path] = {}
        fetch_errors: dict[str, str] = {}
        for pid in partner_structures:
            try:
                partner_cifs[pid] = _fetch(pid, work)
            except Exception as exc:  # noqa: BLE001
                fetch_errors[pid] = f"{type(exc).__name__}: {exc}"
        interface_out = _interface_block(
            partner_cifs, [seqs_by_pid[p] for p in pdb_ids if p in seqs_by_pid],
            target_accession, work,
        )
        if fetch_errors:
            interface_out["fetch_errors"] = fetch_errors
        interface_out.setdefault("classification", "no_partner_structure")
        epitope = interface_out.pop("_epitope", None)
        partner_resnames = interface_out.pop("_partner_resnames", None)
        per_struct: dict[str, dict] = {}
        per_struct_by_rank: dict[str, dict] = {}
        if epitope:
            for pid, byd in vert_by_pid.items():
                for dkey, spec in byd.items():
                    site_rank = spec.get("site_rank")
                    apolar = spec.get("apolar_by_rank") or {}
                    by_rank: dict[str, dict] = {}
                    for rank, vert in sorted(spec.get("verts", {}).items()):
                        try:
                            cls = _classify_site_pocket(
                                prepped_by_pid[pid], chains_by_pid.get(pid, []),
                                vert, epitope, partner_resnames,
                                apolar.get(rank),
                            )
                        except Exception as exc:  # noqa: BLE001
                            cls = {"error": f"{type(exc).__name__}: {exc}"}
                        cls["fpocket_rank"] = rank
                        cls["is_site_pocket"] = rank == site_rank
                        by_rank[str(rank)] = cls
                    # The selected site pocket stays the headline, so an
                    # existing consumer reads the same field it always did.
                    cls = by_rank.get(str(site_rank)) or (
                        next(iter(by_rank.values())) if by_rank else {}
                    )
                    per_struct.setdefault(pid, {})[dkey] = cls
                    # Kept OUT of per_struct: the consensus aggregation below
                    # iterates per_struct[pid].values() and expects every entry
                    # to be one classification.
                    per_struct_by_rank.setdefault(pid, {})[dkey] = by_rank
                    # ---- FOLD THE LATE LABELS ONTO THE POCKET TABLE --------
                    # `interface` and `buried_core` cannot be computed in the
                    # fpocket loop: one needs the partner epitope and the other
                    # needs the enclosure ray-cast, and both exist only here.
                    # They are written back so a reader has ONE per-pocket
                    # record carrying every external label that applies, rather
                    # than having to join two blocks by rank to find out
                    # whether a pocket is anchored on anything at all.
                    # `thresholds` is a block of constants that
                    # `interface_analysis.classify_pocket` returns on EVERY
                    # pocket. Hoisted once, popped everywhere else — same
                    # reasoning as the buried-core cut points.
                    for _c in by_rank.values():
                        th = _c.pop("thresholds", None)
                        if th and "_classification_thresholds" not in interface_out:
                            interface_out["_classification_thresholds"] = th
                    blk = (results[pid].get("by_clustering") or {}).get(dkey) or {}
                    for row in blk.get("pocket_table") or []:
                        c = by_rank.get(str(row.get("rank")))
                        if not c:
                            continue
                        labs = list(row.get("anchors") or [])
                        det = dict(row.get("anchor_detail") or {})
                        ovl = c.get("overlap_fraction")
                        # An overlap computed across an illegal seqid match is
                        # not an anchor. See `_numbering_agreement`.
                        if ovl and c.get(
                            "overlap_unreliable_numbering_mismatch"
                        ) is not True:
                            if "interface" not in labs:
                                labs.append("interface")
                            det["interface_overlap"] = ovl
                        if c.get("buried_core_suspected"):
                            if "buried_core" not in labs:
                                labs.append("buried_core")
                            det["buried_core"] = True
                        row["anchors"] = labs
                        row["anchor_detail"] = det or None
                    summ = blk.get("anchor_summary")
                    if summ is not None:
                        anch = [
                            r for r in (blk.get("pocket_table") or [])
                            if r.get("anchors")
                        ]
                        summ["n_pockets_with_any_anchor"] = len(anch)
                        summ["anchors_seen"] = sorted({
                            lab for r in anch for lab in r["anchors"]
                        })
                        summ["site_hypothesis_basis"] = (
                            "not_established" if not anch
                            else "external_anchor_labels_present"
                        )
                    results[pid].setdefault("pocket_vs_interface", {})[dkey] = {
                        "classification": cls.get(
                            "classification", "no_partner_structure"),
                        "pocket_interface_overlap": cls.get("overlap_fraction"),
                        "enclosure": cls.get("enclosure"),
                        "subunit_enclosure_gain": cls.get(
                            "subunit_enclosure_gain"),
                        "min_distance_to_interface_a": cls.get(
                            "min_distance_to_interface_a"),
                        "adjacent_to_interface": cls.get("adjacent_to_interface"),
                        "also_overlaps_interface": cls.get(
                            "also_overlaps_interface"),
                        "partner_pdb_id": interface_out.get("partner_pdb_id"),
                        "error": cls.get("error"),
                        # Is the seqid match between these two entries legal at
                        # all? Without this the overlap is silently wrong
                        # wherever the numbering conventions differ.
                        "numbering_check": cls.get("numbering_check"),
                        "overlap_unreliable_numbering_mismatch": cls.get(
                            "overlap_unreliable_numbering_mismatch"),
                        # Is this a binding site or the hydrophobic core of a
                        # domain? Geometry only; flag, never filter.
                        "buried_core_suspected": cls.get("buried_core_suspected"),
                        "buried_core_reason": cls.get("buried_core_reason"),
                        # THE CLASSIFICATION IS ONLY AS GOOD AS THE POCKET IT
                        # WAS HANDED. This classifies whichever pocket
                        # `site_pocket_selected_by` chose, so on an apo
                        # homo-oligomer it is classifying the pocket the
                        # residue-number matcher landed on — which on apo
                        # TNF-alpha is 7.7 A off-site. Echoed here so the label
                        # can never be read without its basis.
                        "classified_pocket_selected_by": (
                            results[pid]["by_clustering"][dkey]
                            ["site_pocket_selected_by"]
                        ),
                        # Rule 2b asks for EVERY pocket. This is every pocket
                        # that could be afforded, with the shortfall stated.
                        "by_fpocket_rank": by_rank,
                        "n_pockets_classified": len(by_rank),
                        "n_pockets_detected": spec.get("n_pockets"),
                        "n_pockets_not_classified": max(
                            0, (spec.get("n_pockets") or 0) - len(by_rank)),
                        "not_classified_note": (
                            None
                            if (spec.get("n_pockets") or 0) <= len(by_rank)
                            else (
                                f"{(spec.get('n_pockets') or 0) - len(by_rank)} "
                                "detected pockets were not classified: enclosure "
                                "casts 512 rays per probe point per chain and "
                                f"only the top {MAX_POCKETS_CLASSIFIED} ranks "
                                "plus the selected site pocket are affordable. "
                                "They are still listed with their residues in "
                                "by_clustering.<D>.pockets, so the "
                                "interface-overlap half of rule 2b can be "
                                "computed from interface_residues without "
                                "re-running fpocket."
                            )
                        ),
                    }
            interface_out["per_structure"] = per_struct
            # NOT `per_struct_by_rank` — that is the SAME object already
            # serialised inside every structure, and emitting it twice was the
            # single largest line item in the payload: 156,009 characters of a
            # 529,048-character IL-13 run, against a consumer cap of 180,000.
            # A duplicate is not a second measurement.
            interface_out["per_structure_by_fpocket_rank"] = (
                "MOVED, NOT REMOVED: read "
                "structures.<PDB>.pocket_vs_interface.<D>.by_fpocket_rank, "
                "which is the identical content. It was duplicated here and "
                "cost ~156 kB on a two-structure run."
            )
            interface_out["_buried_core_thresholds"] = BURIED_CORE_THRESHOLDS_NOTE
            interface_out["_numbering_offset_rule"] = (
                "A constant numbering offset between two depositions of one "
                "protein is the normal case, not an anomaly — TL1A carries "
                "three at once (0, +67, +71) and IL-17A carries +23. Recovering "
                "it is a vote over residue names and costs one pass. It is "
                "APPLIED only when it converts an illegal comparison into a "
                "legal one over at least "
                f"{NUMBERING_MIN_COMPARED_FOR_OFFSET} shared positions AND "
                "strictly increases the number of name-agreeing positions, so "
                "an entry already on the partner's numbering is never shifted. "
                "Per-pocket values are in "
                "pocket_vs_interface.<D>.numbering_offset_to_partner."
            )
            # AGGREGATE, NEVER FIRST-WINS. Two symmetry copies of one ligand in
            # one structure can land either side of the overlap boundary:
            # measured on 8DYG ligand U5Q, copy A classified
            # allosteric_candidate at overlap 0.22 and copy B
            # orthosteric_candidate at 0.36, both flagged borderline. A caller
            # that takes whichever came first is tossing a coin between two
            # different mechanistic claims. So the label a caller may quote is
            # the consensus over every classification made, and a disagreement
            # is reported AS a disagreement rather than resolved.
            # AND A CLASSIFICATION BUILT ON AN ILLEGAL SEQID MATCH IS NOT PART
            # OF THE CONSENSUS. `numbering_agrees: false` used to flag the entry
            # and stop there, so the label derived from the flagged overlap
            # travelled up here unmarked and became the field SKILL.md tells
            # callers to quote. Measured on TL1A: `2RE9 -> allosteric_candidate`
            # came ENTIRELY from a 0.227 overlap on A:THR34/PRO35/THR36 —
            # VAL/VAL/ARG in the partner — and `2O0O -> mixed` from the
            # A:HIS118-vs-THR118 artifact. Recoverable offsets are applied
            # upstream in `_classify_site_pocket`; what reaches here still
            # flagged is genuinely uninterpretable and is excluded, listed, and
            # named — never silently dropped and never silently counted.
            def _usable(byd: dict) -> tuple[dict, dict]:
                keep = {
                    k: c for k, c in byd.items()
                    if c.get("overlap_unreliable_numbering_mismatch") is not True
                }
                return keep, {k: c for k, c in byd.items() if k not in keep}

            per_structure_label: dict[str, dict] = {}
            for pid_, byd in per_struct.items():
                keep, skipped = _usable(byd)
                labs = sorted({
                    c.get("classification") for c in keep.values()
                    if c.get("classification")
                })
                skipped_labs = sorted({
                    c.get("classification") for c in skipped.values()
                    if c.get("classification")
                })
                per_structure_label[pid_] = {
                    "classifications_seen": labs,
                    "consensus": (
                        labs[0] if len(labs) == 1
                        else "mixed" if labs
                        else "numbering_mismatch_not_interpretable" if skipped
                        else "no_pocket_to_classify"
                    ),
                    # `overlap_fraction`, NOT `pocket_interface_overlap`. These
                    # entries are the raw classification dicts from
                    # `_classify_site_pocket`; `pocket_interface_overlap` is the
                    # name the field is RENAMED to when it is copied into
                    # results[pid]["pocket_vs_interface"], and reading it here
                    # returned None for every structure at every clustering
                    # value in every run. Not found by the regression — found by
                    # reading the TL1A payload it produced, where the whole map
                    # came back {"1.6": null, "2.4": null} beside overlaps of
                    # 0.267 and 0.133 sitting one level down in per_structure.
                    "overlap_by_clustering": {
                        k: c.get("overlap_fraction")
                        for k, c in byd.items()
                    },
                    # THE FLAG TRAVELS WITH THE LABEL. Both halves of the fix are
                    # here: the excluded entries are named, and the flag is
                    # carried onto the field a caller is told to quote.
                    "numbering_agrees": (
                        None if not byd else
                        all(
                            (c.get("numbering_check") or {}).get(
                                "numbering_agrees") is not False
                            for c in byd.values()
                        )
                    ),
                    "numbering_offsets_applied": sorted({
                        (c.get("numbering_offset_to_partner") or {}).get(
                            "epitope_renumbered_by")
                        for c in byd.values()
                        if (c.get("numbering_offset_to_partner") or {}).get(
                            "applied")
                    }),
                    "excluded_numbering_mismatch": sorted(skipped),
                    "classifications_excluded_numbering_mismatch": skipped_labs,
                    "_numbering_note": (
                        None if not skipped else
                        f"{len(skipped)} classification(s) here were computed "
                        "against non-homologous residues — the two entries are "
                        "not on a common numbering and no constant offset "
                        "recovered it — and are EXCLUDED from `consensus` and "
                        f"from `classifications_seen`: {skipped_labs}. The "
                        "geometric fields (min_distance_to_interface_a, "
                        "enclosure, subunit_enclosure_gain) are unaffected and "
                        "are still in per_structure. Do not quote the excluded "
                        "labels; on TL1A they were allosteric_candidate and "
                        "mixed, both derived from residue-number collisions."
                    ),
                }
                results[pid_]["pocket_vs_interface_consensus"] = (
                    per_structure_label[pid_]
                )
            interface_out["per_structure_consensus"] = per_structure_label
            labels = sorted({
                c.get("classification")
                for byd in per_struct.values()
                for c in _usable(byd)[0].values() if c.get("classification")
            })
            excluded_any = sorted({
                pid_ for pid_, byd in per_struct.items() if _usable(byd)[1]
            })
            interface_out["classifications_seen"] = labels
            interface_out["structures_excluded_numbering_mismatch"] = excluded_any
            interface_out["_aggregation_rule"] = (
                "`classification` here is the CONSENSUS over every pocket "
                "classified in this run, and it is 'mixed' whenever they "
                "disagree. Quote it, or quote per_structure_consensus — never "
                "reach into per_structure and take the first entry. Measured on "
                "8DYG (ligand U5Q): the two symmetry copies gave "
                "allosteric_candidate at overlap 0.22 and orthosteric_candidate "
                "at 0.36, both borderline against the 0.25 boundary. 'mixed' is "
                "the honest answer there and must not be collapsed to one "
                "label; report it as mixed and say the pocket sits on the "
                "boundary. "
                "NUMBERING IS PART OF THIS RULE NOW. A `match_by='seqid'` "
                "comparison between two entries that do not number the protein "
                "the same way matches non-homologous positions, so any "
                "classification whose `numbering_check.numbering_agrees` is "
                "false is EXCLUDED from the consensus here and from "
                "classifications_seen, and is listed in "
                "`excluded_numbering_mismatch` with the label it would have "
                "contributed. Recoverable constant offsets (TL1A +67/+71, "
                "IL-17A +23) are applied upstream before the match, so an entry "
                "reaching this point still flagged is one no constant offset "
                "fixes. A sixth consensus value exists for the case where "
                "EVERY classification on a structure was excluded: "
                "`numbering_mismatch_not_interpretable`. It is not "
                "'no_pocket_to_classify' — pockets were found and classified, "
                "and the classification is what could not be trusted."
            )
            interface_out["classification"] = (
                labels[0] if len(labels) == 1
                else "mixed" if labels
                else "numbering_mismatch_not_interpretable" if excluded_any
                else "no_pocket_to_classify"
            )
            interface_out["_note"] = (
                "destabiliser_candidate is tested FIRST and does not need the "
                "partner at all — burial inside the oligomer is measurable on "
                "the oligomer alone. So a destabiliser call plus an interface "
                "overlap of 0.00, which is what TNF-alpha's SPD304 site gives "
                "against the TNFR2 epitope of 3ALQ, means the pocket is inside "
                "the trimer and nowhere near the receptor epitope. Blocking "
                "TNF/TNFR is not the mechanism; displacing a subunit is."
            )
    interface_out.pop("_epitope", None)
    interface_out.pop("_partner_resnames", None)

    # =======================================================================
    # STAGE 5 — MDPOCKET. The site fixed by construction.
    # =======================================================================
    mdpocket_out: dict = {
        "mdpocket_status": "not_run",
        "mdpocket_reason": "run_mdpocket=False",
    }
    if run_mdpocket:
        try:
            mdpocket_out = _mdpocket_ensemble(
                work,
                {p: prepped_by_pid[p] for p in pdb_ids if p in prepped_by_pid},
                donor["pdb_id"] if donor else None,
                donor["prepped"] if donor else None,
                donor["ligand_xyz"] if donor else None,
                donor["site_chains"] if donor else None,
                target_chains={
                    p: tgt_chains_by_pid[p]
                    for p in pdb_ids if p in tgt_chains_by_pid
                },
            )
            if donor_error:
                mdpocket_out["site_donor_error"] = donor_error
        except Exception as exc:  # noqa: BLE001
            # Belt and braces. _mdpocket_ensemble reports its own failures as
            # data; anything that escapes it still must not cost the run.
            mdpocket_out = {
                "mdpocket_status": "failed",
                "mdpocket_reason": f"{type(exc).__name__}: {exc}",
            }

    # Ensemble spread — volume is the reproducible quantity, druggability is not.
    vols: list[float] = []
    jaccards: dict[str, float] = {}
    n_ligand_confirmed, n_pooled, n_signature_unreliable = 0, 0, 0
    n_signature_foreign = 0
    polymer_ligand_controls = {
        pid: r["polymer_ligand_control"] for pid, r in results.items()
        if r.get("polymer_ligand_control")
    }
    prank_status_counts: dict[str, int] = {}
    # {clustering value: [(pdb_id, centroid, radius), ...]}
    centroids: dict[str, list[tuple[str, list[float], float | None]]] = {}
    for pid, r in results.items():
        for dkey, d in r["by_clustering"].items():
            status = d.get("prank_status")
            if status:
                prank_status_counts[status] = prank_status_counts.get(status, 0) + 1
            if d.get("site_pocket_centroid"):
                centroids.setdefault(dkey, []).append(
                    (
                        pid,
                        d["site_pocket_centroid"],
                        d.get("site_pocket_radius_from_protein_center_a"),
                    )
                )
            if d.get("site_pocket_ligand_site_jaccard") is not None:
                jaccards[f"{pid}@D{dkey}"] = d["site_pocket_ligand_site_jaccard"]
            sp = d["site_pocket"]
            if sp:
                n_pooled += 1
                if d["site_pocket_selected_by"] == "ligand_site_jaccard":
                    n_ligand_confirmed += 1
                elif (
                    d["site_pocket_selected_by"]
                    == "site_signature_unreliable_homooligomer"
                ):
                    n_signature_unreliable += 1
                elif (
                    d["site_pocket_selected_by"]
                    == "site_signature_unreliable_foreign_polymer"
                ):
                    n_signature_foreign += 1
                if sp.get("volume"):
                    vols.append(sp["volume"])
                # AND NOT THE DRUGGABILITY SCORE. There used to be a `drugs`
                # list here, accumulated in exactly the same shape as `vols`
                # one line apart, and that visual symmetry is the whole reason
                # the type error looked like a measurement. Volume is absolute
                # and pools; druggability is normalised inside each file and
                # does not. See ensemble.druggability._removed_pooled_min_max.

    # THE CONTROL. A pocket-matching step is a measurement and needs one: two
    # pockets sharing residue numbers can be 12 A apart and no overlap fraction
    # will say so. A large maximum pairwise distance means the "same site" is
    # not the same site, and the pooled spread above is pooling different
    # pockets.
    # `max_pairwise_centroid_distance_a` IS GONE, DELIBERATELY. It compared
    # pocket centroids across structures that this module does not superpose, so
    # it was the sum of a real site displacement and the two entries' arbitrary
    # rigid-body offsets — the IRAK4 run reported 103.9 A, which is not a
    # measurement of anything. It was presented as THE CONTROL, which made it
    # worse than a stray number.
    #
    # This is not a hypothetical class of error for this project: comparing
    # pockets across structures without a common frame is exactly what produced
    # the 7.7 A off-site tracking that retracted the 651-fold claim. Documenting
    # the caveat beside the number was already tried and the number was quoted
    # anyway. `max_radius_difference_a` — each pocket's distance from its own
    # structure's protein centre — measures the same thing, is frame-invariant,
    # and already existed. It is the control.
    centroid_by_d: dict[str, dict] = {}
    all_radius_diff: list[float] = []
    for dkey, entries in sorted(centroids.items()):
        radii = [e[2] for e in entries if e[2] is not None]
        diff = round(max(radii) - min(radii), 2) if len(radii) > 1 else None
        centroid_by_d[dkey] = {
            "n_structures_with_site_pocket": len(entries),
            # Frame-independent: same quantity measured from each structure's
            # own centre, so it survives the fact that two PDB entries are not
            # deposited in a common coordinate frame. THIS IS THE CONTROL.
            "radius_from_protein_center_a": {e[0]: e[2] for e in entries},
            "max_radius_difference_a": diff,
            # Raw centroids are kept as INPUTS, under a name that says which
            # frame they are in. Do not difference them across entries.
            "centroids_in_own_deposited_frame": {e[0]: e[1] for e in entries},
            "_centroids_frame_warning": (
                "These are in each entry's OWN deposited frame. A difference "
                "between two of them is not a distance between two pockets; it "
                "also contains the two crystals' rigid-body offset. Note also "
                "that a centroid of exactly x=y=z is not an artifact: it is an "
                "on-axis pocket in an assembly whose 3-fold runs along the body "
                "diagonal — 2QE3's assembly operators are literally x,y,z / "
                "z,x,y / y,z,x, so any C3-symmetric cavity has equal "
                "coordinates and zero spread across clustering values. That is "
                "the crystal frame showing through, and it is the same reason "
                "cross-entry centroid distances are meaningless."
            ),
        }
        if diff is not None:
            all_radius_diff.append(diff)

    return {
        "structures": results,
        # The four stages after fpocket, each independently reported. A caller
        # must be able to tell "this did not run" from "this ran and found
        # nothing" from "this died", which is why every one of them carries its
        # own status and reason rather than an absent key.
        "disorder": disorder_out,
        "cryptic": cryptic_out,
        "pocket_vs_interface": interface_out,
        "mdpocket": mdpocket_out,
        "stage_status": {
            "prank": None,  # per structure per D; see prank_status_counts
            "disorder": disorder_out.get("disorder_status"),
            "cryptic": cryptic_out.get("cryptic_status"),
            "interface": interface_out.get("interface_status"),
            "mdpocket": mdpocket_out.get("mdpocket_status"),
        },
        "ensemble": {
            "n_structures": len(pdb_ids),
            "clustering_swept": list(D_VALUES),
            # Which of the pooled pockets are the site and which are a guess.
            "site_pockets_pooled": n_pooled,
            "site_pockets_ligand_confirmed": n_ligand_confirmed,
            "site_pockets_signature_unreliable_homooligomer": n_signature_unreliable,
            "site_pockets_signature_unreliable_foreign_polymer": n_signature_foreign,
            # THE PAIR, PER STRUCTURE. Never one number: a site that exists only
            # while its polymer ligand is present is an induced-fit / occluded
            # site, not an absent one, and `induced_fit_signal` true must force
            # `cryptic_pocket_risk: high`. See `_polymer_ligand_control`.
            "polymer_ligand_control": polymer_ligand_controls,
            "induced_fit_signal_structures": sorted(
                pid for pid, c in polymer_ligand_controls.items()
                if c.get("induced_fit_signal")
            ),
            # The jaccard VALUES, not just the count of structures selected by
            # them. `tractability.ligand_site_jaccard` in the dossier had no
            # source to read: the number was computed per pocket, used to pick
            # the site pocket, and then never surfaced above the individual
            # pocket dict. A basis is not evidence without its value.
            "ligand_site_jaccard_by_structure": jaccards,
            "ligand_site_jaccard": (
                {
                    "min": min(jaccards.values()),
                    "max": max(jaccards.values()),
                    "n": len(jaccards),
                }
                if jaccards
                else None
            ),
            "site_signature": {
                "source": signature_source,
                "n_residues_in": signature_n_residues_in,
                "n_distinct_numbers": len(site_signature),
                # 19 residues collapsing to 11 numbers IS the homotrimer
                # problem, stated as a number rather than as a warning.
                "collapsed_by": (
                    signature_n_residues_in - len(site_signature)
                    if signature_n_residues_in
                    else None
                ),
                "foreign_polymer_residues_dropped": signature_foreign_dropped,
                "foreign_polymer_residues": signature_foreign_residues,
                "_why_foreign": (
                    "collapsed_by counts numbers lost to IDENTICAL protomers. "
                    "This counts numbers imported from a DIFFERENT polymer. "
                    "They are different failures and only the first one was "
                    "guarded. 8QFZ:LFI reported collapsed_by 0 while 9 of its "
                    "13 signature residues belonged to a bicyclic peptide, "
                    "numbered from 1, whose numbers 11-22 also exist on TSLP "
                    "and mean something else. Downstream that produced a "
                    "max_radius_difference_a of 33.52 A and per-structure site "
                    "ranks of 5/7/22/1/36."
                ),
                "donor_homo_oligomer": signature_donor_homo,
                "_warning": (
                    "The signature is a set of residue NUMBERS with chain "
                    "identity discarded. On a homo-oligomer the protomers "
                    "triplicate every number, so a C3-symmetric site cannot be "
                    "resolved in principle and any pocket carrying those "
                    "numbers matches. Structures whose basis is "
                    "site_signature_unreliable_homooligomer or "
                    "site_signature_unreliable_foreign_polymer must not be "
                    "pooled as one site; check site_centroid_control before "
                    "quoting a spread over them."
                ),
            },
            "prank_status_counts": prank_status_counts,
            "site_centroid_control": {
                # THE CONTROL IS max_radius_difference_a AND NOTHING ELSE HERE.
                "max_radius_difference_a": (
                    max(all_radius_diff) if all_radius_diff else None
                ),
                "per_clustering": centroid_by_d,
                "_note": (
                    "A pocket-matching step is a measurement and this is its "
                    "control. Two pockets sharing residue numbers can be 12 A "
                    "apart and an overlap fraction will not tell you. A large "
                    "value here means the 'same site' across the ensemble is "
                    "not the same site and the pooled spread below is pooling "
                    "different pockets. The quantity is each pocket's distance "
                    "from its OWN structure's protein centre, differenced "
                    "across structures — frame-invariant, so it survives the "
                    "fact that two PDB entries are not deposited in a common "
                    "coordinate frame."
                ),
                "_removed_max_pairwise_centroid_distance_a": (
                    "REMOVED, and deliberately not replaced by a null. It "
                    "differenced pocket centroids across structures this module "
                    "does not superpose, so it was a real site displacement "
                    "plus two arbitrary rigid-body offsets — an IRAK4 run "
                    "reported 103.9 A, which is not a measurement of anything, "
                    "under the heading of a control. Comparing pockets across "
                    "structures without a common frame is the exact error that "
                    "retracted this project's 651-fold claim, and a caveat "
                    "printed beside the number did not stop it being quoted. "
                    "Use max_radius_difference_a, or superpose first (the "
                    "mdpocket block does) and quote from there."
                ),
            },
            "_pooling_caveat": (
                "Values are pooled across structures AND clustering values. "
                "For a structure with no drug-like ligand there is no ligand "
                "site to match against, so its 'site pocket' is only the "
                "most druggable pocket anywhere in the chain — it need not be "
                "the site the holo structures point at. Check "
                "site_pocket_selected_by per structure before quoting a "
                "spread as being about one site."
            ),
            "_pooling_caveat_2_trustworthy_basis_is_not_enough": (
                "site_pocket_selected_by == 'ligand_site_jaccard' is a "
                "PER-STRUCTURE guarantee and it does NOT make pooling across "
                "structures safe. Measured on IL-17A: three structures all "
                "selected by ligand_site_jaccard were nonetheless not one site "
                "— 9SQX spans residues 85-142 across both chains, 8DYG spans "
                "A/107-148 plus B/68-104 (a different location), and 8USS is a "
                "MONOMER assembly in which the groove is only half present, so "
                "fpocket buries it at rank 6 of 6 with druggability 0.001. That "
                "0.001 alone produced a 930x pooled range, and "
                "max_radius_difference_a came back at 16.61 A and flagged it. "
                "This is the retracted-651x failure mode recurring WITH a "
                "trustworthy selection basis. So the basis is necessary and not "
                "sufficient: read site_centroid_control.max_radius_difference_a "
                "as well, and do not pool across structures whose assemblies "
                "differ in whether the site is even present. A merge_suspected "
                "or a volume above ~1000 A^3 is the other half of the same "
                "check — at D=2.4 the same IL-17A site came out at 1831 A^3, so "
                "its 0.930 druggability is a merged-site artifact."
            ),
            "volume_a3": {
                "min": min(vols) if vols else None,
                "max": max(vols) if vols else None,
                "spread_pct": (
                    round(100 * (max(vols) - min(vols)) / max(vols), 1)
                    if vols and max(vols)
                    else None
                ),
            },
            "druggability": {
                # min / max / fold_range ARE GONE. They pooled a
                # WITHIN-STRUCTURE-NORMALISED quantity across structures, which
                # is a type error and not a weak measurement — see
                # `_removed_pooled_min_max`. What replaces them is the same
                # information in the only two forms the quantity supports: the
                # rank within each structure, and the range across the D SWEEP
                # WITHIN one structure.
                "_removed_pooled_min_max": (
                    "REMOVED, and deliberately not replaced by nulls. `min`, "
                    "`max` and `fold_range` pooled fpocket's druggability score "
                    "across structures AND across clustering values. That score "
                    "is normalised INSIDE each file: pocket.c:736-756 min-max "
                    "normalises its dominant term, mean_loc_hyd_dens_norm, over "
                    "the CURRENT STRUCTURE'S OWN pocket list whenever "
                    "n_pockets > 1, and the hardcoded (mlhd-8.23)/(24.20-8.23) "
                    "at pocket.c:780 is the single-pocket branch that never "
                    "fires here (4-324 pockets per structure); pscoring.c:325 "
                    "feeds it to the logistic. So the number answers 'how does "
                    "this pocket rank against the others in this file' and "
                    "nothing else, and pooling it across files is a TYPE ERROR "
                    "rather than a noisy measurement. RORgt proves it on one "
                    "protein at one site: 4NB6's site MLHD of 30.722 IS that "
                    "structure's maximum, normalises to 1.0 and scores 0.827, "
                    "while 6C1P's 19.0 against a maximum of 52.767 normalises "
                    "to 0.36 and scores 0.009 — a 90-fold gap from which other "
                    "pockets happened to co-exist in the file. This is also the "
                    "operation that manufactured the withdrawn 651-fold "
                    "TNF-alpha spread, with no matcher error required, and a "
                    "run of this very ensemble emitted fold_range 195.7 while "
                    "the two _warning strings beside it argued against it. A "
                    "caveat next to a number does not stop the number being "
                    "quoted; removing the number does. Read "
                    "`site_pocket_rank_by_structure` (rank + n_pockets + PDB "
                    "ID) and, if you need a range, "
                    "`druggability_range_within_structure_across_d`, which is "
                    "the one range that is legitimate. VOLUME is unaffected and "
                    "is still pooled above: it is an absolute physical quantity "
                    "and it does travel between structures."
                ),
                # THE TWO LEGITIMATE FORMS.
                "site_pocket_rank_by_structure": {
                    f"{pid}@D{dkey}": {
                        "fpocket": (dd.get("site_pocket") or {}).get("rank"),
                        "prank": (dd.get("site_pocket") or {}).get("prank_rank"),
                        "n_pockets": dd.get("n_pockets"),
                        "structure_pdb_id": pid,
                        "selected_by": dd.get("site_pocket_selected_by"),
                        "druggability_score": (
                            (dd.get("site_pocket") or {}).get("druggability_score")
                        ),
                    }
                    for pid, rr in results.items()
                    for dkey, dd in (rr.get("by_clustering") or {}).items()
                    if dd.get("site_pocket")
                },
                "druggability_range_within_structure_across_d": {
                    pid: {
                        "min": min(vv), "max": max(vv),
                        "fold_range": (
                            round(max(vv) / min(vv), 1) if min(vv) > 0 else None
                        ),
                        "n_clustering_values": len(vv),
                    }
                    for pid, vv in (
                        (pid, [
                            (dd.get("site_pocket") or {}).get("druggability_score")
                            for dd in (rr.get("by_clustering") or {}).values()
                            if (dd.get("site_pocket") or {}).get(
                                "druggability_score") is not None
                        ])
                        for pid, rr in results.items()
                    )
                    if len(vv) > 1
                },
                "_why_rank": (
                    "'rank 1 of 30 in 6OIM' is the claim. The value may sit "
                    "beside the rank; the rank is what is asserted. fpocket's "
                    "own rank and PRANK's are TWO WITHIN-STRUCTURE ORDERINGS on "
                    "the same footing — report both, replace neither, and "
                    "report a disagreement as a disagreement."
                ),
                "_warning": (
                    "Druggability is NOT reproducible across structures or "
                    "clustering values. Measured on an apo TNF-alpha ensemble: "
                    "fixing the site BY CONSTRUCTION (one grid definition "
                    "applied to every superposed structure) rather than by "
                    "post-hoc residue matching cut the across-ensemble VOLUME "
                    "CV from ~28% to ~10% (measured 28.1% at D=1.6 against "
                    "9.9%), roughly a 2.8-fold reduction. Both figures carry "
                    "about 1 percentage point of fpocket Monte-Carlo volume "
                    "noise — three identical reruns of one 5-structure ensemble "
                    "gave CVs of 12.1/11.3/10.8% — so quote them to two "
                    "significant figures and do not read a CV difference under "
                    "~1pp as a difference between sites. The improvement is "
                    "real; the third digit is not. NOTE the CV above was "
                    "measured on site_from_density, which is not the ligand "
                    "site; see mdpocket.sites. REPORT DRUGGABILITY AS A RANK "
                    "AMONG THAT STRUCTURE'S POCKETS, WITH THE POCKET COUNT AND "
                    "THE PDB ID — an earlier version of this line said 'report "
                    "druggability as a range', which is void: a range pooled "
                    "across structures is the type error described in "
                    "_removed_pooled_min_max, and only the D-sweep range within "
                    "ONE structure is a range at all. Never drive a verdict "
                    "from the value. Volume is the more reliable number and the "
                    "only one of the two that may be compared across "
                    "structures."
                ),
                "_retracted": (
                    "An earlier version of this warning carried a large "
                    "fold-spread figure across five apo TNF-alpha structures "
                    "'of the same site'. THAT FIGURE IS WITHDRAWN and is "
                    "deliberately not reproduced here, so that it cannot be "
                    "lifted out of this payload. mdpocket showed the "
                    "residue-number matcher that produced it was tracking a "
                    "pocket 7.7 A away from the site it claimed, with 12.2 A of "
                    "internal inconsistency between structures. It was never a "
                    "measurement of one site."
                ),
            },
        },
        "method": {
            "tool": "fpocket 4.2.3 (conda-forge)",
            "tools": {
                "detection": "fpocket 4.2.3 (conda-forge)",
                "rescoring": f"P2Rank {P2RANK_VERSION} rescore (rescore_2024)",
                "site_fixed_by_construction": "mdpocket (ships with fpocket)",
                "cryptic_mechanism": "cryptic_analysis.py (gemmi + numpy)",
                "interface": "interface_analysis.py (gemmi + numpy)",
                "disorder": "metapredict 3.0.2, CPU torch",
            },
            "clustering_swept": list(D_VALUES),
            "ligand_classification": {
                "tool": "ligand_filter.classify_record (stdlib, no RDKit)",
                "records": "RCSB data.rcsb.org/rest/v1/core/chemcomp (SMILES, "
                           "formula, type)",
                "basis": "chemistry of the component, from its SMILES graph",
                "accuracy": "259/262 on ground truth; 61/70 on a blind "
                            "held-out set with zero false positives",
                "replaces": (
                    "a >=18 heavy-atom floor plus two hardcoded comp_id lists. "
                    "Both are deleted. ADP has 27 heavy atoms and so does "
                    "A1IPJ, a genuine inhibitor — no size threshold separates "
                    "them, and no list is ever complete (CHAPS/CPS was simply "
                    "missing). Identity filtering gave 16 holo / 8 apo on "
                    "NLRP3 where a size window gave 19 / 5."
                ),
                "verdicts": list(getattr(_lf_verdicts(), "VERDICTS", ())),
                "undetermined_is_not_apo": (
                    "a component whose record could not be retrieved leaves "
                    "the entry at tier 'undetermined'; see structures.<ID>."
                    "holo_call.undetermined"
                ),
            },
            "ligand_site": "5.0 A heavy-atom shell, single ligand copy, kept chains only",
            "prep": "protein only, altloc A/blank, hydrogens stripped",
            # Provenance: legacy PDB truncates comp_ids to 3 characters and is
            # not issued at all for newer entries, so nothing here is derived
            # from it. The PDB written for fpocket comes out of the mmCIF.
            #
            # DERIVED FROM WHAT WAS ACTUALLY FETCHED, not restated as a
            # constant. The constant said `<ID>.cif`, which is the ASYMMETRIC
            # UNIT, while `_fetch` tries `<ID>-assembly1.cif` first and the
            # per-structure `structure_source` correctly reported `assembly1`.
            # The same payload contradicted itself and the method block was the
            # half that was wrong, so a reader who trusted it concluded we
            # scanned ASUs. That is not cosmetic: ASU-versus-assembly is a
            # documented wrong answer of ours — 9SQX's ASU holds two dimers,
            # scoring all four fused them and moved the real ligand site from
            # rank 1 to rank 9.
            "source_format": (
                "mmCIF, parsed with gemmi. Preferred biological assembly "
                "(files.rcsb.org/download/<ID>-assembly1.cif) with fallback to "
                "the asymmetric unit (<ID>.cif); which one was used is recorded "
                "per structure in structures.<ID>.structure_source, and "
                "summarised in source_used_by_structure below."
            ),
            "source_used_by_structure": {
                pid: r.get("structure_source", "unknown")
                for pid, r in results.items()
            },
            "source_used": sorted(
                {r.get("structure_source", "unknown") for r in results.values()}
            ),
        },
        # PREDICTED ENTRIES IN THIS SCAN, off the default path (CLAUDE.md 4c).
        # Empty labels means an all-experimental scan.
        "predicted_structures_used": {
            "labels": sorted(predicted_ids),
            "uniprot_accession": target_accession if predicted_ids else None,
            "_note": (
                "These labels are ESMFold models, not deposited structures "
                "(structures.<label>.structure_origin == 'esmfold_predicted'). "
                "ANY pocket reported on one is a MODEL pocket — one level less "
                "trustworthy than the rest of the computed axis — and a "
                "predicted entry never donates a site signature and is never a "
                "holo or mdpocket site donor. Reached only when a target has no "
                "experimental structure and no usable homolog; see "
                "structure-select/predicted_structure_fallback.py."
            ),
        },
    }


@app.local_entrypoint()
def main(
    pdb_ids: str = "6OIM,4OBE",
    ligand_codes: str = "",
    uniprot_accession: str = "",
    partner_structures: str = "",
    mdpocket_site_donor: str = "",
    chains: str = "",
    site_residues: str = "",
    run_disorder: bool = True,
    run_cryptic: bool = True,
    run_mdpocket: bool = True,
    out: str = "",
):
    """Smoke test: the KRAS holo/apo pair the calibration was built on.

    Expected: 6OIM's switch-II pocket recovers the MOV site with high Jaccard
    at one D; 4OBE shows the same site collapsed. If 6OIM comes back with a
    low-overlap site, the prep or the parse is broken, not the biology.

    `ligand_codes` DEFAULTS TO EMPTY on purpose. It used to default to "MOV",
    which meant `modal run modal_app.py --pdb-ids 6OIM,4OBE` silently passed a
    ligand code and there was no way to exercise auto-derivation from the CLI at
    all — the one path most likely to be wrong was the one that could not be
    tested. Empty means the function derives MOV from 6OIM itself, which is the
    behaviour that should hold.

    EVERY PARAMETER `pocket_scan` TAKES IS NOW REACHABLE FROM HERE. `chains`,
    `site_residues` and the three stage switches existed on the function and on
    nothing that could call it, which made the single most informative control
    on an oligomer unreachable without editing this file: deleting the third
    protomer is the experiment that separates "the cavity is too small" from "a
    protomer is standing in it", and on TNF-alpha it moves the SPD304 site from
    0.00 A^3 to ~280-550 A^3. TL1A's axial cavity was reported at 49.5-141.1
    A^3 intact and the control was never run, because the CLI could not ask.

        --chains '6OIM=A;1TNF=A,B'
        --site-residues 57,58,59,60,61

    `--out` WRITES THE JSON TO A FILE, and it is the right way to capture it.
    Modal prints its own progress banner and a trailing "Stopping app..." to
    stdout, interleaved with anything printed here, so
    `modal run modal_app.py ... > out.json` produces INVALID JSON and every
    consumer has had to `raw_decode` from the first `{`. With `--out` the
    payload never touches stdout; without it, it goes to stderr, which the
    banner does not share.
    """
    codes = [c.strip() for c in ligand_codes.split(",") if c.strip()]
    partners = [p.strip() for p in partner_structures.split(",") if p.strip()]
    # "6OIM=A;1TNF=A,B" -> {"6OIM": ["A"], "1TNF": ["A", "B"]}
    chain_map: dict[str, list[str]] = {}
    for part in chains.split(";"):
        part = part.strip()
        if not part or "=" not in part:
            continue
        pid, _, cs = part.partition("=")
        picked = [c.strip() for c in cs.split(",") if c.strip()]
        if picked:
            chain_map[pid.strip()] = picked
    resi = [int(r) for r in site_residues.replace(",", " ").split() if r.strip()]
    result = pocket_scan.remote(
        pdb_ids=[p.strip() for p in pdb_ids.split(",")],
        chains=chain_map or None,
        ligand_codes=codes or None,
        site_residues=resi or None,
        uniprot_accession=uniprot_accession.strip() or None,
        partner_structures=partners or None,
        mdpocket_site_donor=mdpocket_site_donor.strip() or None,
        run_disorder=run_disorder,
        run_cryptic=run_cryptic,
        run_mdpocket=run_mdpocket,
    )
    text = json.dumps(result, indent=2)
    if out.strip():
        Path(out.strip()).write_text(text + "\n")
        print(f"wrote {out.strip()}", file=sys.stderr)
    else:
        # stderr, not stdout: Modal owns stdout and puts a banner on either
        # side of whatever is printed there.
        print(text, file=sys.stderr)
