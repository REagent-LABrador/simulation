"""Ground-truth harness for cryptic_analysis.

Reproduces the two hand-measured calibration cases -- KRAS (4OBE apo vs 6OIM
holo, ligand MOV) and TNF-alpha (1TNF apo trimer vs 2AZ5 holo A/B dimer,
ligand comp_id `307`) -- and prints a got-vs-expected table.

    CRYPTIC_FIXTURES=/path/to/structures python3 test_cryptic_analysis.py

THESE ARE CALIBRATION-PROTOCOL EXPECTATIONS. THEY ARE NOT WHAT THE PIPELINE
OUTPUTS, AND THE TWO MUST NOT BE READ AS EACH OTHER.

``run_kras`` and ``run_tnf`` deliberately drive the HAND-CALIBRATION protocol:
``auto_trim=False``, ``match_residue_names=False``, and for KRAS an explicit
switch I / switch II exclusion with a fixed 1-166 fit range. Under that protocol
the displacements are 8.8 A and 1.62 A, which is what the assertions below pin.

The DEPLOYED default -- what ``modal_app.pocket_scan`` actually runs, with
auto-trim on and residue-name matching on, and with no target-specific region
list, because none is hardcoded -- lands 0.1-0.2 A BELOW those figures:
**8.65 A for KRAS and ~1.55 A for TNF-alpha**. A bare ``1.62 +/- 0.03`` in this
file reads as "this is what the pipeline outputs" and it is not; the deployed
run would fail that tolerance, correctly.

Nothing about the science turns on it. ``run_kras_auto`` and ``run_tnf_auto``
run the default and assert the only things that carry a decision -- ``mechanism``
and ``is_cryptic`` -- and those are IDENTICAL under both protocols. The risk this
note exists to close is purely that someone quotes 8.83 or 1.62 as a figure this
pipeline reproduces. ``pocket_scan`` reports the default in
``cryptic.max_backbone_ca_displacement_a`` and re-runs the calibration protocol
separately into ``calibration_protocol``; say which one you are quoting.

Two notes on how the expectations are pinned:

* The two hand calibrations used *different* superposition protocols. KRAS
  excluded switch I (25-40) and switch II (57-75) explicitly and fitted
  residues 1-166; TNF-alpha excluded nothing, because TNF-alpha has no mobile
  region to exclude. Both are reproduced here with the protocol each was
  measured under, and each case is additionally run under the module's
  zero-knowledge default to show the mechanism call does not depend on it.
* Both cases set ``match_residue_names=False`` because the hand fits matched on
  residue number alone. The module defaults to True, which drops the four KRAS
  construct differences (G12C, C51S, C80L, C118S) and the two TNF L143D
  positions from the fit -- changing the fitted-atom count but not the RMSD.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from cryptic_analysis import analyze_cryptic_mechanism  # noqa: E402

# Fixture root must contain pockets/{4OBE,6OIM}.pdb and tnfa/{1TNF,2AZ5}.pdb.
FIXTURES = os.environ.get(
    "CRYPTIC_FIXTURES",
    "/private/tmp/claude-502/-Users-bb-repos/"
    "4b0f8244-20f1-454e-86bd-3f24a4e596c5/scratchpad")
PK = os.path.join(FIXTURES, "pockets")
TN = os.path.join(FIXTURES, "tnfa")

ROWS = []


def check(case, label, got, exp, tol=None, fmt="{}"):
    if tol is None:
        ok = got == exp
    else:
        ok = got is not None and abs(got - exp) <= tol
    show = (lambda v: "None" if v is None else fmt.format(v))
    ROWS.append((case, label, show(exp), show(got), "PASS" if ok else "FAIL"))
    return ok


# ---------------------------------------------------------------- KRAS
def run_kras():
    """4OBE (apo) vs 6OIM (holo), ligand MOV.

    Calibration protocol: fit on G-domain core CA excluding switch I (25-40)
    and switch II (57-75), residues 1-166, no auto-trim. apo_chains restricted
    to A because 4OBE chain B is a second crystallographic copy, not a subunit
    of the biological assembly.
    """
    switch = list(range(25, 41)) + list(range(57, 76))
    r = analyze_cryptic_mechanism(
        os.path.join(PK, "4OBE.pdb"), os.path.join(PK, "6OIM.pdb"), "MOV",
        apo_chains=["A"], exclude_residues=switch, fit_residue_range=(1, 166),
        auto_trim=False, match_residue_names=False)
    c = "KRAS"
    sup, site, con, sc = r["superposition"], r["site"], r["contacts"], r["self_control"]
    check(c, "core CA fitted", sup["n_fitted_ca"], 128)
    check(c, "core CA RMSD (A)", sup["core_ca_rmsd"], 0.86, 0.02, "{:.2f}")
    check(c, "apo protein atoms <2.0 A", con["n_protein_atoms"], 6)
    check(c, "apo protein atoms <3.0 A", r["contacts_wide"]["n_protein_atoms"], 20)
    check(c, "self-control pairs <2.0 A", sc["contact_pairs"], 1)
    check(c, "self-control prot atoms <3.0 A", sc["n_protein_atoms_wide"], 2)
    # Tyr64 points away from the ligand and is outside the 5 A site shell used
    # for classification, so pull the per-residue displacement table from a
    # wider-shell run. The transform is identical; only the reported set grows.
    rw = analyze_cryptic_mechanism(
        os.path.join(PK, "4OBE.pdb"), os.path.join(PK, "6OIM.pdb"), "MOV",
        apo_chains=["A"], exclude_residues=switch, fit_residue_range=(1, 166),
        auto_trim=False, match_residue_names=False, site_radius=8.0,
        compute_free_volume=False)
    disp = {p["resi"]: p["ca_displacement"] for p in rw["site"]["per_residue"]}
    for resi, exp in [(60, 3.1), (62, 4.4), (63, 8.8), (64, 4.3), (68, 3.6)]:
        check(c, f"CA displacement res {resi} (A)", disp.get(resi), exp, 0.1, "{:.2f}")
    # CALIBRATION PROTOCOL, not the deployed default. The default measures
    # 8.65 A on this pair and would fail this tolerance. See the module
    # docstring; `mechanism` and `is_cryptic` below are protocol-independent.
    check(c, "max CA displacement (A) [calibration protocol, NOT the "
             "deployed default, which gives 8.65]",
          site["max_ca_displacement"], 8.8, 0.1, "{:.2f}")
    check(c, "mechanism", r["mechanism"], "loop_or_backbone_motion")
    check(c, "is_cryptic", r["is_cryptic"], True)
    check(c, "self-control passed", sc["passed"], True)
    check(c, "residue-name mismatches found",
          sup["n_residue_name_mismatches"], 4)
    return r


def run_kras_auto():
    """Same pair under the module's zero-knowledge default (auto-trim only)."""
    r = analyze_cryptic_mechanism(
        os.path.join(PK, "4OBE.pdb"), os.path.join(PK, "6OIM.pdb"), "MOV",
        apo_chains=["A"])
    c = "KRAS(auto)"
    check(c, "mechanism", r["mechanism"], "loop_or_backbone_motion")
    check(c, "is_cryptic", r["is_cryptic"], True)
    return r


# ------------------------------------------------------------- TNF-alpha
def run_tnf():
    """1TNF (apo trimer) vs 2AZ5 (holo A/B dimer), ligand comp_id `307`.

    Calibration protocol: fit on every equivalent CA of the two mapped chains,
    no mobile-region exclusion (TNF-alpha has no mobile region to exclude).
    """
    r = analyze_cryptic_mechanism(
        os.path.join(TN, "1TNF.pdb"), os.path.join(TN, "2AZ5.pdb"), "307",
        holo_chains=["A", "B"], ligand_chain="A", auto_trim=False,
        match_residue_names=False)
    c = "TNFa"
    sup, site, con, sc = r["superposition"], r["site"], r["contacts"], r["self_control"]
    fv = r["free_volume"]
    check(c, "CA atoms fitted", sup["n_fitted_ca"], 272)
    check(c, "CA RMSD (A)", sup["core_ca_rmsd"], 1.12, 0.02, "{:.2f}")
    check(c, "displaced apo chain", r["inputs"]["displaced_apo_chains"], ["C"])
    check(c, "total contacts <2.0 A", con["contact_pairs"], 66)
    check(c, "ligand atoms in contact", con["n_ligand_atoms"], 24)
    check(c, "ligand heavy atoms", r["inputs"]["ligand_heavy_atoms"], 40)
    check(c, "min distance (A)", con["min_distance"], 0.55, 0.02, "{:.2f}")
    check(c, "contacts from chain C", con["by_category"]["displaced_chain"], 40)
    check(c, "side-chain contacts", con["by_category"]["sidechain"], 26)
    check(c, "backbone contacts", con["by_category"]["backbone"], 0)
    sc_res = sorted({(x["chain"], x["resi"]) for x in con["residues"]
                     if x["category"] == "sidechain"})
    check(c, "side-chain clash residues", [(str(a), b) for a, b in sc_res],
          [("A", 119), ("B", 119)])
    check(c, "self-control pairs <2.0 A", sc["contact_pairs"], 0)
    check(c, "self-control min dist (A)", sc["min_distance"], 2.60, 0.02, "{:.2f}")
    # CALIBRATION PROTOCOL, not the deployed default. The default measures
    # ~1.55 A on this pair and would fail this +/-0.03 tolerance. That is not a
    # regression; the two protocols are different measurements and only the
    # mechanism call, asserted below, is shared between them.
    check(c, "max site CA displacement (A) [calibration protocol, NOT the "
             "deployed default, which gives ~1.55]",
          site["max_ca_displacement"], 1.62, 0.03, "{:.2f}")
    check(c, "ligand vdW volume (A^3)", fv["ligand_vdw_volume_A3"], 426, 3, "{:.0f}")
    check(c, "free vol, holo (%)", fv["holo"] * 100, 100.0, 0.5, "{:.1f}")
    check(c, "free vol, apo trimer (%)", fv["apo_intact"] * 100, 62.1, 1.0, "{:.1f}")
    check(c, "free vol, -chain C (%)",
          fv["apo_minus_displaced_chains"] * 100, 85.3, 1.0, "{:.1f}")
    check(c, "free vol, -C -Tyr119 (%)",
          fv["apo_minus_displaced_chains_and_sidechains"] * 100, 99.8, 0.5, "{:.1f}")
    check(c, "mechanism", r["mechanism"], "subunit_occlusion")
    check(c, "secondary mechanism", r["secondary_mechanism"], "sidechain_occlusion")
    check(c, "is_cryptic", r["is_cryptic"], False)
    check(c, "self-control passed", sc["passed"], True)
    check(c, "residue-name mismatches found",
          sup["n_residue_name_mismatches"], 2)
    return r


def run_tnf_selfcontrol2():
    """Second TNF self-control: ligand transferred A/B -> C/D within the ASU.

    Uses 2AZ5 as both apo and holo, mapping the A/B dimer onto the C/D dimer.
    """
    r = analyze_cryptic_mechanism(
        os.path.join(TN, "2AZ5.pdb"), os.path.join(TN, "2AZ5.pdb"), "307",
        holo_chains=["A", "B"], apo_chains=["C", "D"], ligand_chain="A",
        auto_trim=False, compute_free_volume=False)
    c = "TNFa ctrl2"
    check(c, "A/B->C/D contacts <2.0 A", r["contacts"]["contact_pairs"], 0)
    check(c, "A/B->C/D min dist (A)", r["contacts"]["min_distance"], 2.67, 0.02, "{:.2f}")
    check(c, "mechanism", r["mechanism"], "none")
    return r


def run_tnf_auto():
    r = analyze_cryptic_mechanism(
        os.path.join(TN, "1TNF.pdb"), os.path.join(TN, "2AZ5.pdb"), "307",
        holo_chains=["A", "B"], ligand_chain="A")
    c = "TNFa(auto)"
    check(c, "mechanism", r["mechanism"], "subunit_occlusion")
    check(c, "is_cryptic", r["is_cryptic"], False)
    return r


def main():
    kras = run_kras()
    run_kras_auto()
    tnf = run_tnf()
    run_tnf_selfcontrol2()
    run_tnf_auto()

    w = [max(len(str(r[i])) for r in ROWS + [("case", "quantity", "expected",
                                              "got", "status")])
         for i in range(5)]
    hdr = ("case", "quantity", "expected", "got", "status")
    line = "  ".join("-" * w[i] for i in range(5))
    print("\n" + "  ".join(hdr[i].ljust(w[i]) for i in range(5)))
    print(line)
    last = None
    for r in ROWS:
        if last is not None and r[0] != last:
            print(line)
        print("  ".join(str(r[i]).ljust(w[i]) for i in range(5)))
        last = r[0]
    print(line)
    nfail = sum(1 for r in ROWS if r[4] == "FAIL")
    print(f"\n{len(ROWS) - nfail}/{len(ROWS)} checks passed"
          + (f", {nfail} FAILED" if nfail else ""))

    print("\n--- KRAS clash attribution at 2.0 A (note: zero backbone) ---")
    for x in kras["contacts"]["residues"]:
        print(f"  {x['category']:16s} {x['chain']}/{x['resn']}{x['resi']}: "
              f"{x['atoms']} min {x['min_distance']}")
    print(f"  by_category: {kras['contacts']['by_category']}")
    print(f"  KRAS rationale: {kras['rationale']}")
    print(f"  KRAS crypticity: {kras['crypticity']['reason']}")
    print("\n--- TNF free-volume separation ---")
    fv = tnf["free_volume"]
    for k in ("holo", "apo_intact", "apo_minus_displaced_chains",
              "apo_minus_clashing_sidechains",
              "apo_minus_displaced_chains_and_sidechains"):
        print(f"  {k:44s} {fv[k]*100:5.1f}%")
    print(f"  TNF rationale: {tnf['rationale']}")
    print(f"  TNF crypticity: {tnf['crypticity']['reason']}")
    return 1 if nfail else 0


if __name__ == "__main__":
    raise SystemExit(main())
