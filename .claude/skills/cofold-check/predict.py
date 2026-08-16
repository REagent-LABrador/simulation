"""cofold-check — the four GPU proto-tools, called in process, each returning
its numbers together with the provenance needed to judge them.

Run this under the proto-tools python (it imports ``proto_tools`` directly —
plain Python import, in-process; not MCP, not a CLI). Modal credentials come
from the environment / ``~/.modal.toml``; nothing here reads a hard-coded
dotenv path.

    from predict import cofold_complex, cofold_affinity, esmfold_predict, bioemu_ensemble

WHAT THIS MODULE WILL AND WILL NOT DO
-------------------------------------
It computes, per run, the things that are measurable per run:

* seed dispersion, and WHICH SITE the seeds converged on (``cofold_complex``);
* a rank ordering plus a positive-control log error (``cofold_affinity``);
* the model's own confidence and the inter-chain contact count
  (``esmfold_predict``);
* the linker that had to be inserted to reach a multimer (``bioemu_ensemble``).

It does NOT ship a calibration constant for any of these tools, and it does not
declare any of them unreliable. Observations about the tools are carried in
``OBSERVATIONS`` with their sample size stated, so a reader can weigh them.

The cross-target benchmark landed on 2026-08-15 and ALL THREE of the former
n=1 observations moved: the sealed-pocket claim was overturned as stated (the
pLDDT family drops on 5 of 5 targets; ipTM/ligand-ipTM are the unreliable ones),
the 1.97-log affinity bias was overturned (no offset detectable over 23 pairs,
and within-target ranking is NOT supported on any of 3 targets), and the ESMFold interface failure
was reproduced exactly and shown to be an input-construction artifact. Every one
had failed in the flattering direction. Nothing downstream needed un-picking,
because no correction was ever applied to a returned number — and that is still
true after the benchmark: a benchmarked observation is a statement about a
population of targets, and your run is one draw from it.
"""

from __future__ import annotations

import math
import os
import statistics
import time
from itertools import combinations
from pathlib import Path
from typing import Any

__all__ = [
    "OBSERVATIONS",
    "bioemu_ensemble",
    "cofold_affinity",
    "cofold_complex",
    "esmfold_predict",
]

# ---------------------------------------------------------------------------
# OBSERVATIONS
#
# Things we have actually seen on our own targets. Every entry states its
# sample size and carries an explicit ``benchmarked`` flag. NONE of these is
# applied to a returned number as a correction, a gate or a downweighting. They
# are attached to the payload so a caller can weigh them, and that is all.
#
# As of the 2026-08-15 cross-target benchmark, three entries carry
# ``benchmarked: True`` with a real n (sealed-pocket confidence, 5 targets /
# 30 folds; affinity, 23 pairs / 3 targets; ESMFold at interfaces, 14 complexes
# / 28 runs). Each of those entries also carries a ``supersedes`` field naming
# the n=1 claim it replaced, because all three of the originals were wrong in
# the flattering direction and a reader who remembers the old figure needs to
# see it withdrawn rather than silently absent.
#
# The distinction that still matters: a benchmarked observation is a statement
# about a POPULATION of targets, and your run is one draw from it. Per-run
# quantities (seed dispersion, the positive-control log error, the contact
# count) are recomputed fresh on every call and are the numbers you should
# actually act on.
# ---------------------------------------------------------------------------
# NOTE ON KEY NAMES: the first three keys below are HISTORICAL. They were named
# for the single case each originally rested on (KRAS, tofacitinib, IL-17A) and
# the names are kept so existing lookups still resolve, but the contents are now
# 5-target, 17-pair and 14-complex benchmarks respectively. The payload labels
# they are returned under were renamed to match the contents.
OBSERVATIONS: dict[str, Any] = {
    "kras_sealed_pocket_confidence": {
        "what_was_measured": (
            "Whether any Boltz-2 output signal reports that a binding site has "
            "been destroyed."
        ),
        "how": (
            "A known ligand site sealed shut with nine phenylalanine "
            "substitutions (physically undruggable by construction) under one "
            "uniform rule — pocket-lining residues within 5.0 A of the cognate "
            "ligand in the holo assembly, the nine with the smallest side "
            "chains mutated to Phe — then wild type and mutant each cofolded "
            "with the cognate ligand and compared. Run on FIVE targets: KRAS "
            "switch-II (6OIM/sotorasib), the JAK1 ATP site (4E4N), the "
            "TNF-alpha trimer-interface site (2AZ5/SPD304), the IL-17A dimer "
            "cavity (5HI3) and the BCL-2 BH3 groove (6O0K)."
        ),
        "sample_size": {
            "targets": 5,
            "target_list": ["KRAS", "JAK1", "TNFa", "IL17A", "BCL2"],
            "states_per_target": 2,
            "seeds_per_state": 3,
            "total_folds": 30,
        },
        "replicated_on_other_targets": True,
        "benchmarked": True,
        "benchmark_date": "2026-08-15",
        "generalises": (
            "MEASURED over 5 targets. The DIRECTION generalises — the pLDDT "
            "family drops on 5 of 5. The MAGNITUDE does not: it ranges from "
            "-0.143 (KRAS, visible) to -0.019 (JAK1, invisible in practice). "
            "Still NOT applied to any returned number."
        ),
        "supersedes": (
            "The previous entry here was ONE target (KRAS) at TWO seeds and "
            "reported the sealed mutant scoring HIGHER than wild type on every "
            "pLDDT-family metric (complex pLDDT 0.940 -> 0.957, confidence "
            "0.919 -> 0.927) with a 0.73 A backbone RMSD to wild type against a "
            "1.02 A wild-type-vs-wild-type baseline. That claim is WITHDRAWN. "
            "Its 1.02 A baseline was seed noise from two seeds; the proper "
            "3-seed baseline on KRAS is 0.23 A."
        ),
        "result": {
            "metric_notices_beyond_2x_seed_sd": {
                "confidence_score": "5 of 5",
                "complex_plddt": "5 of 5",
                "complex_iplddt": "5 of 5",
                "iptm": "3 of 5",
                "ligand_iptm": "3 of 5",
                "ptm": "2 of 5",
                "avg_pae": "2 of 5",
                "complex_pde": "2 of 5",
            },
            "complex_plddt_wild_type_to_sealed": {
                "KRAS": [0.9547, 0.8122, -0.143],
                "JAK1": [0.9666, 0.9480, -0.019],
                "TNFa": [0.8585, 0.8295, -0.029],
                "IL17A": [0.7912, 0.7511, -0.040],
                "BCL2": [0.8094, 0.7632, -0.046],
            },
            "_complex_plddt_note": "[wild type, sealed, delta]; seed sd 0.0009-0.0166",
            "ligand_iptm_ROSE_when_pocket_sealed": {
                "TNFa": [0.8639, 0.9064],
                "IL17A_flat": [0.9262, 0.9248],
            },
            "backbone_ca_rmsd_wt_vs_sealed_against_wt_vs_wt_baseline_a": {
                "KRAS": [1.37, 0.23],
                "JAK1": [0.83, 0.26],
                "IL17A": [5.82, 2.98],
                "BCL2_invisible": [4.64, 4.92],
                "TNFa_invisible": [10.80, 14.85],
            },
            "_backbone_note": (
                "notices on 3 of 5. BCL-2 and TNF-alpha are UNRESOLVED, not "
                "negative: their wild-type seed baselines are as large as the "
                "effect (TNF's baseline sd is 10.09 A), so the test has no "
                "power there."
            ),
        },
        "what_it_shows": (
            "Cofolding confidence DOES register pocket destruction — the pLDDT "
            "family on 5 of 5 — but the ligand-facing metrics are the "
            "treacherous ones: iptm/ligand_iptm move the right way on only 3 of "
            "5, and on TNF-alpha ligand_iptm ROSE when the pocket was sealed "
            "shut. And the magnitudes are mostly unusable: JAK1 fell 0.967 to "
            "0.948, which still reads as an excellent model, so only KRAS moved "
            "enough for an operator to notice unaided."
        ),
        "what_it_does_not_show": (
            "That confidence is a druggability readout. Statistical "
            "detectability is not operational detectability, and nothing tells "
            "you the drop is there without the wild-type control beside it. "
            "NEVER read a high confidence value as evidence that a predicted "
            "pocket is real."
        ),
    },
    "tofacitinib_affinity_error": {
        "what_was_measured": (
            "How far the Boltz-2 affinity value head lands from measured "
            "potency, and whether it can order actives or only separate them "
            "from non-binders."
        ),
        "how": (
            "23 protein-ligand pairs with ChEMBL consensus potencies (assay "
            "confidence 9, relation '='), across JAK1 (12), EGFR (6) and "
            "BCL-2 (5), spanning pChEMBL 5.30-10.36, plus 12 decoys on JAK1. "
            "Diverse chemistry only -- no congeneric series (see "
            "'untested_case')."
        ),
        "sample_size": {
            "targets": 3,
            "target_list": ["JAK1", "EGFR", "BCL2"],
            "pairs": 23,
            "actives_for_triage": 12,
            "decoys": 12,
            "ground_truth": "ChEMBL consensus, not a single cherry-picked value",
        },
        "provenance": {
            "artifact": "out/claim2_{JAK1,EGFR,BCL2}.json",
            "regenerated_by": "analyze.py 2",
            "measured_on": "2026-08-15",
            "remaining_run_failures_on_JAK1": 0,
            "_rule": (
                "Every figure in this entry is regenerated from those "
                "artifacts by analyze.py 2. Quote any of them ONLY with its n "
                "and this date attached -- a number that cannot be traced back "
                "to a file will drift, which is exactly what happened to the "
                "separation figure below."
            ),
        },
        "replicated_on_other_compounds": True,
        "benchmarked": True,
        "benchmark_date": "2026-08-15",
        "generalises": (
            "MEASURED over 23 pairs / 3 targets. There is no offset to "
            "generalise — the mean signed error is indistinguishable from zero, "
            "so nothing is applied as a correction anywhere in this module and "
            "nothing needs to be."
        ),
        "supersedes": (
            "The previous entry here was ONE COMPOUND against ONE literature "
            "value: tofacitinib predicted at 46.4 nM against 0.50 nM measured, "
            "reported as a 1.97-log systematic bias, with a 2.36-log separation "
            "from decoys that was 1 active against 2 decoys and had no n at "
            "all. BOTH figures are withdrawn as stated. ALSO WITHDRAWN: four "
            "later separation figures -- 12x6 -> 2.08 log / AUC 0.972, "
            "12x9 -> 2.32 / 0.981, 12x10 -> 2.36 / 0.983, 12x11 -> 2.27 / "
            "0.977. Each was a read of THIS SAME artifact taken while a repair "
            "pass was still recovering decoys that had failed to run, so each "
            "is the full actives set scored against an incomplete decoy set. "
            "A DECOY THAT FAILED TO RUN IS NOT A DECOY THAT SCORED BADLY: the "
            "six missing decoys were tautomer/CCD-matching failures, not weak "
            "binders, and dropping them shrank the effective n while "
            "flattering the AUC. 12x12 -> 2.13 / 0.958 supersedes all of them. "
            "The verdict -- triage supported -- holds under every one of the "
            "five counts."
        ),
        "result": {
            "absolute_accuracy": {
                "n_pairs": 23,
                "mean_signed_error_log": 0.32,
                "ci95_of_mean_offset": [-0.07, 0.72],
                "p_vs_zero": 0.12,
                "sign_split": "16 predicted too weak / 7 too strong",
                "mae_log": 0.82,
                "rmse_log": 1.01,
                "ground_truth_own_spread_log": 0.76,
                "_ground_truth_note": (
                    "mean ChEMBL sd over the 17 compounds with >=3 "
                    "measurements; ruxolitinib 0.72 over n=38, osimertinib 1.14 "
                    "over n=264. The model's error is now essentially "
                    "INDISTINGUISHABLE from the experimental noise of the data "
                    "scoring it, and 1.97 sits about five standard errors "
                    "outside the confidence interval."
                ),
                "tofacitinib_rescored": {
                    "measured_pchembl_consensus": 8.35,
                    "n_chembl_measurements": 64,
                    "predicted_pchembl": 7.39,
                    "signed_log_error": 0.96,
                    "_note": "the old 1.97 came from one 0.50 nM paper value",
                },
            },
            "within_target_ranking_NOT_SUPPORTED": {
                "JAK1": {"n": 12, "spearman": 0.483, "ci95": [-0.05, 0.77], "p": 0.11},
                "BCL2": {"n": 5, "spearman": 0.600, "p": 0.28},
                "EGFR": {
                    "n": 6,
                    "spearman": 0.314,
                    "p": 0.54,
                    "_provisional": (
                        "The EGFR artifact was STILL BEING REPAIRED at the "
                        "2026-08-15 18:46 read; its n is growing toward 12. "
                        "JAK1 and BCL-2 are final. Re-run analyze.py 2 before "
                        "quoting this row. The actives-vs-decoys figures are "
                        "JAK1-only and are NOT affected."
                    ),
                },
                "_summary": (
                    "Three targets, all three positive, NONE significant. "
                    "Every interval includes zero."
                ),
                "pooled_DO_NOT_USE": {
                    "n": 23,
                    "spearman": 0.564,
                    "p": 0.005,
                    "_why_not": (
                        "pooling targets with different potency offsets "
                        "manufactures rank correlation out of the offset"
                    ),
                },
                "untested_case": (
                    "Measured on DIVERSE CHEMISTRY ONLY. No congeneric series "
                    "could be assembled, because Paperclip's statement timeout "
                    "blocks the GROUP BY assay_id needed to find one. A "
                    "congeneric series is the setting chemists actually rank "
                    "in, and the one where this would most plausibly look "
                    "better -- so read this as NOT SUPPORTED AND NOT YET "
                    "TESTED WHERE IT MATTERS, not as shown to fail. The "
                    "missing series, not the missing compounds, is the real "
                    "limitation of this benchmark."
                ),
            },
            "actives_vs_decoys_CONFIRMED": {
                "n_actives": 12,
                "n_decoys": 12,
                "n_pairs": 144,
                "source_artifact": "out/claim2_JAK1.json",
                "measured_on": "2026-08-15",
                "run_failures_remaining": 0,
                "predicted_pchembl_actives": [7.07, 0.94],
                "predicted_pchembl_decoys": [4.94, 0.83],
                "binder_probability_actives": [0.724, 0.200],
                "binder_probability_decoys": [0.137, 0.075],
                "separation_log_units": 2.13,
                "roc_auc_on_affinity": 0.958,
                "roc_auc_on_binder_probability": 1.000,
                "cohens_d": 2.41,
                "_caveat": (
                    "these decoys (caffeine, metformin) are trivially easy — an "
                    "AUC of 1.0 measures binder/non-binder triage, not potency "
                    "ranking"
                ),
                "_provenance": (
                    "12 actives x 12 decoys on JAK1, from out/claim2_JAK1.json "
                    "as of 2026-08-15 with zero remaining run failures, "
                    "regenerated by analyze.py 2. QUOTE THIS ONLY WITH ITS n. "
                    "See 'supersedes' for the four void mid-repair values."
                ),
            },
        },
        "what_it_shows": (
            "The head TRIAGES. It separates binders from non-binders with a "
            "real n, and it carries no systematic offset to correct for."
        ),
        "what_it_does_not_show": (
            "That it can order actives against each other — within-target "
            "Spearman is +0.48 with a 95% CI including zero, and 'use it to "
            "rank candidates within a target' is WITHDRAWN as a recommendation. "
            "Nor that it measures potency: an 0.85-log MAE is a factor of 7, so "
            "never compare its absolute value against a nanomolar threshold."
        ),
    },
    "il17a_esmfold_dimer": {
        "what_was_measured": (
            "ESMFold's inter-chain geometry across complexes, and whether its "
            "own confidence predicts when it has failed."
        ),
        "how": (
            "14 protein-protein complexes with deposited references, each under "
            "two poly-glycine linker constructions (G25 and G50) = 28 runs, "
            "spanning obligate homo-oligomers, transient hetero-complexes, "
            "antibody-antigen and our own pipeline targets. Input sequences "
            "taken from the reference structure's observed CA residues, giving "
            "1:1 correspondence."
        ),
        "sample_size": {"complexes": 14, "constructions": 2, "runs": 28},
        "replicated_on_other_complexes": True,
        "benchmarked": True,
        "benchmark_date": "2026-08-15",
        "generalises": (
            "MEASURED over 14 complexes. What generalises is the GATE, not a "
            "verdict on the tool: pTM tracks contact recovery at rho +0.79 and "
            "pTM >= 0.80 was 5 of 5 with zero false alarms in 28 runs."
        ),
        "supersedes": (
            "The previous entry here was ONE complex and read as a property of "
            "the tool. It is reproduced EXACTLY and re-attributed: the '1 "
            "contact against 97' is an INPUT-CONSTRUCTION artifact. See "
            "input_construction_confound below."
        ),
        "result": {
            "contact_definition": "CA-CA pairs within 8.0 A between chains",
            "contact_definition_note": (
                "8DYG gives 97 CA-CA PAIRS but only 29 residues-in-contact. "
                "The 97 is a PAIR count and a residue count is not comparable "
                "to it. The harness reproduces the reference numbers exactly "
                "(8DYG: 97 contacts, minCA 4.04 A, COM 12.81 A)."
            ),
            "bimodal": {
                "median_contact_recovery": 0.42,
                "above_50pct_recovery": "6 of 14",
                "exactly_zero_recovery": "6 of 14",
                "successes": {
                    "HIVPR_dimer": 0.90,
                    "SOD1_dimer": 0.90,
                    "KRAS_RAF1RBD": 0.84,
                    "TNFa_trimer": 0.78,
                    "Barnase_barstar": 0.72,
                    "Lysozyme_Fab": 0.54,
                },
            },
            "self_report_tracks_error": {
                "ptm_vs_contact_recovery_spearman": 0.788,
                "ptm_vs_contact_recovery_ci95": [0.57, 0.91],
                "ptm_vs_complex_tm_spearman": 0.943,
                "plddt_vs_contact_recovery_spearman": 0.627,
                "n": 28,
            },
            "ptm_as_a_gate": {
                "no_cut": {"runs": 28, "median_recovery": 0.414, "zero_recovery": 10},
                "ptm_ge_0.60": {"runs": 18, "median_recovery": 0.708, "zero_recovery": 2},
                "ptm_ge_0.80": {"runs": 5, "median_recovery": 0.873, "zero_recovery": 0},
                "false_alarms_in_28_runs": 0,
                "_false_alarm_definition": "a run below pTM 0.60 that recovered >=50% of contacts — there were none",
                "one_false_confidence": (
                    "Trypsin/BPTI, pTM 0.752, recovery 0 — both chains fold "
                    "well (chainTM 0.97) but BPTI docks on the wrong face. This "
                    "is why the usable gate is 0.80 and not 0.70."
                ),
            },
            "input_construction_confound": {
                "_the_point": (
                    "Same tool, same complex, same reference, same contact set. "
                    "Only the input sequence differs."
                ),
                "IL17A_full_uniprot_mature_chain": {
                    "contacts": 1,
                    "contact_recovery": 0.000,
                    "min_inter_chain_ca_a": 7.30,
                    "com_separation_a": 21.84,
                    "ptm": 0.399,
                    "_note": "reproduces the withdrawn claim exactly",
                },
                "IL17A_crystallographically_ordered_core": {
                    "contacts": 55,
                    "contact_recovery": 0.423,
                    "min_inter_chain_ca_a": 4.25,
                    "com_separation_a": 13.51,
                    "complex_tm": 0.861,
                    "ptm": 0.684,
                },
                "not_general": (
                    "TNF-alpha is unaffected by the same swap (0.780 ordered "
                    "core vs 0.807 full mature chain). It bites on chains with "
                    "long disordered termini."
                ),
                "linker_length_is_minor": "paired G50-G25 mean difference -0.032, sd 0.109, n=14",
            },
        },
        "what_it_shows": (
            "ESMFold is BIMODAL at interfaces, not uniformly bad, and its own "
            "pTM tells you which mode you got. Trim the input to the ordered "
            "region, gate at pTM >= 0.80, and treat pTM < 0.6 as 'no answer' "
            "rather than as a negative result."
        ),
        "what_it_does_not_show": (
            "That a low-contact result is a finding about the complex. On a "
            "protein with disordered termini it is a prompt to re-run on the "
            "ordered core. ``esmfold_predict`` still returns the self-report "
            "and the contact count and does not gate, flag or downweight on "
            "this basis — the caller applies the gate."
        ),
    },
    "seed_dispersion_and_site_convergence": {
        "what_was_measured": (
            "How far apart reseeded cofolds land, and whether they land on the "
            "site that was asked about."
        ),
        "how": "Eight seeds of one probe; 24 probe runs across the seed sweep, on KRAS.",
        "sample_size": {"targets": 1, "target": "KRAS", "runs": 24, "seeds_one_probe": 8},
        "replicated_on_other_targets": False,
        "benchmarked": False,
        "generalises": (
            "The seed statistics within this target are sound; whether the "
            "magnitude transfers to other targets is UNKNOWN. This is why "
            "dispersion is RECOMPUTED on every call rather than assumed."
        ),
        "result": {
            "median_pairwise_centroid_dispersion_a": 0.21,
            "seeds_within_0.2_a": "7 of 8",
            "runs_converging_on_one_site": "21 of 24",
            "site_converged_on": "SI/II-P",
            "site_asked_about": "switch-II",
        },
        "what_it_shows": (
            "Tight seed-to-seed agreement coexisted with the runs landing on a "
            "real site that was NOT the site the question was about. So "
            "agreement between seeds is not by itself evidence that the site is "
            "the right one — which is why ``cofold_complex`` returns the "
            "contact residues and makes the caller check them."
        ),
    },
    "bioemu_frame_format": {
        "what_was_measured": "The literal content of BioEmu output frames.",
        "how": "Inspection of a 169-residue KRAS ensemble.",
        "sample_size": {"ensembles": 1, "frames": 16},
        "benchmarked": False,
        "generalises": (
            "This one is a FORMAT property of the tool's output, not a "
            "performance claim, and it is re-checkable on any run: count the "
            "atoms and read the B-factor column."
        ),
        "result": {
            "atoms_per_frame": 835,
            "residues": 169,
            "side_chains_present": False,
            "frames_pre_superposed": True,
            "max_com_spread_a": 0.045,
            "optimal_rotation_from_identity_a": 5e-8,
            "residue_indexing": "zero-indexed",
            "b_factors": "all zero — there is no per-frame confidence",
        },
        "consequence": (
            "Repack side chains before any fpocket / mdpocket run: those tools "
            "define pockets from side-chain atoms, so unrepacked frames inflate "
            "every volume. No alignment step is needed — frames arrive "
            "superposed."
        ),
    },
    "generative_ensembles_on_apo_LITERATURE": {
        "source": "external literature, NOT measured by us",
        "benchmarked_by_us": False,
        "result": {
            "cryptic_pocket_recovery_from_holo_seed": 0.86,
            "cryptic_pocket_recovery_from_apo_seed": 0.56,
            "also_reported": (
                "over-population of partially unfolded and over-extended "
                "conformations; no method reliably predicts the absolute "
                "probability that a pocket is open; all fail below 1% occupancy"
            ),
        },
        "consequence": (
            "Filter frames on radius of gyration, SASA and secondary-structure "
            "sanity before scoring, and do not report a sampled open-state "
            "population as a measurement. Apo is our normal case."
        ),
    },
}

# Contact definition reused everywhere so numbers are comparable across calls
# and comparable to the IL-17A reference figures.
_CA_CONTACT_CUTOFF_A = 8.0
_LIGAND_CONTACT_CUTOFF_A = 4.5

# BioEmu multimer linker. The tool itself rejects >1 chain, so a linker is the
# only route to a multimer ensemble; 5-10 glycines is the working range.
_BIOEMU_LINKER_MIN = 5
_BIOEMU_LINKER_MAX = 10
_BIOEMU_LINKER_DEFAULT = 8

# Rule 12: a predictor that cannot recover a known potent binder within one log
# is uninformative for absolute values on this target. This is a threshold the
# caller applies to a FRESHLY MEASURED control, not a stored correction.
_CONTROL_LOG_THRESHOLD = 1.0


# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
def _require_proto_tools() -> None:
    """Fail loudly and usefully if we are not under the proto-tools python."""
    try:
        import proto_tools  # noqa: F401, PLC0415
    except ImportError as exc:  # pragma: no cover - environment guard
        raise RuntimeError(
            "proto_tools is not importable. This module must run under the "
            "proto-tools python interpreter (the venv that has proto_tools "
            "installed); set PROTO_PY to it and invoke `$PROTO_PY predict.py`. "
            "It is a plain in-process import — there is no MCP server and no "
            "CLI to call instead."
        ) from exc


def _require_modal(device: str) -> None:
    """Check Modal credentials are reachable from the ENVIRONMENT.

    Deliberately does NOT read any dotenv file by path: this module is
    expected to run in sandboxes where no such path exists.
    """
    if device != "modal":
        return
    if os.environ.get("MODAL_TOKEN_ID") and os.environ.get("MODAL_TOKEN_SECRET"):
        return
    if (Path.home() / ".modal.toml").is_file():
        return
    raise RuntimeError(
        "device='modal' but no Modal credentials found. Set MODAL_TOKEN_ID and "
        "MODAL_TOKEN_SECRET in the environment, or set MODAL_PROFILE with a "
        "~/.modal.toml present. Credentials are read from the environment only; "
        "this module never reads a dotenv file by path. Pass device='cpu' to "
        "run locally instead (slow)."
    )


def _prov(
    tool_key: str,
    config: Any,
    *,
    device: str,
    wall_s: float,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Provenance block. Every returned dict carries one."""
    out: dict[str, Any] = {
        "tool": tool_key,
        "invocation": "proto_tools python import, in-process",
        "device": device,
        "wall_clock_s": round(wall_s, 1),
        "config": {
            k: getattr(config, k, None)
            for k in (
                "seed",
                "use_msa",
                "recycling_steps",
                "sampling_steps",
                "diffusion_samples",
                "diffusion_samples_affinity",
                "num_recycles",
                "chain_linker",
                "num_samples",
                "model_name",
                "filter_samples",
                "denoiser_type",
            )
            if hasattr(config, k)
        },
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if extra:
        out.update(extra)
    return out


# ---------------------------------------------------------------------------
# Geometry helpers (gemmi, via the Structure entity's CIF/PDB text)
# ---------------------------------------------------------------------------
def _gemmi_model(structure: Any) -> Any:
    import gemmi  # noqa: PLC0415

    st = gemmi.read_structure_string(structure.structure_cif)
    st.setup_entities()
    return st[0]


def _ca_by_chain(structure: Any) -> dict[str, list[tuple[float, float, float]]]:
    model = _gemmi_model(structure)
    out: dict[str, list[tuple[float, float, float]]] = {}
    for chain in model:
        pts = []
        for res in chain:
            atom = res.find_atom("CA", "*")
            if atom is not None:
                pts.append((atom.pos.x, atom.pos.y, atom.pos.z))
        if pts:
            out[chain.name] = pts
    return out


def _het_atoms(structure: Any) -> dict[str, list[tuple[float, float, float]]]:
    """Non-polymer (ligand) heavy atoms, keyed by chain name."""
    import gemmi  # noqa: PLC0415

    model = _gemmi_model(structure)
    out: dict[str, list[tuple[float, float, float]]] = {}
    for chain in model:
        pts = []
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and info.is_amino_acid():
                continue
            if res.name in ("HOH", "WAT"):
                continue
            for atom in res:
                if atom.element.name != "H":
                    pts.append((atom.pos.x, atom.pos.y, atom.pos.z))
        if pts:
            out[chain.name] = pts
    return out


def _centroid(points: list[tuple[float, float, float]]) -> tuple[float, float, float] | None:
    if not points:
        return None
    n = len(points)
    return (
        sum(p[0] for p in points) / n,
        sum(p[1] for p in points) / n,
        sum(p[2] for p in points) / n,
    )


def _dist(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.dist(a, b)


def _inter_chain_ca_contacts(structure: Any, cutoff: float = _CA_CONTACT_CUTOFF_A) -> dict[str, Any]:
    """Inter-chain CA contacts between the first two protein chains.

    ``contacts`` is the CA-CA PAIR count under ``cutoff``. That is the exact
    definition behind the IL-17A reference figures (1 predicted vs 97
    deposited), verified by re-running the original comparison: 8DYG gives 97
    pairs but only 29 residues-in-contact, so a residue count is NOT the
    documented number. Both are returned; ``contacts`` is the comparable one.
    """
    ca = _ca_by_chain(structure)
    names = [n for n in ca if len(ca[n]) > 1]
    if len(names) < 2:
        return {
            "n_protein_chains": len(names),
            "contacts": None,
            "residues_in_contact": None,
            "min_inter_chain_ca_a": None,
            "com_separation_a": None,
            "definition": None,
            "note": "single protein chain — no interface to measure",
        }
    a, b = ca[names[0]], ca[names[1]]
    dists = [[_dist(pa, pb) for pb in b] for pa in a]
    n_pairs = sum(1 for row in dists for d in row if d < cutoff)
    n_res = sum(1 for row in dists if min(row) < cutoff)
    ca_a, ca_b = _centroid(a), _centroid(b)
    return {
        "n_protein_chains": len(names),
        "chains_compared": [names[0], names[1]],
        "contacts": n_pairs,
        "residues_in_contact": n_res,
        "min_inter_chain_ca_a": round(min(min(row) for row in dists), 2),
        "com_separation_a": round(_dist(ca_a, ca_b), 2) if ca_a and ca_b else None,
        "definition": (
            f"CA-CA pairs within {cutoff} A between chains {names[0]} and "
            f"{names[1]}"
        ),
    }


def _ligand_contact_residues(structure: Any, cutoff: float = _LIGAND_CONTACT_CUTOFF_A) -> list[str]:
    """Protein residues (as 'CHAIN:NUM:RESNAME') within `cutoff` of any ligand heavy atom."""
    import gemmi  # noqa: PLC0415

    het = _het_atoms(structure)
    lig_pts = [p for pts in het.values() for p in pts]
    if not lig_pts:
        return []
    model = _gemmi_model(structure)
    hits: list[str] = []
    for chain in model:
        for res in chain:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_amino_acid():
                continue
            for atom in res:
                if atom.element.name == "H":
                    continue
                p = (atom.pos.x, atom.pos.y, atom.pos.z)
                if any(_dist(p, lp) < cutoff for lp in lig_pts):
                    hits.append(f"{chain.name}:{res.seqid.num}:{res.name}")
                    break
    return hits


def _res_sort_key(label: str) -> tuple[str, int]:
    chain, num, _ = label.split(":", 2)
    return (chain, int(num))


def _pairwise_stats(values: list[float]) -> dict[str, Any]:
    if not values:
        return {"n_pairs": 0, "median_a": None, "min_a": None, "max_a": None}
    return {
        "n_pairs": len(values),
        "median_a": round(statistics.median(values), 3),
        "min_a": round(min(values), 3),
        "max_a": round(max(values), 3),
    }


def _kabsch(mobile: list[tuple], target: list[tuple]) -> tuple[Any, Any, Any]:
    """Rotation + centroids superposing `mobile` onto `target`."""
    import numpy as np  # noqa: PLC0415

    p = np.asarray(mobile, dtype=float)
    q = np.asarray(target, dtype=float)
    pc, qc = p.mean(0), q.mean(0)
    h = (p - pc).T @ (q - qc)
    u, _s, vt = np.linalg.svd(h)
    d = np.diag([1.0, 1.0, float(np.sign(np.linalg.det(vt.T @ u.T)))])
    return vt.T @ d @ u.T, pc, qc


def _ligand_centroids_common_frame(structures: list[Any]) -> list[tuple[float, float, float]]:
    """Ligand centroids after superposing every structure's protein CA onto the first.

    WHY: Boltz-2 emits each sample in its own arbitrary coordinate frame. Taking
    ligand centroids straight off the raw CIFs measured 15.57 A of "dispersion"
    between two seeds whose ligand-contact residue sets were IDENTICAL — a
    number produced entirely by the frames not being aligned. Superpose first,
    always. Verified with a rigid-body control: a rotated and translated copy of
    one structure now returns 0.000 A.
    """
    import numpy as np  # noqa: PLC0415

    ref_ca = _ca_by_chain(structures[0])
    out: list[tuple[float, float, float]] = []
    for i, st in enumerate(structures):
        pts = [p for chain_pts in _het_atoms(st).values() for p in chain_pts]
        c = _centroid(pts)
        if c is None:
            continue
        if i == 0:
            out.append(c)
            continue
        ca = _ca_by_chain(st)
        shared = [k for k in ref_ca if k in ca and len(ref_ca[k]) == len(ca[k])]
        if not shared:
            continue  # cannot superpose -> refuse to emit a number
        mob = [p for k in shared for p in ca[k]]
        tgt = [p for k in shared for p in ref_ca[k]]
        rot, pc, qc = _kabsch(mob, tgt)
        moved = rot @ (np.asarray(c) - pc) + qc
        out.append((float(moved[0]), float(moved[1]), float(moved[2])))
    return out


def _load_reference(reference: Any) -> Any:
    """Accept a Structure, a file path, or raw PDB/CIF text."""
    from proto_tools.entities.structures import Structure  # noqa: PLC0415

    if hasattr(reference, "structure_cif"):
        return reference
    text = str(reference)
    if "\n" not in text and Path(text).is_file():
        return Structure.from_file(text)
    return Structure(structure=text)


def _cofold_control(structures: list[Any], reference: Any) -> dict[str, Any]:
    """CA RMSD of each cofold seed against a supplied reference structure.

    Fills the dossier's ``structure.cofold_control`` block. Scores the cofold
    against the crystal so a reader can see whether cofolding reproduces a
    KNOWN answer for THIS target — the same discipline as the affinity
    positive control, and a per-target measurement rather than a stored claim.

    Residue matching is deliberately strict: if the reference and the
    prediction do not have the same number of CA atoms in the compared chain,
    NO number is emitted and the reason is returned instead. A silently
    mis-paired RMSD is worse than a null.
    """
    import numpy as np  # noqa: PLC0415

    try:
        ref = _load_reference(reference)
    except Exception as exc:  # noqa: BLE001
        return {"reference_loaded": False, "cofold_rmsd_a": None, "reason": f"could not load reference: {exc}"}

    ref_ca = _ca_by_chain(ref)
    if not ref_ca:
        return {"reference_loaded": False, "cofold_rmsd_a": None, "reason": "reference has no CA atoms"}
    ref_chain = max(ref_ca, key=lambda k: len(ref_ca[k]))
    ref_pts = ref_ca[ref_chain]

    rmsds: list[float] = []
    reasons: list[str] = []
    for st in structures:
        ca = _ca_by_chain(st)
        cand = [k for k in ca if len(ca[k]) == len(ref_pts)]
        if not cand:
            reasons.append(
                f"no chain with {len(ref_pts)} CA atoms to match reference chain "
                f"{ref_chain} (prediction has {[(k, len(v)) for k, v in ca.items()]})"
            )
            continue
        mob = ca[cand[0]]
        rot, pc, qc = _kabsch(mob, ref_pts)
        moved = (rot @ (np.asarray(mob) - pc).T).T + qc
        diff = moved - np.asarray(ref_pts)
        rmsds.append(float(np.sqrt((diff**2).sum(1).mean())))

    if not rmsds:
        return {
            "reference_loaded": True,
            "reference_chain": ref_chain,
            "reference_ca_count": len(ref_pts),
            "cofold_rmsd_a": None,
            "reason": (
                "residue counts do not match, so no RMSD was computed. Trim "
                "the reference to the modelled residues and retry. Details: "
                + "; ".join(reasons[:2])
            ),
        }
    return {
        "reference_loaded": True,
        "reference_chain": ref_chain,
        "reference_ca_count": len(ref_pts),
        "cofold_rmsd_a": round(statistics.median(rmsds), 3),
        "cofold_rmsd_a_per_seed": [round(r, 3) for r in rmsds],
        "method": "Kabsch superposition on all CA atoms, 1:1 by index",
        "reproduces_reference": None,
        "trusted": None,
        "_why_null": (
            "reproduces_reference and trusted are judgements, not "
            "measurements, and this module does not make them. There is no "
            "calibrated RMSD threshold here — one would have to come from a "
            "benchmark across targets, and we do not have one. Report the "
            "RMSD and let the reader weigh it."
        ),
    }


def _seed_dispersion(structures: list[Any]) -> dict[str, Any]:
    """Dispersion ACROSS SEEDS — recomputed fresh on every call.

    This is a per-run uncertainty measure, not a stored claim about the method.
    Small dispersion is NOT evidence of correctness: on KRAS, 24 runs at a
    median 0.21 A dispersion converged on a real site that was not the site the
    question was about. Read ``converged_site`` alongside this.
    """
    lig_centroids = _ligand_centroids_common_frame(structures) if structures else []
    lig_pairs = [_dist(a, b) for a, b in combinations(lig_centroids, 2)]

    bb_pairs: list[float] = []
    for a, b in combinations(structures, 2):
        try:
            bb_pairs.append(float(a.backbone_rmsd(b)))
        except Exception:  # noqa: BLE001 - RMSD needs matched chains; report absence, don't crash
            pass

    return {
        "n_seeds": len(structures),
        "ligand_centroid_dispersion_a": {
            **_pairwise_stats(lig_pairs),
            "frame": "protein-CA superposed onto seed 0 before measuring",
            "_why": (
                "Measured off raw un-superposed CIFs this returned 15.57 A for "
                "two seeds whose contact-residue sets were identical. Boltz-2 "
                "emits every sample in its own frame."
            ),
        },
        "backbone_ca_rmsd_across_seeds_a": _pairwise_stats(bb_pairs),
        "interpretation": (
            "Dispersion across seeds is a per-run uncertainty measure and it is "
            "measured, not assumed. It does NOT measure correctness — see "
            "converged_site, and OBSERVATIONS["
            "'seed_dispersion_and_site_convergence'] for the case where tight "
            "agreement accompanied the wrong site."
        ),
    }


def _converged_site(structures: list[Any]) -> dict[str, Any]:
    """Which site did the seeds land on, and did they agree?

    Per-run, and the reason multi-seed is not optional: a caller must be able to
    see that the model answered a different question from the one asked.
    """
    per_seed = [set(_ligand_contact_residues(st)) for st in structures]
    per_seed = [s for s in per_seed if s]
    if not per_seed:
        return {
            "has_ligand": False,
            "consensus_contact_residues": None,
            "seed_agreement_fraction": None,
            "caution": "no ligand in the complex — there is no site to converge on",
        }
    counts: dict[str, int] = {}
    for s in per_seed:
        for r in s:
            counts[r] = counts.get(r, 0) + 1
    n = len(per_seed)
    consensus = sorted([r for r, c in counts.items() if c == n], key=_res_sort_key)
    union = sorted(counts, key=_res_sort_key)
    return {
        "has_ligand": True,
        "n_seeds_with_ligand": n,
        "consensus_contact_residues": consensus,
        "union_contact_residues": union,
        "seed_agreement_fraction": round(len(consensus) / len(union), 3) if union else None,
        "contact_definition": (
            f"protein residue with any heavy atom within {_LIGAND_CONTACT_CUTOFF_A} A "
            f"of a ligand heavy atom"
        ),
        "caution": (
            "THIS IS THE SITE THE MODEL CHOSE, NOT NECESSARILY THE SITE YOU "
            "ASKED ABOUT. On KRAS, 21 of 24 runs converged on SI/II-P when the "
            "question was switch-II — a real site, and the wrong one. Check "
            "these residues against your intended site before using anything "
            "downstream of this call. High seed agreement does not settle it."
        ),
    }


def _metrics_dict(structure: Any) -> dict[str, Any]:
    m = structure.metrics
    d = m.model_dump() if hasattr(m, "model_dump") else dict(m)
    return {k: v for k, v in d.items() if not isinstance(v, list)}


def _as_chain_list(sequences: str | list[str]) -> list[str]:
    return [sequences] if isinstance(sequences, str) else list(sequences)


# ---------------------------------------------------------------------------
# 1. Boltz-2 structure prediction  (proto-tools key: boltz2-prediction)
# ---------------------------------------------------------------------------
def cofold_complex(
    sequences: str | list[str],
    ligand_smiles: str | list[str] | None = None,
    *,
    n_seeds: int = 3,
    reference_structure: Any = None,
    seed: int = 42,
    device: str = "modal",
    use_msa: bool = True,
    recycling_steps: int = 3,
    sampling_steps: int = 200,
    diffusion_samples: int = 1,
    timeout: int = 3600,
    verbose: bool = False,
) -> dict[str, Any]:
    """Fold a protein, a multimer, or a protein-ligand complex with Boltz-2.

    Multi-seed by default. The confidence numbers are returned under
    ``structural_confidence`` and labelled for what they are — metrics about
    the geometry the model drew. They are not gated, corrected or downweighted
    here; ``single_target_observations`` carries the KRAS sealed-pocket result
    with its sample size so a reader can weigh it.

    Args:
        sequences: one protein sequence, or a list for a multimer.
        ligand_smiles: optional ligand SMILES (str or list) to cofold in.
        n_seeds: independent seeds. Seeds are ``seed, seed+1, ...`` — this is
            how the tool itself distinguishes duplicate complexes.

    Returns:
        dict with ``structural_confidence``, ``seed_dispersion``,
        ``converged_site``, ``interface``, ``structures_cif`` and
        ``provenance``.
    """
    _require_proto_tools()
    _require_modal(device)
    if n_seeds < 1:
        raise ValueError("n_seeds must be >= 1")

    from proto_tools.tools.structure_prediction.boltz2.boltz2 import (  # noqa: PLC0415
        Boltz2Config,
        Boltz2Input,
        run_boltz2,
    )

    protein_chains = _as_chain_list(sequences)
    chains: list[dict[str, Any]] = [{"sequence": s, "entity_type": "protein"} for s in protein_chains]
    ligs = [] if ligand_smiles is None else _as_chain_list(ligand_smiles)
    chains += [{"smiles": s, "entity_type": "ligand"} for s in ligs]

    inputs = Boltz2Input(complexes=[{"chains": chains}] * n_seeds)
    config = Boltz2Config(
        device=device,
        seed=seed,
        use_msa=use_msa,
        recycling_steps=recycling_steps,
        sampling_steps=sampling_steps,
        diffusion_samples=diffusion_samples,
        timeout=timeout,
        verbose=verbose,
    )

    t0 = time.time()
    result = run_boltz2(inputs, config)
    wall = time.time() - t0

    per_seed = [_metrics_dict(st) for st in result.structures]

    def _spread(key: str) -> dict[str, Any]:
        vals = [m[key] for m in per_seed if m.get(key) is not None]
        if not vals:
            return {"values": [], "min": None, "max": None}
        return {"values": [round(v, 4) for v in vals], "min": round(min(vals), 4), "max": round(max(vals), 4)}

    payload: dict[str, Any] = {
        "tool": "boltz2-prediction",
        "primary_output": "structures_cif",
        "n_protein_chains": len(protein_chains),
        "n_ligands": len(ligs),
        "is_multimer": len(protein_chains) > 1,
        "seeds": [seed + i for i in range(n_seeds)],
        "structural_confidence": {
            "_what_this_is": (
                "Structural-confidence metrics: how sure the model is about the "
                "geometry it drew. pTM, ipTM and pLDDT are defined as measures "
                "of geometric confidence — they are not binding scores, pocket "
                "scores or druggability scores, and nothing downstream should "
                "read them as such."
            ),
            "per_seed": per_seed,
            "confidence_score": _spread("confidence_score"),
            "complex_plddt": _spread("complex_plddt"),
            "ptm": _spread("ptm"),
            "iptm": _spread("iptm"),
            "ligand_iptm": _spread("ligand_iptm"),
            "avg_pae": _spread("avg_pae"),
        },
        "single_target_observations": {
            "_read_this_first": (
                "Benchmarked 2026-08-15 over 5 targets / 30 folds, and still "
                "NOT applied to anything above. The pLDDT family notices a "
                "sealed pocket on 5 of 5, but usually by 0.02-0.05, which no "
                "reader will see without the wild-type control beside it — and "
                "ligand_iptm ROSE on TNF-alpha. Weigh it yourself."
            ),
            "sealed_pocket_benchmark": OBSERVATIONS["kras_sealed_pocket_confidence"],
        },
        "seed_dispersion": _seed_dispersion(result.structures),
        "converged_site": _converged_site(result.structures),
        "cofold_control": _cofold_control(result.structures, reference_structure)
        if reference_structure is not None
        else {
            "reference_loaded": False,
            "cofold_rmsd_a": None,
            "reason": (
                "no reference_structure supplied. When a crystal structure of "
                "this target exists, pass it — scoring the cofold against the "
                "known answer is the only per-target check available here."
            ),
        },
        "interface": _inter_chain_ca_contacts(result.structures[0]),
        "structures_cif": [st.structure_cif for st in result.structures],
        "provenance": _prov("boltz2-prediction", config, device=device, wall_s=wall),
    }
    if n_seeds == 1:
        payload["seed_dispersion"]["warning"] = (
            "n_seeds=1: no dispersion was measured. The uncertainty on this "
            "structure is UNKNOWN, not zero."
        )
    return payload


# ---------------------------------------------------------------------------
# 2. Boltz-2 affinity  (proto-tools key: boltz2-affinity)
# ---------------------------------------------------------------------------
def _log10_um_from_nm(nm: float) -> float:
    return math.log10(nm / 1000.0)


def _exception_chain_text(exc: BaseException, limit: int = 8) -> str:
    """Full text of an exception and everything it was raised from.

    Needed because proto-tools re-raises the interesting error as a pydantic
    validation TypeError whose message TRUNCATES the payload — the actual cause
    string ('_unphysical.xtc') only survives further down the __cause__ /
    __context__ chain.
    """
    parts: list[str] = []
    seen: set[int] = set()
    cur: BaseException | None = exc
    while cur is not None and len(parts) < limit and id(cur) not in seen:
        seen.add(id(cur))
        parts.append(f"{type(cur).__name__}: {cur}")
        cur = cur.__cause__ or cur.__context__
    return "\n".join(parts)


def cofold_affinity(
    protein: str | list[str],
    ligand_smiles: str | list[str],
    *,
    ligand_names: list[str] | None = None,
    positive_control_smiles: str | None = None,
    positive_control_name: str | None = None,
    positive_control_measured_nm: float | None = None,
    seed: int = 42,
    device: str = "modal",
    use_msa: bool = True,
    diffusion_samples_affinity: int = 5,
    timeout: int = 3600,
    verbose: bool = False,
) -> dict[str, Any]:
    """Triage binders from non-binders against a protein with the Boltz-2
    affinity head.

    The primary output is ``ranking``, and after the 2026-08-15 benchmark it
    should be read as a **binder/non-binder split, not a series ordering**:
    within-target Spearman is +0.483 with a 95% CI of (-0.05, +0.77), p=0.11
    (n=12 on JAK1) and is not significant on any of 3 targets, while
    actives-vs-decoys separation is ROC AUC 0.958 on affinity and 1.000 on
    binder probability, 2.13 log units, on 12 actives x 12 decoys (JAK1,
    out/claim2_JAK1.json, 2026-08-15).

    Absolute values are returned under ``absolute`` marked ``is_a_kd: False``
    and ``benchmarked_against_measured_affinities: False`` — meaning THIS RUN
    was not calibrated against measurements, not that the head is unstudied.
    Benchmarked over 23 pairs it carries no systematic offset (mean signed error
    +0.32 log, p=0.12) but an MAE of 0.82 log, a factor of 6.6, so nothing it
    emits is quotable as a potency. No correction, offset or calibration is
    applied to any returned value.

    Supply ``positive_control_smiles`` + ``positive_control_measured_nm`` to run
    the rule-12 control in the same call. That control is a FRESH per-run
    measurement on your own target, and it is the number to act on.
    """
    _require_proto_tools()
    _require_modal(device)

    from proto_tools.tools.structure_prediction.boltz2.boltz2_affinity import (  # noqa: PLC0415
        Boltz2AffinityConfig,
        Boltz2AffinityInput,
        run_boltz2_affinity,
    )

    ligs = _as_chain_list(ligand_smiles)
    names = list(ligand_names) if ligand_names else [f"ligand_{i}" for i in range(len(ligs))]
    if len(names) != len(ligs):
        raise ValueError("ligand_names must be the same length as ligand_smiles")

    control_idx: int | None = None
    if positive_control_smiles is not None:
        control_idx = 0
        ligs = [positive_control_smiles, *ligs]
        names = [positive_control_name or "positive_control", *names]

    protein_chains = [{"sequence": s, "entity_type": "protein"} for s in _as_chain_list(protein)]
    inputs = Boltz2AffinityInput(
        complexes=[
            {"chains": [*protein_chains, {"smiles": smi, "entity_type": "ligand"}]} for smi in ligs
        ]
    )
    config = Boltz2AffinityConfig(
        device=device,
        seed=seed,
        use_msa=use_msa,
        diffusion_samples_affinity=diffusion_samples_affinity,
        timeout=timeout,
        verbose=verbose,
    )

    t0 = time.time()
    result = run_boltz2_affinity(inputs, config)
    wall = time.time() - t0

    rows: list[dict[str, Any]] = []
    for i, (name, smi, st) in enumerate(zip(names, ligs, result.structures, strict=True)):
        m = _metrics_dict(st)
        v = m.get("affinity_pred_value")
        rows.append(
            {
                "name": name,
                "smiles": smi,
                "is_positive_control": i == control_idx,
                "_sort_key": v,
                "binder_probability": m.get("affinity_probability_binary"),
                "absolute": {
                    "affinity_pred_value_log10_ic50_um": v,
                    "unit": "log10(IC50 in micromolar), lower = predicted stronger",
                    "is_a_kd": False,
                    "is_a_potency_measurement": False,
                    "benchmarked_against_measured_affinities": False,
                    "correction_applied": None,
                    "warning": (
                        "Do NOT report this as a Kd, IC50 or potency and do NOT "
                        "compare it against a nanomolar threshold. Benchmarked "
                        "2026-08-15 over 23 pairs the head carries NO "
                        "systematic offset (mean signed error +0.32 log, 95% CI "
                        "-0.07 to +0.72, p=0.12) but an MAE of 0.82 log — a "
                        "factor of 6.6 — so there is nothing to correct and "
                        "nothing quotable. Use the ranking as a binder/"
                        "non-binder split (NOT as an ordering of actives: "
                        "within-target Spearman +0.48, p=0.11), and use the "
                        "positive control you ran on THIS target."
                    ),
                },
                "raw_metrics": m,
            }
        )

    ranked = sorted([r for r in rows if r["_sort_key"] is not None], key=lambda r: r["_sort_key"])
    for rank, r in enumerate(ranked, start=1):
        r["rank"] = rank
    best = ranked[0]["_sort_key"] if ranked else None
    for r in ranked:
        r["relative_score_log_units_vs_best"] = round(r["_sort_key"] - best, 4) if best is not None else None
    for r in rows:
        r.pop("_sort_key", None)

    span = ranked[-1]["relative_score_log_units_vs_best"] if len(ranked) >= 2 else None

    control: dict[str, Any] = {
        "run": False,
        "reliable": None,
        "note": (
            "Rule 12: a prediction without its control is not a measurement. "
            "Pass positive_control_smiles and positive_control_measured_nm to "
            "measure, on THIS target, how far the predictor lands from a known "
            "binder. That per-run number is worth more than any stored one."
        ),
    }
    if control_idx is not None and positive_control_measured_nm is not None:
        crow = next(r for r in rows if r["is_positive_control"])
        pred = crow["absolute"]["affinity_pred_value_log10_ic50_um"]
        if pred is not None:
            err = pred - _log10_um_from_nm(positive_control_measured_nm)
            control = {
                "run": True,
                "ligand": crow["name"],
                "measured_nm": positive_control_measured_nm,
                "predicted_log10_ic50_um": pred,
                "log_error": round(err, 3),
                "log_error_direction": "predicted weaker than measured"
                if err > 0
                else "predicted stronger than measured",
                "threshold_log_units": _CONTROL_LOG_THRESHOLD,
                "reliable": bool(abs(err) <= _CONTROL_LOG_THRESHOLD),
                "scope": (
                    "ONE compound on THIS target, measured in this run. It is a "
                    "check, not a calibration: it is not applied to the other "
                    "ligands' values and must not be used as an offset."
                ),
                "note": (
                    "reliable=False means the predictor did not recover this "
                    "known binder within one log here, so its absolute values "
                    "for novel chemotypes on this target are uninformative. The "
                    "RANKING may still be usable — read ranking_span_log_units "
                    "and whether the control ranked where you expected."
                ),
            }

    return {
        "tool": "boltz2-affinity",
        "primary_output": "ranking",
        "ranking": [
            {
                "rank": r["rank"],
                "name": r["name"],
                "relative_score_log_units_vs_best": r["relative_score_log_units_vs_best"],
                "binder_probability": r["binder_probability"],
                "is_positive_control": r["is_positive_control"],
            }
            for r in ranked
        ],
        "ranking_span_log_units": span,
        "ranking_usable": len(ranked) >= 2,
        "ranking_caveat": (
            "A ranking of one ligand is not a ranking. Pass more ligands, or a "
            "positive control, so there is something to rank against."
        )
        if len(ranked) < 2
        else (
            "Ranks are within-target only — do not compare ranks or scores "
            "across different protein targets; pooling targets with different "
            "potency offsets manufactures rank correlation out of the offset. "
            "And read this ordering as a BINDER/NON-BINDER split, not as an "
            "ordering of actives: benchmarked 2026-08-15, within-target "
            "Spearman is +0.483, 95% CI (-0.05, +0.77), p=0.11 (n=12, JAK1) "
            "and not significant on any of 3 targets, while actives-vs-decoys "
            "is ROC AUC 0.958 / 1.000 at 2.13 log on 12 actives x 12 decoys."
        ),
        "control": control,
        "per_ligand": rows,
        "single_compound_observations": {
            "_read_this_first": (
                "Benchmarked 2026-08-15 over 23 pairs / 3 targets, and still "
                "NOT applied to anything above. Read it for what it withdraws: "
                "there is no 1.97-log bias to correct for, and ordering actives "
                "against each other is NOT supported (Spearman +0.48, "
                "p=0.11). The ranks above are a binder/non-binder split."
            ),
            "affinity_benchmark": OBSERVATIONS["tofacitinib_affinity_error"],
        },
        "provenance": _prov("boltz2-affinity", config, device=device, wall_s=wall),
    }


# ---------------------------------------------------------------------------
# 3. ESMFold  (proto-tools key: esmfold-prediction)
# ---------------------------------------------------------------------------
def esmfold_predict(
    sequence: str | list[str],
    *,
    seed: int = 42,
    device: str = "modal",
    num_recycles: int = 4,
    chain_linker: str | None = None,
    timeout: int | None = None,
    verbose: bool = False,
) -> dict[str, Any]:
    """Fold with ESMFold and return its own confidence alongside the geometry.

    ``self_report`` carries the model's pTM, average PAE and average pLDDT.
    ``interface`` carries the measured inter-chain CA contact count, the
    residues in contact, the closest inter-chain CA approach and the
    centre-of-mass separation. Both are per-run measurements. This function
    does not gate, flag or downweight its output — the caller judges.
    """
    _require_proto_tools()
    _require_modal(device)

    from proto_tools.tools.structure_prediction.esmfold.esmfold import (  # noqa: PLC0415
        ESMFoldConfig,
        ESMFoldInput,
        run_esmfold,
    )

    chains = _as_chain_list(sequence)
    kwargs: dict[str, Any] = {
        "device": device,
        "seed": seed,
        "num_recycles": num_recycles,
        "verbose": verbose,
    }
    if chain_linker is not None:
        kwargs["chain_linker"] = chain_linker
    if timeout is not None:
        kwargs["timeout"] = timeout
    config = ESMFoldConfig(**kwargs)

    inputs = ESMFoldInput(
        complexes=[{"chains": [{"sequence": s, "entity_type": "protein"} for s in chains]}]
    )

    t0 = time.time()
    result = run_esmfold(inputs, config)
    wall = time.time() - t0

    st = result.structures[0]
    m = _metrics_dict(st)
    is_multimer = len(chains) > 1

    return {
        "tool": "esmfold-prediction",
        "primary_output": "structure_cif",
        "n_chains": len(chains),
        "is_multimer": is_multimer,
        "structure_cif": st.structure_cif,
        "self_report": {
            "avg_plddt": m.get("avg_plddt"),
            "ptm": m.get("ptm"),
            "avg_pae": m.get("avg_pae"),
            "_what_this_is": (
                "The model's own confidence for THIS run. Returned unmodified "
                "and ungated — it is a real per-run signal and the caller "
                "judges it."
            ),
            "for_scale_il17a_case": {
                "monomer": {"ptm": 0.905, "avg_pae": 3.66},
                "dimer": {"ptm": 0.399, "avg_pae": 18.26},
                "caveat": "one complex, for scale only — not a threshold",
            },
        },
        "interface": _inter_chain_ca_contacts(st)
        if is_multimer
        else {"n_protein_chains": 1, "note": "single chain — no interface to measure"},
        "multimer_note": (
            "ESMFold folds multiple chains by joining them with an internal "
            f"linker ({len(config.chain_linker)} glycines by default), which it "
            "strips from the output. Read the contact count and the model's own "
            "pTM/PAE above before relying on the inter-chain geometry."
        )
        if is_multimer
        else None,
        "single_complex_observations": {
            "_read_this_first": (
                "Benchmarked 2026-08-15 over 14 complexes / 28 runs, and still "
                "NOT applied to anything above — the gate is the caller's to "
                "apply. Two things to take from it: pTM >= 0.80 was 5 of 5 with "
                "zero false alarms, and a near-zero contact count on a chain "
                "with disordered termini is a prompt to re-run on the "
                "crystallographically ordered core, not a finding."
            ),
            "esmfold_interface_benchmark": OBSERVATIONS["il17a_esmfold_dimer"],
        },
        "provenance": _prov(
            "esmfold-prediction",
            config,
            device=device,
            wall_s=wall,
            extra={"chain_linker_used": config.chain_linker if is_multimer else None},
        ),
    }


# ---------------------------------------------------------------------------
# 4. BioEmu  (proto-tools key: bioemu-sample)
# ---------------------------------------------------------------------------
def bioemu_ensemble(
    sequence: str | list[str],
    n_samples: int = 32,
    *,
    linker_length: int = _BIOEMU_LINKER_DEFAULT,
    seed: int = 42,
    device: str = "modal",
    model_name: str = "bioemu-v1.1",
    filter_samples: bool = True,
    denoiser_type: str = "dpm",
    timeout: int = 3600,
    verbose: bool = False,
) -> dict[str, Any]:
    """Sample a conformational ensemble with BioEmu.

    BioEmu's own validator rejects any complex with more than one chain, so a
    multimer is only reachable by concatenating the chains through a
    poly-glycine linker. When that happens the returned dict carries
    ``linker.inserted = True`` with the exact residue range of every linker and
    of every original chain, because a linker changes what the ensemble means.

    Frames are backbone + C-beta only, zero-indexed, all B-factors zero — see
    ``frame_caveats``, and check them against the run: they are re-verifiable.
    """
    _require_proto_tools()
    _require_modal(device)
    if not (_BIOEMU_LINKER_MIN <= linker_length <= _BIOEMU_LINKER_MAX):
        raise ValueError(
            f"linker_length must be {_BIOEMU_LINKER_MIN}-{_BIOEMU_LINKER_MAX} "
            f"glycines (got {linker_length})"
        )

    from proto_tools.tools.structure_dynamics.bioemu.bioemu_sample import (  # noqa: PLC0415
        BioEmuConfig,
        BioEmuInput,
        run_bioemu,
    )

    chains = _as_chain_list(sequence)
    linker = "G" * linker_length

    if len(chains) > 1:
        folded = linker.join(chains)
        segments: list[dict[str, Any]] = []
        linkers: list[dict[str, Any]] = []
        pos = 0
        for i, ch in enumerate(chains):
            segments.append(
                {
                    "chain_index": i,
                    "length": len(ch),
                    "residue_range_0indexed": [pos, pos + len(ch) - 1],
                }
            )
            pos += len(ch)
            if i < len(chains) - 1:
                linkers.append(
                    {
                        "after_chain_index": i,
                        "length": linker_length,
                        "sequence": linker,
                        "residue_range_0indexed": [pos, pos + linker_length - 1],
                    }
                )
                pos += linker_length
        linker_record: dict[str, Any] = {
            "inserted": True,
            "reason": (
                "BioEmu's input validator rejects any complex with more than "
                "one chain ('BioEmu only supports single-chain proteins "
                "(monomers)'). Concatenating through a poly-glycine linker is "
                "the only route to a multimer ensemble."
            ),
            "linker_sequence": linker,
            "linker_length": linker_length,
            "n_linkers": len(linkers),
            "linkers": linkers,
            "chain_segments": segments,
            "folded_length": len(folded),
            "what_it_changes": (
                "The ensemble is of a COVALENTLY TETHERED construct, not of the "
                "biological assembly. Inter-chain distances are constrained by "
                "the tether, the relative-orientation distribution is not the "
                "free one, and the linker residues are not part of the protein. "
                "Strip the linker residue ranges above before any pocket "
                "detection or RMSD, and do not report an inter-chain "
                "measurement off these frames as if it were free-solution."
            ),
        }
    else:
        folded = chains[0]
        linker_record = {
            "inserted": False,
            "reason": "single chain — no linker needed",
            "linker_sequence": None,
            "linker_length": None,
            "n_linkers": 0,
            "linkers": [],
            "chain_segments": [
                {"chain_index": 0, "length": len(folded), "residue_range_0indexed": [0, len(folded) - 1]}
            ],
            "folded_length": len(folded),
        }

    config = BioEmuConfig(
        device=device,
        num_samples=n_samples,
        seed=seed,
        model_name=model_name,
        filter_samples=filter_samples,
        denoiser_type=denoiser_type,
        timeout=timeout,
        verbose=verbose,
    )
    inputs = BioEmuInput(complexes=[{"chains": [{"sequence": folded, "entity_type": "protein"}]}])

    # UPSTREAM BUG, reproduced: when the physical-sanity filter actually rejects
    # frames, BioEmu writes a `*_unphysical.xtc` alongside the kept frames and
    # then dies parsing its own filename —
    #   ValueError: Invalid suffix '_unphysical.xtc'
    # It surfaces as a BioEmuOutput validation error ("ensembles Field
    # required"). Hit reproducibly on a glycine-linked 2x60 construct with
    # filter_samples=True; the same call with filter_samples=False succeeds.
    # A linked construct is exactly the input most likely to produce rejectable
    # frames, so the multimer path walks into it. Retry once unfiltered and say
    # so loudly — this CHANGES WHAT THE ENSEMBLE IS.
    filter_fallback: dict[str, Any] | None = None
    t0 = time.time()
    try:
        result = run_bioemu(inputs, config)
    except Exception as exc:  # noqa: BLE001 - narrowed by the message check below
        chain = _exception_chain_text(exc)
        # The proto-tools wrapper re-raises as a BioEmuOutput validation error
        # whose message truncates the payload, so match on the whole chain and
        # fall back to the output-shape signature.
        looks_like_filter_crash = "_unphysical" in chain or (
            "BioEmuOutput" in chain and "ensembles" in chain and "Field required" in chain
        )
        if not (filter_samples and looks_like_filter_crash):
            raise
        config = BioEmuConfig(
            device=device,
            num_samples=n_samples,
            seed=seed,
            model_name=model_name,
            filter_samples=False,
            denoiser_type=denoiser_type,
            timeout=timeout,
            verbose=verbose,
        )
        result = run_bioemu(inputs, config)
        filter_fallback = {
            "triggered": True,
            "requested_filter_samples": True,
            "actual_filter_samples": False,
            "upstream_error": chain[:600],
            "matched_on": "_unphysical" if "_unphysical" in chain else "BioEmuOutput-shape",
            "what_this_means": (
                "The physical-sanity filter REJECTED frames, and upstream then "
                "crashed parsing the '_unphysical.xtc' file it had just "
                "written. The rerun kept ALL frames, including the ones the "
                "filter wanted to discard. These frames are NOT sanity-checked: "
                "steric clashes and chain breaks may be present. Filter them "
                "yourself on radius of gyration, SASA and secondary-structure "
                "sanity before scoring anything."
            ),
        }
    wall = time.time() - t0

    frames = result.ensembles[0].structures

    rgs: list[float] = []
    for f in frames:
        try:
            rgs.append(float(f.gyration_radius()))
        except Exception:  # noqa: BLE001 - a bad frame is a finding, not a crash
            pass
    pair_rmsd: list[float] = []
    for a, b in combinations(frames[: min(len(frames), 8)], 2):
        try:
            pair_rmsd.append(float(a.backbone_rmsd(b)))
        except Exception:  # noqa: BLE001
            pass

    # Re-verify the frame format on THIS run rather than asserting it.
    first_atoms = None
    if frames:
        try:
            first_atoms = sum(
                1 for line in frames[0].structure_pdb.splitlines() if line.startswith(("ATOM", "HETATM"))
            )
        except Exception:  # noqa: BLE001
            pass

    return {
        "tool": "bioemu-sample",
        "primary_output": "frames_pdb",
        "n_frames_returned": len(frames),
        "n_samples_requested": n_samples,
        "frames_filtered_out": n_samples - len(frames),
        "sequence_folded": folded,
        "linker": linker_record,
        "filter_fallback": filter_fallback
        or {"triggered": False, "actual_filter_samples": filter_samples},
        "ensemble_spread": {
            "radius_of_gyration_a": {
                "min": round(min(rgs), 3) if rgs else None,
                "max": round(max(rgs), 3) if rgs else None,
                "median": round(statistics.median(rgs), 3) if rgs else None,
            },
            "pairwise_backbone_rmsd_a": _pairwise_stats(pair_rmsd),
            "pairwise_rmsd_note": "computed on the first 8 frames only (O(n^2))",
        },
        "frame_caveats": {
            "atoms_in_first_frame_THIS_RUN": first_atoms,
            "residues_folded": len(folded),
            "side_chains_present": False,
            "must_repack_before_pocket_detection": True,
            "repack_reason": (
                "Frames are backbone + C-beta only. fpocket and mdpocket define "
                "pockets from side-chain atoms, so scoring these frames "
                "unrepacked inflates every volume. Divide "
                "atoms_in_first_frame_THIS_RUN by residues_folded to confirm — "
                "roughly 5 atoms per residue means no side chains."
            ),
            "residue_indexing": "zero-indexed",
            "per_frame_confidence_available": False,
            "b_factors": "all zero — there is no per-frame confidence to read",
            "pre_superposed": True,
            "alignment_needed_downstream": False,
            "filter_before_scoring_on": ["radius_of_gyration", "SASA", "secondary_structure_sanity"],
        },
        "frames_pdb": [f.structure_pdb for f in frames],
        "observations": {
            "frame_format": OBSERVATIONS["bioemu_frame_format"],
            "apo_degradation_literature": OBSERVATIONS["generative_ensembles_on_apo_LITERATURE"],
        },
        "provenance": _prov("bioemu-sample", config, device=device, wall_s=wall),
    }
