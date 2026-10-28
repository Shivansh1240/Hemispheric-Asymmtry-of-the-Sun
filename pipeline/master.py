"""master.py -- one entry point for every analysis and every figure of the hemispheric-asymmetry project.

The verified scripts stay where they are (one source of truth); this file knows all of them, their order, their inputs
and outputs, and runs any subset in dependency order.

    python master.py list [STAGE]          every step: id, cost, script, outputs, what it is
    python master.py inputs                the external data files (and the catalogue checksum)
    python master.py status                which outputs exist and when they were made
    python master.py order TARGET          the steps TARGET would run, in order
    python master.py docs                  markdown table of every step (README)
    python master.py run TARGET [options]  run TARGET and everything it needs
        TARGET   all | a stage (reanalysis, revision, phase4, figures, paperfigs) | a step id | several ids
        --heavy          also re-run heavy steps (> 10 min); by default a heavy step whose outputs exist is reused
        --legacy         include the legacy / superseded steps (old catalogue, old figure names)
        --verify         snapshot each step's outputs before running it and compare them numerically afterwards
        --dry-run        print what would run
        --only           run only the named steps, not their dependencies

PROJECT (the 'Asymmetry analysis' folder) is found as described below; ASYM_PROJECT overrides it.
Logs go to pipeline/logs/, verification reports to pipeline/logs/verify_*.json.
"""
import argparse, glob, hashlib, json, os, shutil, subprocess, sys, tempfile, time
from dataclasses import dataclass, field

# the 'Asymmetry analysis' folder: $ASYM_PROJECT, else this repository if it holds paper/, else the original location
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT = os.environ.get("ASYM_PROJECT") or (_REPO if os.path.isdir(os.path.join(_REPO, "paper")) else
          r"C:\Users\Shivansh Mishra\Downloads\Quantum Entanglement and Computation - Online LaTeX Editor Overleaf_files"
          r"\Asymmetry analysis")
PAPER = os.path.join(PROJECT, "paper")
HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = os.path.join(HERE, "logs")
CR, REV, P4 = "corrected_reanalysis", "revision", "revision/p4"

# ------------------------------------------------------------------------------------------------ external inputs ----
INPUTS = {
    "RGO_SOON_corrected.txt (sole sunspot source)": (os.path.join(PROJECT, "RGO_SOON_corrected.txt"), "b78c7b08f5465a3d86ae677efcd028b0"),
    "e_PField_MWO_WSO.csv (MWO faculae, WSO; with SEM)": (os.path.join(PROJECT, "e_PField_MWO_WSO.csv"), None),
    "WSO daily filtered polar field": (os.path.join(PROJECT, "Dataset Comparison", "new research", "P3_meridional_flow", "data",
                                                    "Polar Filed data of the sun.txt"), None),
    "combined_sunspot_groups.csv (calibration variants, revision)": (os.path.join(PROJECT, "combined_sunspot_groups.csv"), None),
    "SILSO monthly smoothed SN (epochs, revision)": (r"C:\Users\Shivansh Mishra\Downloads\SN_ms_tot_V2.0.txt", None),
    "SILSO hemispheric catalogue (n6 replication)": (os.path.join(PAPER, REV, "external", "silso_hem", "Catalogue_B.txt"), None),
    "SILSO hemispheric monthly SN (n6 replication)": (os.path.join(PAPER, REV, "external", "silso_hem", "SN_m_hem_V2.0.txt"), None),
    "Tilt workbook, Jun-17 copy (phase 4, addendum B)": (r"C:\Users\Shivansh Mishra\Downloads\scaling polar field\To discuss in meeting"
                                                         r"\Regarding Tilt\tilt_complete_analysis.xlsx", None),
    "Tilt workbook, Jun-22 copy (phase 4 sensitivity)": (r"C:\Users\Shivansh Mishra\Downloads\prdiction cycle26\final sft paper"
                                                         r"\tilt_complete_analysis.xlsx", None),
    "RGO_SOON.txt (old catalogue; legacy steps only)": (os.path.join(PROJECT, "RGO_SOON.txt"), "5b897948aa0213d934f4efa206368fd2"),
}

# ------------------------------------------------------------------------------------------------------ registry ----
@dataclass
class Step:
    id: str
    stage: str
    cwd: str                      # relative to paper/
    script: str
    what: str
    outputs: list = field(default_factory=list)     # relative to cwd
    deps: list = field(default_factory=list)
    args: list = field(default_factory=list)
    cost: str = "fast"            # fast < 1 min, slow 1-10 min, heavy > 10 min
    kind: str = "analysis"        # analysis | figure | check
    legacy: bool = False
    note: str = ""

S = []
def add(*a, **k): S.append(Step(*a, **k))

# --- stage 1: corrected-catalogue reanalysis (paper/corrected_reanalysis) --------------------------------------------
st = "reanalysis"
add("cr.engine_tests", st, CR, "t_engine_tests.py", "unit tests of the SFT engine (bn_sft) used everywhere", kind="check")
add("cr.e1", st, CR, "e1_empirical.py", "empirical channels: A->B (parity-signed X1 vs A_pf), A->C (X1 vs next-cycle area), size classes",
    ["out/e1_empirical.json"])
add("cr.p1", st, CR, "p1_window_scan.py", "window scan s=-5..7, e=s..8; family-wise permutation p; nested LOO (max-r inner criterion)",
    ["out/p1_window_scan.json"])
add("cr.p2", st, CR, "p2_robust.py", "robustness of the two headline correlations", ["out/p2_robust.json"])
add("cr.p3", st, CR, "p3_variants.py", "predictor variants (numerator x latitude normalisation), family-wise p, nested LOO",
    ["out/p3_variants.json"])
add("cr.p4", st, CR, "p4_window_opt.py", "continuous 0.25-yr window optimisation (over-fitting check)",
    ["out/p4_window_opt.json"], cost="slow")
add("cr.p5", st, CR, "p5_window_families.py", "which window search generalises: families F1-F5 with nested LOO",
    ["out/p5_window_families.json"], cost="slow")
add("cr.table1", st, CR, "r1_table1.py", "Table 1: whole-cycle asymmetry by size class (per-region and legacy row-level)",
    ["out/table1_region.csv", "out/table1_rowlevel.csv"], note="recovered 2026-10-24; compares with the stored tables, never overwrites")
add("cr.s0_new", st, CR, "s0_legacy.py", "legacy SFT pipeline on the corrected catalogue (shows what changed)",
    ["out/s0_legacy_new.json"], args=["new"], legacy=True)
add("cr.s0_old", st, CR, "s0_legacy.py", "legacy SFT pipeline on the old catalogue (reproduces the old manuscript numbers)",
    ["out/s0_legacy_old.json"], args=["old"], legacy=True)
add("cr.s1_joy", st, CR, "s1_budget.py", "SFT whole-run budget by population (hi / small / giant), deterministic Joy tilts",
    ["out/s1_new_joy.json", "out/s1_new_joy_series.npz"], args=["tilt=joy"])
add("cr.s1_att", st, CR, "s1_budget.py", "SFT budget, attenuated tilts", ["out/s1_new_att.json"], args=["tilt=att"])
add("cr.s1_scatter", st, CR, "s1_budget.py", "SFT budget, one tilt-scatter realisation", ["out/s1_new_scatter.json"], args=["tilt=scatter"])
add("cr.s1_old_group", st, CR, "s1_budget.py", "SFT budget on the old catalogue, Group key", ["out/s1_old_group_joy.json"],
    args=["catalogue=old_group"], legacy=True)
add("cr.s1_old_cyclegroup", st, CR, "s1_budget.py", "SFT budget on the old catalogue, (Cycle, Group) key",
    ["out/s1_old_cyclegroup_joy.json"], args=["catalogue=old_cyclegroup"], legacy=True)
add("cr.s2", st, CR, "s2_green.py", "per-region SFT contributions via the cap Green's function", ["out/s2_regions_joy.csv"],
    ["cr.s1_joy"], args=["tilt=joy"])
add("cr.s3", st, CR, "s3_attrib.py", "attribution table from the per-region contributions", ["out/s3_attrib_joy.csv"], ["cr.s2"])
add("cr.s4", st, CR, "s4_param_memory.py", "memory of the modelled asymmetry vs transport parameters", ["out/s4_param_memory.json"])
add("cr.s5", st, CR, "s5_tests.py", "giant share, SC21 depletion, area balancing, tilt-scatter seeds", ["out/s5_tests.json"],
    args=["nseeds=20"])
add("cr.x1", st, CR, "x1_run.py", "pre-registered SFT ensemble (101 tilt members x 21 readout configs, per-region contributions)",
    ["out/x1_ensemble.npz"], cost="heavy", note="~33 min; almost everything SFT-based reads this file")
add("cr.x2", st, CR, "x2_compare.py", "model vs observed A_pf for analogs A1-A3 (ratio and difference)", ["out/x2_compare.json"], ["cr.x1"])
add("cr.x3", st, CR, "x3_timeseries.py", "year-by-year |B_N|-|B_S|, model vs observed", ["out/x3_timeseries.json", "out/x3_series.npz"],
    ["cr.x1"], cost="slow")
add("cr.x4", st, CR, "x4_experiments.py", "giant removal, area-matched controls, hemisphere swap, every cycle",
    ["out/x4_experiments.json"], ["cr.x1"], cost="slow")
add("cr.x5", st, CR, "x5_sweeps.py", "size x latitude and cap x readout sweeps", ["out/x5_sweeps.json"], ["cr.x1"], cost="slow",
    note="~10 min")
add("cr.x7", st, CR, "x7_exploratory.py", "exploratory checks on the ensemble", ["out/x7_exploratory.json"], ["cr.x1"])
add("cr.x8", st, CR, "x8_consistency.py", "weighted chi2 consistency of model and predictor with the A_pf errors",
    ["out/x8_consistency.json"], ["cr.x1"])
add("cr.x9", st, CR, "x9_leverage.py", "leverage of the giant class (share of D / share of flux)", ["out/x9_leverage.json"], ["cr.x1"])
add("cr.x10", st, CR, "x10_sc14.py", "out-of-sample Cycle 14 test (A_pf from MWO vs prediction from SC15-24)",
    ["out/x10_sc14.json"], ["cr.x1"], note="recovered 2026-10-24; compares with the stored file, never overwrites")
add("cr.y1", st, CR, "y1_bc.py", "B->C: polar imbalance at minimum vs next-cycle area imbalance (Bhowmik 2019 reproduction), SFT",
    ["out/y1_bc.json", "out/y1_table.csv"], ["cr.x1"])
add("cr.y2", st, CR, "y2_exploratory.py", "parity split and exploratory B->C checks", ["out/y2_exploratory.json"], ["cr.x8"])
add("cr.z1", st, CR, "z1_hpb.py", "algebraic hemispheric polar budget: kernels at eta = 150, 250, 500",
    ["out/z1_hpb.json", "out/z1_hpb_eta250.npz"], ["cr.x1"], cost="slow")
add("cr.z2", st, CR, "z2_hpb_closed.py", "closed-form budget", ["out/z2_hpb_closed.json"], ["cr.z1"])
add("cr.z3", st, CR, "z3_hpb_test.py", "budget vs observations, nine predictors, family-wise p", ["out/z3_hpb_test.json"],
    ["cr.z1", "cr.z2", "cr.x10"])
add("cr.z4", st, CR, "z4_reversals.py", "the two sign changes (SC14, SC21): runs test, giants vs rest, SFT",
    ["out/z4_reversals.json", "out/z4_reversals.npz"], ["cr.x1"])

# --- figures from the reanalysis (written to corrected_reanalysis/figures_pub) --------------------------------------
st = "figures"
add("fig.main", st, CR, "pub_figures.py", "Figs 1-6 and A1: size bars, two channels, window map, size classes, SFT budget, timing, A_pf reductions",
    ["figures_pub/fig1_asymmetry_by_size.pdf", "figures_pub/fig2_two_channels.pdf", "figures_pub/fig3_window_map.pdf",
     "figures_pub/fig4_size_classes.pdf", "figures_pub/fig5_sft_budget.pdf", "figures_pub/fig6_timing.pdf",
     "figures_pub/figA1_apf_reduction.pdf"], ["cr.p1", "cr.s1_joy", "cr.s2", "cr.x8"], kind="figure",
    note="fig5 in the paper is the Phase-4 rebuild (fig.rebuild)")
add("fig.supp", st, CR, "x6_figs.py", "supplementary S1-S7: analogs, time series, giant experiments, depletion, sweeps, sources",
    ["figures_pub/figS1_analogs.pdf", "figures_pub/figS3_experiments.pdf", "figures_pub/figS7_sources.pdf"],
    ["cr.x2", "cr.x3", "cr.x4", "cr.x5", "cr.x9"], kind="figure", note="S3 in the paper is the Phase-4 rebuild")
add("fig.bc", st, CR, "y3_figs.py", "Fig 7 (B->C), Fig 8 (parity split), S8 (predictor comparison)",
    ["figures_pub/fig7_bc.pdf", "figures_pub/fig8_parity.pdf", "figures_pub/figS8_predictors.pdf"], ["cr.y1", "cr.y2", "cr.x8"],
    kind="figure", cost="slow")
add("fig.reversals", st, CR, "z5_figs.py", "Figs 9-13: reversals, butterflies, cumulative imbalance, SFT flip, budget kernel",
    ["figures_pub/fig9_reversals.pdf", "figures_pub/fig13_hpb_kernel.pdf"], ["cr.z1", "cr.z2", "cr.z4"], kind="figure",
    note="figs 9-12 in the paper are the Phase-4 rebuilds (fig.rebuild); fig13 is current")
add("fig.chains", st, CR, "z6_chain_figs.py", "dominance-chain bubble plots A->B and A->C (not in the manuscript)",
    ["figures_pub/fig_chain_AB_AC.pdf", "figures_pub/fig_chain_AC_sc12_25.pdf"], kind="figure")
add("fig.legacy_names", st, CR, "f1_figures.py", "the old figure set under its old names (figures_publication names)",
    ["figures/fig1_asymmetry_by_size.pdf"], ["cr.e1", "cr.p1", "cr.s1_joy", "cr.s2", "cr.s5"], kind="figure", legacy=True,
    note="superseded by fig.main / fig.supp")

# --- stage 2: referee-protocol phases 0-3 (paper/revision; order of revision/run_all.sh) ------------------------------
st = "revision"
add("rev.r0", st, REV, "r0_audit.py", "Phase 0: audit of every number against the corrected outputs", ["out/r0_audit.json"],
    ["cr.x4", "cr.x10", "cr.y1", "fig.main"])
add("rev.t1", st, REV, "t1_speccurve.py", "T1: specification curve over ~2 M definitions; max-statistic null",
    ["out/t1_results.json", "out/t1_r_c1_adopted.npz"], cost="slow", note="~8 min")
add("rev.t1fig", st, REV, "t1_fig.py", "T1 quick plots (superseded by Phase-4 S1, S2)", ["out/t1_speccurve.png"], ["rev.t1"],
    kind="figure", legacy=True)
add("rev.t2", st, REV, "t2_quadrupole.py", "T2: signed quadrupole S form of the same-cycle relation", ["out/t2_results.json"])
add("rev.t3", st, REV, "t3_eiv.py", "T3: errors-in-variables, attainable r under the published sigma", ["out/t3_results.json"], ["cr.x8"])
add("rev.t4", st, REV, "t4_channel2.py", "T4: channel 2 vs persistence, surrogates, rolling-origin forecasts", ["out/t4_results.json"])
add("rev.t5", st, REV, "t5_size_nulls.py", "T5: size nulls, area-matched comparisons, swap percentiles", ["out/t5_results.json"],
    ["cr.x1"], cost="slow", note="~7 min")
add("rev.t5b", st, REV, "t5b_alt.py", "T5b: alternative null", ["out/t5b_alt.json"])
add("rev.t6", st, REV, "t6_splice.py", "T6: RGO->SOON and MWO->WSO splices, nine calibrations", ["out/t6_results.json"])
add("rev.t7", st, REV, "t7_reversals.py", "T7: reversals, region-resampled runs test", ["out/t7_results.json"])
add("rev.t8", st, REV, "t8_multiplicity.py", "T8: multiplicity (BH / Holm) over all tests", ["out/t8_results.json", "out/t8_tests.csv"],
    ["cr.y1", "cr.z3"])
add("rev.t9", st, REV, "t9_large_regions.py", "T9: large-region carrier tests (addenda B2/B3)", ["out/t9_results.json"])
add("rev.t10", st, REV, "t10_polarity.py", "T10: polarity of the parity factor (addendum B2(iii))", ["out/t10_results.json"])
add("rev.n7", st, REV, "n7_prospective.py", "N7: registered prospective test (SC25); PROSPECTIVE.md is hash-stamped",
    ["out/n7_prospective.json"], ["cr.x8"])
add("rev.n1", st, REV, "n1_phaselag.py", "N1: north-south phase lag", ["out/n1_results.json"])
add("rev.n2", st, REV, "n2_reversal.py", "N2: does the first-reversing hemisphere set A_pf", ["out/n2_results.json"], ["rev.n1"])
add("rev.n6", st, REV, "n6_replication.py", "N6: replication on SILSO hemispheric sunspot numbers", ["out/n6_results.json"])
add("rev.n4", st, REV, "n4_kernel.py", "N4: ridge lag kernel, data vs SFT", ["out/n4_results.json"], ["cr.x3"])
add("rev.n3", st, REV, "n3_obsresponse.py", "N3: observation-sampled SFT readout (101 tilt members)", ["out/n3_results.json"],
    ["cr.x1", "cr.x3"], cost="heavy", note="~11 min")
add("rev.n5", st, REV, "n5_sft_variants.py", "N5: SFT transport variants (slow flow etc.)", ["out/n5_results.json"], cost="slow")
add("rev.n5e", st, REV, "n5_check_e.py", "N5 variant e diagnostic", ["out/n5_check_e.json"], ["rev.n5"])
add("rev.p3table4", st, REV, "p3_table4.py", "Phase 3: Table 4 recomputed at n = 12", ["out/p3_table4.json"])
add("rev.manifest", st, REV, "make_manifest.py", "revision/out manifest", ["out/MANIFEST.md"], kind="check")

# --- stage 3: Phase 4 (paper/revision/p4; order of p4/run_all.sh) ------------------------------------------------------
st = "phase4"
add("p4.ws_a", st, P4, "ws_a.py", "WS-A: hemispheric transport asymmetry (delta grid), unit tests, nested test",
    ["out/ws_a_results.json", "out/ws_a_runs.json", "out/ws_a_unittests.txt"], ["cr.x1"], cost="heavy", note="~20 min")
add("p4.ws_b", st, P4, "ws_b.py", "WS-B: measured Joy slope (tilt workbook) vs X1hat and |A_pf|", ["out/ws_b_results.json"])
add("p4.ws_b_extra", st, P4, "ws_b_extra.py", "WS-B exploratory follow-ups (F2)", ["out/ws_b_extra.json"], ["p4.ws_b"])
add("p4.ws_c_scatter", st, P4, "ws_c_scatter.py", "constant-scatter per-region contributions", ["out/ws_c_const_scatter.npz"],
    cost="slow")
add("p4.ws_c", st, P4, "ws_c.py", "WS-C: size x latitude x phase cells, attribution, matched nulls", ["out/ws_c_results.json",
    "out/ws_c_cells.csv"], ["p4.ws_c_scatter", "cr.x1"], cost="heavy", note="~15 min")
add("p4.ws_d", st, P4, "ws_d.py", "WS-D: sign models (what sets the sign of A_pf)", ["out/ws_d_results.json"],
    ["p4.ws_a", "p4.ws_b", "rev.n2"])
add("p4.ws_e", st, P4, "ws_e.py", "WS-E: n-2 seeding, sign patterns", ["out/ws_e_results.json"])
add("p4.ws_f", st, P4, "ws_f.py", "WS-F: the SC14 / SC21 sign changes (bootstrap, matched nulls, leave-k)", ["out/ws_f_results.json",
    "out/ws_f_tables.npz"], ["p4.ws_a", "p4.ws_b", "p4.ws_c_scatter", "rev.n1", "cr.z1"], cost="slow")
add("p4.controls", st, P4, "ws_controls.py", "headline controls: Spearman, drop-one, error-weighted slope", ["out/ws_controls.json"])
add("p4.ws_g", st, P4, "ws_g.py", "Addendum B: best window and best predictor incl. tilt-weighted X2 (nested LOO, permutations)",
    ["out/ws_g_results.json", "out/ws_g_grid.csv"], cost="slow")
add("p4.inner", st, P4, "check_inner_criterion.py", "F9: nested values under both inner criteria (max-r vs LOO error)",
    ["out/check_inner_criterion.json"], ["cr.p5"])
add("fig.wsf", "figures", P4, "ws_f_figs.py", "WS-F figures ws_f1-f6 (180 mm)", ["figures/ws_f1_group_asymmetry.pdf"], ["p4.ws_f"],
    kind="figure")
add("fig.rebuild", "figures", P4, "fig_rebuild.py",
    "the six CURRENT manuscript figures: Fig 5, 9, 10, 11, 15 (S3), 22 -> p4/figures_manuscript (installed by hand into "
    "figures_v2 and paper_figures)", ["figures_manuscript/fig9_reversals.pdf", "figures_manuscript/fig12_sft_flip.pdf"],
    ["p4.ws_f", "cr.x1", "cr.x4", "cr.s1_joy", "rev.t5"], kind="figure")
add("fig.sets", "figures", P4, "folder_figs.py", "analysis figure sets: best window/predictor G1-G6, statistics S1-S9, mechanism M1-M5",
    ["sets/best/G1_map_TA.pdf", "sets/stats/S1_specification_curve.pdf", "sets/mech/M1_transport_asymmetry.pdf"],
    ["p4.ws_g", "p4.ws_a", "p4.ws_b", "p4.ws_c", "p4.ws_d", "p4.ws_e", "p4.controls", "rev.t1", "rev.t3", "rev.t4",
     "rev.n2", "rev.n3", "rev.n4", "rev.n5e"], kind="figure")
add("paperfigs.sets", "paperfigs", P4, "assemble_folders.py", "copy the figure sets into paper_figures/5-7 and verify every file",
    [], ["fig.sets", "fig.wsf"], kind="check")
add("p4.manifest", st, P4, "make_manifest.py", "Phase-4 manifest", ["MANIFEST_P4.md"], kind="check")

# --- stage 4: journal-ready figure folder (TrueType re-render; run by hand originally) ---------------------------------
add("paperfigs.render", "paperfigs", CR, os.path.join("..", "paper_figures", "_build", "render_type42.py"),
    "re-render Figs 1-13, A1, S1-S8 with TrueType fonts into a staging folder", [], ["fig.main", "fig.supp", "fig.bc", "fig.reversals"],
    args=[os.path.join(HERE, "stage_type42")], kind="figure", legacy=True,
    note="the original build of paper_figures/; six of those figures were later replaced by fig.rebuild, so the assembly "
         "step (build_paper_figures.py) now reports differences for them -- run by hand only")

# every Phase-4 script that imports p4lib reads corrected_reanalysis/out/y1_bc.json through p4lib.data()
for s in S:
    if s.cwd == P4 and s.script not in ("check_inner_criterion.py", "assemble_folders.py", "make_manifest.py") \
            and "cr.y1" not in s.deps:
        s.deps.append("cr.y1")
STEPS = {s.id: s for s in S}
STAGES = ["reanalysis", "revision", "phase4", "figures", "paperfigs"]

# ------------------------------------------------------------------------------------------------- comparison ----
def _sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def _same_json(a, b, tol):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same_json(a[k], b[k], tol) for k in a)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(_same_json(x, y, tol) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        return (a != a and b != b) or abs(a - b) <= tol * max(1.0, abs(a), abs(b))
    return a == b

def compare(old, new, tol=1e-9):
    """'identical' | 'same within tol' | 'DIFFERENT' | 'regenerated (figure)'."""
    ext = os.path.splitext(new)[1].lower()
    if _sha(old) == _sha(new):
        return "identical"
    if ext in (".pdf", ".png"):
        return "regenerated (figure; renders carry timestamps)"
    try:
        if ext == ".json":
            return "same within tol" if _same_json(json.load(open(old)), json.load(open(new)), tol) else "DIFFERENT"
        if ext == ".csv":
            import pandas as pd, numpy as np
            a, b = pd.read_csv(old), pd.read_csv(new)
            if a.shape != b.shape or list(a.columns) != list(b.columns):
                return "DIFFERENT"
            num = a.select_dtypes("number").columns
            ok = np.allclose(a[num].values, b[num].values, rtol=tol, atol=tol, equal_nan=True) and \
                 a.drop(columns=num).equals(b.drop(columns=num))
            return "same within tol" if ok else "DIFFERENT"
        if ext == ".npz":
            import numpy as np
            a, b = np.load(old, allow_pickle=True), np.load(new, allow_pickle=True)
            if set(a.files) != set(b.files):
                return "DIFFERENT"
            for k in a.files:
                x, y = a[k], b[k]
                if x.shape != y.shape:
                    return "DIFFERENT"
                if x.dtype.kind in "fc":
                    if not np.allclose(x, y, rtol=tol, atol=tol, equal_nan=True):
                        return "DIFFERENT"
                elif not (x == y).all():
                    return "DIFFERENT"
            return "same within tol"
    except Exception as e:                                           # pragma: no cover
        return f"not comparable ({e.__class__.__name__})"
    return "DIFFERENT (text/binary)"

# ------------------------------------------------------------------------------------------------- execution ----
def outpath(s, o):
    return os.path.normpath(os.path.join(PAPER, s.cwd, o))

def resolve(targets, legacy, only):
    want = []
    for t in targets:
        if t == "all":
            want += [s.id for s in S]
        elif t in STAGES:
            want += [s.id for s in S if s.stage == t]
        elif t in STEPS:
            want.append(t)
        else:
            sys.exit(f"unknown target {t!r}; see `python master.py list`")
    if not legacy:
        explicit = set(t for t in targets if t in STEPS)
        want = [w for w in want if not STEPS[w].legacy or w in explicit]
    order, seen = [], set()
    def visit(i, chain=()):
        if i in seen:
            return
        if i in chain:
            sys.exit("dependency cycle: " + " -> ".join(chain + (i,)))
        if not only:
            for d in STEPS[i].deps:
                visit(d, chain + (i,))
        seen.add(i); order.append(i)
    for w in want:
        visit(w)
    return order

def run(order, heavy, verify, dry):
    os.makedirs(LOGS, exist_ok=True)
    report, t_all = [], time.time()
    for i in order:
        s = STEPS[i]; outs = [outpath(s, o) for o in s.outputs]
        have = outs and all(os.path.exists(o) for o in outs)
        if s.cost == "heavy" and not heavy and have:
            print(f"  reuse  {i:22s} (heavy; outputs exist -- pass --heavy to recompute)"); report.append(dict(step=i, status="reused"))
            continue
        cmd = [sys.executable, s.script] + s.args
        print(f"  run    {i:22s} [{s.cost}] {' '.join(cmd[1:])}   (cwd paper/{s.cwd})", flush=True)
        if dry:
            continue
        snap = None
        if verify:
            snap = tempfile.mkdtemp(prefix="verify_")
            for o in outs:
                if os.path.exists(o):
                    shutil.copy2(o, os.path.join(snap, os.path.basename(o)))
        t0 = time.time()
        with open(os.path.join(LOGS, f"{i}.log"), "w", encoding="utf-8") as log:
            rc = subprocess.call(cmd, cwd=os.path.join(PAPER, s.cwd), stdout=log, stderr=subprocess.STDOUT,
                                 env=dict(os.environ, PYTHONWARNINGS="ignore", MPLBACKEND="Agg", PYTHONIOENCODING="utf-8"))
        dt = time.time() - t0
        rec = dict(step=i, rc=rc, seconds=round(dt, 1))
        missing = [o for o in outs if not os.path.exists(o)]
        if rc != 0 or missing:
            rec["status"] = "FAILED" if rc else "missing outputs"; rec["missing"] = missing
            print(f"         -> {rec['status']} after {dt:.0f}s; see logs/{i}.log"); report.append(rec)
            break
        if verify:
            rec["compare"] = {}
            for o in outs:
                old = os.path.join(snap, os.path.basename(o))
                res = compare(old, o) if os.path.exists(old) else "new file"
                if res.startswith(("DIFFERENT", "not comparable")):
                    # never let a verification run replace a verified result: keep the stored file, park the new one
                    root, ext = os.path.splitext(o); parked = f"{root}.verify_new{ext}"
                    shutil.move(o, parked); shutil.copy2(old, o); res += f" (stored file restored; new result -> {os.path.basename(parked)})"
                rec["compare"][os.path.basename(o)] = res
            shutil.rmtree(snap, ignore_errors=True)
            print("         -> ok %.0fs  %s" % (dt, "; ".join(f"{k}: {v}" for k, v in rec["compare"].items())))
        else:
            print(f"         -> ok {dt:.0f}s")
        rec["status"] = "ok"; report.append(rec)
    if not dry:
        name = os.path.join(LOGS, f"{'verify' if verify else 'run'}_{time.strftime('%Y%m%d_%H%M%S')}.json")
        json.dump(report, open(name, "w"), indent=1)
        print(f"done in {time.time() - t_all:.0f}s; report {name}")
    return report

# ------------------------------------------------------------------------------------------------------ views ----
def cmd_list(stage=None):
    for stg in STAGES:
        if stage and stg != stage:
            continue
        print(f"\n== {stg}")
        for s in S:
            if s.stage != stg:
                continue
            tag = (" [legacy]" if s.legacy else "") + (f" [{s.cost}]" if s.cost != "fast" else "")
            print(f"  {s.id:22s} {s.cwd}/{s.script} {' '.join(s.args)}{tag}\n      {s.what}")
            if s.note:
                print(f"      note: {s.note}")

def cmd_inputs():
    import hashlib as _h
    ok = True
    for name, (p, md5) in INPUTS.items():
        if not os.path.exists(p):
            print(f"  MISSING  {name}\n           {p}"); ok = False; continue
        msg = "present"
        if md5:
            h = _h.md5()
            with open(p, "rb") as f:
                for b in iter(lambda: f.read(1 << 22), b""):
                    h.update(b)
            msg = "md5 OK" if h.hexdigest() == md5 else f"md5 MISMATCH ({h.hexdigest()})"; ok &= msg == "md5 OK"
        print(f"  {msg:9s} {name}")
    return ok

def cmd_status():
    for stg in STAGES:
        print(f"\n== {stg}")
        for s in (x for x in S if x.stage == stg):
            outs = [outpath(s, o) for o in s.outputs]
            if not outs:
                print(f"  {s.id:22s} (no file outputs)"); continue
            ex = [o for o in outs if os.path.exists(o)]
            when = time.strftime("%Y-%m-%d %H:%M", time.localtime(max(os.path.getmtime(o) for o in ex))) if ex else "-"
            print(f"  {s.id:22s} {len(ex)}/{len(outs)} outputs  last {when}")

def cmd_docs():
    """Markdown table of every step, generated from the registry (used for README.md)."""
    for stg in STAGES:
        print(f"\n### {stg}\n\n| Step | Script | Cost | What it does | Main outputs |\n|---|---|---|---|---|")
        for s in (x for x in S if x.stage == stg):
            outs = ", ".join(f"`{o}`" for o in s.outputs[:3]) + (" ..." if len(s.outputs) > 3 else "")
            what = s.what + (f" (*{s.note}*)" if s.note else "") + (" **[legacy]**" if s.legacy else "")
            print(f"| `{s.id}` | `{s.cwd}/{s.script}{(' ' + ' '.join(s.args)) if s.args else ''}` | {s.cost} | {what} | {outs} |")

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["list", "inputs", "status", "order", "run", "docs"])
    ap.add_argument("targets", nargs="*")
    ap.add_argument("--heavy", action="store_true"); ap.add_argument("--legacy", action="store_true")
    ap.add_argument("--verify", action="store_true"); ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", action="store_true")
    a = ap.parse_args()
    if a.command not in ("list", "docs") and not os.path.isdir(PAPER):     # list and docs only read the registry
        sys.exit(f"project not found: {PAPER}\nset ASYM_PROJECT to the 'Asymmetry analysis' folder")
    if a.command == "list":
        cmd_list(a.targets[0] if a.targets else None)
    elif a.command == "inputs":
        sys.exit(0 if cmd_inputs() else 1)
    elif a.command == "status":
        cmd_status()
    elif a.command == "docs":
        cmd_docs()
    else:
        if not a.targets:
            sys.exit("give a target: all | stage | step id")
        order = resolve(a.targets, a.legacy, a.only)
        if a.command == "order":
            for i in order:
                print(f"  {i:22s} [{STEPS[i].cost}] {STEPS[i].cwd}/{STEPS[i].script} {' '.join(STEPS[i].args)}")
        else:
            run(order, a.heavy, a.verify, a.dry_run)

if __name__ == "__main__":
    main()
