# pipeline/ — one entry point for the whole analysis

Two files:

- **`master.py`** knows every analysis and figure script of the project (about 100 scripts in four stages): their inputs, outputs, dependencies and run time. It runs any subset in dependency order and can check that a rerun reproduces the stored result.
- **`asymcore.py`** puts the reusable pieces (data layer, predictor, polar-field reductions, statistics, the SFT model, plotting style) behind one import, so new work uses exactly the data, conventions and engine that produced the published numbers.

Neither file re-implements anything. They call the verified project modules (`lib_asym`, `common`, `p4lib`, `x_engine`, `bn_sft`) and the stage scripts under `paper/`. That tree has to be present for anything except `list` and `docs` to work.

## Where the project is found

`master.py` and `asymcore.py` look for the project folder (the one that holds `paper/` and the catalogue) in this order:

1. the environment variable `ASYM_PROJECT`;
2. the root of this repository, if it contains `paper/`;
3. the original location on the author's machine.

```bash
export ASYM_PROJECT="/path/to/Asymmetry analysis"
```

Two inputs still have absolute paths inside the project code: the SILSO monthly SN (`revision/common.py`) and the tilt workbooks (`revision/p4/p4lib.py`). Change those two lines if you move the data.

## Quick start

```bash
python master.py inputs          # external data: catalogue md5, polar, SILSO and tilt files
python master.py list            # every step: id, script, cost, what it does
python master.py status          # which outputs exist and when they were made
python asymcore.py --sft         # self-test of the core library, ~30 s
python master.py run figures     # remake all figures from the stored analysis outputs
python master.py run all         # the whole project; heavy steps are reused if their outputs exist
python master.py run all --heavy # really everything, ~2 h
python master.py run cr.x1 --heavy --only --verify   # rerun one step and compare with the stored result
```

A target is `all`, a stage (`reanalysis`, `revision`, `phase4`, `figures`, `paperfigs`) or one or more step ids (`cr.p1 rev.t4`).

| Option | Effect |
|---|---|
| `--heavy` | also recompute steps over 10 min; by default they are reused when their outputs exist |
| `--only` | skip dependencies |
| `--legacy` | include superseded steps (old catalogue, old figure names) |
| `--dry-run` | print the plan, run nothing |
| `--verify` | snapshot each step's outputs, rerun, compare: JSON to 1e-9, CSV/NPZ with `allclose`, anything else by checksum |

`--verify` never overwrites a stored result. If a rerun differs, the stored file is put back and the new one is kept beside it as `*.verify_new.*`. Each step logs to `logs/<step>.log` and each run writes a JSON report to `logs/`.

## The four stages

| Stage | Folder | What it covers |
|---|---|---|
| `reanalysis` | `paper/corrected_reanalysis/` | the corrected-catalogue reanalysis: both channels, window scans, SFT budget and the 101-member ensemble (`cr.x1`, ~33 min), experiments, B→C, the algebraic polar budget, the two sign changes |
| `revision` | `paper/revision/` | the referee protocol, phases 0–3: specification curve, errors-in-variables, channel 2 vs persistence, size nulls, splices, multiplicity, phase lag, lag kernel, SILSO replication, the prospective SC25 test |
| `phase4` | `paper/revision/p4/` | what sets the amplitude and sign of A_pf (WS-A to WS-F), the headline controls, best window and predictor (incl. the tilt-weighted X2) |
| `figures` | the three folders above | every figure |
| `paperfigs` | `paper/paper_figures/` | copying the figure sets into the journal folder and checking them |

Run times: heavy (>10 min) are `cr.x1` ≈ 33 min, `p4.ws_a` ≈ 20, `p4.ws_c` ≈ 15, `rev.n3` ≈ 11. Slow (1–10 min) include `cr.x5` ≈ 10, `rev.t1` ≈ 8, `rev.t5` ≈ 7, `cr.z1` ≈ 2, `fig.bc` ≈ 1.5. Everything else takes seconds.

## Which script makes each manuscript figure

The manuscript compiles from `paper/figures_v2/`; journal copies are in `paper/paper_figures/`. Six figures (5, 9, 10, 11, 15, 22) come from the Phase-4 rebuild `fig.rebuild` and were installed by hand. The older versions that `fig.main`, `fig.supp` and `fig.reversals` still write to `corrected_reanalysis/figures_pub/` are superseded.

| Fig. | File | Made by |
|---|---|---|
| 1 | `fig1_asymmetry_by_size` | `fig.main` (pub_figures.py) |
| 2 | `fig2_two_channels` | `fig.main` |
| 3 | `fig3_window_map` | `fig.main` |
| 4 | `fig4_size_classes` | `fig.main` |
| 5 | `fig5_sft_budget` | **`fig.rebuild`** |
| 6 | `fig6_timing` | `fig.main` |
| 7 | `fig8_parity` | `fig.bc` (y3_figs.py) |
| 8 | `fig7_bc` | `fig.bc` |
| 9 | `fig9_reversals` | **`fig.rebuild`** |
| 10 | `fig11_flip_imbalance` | **`fig.rebuild`** |
| 11 | `fig12_sft_flip` | **`fig.rebuild`** |
| 12 (App.) | `figA1_apf_reduction` | `fig.main` |
| 13–19 | `figS1_analogs` … `figS7_sources` | `fig.supp` (x6_figs.py), except 15 |
| 15 | `figS3_experiments` | **`fig.rebuild`** |
| 20 | `figS8_predictors` | `fig.bc` |
| 21 | `fig13_hpb_kernel` | `fig.reversals` (z5_figs.py) |
| 22 | `fig10_flip_butterflies` | **`fig.rebuild`** |

Not in the manuscript: chain plots (`fig.chains`), WS-F figures (`fig.wsf`), and the analysis sets G1–G6 (best window/predictor), S1–S9 (statistics checks), M1–M5 (mechanism tests) from `fig.sets`.

## Building on it: `asymcore`

```python
import asymcore as A
d   = A.daily()                                    # daily rows, unique region ids, peak area per region
T   = A.targets()                                  # A_pf, S, A^n, A^(n+1), signed polar field, errors
ann = A.annual(d, mask=A.size_mask(d, lo=500, latmax=20))      # any sub-population
x   = A.x1_vector(A.C1, 0, 1, signed=True, ann=ann)            # X1hat over t_R+0..+1
A.corr(x, [T["Apf"][n] for n in A.C1])
pred, chosen = A.nested_window_selection(X_family, y)          # the selection rule the manuscript uses
rmax, j, p = A.perm_max_r(X_family, y)                          # family-wise p of a search maximum
m = A.SFTModel(); r = m.readout(seed=None); A.apf_of(r["full"][21])    # the SFT, read out like the data
A.ensemble_apf("A1")                               # modelled A_pf per cycle from the 101-member ensemble
A.save(fig + A.theme(), "name", 85, 65, "out")     # paper style: TrueType PDF + 600-dpi PNG
```

Any new predictor should go through `nested_window_selection` and `perm_max_r`, as the paper's own predictors did. A correlation picked from a scan is not a result until it has a nested value and a family-wise p.

### Conventions that must not change

Each of these was once the source of a real error.

- **Parity factor** (−1)^(n+1), applied only to polar-field targets. It is an empirical regularity (it equals −sgn B_N at every minimum in the record); do not call it a Hale-polarity correction.
- **Area targets** are unsigned.
- **Size classes** use the peak area of each region, which needs the unique region ids of `RGO_SOON_corrected.txt`.
- **Rows at Lat = 0** are dropped from hemispheric sums.
- **A_pf** is taken ±1.6 yr about the minimum that ends the cycle: MWO faculae up to SC19, WSO daily filtered field from SC20.

## Last verification (24 Oct 2026)

`python asymcore.py --sft` passed in 32 s: catalogue md5; channel 1 r = 0.9273, signs 10/10; channel 2 r = 0.7707, signs 10/12; Table 1 SC12 large class −0.2181; ensemble A1 modelled A_pf(21) −0.0594; a fresh SFT readout matching ensemble member 0 to 1e-6.

`master.py run … --verify` gave byte-identical data outputs for `cr.p1`, `cr.table1`, `cr.x10`, `cr.y1`, `cr.e1`, `cr.s1_joy`, `cr.z1`, `rev.t2`, `rev.n1`, `rev.t4`, `p4.controls` and `p4.ws_b`. `fig.chains` and `fig.bc` re-rendered without error (PDFs carry timestamps, so they are not compared byte by byte).

Not rerun that day: the heavy steps (`cr.x1`, `rev.n3`, `p4.ws_a`, `p4.ws_c`) and the remaining slow ones. Their stored outputs are what later steps read.

Two outputs had been made by one-off commands and now have scripts: `out/table1_region.csv` / `out/table1_rowlevel.csv` (`cr.table1`) and `out/x10_sc14.json` (`cr.x10`). Both reproduce the stored files exactly and never overwrite them. Still documented only, with no script: `out/x11_predwin_analogs.json` and `out/s6_cap_robustness.json` (no step reads them).

## Deliberately outside the pipeline

- Scripts that edited the manuscript or the referee reply (`revision/p3_apply_table4.py`, `revision/p3_stub.py`, `revision/p4/apply_p4_edits.py`, `revision/p4/apply_f9_edits.py`, `revision/p4/remap_reply_lines.py`). They were run once and must not be rerun blindly.
- `revision/p4/ws_g_build.py`, a one-off code split, already applied.
- `paper/sft_symmetrization/`, the legacy pipeline. It carries the old-catalogue numbers.
- `paper_figures/_build/build_paper_figures.py`, which built `paper_figures/` from the TrueType re-render. It now reports differences for the six rebuilt figures, so it is run by hand only.

## External inputs (`python master.py inputs`)

| File | Location in the project |
|---|---|
| `RGO_SOON_corrected.txt` (md5 `b78c7b08…`) | project root |
| `e_PField_MWO_WSO.csv` | project root |
| WSO daily polar file | `Dataset Comparison/new research/P3_meridional_flow/data/` |
| `combined_sunspot_groups.csv` | project root |
| SILSO `SN_ms_tot_V2.0.txt` | absolute path in `revision/common.py` |
| SILSO hemispheric catalogue and monthly SN | `revision/external/silso_hem/` |
| tilt workbook, Jun-17 and Jun-22 copies | absolute paths in `revision/p4/p4lib.py` |
| `RGO_SOON.txt` (legacy steps only) | project root |

## All steps

The table below is generated from the registry by `python master.py docs`.

### reanalysis

| Step | Script | Cost | What it does | Main outputs |
|---|---|---|---|---|
| `cr.engine_tests` | `corrected_reanalysis/t_engine_tests.py` | fast | unit tests of the SFT engine (bn_sft) used everywhere |  |
| `cr.e1` | `corrected_reanalysis/e1_empirical.py` | fast | empirical channels: A->B (parity-signed X1 vs A_pf), A->C (X1 vs next-cycle area), size classes | `out/e1_empirical.json` |
| `cr.p1` | `corrected_reanalysis/p1_window_scan.py` | fast | window scan s=-5..7, e=s..8; family-wise permutation p; nested LOO (max-r inner criterion) | `out/p1_window_scan.json` |
| `cr.p2` | `corrected_reanalysis/p2_robust.py` | fast | robustness of the two headline correlations | `out/p2_robust.json` |
| `cr.p3` | `corrected_reanalysis/p3_variants.py` | fast | predictor variants (numerator x latitude normalisation), family-wise p, nested LOO | `out/p3_variants.json` |
| `cr.p4` | `corrected_reanalysis/p4_window_opt.py` | slow | continuous 0.25-yr window optimisation (over-fitting check) | `out/p4_window_opt.json` |
| `cr.p5` | `corrected_reanalysis/p5_window_families.py` | slow | which window search generalises: families F1-F5 with nested LOO | `out/p5_window_families.json` |
| `cr.table1` | `corrected_reanalysis/r1_table1.py` | fast | Table 1: whole-cycle asymmetry by size class (per-region and legacy row-level) (*recovered 2026-10-24; compares with the stored tables, never overwrites*) | `out/table1_region.csv`, `out/table1_rowlevel.csv` |
| `cr.s0_new` | `corrected_reanalysis/s0_legacy.py new` | fast | legacy SFT pipeline on the corrected catalogue (shows what changed) **[legacy]** | `out/s0_legacy_new.json` |
| `cr.s0_old` | `corrected_reanalysis/s0_legacy.py old` | fast | legacy SFT pipeline on the old catalogue (reproduces the old manuscript numbers) **[legacy]** | `out/s0_legacy_old.json` |
| `cr.s1_joy` | `corrected_reanalysis/s1_budget.py tilt=joy` | fast | SFT whole-run budget by population (hi / small / giant), deterministic Joy tilts | `out/s1_new_joy.json`, `out/s1_new_joy_series.npz` |
| `cr.s1_att` | `corrected_reanalysis/s1_budget.py tilt=att` | fast | SFT budget, attenuated tilts | `out/s1_new_att.json` |
| `cr.s1_scatter` | `corrected_reanalysis/s1_budget.py tilt=scatter` | fast | SFT budget, one tilt-scatter realisation | `out/s1_new_scatter.json` |
| `cr.s1_old_group` | `corrected_reanalysis/s1_budget.py catalogue=old_group` | fast | SFT budget on the old catalogue, Group key **[legacy]** | `out/s1_old_group_joy.json` |
| `cr.s1_old_cyclegroup` | `corrected_reanalysis/s1_budget.py catalogue=old_cyclegroup` | fast | SFT budget on the old catalogue, (Cycle, Group) key **[legacy]** | `out/s1_old_cyclegroup_joy.json` |
| `cr.s2` | `corrected_reanalysis/s2_green.py tilt=joy` | fast | per-region SFT contributions via the cap Green's function | `out/s2_regions_joy.csv` |
| `cr.s3` | `corrected_reanalysis/s3_attrib.py` | fast | attribution table from the per-region contributions | `out/s3_attrib_joy.csv` |
| `cr.s4` | `corrected_reanalysis/s4_param_memory.py` | fast | memory of the modelled asymmetry vs transport parameters | `out/s4_param_memory.json` |
| `cr.s5` | `corrected_reanalysis/s5_tests.py nseeds=20` | fast | giant share, SC21 depletion, area balancing, tilt-scatter seeds | `out/s5_tests.json` |
| `cr.x1` | `corrected_reanalysis/x1_run.py` | heavy | pre-registered SFT ensemble (101 tilt members x 21 readout configs, per-region contributions) (*~33 min; almost everything SFT-based reads this file*) | `out/x1_ensemble.npz` |
| `cr.x2` | `corrected_reanalysis/x2_compare.py` | fast | model vs observed A_pf for analogs A1-A3 (ratio and difference) | `out/x2_compare.json` |
| `cr.x3` | `corrected_reanalysis/x3_timeseries.py` | slow | year-by-year |B_N|-|B_S|, model vs observed | `out/x3_timeseries.json`, `out/x3_series.npz` |
| `cr.x4` | `corrected_reanalysis/x4_experiments.py` | slow | giant removal, area-matched controls, hemisphere swap, every cycle | `out/x4_experiments.json` |
| `cr.x5` | `corrected_reanalysis/x5_sweeps.py` | slow | size x latitude and cap x readout sweeps (*~10 min*) | `out/x5_sweeps.json` |
| `cr.x7` | `corrected_reanalysis/x7_exploratory.py` | fast | exploratory checks on the ensemble | `out/x7_exploratory.json` |
| `cr.x8` | `corrected_reanalysis/x8_consistency.py` | fast | weighted chi2 consistency of model and predictor with the A_pf errors | `out/x8_consistency.json` |
| `cr.x9` | `corrected_reanalysis/x9_leverage.py` | fast | leverage of the giant class (share of D / share of flux) | `out/x9_leverage.json` |
| `cr.x10` | `corrected_reanalysis/x10_sc14.py` | fast | out-of-sample Cycle 14 test (A_pf from MWO vs prediction from SC15-24) (*recovered 2026-10-24; compares with the stored file, never overwrites*) | `out/x10_sc14.json` |
| `cr.y1` | `corrected_reanalysis/y1_bc.py` | fast | B->C: polar imbalance at minimum vs next-cycle area imbalance (Bhowmik 2019 reproduction), SFT | `out/y1_bc.json`, `out/y1_table.csv` |
| `cr.y2` | `corrected_reanalysis/y2_exploratory.py` | fast | parity split and exploratory B->C checks | `out/y2_exploratory.json` |
| `cr.z1` | `corrected_reanalysis/z1_hpb.py` | slow | algebraic hemispheric polar budget: kernels at eta = 150, 250, 500 | `out/z1_hpb.json`, `out/z1_hpb_eta250.npz` |
| `cr.z2` | `corrected_reanalysis/z2_hpb_closed.py` | fast | closed-form budget | `out/z2_hpb_closed.json` |
| `cr.z3` | `corrected_reanalysis/z3_hpb_test.py` | fast | budget vs observations, nine predictors, family-wise p | `out/z3_hpb_test.json` |
| `cr.z4` | `corrected_reanalysis/z4_reversals.py` | fast | the two sign changes (SC14, SC21): runs test, giants vs rest, SFT | `out/z4_reversals.json`, `out/z4_reversals.npz` |

### revision

| Step | Script | Cost | What it does | Main outputs |
|---|---|---|---|---|
| `rev.r0` | `revision/r0_audit.py` | fast | Phase 0: audit of every number against the corrected outputs | `out/r0_audit.json` |
| `rev.t1` | `revision/t1_speccurve.py` | slow | T1: specification curve over ~2 M definitions; max-statistic null (*~8 min*) | `out/t1_results.json`, `out/t1_r_c1_adopted.npz` |
| `rev.t1fig` | `revision/t1_fig.py` | fast | T1 quick plots (superseded by Phase-4 S1, S2) **[legacy]** | `out/t1_speccurve.png` |
| `rev.t2` | `revision/t2_quadrupole.py` | fast | T2: signed quadrupole S form of the same-cycle relation | `out/t2_results.json` |
| `rev.t3` | `revision/t3_eiv.py` | fast | T3: errors-in-variables, attainable r under the published sigma | `out/t3_results.json` |
| `rev.t4` | `revision/t4_channel2.py` | fast | T4: channel 2 vs persistence, surrogates, rolling-origin forecasts | `out/t4_results.json` |
| `rev.t5` | `revision/t5_size_nulls.py` | slow | T5: size nulls, area-matched comparisons, swap percentiles (*~7 min*) | `out/t5_results.json` |
| `rev.t5b` | `revision/t5b_alt.py` | fast | T5b: alternative null | `out/t5b_alt.json` |
| `rev.t6` | `revision/t6_splice.py` | fast | T6: RGO->SOON and MWO->WSO splices, nine calibrations | `out/t6_results.json` |
| `rev.t7` | `revision/t7_reversals.py` | fast | T7: reversals, region-resampled runs test | `out/t7_results.json` |
| `rev.t8` | `revision/t8_multiplicity.py` | fast | T8: multiplicity (BH / Holm) over all tests | `out/t8_results.json`, `out/t8_tests.csv` |
| `rev.t9` | `revision/t9_large_regions.py` | fast | T9: large-region carrier tests (addenda B2/B3) | `out/t9_results.json` |
| `rev.t10` | `revision/t10_polarity.py` | fast | T10: polarity of the parity factor (addendum B2(iii)) | `out/t10_results.json` |
| `rev.n7` | `revision/n7_prospective.py` | fast | N7: registered prospective test (SC25); PROSPECTIVE.md is hash-stamped | `out/n7_prospective.json` |
| `rev.n1` | `revision/n1_phaselag.py` | fast | N1: north-south phase lag | `out/n1_results.json` |
| `rev.n2` | `revision/n2_reversal.py` | fast | N2: does the first-reversing hemisphere set A_pf | `out/n2_results.json` |
| `rev.n6` | `revision/n6_replication.py` | fast | N6: replication on SILSO hemispheric sunspot numbers | `out/n6_results.json` |
| `rev.n4` | `revision/n4_kernel.py` | fast | N4: ridge lag kernel, data vs SFT | `out/n4_results.json` |
| `rev.n3` | `revision/n3_obsresponse.py` | heavy | N3: observation-sampled SFT readout (101 tilt members) (*~11 min*) | `out/n3_results.json` |
| `rev.n5` | `revision/n5_sft_variants.py` | slow | N5: SFT transport variants (slow flow etc.) | `out/n5_results.json` |
| `rev.n5e` | `revision/n5_check_e.py` | fast | N5 variant e diagnostic | `out/n5_check_e.json` |
| `rev.p3table4` | `revision/p3_table4.py` | fast | Phase 3: Table 4 recomputed at n = 12 | `out/p3_table4.json` |
| `rev.manifest` | `revision/make_manifest.py` | fast | revision/out manifest | `out/MANIFEST.md` |

### phase4

| Step | Script | Cost | What it does | Main outputs |
|---|---|---|---|---|
| `p4.ws_a` | `revision/p4/ws_a.py` | heavy | WS-A: hemispheric transport asymmetry (delta grid), unit tests, nested test (*~20 min*) | `out/ws_a_results.json`, `out/ws_a_runs.json`, `out/ws_a_unittests.txt` |
| `p4.ws_b` | `revision/p4/ws_b.py` | fast | WS-B: measured Joy slope (tilt workbook) vs X1hat and |A_pf| | `out/ws_b_results.json` |
| `p4.ws_b_extra` | `revision/p4/ws_b_extra.py` | fast | WS-B exploratory follow-ups (F2) | `out/ws_b_extra.json` |
| `p4.ws_c_scatter` | `revision/p4/ws_c_scatter.py` | slow | constant-scatter per-region contributions | `out/ws_c_const_scatter.npz` |
| `p4.ws_c` | `revision/p4/ws_c.py` | heavy | WS-C: size x latitude x phase cells, attribution, matched nulls (*~15 min*) | `out/ws_c_results.json`, `out/ws_c_cells.csv` |
| `p4.ws_d` | `revision/p4/ws_d.py` | fast | WS-D: sign models (what sets the sign of A_pf) | `out/ws_d_results.json` |
| `p4.ws_e` | `revision/p4/ws_e.py` | fast | WS-E: n-2 seeding, sign patterns | `out/ws_e_results.json` |
| `p4.ws_f` | `revision/p4/ws_f.py` | slow | WS-F: the SC14 / SC21 sign changes (bootstrap, matched nulls, leave-k) | `out/ws_f_results.json`, `out/ws_f_tables.npz` |
| `p4.controls` | `revision/p4/ws_controls.py` | fast | headline controls: Spearman, drop-one, error-weighted slope | `out/ws_controls.json` |
| `p4.ws_g` | `revision/p4/ws_g.py` | slow | Addendum B: best window and best predictor incl. tilt-weighted X2 (nested LOO, permutations) | `out/ws_g_results.json`, `out/ws_g_grid.csv` |
| `p4.inner` | `revision/p4/check_inner_criterion.py` | fast | F9: nested values under both inner criteria (max-r vs LOO error) | `out/check_inner_criterion.json` |
| `p4.manifest` | `revision/p4/make_manifest.py` | fast | Phase-4 manifest | `MANIFEST_P4.md` |

### figures

| Step | Script | Cost | What it does | Main outputs |
|---|---|---|---|---|
| `fig.main` | `corrected_reanalysis/pub_figures.py` | fast | Figs 1-6 and A1: size bars, two channels, window map, size classes, SFT budget, timing, A_pf reductions (*fig5 in the paper is the Phase-4 rebuild (fig.rebuild)*) | `figures_pub/fig1_asymmetry_by_size.pdf`, `figures_pub/fig2_two_channels.pdf`, `figures_pub/fig3_window_map.pdf` ... |
| `fig.supp` | `corrected_reanalysis/x6_figs.py` | fast | supplementary S1-S7: analogs, time series, giant experiments, depletion, sweeps, sources (*S3 in the paper is the Phase-4 rebuild*) | `figures_pub/figS1_analogs.pdf`, `figures_pub/figS3_experiments.pdf`, `figures_pub/figS7_sources.pdf` |
| `fig.bc` | `corrected_reanalysis/y3_figs.py` | slow | Fig 7 (B->C), Fig 8 (parity split), S8 (predictor comparison) | `figures_pub/fig7_bc.pdf`, `figures_pub/fig8_parity.pdf`, `figures_pub/figS8_predictors.pdf` |
| `fig.reversals` | `corrected_reanalysis/z5_figs.py` | fast | Figs 9-13: reversals, butterflies, cumulative imbalance, SFT flip, budget kernel (*figs 9-12 in the paper are the Phase-4 rebuilds (fig.rebuild); fig13 is current*) | `figures_pub/fig9_reversals.pdf`, `figures_pub/fig13_hpb_kernel.pdf` |
| `fig.chains` | `corrected_reanalysis/z6_chain_figs.py` | fast | dominance-chain bubble plots A->B and A->C (not in the manuscript) | `figures_pub/fig_chain_AB_AC.pdf`, `figures_pub/fig_chain_AC_sc12_25.pdf` |
| `fig.legacy_names` | `corrected_reanalysis/f1_figures.py` | fast | the old figure set under its old names (figures_publication names) (*superseded by fig.main / fig.supp*) **[legacy]** | `figures/fig1_asymmetry_by_size.pdf` |
| `fig.wsf` | `revision/p4/ws_f_figs.py` | fast | WS-F figures ws_f1-f6 (180 mm) | `figures/ws_f1_group_asymmetry.pdf` |
| `fig.rebuild` | `revision/p4/fig_rebuild.py` | fast | the six CURRENT manuscript figures: Fig 5, 9, 10, 11, 15 (S3), 22 -> p4/figures_manuscript (installed by hand into figures_v2 and paper_figures) | `figures_manuscript/fig9_reversals.pdf`, `figures_manuscript/fig12_sft_flip.pdf` |
| `fig.sets` | `revision/p4/folder_figs.py` | fast | analysis figure sets: best window/predictor G1-G6, statistics S1-S9, mechanism M1-M5 | `sets/best/G1_map_TA.pdf`, `sets/stats/S1_specification_curve.pdf`, `sets/mech/M1_transport_asymmetry.pdf` |

### paperfigs

| Step | Script | Cost | What it does | Main outputs |
|---|---|---|---|---|
| `paperfigs.sets` | `revision/p4/assemble_folders.py` | fast | copy the figure sets into paper_figures/5-7 and verify every file |  |
| `paperfigs.render` | `corrected_reanalysis/../paper_figures/_build/render_type42.py pipeline/stage_type42` | fast | re-render Figs 1-13, A1, S1-S8 with TrueType fonts into a staging folder (*the original build of paper_figures/; six of those figures were later replaced by fig.rebuild, so the assembly step (build_paper_figures.py) now reports differences for them -- run by hand only*) **[legacy]** |  |
