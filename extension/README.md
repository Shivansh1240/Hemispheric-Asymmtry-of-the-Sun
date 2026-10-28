# extension/ — beyond the paper: mechanism, persistence and models

The paper (see `pipeline/`) establishes the two channels. This second stage asks why they work and how far models can reproduce them. It is organised around four questions:

| | Question |
|---|---|
| **Q1** | What sets the polar-field asymmetry at a minimum? |
| **Q2** | Why that window, where does the parity sign come from, and what carries over to the next cycle? |
| **Q3** | Why does one hemisphere stay dominant for several cycles, and how does the dominance flip? |
| **Q4** | Is the asymmetry predetermined by the emergence record, or renewed stochastically, and can a dynamo reproduce it? |

The work was done in rounds, each written up as a report that only adds to the previous ones. Every number in a report is formatted from the stored JSON outputs by a script; none is typed by hand. That is why the reports are kept here exactly as generated.

| Report | Content |
|---|---|
| `REPORT.md` (v3, not in the repo yet) | the base answers to Q1–Q4 with the project's 1-D SFT and its own stochastic 2-D dynamo |
| [`reports/REPORT_v4_closure.md`](reports/REPORT_v4_closure.md) | closure tests T1–T5 on persistence, and a claim-by-claim ledger of every figure |
| [`reports/REPORT_v5_additions.md`](reports/REPORT_v5_additions.md) | point-by-point comparison with Bhowmik (2019, A&A 632, A117) |
| [`reports/REPORT_v7_surya.md`](reports/REPORT_v7_surya.md) | Bhowmik's dynamo (Passos et al. 2014 on SURYA, with pumping) run as a 6-seed ensemble |
| [`reports/REPORT_v8_bhowmik.md`](reports/REPORT_v8_bhowmik.md) | every earlier simulation repeated in Bhowmik's SFT and dynamo settings |

## Main findings so far

Each result is labelled as in the reports: a *Sun result* comes from data, a *model result* from a simulation.

**Q1 — what sets A_pf**
- The parity-signed area imbalance near maximum tracks A_pf: r = +0.93, n = 10, 10/10 signs. Redrawing A_pf within its errors brings this to 0.74 (0.55–0.85). The relation is robust to area noise but fragile to polar-field errors (Sun result).
- It lives in the large regions: r = 0.97 for regions ≥ 1000 MSH, 0.94 for ≥ 500, 0.58 for < 100 (Sun result).
- The SFT reproduces A_pf only in odd cycles: r = +0.97 odd, +0.01 even, +0.61 overall, in both the project's SFT and Bhowmik's (model result). The even-cycle failure is the main limit on every model statement.

**Q2 — window, sign, next cycle**
- In the 1-D SFT the polar imbalance is built from emergence at t_R + 2.3 to + 5.1 yr (16–84 %) (model result). The paper's windows are the top of a broad plateau in the scan. The window the SFT points to ranks only 74th and 51st of 104 (Sun result).
- Even cycles switch their dominant hemisphere after maximum (5 of 7), odd cycles never do (0 of 6). With declining-phase emergence, the A_pf sign comes out right in 9 of 10 cycles (r = +0.68) *without* any parity factor. The switch replicates on the Debrecen (DPD) and Temmer/SILSO records over the cycles they cover. This is a candidate explanation of the parity factor, but it is still **exploratory**: starting the window at t_R + 2 gives only 6/10.
- The magnitude of the polar asymmetry, the area imbalance and the N–S timing lag each track the next cycle (r = +0.62, +0.77, +0.63) (Sun result). The area link is mostly persistence, and the lag link leans on SC19/20.

**Q3 — persistence and flips**
- The area imbalance changes slowly from cycle to cycle: lag-1 r = +0.57 (n = 12), sign runs [2, 7, 4] against 7.0 expected for independent cycles (Sun result).
- Wherever two consecutive cycles both have a firm sign, the sign is kept (4 of 4 pairs). Every sign change goes through a near-balanced cycle whose own sign a region bootstrap does not fix (SC14: 67 % of draws, SC21: 61 %) (Sun result).
- The 1-D SFT carries no memory of the imbalance from one cycle to the next, with or without decay (model result). A quasi-periodic slow mode (AR(2), ~6 cycles) is suggested but not established.

**Q4 — stochastic or predetermined, and the dynamo**
- Hemispheric totals do not predetermine A_pf; the timing and latitude pattern of emergence does (model result).
- The project's own 2-D dynamo renews the asymmetry every cycle. It matches the Sun's persistence only if the source fluctuations stay coherent for ~30 yr (lag-1 r +0.58 at τ_c = 30 yr). A constant 1–10 % hemispheric bias does not create persistence (model result).
- **Bhowmik's dynamo on SURYA does carry persistence** at its standard noise: pooled lag-1 r = +0.45 (seed range +0.26 to +0.56), against +0.57 for the Sun (model result). So the lack of persistence was a property of the project's own model, not of flux-transport dynamos in general. Even there, the same-sign runs (4.0–5.5 per 13 cycles) are shorter than the Sun's 7.

**Against Bhowmik (2019)**
- Her polar-flux → next-cycle link (Pearson 0.73) comes out at +0.59 (n = 8) on her pairs and definition with the project's polar record, and +0.74 in signed form: consistent in sign, different in size.
- On her cycles 17–24 both SFTs reproduce the polar asymmetry at r ≈ 0.9; adding SC15–16 lowers it to +0.61.
- They differ on the sources. In the project's SFT, tilt scatter gives a median polar asymmetry of 3.7 % (hers: up to 36 %). A hemispheric flow difference gives 15–37 % (hers: about 3 %). A 2.5-yr hemispheric lag changes the final polar asymmetry by 30.6 % (hers: unaffected).

## Status

- **The v8 batch was stopped early** (27 Oct 2026, laptop overheating). Not run: τ_c = 3 and 5 yr at production noise, κ = 0.3 and 3, flow seeds 2–3, the T3 production runs (T3 in v8 uses the calibration runs) and all of T4 in SURYA. The v4 conclusions on T3/T4 are therefore not yet replaced.
- One claim from v3 does not survive v8: "placement, not amount, sets the modelled A_pf". In Bhowmik's SFT the amount of flux alone reaches r = 0.50.
- In Bhowmik's SFT, normalised A_pf amplitudes depend on the dipole carried over from the start. Quote D-based agreement (0.66) instead.

## Files here

| File | Does |
|---|---|
| `closure_tests.py` | T1 (timing-lag persistence), T2 (firm pairs, region bootstrap), T5 (AR(2) check), T3/T4 summaries |
| `closure_dynamo.py` | T3 (forcing coherence time) and T4 (constant hemispheric bias) runs of the 2-D dynamo |
| `closure_v5.py` | v5: pure-lag control, parity switch on independent records, bridge to Bhowmik's 0.91, her sensitivities, her link |
| `surya_runs.py` | v7: launches and analyses the SURYA runs of Bhowmik's dynamo (validation, mirror pair, 6-seed ensemble) |
| `surya_v8.py` | v8: the earlier dynamo experiments repeated in SURYA (ρ, κ, flow, τ_c, bias) |
| `sft_bhowmik.py` | v8: reruns every SFT analysis with Bhowmik's settings (no decay, her tilt law) by patching the shared core |

## What is needed to rerun this

These scripts import modules and read outputs that are **not in the repository yet**:

- `core.py` (catalogue, 1-D SFT model, mirror/symmetrised catalogues), `q4_dynamo.py`, `q4_dynamo_analyse.py`;
- the analysis scripts that `sft_bhowmik.py` reruns: `q1_sources`, `q2_window`, `q2b_phase_sign`, `q3_persistence`, `q4_stochastic`, `ens_verify`, `ens_extra`, `ens_transport_series`, `explore`, `figs_main`, `v0_checks`, `plots_src/`;
- `surya/surya_bhowmik.f` and `surya/surya_bhowmik_ext.f` (see `../dynamo/README.md`);
- `REPORT.md` (v3) and the stored outputs under `out/` that later steps read.

The scripts expect to sit in a folder next to the pipeline folder; `surya_runs.py` reads `../master_code/init.dat`, which is `dynamo/init.dat` here.
