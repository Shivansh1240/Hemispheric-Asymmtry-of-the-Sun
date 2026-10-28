# extension/ — beyond the paper: mechanism, persistence and models

The paper (see `pipeline/`) establishes the two channels. This second stage asks why they work and how far models can reproduce them. It is organised around four questions:

| | Question |
|---|---|
| **Q0** | How do the different written forms of channel 1 relate, and which earlier numbers have a stored source? |
| **Q1** | What sets the polar-field asymmetry at a minimum? |
| **Q2** | Why that window, where does the parity sign come from, and what carries over to the next cycle? |
| **Q3** | Why does one hemisphere stay dominant for several cycles, and how does the dominance flip? |
| **Q4** | Is the asymmetry predetermined by the emergence record, or renewed stochastically, and can a dynamo reproduce it? |

The work was done in rounds, each written up as a report that only adds to the previous ones. Every number in a report is formatted from the stored JSON outputs by a script; none is typed by hand. That is why the reports are kept here exactly as generated.

| Report | Content |
|---|---|
| [`reports/REPORT_v3.md`](reports/REPORT_v3.md) | the base answers to Q0–Q4 with the project's 1-D SFT and its own stochastic 2-D dynamo, ensembles and the exploratory track A–E (`REPORT.md` in the original folder) |
| [`reports/REPORT_v4_closure.md`](reports/REPORT_v4_closure.md) | closure tests T1–T5 on persistence, and a claim-by-claim ledger of every figure |
| [`reports/REPORT_v5_additions.md`](reports/REPORT_v5_additions.md) | point-by-point comparison with Bhowmik (2019, A&A 632, A117) |
| [`reports/REPORT_v7_surya.md`](reports/REPORT_v7_surya.md) | Bhowmik's dynamo (Passos et al. 2014 on SURYA, with pumping) run as a 6-seed ensemble |
| [`reports/REPORT_v8_bhowmik.md`](reports/REPORT_v8_bhowmik.md) | every earlier simulation repeated in Bhowmik's SFT and dynamo settings |

## Main findings so far

Each result is labelled as in the reports: a *Sun result* comes from data, a *model result* from a simulation.

**Q0 — one relation, two forms**
- Channel 1 needs no separate parity factor when it is written against the *signed* polar sum S = (B_N + B_S)/(|B_N| + |B_S|). Because (−1)^(n+1) = −sgn B_N at every minimum, X̂1 · A_pf = −X1 · S identically, so r(X̂1, A_pf) and −r(X1, S) are the same number (`q0_reconcile.py`). The parity factor is therefore fixed by the observed polarity, not fitted. What still needs a physical explanation is why the sign comes out that way (see the even-cycle switch under Q2).
- The "r ≈ −0.98, sign reversal" value quoted for channel 1 in an earlier brief is this same relation in first differences: r(X1, S) = −0.94 in levels (n = 10, all 10 signs reversed) and −0.98 in first differences (n = 9). It is not a separate or stronger result and should not be quoted as one; channel 1 is quoted once, as +0.93. `q0_candidates.py` lists every simple candidate relation with its n.

**Q1 — what sets A_pf**
- The parity-signed area imbalance near maximum tracks A_pf: r = +0.93, n = 10, 10/10 signs. Redrawing A_pf within its errors brings this to 0.74 (0.55–0.85). The relation is robust to area noise but fragile to polar-field errors (Sun result).
- It lives in the large regions: r = 0.97 for regions ≥ 1000 MSH, 0.94 for ≥ 500, 0.58 for < 100 (Sun result).
- The SFT reproduces A_pf only in odd cycles: r = +0.97 odd, +0.01 even, +0.61 overall, in both the project's SFT and Bhowmik's (model result). The even-cycle failure is the main limit on every model statement.

**Q2 — window, sign, next cycle**
- In the 1-D SFT the polar imbalance is built from emergence at t_R + 2.3 to + 5.1 yr (16–84 %) (model result). The paper's windows are the top of a broad plateau in the scan. The window the SFT points to ranks only 74th and 51st of 104 (Sun result).
- Even cycles switch their dominant hemisphere after maximum (5 of 7), odd cycles never do (0 of 6). With declining-phase emergence, the A_pf sign comes out right in 9 of 10 cycles (r = +0.68) *without* any parity factor. The switch replicates on the Debrecen (DPD) and Temmer/SILSO records over the cycles they cover. This is a candidate explanation of the parity factor, but it is still **exploratory**: starting the window at t_R + 2 gives only 6/10.
- The magnitude of the polar asymmetry, the area imbalance and the N–S timing lag each track the next cycle (r = +0.62, +0.77, +0.63) (Sun result). The area link is mostly persistence, and the lag link leans on SC19/20.

**Q3 — persistence and flips**
- Each pole is built by both hemispheres (about 60 % own, 40 % opposite, through cross-equatorial leading flux), and the imbalance |D| is a small residual of large opposing contributions (model result).
- The area imbalance changes slowly from cycle to cycle: lag-1 r = +0.57 (n = 12), sign runs [2, 7, 4] against 7.0 expected for independent cycles (Sun result).
- Wherever two consecutive cycles both have a firm sign, the sign is kept (4 of 4 pairs). Every sign change goes through a near-balanced cycle whose own sign a region bootstrap does not fix (SC14: 67 % of draws, SC21: 61 %) (Sun result).
- The 1-D SFT carries no memory of the imbalance from one cycle to the next, with or without decay (model result). A quasi-periodic slow mode (AR(2), ~6 cycles) is suggested but not established.

**Q4 — stochastic or predetermined, and the dynamo**
- Hemispheric totals do not predetermine A_pf; the timing and latitude pattern of emergence does (model result).
- The three headline correlations survive removing the largest regions of each cycle (X̂1 → A_pf stays 0.93 at k = 1), but the 7-cycle northern run and the SC21 southern sign do not: removing the single largest region per cycle breaks the run, and removing the top 10 turns A(21) positive (Sun result).
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

| Script | Question | Does |
|---|---|---|
| `core.py` | shared | the layer every script imports: catalogue, A_pf, the bn_sft engine through `pipeline/asymcore.py`, a cap-70 Green's function rescalable in decay time, a readout that keeps every region's contribution to every minimum, and catalogue operations (mirror, symmetrise, rescale, move a region) |
| `v0_checks.py` | engine | checks before any number is quoted: the cap-70 Green's function reproduces ensemble member 0; decay-time rescaling equals a direct build; mirror checks (N ↔ S) of model and data, also with hemispheric transport |
| `q0_reconcile.py` | Q0 | channel 1 in its parity-signed and signed-field (S) forms, first differences, and the magnitude/sign pieces that link it to the next cycle |
| `q0_candidates.py` | Q0 | every simple candidate relation (signed/unsigned, odd/even subsets) against A_pf, S and A(n+1), with n |
| `q1_sources.py` | Q1 | from a symmetrised catalogue (identical hemispheres, A_pf = 0), switches on one difference at a time: area, amount only, placement only, tilt scatter, transport, all; plus X̂1 × T_lin and an algebraic hemispheric dipole |
| `q2_window.py` | Q2 | exact per-region attribution of the modelled polar imbalance, the "physics window" chosen by the model alone, the full data scan, size × latitude table, the effective sign and next-cycle sign rules |
| `q2b_phase_sign.py` | Q2 | the parity sign modelled explicitly: X1 at maximum vs in the declining-phase window t_R+3 → +4 (chosen from the model only), the even/odd switch, A_pf against the unsigned decline-phase X1, neighbouring-window sensitivity |
| `q3_persistence.py` | Q3 | run lengths, n vs n+2 seeding, the flip cycles with and without low-latitude large regions, the SC21 flip on DPD and Temmer et al. (2006) records, model memory vs decay time |
| `q4_stochastic.py` | Q4 | randomised emergence in the SFT (sizes re-dealt within a hemisphere; time/latitude slots re-dealt between hemispheres; tilt only), and the key correlations without the k largest regions |
| `q4_dynamo.py` | Q4 | the 2-D Babcock–Leighton flux-transport dynamo on the full sphere, with stochastic forcing per hemisphere (σ, τ_c, ρ) and a cross-equatorial coupling knob κ; transport operator imported read-only from `ftd_core.py` |
| `q4_dynamo_analyse.py` | Q4 | per-cycle statistics of a dynamo run, defined exactly as for the Sun (A, A_pf, S, X at maximum and in the decline) |
| `q4_dynamo_checks.py` | Q4 | symmetry without noise, exact mirror with swapped forcing, and resolution (97 × 181 vs 65 × 121) |
| `q4_dynamo_ensemble.py` | Q4 | the dynamo table over seeds, with pooled r over all adjacent cycle pairs |
| `q4_ftd_verify.py` | Q4 | runs the transport solver's own eight verification checks, unchanged |
| `explore.py` | A–E | exploratory: move the largest low-latitude region to 35°, a Green's-function predictor X_G, the fewest regions that make up each cycle's imbalance, an AR(2) memory model, N–S timing |
| `ens_verify.py` | H1–H3 | 100-member emergence ensemble (area noise, A_pf within 1σ, tilt scatter) and 27 transport configurations for the three headline results |
| `ens_extra.py` | H1, Q3 | splits the H1 spread into area and A_pf parts; region bootstrap of the flip-cycle asymmetries |
| `ens_transport_series.py` | H1 | per-cycle modelled A_pf for the 27 transport configurations (band in the ensemble figure) |
| `figs_main.py` | all | plot types P1–P7 and the Q2 table |
| `figs_q.py` | all | one figure per question (Q0–Q3), the ensemble figure, P6 and the exploratory figures |
| `q4_figs.py` | Q4 | randomised emergence, the correlations without the largest regions, and the 2-D dynamo figures |
| `closure_tests.py` | v4 | T1 (timing-lag persistence), T2 (firm pairs), T5 (AR(2) check), T3/T4 summaries |
| `closure_dynamo.py` | v4 | T3 (forcing coherence time) and T4 (constant hemispheric bias) runs of the 2-D dynamo |
| `closure_v5.py` | v5 | pure-lag control, parity switch on independent records, bridge to Bhowmik's 0.91, her sensitivities, her link |
| `surya_runs.py` | v7 | SURYA runs of Bhowmik's dynamo: validation, mirror pair, 6-seed ensemble |
| `surya_v8.py` | v8 | the earlier dynamo experiments repeated in SURYA (ρ, κ, flow, τ_c, bias) |
| `sft_bhowmik.py` | v8 | reruns every SFT analysis above with Bhowmik's settings (no decay, her tilt law) by patching the shared core |

Run order for v3 is the table at the end of [`reports/REPORT_v3.md`](reports/REPORT_v3.md) (`v0_checks` first; `q2b_phase_sign` needs `explore.json`, so run it after `explore.py`), then the closure, v5, v7 and v8 scripts.

## What is needed to rerun this

All analysis scripts of the extension are here, including the shared `core.py` and the 2-D dynamo. What is still **not in the repository**:

- `pipeline/asymcore.py` itself works, but it needs the project's `paper/` tree (`lib_asym`, `common`, `p4lib`, `x_engine`, `bn_sft`, and `sft_asym`, which `core.py` imports directly); see `../pipeline/README.md`;
- `ftd_core.py` and `ftd_verify.py`, the 2-D transport solver from the separate declining-phase project (`sc24_final/`), which `q4_dynamo.py` imports. Set `FTD_DIR` to its folder;
- `surya/surya_bhowmik.f` and `surya/surya_bhowmik_ext.f` (see `../dynamo/README.md`);
- `plots_src/` (the figure and report writers of v4–v8) and the stored outputs under `out/` that later steps read;
- data: `../Datasets Used/tilt_complete_analysis.xlsx` (Q1, v5) and the Debrecen DPD file, whose absolute Windows path is written into `q3_persistence.py` (line 95).

The scripts expect to sit in a folder next to the pipeline folder. `core.py` looks for `asymcore.py` in `../master_code/` (the original layout) and `../pipeline/` (this repository). `surya_runs.py` reads `../master_code/init.dat`, which is `dynamo/init.dat` here.
