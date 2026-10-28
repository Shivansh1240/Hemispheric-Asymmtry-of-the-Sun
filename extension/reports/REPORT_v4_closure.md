# Closure: persistence of hemispheric dominance (addendum to REPORT.md)

REPORT.md is unchanged; this file adds the closure tests T1–T5, each run once (stop rule), and a claim-by-claim ledger of every figure. Plain statistics only: Pearson r with n, sign counts, run lengths, 16–84 % bands. Every number below is formatted from `out/closure/*.json` or `plots/registry.json` by `plots_src/write_closure.py`; none is typed.

## Answer

**Why does one hemisphere's dominance persist for several cycles, yet not hold constantly?** *Sun result:* the area imbalance changes slowly from one cycle to the next (lag-1 r +0.57, n = 12; runs [2, 7, 4] against 7.0 for independent cycles), and wherever two consecutive cycles both have a firm sign, the sign is kept (4/4 pairs, T2). *Sun result:* it does not hold constantly, and where it changes the imbalance passes through zero: every sign change in SC12–24 goes through a near-balanced cycle whose own sign is not established (SC14, 67 % of draws; SC21, 61 % of draws), and 8 of the 12 consecutive pairs involve such a cycle (T2), so the persistence record itself rests on cycles whose sign a resampling of their regions does not fix. *Sun result:* the hemisphere that peaks later keeps doing so in the same runs (timing-lag sign runs [4, 4, 5], lag-1 r +0.34, n = 12, T1), but the lag agrees in sign with the area imbalance in 10/13 cycles, so it is the same phenomenon seen in timing, not an independent carrier. *Not established:* a quasi-periodic slow mode. The AR(2) fit implies a damped oscillation (|z| = 0.765, period 6.0 cycles), but an AR(1) record of 13 cycles gives a lag-2 coefficient below −0.47 in about one case in six (T5). *Model result:* neither model supplies the memory by itself. The 1-D SFT has none (Q3); the stochastic 2-D dynamo with production forcing renews both the area and the timing asymmetry every cycle (T1: lag-1 r of the lag −0.06 to +0.15); of the forcing coherence times run (1, 3, 5, 10 and 30 yr), only 30 yr, about 2 of its 14-yr cycles, gives the Sun's persistence: pooled lag-1 r is +0.58 (n = 204), the mean longest run 6.1, rms A_pf 0.12 and 94 % of cycles Sun-like, all close to the Sun, whereas coherence times of 1–10 yr give lag-1 r +0.01 to +0.19 (T3). In the model the persistence and the Sun-like polar behaviour appear together, not at each other's expense; the dynamo adds little memory of its own, and what would keep a hemisphere's source efficiency coherent for decades on the Sun is not established. *Model result:* a constant 1–10 % hemispheric difference in α0, v0 or η does not create cycle-to-cycle persistence (pooled lag-1 r −0.10 to +0.01); the 10 % α0 difference lengthens runs (mean longest 5.7) without raising lag-1 r (−0.02); the longer runs come with the north dominating 70 % of cycles (the Sun: 7/13), and no bias gives a same-sign polar offset in every cycle (A_pf > 0 in 47 %–56 % of cycles) (T4). In short: in the record, dominance persists while the imbalance stays large and changes sign only through a near-balanced cycle (Sun result); in the dynamo, a hemispheric source fluctuation that stays coherent for about 2 cycles gives the Sun's persistence, run length, polar asymmetry and Sun-like share together (model result); what would supply that coherence on the Sun is not established.

## T1. Persistence of the N–S timing lag

**Figure** `plots/Q3_persistence_and_flips/T1_phase_lag_persistence.png`.
**Number** Sun: lag-1 r of the per-cycle lag t_N − t_S = **+0.34 (n = 12)**; sign runs **[4, 4, 5]** (3 runs; independent cycles 7.0); longest 5 (independent 4.0). Dynamo, same definition, every production setting: pooled lag-1 r −0.06 to +0.15.
**Sentence** The Sun's timing lag keeps its sign for 4–5 cycles at a time, in as few runs as the area asymmetry; no dynamo setting carries it from one cycle to the next.

- The lag has the sign of A(n) in 10/13 cycles and of A(n+1) in 9/12, so it is not independent evidence.
- In the dynamo, the mean longest run of the lag sign is 3.4–5.2 cycles per 13; the highest is the single coarse-grid κ = 10 run.

## T2. Firm pairs

**Figure** `plots/Q3_persistence_and_flips/T2_firm_pairs.png`.
**Number** Region bootstrap of A(n), 1000 draws, all 13 cycles: **8/13 firm** (observed sign in ≥ 90 % of draws). Consecutive firm–firm pairs keeping the sign: **4/4**. Pairs involving a non-firm cycle: **8/12**.
**Sentence** Where both cycles are firm, dominance always persists, but most of the record's pairs involve a near-balanced cycle whose sign is not established.

| Cycle | A(n) | draws with the observed sign | firm |
|---|---|---|---|
| SC12 | −0.175 | 100 % | yes |
| SC13 | −0.102 | 99 % | yes |
| SC14 | +0.025 | 67 % | no |
| SC15 | +0.097 | 98 % | yes |
| SC16 | +0.088 | 97 % | yes |
| SC17 | +0.004 | 53 % | no |
| SC18 | +0.028 | 75 % | no |
| SC19 | +0.176 | 100 % | yes |
| SC20 | +0.170 | 100 % | yes |
| SC21 | −0.013 | 61 % | no |
| SC22 | −0.067 | 96 % | yes |
| SC23 | −0.066 | 96 % | yes |
| SC24 | −0.007 | 57 % | no |

- Firm–firm pairs: SC12–13, SC15–16, SC19–20, SC22–23. Sign changes: SC13→14, SC20→21; each passes through a non-firm cycle.
- The six cycles bootstrapped before reproduce within a few points: SC13 99 % → 99 %, SC14 69 % → 67 %, SC17 56 % → 53 %, SC20 100 % → 100 %, SC21 64 % → 61 %, SC24 55 % → 57 %.
- The bootstrap measures sampling spread of the regions, not measurement error.

## T3. How much forcing memory the dynamo needs

**Figure** `plots/Q4_stochastic_vs_predetermined/T3_coherence_time_sweep.png`.
**Number** The most persistent setting is τ_c = 30 yr: pooled lag-1 r +0.58 (n = 204), mean longest run 6.1 in 13 cycles, Sun-like cycles 94 %, rms A_pf 0.12. The Sun: lag-1 r +0.57, longest run 7, rms A_pf 0.10.
**Sentence** Of the coherence times run, only τ_c = 30 yr (about 2 cycles) lets the dynamo carry hemispheric dominance as the Sun does; at that setting the polar asymmetry and the Sun-like share are also the closest to the Sun of any setting run. With coherence of 10 yr or less it renews the asymmetry almost every cycle.

| τ_c (yr) | σ | rms A achieved | pooled lag-1 r of A (seed range; n pairs) | mean longest run / 13 (range) | Sun-like cycles (range) | rms A_pf |
|---|---|---|---|---|---|---|
| 1 | 0.600 | 0.105 | +0.01 (−0.13 to +0.16; 204) | 4.1 (3.9–4.9) | 13 % (0 %–33 %) | 0.50 |
| 3 | 0.297 | 0.130 | +0.05 (−0.09 to +0.16; 203) | 4.4 (4.0–4.4) | 53 % (26 %–78 %) | 0.34 |
| 5 | 0.242 | 0.129 | +0.12 (+0.06 to +0.15; 203) | 3.8 (3.5–3.9) | 4 % (1 %–41 %) | 0.42 |
| 10 | 0.161 | 0.094 | +0.19 (−0.00 to +0.31; 204) | 4.9 (3.8–5.6) | 74 % (1 %–77 %) | 0.30 |
| 30 | 0.140 | 0.088 | +0.58 (+0.53 to +0.62; 204) | 6.1 (5.3–7.0) | 94 % (88 %–94 %) | 0.12 |

- κ = 1, ρ = 0, 3 seeds per τ_c, 1100 yr each (first 100 yr dropped). τ_c = 1 yr reuses the production seeds 1–3 (σ = 0.6).
- σ was re-tuned once per τ_c from two 600-yr calibration runs: τ_c = 3: σ 0.25 gave rms A 0.084 → σ 0.297; τ_c = 5: σ 0.18 gave rms A 0.074 → σ 0.242; τ_c = 10: σ 0.12 gave rms A 0.074 → σ 0.161; τ_c = 30: σ 0.1 gave rms A 0.071 → σ 0.140.
- rms A_pf here is over all cycles of the three seeds pooled.

## T4. A constant hemispheric difference

**Figure** `plots/Q4_stochastic_vs_predetermined/T4_hemispheric_bias.png`.
**Number** With the noise on, a 1–10 % N−S difference in α0, v0 or η gives pooled lag-1 r of A **−0.10 to +0.01** (unbiased +0.01; Sun +0.57). The 10 % α0 difference makes the north dominate **70 %** of cycles (mean longest run 5.7); A_pf > 0 in 47 %–56 % of cycles across all biases.
**Sentence** A constant hemispheric bias does not make the area asymmetry persist from cycle to cycle; at 10 % in α0 it only tilts the odds toward one hemisphere, and it never gives a same-sign polar offset in every cycle. It fails the test.

| Bias (larger in N) | pooled lag-1 r of A (n pairs) | mean longest run / 13 | cycles with A > 0 | mean A | cycles with A_pf > 0 | mean A_pf | rms A_pf |
|---|---|---|---|---|---|---|---|
| reference (delta 0) | +0.01 (204) | 4.1 | 50 % | −0.00 | 52 % | +0.02 | 0.50 |
| alpha 1 % | −0.04 (203) | 4.1 | 51 % | +0.00 | 51 % | +0.02 | 0.48 |
| alpha 3 % | −0.10 (203) | 3.9 | 55 % | +0.01 | 52 % | −0.00 | 0.48 |
| alpha 10 % | −0.02 (208) | 5.7 | 70 % | +0.05 | 56 % | +0.05 | 0.45 |
| v0 1 % | −0.07 (207) | 3.9 | 50 % | +0.00 | 47 % | −0.05 | 0.43 |
| v0 3 % | −0.04 (207) | 3.6 | 50 % | +0.01 | 49 % | −0.01 | 0.43 |
| v0 10 % | −0.04 (206) | 4.4 | 47 % | −0.00 | 50 % | +0.02 | 0.45 |
| eta 1 % | −0.09 (203) | 4.5 | 49 % | −0.00 | 51 % | −0.02 | 0.49 |
| eta 3 % | −0.10 (207) | 4.1 | 49 % | −0.00 | 55 % | +0.02 | 0.47 |
| eta 10 % | +0.01 (205) | 3.7 | 50 % | +0.01 | 48 % | −0.04 | 0.44 |

- The Sun: A > 0 in 7/13 cycles (SC12–24); A_pf > 0 in 60 % of SC15–24.
- κ = 1, σ = 0.6, τ_c = 1 yr, ρ = 0; 3 seeds (1–3) per bias, the same noise realisations as the reference. A zero bias reproduces the production model (largest relative difference 6e-16; `out/closure/dynamo_checks.json`); a v0 bias of δ gives surface speeds v0 (1 ± δ/2) at 45° N / S.

## T5. The AR(2) memory fit

**Figure** `plots/Q3_persistence_and_flips/T5_ar2_check.png`.
**Number** From a = +0.77, b = −0.59: complex roots, **modulus 0.765, period 6.0 cycles**. Fitting AR(2) to 20000 13-cycle series simulated from the fitted AR(1) (a = 0.39) gives b in **−0.47 to +0.12** (16–84 %; median −0.18). The Sun's −0.59 lies below that band.
**Sentence** The fit describes a damped oscillation of about 6 cycles, but a record of 13 cycles cannot establish it: short AR(1) records already bias b negative, and about one in six gives b below the band's lower edge.

## Claim ledger (every figure, by question)

Each row is one figure's claim, as registered when the figure was built. Numbers come from `out/*.json` and were checked against REPORT.md (`plots/number_checks.md`, `plots/number_checks_T.md`; 0 mismatches).

| Question | Claim | Key numbers | Result type | Caveats | Figure |
|---|---|---|---|---|---|
| Q1 | Area imbalance near maximum, parity-signed, tracks A_pf | r=+0.93, n=10, 10/10; A_pf-redraw 0.74 (0.55-0.85) | Sun result | scan-chosen window; fragile to A_pf errors; sign needs the parity label | `plots/Q1_what_sets_polar_asymmetry/F1_x1hat_vs_apf.png` |
| Q1 | The A_pf relation lives in large regions; the next-cycle relation in all sizes | >=1000 MSH 0.97, >=500 0.94, <100 0.58; next cycle 0.65-0.85 | Sun result | scan-chosen windows; small cells noisy (region counts shown) | `plots/Q1_what_sets_polar_asymmetry/F2_size_latitude_tables.png` |
| Q1 | Where/when flux emerges, not how much, sets the modelled A_pf; agreement only in odd cycles | area 0.61 (odd 0.97, even 0.01); amount 0.20; placement 0.60; tilt 0.10; all 0.55 | model result | SFT matches odd cycles only; transport-only gives a constant offset | `plots/Q1_what_sets_polar_asymmetry/F3_sft_ablation.png` |
| Q1 | The SFT reproduces A_pf only in odd cycles | odd +0.97, even +0.01, all +0.61 | model result vs Sun | even-cycle failure is the main limit of every model statement | `plots/Q1_what_sets_polar_asymmetry/F4_model_vs_observed_per_cycle.png` |
| Q1 | H1 robust to area noise, fragile to A_pf errors; model agreement transport-independent | area 0.92-0.93; A_pf redraw 0.74 (0.55-0.85); transport r 0.56-0.63 | Sun result (a); model result (b) | redraw adds a second error realisation: a sensitivity, not an estimate | `plots/Q1_what_sets_polar_asymmetry/S1_robustness.png` |
| Q2 | The polar imbalance at minimum is built by declining-phase emergence because the imbalance kernel decays | imbalance t_R+2.3..+5.1 yr; kernel 1.00/0.27/0.008 | model result (1-D SFT, holds over 27 transport configs) | window edges shift +-0.7 yr with v0; 1-D model | `plots/Q2_window_sign_next_cycle/F1_green_function_epochs.png` |
| Q2 | The paper windows are maxima of a plateau; the physics window does not reproduce them | paper ranks 1, 1; physics ranks 74, 51 of 104 | Sun result | scan-chosen windows; n = 10 and 12 | `plots/Q2_window_sign_next_cycle/F2_window_scans.png` |
| Q2 | Even cycles switch dominant hemisphere after maximum; declining-phase emergence gives the A_pf sign without parity | switch even 5/7, odd 0/6; r = +0.68, 9/10 signs | exploratory (Sun data; window from the model) | window-sensitive (start at t_R+2 gives 6/10); n = 10; r below the parity-signed 0.93 | `plots/Q2_window_sign_next_cycle/F3_parity_as_phase_effect.png` |
| Q2 | Polar-asymmetry magnitude, area imbalance and timing lag each track the next cycle | \|A_pf\| +0.62 (n=9); X1 +0.77 (n=12); lag +0.63 (n=12; +0.53 w/o SC19/20) | Sun result | n = 9-12; (b) scan-chosen window, mostly persistence; (c) leans on SC19/20 | `plots/Q2_window_sign_next_cycle/F4_next_cycle_links.png` |
| Q2 | The modelled polar imbalance is essentially which hemisphere peaks later | +0.87, +0.87 (n = 13) | model result (a); Sun result (b); exploratory | exploratory; n = 13 | `plots/Q2_window_sign_next_cycle/F5_timing_lag.png` |
| Q3 | Area dominance persists (lives in >=500 MSH); A_pf does not | runs [2,7,4] vs 7.0 independent; A_pf 7 runs in 10 | Sun result | the 7-cycle run depends on single regions in SC17/SC24 (Q4 F2) | `plots/Q3_persistence_and_flips/F1_dominance_and_cumulative.png` |
| Q3 | Lag-1 persistence, no n+2 seeding; the next-cycle sign follows persistence | lag-1 +0.57 (n=12); lag-2 -0.28 (n=11); persistence 10/12 | Sun result | n = 11-12 | `plots/Q3_persistence_and_flips/F2_return_maps_sign_rules.png` |
| Q3 | Large low-latitude regions decide near-balanced cycles in both directions; flip signs are not established | bootstrap 64/69/56/55 % vs 99/100 % | Sun result | SC18 not bootstrapped; large-region carrier not seen in DPD (F4) | `plots/Q3_persistence_and_flips/F3_flip_anatomy.png` |
| Q3 | SC21 sign replicates on independent records; SC14 is a negative control | SC21: -0.013/-0.014/-0.018; SC14 model -0.032, MWO -0.147+-0.143 | Sun result (a); model result vs Sun (b) | magnitudes at noise level (bootstrap 64 %); carrier not seen in DPD; SC13 has no polar data | `plots/Q3_persistence_and_flips/F4_sc21_records_sc14_control.png` |
| Q3 | The 1-D SFT has no inter-cycle memory of the imbalance at any decay time | lag-1 0.37-0.44; old share 2e-06-5e-05 | model result | 1-D model; deterministic tilts in (b) | `plots/Q3_persistence_and_flips/F5_sft_memory_tau.png` |
| Q4 | Hemispheric totals do not predetermine A_pf; the emergence timing/latitude pattern does | 0.54 / 0.51 / 0.11; shares 0.77 / 0.67 / 0.28 | model result | SFT agreement is in odd cycles only | `plots/Q4_stochastic_vs_predetermined/F1_randomised_emergence.png` |
| Q4 | Correlations are robust to removing the largest regions; the persistence record and the SC21 sign are not | H1 0.93->0.83 (k=20); runs 3->6 at k=1; A(21) -0.013 -> +0.003 at k=10 | Sun result | 13 cycles; one realisation of the record | `plots/Q4_stochastic_vs_predetermined/F2_without_largest_regions.png` |
| Q4 | The dynamo renews hemispheric asymmetry every cycle, unlike the Sun's runs | Sun lag-1 +0.57; model whole run +0.00 | model result vs Sun | one seed shown; amplitudes resolution-sensitive; model polar asymmetry too large (Q4 F4) | `plots/Q4_stochastic_vs_predetermined/F3_dynamo_example.png` |
| Q4 | No dynamo setting carries hemispheric dominance; the model's polar asymmetry is too large | pooled lag-1 -0.05..+0.13 vs Sun +0.57; longest run 3.4-4.5 vs 7; rms A_pf 0.27-0.45 vs 0.10 | model result vs Sun | model not a faithful Sun analogue (polar asymmetry too large; amplitudes resolution-sensitive) | `plots/Q4_stochastic_vs_predetermined/F4_dynamo_seed_ensembles.png` |
| T | The Sun's N-S timing lag keeps its sign for 4-5 cycles at a time; the dynamo's lag renews every cycle | Sun lag-1 r +0.34 (n = 12), runs [4, 4, 5]; model pooled −0.06 to +0.15 | Sun result; model result | 13 cycles; lag and A(n) share sign in 10/13, so not independent evidence | `plots/Q3_persistence_and_flips/T1_phase_lag_persistence.png` |
| T | Where both cycles are firm, dominance persists (4/4); most pairs involve a near-balanced cycle whose sign is not established | 8/13 firm; firm-firm 4/4; 8/12 pairs with a non-firm cycle | Sun result | only 4 firm-firm pairs; the bootstrap is sampling spread, not measurement error | `plots/Q3_persistence_and_flips/T2_firm_pairs.png` |
| T | The negative lag-2 term is outside the AR(1) 16-84 % band but the record is too short to establish an oscillation | \|z\| = 0.765, period 6.0 cycles; AR(1) b band [−0.47, +0.12] vs −0.59 | Sun result (fit); not established | 11 triples; band is 16-84 %, so ~1 in 6 AR(1) series fall below it | `plots/Q3_persistence_and_flips/T5_ar2_check.png` |
| T | Persistence of A in the dynamo as a function of the forcing coherence time, with the polar asymmetry and Sun-like share alongside | tau_c 1: lag-1 +0.01, run 4.1, rms A_pf 0.50; tau_c 3: lag-1 +0.05, run 4.4, rms A_pf 0.34; tau_c 5: lag-1 +0.12, run 3.8, rms A_pf 0.42; tau_c 10: lag-1 +0.19, run 4.9, rms A_pf 0.30; tau_c 30: lag-1 +0.58, run 6.1, rms A_pf 0.12 | model result | 3 seeds; one calibration step for sigma (rms A 0.13 and 0.13 at tau_c 3 and 5, not 0.10); polar asymmetry too large except at tau_c 30 (0.12); 10-30 yr and > 30 yr not run | `plots/Q4_stochastic_vs_predetermined/T3_coherence_time_sweep.png` |
| T | Whether a constant hemispheric bias makes area dominance persist without a same-sign polar-field offset | alpha 10 %: lag-1 −0.02, A>0 70 %, A_pf>0 56 %; v0 10 %: lag-1 −0.04, A>0 47 %, A_pf>0 50 %; eta 10 %: lag-1 +0.01, A>0 50 %, A_pf>0 48 % | model result | 3 seeds; one model; dynamo polar asymmetry too large | `plots/Q4_stochastic_vs_predetermined/T4_hemispheric_bias.png` |

## What does not hold (closure tests)

- The dynamo does not carry dominance with forcing coherence of 10 yr or less (pooled lag-1 r +0.01 to +0.19) or with a constant 1–10 % hemispheric bias (T4). The σ re-tuning missed at τ_c = 3 and 5 yr (rms A 0.13 and 0.13, not 0.10). τ_c = 30 yr is the end of the sweep: coherence times between 10 and 30 yr and beyond 30 yr were not run (stop rule).
- The AR(2) oscillation (period about 6 cycles) is not established (T5).
- The persistence record leans on near-balanced cycles: 5/13 cycles are not firm (T2).
- The timing lag is not independent evidence for persistence (T1).

## Files

- `closure_tests.py` (T1, T2, T5, and the T3/T4 summaries), `closure_dynamo.py` (T3/T4 runs; hemispheric-bias subclass of the Q4 dynamo).
- `out/closure/`: T1–T5 JSON, T3 calibration, run log, dynamo runs, checks.
- `plots_src/qt.py`, `plots_src/qt_dyn.py` (figures), `plots_src/write_closure.py` (this file and `plots/CLOSURE.md`).
