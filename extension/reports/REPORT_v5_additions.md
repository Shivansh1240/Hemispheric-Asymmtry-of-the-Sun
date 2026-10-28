# v5 additions: comparison with Bhowmik (2019)

Addendum to REPORT.md and REPORT_v4_closure.md (both unchanged). Plain statistics: Pearson r with n, sign counts, differences, 16–84 % bands. Every number is formatted from `out/v5/*.json` by `plots_src/write_v5.py`; code in `closure_v5.py`, figures in `plots_src/qv5.py`.

**Source.** Bhowmik, P. 2019, A&A 632, A117 (`Reference Papers/aa34425-18 (1).pdf`; verified against the text): Pearson 0.73 and Spearman 0.81 for |polar-flux asymmetry| at the minimum of cycle n vs |total-area asymmetry| of cycle n+1, n = 15–22 (**flag:** her Fig. 2 caption gives the Spearman coefficient as 0.73); one-factor SFT sensitivities 36 % (tilt scatter), 14 % (+70 % flux), about 3 % (50 % flow); lags of 6–30 months leave the final polar flux unaffected and move the reversals (r = 0.95); data-driven SFT vs MWO, cycles 17–24: polar flux 0.84 (0.94 without the 18/19 minimum), N−S asymmetry 0.91. Earlier work: Goel & Choudhuri (2009, RAA 9, 115) built the poloidal field at each cycle start from the polar faculae record and fed it to a dynamo to study hemispheric asymmetry.

## N1. Pure-lag control

**Figure** `plots/Q2_window_sign_next_cycle/F6_pure_lag_control.png`.
**Number** A 1-yr shift of one hemisphere's emergence changes D by a median **1.18 rms D** (max 2.44; cycles 13–23, N and S, ±1 yr). Label: **lag matters in this SFT; differs from Bhowmik**.
**Sentence** In this SFT the timing of emergence, not only the amount of late flux, sets the polar imbalance; Bhowmik's SFT does not show this.

- a) Symmetrised catalogue, S delayed: final polar-field asymmetry median 12.1 % at 1 yr and 30.6 % at 2.5 yr (hers: unaffected); reversal-time difference +0.65 yr and +1.49 yr (hers about 1 yr at 2.5 yr, r = 0.95 with the lag).
- b) Real catalogue, ±0.5 and ±1 yr, N and S: see the figure; SC12 and SC24 are edges (shifted flux leaves the record) and are excluded from the label.
- **Checks.** Green's-function shortcut vs one direct time-stepped SFT run (S +1 yr): max relative difference 2.1e-04. Mirror check, D(N shifted) = −D(S shifted, mirrored catalogue): 9.1e-14 of rms D.
- **Explanation tested.** Not the decay term: with no decay (τ = 1e6 yr) a 1-yr shift gives 1.10 rms D and a 2.5-yr delay 73.1 % (τ = 10 yr: 1.16). The remaining differences are hers: a 2-D SFT driven by synthetic single cycles, ours: a 1-D SFT read at each minimum of overlapping real cycles. Which of these matters is not established.

## N4. The even-cycle switch on independent records

**Figure** `plots/Q2_window_sign_next_cycle/F7_parity_independent_records.png`.

| Record | Cycles | switch, even | switch, odd | decline sign = A_pf sign | switch agrees with RGO/USAF |
|---|---|---|---|---|---|
| RGO/USAF (x1) | 12–24 | 5/7 | 0/6 | 9/10 | 13/13 |
| RGO/USAF (plain sums) | 12–24 | 5/7 | 0/6 | 9/10 | 13/13 |
| DPD (x1) | 21–23 | 1/1 | 0/2 | 3/3 | 3/3 |
| Temmer/SILSO (plain sums) | 18–24 | 4/4 | 0/3 | 7/7 | 7/7 |

- Windows: maximum t_R+0→+1, decline t_R+3→+4, as in Q2. DPD runs through the project's own annual()/x1() (same latitude normalisation); the sunspot-number records have no latitudes, so they use plain yearly N−S sums; RGO with plain sums gives the same counts as with x1 (the control).
- The records overlap RGO/USAF in time: they test the measurement, not new cycles. DPD covers only three cycles.

## N2. Bridge to her 0.91

**Figure** `plots/Q1_what_sets_polar_asymmetry/S2_bridge_to_bhowmik.png`.

| Step (cumulative) | asymmetry r (n) | odd / even | polar field r (n) |
|---|---|---|---|
| 0 our area-only run, A_pf, SC15-24 | +0.61 (10) | +0.97 / +0.01 | +0.78 (20) |
| 1 cycles 17-24 | +0.90 (8) | +0.97 / +0.71 | +0.85 (16) |
| 2 unnormalised |B_N|-|B_S| | +0.96 (8) | +0.99 / +0.90 | +0.85 (16) |
| 3 her tilt law 0.7 T sqrt|lat| | +0.95 (8) | +0.98 / +0.85 | +0.77 (16) |
| 4 MWO only (minima 17-22; 23-24 not in the project's MWO file) | +0.95 (6) | +1.00 / +0.97 | +0.82 (12) |
| 5 drop the 18/19 minimum (polar flux only) | — | — | +0.84 (10) |

- **Largest move:** restricting to her cycles 17–24 (+0.30 in the asymmetry r). Our area-only SFT already reaches about 0.9 on her cycles; the 0.61 over SC15–24 comes from the two extra cycles SC15–16 (SC16 is the main miss).
- Her tilt law uses the workbook's per-cycle T (MWO, KK, DPD) in place of Jiang et al. (2011) T_n, which the project does not hold; it lowers r slightly. Polar fields: WSO is put on the MWO scale by one factor fitted on minima 20–22 (0.0133); the project's MWO file covers minima 17–22 only.

## N3. Her sensitivities in her metric

**Figure** `plots/Q1_what_sets_polar_asymmetry/S3_bhowmik_sensitivities.png`. Metric: (|B_N|−|B_S|)/|B_sym| at each minimum, cycles 13–23, symmetrised catalogue.

| Change | hers | ours: median (max) |
|---|---|---|
| +70 % southern flux | 14 % | 17.5 % (22.0 %) |
| +50 % southern peak flow | about 3 % | 15.4 % (20.2 %) |
| −50 % southern peak flow | about 3 % | 37.1 % (42.9 %) |
| tilt scatter, 50 realisations (Jiang 2014 σ) | up to 36 % | 3.7 % (19.3 %) |
| southern delay 2.5 yr (N1a) | unaffected | 30.6 % |
| REPORT's 10 % flow difference (v0 16.5 / 13.5 m s⁻¹) | — | 8.9 % (A_pf 0.044) |

- Her +70 % run also lowers the tilt with flux; ours scales flux only. Mirror check: 3.2e-14.

## N5. Positioning against Bhowmik (2019)

| Her claim | Her number and definition | Our number and definition | Verdict | Usable sentence |
|---|---|---|---|---|
| The polar-flux imbalance at minimum sets the next cycle's hemispheric asymmetry | Pearson 0.73 (Spearman 0.81 in the text, 0.73 in the Fig. 2 caption); \|ΔΦ\| at the minimum of cycle n (MWO) vs \|ΔA_total\| of cycle n+1, n = 15–22 | her definition, project MWO and RGO/USAF: +0.59 (n = 8); signed +0.74 (6/8 signs); ours, \|A_pf(n)\| vs \|A(n+1)\|: +0.62 (n = 9) | different in size, consistent in sign | On her pairs and definition the project's polar record gives r = +0.59 (n = 8) against her 0.73, and +0.74 in signed form; the link holds, but its size depends on the polar-field record (Bhowmik 2019; Goel & Choudhuri 2009). |
| A data-driven SFT reproduces the observed polar asymmetry | N−S asymmetry r = 0.91 (n = 8), polar flux 0.84 (0.94 without 18/19), cycles 17–24 | area-only SFT on cycles 17–24: +0.90 normalised, +0.96 unnormalised (n = 8); polar field +0.85 (n = 16); over SC15–24 +0.61 | consistent on her cycles | Over her cycles 17–24 a deterministic SFT reproduces the polar asymmetry at r ≈ 0.9 in both studies; adding SC15–16 lowers it to +0.61. |
| A 0.5–2.5-yr hemispheric lag leaves the final polar flux unaffected and moves the reversals | final flux unaffected; reversal difference vs lag r = 0.95 | symmetrised catalogue: 12.1 % at 1 yr, 30.6 % at 2.5 yr; reversal difference +1.49 yr at 2.5 yr | different (amplitude), consistent (reversal) | In our 1-D SFT a hemispheric lag changes the final polar asymmetry (30.6 % at 2.5 yr), unlike in Bhowmik (2019); both models shift the reversals by about a year. |
| Tilt scatter is the largest source of asymmetry | up to 36 % | median 3.7 %, max 19.3 % | different | With the Jiang et al. (2014) scatter on the real regions, tilt scatter gives a median 3.7 % polar asymmetry, well below her 36 % maximum. |
| More flux in one hemisphere gives a modest asymmetry | +70 % flux: 14 % | 17.5 % | consistent | A 70 % flux excess gives a 17.5 % polar asymmetry, close to her 14 %. |
| The meridional flow hardly matters | ±50 % flow: about 3 % | 15.4 % and 37.1 % | different | Our 1-D SFT is far more sensitive to a hemispheric flow difference than her 2-D SFT. |

## N6. Claim wording

| Claim | Current sentence | Sentence N1 / N4 support |
|---|---|---|
| 8 (parity reading) | In the data that late flux gets the observed sign right in 9 of 10 cycles without any parity factor, and the parity rule at maximum appears to follow from even cycles switching their dominant hemisphere between maximum and decline (5/7 even, 0/6 odd). | The even-cycle switch replicates on the two other records over the cycles they cover: Debrecen DPD (cycles 21–23: even 1/1, odd 0/2) and the Temmer/SILSO hemispheric sunspot numbers (cycles 18–24: even 4/4, odd 0/3), and the declining-phase sign matches the observed A_pf in 3/3 and 7/7 cycles. These records overlap RGO/USAF in time, so they confirm the measurement, not new cycles; the reading stays exploratory. |
| 10 (timing lag) | The SFT's polar imbalance is essentially "which hemisphere peaks later" (r = 0.87, n = 13). | In the SFT the emergence timing itself moves the polar imbalance: shifting one hemisphere's emergence by 1 yr changes D by a median 1.18 rms D (cycles 13–23), so the 0.87 link is not only a proxy for late flux (model result). This differs from Bhowmik (2019), whose lags of 0.5–2.5 yr left the final polar flux unaffected; our decay term is not the reason (without decay the effect is larger), and the reason is not established. |

## Numbers that differ from REPORT.md

- N3 REPORT 10% flow A_pf: 0.044 against REPORT's 0.05. REPORT.md states the 10 % flow-difference offset as "≈ 0.05"; the stored value is 0.044.

## Steps changed or skipped

- N2 step 3: Jiang et al. (2011) T_n are not in the project; the workbook's per-cycle T was used instead (flagged).
- N2 step 4: her MWO dataset is only partly available (project MWO file: minima 15–22); minima 23–24 are skipped in that step.
- N5: her Spearman coefficient is not recomputed (rank tests are outside this project's plain statistics).

## Files

`closure_v5.py` (n1, n1x, n4, n2, n3, n3b, n5); `out/v5/` (JSON, progress.md); `plots_src/qv5.py`, `make_v5.py`, `write_v5.py`; figures `plots/Q2_window_sign_next_cycle/F6_*`, `F7_*`, `plots/Q1_what_sets_polar_asymmetry/S2_*`, `S3_*`; checks `plots/number_checks_v5.md`.
