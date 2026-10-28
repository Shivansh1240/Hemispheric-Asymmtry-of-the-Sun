# v8: everything redone in Bhowmik's models

Earlier reports and numbers are unchanged. Here every simulation of v3-v7 is repeated with Bhowmik's own settings: her SFT (Bhowmik 2019, Sects. 3.1 and 6; the data-driven SFT of Bhowmik & Nandy 2018) and her dynamo (Passos et al. 2014 on SURYA with pumping, v7). Code: `sft_bhowmik.py`, `sft_bhowmik_checks.py`, `surya/surya_bhowmik_ext.f`, `surya_v8.py`; numbers `out/v8_sft/`, `out/v8_surya/`; figures `plots_v8/` (`plots_src/qv8_sft.py`, `qv8_dyn.py`, `make_v8.py`).

## 1. Settings

| Item | our SFT (v3-v5) | Bhowmik | v8 |
|---|---|---|---|
| transport | eta 250 km2/s, v0 15 m/s, lambda0 75 deg | same | same |
| flux | 7e19 A Mx | same | same |
| decay | tau = 5 yr | none | none (150-yr Green's function) |
| tilt | T_n(cat) sqrt(\|lat\|), project T_n 0.51-0.95 | 0.7 T_n sqrt(\|lat\|), T_n = 1.73 - 0.0035 S_n (Jiang et al. 2011 Eq. 15) | hers; S_n from Jiang Table 1, S_24 76.7, S_25 105.2 from SILSO V2 scaled to S_23 |
| tilt scatter | Jiang et al. 2014 | same (synthetic runs); none in her data-driven run | ensembles around her mean tilt |
| run | continuous from SC12 | continuous from SC15, dipole start | continuous from SC12 (start effect checked) |

Dynamo: v7 settings (Passos et al. 2014 on the SURYA guide standard, 100 % / 200 % uniform noise renewed every 0.5 yr, pumping 3.6 m/s). The earlier experiments' knobs are added in `surya_bhowmik_ext.f`: noise coherence tau_c, N-S noise correlation rho (S draws set equal to N with probability rho), equatorial latitudinal diffusivity factor kappa (flux form, A and B), flow speed, and a constant N-S bias in alpha_BL, the flow or eta.

## 2. Checks

- Our catalogue tilt is exactly T_n sqrt|lat| per cycle (constant to 4 decimals within each cycle), with project T_n 0.51-0.95.
- Without decay and with our T_n the polar field does not reverse at the SC19, SC21 and SC23 minima, so her T_n is required.
- Her 2019 hemispheric law T = 1.71 - 0.12 m F gives T < 0 for every observed cycle (calibrated on synthetic cycles about 20 times weaker), so it is not used.
- With her T_n every minimum reverses.
- Start effect. Without decay the 1-D SFT integrates: flux from earlier cycles ends up on the caps as an exactly antisymmetric dipole. That dipole adds no N-S asymmetry (earlier-cycle share of D at most 8e-05) but sets |B_N| + |B_S|. In our run the carried dipole is 33.1 G, and all poles reverse for 17.5-34.0 G. Over that range r(D, observed A_pf) = 0.66 with 7/10 signs, while r(A_pf, observed) = 0.34-0.72 and the lag-1 r of the modelled A_pf = −0.15 to +0.26. Normalised A_pf results of her SFT therefore depend on a start choice the data do not fix; D-based results do not.
- Dynamo: with the new knobs off the extended code is bit-identical to v7 (2-yr run, `out/v8_surya/_check`).
- Dynamo: the flux-form latitudinal operator matches the analytic one to 0.3 % (grid error).
- Dynamo: every new knob is mirror-symmetric to about 1e-6 (`out/v8_surya/validation.json`).

## 3. SFT results: before (our SFT) and after (Bhowmik's SFT)

| Figure | Result | our SFT | Bhowmik's SFT | Verdict |
|---|---|---|---|---|
| Q1 F3/F4 | area-only SFT vs observed A_pf, SC15-24 (odd, even) | 0.61 (0.97, 0.01) | 0.66 (0.95, 0.01); unnormalised D 0.66 | holds: agreement only in odd cycles |
| Q1 F3 | amount only / placement only | 0.20 / 0.60 | 0.50 / 0.58 | changes: the amount of flux carries agreement too; 'placement, not amount' no longer holds |
| Q1 F3 | rms modelled A_pf (observed 0.102) | 0.044 | 0.224 | changes: amplitude too large (carried dipole, see checks) |
| Q1 S1 | r(model, observed) over transport configurations | 0.56-0.63 (27 configs) | 0.51-0.73 (9 configs, no decay) | holds |
| Q2 F1 | imbalance built at t_R + (16-84 %) | 2.3 to 5.1 yr | 2.3 to 5.0 yr | holds; the imbalance kernel falls off without any decay term |
| Q2 F1 | N pole built at t_R + (16-84 %) | 0.3 to 3.6 yr | 3.6 to 4.8 yr | shifts later |
| Q2 F5 | r(timing lag, model A_pf), n = 13 | 0.87 | 0.53 | weaker |
| Q3 F4 | SC14 negative control: model A_pf (MWO -0.147) | -0.032 | -0.019 | holds (stays south) |
| Q3 F4 | SC21 model A_pf (observed -0.077) | -0.059 | -0.557 | same sign; amplitude inflated |
| Q3 F5 | share of the imbalance from earlier cycles | about 1e-5 | at most 8e-05, also without decay | holds: no inter-cycle memory of the imbalance |
| Q3 F5 | lag-1 r of modelled A_pf | 0.37-0.44 (decay 2-17 yr) | −0.04 without decay (−0.15 to +0.26 over carried dipoles) | changes: no A_pf persistence in her SFT |
| Q4 F1 | r(member, obs): real / sizes re-dealt / timing re-dealt | 0.54 / 0.51 / 0.11 | 0.55 / 0.47 / 0.20 | holds: timing/latitude matters most |
| Q4 F1 | predetermined share of the variance | 0.77 / 0.67 / 0.28 | 0.10 / 0.09 / 0.03 | changes: tilt scatter dominates every member |
| v5 N1 | 1-yr shift of one hemisphere: median \|dD\| / rms D | 1.18 | 1.08 | holds: still unlike her 'unaffected' |
| v5 N2 | her cycles 17-24, unnormalised \|B_N\|-\|B_S\| (hers 0.91) | 0.96 | 0.94 | holds |
| v5 N2 | her cycles 17-24, polar \|B\| (hers 0.84) | 0.85 | 0.49 | worse (field level set by the carried dipole) |
| v5 N3 | +70 % S flux / +-50 % S flow / tilt scatter (hers 14 / 3 / up to 36 %) | 18 / 15-37 / 4 % | 16 / 19-40 / 4 % | holds |

## 4. Dynamo results (Bhowmik's dynamo on SURYA)

Sun: lag-1 r of A = +0.57, longest run 7, rms A_pf 0.10. Earlier model in brackets (same setting; flow ×2/3 = 10 of 15 m/s there).

| Setting | runs | cycles | pooled lag-1 r of A | mean longest run / 13 | rms A_pf | Sun-like | timing-lag lag-1 r |
|---|---|---|---|---|---|---|---|
| kappa=1 | 6 | 270 | **+0.45** (+0.02) | 5.5 (4.1) | 0.27 (0.44) | 78 % | +0.14 (+0.05) |
| rho=0.8 | 3 | 107 | **+0.50** (−0.05) | 4.0 (4.2) | 0.12 (0.27) | 86 % | +0.13 (+0.01) |
| v0=20 m/s | 1 | 27 | **+0.29** (−0.01) | 4.9 (4.3) | 0.57 (0.41) | 33 % | +0.66 (+0.04) |
| rho=0.5 | 1 | 35 | **+0.37** (−0.14) | 4.8 (3.9) | 0.25 (0.41) | 89 % | +0.20 (+0.02) |

Detector check (cycle minima from the 2-yr smoothed eruption series, at least 0.7 period apart, instead of the Q4 detector; in SURYA eruptions never stop, so the Q4 detector puts 3-7 of 35-45 boundaries per run early or late):

| Setting | lag-1 r of A (Q4 detector → check) | longest / 13 | rms A_pf | Sun-like |
|---|---|---|---|---|
| kappa=1 | +0.45 → +0.46 | 5.5 → 5.2 | 0.27 → 0.31 | 78 → 86 % |
| rho=0.8 | +0.50 → +0.51 | 4.0 → 4.8 | 0.12 → 0.13 | 86 → 94 % |
| v0=20 m/s | +0.29 → +0.35 | 4.9 → 4.9 | 0.57 → 0.59 | 33 → 36 % |
| rho=0.5 | +0.37 → +0.34 | 4.8 → 5.1 | 0.25 → 0.26 | 89 → 95 % |

**Not run** (batch stopped on 2026-10-27 at the user's request, laptop overheating): tau_c 3 and 5 yr at production noise, kappa 0.3 and 3, flow seeds 2-3, the T3 production runs (T3 below uses the calibration runs) and all of T4.


### T3, noise coherence time (two calibration runs per tau_c, noise scaled by a first-guess factor)

| tau_c (yr) | scale s | rms A | pooled lag-1 r | longest / 13 | rms A_pf | Sun-like | earlier model lag-1 |
|---|---|---|---|---|---|---|---|
| 0.5 (Bhowmik production, not re-tuned) | 1.000 | 0.18 | **+0.45** | 5.5 | 0.29 | 78 % | — |
| 0.5 | 0.550 | 0.10 | **+0.38** | 4.1 | 0.17 | 88 % | — |
| 1 | 0.400 | 0.07 | **+0.20** | 3.7 | 0.17 | 74 % | +0.01 |
| 3 | 0.250 | 0.10 | **+0.31** | 5.5 | 0.15 | 97 % | +0.05 |
| 5 | 0.190 | 0.09 | **+0.28** | 4.0 | 0.11 | 100 % | +0.12 |
| 10 | 0.140 | 0.07 | **+0.18** | 5.0 | 0.08 | 100 % | +0.19 |
| 30 | 0.110 | 0.08 | **+0.72** | 5.2 | 0.08 | 100 % | +0.58 |

### T4, constant hemispheric bias (production noise)

| bias | pooled lag-1 r | longest / 13 | A > 0 | A_pf > 0 | mean A |
|---|---|---|---|---|---|
| reference (delta 0) | +0.45 | 5.5 | 46 % | 51 % | -0.023 |

## 5. Bottom line

- Her SFT: the observational links and the SFT's odd-cycle-only agreement hold. 'Placement, not amount' does not hold (amount alone now reaches 0.50). Normalised A_pf amplitudes depend on the carried dipole, so quote D-based agreement (0.66).
- Her SFT carries no inter-cycle memory of the imbalance, with or without decay.
- Her dynamo carries hemispheric dominance from cycle to cycle in every setting run (pooled lag-1 r +0.29 to +0.50; Sun +0.57). Our earlier dynamo gave −0.14 to +0.02. The lack of persistence was a property of our model, not of flux-transport dynamos.
- Even in her dynamo the dominant hemisphere's runs are near independent-cycle length (4.0-5.5 per 13 against the Sun's 7). Correlated amplitudes are reproduced; the Sun's long same-sign runs are not.
- Dynamo statements are model results from a stopped batch (see 'Not run'). T3 and T4 conclusions from v4 are not replaced here.

