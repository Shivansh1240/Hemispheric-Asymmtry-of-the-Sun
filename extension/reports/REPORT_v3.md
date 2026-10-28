# Hemispheric asymmetry: analysis v3 (2026-10-24)

Gupta, Mishra & Karak, RAA. Plain statistics only: Pearson r with n, sign counts and run lengths. Ensembles are shown as 16–84 % bands. There are no p-values.

**Data and engine.** Everything runs on the project's verified layer through `master_code/asymcore.py`:
- `RGO_SOON_corrected.txt` (md5 checked);
- the adopted A_pf (MWO faculae SC15–19, WSO SC20–24, ±1.6 yr about the minimum ending cycle n);
- the `bn_sft` 1-D SFT: η = 250 km² s⁻¹, v0 = 15 m s⁻¹, λ0 = 75°, τ = 5 yr; cap ≥ 70°, area-weighted.

**Checks passed before any number was used** (`v0_checks.py`):
- The new Green's-function readout reproduces the project's deterministic run to 3e-16.
- Decay-time rescaling of the kernel equals a direct build to 7e-15.
- **Mirror check** (N↔S swapped in the catalogue, in the daily areas, and in η_N↔η_S): the model returns N' = −S, S' = −N to 4e-15 (1e-14 with hemispheric transport). The data return X1' = −X1 and A' = −A exactly. Both headline r values are unchanged. The parity sign therefore enters only where it should.

---

## 0. First job: C1 (r ≈ −0.98, sign reversal) versus unsigned X1 (r ≈ +0.77)

**What C1 is.** No stored file holds "−0.98", so I reproduced it.
- It is the same-cycle relation written against the *signed* polar field S = (B_N+B_S)/(|B_N|+|B_S|): r(X1, S) = **−0.94 (n = 10), with all 10 signs reversed**.
- In first differences it is **−0.98 (n = 9)**.
- Because (−1)^(n+1) = −sgn B_N at all 10 minima (checked), X̂1·A_pf ≡ −X1·S. The −0.94/−0.98 and the +0.93 are **one relation**: the "sign reversal" is the polarity of B, not physics.

**Why it is not in conflict with +0.77.** The +0.77 is a different relation:
- the target is the next cycle's area, not the polar field;
- the window is t_R−1→+4, not t_R+0→+1;
- it is unsigned.

What connects them is the magnitude:
- |A_pf(n)| vs |A(n+1)|: r = 0.62 (n = 9).
- |X̂1(t_R+0→+1)| vs |A(n+1)|: r = 0.63 (n = 12).
- Persistence of magnitude, |A(n)| vs |A(n+1)|, is only r = 0.09 (n = 12).
- The sign of A(n+1) follows persistence (10/12), not alternation (2/12) or the n+2 rule (7/11).

**Primary result: C1, quoted once, as r(X̂1, A_pf) = +0.93 (n = 10, 10/10 signs).** Why:
1. It is the measured first link of the hypothesis chain.
2. The signed-field −0.94 is the same number, and the −0.98 is a first-difference value on n = 9 that should not be quoted as a separate or stronger result.
3. The +0.77 adds no sign information beyond persistence. A(n) alone gives 10/12 signs, the same as X1, and r = 0.57.

Two caveats travel with C1:
- Its magnitude is fragile to the published A_pf errors (Section E).
- Its sign needs the parity label. Section Q2 now offers a physical reading of that label.

Figure: `figs/Q0_reconcile.png`.

---

## Q1. What sets the polar-field asymmetry?

**Figure** `figs/Q1_what_sets_Apf.png`.
**Number** Area only: r(model, observed A_pf) = **0.61 (n = 10)**. The full model gives 0.55 (16–84 %: 0.42–0.66).
**Sentence** In the SFT, where and when the flux emerges recovers all of the agreement the full model has with the observed A_pf (r = 0.61 against 0.55) and 75 % of the full model's variance; the amount of area alone does not (r = 0.20), and the agreement is entirely in odd cycles (r = 0.97; even cycles 0.01).

Each run starts from a symmetrised catalogue (every region split into half-flux copies at ±|lat|), which gives A_pf = 2e-15. One difference at a time is then switched on.

| Experiment | r with observed A_pf (n = 10) | rms A_pf (observed 0.102) |
|---|---|---|
| (a) area only: real emergence, Joy tilts, symmetric transport | **0.61** (odd 0.97, even 0.01; signs 7/10) | 0.044 |
| (a1) amount only: hemispheric flux share, symmetric placement | 0.20 | 0.021 |
| (a2) placement only: real latitudes and timing, flux equalised | **0.60** | 0.039 |
| (b) tilt scatter only (100 members) | 0.10 (−0.32 to 0.42) | 0.017 |
| (c) transport only: η or v0 differing by ±10–20 % | ±0.18 to ±0.24 (a constant offset of the same sign every cycle) | 0.027–0.095 |
| area + tilt scatter (25 members) | 0.55 (0.44–0.66) | 0.052 |
| (d) all: area + scatter + 8 transport asymmetries (200 members) | 0.55 (0.42–0.66) | 0.072 |

- **X̂1·T_lin** (cycle Joy slope from the tilt workbook): r = 0.89 against 0.93 for X̂1 alone. No gain. T_lin vs |A_pf|: r = −0.83 (n = 10).
- **Algebraic hemispheric dipole with η_N ≠ η_S** (λ_R from each hemisphere's η): r peaks at δ = 0 (0.45) and falls to 0.30 at δ = ±0.6. A diffusivity difference does not help.
- **Does not hold:** a hemispheric transport difference produces a same-sign offset every cycle. It cannot produce the cycle-to-cycle signal. A 10 % flow difference does, however, give an offset (≈ 0.05) comparable to the observed amplitude.

---

## Q2. Which emergence epoch feeds each pole, and the next-cycle asymmetry

**Figures**
- `figs/P3_green_epoch_latitude.png` (Green's function);
- `figs/P2_window_heatmaps.png` (whole window table);
- `figs/Q2_size_latitude_table.png` (whole size × latitude table);
- `figs/Q2b_phase_sign.png` (sign model).

**Number** At the physics window t_R+3→+4, the unsigned X1 vs A_pf is **r = 0.68 (n = 10), with 9/10 signs and no parity factor**. The parity-signed X̂1 at t_R+0→+1 gives 0.93.
**Sentence** The SFT builds each minimum's polar imbalance from flux emerging 2.3–5.1 yr after t_R. In the data that late flux gets the observed sign right in 9 of 10 cycles without any parity factor, and the parity rule at maximum appears to follow from even cycles switching their dominant hemisphere between maximum and decline (5/7 even, 0/6 odd).

**Green's function** (exact per-region attribution, cycles 12–24):
- **Each pole is built by both hemispheres.** 60 % of the N pole comes from northern regions and 40 % from southern ones, through cross-equatorial leading flux.
- **The poles are built at t_R+0.3 to +3.6 yr** (16–84 %).
- **The imbalance |B_N| − |B_S| is built later, at t_R+2.3 to +5.1 yr (median +3.7).** A region's contribution to the imbalance decays within about 5–8 yr. By the next minimum, low-latitude flux feeds both poles equally. The kernel ratio (imbalance per own-pole contribution) at 15° is 1.00 at 2 yr, 0.27 at 5 yr and 0.008 at 8 yr.

**Physics window.** It is chosen from the model only: the shortest window whose emergence carries the model's own imbalance at r ≥ 0.9 over 13 cycles. The result is **t_R+3→+4**.
- In the data's whole scan (104 windows, `P2` panels a–b), this window ranks 74th for the parity-signed A_pf relation (0.53) and 51st for the next cycle (0.67).
- **Does not hold:** the physics window does not reproduce the paper's best parity-signed correlations (0.93 and 0.77).

**Size × latitude, the whole table** (`Q2_size_latitude_table.png`).
- At t_R+0→+1 the A_pf relation is carried by the large regions: ≥1000 MSH gives 0.97, ≥500 MSH 0.94, and <100 MSH only 0.58.
- The next-cycle relation is spread over all size classes (0.65–0.85 at all latitudes). No cell stands out. Cells where the predictor is zero in more than two cycles are greyed as degenerate.

**Effective sign, modelled explicitly** (`Q2b`):
- At maximum, s_eff = sgn(A_pf·X1) equals the parity in 10/10 cycles. The polarity-blind SFT matches it in 7/10.
- At the physics window, s_eff = +1 in 9/10 cycles. The one failure is **SC16, which is also the SFT's sign failure**.
- Over SC12–24 the dominant hemisphere at maximum and in the decline differ in **5/7 even and 0/6 odd cycles**. SC12–14 have no polar data, so they are an out-of-sample check: 2 of 3 cycles are consistent.
- A reading that follows: the parity factor in C1 may not be a polarity effect. The polar imbalance follows the declining-phase emergence (as the SFT says), and even cycles tend to swap their dominant hemisphere after maximum.
- Caveats: this is exploratory, the r is lower (0.68 against 0.93; odd 0.95, even 0.38), and it holds only for decline windows starting at t_R+3 or later. Starting at +2 gives 6/10 signs.

**Next-cycle sign rules** (hits / n):

| Rule | Hits |
|---|---|
| persistence, sgn A(n) | 10/12 |
| sgn X1(t_R−1→+4) | 10/12 |
| sgn X1, physics window | 9/12 |
| n+2 rule, sgn A(n−1) | 7/11 |
| alternation | 2/12 |
| sgn A_pf(n) | 7/9 |

A two-part model, sign from A(n) and magnitude from |A_pf(n)|, gives r = 0.80 (n = 9), against 0.79 for the unsigned X1 on the same 9 cycles.

---

## Q3. Persistence and flips

**Figures**
- `figs/Q3c_persistence_tau.png` (runs, τ sweep);
- `figs/Q3a_flip_butterflies.png`;
- `figs/Q3b_sft_and_datacheck.png`;
- `figs/P4_cumulative_NS.png`;
- `figs/P5_return_map.png`.

**Number** Dominant-hemisphere run lengths SC12–24 are **[2, 7, 4], i.e. 3 runs**. Independent cycles would give 7.0 runs on average, with a longest run of 4.0.
**Sentence** Area dominance is persistent and the persistence lives in the regions ≥500 MSH, but the polar asymmetry itself is not persistent, there is no n+2 seeding, and the SFT carries no cycle-to-cycle memory of the imbalance at any decay time from 2 to 17 yr.

**Persistence by class**

| Series | Run lengths | Independent reference |
|---|---|---|
| ≥500 MSH | [2, 7, 4] | — |
| low-latitude ≥500 MSH | [2, 7, 4] | — |
| 100–499 MSH | 6 runs | — |
| <100 MSH | 6 runs | — |
| A_pf (SC15–24) | 7 runs in 10 | 5.5 |

**n vs n+2**
- r(A(n), A(n+1)) = +0.57 (n = 12).
- r(A(n), A(n+2)) = −0.28 (n = 11).
- The residual of A(n+1) after A(n), against A(n−1): −0.58 (n = 11).
- r(A_pf(n), A(n+2)) = +0.39 (n = 9), against r(A_pf(n), A(n+1)) = +0.68.
- **No n+2 seeding.** The lag-2 term is negative.

**Flips and near-flips** (whole-cycle A: all / low-latitude ≥500 MSH / rest)

| Cycle | All | Low-lat ≥500 | Rest | Effect of the low-latitude large regions |
|---|---|---|---|---|
| SC14 flip | +0.025 | +0.049 | +0.002 | tip it |
| SC21 flip | −0.013 | −0.038 | +0.003 | tip it |
| SC17 near flip | +0.004 | +0.076 | −0.053 | **hold the regime** |
| SC24 near flip | −0.007 | −0.123 | +0.034 | **hold the regime** |
| SC18 | +0.028 | +0.041 | +0.018 | — |

These regions are the decider in any near-balanced cycle, in both directions. That is not specific to flips.

**Region bootstrap** (1000 draws). The flip cycles' signs are not established:

| Cycle | Draws with the observed sign |
|---|---|
| SC21 | 64 % |
| SC14 | 69 % |
| SC17 | 56 % |
| SC24 | 55 % |
| SC13 | 99 % |
| SC20 | 100 % |

The two "flips" are moves from a firm regime into a balanced cycle.

**DATA CHECK: the SC21 flip on independent records**

| Record | SC20 | SC21 | Low-lat ≥500, SC21 |
|---|---|---|---|
| RGO / USAF-NOAA | +0.170 | −0.013 | −0.038 |
| Debrecen DPD (own plates, 1974–2013) | (partial) | **−0.014** | **−0.007** |
| Temmer et al. 2006 hemispheric SN | +0.124 | **−0.018** | — |

- DPD vs RGO_SOON yearly N−S: r = 0.99 (n = 40 yr; 1977–86 alone 0.99).
- **The sign of SC21 replicates in all three records.**
- **The "carried by the low-latitude large regions" part does not replicate in DPD.** There the class is −0.007 and the rest −0.020. It stays a catalogue-dependent statement, not a headline.

**SFT, SC21 vs SC14** (`Q3b` a)

| | Model | Observed |
|---|---|---|
| SC20 polar imbalance | S (A_pf −0.035) | S (−0.057) |
| SC21 polar imbalance | S (A_pf −0.059; 100 % of tilt members) | S (−0.077) |
| SC14 polar imbalance | S (−0.032) | S (MWO −0.147 ± 0.143) |

- In SC21 the area flip and the polar imbalance agree (both S). In SC14 the area flipped to N but the polar imbalance stayed S, in both model and data. **SC14 is the negative control**: the area flip did not reach the polar field.
- In SC21 the low-latitude ≥500 MSH regions supply 54 % of the model D. **Without them D(21) stays southern (−0.61 of −1.33), so they are not needed for the sign.**

**Decay-time sweep, τ = 2–17 yr** (exact kernel rescaling)
- Lag-1 r of the modelled A_pf stays at 0.37–0.44 (0.21–0.29 for the tilt-ensemble medians) at every τ.
- r(model, observed) stays at 0.60–0.66.
- **Flux from earlier cycles adds large but equal fields to both poles.** Its share of the imbalance D is about 1e-5 at every τ: by the next minimum the old flux is a pure dipole.
- **The 1-D SFT cannot carry multi-cycle persistence of the imbalance at any τ.** The model's lag-1 memory comes from the emergence pattern itself. (This points to Q4: the 2-D flux-transport dynamo.)

---

## E. Ensemble verification

**Figure** `figs/ENS_headlines.png`.

**H1, X̂1 vs A_pf, 0.93 (n = 10)**
- Area noise (σ_ln = 0.20 per daily row, 100 members): 0.92–0.93, i.e. it does not move.
- A_pf redrawn within its published 1σ: **0.74 (0.55–0.85)**. H1 is fragile to polar-field errors of the published size. This ensemble adds a second error realisation, so it is a sensitivity, not an unbiased estimate.

**H2, X1 vs A(n+1), 0.77 (n = 12)**
- Area noise: 0.76–0.78. Robust.

**H3, the Cycle 21 flip**
- Area noise: A(21) < 0 in 100 % of members, A(20) > 0 in 100 %.
- Region bootstrap: A(21) < 0 in only 64 % of draws (see Q3).
- Model D(21) < 0 in 100 % of the emergence members (tilt scatter + area noise). Model D(20) > 0 in only 5 %, so the model has no polar "flip" between 20 and 21; neither do the observations, where both are S.

**Longitudes.** Neither X1 nor the 1-D longitude-averaged SFT uses longitude, so randomising it is exactly a no-op. A test needs a 2-D model.

**Transport ensemble** (η 150/250/400 × v0 10/15/20 × τ 3/5/10, 27 configs). Range over the configs:

| Quantity | Range |
|---|---|
| r(model, observed A_pf) | 0.56–0.63 |
| r(X̂1, model A_pf) | 0.64–0.72 |
| D(21) southern | 27/27 configs; 100 % of the tilt members in every config |
| low-latitude ≥500 MSH share of D | 0.55–0.64 |
| imbalance epoch | t_R+3.0 (v0 = 10) to t_R+4.4 (v0 = 20) |
| r(\|model A_pf(n)\|, \|A(n+1)\|) | −0.29 to +0.57 (median 0.09) |

The ensemble band is drawn on `ENS_headlines` panel c.

**Model results vs Sun results**

| Conclusion | Holds across the 27 transport configs? | Status |
|---|---|---|
| X̂1 → A_pf 0.93; X1 → A(n+1) 0.77 | (data only) | **Sun result**, though H1 is fragile to A_pf errors |
| Area dominance persists (3 runs vs 7 independent); A_pf does not | (data only) | **Sun result** |
| SC21's southern sign in three independent records | (data only) | **Sun result**, but only 64 % in the region bootstrap |
| Decline-phase X1 gives the A_pf sign 9/10 without parity; even cycles switch (5/7 vs 0/6) | the window start t_R+3 to +4.4 holds in every config | **Sun result in sign** (data); window **set by a model-class property** |
| Imbalance built at t_R+2.3 to +5.1 yr; no memory beyond about 5–8 yr; older cycles contribute about 1e-5 of D | yes, all configs and all τ = 2–17 | **model-class result** (1-D SFT); not observable directly |
| Low-latitude ≥500 MSH supply about 57–60 % of D | yes (0.55–0.64) | **model result**, robust within the model |
| Model D(21) southern | yes (27/27) | **model result**, which agrees with the observed A_pf(21) |
| Exact physics-window boundaries (t_R+3 → +4) | **no**: they shift by about ±0.7 yr with v0 | **model result for one parameter set**; quote it as "declining phase" |
| \|model A_pf\| sets \|A(n+1)\| | **no** (−0.29 to +0.57) | not supported in the model; only in the data (0.62, n = 9) |

---

## Exploratory track A–E

Figures: `figs/EXP_A_B_D_E.png` and `figs/P6_concentration.png`.

**A. Counterfactual: remove, or move to 35°, the largest low-latitude region of each cycle.**
- The modelled D changes by a median of 0.05 % (removed) and 0.03 % (moved). The maximum is SC18: −23 % removed, +20 % moved.
- **No sign flips in 13 cycles.**
- The largest regions emerge near maximum, and by the minimum their imbalance contribution has decayed. The top individual contributors emerge at t_R+3.7–4.6 yr.
- **Dead end for "one rogue giant decides the cycle"**, as far as the SFT can tell.

**B. X_G**, the Green's-function kernel × flux × separation × sin(tilt), signed by hemisphere.
- It is effectively the linearised SFT: r(X_G, model D) = 0.98 (n = 13).
- r(X_G, A_pf) = 0.56 (n = 10; odd 0.95, even 0.21). The parity-signed version gives 0.48, against 0.93 for X̂1.
- Next cycle: 0.67, against 0.77 for X1.
- **The physics weighting does not improve the predictor. Dead end** as a better predictor. It is useful only as the reason for the window result in Q2.

**C. Sparse attribution.**
- A median of 17 regions (range 1–76), about 0.15 % of the regions reaching a minimum, reach each cycle's |D|.
- But the cumulative positive contributions overshoot |D| by 1.7–6.8×. **|D| is a small residual of large opposing contributions.** That is why single-region counterfactuals barely move it, and why the concentration curve (`P6`) is not a simple "few regions do it".

**D. Memory model** A(n+1) = c + a A(n) + b A(n−1), n = 11.
- a = +0.77, b = **−0.59**, r(fit) = 0.76.
- The n+2 effect does **not** appear; the lag-2 term is negative (anti-persistent).
- The longest run of a simulated 13-cycle series averages 5.3 (AR(2)) and 6.2 (AR(1)). The observed run is 7; for independent fair-coin cycles it is 4.0.
- **A slow mode plus noise kicks is consistent with the record. It is not demanded by it.**

**E. Hemispheric timing**, the N−S lag of the area-weighted centroid per cycle.
- r(lag, A(n)) = 0.32 (n = 13).
- **r(lag, modelled A_pf) = 0.87 (n = 13), and r(lag, decline-phase X1) = 0.87.** The SFT's polar imbalance is essentially "which hemisphere peaks later" (`out/q2b_phase_sign.json`).
- r(lag(n), A(n+1)) = 0.63 (n = 12; 9/12 signs; 0.53 without SC19/20).

**Looks genuinely new**
1. **(Q2/E, the strongest) The parity sign as a phase effect.** The polar imbalance follows the declining-phase emergence. Even cycles swap their dominant hemisphere after maximum in 5/7 cases (odd 0/6), which turns the maximum-phase X1 into a parity-signed predictor. This would give the "label" a mechanism consistent with SFT physics. It needs a check on more cycles and on DPD positions.
2. **(E) The N−S timing lag is what the SFT's imbalance measures** (0.87). The later-peaking hemisphere tends to dominate the next cycle (0.63). This is a candidate observable link from polar imbalance to the next cycle.
3. **(C/Q3) The imbalance is a small residual of large opposing contributions** and has no inter-cycle memory in 1-D SFT. Each pole is fed 40 % by the opposite hemisphere.

**Dead ends**
- A: single-giant counterfactuals.
- B: X_G as a better predictor.
- D: the n+2 seeding, which is negative.
- Hemispheric diffusivity asymmetry.
- Tilt-weighting (X̂1·T_lin).

---

## Q4. Stochastic or predetermined?

### Q4a. Randomised emergence with the hemispheric totals kept (SFT, 100 members each)

**Figure** `figs/Q4a_randomised_emergence.png`.
**Number** Re-dealing *when and where* flux emerges, with each hemisphere keeping its own regions and total flux, gives r(member, observed A_pf) = **0.11 (16–84 %: −0.12 to 0.36; n = 10)**. With the real emergence pattern it is 0.54.
**Sentence** The hemispheric totals do not predetermine the polar asymmetry; the hemispheric pattern of emergence timing and latitude does, and which particular regions are large hardly matters.

| Experiment | r(member, observed A_pf) | Predetermined share of variance | Cycles keeping the observed sign in ≥84 % of members |
|---|---|---|---|
| real catalogue, tilt scatter only | 0.54 (0.41–0.64) | 0.77 | 6/10 |
| R1: sizes re-dealt within each hemisphere's own slots | 0.51 (0.36–0.60) | 0.67 | 5/10 |
| R2: slots (time, \|lat\|) re-dealt across hemispheres, totals kept | **0.11 (−0.12 to 0.36)** | 0.28 | 1/10 |

"Predetermined share" means the variance of the ensemble mean across cycles over the total variance. Given the real emergence pattern, about three quarters of the modelled A_pf is fixed and about a quarter is tilt noise.

`StochasticJoyEnsemble` (Downloads/files (24)/sft_extensions_r1.py) drives a synthetic Hathaway source, not the catalogue. The same idea, per-region tilt scatter, is applied here to the real regions.

### Q4b. The key correlations with and without the largest regions

**Figure** `figs/Q4b_without_largest.png`.
**Number** Removing the single largest region of each cycle (15 regions, 1.1 % of the area) changes the dominant-hemisphere runs from **[2, 7, 4] to [2, 3, 1, 3, 3, 1]**. It leaves X̂1 → A_pf at 0.93 and X1 → A(n+1) at 0.75.
**Sentence** The correlations are not driven by a handful of rogue regions, but the persistence record and the SC21 sign are.

**Data** (top k regions per cycle, by peak area, removed from the daily rows):

| k removed | X̂1 → A_pf (n = 10) | X1 → A(n+1) (n = 12) | decline X1 → A_pf, unsigned | lag-1 r of A | runs | A(21) |
|---|---|---|---|---|---|---|
| 0 | 0.93 | 0.77 | 0.68 (9/10) | 0.57 | [2, 7, 4] | −0.013 |
| 1 | 0.93 | 0.75 | 0.68 (9/10) | 0.56 | [2, 3, 1, 3, 3, 1] | −0.020 |
| 5 | 0.93 | 0.70 | 0.69 (9/10) | 0.49 | 6 runs | −0.008 |
| 10 (7 % of area) | 0.86 | 0.69 | 0.71 (9/10) | 0.51 | 6 runs | **+0.003** |
| 20 (11 % of area) | 0.83 | 0.68 | 0.72 (9/10) | 0.48 | 6 runs | **+0.015** |
| *only* the top 10 | 0.82 | 0.07 | 0.78 (4/10) | 0.05 | 8 runs | −0.300 |

- The k = 1 change in the runs comes from SC17 and SC24, the two near-flip cycles, changing sign.
- **Rogue-region dependence, said plainly:**
  - **The 7-cycle northern run depends on single regions** in near-balanced cycles.
  - **SC21's southern sign depends on its 10 largest regions.**
  - The three correlations do not depend on a handful.

**SFT.** Removing the top 1–20 regions per cycle leaves r(model, observed) at 0.59–0.65. It leaves r with the full model at 0.97–1.00. In the model the largest regions are not special (consistent with Q3 and exploratory A).

### Q4c. The 2-D flux-transport dynamo with hemispheric coupling

**Figures**
- `figs/Q4c_dynamo_coupling.png`: one run per setting, with a butterfly;
- `figs/Q4d_dynamo_seed_ensembles.png`: seed ranges and pooled values.

**Number** Over 5 seeds (339–342 cycle pairs per setting), pooled lag-1 r of the cycle asymmetry is **+0.02 at κ = 1 and −0.05 at ρ = 0.8**. The Sun gives +0.57 (n = 12). The best of all settings is +0.13 (forcing coherence 5 yr, n = 203).
**Sentence** A stochastically forced 2-D Babcock–Leighton dynamo renews its hemispheric asymmetry every cycle under every coupling, seed, forcing coherence time and flow speed tried, so it does not produce the observed persistence.

**Why SFT was not enough (Q3).** In the 1-D SFT, earlier cycles contribute about 1e-5 of the imbalance at every τ. Hence the move to the dynamo, as the brief asked.

**The model**
- The transport is the declining-phase solver (`sc24_final/ftd_core.py`, read-only).
- Its eight checks pass unchanged (`q4_ftd_verify.py`): mass conservation to 5e-16, second-order operator, exact superposition, and the others.
- Closed into a self-excited loop (`q4_dynamo.py`):
  - a non-local Babcock–Leighton source α cosθ sinθ B̄_tach/(1+(B̄/B0)²), with emergence limited to |lat| < 35°;
  - α0 = 40 cm s⁻¹, v0 = 15 m s⁻¹;
  - full sphere, 65 × 121, period 13–14 yr.
- **Stochastic forcing per hemisphere:** α_h = α0(1 + σ ξ_h).
  - ξ is renewed every τ_c (1 yr by default) and σ = 0.6, chosen so the rms cycle asymmetry matches the Sun (0.10) at κ = 1.
  - For τ_c = 3 and 5 yr, σ is scaled as 1/√τ_c (0.35, 0.27). That under-corrects: rms A comes out at 0.14–0.15.
- **Hemispheric coupling, two knobs:**
  - κ multiplies the latitudinal diffusivity within about 10° of the equator;
  - ρ is the correlation between the hemispheres' forcing.
- **Slower flow:** v0 = 10 m s⁻¹ gives a 21-yr period under noise (15.8 yr without noise). All four indicators agree on that period: both polar reversals, the dipole proxy and the toroidal field.

**Checks of the closed model** (`q4_dynamo_checks.py`)
- With no noise the hemispheres stay identical to 6e-16.
- With the forcing swapped between hemispheres the run is the exact mirror (1e-15).
- At 97 × 181:
  - the period is unchanged (12.98 yr);
  - lag-1 r stays non-positive (−0.13, against −0.27 at 65 × 121 over the same 400 yr);
  - **but the asymmetry amplitude differs** (rms 0.17 vs 0.11). Dynamo amplitudes are therefore resolution-sensitive; only the qualitative persistence result is used.
- At 49 × 91 the period is 15 yr, so that grid is used only for the κ = 10 point, against its own κ = 1 partner.
- **Cycle detection fixed:** when the hemispheres drift out of phase, the dipole proxy P_N − P_S can cross zero without a reversal. One seed gave 113 crossings for 69 reversals. Crossings closer than 0.6 of the run's own polar-reversal period are now merged, and every run was re-analysed. The earlier "11-yr period on the coarse grid" was this artefact; it is 15 yr.

**Table** (generated by `q4_dynamo_ensemble.py` into `out/dynamo/ensemble_table.md`)
- Multi-seed settings: the median over seeds, the seed range in brackets, and r pooled over all adjacent cycle pairs of all seeds.
- Each run: 1100 yr (1250 yr for 10 m/s), minus a 100-yr (150-yr) spin-up.

| Setting | seeds (cycles) | Sun-like cycles | rms A | rms A_pf | rms N–S lag (yr) | lag-1 r of A, pooled (seed range) | mean longest run / 13 (range) | r(\|A_pf(n)\|, \|A(n+1)\|), pooled (range) | r(A_pf(n), A(n+1)), pooled (range) | period (yr) |
|---|---|---|---|---|---|---|---|---|---|---|
| **Sun** | — | — | 0.10 | 0.10 | 0.53 | **+0.57** (n = 12) | **7** | **+0.62** (n = 9) | **+0.68** (n = 9) | 11 |
| κ = 1 | 5 (347) | 21 % | 0.10 | 0.44 | 0.73 | +0.02 (−0.13 to +0.16) (n = 342) | 4.1 (3.5–4.9) | +0.02 (−0.18 to +0.28) | −0.07 (−0.16 to −0.03) | 13.9 |
| ρ = 0.8 | 5 (344) | 81 % | 0.06 | 0.27 | 0.38 | −0.05 (−0.18 to +0.00) (n = 339) | 4.2 (3.5–4.4) | +0.08 (−0.16 to +0.33) | +0.08 (−0.14 to +0.25) | 13.9 |
| τ_c = 3 yr | 3 (207) | 41 % | 0.14 | 0.33 | 0.96 | −0.01 (−0.13 to +0.08) (n = 204) | 4.5 (3.7–4.9) | +0.07 (−0.09 to +0.16) | −0.02 (−0.07 to +0.02) | 14.3 |
| τ_c = 5 yr | 3 (206) | 4 % | 0.15 | 0.42 | 0.87 | +0.13 (+0.08 to +0.16) (n = 203) | 3.9 (3.8–4.1) | +0.24 (+0.15 to +0.39) | −0.16 (−0.16 to −0.06) | 14.7 |
| v0 = 10 m/s | 3 (154) | 43 % | 0.09 | 0.41 | 0.89 | −0.01 (−0.05 to +0.04) (n = 151) | 4.3 (2.4–4.4) | +0.21 (+0.01 to +0.42) | −0.07 (−0.24 to +0.00) | 20.9 |
| κ = 0.3 | 1 (68) | 34 % | 0.15 | 0.39 | 0.85 | −0.04 (n = 67) | 3.8 | +0.15 | −0.18 | 13.0 |
| κ = 3 | 1 (69) | 29 % | 0.11 | 0.45 | 0.91 | −0.26 (n = 68) | 3.4 | +0.09 | −0.23 | 13.0 |
| ρ = 0.5 | 1 (69) | 54 % | 0.08 | 0.41 | 0.68 | −0.14 (n = 68) | 3.9 | −0.03 | −0.19 | 12.9 |
| κ = 1, 49 × 91 | 1 (69) | 30 % | 0.13 | 0.43 | 0.78 | +0.05 (n = 68) | 4.3 | −0.12 | +0.02 | 15.0 |
| κ = 10, 49 × 91 | 1 (65) | 55 % | 0.11 | 0.37 | 0.87 | −0.13 (n = 64) | 3.8 | +0.00 | −0.23 | 15.0 |

"Sun-like" means both poles reverse inside the cycle within 2 yr of each other and the polarity alternates.

**What holds**
- **No setting persists.**
  - Pooled lag-1 r is −0.05 to +0.13 across the multi-seed settings (151–342 pairs), and −0.26 to +0.05 in the single runs. The Sun gives +0.57.
  - The longest run averages 3.4–4.5 per 13 cycles: independent cycles give 4.0, the Sun 7.
- **The longer forcing coherence (5 yr) is the only setting with a consistent positive memory**: +0.13, positive in all three seeds.
  - It also gives the largest magnitude link (+0.24).
  - But it is the least Sun-like (4 % Sun-like cycles; polar asymmetry 0.42).
  - The memory comes from the forcing itself: 5 yr is about a third of a cycle.
- **The slower flow does not help.** At 10 m/s, lag-1 r is −0.01 (n = 151) and the period becomes 21 yr.
- **The polar-to-next-cycle link is weak or absent.**
  - Magnitude: pooled +0.02 to +0.24.
  - Signed r(A_pf(n), A(n+1)): +0.08 at best and negative in 4 of 5 multi-seed settings. The Sun gives +0.68.
- **The earlier single-run magnitude link was a fluke.** ρ = 0.8 gave +0.33, but the seed range is −0.16 to +0.33.
- **Seed scatter is large.** Over about 69 cycles, a single run's lag-1 r ranges ±0.15 between seeds. A 13-cycle record such as the Sun's is far noisier than that.

**What does not reproduce**
- The model's polar asymmetry is 2.7–4.5× too large: its hemispheres drift apart too easily. Even at ρ = 0.8, which is 81 % Sun-like, it is 0.27 against 0.10.
- The even/odd switch of Q2 is absent: the switch fraction is about 0.5 for both polarities.
- The decline-window sign rule is absent: 38/69 cycles in the κ = 1 run, and r(X_dec, A_pf) is −0.08 to +0.27 across settings.

**Coupling knobs**
- Equatorial diffusion helps synchronisation only when strong: κ = 10 roughly doubles the Sun-like share on the coarse grid (30 → 55 %).
- Correlated forcing does more: ρ = 0.8 gives 81 %.
- Neither creates persistence.

**Reading.** In this dynamo, hemispheric asymmetry is stochastic cycle by cycle, robustly across seeds, coupling, coherence time and flow speed. The Sun's persistence is therefore either a property this model lacks, or it is weaker than the 13-cycle record suggests. Q4b shows the 7-cycle run depends on single regions in SC17 and SC24.

**Not done**
- forcing memory beyond 5 yr;
- a source that couples the hemispheres directly (e.g. cross-equatorial regions);
- production runs at 97 × 181.

### Q4 in one line

Given where and when flux emerges, the polar asymmetry is largely set (77 % of model variance). That pattern is not fixed by hemispheric totals or by a handful of giants. Neither the SFT nor the coupled 2-D dynamo carries a hemisphere's dominance into the next cycle (dynamo: 5 seeds, coupling, coherence time 1–5 yr, flow 10–15 m/s). The data's persistence beyond lag 1 rests on a few regions in near-balanced cycles. The evidence leans to **stochastic**.

---

## What does not hold (plain list)

- The SFT reproduces the observed A_pf only in odd cycles (0.97; even 0.01). Every model statement about even cycles is weak.
- The physics window does not reproduce the 0.93 / 0.77 correlations. Data and model disagree on *when* the polar imbalance is set, unless the even-cycle switch reading above is right.
- H1 drops to about 0.74 when A_pf is redrawn within its published errors.
- The signs of the SC21 and SC14 flips are not established by region resampling (64 %, 69 %).
- The large-region carrier of the SC21 flip is not seen in DPD.
- There is no n+2 seeding, in the data or in the AR(2) fit.
- The 1-D SFT has no inter-cycle memory of the imbalance at any τ.

## Files

Code and outputs are in `analysis_v3/`. Run in this order:

| Script | What it does |
|---|---|
| `v0_checks.py` | engine and mirror checks |
| `q0_reconcile.py` | C1 vs +0.77 |
| `q1_sources.py` | SFT switches |
| `q2_window.py` | Green's function, scans, signs |
| `q2b_phase_sign.py` | the explicit sign model and the timing-lag links (needs `explore.json`, so run it after `explore.py`) |
| `q3_persistence.py` | runs, flips, DPD/Temmer, τ sweep |
| `ens_verify.py` | ensembles |
| `ens_extra.py` | H1 split, region bootstrap |
| `ens_transport_series.py` | per-cycle A_pf of the 27 transport configs (band in `ENS_headlines`) |
| `explore.py` | A–E |
| `figs_main.py` | P1–P7, Q2 table |
| `figs_q.py` | Q0–Q3, ensembles, P6, exploratory |
| `q4_stochastic.py` | Q4a randomised emergence, Q4b with/without the largest regions |
| `q4_ftd_verify.py` | the 2-D solver's own 8 checks (read-only import from the declining-phase project) |
| `q4_dynamo.py` | self-excited 2-D dynamo with per-hemisphere forcing and coupling (κ, ρ); `run NAME --kappa ... --rho ...` |
| `q4_dynamo_analyse.py` | cycle statistics of a run, as for the Sun (writes `out/dynamo/*_stats.json`) |
| `q4_dynamo_checks.py` | symmetry, mirror and resolution checks of the closed dynamo |
| `q4_dynamo_ensemble.py` | the Q4c table over seeds (`out/dynamo/ensemble_table.json`, `.md`) |
| `q4_figs.py` | Q4a–c figures |

The production dynamo runs were launched with:

```
q4_dynamo.py run prod_<name> --years 1100 --alpha0 40 --sigma 0.6 --tau_c 1 --wlat 35 [--kappa K] [--rho R] [--seed S] [--nr 49 --nth 91]
```

- Numbers are in `out/*.json`.
- Figures are in `figs/`, as PDF with TrueType fonts and 300-dpi PNG.
- External inputs (read only):
  - DPD: `Downloads/scaling polar field/.../Tilt Datasets/dpdtiltdata.txt`. Cite Baranyi 2015 and Győri et al. 2011, per the file's data policy.
  - Temmer/SILSO hemispheric SN: `revision/external/silso_hem/Catalogue_B.txt`.
  - Tilt workbook: `Datasets Used/tilt_complete_analysis.xlsx`.
