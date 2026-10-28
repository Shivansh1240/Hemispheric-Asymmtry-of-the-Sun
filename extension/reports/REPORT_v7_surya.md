# v7: Bhowmik's dynamo run on the SURYA code

Earlier reports are unchanged. Code: `surya/surya_bhowmik.f` (a copy of `master_code/surya_2.f`, which is untouched), driver `surya_runs.py`, numbers `out/v7_surya/`, figure `plots/Q4_stochastic_vs_predetermined/F6_bhowmik_dynamo_surya.png`. Compiled with gfortran 15.2 (MSYS2 UCRT64).

## Set-up

| Item | Bhowmik / Passos / guide | supplied surya_2.f | used here |
|---|---|---|---|
| code | public SURYA (Choudhuri 2005 guide) | supplied surya_2.f, a modified SURYA | SURYA, reset to the guide's standard |
| meridional flow | v0 29 m/s, R_p 0.61 R, beta 1.5 / 1.8, Gamma 3.47, r0 = (R - R_b)/4 | v0 24, R_p 0.635 R, beta 1.5 / 1.3, Gamma 3.1, r0 = (R - R_b)/3.5 | guide; normalised to 29 m/s at 45 deg |
| diffusivity | eta_SCZ 2.6e12, eta_SCZ1 4e10, widths 0.025 R, steps at 0.70 / 0.72 / 0.975 R | eta_SCZ 3.0e12, widths 0.030 R, step at 0.725 R | guide |
| rotation | CNC fit, width 0.025 R | same fit, width 0.030 R | guide |
| poloidal source | Passos: alpha_BL 27 m/s on 1e3-1e5 G (Eq. 4) + alpha_MF 0.4 m/s quenched at 1e4 G (Eq. 5) | alpha 30 m/s, no window, no alpha_MF | Passos |
| fluctuations | Passos: uniform +-1, every 6 months, BL 100 %, MF 200 %, hemispheres independent | none | Passos (Bhowmik's own runs are data-driven) |
| pumping | Bhowmik Eq. 7, 3.6 m/s, poloidal field only | none | Bhowmik |
| buoyancy | B_c 8e4 G (Passos) every 8.8e5 s (guide), half erupts | B_c 8e4 G every 8e5 s, deposited at r = 0.989 R | Passos / guide; deposition as in the supplied code |
| grid, start | 256 x 256 | 256 x 256, init.dat | same; first 100 yr dropped |

## Checks

- Noise-free run: period 11.0 yr (10 minima after 30 yr), eruptions in 100 % of 0.05-yr bins, base field up to 3.51 x 1e5 G; north and south identical to 6e-06.
- Mirror check (noise swapped between hemispheres): polar field 4e-08, eruptions 0e+00.

## Persistence in Bhowmik's dynamo

| Model | cycles | pooled lag-1 r of A (seed range) | mean longest run / 13 | rms A | rms A_pf | Sun-like cycles |
|---|---|---|---|---|---|---|
| **Sun** | 13 | **+0.57** | **7** | 0.10 | 0.10 | — |
| our Q4 dynamo, kappa = 1 | 347 | +0.02 (−0.13 to +0.16) | 4.1 | 0.10 | 0.44 | 21 % |
| our Q4 dynamo, tau_c = 30 yr (T3) | 207 | +0.58 (+0.53 to +0.62) | 6.1 | 0.09 | 0.12 | 94 % |
| **Bhowmik's dynamo (SURYA)** | 270 | **+0.45** (+0.26 to +0.56) | **5.5** | 0.18 | **0.29** | 78 % |

- Period 10.7 yr; timing-lag lag-1 r +0.14 (Sun +0.34). A is the erupted toroidal flux asymmetry per cycle (Bhowmik's B_Dyn proxy), A_pf the polar-cap (>= 70 deg) asymmetry at each minimum; same per-cycle code as Q4.

| seed | cycles | period (yr) | lag-1 r | longest run | rms A | rms A_pf | Sun-like |
|---|---|---|---|---|---|---|---|
| 1 | 45 | 10.7 | +0.50 | 4.2 | 0.20 | 0.24 | 84 % |
| 2 | 45 | 10.7 | +0.56 | 6.4 | 0.20 | 0.40 | 56 % |
| 3 | 45 | 10.6 | +0.31 | 5.6 | 0.16 | 0.27 | 84 % |
| 4 | 45 | 10.5 | +0.37 | 6.6 | 0.15 | 0.27 | 53 % |
| 5 | 45 | 10.6 | +0.39 | 5.2 | 0.17 | 0.20 | 84 % |
| 6 | 45 | 10.8 | +0.26 | 5.5 | 0.20 | 0.32 | 71 % |
