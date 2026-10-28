# dynamo/ — Surya flux-transport dynamo (extension)

The observational analysis says *which* emergence pattern sets the polar asymmetry. A dynamo model is needed to ask whether that relation, and the cycle-to-cycle persistence of the dominant hemisphere, survive once the field is carried through the convection zone. This folder holds the model this extension will use: a modified copy of **Surya**, the 2-D axisymmetric kinematic flux-transport dynamo of A. R. Choudhuri's group (IISc), plus a relaxed initial state.

> **Credit.** Surya was written by Choudhuri and collaborators (Choudhuri, Schüssler & Dikpati 1995; Nandy & Choudhuri 2002; Chatterjee, Nandy & Choudhuri 2004; guide: Choudhuri 2005). The group gives the code and its guide to anyone who asks, but it has no published open-source licence. `surya_2.f` is a modified version used in this project. Cite the papers above if you use it, and ask the Surya authors before redistributing it further.

## Files

| File | What it is |
|---|---|
| `surya_2.f` | the dynamo code (fixed-form Fortran 77) |
| `init.dat` | relaxed starting state on the 257 × 257 grid: θ, r, r sinθ·A, B_φ for every grid point (66 049 rows) |
| `Makefile` | builds `surya` with gfortran |
| `run.sh` | one run in a fresh directory, then the plot |
| `plot_run.py` | butterfly diagrams + a north–south check of a finished run |
| `figures/butterfly.png` | output of the reference run below |

The asymmetry experiments run on extended copies of this code; see the section on asymmetry below.

## Running

```bash
cd dynamo
./run.sh            # builds, runs into run/, writes run/butterfly.png
```

The reference run takes about 13 minutes on one core (gfortran 13, `-O2`). Always run into an empty directory. The code opens `rad.dat`, `butbot.dat`, `series_b`, `ber.dat` and `psi.dat` with `access='append'`, so a second run in the same place silently adds to the first. `run.sh` refuses a non-empty directory for this reason.

## What this copy is set to

Code units are 10¹⁰ cm and 10⁸ s, so η is in 10¹² cm² s⁻¹ and speeds in m s⁻¹. The grid runs from 0.55 R☉ to R☉ in radius and pole to pole in θ (`nmax` = 257).

| Parameter | Value | Meaning |
|---|---|---|
| `irelax` | 1 | start from `init.dat` (0 = start from a seed toroidal field) |
| `tmax`, `dt` | 10, 2×10⁻⁴ | 10 units = 31.7 yr, 50 000 steps |
| `v0` | −24 | meridional-flow amplitude |
| `et0`, `et1` | 3.0, 0.04 | poloidal diffusivity in the convection zone; toroidal diffusivity in the bulk (it rises to `et0` above 0.975 R☉) |
| `al0` | 30 | α amplitude (Babcock–Leighton α confined to r > 0.95 R☉) |
| rotation | Ω = 271.9 (core) and 289.5 − 39.4 cos²θ − 42.2 cos⁴θ (envelope), in rad per 10⁸ s | helioseismic fit, i.e. 432.8 nHz in the core and 460.7 nHz at the equator; tachocline at 0.7 R☉ |

**Magnetic buoyancy.** Every 40 steps (8×10⁵ s ≈ 9 days), wherever |B_φ| > 0.8 (8×10⁴ G; the field unit is 10⁵ G) between 0.712 and 0.773 R☉, half of the toroidal field is removed there and added, scaled by the radius ratio, at r = 0.988 R☉ (grid point n − 6). Each event is logged to `ber.dat` as (t, θ, B_φ).

## Outputs

| File | Content |
|---|---|
| `rad.dat` | t, θ, B_r near the surface (r = 0.996 R☉), every π/20 time units (~0.5 yr) |
| `butbot.dat` | t, θ, B_φ at r = 0.706 R☉, same times |
| `series_b` | t (yr), \|B_φ\| at r = 0.70 R☉, latitude 14.8° N, every 40 steps |
| `ber.dat` | buoyant eruptions: t, θ, B_φ |
| `final.dat` | final state, same format as `init.dat`: copy it to `init.dat` to continue a run |
| `psi.dat`, `v_theta.dat` | stream function and v_θ of the meridional flow (11 MB each) |

## Reference run

![Butterfly diagrams of the reference run](figures/butterfly.png)

Output of `plot_run.py` for the run above (31.4 yr, 63 snapshots):

```
B_phi(0.706 R)  max |A_NS| = 3.74e-09   dipolar-parity residual = 9.90e-08
B_r(surface)    max |A_NS| = 3.83e-07   dipolar-parity residual = 9.01e-06
|B_phi| at r = 0.70 R, lat 14.8 N: mean 0.730, max 1.142 (code units)
buoyant eruptions: 5122 (2561 north, 2561 south)
```

The polar field reverses about every 11 yr (a 22-yr magnetic cycle). The north–south asymmetry is zero to the precision of the output files, and the eruptions split exactly evenly.

## On its own this copy has no hemispheric asymmetry

`init.dat` is exactly symmetric about the equator: A is even and B_φ odd under θ → π − θ, to 10⁻⁹ (the print precision). The α effect, flows, diffusivity and the buoyancy rule are symmetric too, so the run stays a pure dipole (see the numbers above). That is the right starting point for a controlled experiment, but nothing about asymmetry can be learnt from `surya_2.f` alone.

The project's asymmetry runs therefore use two extended copies, built from this file and kept separate so that `surya_2.f` stays untouched:

- **`surya_bhowmik.f`** (report v7): SURYA reset to the standard of the public guide (Choudhuri 2005), with the poloidal source and fluctuations of Passos et al. (2014): α_BL 27 m s⁻¹ on 10³–10⁵ G plus a quenched mean-field α of 0.4 m s⁻¹, uniform ±1 noise renewed every 6 months (100 % and 200 %), independent in the two hemispheres, and the 3.6 m s⁻¹ pumping of Bhowmik (2019). This is the dynamo Bhowmik (2019) used.
- **`surya_bhowmik_ext.f`** (report v8): the same plus switches for the noise coherence time τ_c, the N–S noise correlation ρ, an equatorial diffusivity factor κ, the flow speed, and a constant N–S bias in α_BL, the flow or η. With every switch off it reproduces `surya_bhowmik.f` bit for bit.

Both read a one-line `run.in` and write `polar.dat` (t, polar fields, eruptions and base fields per hemisphere). They are driven by `extension/surya_runs.py` and `extension/surya_v8.py`; the results are in `extension/reports/`. **The two Fortran files are not in the repository yet.** Add them to `dynamo/` before anyone tries to rerun v7 or v8.

### How `surya_2.f` differs from the public guide standard

From the comparison in report v7:

| Item | public SURYA (guide) | `surya_2.f` |
|---|---|---|
| meridional flow | v0 29 m s⁻¹, R_p 0.61 R☉, β 1.5 / 1.8, Γ 3.47, r₀ = (R☉ − R_b)/4 | v0 24, R_p 0.635 R☉, β 1.5 / 1.3, Γ 3.1, r₀ = (R☉ − R_b)/3.5 |
| diffusivity | η_SCZ 2.6×10¹² cm² s⁻¹, widths 0.025 R☉, steps at 0.70 / 0.72 / 0.975 R☉ | η_SCZ 3.0×10¹², widths 0.030 R☉, step at 0.725 R☉ |
| rotation | same fit, width 0.025 R☉ | width 0.030 R☉ |
| buoyancy | checked every 8.8×10⁵ s | B_c = 8×10⁴ G checked every 8×10⁵ s, half deposited at 0.988 R☉ |

## Known quirks of the code

- Output files are appended, not overwritten (see above).
- `rad.dat` and `butbot.dat` are written only when t lands within one step of a multiple of π/20, so the output cadence is ~0.5 yr. That is too coarse for the ±1.6 yr polar window. Write every few hundred steps before using it for A_pf.
- `PAUSE` statements (in `tridag` and `plgndr`) are deleted features in Fortran 95. gfortran accepts them with `-std=legacy`.
- When reading `init.dat`, A is divided by r sinθ; at θ = π this is 0/sin(π − 10⁻⁷). It works only because A = 0 there in the file.
