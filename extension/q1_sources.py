"""q1_sources.py -- Q1: what sets the polar-field asymmetry?

Bhowmik-parameter SFT (eta 250 km^2/s, v0 15 m/s, lambda0 75 deg; tau 5 yr is the project's addition), cap-70
readout +-1.6 yr about the minimum that ends each cycle, observed A_pf for SC15-24 (n = 10).

Reference: the SYMMETRISED catalogue (every region split into two half-flux copies at +|lat| and -|lat|) gives
identical hemispheres and A_pf = 0.  Hemispheric differences are then switched on one at a time:
  (a)  area only       real catalogue (amount, latitude and timing of emergence), Joy tilts, symmetric transport
       (a1) amount only    symmetric copies, each hemisphere scaled to its observed flux share in that cycle
       (a2) placement only real catalogue, the two hemispheres scaled to equal flux in each cycle
  (b)  tilt only       symmetric copies + independent tilt-scatter draws sigma(A) per copy (100 members)
  (c)  transport only  symmetric copies + eta_N != eta_S (+-10, 20 %) or v0_N != v0_S (+-10, 20 %)
  (d)  all             real catalogue + tilt scatter + each transport asymmetry (8 configs x 25 members)
Also: X1hat * T_lin (cycle Joy slope from the tilt workbook) and the algebraic hemispheric dipole predictor with
lambda_R from a different diffusivity in each hemisphere.
"""
import numpy as np
import pandas as pd
from core import *

obs = np.array([APF[n] for n in C1]); odd = np.array([n % 2 == 1 for n in C1])
out = {}

def met(a):
    a = np.asarray(a, float)
    return dict(rms=float(np.sqrt(np.mean(a ** 2))), r_obs=r(a, obs), signs=signs(a, obs),
                r_odd=r(a[odd], obs[odd]), r_even=r(a[~odd], obs[~odd]))

def flux_share(cat):
    f = {}
    for n in range(12, 26):
        m = cat.cycle.values == n
        fN = cat.flux.values[m & (cat.latc.values > 0)].sum(); fS = cat.flux.values[m & (cat.latc.values < 0)].sum()
        f[n] = (fN, fS)
    return f

tic("base model"); M = Model()
sym = symmetrised(CAT); FS = flux_share(CAT)
res0 = M.run(sym); out["reference_symmetric"] = met(apf_vec(res0)); print("reference", out["reference_symmetric"])

tic("(a) area only"); ra = M.run(CAT); a_area = apf_vec(ra); out["a_area"] = met(a_area)
# (a1) amount only
fac = np.empty(len(sym)); north = sym.latc.values > 0
for n in range(12, 26):
    m = sym.cycle.values == n; fN, fS = FS[n]; tot = fN + fS
    if tot > 0:
        fac[m & north] = 2 * fN / tot; fac[m & ~north] = 2 * fS / tot
fac[~np.isfinite(fac)] = 1.0
a_amount = apf_vec(M.run(scaled(sym, fac))); out["a1_amount_only"] = met(a_amount)
# (a2) placement only
fac2 = np.ones(len(CAT)); nN = CAT.latc.values > 0
for n in range(12, 26):
    m = CAT.cycle.values == n; fN, fS = FS[n]
    if fN > 0 and fS > 0:
        fac2[m & nN] = 0.5 * (fN + fS) / fN; fac2[m & ~nN] = 0.5 * (fN + fS) / fS
a_place = apf_vec(M.run(scaled(CAT, fac2))); out["a2_placement_only"] = met(a_place)
for k in ("a_area", "a1_amount_only", "a2_placement_only"):
    print(k, out[k])

tic("(b) tilt only, 100 members")
sig_sym = scatter_sigma(2 * sym.area.values)
B_tilt = []
for k in range(100):
    t = sym.alpha.values + np.random.default_rng(9100 + k).normal(0, 1, len(sym)) * sig_sym
    B_tilt.append(apf_vec(M.run(sym, t)))
B_tilt = np.array(B_tilt)
out["b_tilt_only"] = dict(rms_median=float(np.median(np.sqrt(np.mean(B_tilt ** 2, 1)))),
                          r_obs_median=float(np.median([r(b, obs) for b in B_tilt])),
                          r_obs_16_84=np.percentile([r(b, obs) for b in B_tilt], [16, 84]).tolist(),
                          per_cycle_sd=dict(zip(C1, B_tilt.std(0))))
print("b", out["b_tilt_only"]["rms_median"], out["b_tilt_only"]["r_obs_median"], out["b_tilt_only"]["r_obs_16_84"])

tic("(c) transport only and (d) all")
TCFG = {f"eta{d:+.1f}": dict(etaN=250 * (1 + d), etaS=250 * (1 - d)) for d in (-0.2, -0.1, 0.1, 0.2)}
TCFG.update({f"v0{d:+.1f}": dict(v0N=15 * (1 + d), v0S=15 * (1 - d)) for d in (-0.2, -0.1, 0.1, 0.2)})
C_tr, D_all = {}, {}
tilt_members = [tilt_draw(CAT, 7000 + k) for k in range(25)]
for name, kw in TCFG.items():
    Mt = Model(**kw)
    C_tr[name] = apf_vec(Mt.run(sym))
    D_all[name] = np.array([apf_vec(Mt.run(CAT, t)) for t in tilt_members])
    print(name, "transport only rms", np.sqrt(np.mean(C_tr[name] ** 2)).round(4), "all: median r_obs",
          np.median([r(x, obs) for x in D_all[name]]).round(3))
out["c_transport_only"] = {k: dict(met(v), apf=dict(zip(C1, v))) for k, v in C_tr.items()}
Dstack = np.concatenate(list(D_all.values()))
rD = [r(x, obs) for x in Dstack]
out["d_all"] = dict(n_members=len(Dstack), rms_median=float(np.median(np.sqrt(np.mean(Dstack ** 2, 1)))),
                    r_obs_median=float(np.median(rD)), r_obs_16_84=np.percentile(rD, [16, 84]).tolist(),
                    r_with_area_only_median=float(np.median([r(x, a_area) for x in Dstack])))
# Joy tilts + scatter, symmetric transport (the project's own ensemble configuration)
Esym = np.array([apf_vec(M.run(CAT, t)) for t in tilt_members])
rE = [r(x, obs) for x in Esym]
out["area_plus_tilt_scatter"] = dict(rms_median=float(np.median(np.sqrt(np.mean(Esym ** 2, 1)))),
                                     r_obs_median=float(np.median(rE)), r_obs_16_84=np.percentile(rE, [16, 84]).tolist())
# how much of the full model is the area term:  1 - var(d - a)/var(d), member by member
out["area_share_of_full_model_variance_median"] = float(np.median([1 - np.var(x - a_area) / np.var(x) for x in Dstack]))
print("d", out["d_all"]); print("area+scatter", out["area_plus_tilt_scatter"])
print("area share of full-model variance", out["area_share_of_full_model_variance_median"])

tic("X1hat * T_lin")
tw = pd.read_excel(os.path.join(HERE, "..", "Datasets Used", "tilt_complete_analysis.xlsx"), sheet_name="Phase_Summary")
tw = tw[(tw.Filter == "Ds>2.5") & (tw.Scope == "Full")]
Tlin = {n: float(tw[(tw.Cycle == n) & (tw.Dataset == ("KK" if n <= 21 else "DPD"))].T_lin.iloc[0]) for n in C1}
xh = x1v(C1, 0, 1, True)
out["Tlin"] = Tlin
out["X1hat"] = dict(r=r(xh, obs), n=10)
out["X1hat_times_Tlin"] = dict(r=r(xh * np.array([Tlin[n] for n in C1]), obs), n=10)
out["Tlin_vs_|Apf|"] = dict(r=r([Tlin[n] for n in C1], np.abs(obs)), n=10)
print("X1hat", out["X1hat"], "X1hat*Tlin", out["X1hat_times_Tlin"], "Tlin vs |Apf|", out["Tlin_vs_|Apf|"])

tic("algebraic hemispheric dipole with eta_N != eta_S")
R_SUN, TAU, U0, LAM0, SIG_SPOT = 6.96e10, 5.0, 1500.0, 75.0, 0.56
lamR = lambda eta: float(np.sqrt(np.radians(SIG_SPOT) ** 2 + (eta * 1e10 / R_SUN ** 2) / ((U0 * np.pi / np.radians(LAM0)) / R_SUN)))
h = 0.5 * CAT.ddeg.values * np.sin(np.radians(CAT.alpha.values))
base = CAT.parity.values * (3.0 / (4 * np.pi * R_SUN ** 2)) * 0.5 * CAT.flux.values * \
    (np.sin(np.radians(CAT.latc.values + h)) - np.sin(np.radians(CAT.latc.values - h)))
yr, cy, la = CAT.year.values, CAT.cycle.values, CAT.latc.values
def adp(etaN, etaS):
    qN = np.exp(-np.radians(np.abs(la)) ** 2 / (2 * lamR(etaN) ** 2)); qS = np.exp(-np.radians(np.abs(la)) ** 2 / (2 * lamR(etaS) ** 2))
    v = []
    for n in C1:
        m = (cy == n) & (yr >= T0[n]) & (yr <= T0[n + 1]); w = base * np.exp((yr - T0[n + 1]) / TAU)
        dN = (w * qN)[m & (la > 0)].sum(); dS = (w * qS)[m & (la < 0)].sum()
        v.append((abs(dN) - abs(dS)) / (abs(dN) + abs(dS)))
    return np.array(v)
alg = {}
for d in np.round(np.arange(-0.6, 0.61, 0.1), 1):
    a = adp(250 * (1 + d), 250 * (1 - d)); alg[float(d)] = dict(r_obs=r(a, obs), rms=float(np.sqrt(np.mean(a ** 2))), lamR_N_deg=np.degrees(lamR(250 * (1 + d))), lamR_S_deg=np.degrees(lamR(250 * (1 - d))))
out["algebraic_hemi_eta"] = alg
print({k: round(v["r_obs"], 3) for k, v in alg.items()})

out["series"] = dict(obs=dict(zip(C1, obs)), area=dict(zip(C1, a_area)), amount=dict(zip(C1, a_amount)),
                     placement=dict(zip(C1, a_place)), tilt_only_p16=dict(zip(C1, np.percentile(B_tilt, 16, 0))),
                     tilt_only_p84=dict(zip(C1, np.percentile(B_tilt, 84, 0))),
                     all_p16=dict(zip(C1, np.percentile(Dstack, 16, 0))), all_p84=dict(zip(C1, np.percentile(Dstack, 84, 0))),
                     all_median=dict(zip(C1, np.median(Dstack, 0))))
np.savez(os.path.join(OUT, "q1_ensembles.npz"), tilt_only=B_tilt, all=Dstack, area_scatter=Esym)
dump("q1_sources.json", out)
tic("done")
