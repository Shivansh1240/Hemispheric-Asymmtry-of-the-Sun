"""explore.py -- exploratory track A-E (does not replace Q1-Q3).

A  counterfactual: remove the largest low-latitude (|lat| <= 20) region of each cycle; then put it back at 35 deg
   (same hemisphere, time, area; Joy tilt for the new latitude).  Change in modelled D, A_pf and in X1.
B  X_G = sum over regions of K(|lat|, lag) * Phi * sin(alpha) * (+1 N, -1 S), K the SFT imbalance kernel of a
   small tilted bipole (Green's function), read at the minimum ending the cycle.  Compared with X1.
C  sparse attribution: the fewest regions whose (exactly additive) contributions reach each cycle's |D|.
D  memory model A(n+1) = a A(n) + b A(n-1) + noise; simulated run lengths and lag-2 correlation.
E  hemispheric timing: area-weighted centroid time of each hemisphere per cycle, N - S lag vs dominance.
"""
import itertools
import numpy as np
from core import *

out = {}
M = Model(); base = M.run(CAT, keep=True)
apf_obs = np.array([APF[n] for n in C1])
yr = CAT.year.values; la = CAT.latc.values; ar = CAT.area.values; cy = CAT.cycle.values; gid = CAT.gid.values

def dD_of(rr, n):
    sel, c = rr["reg"][n]; N, S = rr["full"][n]
    return sel, np.sign(N) * c[:, 0] - np.sign(S) * c[:, 1], abs(N) - abs(S)

tic("A: counterfactuals")
cf = {}
for n in CYC:
    m = np.flatnonzero((cy == n) & (np.abs(la) <= 20))
    i = m[np.argmax(ar[m])]
    keep = np.ones(len(CAT), bool); keep[i] = False
    rem = M.run(CAT[keep].reset_index(drop=True))
    mv = CAT.copy(); Tn = mv.alpha.values[i] / np.sqrt(abs(la[i])); newlat = 35.0 * np.sign(la[i])
    mv.loc[i, "latc"] = newlat; mv.loc[i, "alpha"] = Tn * np.sqrt(35.0)
    mov = M.run(mv)
    # data side: X1 without that region / with its rows moved to 35 deg
    g = gid[i]; dmask = (D.gid != g).values
    ann_rm = A.annual(D[dmask])
    Dm = D.copy(); Dm.loc[~dmask, "abslat"] = 35.0; ann_mv = A.annual(Dm)
    D0, Dr, Dmv = d_of(base["full"][n]), d_of(rem["full"][n]), d_of(mov["full"][n])
    cf[n] = dict(region_area=float(ar[i]), region_lat=float(la[i]), region_year=float(yr[i]),
                 D=D0, D_removed=Dr, D_moved35=Dmv, dD_removed_frac=(Dr - D0) / abs(D0), dD_moved_frac=(Dmv - D0) / abs(D0),
                 Apf=apf_of(base["full"][n]), Apf_removed=apf_of(rem["full"][n]), Apf_moved35=apf_of(mov["full"][n]),
                 sign_flip_removed=bool(np.sign(Dr) != np.sign(D0)), sign_flip_moved=bool(np.sign(Dmv) != np.sign(D0)),
                 X1_01=x1(n, 0, 1), X1_01_removed=x1(n, 0, 1, ann=ann_rm), X1_01_moved=x1(n, 0, 1, ann=ann_mv),
                 X1_m14=x1(n, -1, 4), X1_m14_removed=x1(n, -1, 4, ann=ann_rm),
                 A=AN[n], A_removed=cycA(n, D.gid != g))
    v = cf[n]
    print(f"SC{n}: {v['region_area']:5.0f} MSH at {v['region_lat']:+5.1f}  D {v['D']:+.2e} -> removed {v['dD_removed_frac']:+.2f}, moved {v['dD_moved_frac']:+.2f} "
          f"| flip {v['sign_flip_removed']}/{v['sign_flip_moved']} | X1(0..1) {v['X1_01']:+.0f} -> {v['X1_01_removed']:+.0f} / moved {v['X1_01_moved']:+.0f}")
out["A_counterfactual"] = cf
c1 = [n for n in C1]
out["A_summary"] = dict(median_dD_removed=float(np.median([cf[n]["dD_removed_frac"] for n in CYC])),
                        median_dD_moved=float(np.median([cf[n]["dD_moved_frac"] for n in CYC])),
                        n_sign_flips_removed=int(sum(cf[n]["sign_flip_removed"] for n in CYC)),
                        n_sign_flips_moved=int(sum(cf[n]["sign_flip_moved"] for n in CYC)), n=len(CYC),
                        r_X1hat_obs=r(x1v(C1, 0, 1, True), apf_obs),
                        r_X1hat_removed_obs=r([parity(n) * cf[n]["X1_01_removed"] for n in C1], apf_obs))
print(out["A_summary"])

tic("B: Green's-function weighted predictor X_G")
# kernel: unit-flux bipole, separation 6 deg, tilt 2 deg, N hemisphere, normal orientation (+ at the northern spot)
s = M.s; Kl = np.arange(0.5, 45.0, 1.0)
lagK = np.arange(M.nk)
Ksrc = []
for l0 in Kl:
    h = 0.5 * 6.0 * np.sin(np.radians(2.0)); lp, lm = l0 + h, l0 - h
    sw = s.sigma_src(np.array([lm]), np.array([lp]), np.array([6.0 * np.cos(np.radians(2.0))]))[0]
    Ksrc.append(s.bmr_gauss(lm, lp, 1.0, sw) / (6.0 * np.sin(np.radians(2.0))))
Ksrc = np.array(Ksrc)                                         # (nl0, nlat), per unit (Phi * d * sin alpha)
KN = np.einsum("aj,jk->ak", Ksrc, M.G[:, 0, :]); KS = np.einsum("aj,jk->ak", Ksrc, M.G[:, 1, :])
Kimb = KN + KS            # dD = s_N (c_N + c_S) for an N-hemisphere region (s_S = -s_N at a minimum)
np.savez(os.path.join(OUT, "explore_kernel.npz"), lat=Kl, KN=KN, KS=KS, Kimb=Kimb, dt_yr=s.dt / SEC)
def XG(n, cat=CAT, own_only=True, tilt=None):
    tl = cat.alpha.values if tilt is None else tilt
    ks = np.floor((cat.year.values - T0[12]) * SEC / s.dt).astype(int)
    tk = np.floor((np.array([X.READ[n] - 1.6, X.READ[n] + 1.6]) - T0[12]) * SEC / s.dt).astype(int)
    m = (cat.cycle.values == n) if own_only else np.ones(len(cat), bool)
    m &= ks < tk[0]
    il = np.clip(np.round(np.abs(cat.latc.values[m]) - 0.5).astype(int), 0, len(Kl) - 1)
    lags = np.arange(tk[0], tk[1] + 1)[None, :] - ks[m][:, None]
    k = Kimb[il[:, None], np.clip(lags, 0, M.nk - 1)].mean(1)
    w = k * cat.flux.values[m] * cat.ddeg.values[m] * np.sin(np.radians(tl[m])) * np.sign(cat.latc.values[m])
    return float(w.sum())
xg = np.array([XG(n) for n in CYC]); dmod = np.array([d_of(base["full"][n]) for n in CYC])
xg10 = xg[3:]
out["B_XG"] = dict(r_XG_vs_modelD=r(xg, dmod), n_model=len(CYC),
                   r_XG_vs_obsApf=r(xg10, apf_obs), r_parityXG_vs_obsApf=r([parity(n) * v for n, v in zip(C1, xg10)], apf_obs),
                   r_X1hat_vs_obsApf=r(x1v(C1, 0, 1, True), apf_obs), r_X1_unsigned_vs_obsApf=r(x1v(C1, 0, 1), apf_obs),
                   r_XG_vs_A_n1=r(xg[:-1], [AN[n + 1] for n in CYC[:-1]]), r_X1m14_vs_A_n1=r(x1v(C2, -1, 4), [AN[n + 1] for n in C2]),
                   r_absXG_vs_absA_n1=r(np.abs(xg[:-1]), [abs(AN[n + 1]) for n in CYC[:-1]]),
                   r_XG_odd=r(xg10[::2], apf_obs[::2]), r_XG_even=r(xg10[1::2], apf_obs[1::2]), XG=dict(zip(CYC, xg)))
print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in out["B_XG"].items() if k != "XG"})

tic("C: sparse attribution / concentration")
conc = {}; sparse = {}
for n in CYC:
    sel, dD, Dn = dD_of(base, n)
    a = np.sign(Dn) * dD; o = np.argsort(-a); cs = np.cumsum(a[o]) / abs(Dn)
    k = int(np.argmax(cs >= 1.0) + 1) if (cs >= 1.0).any() else None
    nonzero = int((np.abs(dD) > 1e-6 * np.abs(dD).max()).sum())
    sparse[n] = dict(k_reach_D=k, n_contributing=len(sel), frac_regions=(k / len(sel)) if k else None,
                     top_areas=ar[sel][o][:k if k else 5].tolist()[:10], top_lat=la[sel][o][:10].tolist(),
                     top_own=bool(np.all(cy[sel][o][: (k or 1)] == n)), share_top1=float(cs[0]), share_top5=float(cs[min(4, len(cs) - 1)]),
                     share_top10=float(cs[min(9, len(cs) - 1)]), max_cum=float(cs.max()))
    conc[n] = cs[:200]
    print(f"SC{n}: k = {k} of {len(sel)} regions reach |D|; top1 {cs[0]:.2f} top5 {cs[4]:.2f} top10 {cs[9]:.2f} max cum {cs.max():.2f}")
out["C_sparse"] = sparse
ks_ = [v["k_reach_D"] for v in sparse.values() if v["k_reach_D"]]
out["C_summary"] = dict(k_median=float(np.median(ks_)), k_min=int(min(ks_)), k_max=int(max(ks_)),
                        frac_median=float(np.median([v["frac_regions"] for v in sparse.values() if v["frac_regions"]])))
np.savez(os.path.join(OUT, "explore_concentration.npz"), **{f"c{n}": conc[n] for n in CYC})
print(out["C_summary"])

tic("D: memory model")
nn = list(range(13, 24)); y = np.array([AN[n + 1] for n in nn]); x_1 = np.array([AN[n] for n in nn]); x_2 = np.array([AN[n - 1] for n in nn])
Xd = np.column_stack([np.ones(len(nn)), x_1, x_2]); coef, *_ = np.linalg.lstsq(Xd, y, rcond=None)
fit = Xd @ coef; sd = float(np.std(y - fit, ddof=3))
X1d = np.column_stack([np.ones(len(nn)), x_1]); c1_, *_ = np.linalg.lstsq(X1d, y, rcond=None); sd1 = float(np.std(y - X1d @ c1_, ddof=2))
out["D_memory"] = dict(n=len(nn), a=float(coef[1]), b=float(coef[2]), c=float(coef[0]), r_fit=r(fit, y), resid_sd=sd,
                       ar1_a=float(c1_[1]), ar1_r=r(X1d @ c1_, y), ar1_resid_sd=sd1)
rng = np.random.default_rng(42)
def simulate(a, b, c, sd, ncyc=13, nsim=20000):
    L, lag2, flips_big = [], [], []
    for _ in range(nsim):
        z = [AN[12], AN[13]]
        for k in range(200 + ncyc):
            z.append(c + a * z[-1] + b * z[-2] + rng.normal(0, sd))
        z = np.array(z[-ncyc:]); sg = np.sign(z)
        L.append(max(len(list(g)) for _, g in itertools.groupby(sg)))
        lag2.append(np.corrcoef(z[:-2], z[2:])[0, 1])
    return dict(mean_longest_run=float(np.mean(L)), frac_longest_ge7=float(np.mean(np.array(L) >= 7)), mean_lag2_r=float(np.nanmean(lag2)),
                mean_lag1_implied=float(a / (1 - b)) if b != 1 else None)
out["D_sim_ar2"] = simulate(coef[1], coef[2], coef[0], sd)
out["D_sim_ar1"] = simulate(c1_[1], 0.0, c1_[0], sd1)
out["D_sim_white"] = simulate(0.0, 0.0, float(np.mean(y)), float(np.std(y)))
print(out["D_memory"]); print("AR2 sim", out["D_sim_ar2"]); print("AR1 sim", out["D_sim_ar1"]); print("white", out["D_sim_white"])

tic("E: hemispheric timing")
tim = {}
for n in range(12, 25):
    m = (D.dec >= SCB[n]) & (D.dec < SCB[n + 1])
    for h in "NS":
        mh = m & (D.hemi == h); tim.setdefault(n, {})[h] = float(np.average(D.dec[mh], weights=D.area[mh]))
    tim[n]["lag_N_minus_S"] = tim[n]["N"] - tim[n]["S"]
lagv = np.array([tim[n]["lag_N_minus_S"] for n in range(12, 25)])
out["E_timing"] = dict(per_cycle=tim, r_lag_vs_A=r(lagv, [AN[n] for n in range(12, 25)]), n_A=13,
                       r_lag_vs_Apf=r(lagv[3:], apf_obs), n_Apf=10, r_lag_vs_X1hat=r(lagv[3:], x1v(C1, 0, 1, True)),
                       r_lag_vs_A_next=r(lagv[:-1], [AN[n + 1] for n in range(12, 24)]), n_next=12,
                       r_abs_lag_vs_abs_A=r(np.abs(lagv), [abs(AN[n]) for n in range(12, 25)]))
print({k: v for k, v in out["E_timing"].items() if k != "per_cycle"})
print("lags", np.round(lagv, 2))
dump("explore.json", out)
tic("done")
