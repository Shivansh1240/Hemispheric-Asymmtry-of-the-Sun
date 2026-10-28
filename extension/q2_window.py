"""q2_window.py -- Q2: which emergence epoch feeds each pole, the window that follows from it, the whole scan table,
and an explicit model of the effective (parity) sign.

1. Exact attribution (deterministic Joy SFT, cap-70 readout +-1.6 yr about the minimum ending cycle n): every
   region's contribution c_i = (c_N, c_S).  With s_N = sgn N(n), s_S = sgn S(n) the pole-building parts are s_N c_N
   and s_S c_S, and the imbalance D = |N| - |S| = s_N N - s_S S splits exactly into dD_i = s_N c_N,i - s_S c_S,i.
2. Physics window: for every window [t_R+s, t_R+e+1) the model's windowed imbalance D_w(n) = sum of dD_i over
   regions emerging in it; r(D_w, D) over cycles 12-24 (n = 13) says how well that window carries the model's own
   polar imbalance.  No observation enters this choice.
3. Data scans (the whole table): r over window start x end for both targets; r over size class x latitude band.
4. Sign: effective sign s_eff(n) = sgn(A_pf X1) against parity; sign rules for the next cycle; magnitude links.
"""
import numpy as np
from core import *

out = {}
tic("model with per-region contributions")
M = Model(); res = M.run(CAT, keep=True)
yr = CAT.year.values; la = CAT.latc.values; ar = CAT.area.values; cy = CAT.cycle.values

E_EDGES = np.arange(-12.0, 12.01, 0.5); L_EDGES = np.arange(-45.0, 45.01, 2.5)
HN = np.zeros((len(E_EDGES) - 1, len(L_EDGES) - 1)); HS = np.zeros_like(HN); HD = np.zeros_like(HN)
normN = normS = normD = 0.0
rows = []
for n in CYC:
    sel, c = res["reg"][n]; N, S = res["full"][n]; sN, sS = np.sign(N), np.sign(S)
    trel = yr[sel] - TR[n]; dD = sN * c[:, 0] - sS * c[:, 1]; Dn = abs(N) - abs(S); sD = np.sign(Dn)
    HN += np.histogram2d(trel, la[sel], [E_EDGES, L_EDGES], weights=sN * c[:, 0])[0]; normN += abs(N)
    HS += np.histogram2d(trel, la[sel], [E_EDGES, L_EDGES], weights=sS * c[:, 1])[0]; normS += abs(S)
    HD += np.histogram2d(trel, la[sel], [E_EDGES, L_EDGES], weights=sD * dD)[0]; normD += abs(Dn)
    rows.append(dict(n=n, trel=trel, dD=dD, D=Dn, own=(cy[sel] == n), lat=la[sel], area=ar[sel], cN=c[:, 0] * sN, cS=c[:, 1] * sS))
np.savez(os.path.join(OUT, "q2_green_maps.npz"), HN=HN / normN, HS=HS / normS, HD=HD / normD, E=E_EDGES, L=L_EDGES)

# epoch profile of the aligned imbalance, summed over cycles
prof = HD.sum(1) / normD; ec = 0.5 * (E_EDGES[1:] + E_EDGES[:-1])
out["epoch_profile_imbalance"] = dict(zip(np.round(ec, 2), prof))
cum = np.cumsum(prof) / prof.sum()
out["imbalance_epoch_quantiles"] = {q: float(np.interp(q, cum, E_EDGES[1:])) for q in (0.16, 0.5, 0.84)}
profN = HN.sum(1) / normN; profS = HS.sum(1) / normS
for lab, pr in (("N", profN), ("S", profS)):
    cm = np.cumsum(pr) / pr.sum()
    out[f"pole_{lab}_epoch_quantiles"] = {q: float(np.interp(q, cm, E_EDGES[1:])) for q in (0.16, 0.5, 0.84)}
print("imbalance epoch quantiles", out["imbalance_epoch_quantiles"], "N", out["pole_N_epoch_quantiles"], "S", out["pole_S_epoch_quantiles"])

# per-hemisphere: which hemisphere's regions feed which pole (fractions)
fr = {}
for lab, col in (("N", "cN"), ("S", "cS")):
    tot = sum(rw[col].sum() for rw in rows)
    fr[f"pole_{lab}_from_N_regions"] = float(sum(rw[col][rw["lat"] > 0].sum() for rw in rows) / tot)
    fr[f"pole_{lab}_from_own_cycle"] = float(sum(rw[col][rw["own"]].sum() for rw in rows) / tot)
out["pole_feeding"] = fr; print(fr)

tic("physics window: r(D_window, D) over the model's own cycles")
WIN = [(s, e) for s in range(-5, 8) for e in range(s, 9)]
Dtot = np.array([rw["D"] for rw in rows]); phys = {}
for s, e in WIN:
    Dw = np.array([rw["dD"][(rw["trel"] >= s) & (rw["trel"] < e + 1)].sum() for rw in rows])
    phys[(s, e)] = r(Dw, Dtot) if np.std(Dw) > 0 else np.nan
out["physics_window_r"] = {f"{s},{e}": v for (s, e), v in phys.items()}
cands = sorted([(e - s, -phys[(s, e)], s, e) for (s, e) in WIN if phys[(s, e)] >= 0.9])
s_ph, e_ph = cands[0][2], cands[0][3]
best = max(WIN, key=lambda w: phys[w])
out["physics_window"] = dict(rule="shortest window with r(D_window, D) >= 0.9 in the model (n = 13 cycles)",
                             s=s_ph, e=e_ph, r_model=phys[(s_ph, e_ph)], best_any=dict(s=best[0], e=best[1], r=phys[best]))
print("physics window", out["physics_window"])

tic("data window scan, both targets")
apf = np.array([APF[n] for n in C1]); an1 = np.array([AN[n + 1] for n in C2])
scan = {}
for s, e in WIN:
    scan[(s, e)] = (r(x1v(C1, s, e, True), apf), r(x1v(C2, s, e, False), an1))
out["scan_window"] = {f"{s},{e}": dict(r_apf_n10=a, r_next_n12=b) for (s, e), (a, b) in scan.items()}
out["data_at_physics_window"] = dict(r_apf=scan[(s_ph, e_ph)][0], n_apf=10, r_next=scan[(s_ph, e_ph)][1], n_next=12)
out["data_at_paper_windows"] = dict(r_apf_0_1=scan[(0, 1)][0], r_next_m1_4=scan[(-1, 4)][1])
out["scan_rank_of_physics_window"] = dict(
    apf=int(sum(v[0] > scan[(s_ph, e_ph)][0] for v in scan.values())) + 1,
    next=int(sum(v[1] > scan[(s_ph, e_ph)][1] for v in scan.values())) + 1, of=len(scan))
print("data at physics window", out["data_at_physics_window"], "rank", out["scan_rank_of_physics_window"])

tic("size x latitude table")
SIZES = [("all", 0, np.inf), ("<100", 0, 100), ("100-499", 100, 500), ("500-999", 500, 1000), (">=1000", 1000, np.inf), (">=500", 500, np.inf)]
LATS = [("all", 0, 90), ("0-10", 0, 10), ("10-20", 10, 20), ("20-30", 20, 30), (">=30", 30, 90), ("0-20", 0, 20)]
tab = []
for sl, lo, hi in SIZES:
    for ll, a0, a1 in LATS:
        msk = band_mask(D, lo, hi, a0, a1); nreg = int(D.gid[msk].nunique())
        an = ann_of(msk)
        row = dict(size=sl, lat=ll, n_regions=nreg)
        for lab, (s, e) in (("phys", (s_ph, e_ph)), ("w01", (0, 1)), ("wm14", (-1, 4))):
            row[f"r_apf_{lab}"] = r(x1v(C1, s, e, True, an), apf)
            row[f"r_next_{lab}"] = r(x1v(C2, s, e, False, an), an1)
        tab.append(row)
out["scan_size_lat"] = tab
for t in tab:
    print(f"{t['size']:8s} {t['lat']:6s} n={t['n_regions']:6d}  apf phys {t['r_apf_phys']:+.2f} 0-1 {t['r_apf_w01']:+.2f} | next phys {t['r_next_phys']:+.2f} -1..4 {t['r_next_wm14']:+.2f}")

tic("sign model")
x01 = x1v(C1, 0, 1); xph = x1v(C1, s_ph, e_ph)
mod = apf_vec(res, C1)
seff = np.sign(apf * x01); seff_ph = np.sign(apf * xph); seff_mod = np.sign(mod * x01)
par = np.array([parity(n) for n in C1])
out["sign_polar"] = dict(cycles=C1, s_eff_obs=seff, s_eff_obs_physwin=seff_ph, s_eff_model=seff_mod, parity=par,
                         obs_matches_parity=int((seff == par).sum()), model_matches_parity=int((seff_mod == par).sum()),
                         model_matches_obs_s_eff=int((seff_mod == seff).sum()))
# A_pf = b * s(n) * X1 with s(n) = parity vs s(n) = +1 (no sign rule)
for lab, s_ in (("s=parity", par), ("s=+1", np.ones(10)), ("s=model_sign", seff_mod)):
    xx = s_ * x01; b = np.polyfit(xx, apf, 1)
    out[f"polar_sign_model_{lab}"] = dict(r=r(xx, apf), n=10, slope=float(b[0]), signs=signs(xx, apf))
print({k: v for k, v in out.items() if k.startswith("polar_sign_model")})
print("s_eff obs vs parity", out["sign_polar"]["obs_matches_parity"], "/10; model vs parity", out["sign_polar"]["model_matches_parity"], "/10")

# next cycle: sign rules and magnitude links
nn = C2; A1 = np.array([AN[n + 1] for n in nn])
rules = {"persistence sgn A(n)": [AN[n] for n in nn], "sgn X1(-1..+4)": x1v(nn, -1, 4),
         "sgn X1(physics window)": x1v(nn, s_ph, e_ph), "n+2 rule sgn A(n-1)": [AN[n - 1] if n - 1 >= 12 else np.nan for n in nn],
         "alternation -sgn A(n)": [-AN[n] for n in nn]}
c9 = list(range(15, 24)); A9 = np.array([AN[n + 1] for n in c9])
rules9 = {"sgn A_pf(n)": [APF[n] for n in c9], "sgn parity*A_pf(n)": [parity(n) * APF[n] for n in c9],
          "sgn S(n)": [S_SIGNED[n] for n in c9], "persistence sgn A(n)": [AN[n] for n in c9]}
sr = {}
for k, v in rules.items():
    v = np.array(v, float); m = np.isfinite(v); sr[k] = dict(hits=int((np.sign(v[m]) == np.sign(A1[m])).sum()), n=int(m.sum()))
for k, v in rules9.items():
    sr[k + " [15-23]"] = dict(hits=int((np.sign(v) == np.sign(A9)).sum()), n=9)
out["next_sign_rules"] = sr
mags = {"|A_pf(n)|": (c9, [abs(APF[n]) for n in c9]), "|X1(phys)|": (nn, np.abs(x1v(nn, s_ph, e_ph))),
        "|X1(-1..+4)|": (nn, np.abs(x1v(nn, -1, 4))), "|X1(0..+1)|": (nn, np.abs(x1v(nn, 0, 1))), "|A(n)|": (nn, [abs(AN[n]) for n in nn]),
        "|model A_pf(n)|": (c9, np.abs(apf_vec(res, c9)))}
out["next_magnitude_links"] = {k: dict(r=r(v, [abs(AN[n + 1]) for n in cc]), n=len(cc)) for k, (cc, v) in mags.items()}
# explicit two-part model: magnitude from |A_pf(n)|, sign from persistence  ->  A(n+1)_hat = sgn A(n) * |A_pf(n)|
twop = np.array([np.sign(AN[n]) * abs(APF[n]) for n in c9])
out["two_part_model_sgnA(n)*|Apf(n)|"] = dict(r=r(twop, A9), n=9)
out["two_part_model_sgnX1*|Apf(n)|"] = dict(r=r(np.sign(x1v(c9, -1, 4)) * np.abs([APF[n] for n in c9]), A9), n=9)
out["unsigned_X1_-1..4_same_9"] = dict(r=r(x1v(c9, -1, 4), A9), n=9)
print(sr); print(out["next_magnitude_links"])
print("two-part", out["two_part_model_sgnA(n)*|Apf(n)|"], out["two_part_model_sgnX1*|Apf(n)|"], "X1 same 9", out["unsigned_X1_-1..4_same_9"])
dump("q2_window.json", out)
tic("done")
