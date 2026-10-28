"""ens_verify.py -- ensemble verification of the three headline results.

H1  X1hat (t_R+0 -> +1, parity-signed) vs A_pf, n = 10
H2  X1 (t_R-1 -> +4, unsigned) vs next-cycle area asymmetry A(n+1), n = 12
H3  the Cycle 21 flip: A(20) > 0, A(21) < 0 in area; the modelled polar imbalance D(21) < 0

Emergence ensemble (100 members)
  * every daily area row x lognormal factor, sigma_ln = 0.20 (per-observation measurement error, assumed);
    X1, A(n) and A(n+1) are all recomputed from the perturbed rows;
  * A_pf perturbed by its published 1-sigma (independent-error reduction);
  * model: tilt scatter sigma(A) per region (as the project's ensemble) and the same area factors per region;
  * longitudes: neither X1 nor the 1-D (longitude-averaged) SFT uses longitude, so randomising it changes nothing.
Transport ensemble: eta in {150, 250, 400} km^2/s x v0 in {10, 15, 20} m/s x tau in {3, 5, 10} yr (27 configs,
tau by exact kernel rescaling), deterministic Joy tilts and 10 tilt-scatter members each.
"""
import numpy as np
from core import *

out = {}
rng = np.random.default_rng(20261024)
yrs = D.yr.values.astype(int); y0 = yrs.min(); ny = yrs.max() - y0 + 1; iy = yrs - y0
north = (D.hemi == "N").values; area0 = D.area.values; abslat = D.abslat.values; dec = D.dec.values
cnt = np.bincount(iy, minlength=ny); ml = np.bincount(iy, weights=abslat, minlength=ny) / np.maximum(cnt, 1)
yr_axis = np.arange(y0, y0 + ny)
cyc_idx = {n: (dec >= SCB[n]) & (dec < SCB.get(n + 1, 9999)) for n in range(12, 26)}

def fast(area):
    AN_ = np.bincount(iy, weights=np.where(north, area, 0.0), minlength=ny)
    AS_ = np.bincount(iy, weights=np.where(north, 0.0, area), minlength=ny)
    x1y = np.where(cnt > 0, (AN_ - AS_) / np.where(ml > 0, ml, 1), 0.0)
    def X(n, s, e, signed):
        m = (yr_axis + 0.5 >= TR[n] + s) & (yr_axis + 0.5 < TR[n] + e + 1)
        v = x1y[m].sum(); return v * parity(n) if signed else v
    A_ = {}
    for n, m in cyc_idx.items():
        N = area[m & north].sum(); S = area[m & ~north].sum(); A_[n] = (N - S) / (N + S)
    return X, A_

X0, A0 = fast(area0)
chk = dict(x1hat=float(max(abs(X0(n, 0, 1, True) - x1(n, 0, 1, True)) for n in C1)),
           x1=float(max(abs(X0(n, -1, 4, False) - x1(n, -1, 4)) for n in C2)),
           A=float(max(abs(A0[n] - AN[n]) for n in range(12, 25))))
out["fast_reimplementation_max_abs_diff"] = chk; print("fast path check", chk)

apf = np.array([APF[n] for n in C1]); sig = np.array([SIG[n] for n in C1])
H1, H2, H3, X1E, XHE, AE = [], [], [], [], [], []
tic("emergence ensemble (data side)")
for k in range(100):
    f = np.exp(rng.normal(0.0, 0.20, len(area0)))
    Xk, Ak = fast(area0 * f)
    xh = np.array([Xk(n, 0, 1, True) for n in C1]); x14 = np.array([Xk(n, -1, 4, False) for n in C2])
    an1 = np.array([Ak[n + 1] for n in C2]); ap = apf + rng.normal(0, 1, 10) * sig
    H1.append(r(xh, ap)); H2.append(r(x14, an1)); H3.append((Ak[20], Ak[21]))
    XHE.append(xh); X1E.append(x14); AE.append(an1)
H1, H2, H3 = np.array(H1), np.array(H2), np.array(H3)
out["H1_emergence"] = dict(r_obs=r(x1v(C1, 0, 1, True), apf), median=float(np.median(H1)), p16_84=np.percentile(H1, [16, 84]).tolist(),
                           p2_98=np.percentile(H1, [2.5, 97.5]).tolist(), n=10)
out["H2_emergence"] = dict(r_obs=r(x1v(C2, -1, 4), [AN[n + 1] for n in C2]), median=float(np.median(H2)),
                           p16_84=np.percentile(H2, [16, 84]).tolist(), p2_98=np.percentile(H2, [2.5, 97.5]).tolist(), n=12)
out["H3_area_emergence"] = dict(A20_median=float(np.median(H3[:, 0])), A21_median=float(np.median(H3[:, 1])),
                                A21_p16_84=np.percentile(H3[:, 1], [16, 84]).tolist(), frac_A21_negative=float((H3[:, 1] < 0).mean()),
                                frac_A20_positive=float((H3[:, 0] > 0).mean()))
print(out["H1_emergence"]); print(out["H2_emergence"]); print(out["H3_area_emergence"])
np.savez(os.path.join(OUT, "ens_data.npz"), H1=H1, H2=H2, H3=H3, XHE=np.array(XHE), X1E=np.array(X1E), AE=np.array(AE))

tic("emergence ensemble (model side): tilt scatter + area factors per region")
M = Model()
gL = (np.abs(CAT.latc.values) <= 20) & (CAT.area.values >= 500)
def giant_share(rr):
    sh = []
    for n in C1:
        sel, c = rr["reg"][n]; N, S = rr["full"][n]; dD = np.sign(N) * c[:, 0] - np.sign(S) * c[:, 1]
        own = CAT.cycle.values[sel] == n; Dn = abs(N) - abs(S)
        sh.append(float(dD[own & gL[sel]].sum() / Dn))
    return sh
base = M.run(CAT, keep=True)
out["giant_share_base"] = dict(per_cycle=dict(zip(C1, giant_share(base))), median=float(np.median(giant_share(base))))
print("low-lat large share of D, base:", round(out["giant_share_base"]["median"], 3))
ME = []
for k in range(100):
    fac = np.exp(rng.normal(0.0, 0.20 / np.sqrt(8), len(CAT)))     # region totals average ~ several daily rows
    t = tilt_draw(CAT, 7001 + k)
    rr = M.run(scaled(CAT, fac), t)
    ME.append([apf_of(rr["full"][n]) for n in CYC] + [d_of(rr["full"][n]) for n in CYC])
ME = np.array(ME); nC = len(CYC)
mA, mD = ME[:, :nC], ME[:, nC:]
i21 = CYC.index(21); i20 = CYC.index(20)
xh = x1v(C1, 0, 1, True)
out["model_emergence"] = dict(frac_D21_negative=float((mD[:, i21] < 0).mean()), frac_D20_positive=float((mD[:, i20] > 0).mean()),
                              r_model_obs_median=float(np.median([r(a[3:], apf) for a in mA])),
                              r_model_obs_16_84=np.percentile([r(a[3:], apf) for a in mA], [16, 84]).tolist(),
                              r_X1hat_model_median=float(np.median([r(xh, a[3:]) for a in mA])))
print(out["model_emergence"])
np.savez(os.path.join(OUT, "ens_model_emergence.npz"), A=mA, D=mD, cycles=CYC)

tic("transport ensemble (27 configs)")
TR_ = {}
obsA1 = [abs(AN[n + 1]) for n in range(15, 24)]
for eta in (150.0, 250.0, 400.0):
    for v0 in (10.0, 15.0, 20.0):
        Mb = Model(eta=eta, v0=v0, tau=5.0)
        for tau in (3.0, 5.0, 10.0):
            Mt = Mb if tau == 5.0 else Model(eta=eta, v0=v0, tau=tau, base=Mb)
            rr = Mt.run(CAT, keep=True)
            a = np.array([apf_of(rr["full"][n]) for n in CYC]); dd = np.array([d_of(rr["full"][n]) for n in CYC])
            # imbalance epoch median (physics window centre) for this transport
            num = []; ep = []
            for n in CYC:
                sel, c = rr["reg"][n]; N, S = rr["full"][n]; dD = np.sign(N) * c[:, 0] - np.sign(S) * c[:, 1]
                num.append(np.sign(abs(N) - abs(S)) * dD); ep.append(CAT.year.values[sel] - TR[n])
            num = np.concatenate(num); ep = np.concatenate(ep); o = np.argsort(ep); cw = np.cumsum(num[o]) / num.sum()
            ens = [apf_vec(Mt.run(CAT, tilt_draw(CAT, 7001 + k)), C1) for k in range(10)]
            TR_[f"eta{eta:.0f}_v{v0:.0f}_tau{tau:.0f}"] = dict(
                r_model_obs=r(a[3:], apf), r_model_obs_tiltens_median=float(np.median([r(x, apf) for x in ens])),
                r_X1hat_model=r(xh, a[3:]), r_abs_model_vs_abs_A_n1=r(np.abs(a[3:12]), obsA1),
                sign_D21=int(np.sign(dd[i21])), sign_D20=int(np.sign(dd[i20])),
                frac_D21_neg_tiltens=float(np.mean([e[C1.index(21)] < 0 for e in ens])),     # sgn A_pf = sgn D
                giant_share_median=float(np.median(giant_share(rr))), imbalance_epoch_median=float(np.interp(0.5, cw, ep[o])),
                rms_model_apf=float(np.sqrt(np.mean(a[3:] ** 2))))
            v = TR_[f"eta{eta:.0f}_v{v0:.0f}_tau{tau:.0f}"]
            print(f"eta {eta:3.0f} v0 {v0:2.0f} tau {tau:2.0f}: r(mod,obs) {v['r_model_obs']:+.2f} r(X1hat,mod) {v['r_X1hat_model']:+.2f} "
                  f"D21 {v['sign_D21']:+d} (ens neg {v['frac_D21_neg_tiltens']:.1f}) giant {v['giant_share_median']:.2f} epoch {v['imbalance_epoch_median']:.1f}")
out["transport"] = TR_
vals = lambda k: np.array([v[k] for v in TR_.values()])
out["transport_summary"] = {k: dict(min=float(vals(k).min()), median=float(np.median(vals(k))), max=float(vals(k).max()))
                            for k in ("r_model_obs", "r_model_obs_tiltens_median", "r_X1hat_model", "r_abs_model_vs_abs_A_n1",
                                      "giant_share_median", "imbalance_epoch_median", "frac_D21_neg_tiltens", "rms_model_apf")}
out["transport_summary"]["D21_negative_configs"] = int((vals("sign_D21") < 0).sum())
out["transport_summary"]["n_configs"] = len(TR_)
print(json.dumps(out["transport_summary"], indent=1))
dump("ens_verify.json", out)
tic("done")
