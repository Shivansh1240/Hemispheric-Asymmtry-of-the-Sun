"""v0_checks.py -- engine checks before any result is quoted.

1. The cap-70-only Green's function reproduces the project's deterministic readout (ensemble member 0, analog A1).
2. Decay-time rescaling of the kernel equals a direct build at that tau.
3. MIRROR CHECK: swap N and S in the catalogue (lat -> -lat, parity column unchanged) and in the daily area rows.
   Expected exactly: model N' = -S, S' = -N (so A_pf' = -A_pf, D' = -D, parity-signed quantities flip, r unchanged);
   data X1' = -X1, A' = -A.  Also with hemispheric transport (eta_N != eta_S): mirror + swap eta_N <-> eta_S.
"""
import numpy as np
from core import *

out = {}
tic("parity column: unique values per (cycle, hemisphere)")
g = CAT.assign(h=np.sign(CAT.latc)).groupby(["cycle", "h"]).parity.agg(lambda v: sorted(set(v)))
out["parity_column"] = {f"{c}{'N' if h > 0 else 'S'}": v for (c, h), v in g.items() if 12 <= c <= 25}
print({k: v for k, v in list(out["parity_column"].items())[:8]})

tic("building default model (eta 250, v0 15, tau 5)")
M = Model()
res = M.run(CAT)
a = apf_vec(res, C1)
ens = A.ensemble_apf("A1", C1, member=0)
out["check1_max_abs_diff_vs_ensemble_member0"] = float(np.max(np.abs(a - ens)))
out["model_apf_det"] = dict(zip(C1, a))
print("max |A_pf - ensemble member 0| =", out["check1_max_abs_diff_vs_ensemble_member0"])

tic("tau rescaling vs direct build (tau = 8 yr)")
M8r = Model(tau=8.0, base=M); M8d = Model(tau=8.0)
r8r, r8d = M8r.run(CAT), M8d.run(CAT)
out["check2_tau8_rescaled_vs_direct_max_rel"] = float(max(np.max(np.abs(r8r["full"][n] - r8d["full"][n]) / np.abs(r8d["full"][n]).max()) for n in CYC))
print("tau rescale max rel diff", out["check2_tau8_rescaled_vs_direct_max_rel"])

tic("mirror check: model")
rm = M.run(mirror(CAT))
dN = max(abs(rm["full"][n][0] + res["full"][n][1]) for n in CYC)
dS = max(abs(rm["full"][n][1] + res["full"][n][0]) for n in CYC)
scale = max(np.abs(res["full"][n]).max() for n in CYC)
out["check3_model_mirror_max_rel"] = float(max(dN, dS) / scale)
am = apf_vec(rm, C1)
out["check3_apf_mirror_max_abs(a'+a)"] = float(np.max(np.abs(am + a)))
print("model mirror rel residual", out["check3_model_mirror_max_rel"], " A_pf' + A_pf max", out["check3_apf_mirror_max_abs(a'+a)"])

tic("mirror check: hemispheric transport")
Mh = Model(etaN=300.0, etaS=200.0); Mhm = Model(etaN=200.0, etaS=300.0)
rh, rhm = Mh.run(CAT), Mhm.run(mirror(CAT))
out["check3b_hemi_transport_mirror_max_rel"] = float(max(max(abs(rhm["full"][n][0] + rh["full"][n][1]), abs(rhm["full"][n][1] + rh["full"][n][0])) for n in CYC) / scale)
print("hemi transport mirror rel residual", out["check3b_hemi_transport_mirror_max_rel"])

tic("mirror check: data")
Dm = D.copy(); Dm["hemi"] = np.where(D.hemi == "N", "S", "N"); Dm["lat"] = -D.lat
annm = A.annual(Dm)
x = x1v(C2, -1, 4); xm = x1v(C2, -1, 4, ann=annm)
xs = x1v(C1, 0, 1, True); xsm = x1v(C1, 0, 1, True, ann=annm)
Am = np.array([A.cycle_asymmetry(n, None, Dm) for n in range(12, 25)]); Ao = np.array([AN[n] for n in range(12, 25)])
out["check3c_data_mirror"] = dict(x1_max=float(np.max(np.abs(x + xm))), x1hat_max=float(np.max(np.abs(xs + xsm))),
                                  A_max=float(np.max(np.abs(Ao + Am))))
apf_obs = np.array([APF[n] for n in C1])
out["check3c_r_unchanged"] = dict(ch1=r(xs, apf_obs), ch1_mirror=r(xsm, -apf_obs),
                                  ch2=r(x, [AN[n + 1] for n in C2]), ch2_mirror=r(xm, [-AN[n + 1] for n in C2]))
print(out["check3c_data_mirror"], out["check3c_r_unchanged"])
dump("v0_checks.json", out)
tic("done")
