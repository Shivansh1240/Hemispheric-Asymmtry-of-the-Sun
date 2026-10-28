"""ens_extra.py -- splits the H1 ensemble spread into its two sources, and adds the sampling (region-bootstrap)
spread of the flip-cycle asymmetries, which measurement-noise ensembles cannot show.

H1 (a) area factors only (A_pf fixed);  (b) A_pf redrawn within its published 1 sigma only (areas fixed).
Note that (b) adds a second error realisation to values that already carry one, so it shows how fragile r is to
errors of the published size; it is not an unbiased estimate of the true r.
Region bootstrap: the regions of a cycle resampled with replacement (all daily rows of a drawn region kept),
1000 draws, for SC13, 14, 20, 21, 17, 24.
"""
import numpy as np
from core import *

out = {}
rng = np.random.default_rng(7)
apf = np.array([APF[n] for n in C1]); sig = np.array([SIG[n] for n in C1]); xh = x1v(C1, 0, 1, True)
E = np.load(os.path.join(OUT, "ens_data.npz"))
ra = [r(x, apf) for x in E["XHE"]]
rb = [r(xh, apf + rng.normal(0, 1, 10) * sig) for _ in range(1000)]
out["H1_area_only"] = dict(median=float(np.median(ra)), p16_84=np.percentile(ra, [16, 84]).tolist(), n=10)
out["H1_apf_error_only"] = dict(median=float(np.median(rb)), p16_84=np.percentile(rb, [16, 84]).tolist(), n=10)
print(out)

tic("region bootstrap")
north = (D.hemi == "N").values; area = D.area.values; g = D.gid.values; dec = D.dec.values
for n in (13, 14, 17, 20, 21, 24):
    m = (dec >= SCB[n]) & (dec < SCB[n + 1]); gi = g[m]; a = area[m]; nh = north[m]
    ids, inv = np.unique(gi, return_inverse=True)
    aN = np.bincount(inv, weights=np.where(nh, a, 0.0)); aS = np.bincount(inv, weights=np.where(nh, 0.0, a))
    LL = band_mask(D[m], 500, np.inf, 0, 20).values; isLL = np.bincount(inv, weights=LL.astype(float)) > 0
    vals, valsLL = [], []
    for _ in range(1000):
        k = rng.integers(0, len(ids), len(ids)); N_, S_ = aN[k].sum(), aS[k].sum(); vals.append((N_ - S_) / (N_ + S_))
        kk = k[isLL[k]]; NL, SL = aN[kk].sum(), aS[kk].sum(); valsLL.append((NL - SL) / (NL + SL))
    vals = np.array(vals); valsLL = np.array(valsLL)
    out[f"boot_SC{n}"] = dict(A=AN[n], p16_84=np.percentile(vals, [16, 84]).tolist(), p2_98=np.percentile(vals, [2.5, 97.5]).tolist(),
                              frac_same_sign=float(np.mean(np.sign(vals) == np.sign(AN[n]))), n_regions=int(len(ids)),
                              lowlat_large_p16_84=np.percentile(valsLL, [16, 84]).tolist(), lowlat_large_frac_same_sign_as_total=float(np.mean(np.sign(valsLL) == np.sign(AN[n]))))
    print(n, out[f"boot_SC{n}"])
dump("ens_extra.json", out)
