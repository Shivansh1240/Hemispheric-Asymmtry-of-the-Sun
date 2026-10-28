"""q4_dynamo_checks.py -- checks of the self-excited dynamo before its statistics are used.

0  the declining-phase solver's own eight checks (q4_ftd_verify.py): 8/8 must pass
1  SYMMETRY   no noise, dipole seed: the hemispheres must stay identical (|E_N - E_S| / E at round-off)
2  MIRROR     the same noise sequence with the hemispheres' forcing swapped must return the exact mirror image:
              E_N(b) = E_S(a), E_S(b) = E_N(a), |P_N(b)| = |P_S(a)|
3  RESOLUTION 97 x 181 against 65 x 121 (same parameters and seed): the cycle statistics must agree within their
              own sampling scatter (the trajectories need not, the system is chaotic under noise)
"""
import os, json
import numpy as np
from q4_dynamo import OUT
from q4_dynamo_analyse import per_cycle, stats

L = lambda n: dict(np.load(os.path.join(OUT, f"{n}.npz"), allow_pickle=True))
out = {}
s = L("t3_w35_k1_nonoise"); E = s["EN"] + s["ES"]
out["symmetry_max_rel"] = float(np.max(np.abs(s["EN"] - s["ES"])) / E.max())
a, b = L("chk_mirror_a"), L("chk_mirror_b")
out["mirror_E_max_rel"] = float(max(np.max(np.abs(a["EN"] - b["ES"])), np.max(np.abs(a["ES"] - b["EN"]))) / np.max(a["EN"] + a["ES"]))
out["mirror_P_max_rel"] = float(max(np.max(np.abs(np.abs(a["PN"]) - np.abs(b["PS"]))), np.max(np.abs(np.abs(a["PS"]) - np.abs(b["PN"])))) / np.max(np.abs(a["PN"])))
out["mirror_signs_PN_b_vs_PS_a"] = float(np.mean(np.sign(b["PN"]) == np.sign(a["PS"])))
hi = stats(per_cycle(L("chk_res97"))); lo_res = L("prod_k1")
m = lo_res["t"] <= 400
lo = stats(per_cycle({k: (v[m] if hasattr(v, "shape") and v.shape[:1] == lo_res["t"].shape else v) for k, v in lo_res.items()}))
keys = ("n_cycles", "period", "rms_A", "rms_Apf", "rms_lag", "lag1_A", "lag1_Apf")
out["resolution"] = {k: dict(r97x181=hi[k], r65x121=lo[k]) for k in keys}
for k, v in out.items():
    print(k, v)
json.dump(out, open(os.path.join(OUT, "dynamo_checks.json"), "w"), indent=1, default=float)
