"""q4_dynamo_ensemble.py -- the Q4c table over seeds.

For every setting: the per-seed values (median and range over seeds) of the same columns as the single-run table, and
pooled Pearson r over all adjacent cycle pairs of all seeds (one r, with its n).  Plain statistics only."""
import os, json
import numpy as np
from q4_dynamo import OUT
from q4_dynamo_analyse import per_cycle, stats, stats_sunlike, r

SETTINGS = {
    "kappa=1": ["prod_k1", "prod_k1s2", "prod_k1_s3", "prod_k1_s4", "prod_k1_s5"],
    "rho=0.8": ["prod_rho08", "prod_rho08_s2", "prod_rho08_s3", "prod_rho08_s4", "prod_rho08_s5"],
    "tau_c=3 yr": ["prod_tc3_s1", "prod_tc3_s2", "prod_tc3_s3"],
    "tau_c=5 yr": ["prod_tc5_s1", "prod_tc5_s2", "prod_tc5_s3"],
    "v0=10 m/s": ["prod_v10_s1", "prod_v10_s2", "prod_v10_s3"],
    "kappa=0.3": ["prod_k03"], "kappa=3": ["prod_k3"], "rho=0.5": ["prod_rho05"],
    "kappa=1 (49x91)": ["prod_k1lr"], "kappa=10 (49x91)": ["prod_k10lr"],
}
COLS = ["period", "rms_A", "rms_Apf", "rms_lag", "lag1_A", "lag2_A", "longest_run_mean_per_13", "r_absApf_absnextA", "r_Apf_next_A",
        "r_Xmax_Apf_paritysigned", "r_Xdec_Apf_unsigned"]

out = {}
for lab, names in SETTINGS.items():
    have = [n for n in names if os.path.exists(os.path.join(OUT, f"{n}.npz"))]
    if not have:
        continue
    per, pool = [], dict(A0=[], A1=[], A2=[], P0=[], Aabs1=[], n_cycles=0)
    sun = []
    for n in have:
        res = dict(np.load(os.path.join(OUT, f"{n}.npz"), allow_pickle=True))
        spin = 150.0 if "v10" in n else 100.0               # the 10 m/s dynamo saturates more slowly (tune_v10)
        rows = per_cycle(res, spin=spin); st = stats(rows); per.append(st); sun.append(stats_sunlike(rows)["frac_sunlike"])
        A = np.array([c["A"] for c in rows]); Apf = np.array([c["Apf"] for c in rows])
        pool["A0"] += list(A[:-1]); pool["A1"] += list(A[1:]); pool["P0"] += list(Apf[:-1]); pool["n_cycles"] += len(A)
        pool["A2"] += list(zip(A[:-2], A[2:]))
    pooled = dict(n_seeds=len(have), n_cycles=pool["n_cycles"], n_pairs=len(pool["A0"]),
                  lag1_A=r(pool["A0"], pool["A1"]), lag2_A=r([a for a, _ in pool["A2"]], [b for _, b in pool["A2"]]),
                  r_absApf_absnextA=r(np.abs(pool["P0"]), np.abs(pool["A1"])), r_Apf_next_A=r(pool["P0"], pool["A1"]))
    summ = {c: dict(median=float(np.median([p[c] for p in per])), min=float(np.min([p[c] for p in per])), max=float(np.max([p[c] for p in per])))
            for c in COLS}
    summ["frac_sunlike"] = dict(median=float(np.median(sun)), min=float(np.min(sun)), max=float(np.max(sun)))
    out[lab] = dict(runs=have, per_seed=summ, pooled=pooled)
    print(f"\n{lab}: {len(have)} seeds, {pooled['n_cycles']} cycles")
    for c in ["frac_sunlike"] + COLS:
        v = summ[c]; print(f"   {c:28s} median {v['median']:+.2f}  range [{v['min']:+.2f}, {v['max']:+.2f}]")
    print(f"   POOLED lag-1 r(A) {pooled['lag1_A']:+.2f} (n = {pooled['n_pairs']}), lag-2 {pooled['lag2_A']:+.2f}, "
          f"r(|Apf|,|A(n+1)|) {pooled['r_absApf_absnextA']:+.2f}, r(Apf, A(n+1)) {pooled['r_Apf_next_A']:+.2f}")
json.dump(out, open(os.path.join(OUT, "ensemble_table.json"), "w"), indent=1)

# markdown table, same columns as REPORT Q4c (written from the numbers above, nothing typed by hand)
lines = ["| Setting | seeds (cycles) | Sun-like cycles | rms A | rms A_pf | rms N-S lag (yr) | lag-1 r of A, pooled (seed range) | mean longest run / 13 (range) | r(\\|A_pf(n)\\|, \\|A(n+1)\\|), pooled (range) | r(A_pf(n), A(n+1)), pooled (range) | period (yr) |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
for lab, v in out.items():
    ps, po = v["per_seed"], v["pooled"]; one = po["n_seeds"] == 1
    rng = lambda c: "" if one else f" ({ps[c]['min']:+.2f} to {ps[c]['max']:+.2f})"
    lines.append(f"| {lab} | {po['n_seeds']} ({po['n_cycles']}) | {ps['frac_sunlike']['median']*100:.0f} % | {ps['rms_A']['median']:.2f} | {ps['rms_Apf']['median']:.2f} | "
                 f"{ps['rms_lag']['median']:.2f} | {po['lag1_A']:+.2f}{rng('lag1_A')} (n = {po['n_pairs']}) | {ps['longest_run_mean_per_13']['median']:.1f}"
                 + ("" if one else f" ({ps['longest_run_mean_per_13']['min']:.1f}-{ps['longest_run_mean_per_13']['max']:.1f})")
                 + f" | {po['r_absApf_absnextA']:+.2f}{rng('r_absApf_absnextA')} | {po['r_Apf_next_A']:+.2f}{rng('r_Apf_next_A')} | {ps['period']['median']:.1f} |")
open(os.path.join(OUT, "ensemble_table.md"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
print("\n".join(lines))
