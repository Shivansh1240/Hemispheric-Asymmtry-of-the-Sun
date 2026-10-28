"""closure_tests.py -- closure tests T1, T2, T5 (data and stored runs) and the T3 / T4 summaries (closure_dynamo.py runs).
New file; nothing existing is changed.  Plain statistics: Pearson r with n, sign counts, run lengths, 16-84 % bands.

T1  N-S phase lag persistence.  lag(n) = t_N - t_S, the emergence-weighted mean time of each hemisphere in cycle n
    (Sun: out/explore.json E_timing, n = 13; dynamo: q4_dynamo_analyse.per_cycle 'lag', the same definition, for every
    existing production setting).  Lag-1 r of the lag series and the runs of its sign (13-cycle segments in the model).
T2  firm pairs.  Region bootstrap of A(n) for all 13 cycles SC12-24 (regions resampled with replacement, all daily rows
    of a drawn region kept, 1000 draws; the method of ens_extra.py).  Firm: the observed sign in >= 90 % of draws.
    Persistence hits / pairs over consecutive firm-firm pairs; the number of the 12 pairs that involve a non-firm cycle.
T5  AR(2) check.  From the stored fit (a, b): root modulus and period in cycles.  Then 13-cycle series simulated from the
    stored AR(1) fit, AR(2) refitted exactly as explore.py does (11 triples), 16-84 % band of the fitted b.

    python closure_tests.py t1 t2 t5      -> out/closure/T1.json, T2.json, T5.json
    python closure_tests.py t3t4          -> out/closure/T3.json, T4.json  (after closure_dynamo.py launch)
"""
import os, sys, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
OUT = os.path.join(HERE, "out"); OUTC = os.path.join(OUT, "closure"); os.makedirs(OUTC, exist_ok=True)

def r(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[m], y[m])[0, 1]) if m.sum() > 2 else float("nan")

def runs(seq):
    return [len(list(g)) for _, g in itertools.groupby(seq)]

def seg_runs(sg, L=13):
    seg = [runs(list(sg[i:i + L])) for i in range(0, len(sg) - L + 1)]
    return dict(mean_runs_per_13=float(np.mean([len(s) for s in seg])), mean_longest_per_13=float(np.mean([max(s) for s in seg])))

def J(name):
    return json.load(open(os.path.join(OUT, name)))

def dump(name, d):
    json.dump(d, open(os.path.join(OUTC, name), "w"), indent=1, default=float)
    print(f"wrote out/closure/{name}")

# existing production settings, as in q4_dynamo_ensemble.py (copied: importing that script would re-run it)
SETTINGS = {
    "kappa=1": ["prod_k1", "prod_k1s2", "prod_k1_s3", "prod_k1_s4", "prod_k1_s5"],
    "rho=0.8": ["prod_rho08", "prod_rho08_s2", "prod_rho08_s3", "prod_rho08_s4", "prod_rho08_s5"],
    "tau_c=3 yr": ["prod_tc3_s1", "prod_tc3_s2", "prod_tc3_s3"],
    "tau_c=5 yr": ["prod_tc5_s1", "prod_tc5_s2", "prod_tc5_s3"],
    "v0=10 m/s": ["prod_v10_s1", "prod_v10_s2", "prod_v10_s3"],
    "kappa=0.3": ["prod_k03"], "kappa=3": ["prod_k3"], "rho=0.5": ["prod_rho05"],
    "kappa=1 (49x91)": ["prod_k1lr"], "kappa=10 (49x91)": ["prod_k10lr"],
}

def model_rows(name, spin=None):
    from q4_dynamo_analyse import per_cycle
    from closure_dynamo import load
    if spin is None:
        spin = 150.0 if "v10" in name else 100.0
    return per_cycle(load(name), spin=spin)

def pooled(names, key):
    """pooled lag-1 r of a per-cycle quantity over all adjacent pairs of all seeds, plus per-seed run statistics."""
    x0, x1, per = [], [], []
    for n in names:
        v = np.array([c[key] for c in model_rows(n)], float)
        x0 += list(v[:-1]); x1 += list(v[1:])
        per.append(dict(run=n, n_cycles=len(v), lag1=r(v[:-1], v[1:]), **seg_runs(np.sign(v))))
    return dict(n_seeds=len(names), n_pairs=len(x0), lag1_pooled=r(x0, x1),
                lag1_seed_range=[min(p["lag1"] for p in per), max(p["lag1"] for p in per)],
                mean_longest_per_13=float(np.median([p["mean_longest_per_13"] for p in per])),
                mean_longest_seed_range=[min(p["mean_longest_per_13"] for p in per), max(p["mean_longest_per_13"] for p in per)],
                mean_runs_per_13=float(np.median([p["mean_runs_per_13"] for p in per])), per_seed=per)

# ------------------------------------------------------------------------------------------------------------ T1
def t1():
    E = J("explore.json")["E_timing"]["per_cycle"]; cyc = sorted(int(k) for k in E)
    lag = np.array([E[str(n)]["lag_N_minus_S"] for n in cyc]); s = J("plot_series.json")
    A = np.array([s["A"][str(n)] for n in cyc]); sg = np.sign(lag)
    ind = J("q3_persistence.json")["runs"]["area_all_SC12-24"]["independent_ref"]
    sun = dict(cycles=cyc, lag=lag.tolist(), n=len(lag), lag1_r=r(lag[:-1], lag[1:]), n_pairs=len(lag) - 1,
               signs=sg.astype(int).tolist(), run_lengths=runs(list(sg)), n_runs=len(runs(list(sg))), longest=max(runs(list(sg))),
               lag1_r_A=r(A[:-1], A[1:]), A_run_lengths=runs(list(np.sign(A))),
               same_sign_lag_A=int((np.sign(lag) == np.sign(A)).sum()), same_sign_lag_nextA=int((np.sign(lag[:-1]) == np.sign(A[1:])).sum()),
               independent_ref=ind)
    print("T1 Sun:", {k: v for k, v in sun.items() if k not in ("lag", "signs", "cycles")})
    mod = {}
    for lab, names in SETTINGS.items():
        mod[lab] = dict(lag=pooled(names, "lag"), A=pooled(names, "A"))
        print(f"T1 {lab}: lag lag-1 {mod[lab]['lag']['lag1_pooled']:+.2f} (n = {mod[lab]['lag']['n_pairs']}), longest {mod[lab]['lag']['mean_longest_per_13']:.1f}; "
              f"A lag-1 {mod[lab]['A']['lag1_pooled']:+.2f}")
    dump("T1.json", dict(sun=sun, model=mod))

# ------------------------------------------------------------------------------------------------------------ T2
def t2(ndraw=1000, firm=0.90):
    from core import D, SCB, AN
    rng = np.random.default_rng(2026)
    north = (D.hemi == "N").values; area = D.area.values; g = D.gid.values; dec = D.dec.values
    cyc = list(range(12, 25)); rows = {}
    for n in cyc:
        m = (dec >= SCB[n]) & (dec < SCB[n + 1]); gi = g[m]; a = area[m]; nh = north[m]
        ids, inv = np.unique(gi, return_inverse=True)
        aN = np.bincount(inv, weights=np.where(nh, a, 0.0)); aS = np.bincount(inv, weights=np.where(nh, 0.0, a))
        vals = []
        for _ in range(ndraw):
            k = rng.integers(0, len(ids), len(ids)); N_, S_ = aN[k].sum(), aS[k].sum(); vals.append((N_ - S_) / (N_ + S_))
        vals = np.array(vals); A0 = float((aN.sum() - aS.sum()) / (aN.sum() + aS.sum()))
        assert abs(A0 - AN[n]) < 1e-9, (n, A0, AN[n])
        f = float(np.mean(np.sign(vals) == np.sign(A0)))
        rows[n] = dict(A=A0, frac_same_sign=f, firm=bool(f >= firm), p16_84=np.percentile(vals, [16, 84]).tolist(), n_regions=int(len(ids)))
        print(f"T2 SC{n}: A {A0:+.3f}, observed sign in {f * 100:.1f} % of draws, {'firm' if f >= firm else 'not firm'}")
    pairs = [(n, n + 1) for n in cyc[:-1]]
    ff = [(a, b) for a, b in pairs if rows[a]["firm"] and rows[b]["firm"]]
    hits = [(a, b) for a, b in ff if np.sign(rows[a]["A"]) == np.sign(rows[b]["A"])]
    allhits = [(a, b) for a, b in pairs if np.sign(rows[a]["A"]) == np.sign(rows[b]["A"])]
    old = J("ens_extra.json")
    repro = {n: dict(previous=old[f"boot_SC{n}"]["frac_same_sign"], now=rows[n]["frac_same_sign"]) for n in (13, 14, 17, 20, 21, 24)}
    out = dict(n_draws=ndraw, firm_threshold=firm, cycles=rows, n_firm=int(sum(v["firm"] for v in rows.values())),
               firm_cycles=[n for n in cyc if rows[n]["firm"]], nonfirm_cycles=[n for n in cyc if not rows[n]["firm"]],
               firm_firm_pairs=ff, firm_firm_hits=len(hits), firm_firm_n=len(ff), firm_firm_breaks=[p for p in ff if p not in hits],
               pairs_with_nonfirm=len(pairs) - len(ff), all_pairs_hits=len(allhits), all_pairs_n=len(pairs),
               reproduces_previous_bootstrap=repro)
    print(f"T2: firm cycles {out['firm_cycles']}; firm-firm pairs {len(hits)}/{len(ff)} persist; {out['pairs_with_nonfirm']} of 12 pairs involve a non-firm cycle")
    dump("T2.json", out)

# ------------------------------------------------------------------------------------------------------------ T5
def t5(nsim=20000):
    Dm = J("explore.json")["D_memory"]; a, b = Dm["a"], Dm["b"]
    roots = np.roots([1.0, -a, -b]); mod = float(np.abs(roots[0])); ang = float(np.abs(np.angle(roots[0])))
    period = float(2 * np.pi / ang) if ang > 0 else float("inf")
    s = J("plot_series.json"); AN = {int(k): v for k, v in s["A"].items()}
    nn = list(range(13, 24)); y = np.array([AN[n + 1] for n in nn]); x1 = np.array([AN[n] for n in nn]); x2 = np.array([AN[n - 1] for n in nn])
    X1d = np.column_stack([np.ones(len(nn)), x1]); c1, *_ = np.linalg.lstsq(X1d, y, rcond=None); sd1 = float(np.std(y - X1d @ c1, ddof=2))
    Xd = np.column_stack([np.ones(len(nn)), x1, x2]); c2, *_ = np.linalg.lstsq(Xd, y, rcond=None)
    assert abs(c1[1] - Dm["ar1_a"]) < 1e-12 and abs(sd1 - Dm["ar1_resid_sd"]) < 1e-12 and abs(c2[2] - b) < 1e-12   # same fits as explore.py
    rng = np.random.default_rng(5)
    bs, as_ = [], []
    for _ in range(nsim):
        z = [AN[12], AN[13]]
        for k in range(200 + 13):
            z.append(c1[0] + c1[1] * z[-1] + rng.normal(0, sd1))
        z = np.array(z[-13:])
        yy = z[2:]; X = np.column_stack([np.ones(11), z[1:-1], z[:-2]])
        cf, *_ = np.linalg.lstsq(X, yy, rcond=None); as_.append(cf[1]); bs.append(cf[2])
    bs = np.array(bs); as_ = np.array(as_)
    p16, p50, p84 = np.percentile(bs, [16, 50, 84])
    out = dict(a=a, b=b, roots=[[float(z.real), float(z.imag)] for z in roots], root_modulus=mod, angle_deg=float(np.degrees(ang)), period_cycles=period,
               ar1=dict(c=float(c1[0]), a=float(c1[1]), resid_sd=sd1), nsim=nsim, n_triples=11,
               b_fit_p16_50_84=[float(p16), float(p50), float(p84)], a_fit_p16_50_84=np.percentile(as_, [16, 50, 84]).tolist(),
               b_obs_inside_16_84=bool(p16 <= b <= p84), b_fit_p2_98=np.percentile(bs, [2.5, 97.5]).tolist())
    print(f"T5: roots {roots}, modulus {mod:.3f}, period {period:.2f} cycles; AR(1)-simulated b 16-84 % [{p16:+.2f}, {p84:+.2f}] (median {p50:+.2f}); "
          f"observed {b:+.2f} inside: {out['b_obs_inside_16_84']}")
    old = os.path.join(OUTC, "T5.json")
    if os.path.exists(old):                                   # deterministic (seed 5): a rerun must give the same band
        assert np.allclose(json.load(open(old))["b_fit_p16_50_84"], out["b_fit_p16_50_84"], atol=1e-12)
    np.savez_compressed(os.path.join(OUTC, "T5_samples.npz"), b_fit=bs, a_fit=as_)
    dump("T5.json", out)

# ------------------------------------------------------------------------------------------------------- T3 / T4
def summary(names, spin=100.0):
    from q4_dynamo_analyse import per_cycle, stats_sunlike
    A0, A1, per, allA, allApf, sl, sgn_apf, sgn_A = [], [], [], [], [], [], [], []
    for n in names:
        rows = per_cycle(__import__("closure_dynamo").load(n), spin=spin)
        A = np.array([c["A"] for c in rows]); Apf = np.array([c["Apf"] for c in rows])
        A0 += list(A[:-1]); A1 += list(A[1:]); allA += list(A); allApf += list(Apf)
        sl.append(stats_sunlike(rows)["frac_sunlike"])
        per.append(dict(run=n, n_cycles=len(A), lag1=r(A[:-1], A[1:]), rms_A=float(np.sqrt(np.mean(A ** 2))),
                        rms_Apf=float(np.sqrt(np.mean(Apf ** 2))), **seg_runs(np.sign(A))))
    allA = np.array(allA); allApf = np.array(allApf)
    return dict(runs=names, n_cycles=len(allA), n_pairs=len(A0), lag1_A_pooled=r(A0, A1),
                lag1_seed_range=[min(p["lag1"] for p in per), max(p["lag1"] for p in per)],
                mean_longest_per_13=float(np.median([p["mean_longest_per_13"] for p in per])),
                mean_longest_seed_range=[min(p["mean_longest_per_13"] for p in per), max(p["mean_longest_per_13"] for p in per)],
                rms_A=float(np.sqrt(np.mean(allA ** 2))), rms_Apf=float(np.sqrt(np.mean(allApf ** 2))),
                rms_Apf_seed_median=float(np.median([p["rms_Apf"] for p in per])),
                frac_sunlike=float(np.median(sl)), frac_sunlike_range=[min(sl), max(sl)],
                mean_A=float(np.mean(allA)), mean_Apf=float(np.mean(allApf)),
                frac_A_pos=float(np.mean(allA > 0)), frac_Apf_pos=float(np.mean(allApf > 0)), per_seed=per)

def t3t4():
    t3(); t4()

def t3():
    import closure_dynamo as CD
    cal = json.load(open(os.path.join(OUTC, "T3_calibration.json")))
    t3 = {"1": dict(sigma=0.6, reused=CD.REF_TAU1, **summary(CD.REF_TAU1))}
    for tc in (3.0, 5.0, 10.0, 30.0):
        names = [f"T3_tc{tc:g}_s{s}" for s in CD.SEEDS]
        t3[f"{tc:g}"] = dict(sigma=cal[f"{tc:g}"]["sigma_prod"], calibration=cal[f"{tc:g}"], **summary(names))
    for k, v in t3.items():
        print(f"T3 tau_c {k}: sigma {v['sigma']:.3f}, rms A {v['rms_A']:.3f}, lag-1 {v['lag1_A_pooled']:+.2f} (n = {v['n_pairs']}), "
              f"longest {v['mean_longest_per_13']:.1f}, Sun-like {v['frac_sunlike'] * 100:.0f} %, rms Apf {v['rms_Apf']:.2f}")
    dump("T3.json", t3)

def t4():
    import closure_dynamo as CD
    t4 = {"reference (delta 0)": summary(CD.REF_TAU1)}
    for kind in CD.KINDS:
        for dl in CD.DELTAS:
            t4[f"{kind} {int(round(dl * 100))} %"] = summary([f"T4_{kind}_d{int(round(dl * 100))}_s{s}" for s in CD.SEEDS])
    for k, v in t4.items():
        print(f"T4 {k}: lag-1 {v['lag1_A_pooled']:+.2f} (n = {v['n_pairs']}), longest {v['mean_longest_per_13']:.1f}, mean A {v['mean_A']:+.3f}, "
              f"A>0 {v['frac_A_pos'] * 100:.0f} %, mean Apf {v['mean_Apf']:+.3f}, Apf>0 {v['frac_Apf_pos'] * 100:.0f} %, rms Apf {v['rms_Apf']:.2f}")
    s = J("plot_series.json"); apf = np.array([s["Apf"][str(n)] for n in range(15, 25)])
    t4["_sun"] = dict(mean_Apf=float(np.mean(apf)), frac_Apf_pos=float(np.mean(apf > 0)), n=len(apf))
    dump("T4.json", t4)

if __name__ == "__main__":
    for a in sys.argv[1:]:
        globals()[a]()
