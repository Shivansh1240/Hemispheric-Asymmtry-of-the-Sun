"""surya_v8.py -- the earlier dynamo experiments (Q4 F3/F4, T1, T3, T4) redone in Bhowmik's dynamo
(surya/surya_bhowmik_ext.f: Passos et al. 2014 on the public SURYA, plus pumping; v7) instead of our 2-D model.
New file; surya_runs.py and the v7 numbers are not changed.

Settings mirror the earlier runs (q4_dynamo.py, closure_dynamo.py), translated to Surya:
  production   Passos noise (BL 100 %, MF 200 %, uniform, renewed every 0.5 yr), rho 0, kappa 1, v0 29 m/s.
               The kappa = 1 reference is the v7 ensemble (out/v7_surya, 6 seeds).
  Q4 F4        rho 0.8 (3 seeds), rho 0.5 (1), kappa 0.3 and 3 (1 each), v0 20 m/s (2/3 of 29, as 10 of 15 m/s; 3 seeds),
               tau_c 3 and 5 yr at the production noise (3 seeds each).
  Q4 F3        one production run with the butterfly written (but.dat).
  T3           tau_c = 0.5, 1, 3, 5, 10, 30 yr; both noise levels scaled by one factor s per tau_c so that rms A is
               0.10 (the Sun's).  One calibration step: two 250-yr runs at a first-guess s, then
               s_prod = s_cal * 0.10 / rms A_cal.  3 seeds each.
  T4           constant N-S bias delta = 1, 3, 10 % in alpha_BL, the flow or eta (definitions in the Fortran header),
               production noise, 3 seeds each.
Runs start from the noise-free saturated state (out/v7_surya/det/final.dat); the first `spin` years are dropped.

    python surya_v8.py validate     8 short mirror runs -> out/v8_surya/validation.json
    python surya_v8.py launch       validation, then the whole batch (9 at a time; T3 production after calibration)
    python surya_v8.py analyse      out/v8_surya/summary.json
"""
import os, sys, json, shutil, subprocess, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
EXE = os.path.join(HERE, "surya", "surya_bhowmik_ext.exe")
INIT = os.path.join(HERE, "out", "v7_surya", "det", "final.dat")
OUT = os.path.join(HERE, "out", "v8_surya"); os.makedirs(OUT, exist_ok=True)
YR = 3600.0 * 24 * 365 / 1e8
NCONC = 9
PROD = dict(f0r=3.6, flvBL=1.0, flvMF=2.0, iswap=0, tauc=0.5, rho=0.0, xkap=1.0, v0abs=29.0, kbias=0, delta=0.0, ibut=0)
YEARS, SPIN = 450, 50                       # 400 analysed years, about 37 cycles per run
YEARS_SLOW, SPIN_SLOW = 500, 100            # flow and coupling changes move the period: longer transient dropped
SEEDS = [1, 2, 3]
TAUS = [0.5, 1.0, 3.0, 5.0, 10.0, 30.0]
S_GUESS = {0.5: 0.55, 1.0: 0.40, 3.0: 0.25, 5.0: 0.19, 10.0: 0.14, 30.0: 0.11}
CAL_SEEDS, CAL_YEARS = [11, 12], 250
TARGET_RMS_A = 0.10
DELTAS = [0.01, 0.03, 0.10]
KINDS = {"alpha": 1, "v0": 2, "eta": 3}
LOG = os.path.join(OUT, "launch.log")

def cfg(**kw):
    c = dict(PROD, years=YEARS, spin=SPIN); c.update(kw); return c

def val_runs():
    R = {}
    knobs = dict(tauc=1.0, rho=0.5, xkap=3.0, v0abs=20.0)
    R["val_knobs_0"] = cfg(years=5, spin=0, seed=7, **knobs)
    R["val_knobs_1"] = cfg(years=5, spin=0, seed=7, iswap=1, **knobs)
    for k, kb in KINDS.items():
        R[f"val_{k}_p"] = cfg(years=5, spin=0, seed=7, kbias=kb, delta=0.10)
        R[f"val_{k}_m"] = cfg(years=5, spin=0, seed=7, kbias=kb, delta=-0.10, iswap=1)
    return R

def batch_runs():
    """everything except the T3 production runs, in queue order"""
    R = {}
    for tc in TAUS:
        for s in CAL_SEEDS:
            g = S_GUESS[tc]
            R[f"T3cal_tc{tc:g}_s{s}"] = cfg(years=CAL_YEARS, seed=s, tauc=tc, flvBL=g, flvMF=2 * g)
    R["prod_but_s1"] = cfg(seed=1, ibut=1)
    for s in SEEDS:
        R[f"rho08_s{s}"] = cfg(seed=s, rho=0.8)
    R["rho05_s1"] = cfg(seed=1, rho=0.5)
    for s in SEEDS:
        R[f"v20_s{s}"] = cfg(seed=s, v0abs=20.0, years=YEARS_SLOW, spin=SPIN_SLOW)
    R["kap03_s1"] = cfg(seed=1, xkap=0.3, years=YEARS_SLOW, spin=SPIN_SLOW)
    R["kap3_s1"] = cfg(seed=1, xkap=3.0, years=YEARS_SLOW, spin=SPIN_SLOW)
    for tc in (3.0, 5.0):
        for s in SEEDS:
            R[f"tc{tc:g}_s{s}"] = cfg(seed=s, tauc=tc)
    return R

def t4_runs():
    return {f"T4_{k}_d{int(round(d * 100))}_s{s}": cfg(seed=s, kbias=kb, delta=d)
            for k, kb in KINDS.items() for d in DELTAS for s in SEEDS}

def t3_runs(cal):
    R = {}
    for tc in TAUS:
        sp = cal[f"{tc:g}"]["s_prod"]
        for s in SEEDS:
            R[f"T3_tc{tc:g}_s{s}"] = cfg(seed=s, tauc=tc, flvBL=sp, flvMF=2 * sp)
    return R

def prepare(name, c):
    d = os.path.join(OUT, name); os.makedirs(d, exist_ok=True)
    shutil.copy(INIT, os.path.join(d, "init.dat"))
    with open(os.path.join(d, "run.in"), "w") as f:
        f.write(f"{c['years'] * YR:.6f} {c['seed']} {c['f0r']} {c['flvBL']:.6f} {c['flvMF']:.6f} {c['iswap']} "
                f"{c['tauc']} {c['rho']} {c['xkap']} {c['v0abs']} {c['kbias']} {c['delta']} {c['ibut']}\n")
    json.dump(c, open(os.path.join(d, "config.json"), "w"))
    return d

def done(name):
    return os.path.exists(os.path.join(OUT, name, "final.dat"))

def log(msg):
    line = f"{time.strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    open(LOG, "a").write(line + "\n")

def run_queue(queue, on_finish=None):
    """queue: list of (name, config); on_finish(name) may return more (name, config) items to put at the front."""
    procs = {}
    queue = [q for q in queue if not done(q[0])]
    while queue or procs:
        while queue and len(procs) < NCONC:
            name, c = queue.pop(0)
            d = prepare(name, c)
            procs[name] = subprocess.Popen([EXE], cwd=d, stdout=open(os.path.join(d, "stdout.txt"), "w"), stderr=subprocess.STDOUT)
            log(f"start {name} ({c['years']} yr)")
        time.sleep(20)
        for name in [n for n, p in procs.items() if p.poll() is not None]:
            rc = procs.pop(name).returncode
            log(f"finish {name}: exit {rc}, {len(procs)} running, {len(queue)} queued")
            if on_finish:
                queue[0:0] = [q for q in (on_finish(name) or []) if not done(q[0])]

def load(name):
    a = np.loadtxt(os.path.join(OUT, name, "polar.dat"))
    return dict(t=a[:, 0], PN=a[:, 1], PS=a[:, 2], EN=a[:, 3], ES=a[:, 4], BN=a[:, 5], BS=a[:, 6])

def check_validation():
    out = {}
    a, b = load("val_knobs_0"), load("val_knobs_1")
    out["knobs mirror (tau_c 1, rho 0.5, kappa 3, v0 20)"] = dict(
        PN_vs_swapped_PS=float(np.max(np.abs(a["PN"] + b["PS"])) / np.max(np.abs(a["PN"]))),
        EN_vs_swapped_ES=float(np.max(np.abs(a["EN"] - b["ES"])) / max(np.max(a["EN"]), 1e-30)))
    for k in KINDS:
        a, b = load(f"val_{k}_p"), load(f"val_{k}_m")
        out[f"{k} bias +10 % vs -10 % swapped"] = dict(
            PN_vs_swapped_PS=float(np.max(np.abs(a["PN"] + b["PS"])) / np.max(np.abs(a["PN"]))),
            EN_vs_swapped_ES=float(np.max(np.abs(a["EN"] - b["ES"])) / max(np.max(a["EN"]), 1e-30)))
    out["zero-knob identity"] = "surya_bhowmik_ext.exe with tauc 0.5, rho 0, xkap 1, v0 29, kbias 0 reproduces surya_bhowmik.exe bit for bit (2-yr run, seed 7: polar.dat and final.dat identical; out/v8_surya/_check)"
    json.dump(out, open(os.path.join(OUT, "validation.json"), "w"), indent=1)
    return out

def validate():
    R = val_runs()
    run_queue(list(R.items()))
    out = check_validation()
    log("validation " + json.dumps({k: v for k, v in out.items() if isinstance(v, dict)}))
    return out

def rms_A(names):
    from q4_dynamo_analyse import per_cycle
    A = []
    for n in names:
        c = json.load(open(os.path.join(OUT, n, "config.json")))
        A += [x["A"] for x in per_cycle(load(n), spin=c["spin"])]
    A = np.array(A)
    return float(np.sqrt(np.mean(A ** 2))), int(len(A))

def calibrate():
    cal = {}
    for tc in TAUS:
        names = [f"T3cal_tc{tc:g}_s{s}" for s in CAL_SEEDS]
        r_, n_ = rms_A(names)
        g = S_GUESS[tc]
        cal[f"{tc:g}"] = dict(s_cal=g, rms_A_cal=r_, n_cycles_cal=n_, s_prod=g * TARGET_RMS_A / r_, runs=names)
    json.dump(cal, open(os.path.join(OUT, "T3_calibration.json"), "w"), indent=1)
    log("calibration " + json.dumps({k: round(v["s_prod"], 3) for k, v in cal.items()}))
    return cal

def launch():
    if not os.path.exists(os.path.join(OUT, "validation.json")):
        validate()
    B = batch_runs(); T4 = t4_runs()
    cal_names = {n for n in B if n.startswith("T3cal")}
    state = dict(t3_added=False)

    def on_finish(name):
        if name in cal_names and not state["t3_added"] and all(done(n) for n in cal_names):
            state["t3_added"] = True
            return list(t3_runs(calibrate()).items())
        return []

    queue = list(B.items()) + list(T4.items())
    if all(done(n) for n in cal_names):
        state["t3_added"] = True
        cal = json.load(open(os.path.join(OUT, "T3_calibration.json"))) if os.path.exists(os.path.join(OUT, "T3_calibration.json")) else calibrate()
        queue = list(B.items()) + list(t3_runs(cal).items()) + list(T4.items())
    log(f"launch: {sum(not done(n) for n, _ in queue)} runs to do")
    run_queue(queue, on_finish)
    log("all runs finished")

# ----------------------------------------------------------------------------------------------- analysis ----
V7 = os.path.join(HERE, "out", "v7_surya")
REF = [f"v7:ens_s{s}" for s in range(1, 7)]                  # Bhowmik production, v7 ensemble (spin 100 yr)
GROUPS = {"kappa=1": REF, "rho=0.8": [f"rho08_s{s}" for s in SEEDS], "tau_c=3 yr": [f"tc3_s{s}" for s in SEEDS],
          "tau_c=5 yr": [f"tc5_s{s}" for s in SEEDS], "v0=20 m/s": [f"v20_s{s}" for s in SEEDS],
          "kappa=0.3": ["kap03_s1"], "kappa=3": ["kap3_s1"], "rho=0.5": ["rho05_s1"], "production (butterfly run)": ["prod_but_s1"]}

def _dir(name):
    return os.path.join(V7, name[3:]) if name.startswith("v7:") else os.path.join(OUT, name)

def _spin(name):
    return 100.0 if name.startswith("v7:") else json.load(open(os.path.join(_dir(name), "config.json")))["spin"]

def _finished(name):
    return os.path.exists(os.path.join(_dir(name), "final.dat"))

def _load(name):
    a = np.loadtxt(os.path.join(_dir(name), "polar.dat"))
    return dict(t=a[:, 0], PN=a[:, 1], PS=a[:, 2], EN=a[:, 3], ES=a[:, 4], BN=a[:, 5], BS=a[:, 6])

def _runs(seq):
    import itertools
    return [len(list(g)) for _, g in itertools.groupby(seq)]

def _seg(sg, L=13):
    seg = [_runs(list(sg[i:i + L])) for i in range(0, len(sg) - L + 1)]
    return (float(np.mean([max(s) for s in seg])), float(np.mean([len(s) for s in seg]))) if seg else (float("nan"), float("nan"))

def summary(names):
    """the Q4 / closure statistics for a list of runs (pooled over adjacent pairs within each run)."""
    from q4_dynamo_analyse import per_cycle, stats, stats_sunlike, r as rr
    names = [n for n in names if _finished(n)]
    if not names:
        return None
    P = {k: ([], []) for k in ("A", "lag", "Apf")}; allA, allApf, allLag, apf0, A1 = [], [], [], [], []
    per = []
    for n in names:
        rows = per_cycle(_load(n), spin=_spin(n)); st = stats(rows); sl = stats_sunlike(rows)
        A = np.array([c["A"] for c in rows]); Apf = np.array([c["Apf"] for c in rows]); lag = np.array([c["lag"] for c in rows])
        for k, v in (("A", A), ("lag", lag), ("Apf", Apf)):
            P[k][0].extend(v[:-1]); P[k][1].extend(v[1:])
        allA += list(A); allApf += list(Apf); allLag += list(lag); apf0 += list(Apf[:-1]); A1 += list(A[1:])
        lA, nA = _seg(np.sign(A)); lL, _ = _seg(np.sign(lag))
        per.append(dict(run=n, n_cycles=len(A), period=st["period"], lag1=st["lag1_A"], longest=lA, runs_per_13=nA,
                        lag1_timing=rr(lag[:-1], lag[1:]), longest_timing=lL, rms_A=st["rms_A"], rms_Apf=st["rms_Apf"],
                        mean_A=float(np.mean(A)), r_Apf_next_A=st["r_Apf_next_A"], frac_sunlike=sl["frac_sunlike"],
                        switch_plus=st["switch_frac_par_plus"], switch_minus=st["switch_frac_par_minus"],
                        A=A.tolist(), Apf=Apf.tolist(), lag=lag.tolist()))
    allA = np.array(allA); allApf = np.array(allApf); allLag = np.array(allLag)
    med = lambda k: float(np.median([p[k] for p in per])); rng_ = lambda k: [min(p[k] for p in per), max(p[k] for p in per)]
    return dict(runs=names, n_runs=len(names), n_cycles=int(len(allA)), n_pairs=len(P["A"][0]),
                lag1_A_pooled=rr(*P["A"]), lag1_seed_range=rng_("lag1"), mean_longest_per_13=med("longest"),
                mean_longest_seed_range=rng_("longest"), rms_A=float(np.sqrt(np.mean(allA ** 2))),
                rms_Apf=float(np.sqrt(np.mean(allApf ** 2))), rms_Apf_seed_median=med("rms_Apf"),
                frac_sunlike=med("frac_sunlike"), frac_sunlike_range=rng_("frac_sunlike"), period=med("period"),
                mean_A=float(np.mean(allA)), mean_Apf=float(np.mean(allApf)), frac_A_pos=float(np.mean(allA > 0)),
                frac_Apf_pos=float(np.mean(allApf > 0)), r_Apf_next_A_pooled=rr(apf0, A1),
                lag1_Apf_pooled=rr(*P["Apf"]), timing_lag1_pooled=rr(*P["lag"]), timing_longest_per_13=med("longest_timing"),
                rms_timing_lag=float(np.sqrt(np.mean(allLag ** 2))), per_seed=per)

def cycles_alt(res, smooth_yr=2.0):
    """detector check: minima of the 2-yr running-mean total eruption at least 0.7 period apart (period from the N-polar
    sign changes).  In SURYA eruptions never stop, so the standard detector's emergence minimum between two reversals is
    shallow and a few boundaries land early or late (3-7 of 35-45 cycles per run)."""
    from scipy.signal import find_peaks
    t, E = res["t"], res["EN"] + res["ES"]; dt = t[1] - t[0]; w = max(1, int(smooth_yr / dt))
    Es = np.convolve(E, np.ones(w) / w, mode="same"); k1 = max(1, int(1.0 / dt))
    zN = np.flatnonzero(np.diff(np.sign(np.convolve(res["PN"], np.ones(k1) / k1, mode="same"))) != 0)
    P = (t[-1] - t[0]) / max(len(zN), 1)
    pk, _ = find_peaks(-Es[w:-w], distance=max(1, int(0.7 * P / dt)))
    return pk + w, Es

def detector_check(groups):
    import q4_dynamo_analyse as QA
    orig = QA.cycles; QA.cycles = cycles_alt
    try:
        out = {}
        for k, names in groups.items():
            s = summary(names)
            if s:
                out[k] = {f: s[f] for f in ("n_cycles", "lag1_A_pooled", "mean_longest_per_13", "rms_A", "rms_Apf", "rms_Apf_seed_median", "frac_sunlike",
                                            "frac_A_pos", "frac_Apf_pos", "timing_lag1_pooled", "period")}
    finally:
        QA.cycles = orig
    return out

def analyse():
    out = dict(validation=json.load(open(os.path.join(OUT, "validation.json"))))
    out["settings"] = {k: summary(v) for k, v in GROUPS.items()}
    cal = json.load(open(os.path.join(OUT, "T3_calibration.json"))) if os.path.exists(os.path.join(OUT, "T3_calibration.json")) else {}
    t3 = {"0.5 (Bhowmik production, not re-tuned)": dict(s=1.0, **(summary(REF) or {}))}
    for tc in TAUS:
        k = f"{tc:g}"
        if k in cal:
            prod = summary([f"T3_tc{k}_s{s}" for s in SEEDS])
            if prod:
                t3[k] = dict(s=cal[k]["s_prod"], calibration=cal[k], source="production", **prod)
            else:                       # batch stopped early: the two calibration runs (200 analysed yr each) stand in
                t3[k] = dict(s=cal[k]["s_cal"], calibration=cal[k], source="calibration runs", **(summary(cal[k]["runs"]) or {}))
    out["T3"] = t3
    t4 = {"reference (delta 0)": summary(REF)}
    for kind in KINDS:
        for d in DELTAS:
            t4[f"{kind} {int(round(d * 100))} %"] = summary([f"T4_{kind}_d{int(round(d * 100))}_s{s}" for s in SEEDS])
    out["T4"] = t4
    g = dict(GROUPS)
    g.update({f"T3 tau_c {tc:g}": [f"T3_tc{tc:g}_s{s}" for s in SEEDS] for tc in TAUS})
    g.update({f"T4 {k} {int(round(d * 100))} %": [f"T4_{k}_d{int(round(d * 100))}_s{s}" for s in SEEDS] for k in KINDS for d in DELTAS})
    out["detector_check"] = detector_check(g)
    json.dump(out, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)
    for sec in ("settings", "T4"):
        for k, v in out[sec].items():
            if v:
                print(f"{sec} {k}: runs {v['n_runs']}, cycles {v['n_cycles']}, lag-1 {v['lag1_A_pooled']:+.2f}, longest {v['mean_longest_per_13']:.1f}, "
                      f"rms A {v['rms_A']:.2f}, rms Apf {v['rms_Apf']:.2f}, Sun-like {v['frac_sunlike'] * 100:.0f} %, mean A {v['mean_A']:+.3f}, timing lag-1 {v['timing_lag1_pooled']:+.2f}")
    for k, v in out["T3"].items():
        if v.get("n_runs"):
            print(f"T3 tau_c {k}: s {v['s']:.3f}, rms A {v['rms_A']:.3f}, lag-1 {v['lag1_A_pooled']:+.2f}, longest {v['mean_longest_per_13']:.1f}, rms Apf {v['rms_Apf']:.2f}")

if __name__ == "__main__":
    {"validate": validate, "launch": launch, "analyse": analyse}[sys.argv[1]]()
