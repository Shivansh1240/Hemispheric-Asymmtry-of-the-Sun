"""surya_runs.py -- runs of surya/surya_bhowmik.f (Bhowmik 2019's dynamo: Passos et al. 2014 on the public SURYA, plus
pumping) and their analysis with the same per-cycle statistics as the Q4 dynamo (q4_dynamo_analyse.py).  New file.

    python surya_runs.py launch      validation (noise-free, mirror pair) + 6-seed ensemble, all in parallel
    python surya_runs.py analyse     out/v7_surya/summary.json
"""
import os, sys, json, shutil, subprocess, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
EXE = os.path.join(HERE, "surya", "surya_bhowmik.exe"); INIT = os.path.join(HERE, "..", "master_code", "init.dat")
OUT = os.path.join(HERE, "out", "v7_surya"); os.makedirs(OUT, exist_ok=True)
YR = 3600.0 * 24 * 365 / 1e8                                   # code time unit is 1e8 s; the code's own years are 365 d
BHOWMIK = dict(f0r=3.6, flvBL=1.0, flvMF=2.0)                     # Passos Fig. 8 noise levels, Bhowmik pumping
RUNS = dict(det=dict(years=150, seed=1, f0r=3.6, flvBL=0.0, flvMF=0.0, iswap=0),
            mir0=dict(years=20, seed=7, iswap=0, **BHOWMIK), mir1=dict(years=20, seed=7, iswap=1, **BHOWMIK),
            **{f"ens_s{s}": dict(years=600, seed=s, iswap=0, **BHOWMIK) for s in range(1, 7)})

def prepare(name, c):
    d = os.path.join(OUT, name); os.makedirs(d, exist_ok=True)
    shutil.copy(INIT, os.path.join(d, "init.dat"))
    with open(os.path.join(d, "run.in"), "w") as f:
        f.write(f"{c['years'] * YR:.6f} {c['seed']} {c['f0r']} {c['flvBL']} {c['flvMF']} {c['iswap']}\n")
    json.dump(c, open(os.path.join(d, "config.json"), "w"))
    return d

def launch():
    procs = {}
    for name, c in RUNS.items():
        d = prepare(name, c)
        if os.path.exists(os.path.join(d, "final.dat")):
            continue
        procs[name] = subprocess.Popen([EXE], cwd=d, stdout=open(os.path.join(d, "stdout.txt"), "w"), stderr=subprocess.STDOUT)
    print(f"{time.strftime('%H:%M:%S')} started {list(procs)}", flush=True)
    while procs:
        time.sleep(30)
        for name in [n for n, p in procs.items() if p.poll() is not None]:
            print(f"{time.strftime('%H:%M:%S')} finished {name}: exit {procs.pop(name).returncode}", flush=True)
    print("all runs finished", flush=True)

def load(name):
    a = np.loadtxt(os.path.join(OUT, name, "polar.dat"))
    return dict(t=a[:, 0], PN=a[:, 1], PS=a[:, 2], EN=a[:, 3], ES=a[:, 4], BN=a[:, 5], BS=a[:, 6])

def analyse(spin=100.0):
    from q4_dynamo import cycles
    from q4_dynamo_analyse import per_cycle, stats, stats_sunlike, r as rr
    out = {}
    det = load("det"); mins, _ = cycles(det); tm = det["t"][mins]; tm = tm[tm > 30]
    late = det["t"] > 30; E = det["EN"] + det["ES"]
    out["deterministic"] = dict(period_yr=float(np.median(np.diff(tm))) if len(tm) > 2 else None, n_minima=int(len(tm)),
                                erupting_fraction=float(np.mean(E[late] > 0)), base_field_max_1e5G=float(det["BN"][late].max()),
                                NS_max_rel_diff=float(np.max(np.abs(det["PN"] + det["PS"])) / np.max(np.abs(det["PN"]))))
    a, b = load("mir0"), load("mir1")
    out["mirror"] = dict(PN_vs_swapped_PS=float(np.max(np.abs(a["PN"] + b["PS"])) / np.max(np.abs(a["PN"]))),
                         EN_vs_swapped_ES=float(np.max(np.abs(a["EN"] - b["ES"])) / max(np.max(a["EN"]), 1e-30)))
    A0, A1, allA, allApf, per, sl, lg0, lg1 = [], [], [], [], [], [], [], []
    for s in range(1, 7):
        name = f"ens_s{s}"
        if not os.path.exists(os.path.join(OUT, name, "final.dat")):
            continue
        rows = per_cycle(load(name), spin=spin); st = stats(rows)
        A = np.array([c["A"] for c in rows]); Apf = np.array([c["Apf"] for c in rows]); lg = np.array([c["lag"] for c in rows])
        A0 += list(A[:-1]); A1 += list(A[1:]); allA += list(A); allApf += list(Apf); lg0 += list(lg[:-1]); lg1 += list(lg[1:])
        sl.append(stats_sunlike(rows)["frac_sunlike"])
        per.append(dict(seed=s, n_cycles=len(A), period=st["period"], lag1=st["lag1_A"], longest=st["longest_run_mean_per_13"],
                        rms_A=st["rms_A"], rms_Apf=st["rms_Apf"], r_Apf_next_A=st["r_Apf_next_A"],
                        switch_plus=st["switch_frac_par_plus"], switch_minus=st["switch_frac_par_minus"],
                        frac_sunlike=sl[-1], A=A.tolist(), Apf=Apf.tolist()))
    if per:
        allA = np.array(allA); allApf = np.array(allApf)
        out["ensemble"] = dict(per_seed=per, n_cycles=len(allA), n_pairs=len(A0), lag1_A_pooled=rr(A0, A1),
                               lag1_seed_range=[min(p["lag1"] for p in per), max(p["lag1"] for p in per)],
                               longest_median=float(np.median([p["longest"] for p in per])), rms_A=float(np.sqrt(np.mean(allA ** 2))),
                               rms_Apf=float(np.sqrt(np.mean(allApf ** 2))), frac_sunlike=float(np.median(sl)),
                               period=float(np.median([p["period"] for p in per])), lag1_timing_lag=rr(lg0, lg1))
    json.dump(out, open(os.path.join(OUT, "summary.json"), "w"), indent=1, default=float)
    print(json.dumps(out, default=float)[:1500])

if __name__ == "__main__":
    {"launch": launch, "analyse": analyse}[sys.argv[1]]()
