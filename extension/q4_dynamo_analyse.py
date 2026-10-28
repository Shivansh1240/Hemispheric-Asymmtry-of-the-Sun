"""q4_dynamo_analyse.py -- cycle-by-cycle statistics of a self-excited dynamo run, computed exactly as for the Sun.

Per cycle n (minimum to minimum of the 1-yr running mean of the total BL emergence):
  A(n)       (E_N - E_S)/(E_N + E_S), E = BL flux deposited (the model's emergence)
  A_pf(n)    (|P_N|-|P_S|)/(|P_N|+|P_S|), P = cap-mean surface B_r over |lat| >= 70, +-1.6 yr about the minimum
             that ends n; S(n) its signed form; sgnBN its polarity
  Xmax(n)    E_N - E_S summed over [t_max, t_max + 2) yr (t_max = peak of the smoothed emergence; the t_R analogue)
  Xdec(n)    E_N - E_S summed over [t_max + 3, t_max + 5) yr (the declining-phase window of Q2)
Statistics are the ones used for the Sun: Pearson r with n, signs, run lengths (also in 13-cycle segments).
"""
import os, sys, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from q4_dynamo import cycles, OUT                               # noqa: E402

def r(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[m], y[m])[0, 1]) if m.sum() > 2 else float("nan")

def runs(seq):
    return [len(list(g)) for _, g in itertools.groupby(seq)]

def per_cycle(res, spin=100.0):
    t = res["t"]; EN, ES, PN, PS = res["EN"], res["ES"], res["PN"], res["PS"]
    mins, Es = cycles(res)
    mins = [i for i in mins if t[i] >= spin]
    rows = []
    for a, b in zip(mins[:-1], mins[1:]):
        sl = slice(a, b); eN, eS = EN[sl].sum(), ES[sl].sum()
        imax = a + int(np.argmax(Es[sl])); tmax = t[imax]
        def win(t0, t1):
            m = (t >= t0) & (t < t1); return float(EN[m].sum() - ES[m].sum())
        m = np.abs(t - t[b]) <= 1.6
        pN, pS = PN[m].mean(), PS[m].mean()
        tN = float(np.average(t[sl], weights=EN[sl])); tS = float(np.average(t[sl], weights=ES[sl]))   # data's E definition
        # polar reversal of each hemisphere inside the cycle (1-yr running mean); offset S - N
        mm = (t >= t[a]) & (t < t[b]); kk = np.ones(20) / 20
        def zc(P):
            Ps = np.convolve(P, kk, mode="same")[mm]; z = np.flatnonzero(np.diff(np.sign(Ps)) != 0)
            return t[mm][z[0]] if len(z) else np.nan
        rev_off = zc(PS) - zc(PN)
        rows.append(dict(t0=t[a], t1=t[b], tmax=tmax, A=(eN - eS) / (eN + eS), amp=float(Es[sl].max()),
                         Apf=(abs(pN) - abs(pS)) / (abs(pN) + abs(pS)), S=(pN + pS) / (abs(pN) + abs(pS)),
                         sgnBN=float(np.sign(pN)), Xmax=win(tmax, tmax + 2), Xdec=win(tmax + 3, tmax + 5),
                         Xm14=win(tmax - 1, tmax + 5), lag=tN - tS, rev_off=float(rev_off)))
    return rows

def stats(rows):
    A = np.array([c["A"] for c in rows]); Apf = np.array([c["Apf"] for c in rows]); S = np.array([c["S"] for c in rows])
    par = -np.array([c["sgnBN"] for c in rows])                     # the data's (-1)^(n+1) = -sgn B_N
    Xm = np.array([c["Xmax"] for c in rows]); Xd = np.array([c["Xdec"] for c in rows]); X14 = np.array([c["Xm14"] for c in rows])
    n = len(A); seq = list(np.sign(A))
    seg = [runs(seq[i:i + 13]) for i in range(0, n - 12)]
    sw = np.sign(Xm) != np.sign(Xd)
    lag = np.array([c["lag"] for c in rows])
    out = dict(n_cycles=n, period=float(np.median([c["t1"] - c["t0"] for c in rows])),
               rms_lag=float(np.sqrt(np.mean(lag ** 2))), r_lag_A=r(lag, A), r_lag_nextA=r(lag[:-1], A[1:]),
               rms_A=float(np.sqrt(np.mean(A ** 2))), rms_Apf=float(np.sqrt(np.mean(Apf ** 2))),
               lag1_A=r(A[:-1], A[1:]), lag2_A=r(A[:-2], A[2:]), lag1_Apf=r(Apf[:-1], Apf[1:]),
               runs_mean_per_13=float(np.mean([len(s) for s in seg])) if seg else None,
               longest_run_mean_per_13=float(np.mean([max(s) for s in seg])) if seg else None,
               longest_run_all=int(max(runs(seq))),
               r_Apf_next_A=r(Apf[:-1], A[1:]), r_absApf_absnextA=r(np.abs(Apf[:-1]), np.abs(A[1:])),
               r_A_Apf_same=r(A, Apf),
               r_Xmax_Apf_unsigned=r(Xm, Apf), r_Xmax_Apf_paritysigned=r(par * Xm, Apf), r_Xmax_S=r(Xm, S),
               r_Xdec_Apf_unsigned=r(Xd, Apf), signs_Xdec_Apf=f"{int((np.sign(Xd) == np.sign(Apf)).sum())}/{n}",
               signs_Xmax_parity_Apf=f"{int((np.sign(par * Xm) == np.sign(Apf)).sum())}/{n}",
               r_X14_nextA=r(X14[:-1], A[1:]),
               switch_frac_par_plus=float(sw[par > 0].mean()) if (par > 0).any() else None,
               switch_frac_par_minus=float(sw[par < 0].mean()) if (par < 0).any() else None,
               polarity_alternates=float(np.mean(par[1:] != par[:-1])))
    return out

def sunlike(rows, max_off=2.0):
    """cycles whose two poles both reverse inside the cycle, within max_off yr of each other, and whose polarity
    alternates with the previous cycle -- the Sun-like part of a noisy run (consecutive pairs only kept by stats via
    the lag terms on the filtered list, so lags here mix non-adjacent cycles; quoted separately)."""
    keep = []
    for i, c in enumerate(rows):
        ok = np.isfinite(c["rev_off"]) and abs(c["rev_off"]) <= max_off and (i == 0 or c["sgnBN"] != rows[i - 1]["sgnBN"])
        keep.append(bool(ok))
    return keep

def stats_sunlike(rows):
    """lag statistics restricted to ADJACENT pairs of Sun-like cycles."""
    k = sunlike(rows); A = np.array([c["A"] for c in rows]); Apf = np.array([c["Apf"] for c in rows])
    pair = np.array([k[i] and k[i + 1] for i in range(len(rows) - 1)])
    sub = [c for c, kk in zip(rows, k) if kk]
    return dict(n_sunlike=int(sum(k)), frac_sunlike=float(np.mean(k)), n_adjacent_pairs=int(pair.sum()),
                rms_Apf=float(np.sqrt(np.mean(np.array([c["Apf"] for c in sub]) ** 2))) if sub else None,
                lag1_A=r(A[:-1][pair], A[1:][pair]), r_absApf_absnextA=r(np.abs(Apf[:-1][pair]), np.abs(A[1:][pair])),
                r_Apf_nextA=r(Apf[:-1][pair], A[1:][pair]),
                r_Xmax_Apf_unsigned=r([c["Xmax"] for c in sub], [c["Apf"] for c in sub]),
                r_Xdec_Apf_unsigned=r([c["Xdec"] for c in sub], [c["Apf"] for c in sub]))

if __name__ == "__main__":
    for name in sys.argv[1:]:
        res = dict(np.load(os.path.join(OUT, f"{name}.npz"), allow_pickle=True))
        rows = per_cycle(res); st = stats(rows); sl = stats_sunlike(rows)
        json.dump(dict(stats=st, sunlike=sl, cycles=rows, meta=json.loads(str(res["meta"]))), open(os.path.join(OUT, f"{name}_stats.json"), "w"), indent=1, default=float)
        print(name, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in st.items()}))
        print(name, "SUN-LIKE", json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in sl.items()}))
