"""closure_v5.py -- v5 additions against Bhowmik (2019, A&A 632, A117).  New file; existing code is imported, not changed.

    python closure_v5.py n1|n4|n2|n3|n5      one item per call; results in out/v5/<item>.json, progress in out/v5/progress.md

N1  pure-lag control.  Shifting one hemisphere's emergence times is a re-evaluation of the exact Green's function
    (core.Model.run on a catalogue with shifted years); one direct time-stepped SFT run checks that shortcut, and the
    mirror check requires D(N shifted) = -D(S shifted, mirrored catalogue).
    a) symmetrised catalogue (core.symmetrised: identical hemispheres), S delayed 0.5-2.5 yr: change in the final polar
       field asymmetry (|B_N|-|B_S|)/|B_sym| and in the reversal-time difference, per cycle.
    b) real catalogue, N or S shifted by +-0.5, +-1.0 yr: per-cycle change in D as a fraction of rms D.
N4  parity switch on independent records: DPD (1974-) through the project's own annual()/x1() (same windows and
    latitude normalisation as RGO/USAF); Temmer (1945-1991) / SILSO (1992-) hemispheric sunspot numbers, which carry no
    latitudes, as plain yearly N-S sums over the same windows (RGO with plain sums is the control).
N2  bridge from our area-only SFT (r 0.61) to Bhowmik's data-driven 0.91, one change at a time (cumulative).
N3  her one-factor sensitivities in her metric (|B_N|-|B_S|)/|B_sym| on the symmetrised catalogue.
N5  her observational link recomputed on her pairs and definition.
"""
import os, sys, json, time
import numpy as np
import pandas as pd
from scipy.signal import fftconvolve

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from core import *                                                   # noqa: F401,F403  (CAT, Model, mirror, symmetrised, ...)
import core as K

OV = os.path.join(OUT, "v5"); os.makedirs(OV, exist_ok=True)
EDGE = (12, 24)                                                      # first/last readout cycles: shifted flux can leave the record
INNER = [n for n in CYC if n not in EDGE]

def save(name, obj):
    json.dump(obj, open(os.path.join(OV, name), "w"), indent=1, default=lambda v: v.tolist() if hasattr(v, "tolist") else float(v))

def progress(item, status, files):
    with open(os.path.join(OV, "progress.md"), "a", encoding="utf-8") as f:
        f.write(f"- {time.strftime('%Y-%m-%d %H:%M')} {item}: {status} -- {files}\n")

def shifted(cat, hemi, dt):
    c = cat.copy(); m = (c.latc.values > 0) if hemi == "N" else (c.latc.values < 0)
    y = c.year.values.copy(); y[m] = y[m] + dt; c["year"] = y
    return c

def dvals(res, cyc=CYC):
    return np.array([d_of(res["full"][n]) for n in cyc])

def absNS(res, cyc=CYC):
    return np.array([np.abs(np.asarray(res["full"][n], float)) for n in cyc])      # (ncyc, 2): |B_N|, |B_S|

# ---------------------------------------------------------------------------------------------------------- N1
def direct_run(m, cat):
    """time-stepped SFT with the sources injected at their emergence steps: the check of the Green's-function shortcut."""
    s = m.s; Ss = X.sources(cat, cat.alpha.values).astype(np.float64)
    ks = np.floor((cat.year.values - T0[12]) * SEC / s.dt).astype(int)
    ns = int(np.ceil((max(X.READ.values()) + 1.6 + 0.1 - T0[12]) * SEC / s.dt)) + 1
    inj = np.zeros((ns, Ss.shape[1])); ok = (ks >= 0) & (ks < ns)
    np.add.at(inj, ks[ok], Ss[ok])
    B = -np.sin(s.latr); vals = np.empty((ns, 2))
    for k in range(ns):
        B = B + inj[k]; vals[k] = B @ WCAP; B = s.step(B)
    tk = T0[12] + np.arange(ns) * s.dt / SEC
    return {n: vals[np.abs(tk - X.READ[n]) <= 1.6].mean(0) for n in CYC}, int((~ok).sum())

def cap_series(m, cat, ic_series):
    """N/S cap field at every step, by convolving the per-cell source series with the Green's function."""
    s = m.s; Ss = X.sources(cat, cat.alpha.values).astype(np.float64)
    ks = np.floor((cat.year.values - T0[12]) * SEC / s.dt).astype(int); ns = ic_series.shape[0]
    src = np.zeros((ns, Ss.shape[1])); ok = (ks >= 0) & (ks < ns); np.add.at(src, ks[ok], Ss[ok])
    out = ic_series.copy()
    for h in (0, 1):
        conv = fftconvolve(src, m.G[:, h, :].T, axes=0)[:ns]                   # (ns, nlat): per-cell contribution
        out[:, h] += conv.sum(1)
    return out

def ic_series(m):
    s = m.s; ns = int(np.ceil((max(X.READ.values()) + 1.6 + 0.1 - T0[12]) * SEC / s.dt)) + 1
    B = -np.sin(s.latr); vals = np.empty((ns, 2))
    for k in range(ns):
        vals[k] = B @ WCAP; B = s.step(B)
    return vals, T0[12] + np.arange(ns) * s.dt / SEC

def reversals(series, tk):
    """per cycle and hemisphere: zero crossing of the 1-yr running-mean cap field nearest the observed reversal t_R."""
    w = max(1, int(round(1.0 / (tk[1] - tk[0])))); ker = np.ones(w) / w; out = {}
    for h in (0, 1):
        sm = np.convolve(series[:, h], ker, mode="same"); z = np.flatnonzero(np.diff(np.sign(sm)) != 0)
        zt = tk[z]
        for n in CYC:
            c = zt[(zt > SCB[n]) & (zt < SCB[n + 1])]
            out.setdefault(n, [np.nan, np.nan])[h] = float(c[np.argmin(np.abs(c - TR[n]))]) if len(c) else np.nan
    return out

def n1():
    t0 = time.time(); m = Model(); res0 = m.run(CAT); D0 = dvals(res0); rmsD = float(np.sqrt(np.mean(D0 ** 2)))
    # --- shortcut check: one direct time-stepped SFT run (S shifted by +1 yr)
    cs = shifted(CAT, "S", 1.0); g = m.run(cs); dr, dropped = direct_run(m, cs)
    chk = float(max(np.max(np.abs(np.asarray(g["full"][n]) - dr[n])) / np.max(np.abs(dr[n])) for n in CYC))
    # --- mirror check: N shifted on the catalogue = minus (S shifted on the mirrored catalogue)
    dN = dvals(m.run(shifted(CAT, "N", 1.0))); dSm = dvals(m.run(shifted(mirror(CAT), "S", 1.0)))
    mir = float(np.max(np.abs(dN + dSm)) / rmsD)
    # --- b) real catalogue
    b = {}
    for hemi in ("N", "S"):
        for dt in (-1.0, -0.5, 0.5, 1.0):
            dd = dvals(m.run(shifted(CAT, hemi, dt))) - D0
            b[f"{hemi}{dt:+.1f}"] = dict(dD=dd.tolist(), frac=(np.abs(dd) / rmsD).tolist())
    f1 = np.array([[b[f"{h}{d:+.1f}"]["frac"][CYC.index(n)] for n in INNER] for h in "NS" for d in (-1.0, 1.0)])
    med1, max1 = float(np.median(f1)), float(f1.max())
    label = "proxy for late flux" if med1 < 0.2 else ("partly lag" if med1 < 0.5 else "lag matters in this SFT; differs from Bhowmik")
    # --- a) symmetrised catalogue, S delayed
    sym = symmetrised(CAT); r0 = m.run(sym); B0 = absNS(r0); Bsym = B0[:, 0]
    assert np.max(np.abs(B0[:, 0] - B0[:, 1]) / Bsym) < 1e-9                       # identical hemispheres
    ics, tk = ic_series(m); rev0 = reversals(cap_series(m, sym, ics), tk)
    a = {}
    for dt in (0.5, 1.0, 1.5, 2.0, 2.5):
        cc = shifted(sym, "S", dt); B = absNS(m.run(cc)); rv = reversals(cap_series(m, cc, ics), tk)
        asym = (B[:, 0] - B[:, 1]) / Bsym; dS = (B[:, 1] - Bsym) / Bsym; dN = (B[:, 0] - Bsym) / Bsym
        drev = [(rv[n][1] - rv[n][0]) - (rev0[n][1] - rev0[n][0]) for n in CYC]
        ii = [CYC.index(n) for n in INNER]
        a[f"{dt:.1f}"] = dict(asym=asym.tolist(), dN=dN.tolist(), dS=dS.tolist(), drev=drev,
                              asym_med_abs=float(np.median(np.abs(asym[ii]))), asym_p16_84=np.percentile(asym[ii], [16, 84]).tolist(),
                              drev_med=float(np.nanmedian(np.array(drev)[ii])), drev_p16_84=np.nanpercentile(np.array(drev)[ii], [16, 84]).tolist())
    out = dict(cycles=CYC, edge=list(EDGE), rmsD=rmsD, D0=D0.tolist(), shortcut_vs_direct_max_rel=chk, direct_dropped_sources=dropped,
               mirror_max_rel=mir, b=b, b_1yr_median_frac=med1, b_1yr_max_frac=max1, label=label, a=a,
               bhowmik=dict(lag_yr=[0.5, 2.5], final_flux_change="unaffected (her Sect. 4.1)", reversal_diff_at_2p5yr_approx=1.0,
                            r_lag_reversal=0.9507, source="Bhowmik 2019 A&A 632 A117, Sect. 4.1 and Fig. 3c (PDF in Reference Papers)"))
    save("N1.json", out)
    print(f"N1 shortcut vs direct SFT max rel diff {chk:.1e} (sources outside the record {dropped}); mirror {mir:.1e}")
    print(f"N1b 1-yr shift: median |dD|/rmsD {med1:.2f}, max {max1:.2f} (cycles 13-23) -> {label}")
    for dt in ("1.0", "2.5"):
        print(f"N1a S delay {dt} yr: median |asym| {a[dt]['asym_med_abs']*100:.1f} %, reversal diff {a[dt]['drev_med']:+.2f} yr")
    print(f"N1 done in {time.time() - t0:.0f} s")
    progress("N1", f"done; label '{label}'", "out/v5/N1.json")

def n1x():
    """explanation of N1: Bhowmik's SFT has no surface decay term; ours decays with tau = 5 yr.  The same N1 measures
    with the kernel rescaled exactly to tau = 10 yr and to no decay (tau = 1e6 yr)."""
    base = Model(); out = {}
    for tau in (5.0, 10.0, 1e6):
        m = base if tau == 5.0 else Model(tau=tau, base=base)
        D0 = dvals(m.run(CAT)); rmsD = float(np.sqrt(np.mean(D0 ** 2)))
        f = [np.abs(dvals(m.run(shifted(CAT, h, d))) - D0)[[CYC.index(n) for n in INNER]] / rmsD for h in "NS" for d in (-1.0, 1.0)]
        sym = symmetrised(CAT); Bs = absNS(m.run(sym))[:, 0]; ii = [CYC.index(n) for n in INNER]
        a = {}
        for dt in (1.0, 2.5):
            B = absNS(m.run(shifted(sym, "S", dt))); a[f"{dt:.1f}"] = float(np.median(np.abs((B[:, 0] - B[:, 1]) / Bs)[ii]))
        out[f"{tau:g}"] = dict(b_1yr_median_frac=float(np.median(f)), a_asym_med_abs=a)
        print(f"N1x tau {tau:g} yr: 1-yr shift median |dD|/rmsD {np.median(f):.2f}; S delay 1 / 2.5 yr median |asym| {a['1.0']*100:.1f} / {a['2.5']*100:.1f} %")
    save("N1x_decay.json", out)
    progress("N1x", "done (decay-time explanation)", "out/v5/N1x_decay.json")

# ---------------------------------------------------------------------------------------------------------- N4
def win_years(n, s, e):
    return [y for y in range(1870, 2030) if TR[n] + s <= y + 0.5 < TR[n] + e + 1]

def plain_x(yearly, n, s, e):
    """sum of yearly N-S over the same whole-year window as x1; None when a year is missing from the record."""
    ys = win_years(n, s, e)
    return float(sum(yearly[y] for y in ys)) if all(y in yearly for y in ys) else None

def n4():
    q2b = load("q2b_phase_sign.json")["switch_max_vs_decline"]; pw = load("q2_window.json")["physics_window"]
    PS, PE = pw["s"], pw["e"]; assert (PS, PE) == (3, 4)
    recs = {}
    # RGO/USAF, x1 (reference, as in q2b) and plain sums (control for the sunspot-number records)
    rgo_x = {n: (x1(n, 0, 1), x1(n, PS, PE)) for n in range(12, 25)}
    yr = {int(r.yr): float(r.AN - r.AS) for r in ANN.itertuples()}
    rgo_p = {n: (plain_x(yr, n, 0, 1), plain_x(yr, n, PS, PE)) for n in range(12, 25)}
    # DPD through the project's annual()/x1(): same windows and latitude normalisation
    dp = pd.read_pickle(os.path.join(OUT, "dpd_groups.pkl"))
    full_years = [y for y in range(1974, 2018) if ((dp.dec >= y) & (dp.dec < y + 1)).any() and y + 1 <= dp.dec.max()]
    dd = pd.DataFrame(dict(yr=np.floor(dp.dec.values).astype(int), dec=dp.dec.values, hemi=np.where(dp.north.values, "N", "S"),
                           area=dp.area.values, abslat=np.abs(dp.lat.values), lat=dp.lat.values, gid=dp.gid.values,
                           cyc=np.searchsorted([SCB[n] for n in range(12, 27)], dp.dec.values, side="right") + 11,
                           rmax=dp.rmax.values, rlat=dp.rlat.values))
    ann_d = A.annual(dd)
    cov = lambda n: all(y in full_years for y in win_years(n, 0, 1) + win_years(n, PS, PE))
    dpd_x = {n: (x1(n, 0, 1, ann=ann_d), x1(n, PS, PE, ann=ann_d)) for n in range(12, 25) if cov(n)}
    # Temmer et al. 2006 (1945-1991) then SILSO (1992-) hemispheric sunspot numbers, full years only
    EXT = os.path.join(P.REV, "external", "silso_hem")
    CBN = ["date", "n_comb", "s_comb", "n_comb_sm", "s_comb_sm", "n_area", "s_area", "n_area_sm", "s_area_sm",
           "n_t06", "s_t06", "n_t06_sm", "s_t06_sm", "n_silso", "s_silso", "n_silso_sm", "s_silso_sm"]
    cb = pd.read_csv(os.path.join(EXT, "Catalogue_B.txt"), sep=r"\s+", header=None, names=CBN)
    cb["y"] = cb.date.str[:4].astype(int)
    t06 = cb[(cb.y <= 1991) & (cb.n_t06 >= 0) & (cb.s_t06 >= 0)]; sil = cb[(cb.y >= 1992) & (cb.n_silso >= 0) & (cb.s_silso >= 0)]
    yt = {}
    for sub, a, b in ((t06, "n_t06", "s_t06"), (sil, "n_silso", "s_silso")):
        for y, g in sub.groupby("y"):
            if len(g) == 12:
                yt[int(y)] = float((g[a] - g[b]).sum())
    tem_p = {n: v for n in range(12, 25) if None not in (v := (plain_x(yt, n, 0, 1), plain_x(yt, n, PS, PE)))}
    def table(xs, ref_switch=None):
        rows = {}
        for n, (m_, d_) in xs.items():
            if m_ is None or d_ is None:
                continue
            rows[n] = dict(max_sign=int(np.sign(m_)), decline_sign=int(np.sign(d_)), switch=bool(np.sign(m_) != np.sign(d_)),
                           decline_vs_Apf=(bool(np.sign(d_) == np.sign(APF[n])) if n in APF else None))
        ev = [n for n in rows if n % 2 == 0]; od = [n for n in rows if n % 2]
        ap = [n for n in rows if rows[n]["decline_vs_Apf"] is not None]
        return dict(cycles=rows, even=f"{sum(rows[n]['switch'] for n in ev)}/{len(ev)}", odd=f"{sum(rows[n]['switch'] for n in od)}/{len(od)}",
                    decline_sign_vs_Apf=f"{sum(rows[n]['decline_vs_Apf'] for n in ap)}/{len(ap)}")
    recs = {"RGO/USAF (x1)": table(rgo_x), "RGO/USAF (plain sums)": table(rgo_p), "DPD (x1)": table(dpd_x), "Temmer/SILSO (plain sums)": table(tem_p)}
    ref = recs["RGO/USAF (x1)"]; assert (ref["even"], ref["odd"]) == (q2b["even"], q2b["odd"]) == ("5/7", "0/6")
    assert all(ref["cycles"][n]["switch"] == q2b["per_cycle"][str(n)] for n in ref["cycles"])
    # cycle-by-cycle agreement with RGO on the switch, where both records cover the cycle
    for k, v in recs.items():
        common = [n for n in v["cycles"] if n in ref["cycles"]]
        v["switch_agrees_with_RGO"] = f"{sum(v['cycles'][n]['switch'] == ref['cycles'][n]['switch'] for n in common)}/{len(common)}"
    save("N4.json", dict(windows=dict(max=[0, 1], decline=[PS, PE]), dpd_full_years=[min(full_years), max(full_years)],
                         temmer_silso_years=[min(yt), max(yt)], records=recs))
    for k, v in recs.items():
        print(f"N4 {k}: cycles {min(v['cycles'])}-{max(v['cycles'])}; switch even {v['even']}, odd {v['odd']}; "
              f"decline sign = A_pf sign {v['decline_sign_vs_Apf']}; switch agrees with RGO {v['switch_agrees_with_RGO']}")
    progress("N4", "done", "out/v5/N4.json")

# ---------------------------------------------------------------------------------------------------------- N2
def her_T():
    """per-cycle sqrt-law tilt coefficient T_sqr from the project's tilt workbook (No_Filter, Full phase), averaged over
    the datasets that cover the cycle (MWO, KK 15-21; DPD 21-24); other cycles take the mean.  Stands in for Jiang et al.
    (2011) T_n, whose values are not in the project (Bhowmik cites their Eq. 15 without listing them)."""
    tw = pd.read_excel(os.path.join(HERE, "..", "Datasets Used", "tilt_complete_analysis.xlsx"), sheet_name="Phase_Summary")
    tw = tw[(tw.Filter == "No_Filter") & (tw.Scope == "Full")]
    T = tw.groupby("Cycle").T_sqr.mean().to_dict(); mean = float(np.mean(list(T.values())))
    return {int(k): float(v) for k, v in T.items()}, mean

def n2():
    m = Model(); obsA = A.polar_signed("adopted"); obsM = A.polar_signed("homogeneous")          # MWO <= 19 then WSO; MWO only
    f = sum(np.abs(obsM[n]).sum() for n in (20, 21, 22)) / sum(np.abs(obsA[n]).sum() for n in (20, 21, 22))
    harm = {n: np.asarray(obsA[n], float) * (f if n >= 20 else 1.0) for n in C1}                  # WSO put on the MWO scale
    def stats(res, cyc, obs, kind):
        mod = {n: np.asarray(res["full"][n], float) for n in cyc}
        if kind == "norm":
            x = [apf_of(mod[n]) for n in cyc]; y = [APF[n] for n in cyc] if obs is None else [apf_of(obs[n]) for n in cyc]
        else:
            x = [abs(mod[n][0]) - abs(mod[n][1]) for n in cyc]; y = [abs(obs[n][0]) - abs(obs[n][1]) for n in cyc]
        od = [i for i, n in enumerate(cyc) if n % 2]; ev = [i for i, n in enumerate(cyc) if n % 2 == 0]
        sub = lambda ii: r(np.array(x)[ii], np.array(y)[ii]) if len(ii) > 2 else float("nan")
        return dict(r=r(x, y), n=len(cyc), r_odd=sub(od), n_odd=len(od), r_even=sub(ev), n_even=len(ev))
    def flux(res, cyc, obs):
        x = np.ravel([np.abs(np.asarray(res["full"][n], float)) for n in cyc]); y = np.ravel([np.abs(obs[n]) for n in cyc])
        return dict(r=r(x, y), n=len(x))
    T, Tm = her_T()
    tilt = 0.7 * np.array([T.get(int(c), Tm) for c in CAT.cycle.values]) * np.sqrt(np.abs(CAT.latc.values)) * np.sign(CAT.alpha.values)
    unit = float(np.median(np.abs(tilt)) / np.median(np.abs(CAT.alpha.values)))
    base, hers = m.run(CAT), m.run(CAT, tilt=tilt)
    C17 = list(range(17, 25)); C17m = list(range(17, 23)); C17x = [n for n in C17m if n != 18]
    steps = [("0 our area-only run, A_pf, SC15-24", stats(base, C1, None, "norm"), flux(base, C1, harm)),
             ("1 cycles 17-24", stats(base, C17, None, "norm"), flux(base, C17, harm)),
             ("2 unnormalised |B_N|-|B_S|", stats(base, C17, harm, "raw"), flux(base, C17, harm)),
             ("3 her tilt law 0.7 T sqrt|lat|", stats(hers, C17, harm, "raw"), flux(hers, C17, harm)),
             ("4 MWO only (minima 17-22; 23-24 not in the project's MWO file)", stats(hers, C17m, obsM, "raw"), flux(hers, C17m, obsM)),
             ("5 drop the 18/19 minimum (polar flux only)", None, flux(hers, C17x, obsM))]
    s0 = steps[0][1]; assert abs(s0["r"] - 0.61) < 0.005 and abs(s0["r_odd"] - 0.97) < 0.005 and abs(s0["r_even"] - 0.01) < 0.005
    rows = [dict(step=a, asym=b, flux=c) for a, b, c in steps]
    da = [(rows[i]["step"], rows[i]["asym"]["r"] - rows[i - 1]["asym"]["r"]) for i in range(1, 5)]
    dfl = [(rows[i]["step"], rows[i]["flux"]["r"] - rows[i - 1]["flux"]["r"]) for i in range(1, 6)]
    big_a = max(da, key=lambda t: abs(t[1])); big_f = max(dfl, key=lambda t: abs(t[1]))
    reach = any((rw["asym"] or {}).get("r", 0) >= 0.85 for rw in rows) or any(rw["flux"]["r"] >= 0.85 for rw in rows)
    save("N2.json", dict(steps=rows, wso_to_mwo_factor=f, her_T=T, her_T_mean=Tm, tilt_median_ratio_hers_to_ours=unit,
                         largest_step_asym=big_a, largest_step_flux=big_f, reaches_0p9=bool(reach),
                         bhowmik=dict(asym_r=0.91, asym_n=8, flux_r=0.84, flux_r_excl_18_19=0.94, cycles="17-24")))
    for rw in rows:
        a = rw["asym"]; fl = rw["flux"]
        print(f"N2 {rw['step'][:34]:34s} asym " + (f"r {a['r']:+.2f} n {a['n']} (odd {a['r_odd']:+.2f}, even {a['r_even']:+.2f})" if a else "--")
              + f" | flux r {fl['r']:+.2f} n {fl['n']}")
    print(f"N2 largest move: asym {big_a[0][:2]} {big_a[1]:+.2f}; flux {big_f[0][:2]} {big_f[1]:+.2f}; reaches ~0.9: {reach}; "
          f"WSO->MWO factor {f:.4f}; her/our median tilt {unit:.2f}")
    progress("N2", "done", "out/v5/N2.json")

# ---------------------------------------------------------------------------------------------------------- N3
def boostS(cat, f):
    """southern rows' flux x f (separation x sqrt f, area x f), as core.scaled does for a whole catalogue."""
    c = cat.copy(); s = c.latc.values < 0; g = np.where(s, f, 1.0)
    c["flux"] = c.flux.values * g; c["ddeg"] = c.ddeg.values * np.sqrt(g); c["area"] = c.area.values * g
    return c

def n3(nreal=50):
    base = Model(); sym = symmetrised(CAT); B0 = absNS(base.run(sym)); Bs = B0[:, 0]; ii = [CYC.index(n) for n in INNER]
    def metric(res):
        B = absNS(res); return (B[:, 0] - B[:, 1]) / Bs
    def summ(v):
        v = np.abs(np.asarray(v)); return dict(median=float(np.median(v)), p84=float(np.percentile(v, 84)), max=float(v.max()))
    out = {}
    m70 = metric(base.run(boostS(sym, 1.7)))[ii]
    mir = float(np.max(np.abs(metric(base.run(boostS(sym, 1.7))) + metric(base.run(mirror(boostS(mirror(sym), 1.7)))))))
    out["+70% S flux"] = summ(m70)
    for lab, v0s in (("+50% S peak flow", 22.5), ("-50% S peak flow", 7.5), ("-10% S peak flow", 13.5)):
        mm = Model(v0S=v0s); res = mm.run(sym); out[lab] = summ(metric(res)[ii])
        if v0s == 13.5:
            ap = np.array([apf_of(res["full"][n]) for n in INNER]); out[lab]["A_pf_median_abs"] = float(np.median(np.abs(ap)))
    ts = np.array([metric(base.run(sym, tilt=tilt_draw(sym, seed)))[ii] for seed in range(1, nreal + 1)])    # (nreal, ncyc)
    out["tilt scatter (Jiang 2014), 50 realisations"] = dict(summ(ts), max_per_realisation_median=float(np.median(np.abs(ts).max(1))))
    her = {"+70% S flux": 0.14, "+50% S peak flow": 0.03, "-50% S peak flow": 0.03, "tilt scatter (Jiang 2014), 50 realisations": 0.36}
    save("N3.json", dict(metric="(|B_N|-|B_S|)/|B_sym| at each minimum, cycles 13-23, symmetrised catalogue", ours=out, hers=her,
                         mirror_max=mir, note_flux="her +70% run also lowers the tilt with flux (T = 1.71 - 0.12 Phi); ours scales flux only"))
    for k, v in out.items():
        print(f"N3 {k:44s} median {v['median']*100:5.1f} %, p84 {v['p84']*100:5.1f} %, max {v['max']*100:5.1f} %  (hers {her.get(k, float('nan'))*100:.0f} %)"
              + (f"; A_pf {v['A_pf_median_abs']:.3f}" if "A_pf_median_abs" in v else ""))
    print(f"N3 mirror check max {mir:.1e}")
    progress("N3", "done", "out/v5/N3.json")

def n3b():
    """REPORT's '10 % flow difference -> offset 0.05' is v0N = 15 x 1.1, v0S = 15 x 0.9 (q1_sources TCFG 'v0+0.1');
    restated in Bhowmik's metric on the same symmetrised catalogue."""
    base = Model(); sym = symmetrised(CAT); Bs = absNS(base.run(sym))[:, 0]; ii = [CYC.index(n) for n in INNER]
    res = Model(v0N=16.5, v0S=13.5).run(sym); B = absNS(res)
    met = ((B[:, 0] - B[:, 1]) / Bs)[ii]; ap = np.array([apf_of(res["full"][n]) for n in INNER])
    q1 = load("q1_sources.json")["c_transport_only"]["v0+0.1"]["apf"]
    o = dict(metric_median=float(np.median(met)), metric_max=float(np.max(np.abs(met))), A_pf_median=float(np.median(ap)),
             A_pf_median_q1_sources=float(np.median(list(q1.values()))))
    d = json.load(open(os.path.join(OV, "N3.json"))); d["ours"]["REPORT 10% flow difference (v0N 16.5, v0S 13.5)"] = o; save("N3.json", d)
    print(f"N3b REPORT setting: A_pf median {o['A_pf_median']:+.3f} (q1_sources {o['A_pf_median_q1_sources']:+.3f}); "
          f"her metric median {o['metric_median']*100:+.1f} %, max {o['metric_max']*100:.1f} %")
    progress("N3b", "done (REPORT's 10% flow difference restated)", "out/v5/N3.json")

# ---------------------------------------------------------------------------------------------------------- N5
def n5():
    """Bhowmik's observational link on her pairs and definition: | |Phi_N| - |Phi_S| | at the minimum of cycle n (MWO,
    n = 15-22) against | A_N - A_S | of the total sunspot area of cycle n+1 (unnormalised).  Pearson r only (her
    Spearman 0.81 is not recomputed: rank tests are outside this project's plain statistics)."""
    obsM = A.polar_signed("homogeneous")                                   # MWO polar faculae, minima ending 15..22
    tot = {}
    for n in range(15, 25):
        mm = (D.dec.values >= SCB[n]) & (D.dec.values < SCB[n + 1]); h = D.hemi.values
        tot[n] = (float(D.area.values[mm & (h == "N")].sum()), float(D.area.values[mm & (h == "S")].sum()))
    pr = list(range(15, 23))
    xs = [abs(obsM[n][0]) - abs(obsM[n][1]) for n in pr]; ys = [tot[n + 1][0] - tot[n + 1][1] for n in pr]
    her_abs = dict(r=r(np.abs(xs), np.abs(ys)), n=len(pr)); her_signed = dict(r=r(xs, ys), n=len(pr), signs=int(signs(xs, ys)))
    ours = list(range(15, 24))
    our_abs = dict(r=r([abs(APF[n]) for n in ours], [abs(AN[n + 1]) for n in ours]), n=len(ours))
    our_signed = dict(r=r([APF[n] for n in ours], [AN[n + 1] for n in ours]), n=len(ours))
    assert abs(our_abs["r"] - 0.62) < 0.005 and abs(our_signed["r"] - 0.68) < 0.005               # REPORT Q0 / Q3
    # same normalised definition restricted to her pairs (n = 15-22)
    our_her_pairs = dict(r=r([abs(APF[n]) for n in pr], [abs(AN[n + 1]) for n in pr]), n=len(pr))
    out = dict(her_definition_abs=her_abs, her_definition_signed=her_signed, ours_abs=our_abs, ours_signed=our_signed,
               ours_abs_on_her_pairs=our_her_pairs, pairs=[[n, n + 1] for n in pr], data="MWO polar faculae (homogeneous), RGO/USAF areas",
               bhowmik=dict(pearson=0.73, spearman_text=0.81, fig2_caption="Spearman given as 0.73 in the caption (flag)"))
    save("N5.json", out)
    print(f"N5 her definition (|dPhi| vs |dA_total|, n = 15-22, MWO): r {her_abs['r']:+.2f} (n {her_abs['n']}); signed r {her_signed['r']:+.2f}, "
          f"signs {her_signed['signs']}/{her_signed['n']}")
    print(f"N5 ours (|A_pf| vs |A(n+1)|): r {our_abs['r']:+.2f} (n {our_abs['n']}); on her pairs {our_her_pairs['r']:+.2f} (n {our_her_pairs['n']}); signed {our_signed['r']:+.2f}")
    progress("N5", "computed (table in REPORT_v5_additions.md)", "out/v5/N5.json")

if __name__ == "__main__":
    globals()[sys.argv[1]]()
