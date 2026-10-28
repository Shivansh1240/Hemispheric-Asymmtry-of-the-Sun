"""q3_persistence.py -- Q3: persistence and flips.

1. Run lengths of the dominant hemisphere (whole-cycle area sign), against what independent cycles would give.
2. n vs n+2 seeding: lag-1 and lag-2 correlations, plain r with n.
3. Flip and near-flip cycles: asymmetry of low-latitude large regions (>= 500 MSH, |lat| <= 20 at peak) and of
   all other regions; the SFT readout for SC13/14 and SC20/21 with the low-latitude large regions removed.
4. DATA CHECK of the SC21 flip on records independent of the RGO -> USAF/NOAA splice: Debrecen DPD (1974-2013, own
   plates, corrected whole-spot area) and the Temmer et al. (2006) hemispheric sunspot numbers (1945-1991).
5. Surface decay time tau = 2-17 yr: how the model's cycle-to-cycle memory changes (exact kernel rescaling).
"""
import itertools
import numpy as np
import pandas as pd
from core import *

out = {}
cyc13 = list(range(12, 25)); sgn = {n: int(np.sign(AN[n])) for n in range(12, 26)}

def runs(seq):
    return [len(list(g)) for _, g in itertools.groupby(seq)]

tic("run lengths")
seqs = {"area_all_SC12-24": [sgn[n] for n in cyc13], "area_all_SC12-25prov": [sgn[n] for n in range(12, 26)]}
d = D
for lab, lo, hi in (("large>=500", 500, np.inf), ("small<100", 0, 100), ("medium100-499", 100, 500)):
    seqs[f"area_{lab}"] = [int(np.sign(cycA(n, band_mask(d, lo, hi)))) for n in cyc13]
seqs["lowlat_large>=500_|lat|<=20"] = [int(np.sign(cycA(n, band_mask(d, 500, np.inf, 0, 20)))) for n in cyc13]
seqs["Apf_SC15-24"] = [int(np.sign(APF[n])) for n in C1]
seqs["X1hat_0_1_SC15-24"] = [int(np.sign(x1(n, 0, 1, True))) for n in C1]
# independence reference: all 2^n equally likely sign sequences (fair coin, cycles independent)
def indep(nc):
    L = [runs(s) for s in itertools.product([-1, 1], repeat=nc)]
    return dict(mean_runs=float(np.mean([len(x) for x in L])), mean_longest=float(np.mean([max(x) for x in L])),
                frac_longest_ge7=float(np.mean([max(x) >= 7 for x in L])))
ref = {nc: indep(nc) for nc in (10, 13, 14)}
out["runs"] = {k: dict(seq=v, run_lengths=runs(v), n_runs=len(runs(v)), longest=max(runs(v)),
                       independent_ref=ref[len(v)]) for k, v in seqs.items()}
for k, v in out["runs"].items():
    print(f"{k:32s} runs {v['run_lengths']}  (independent: {v['independent_ref']['mean_runs']:.1f} runs, longest {v['independent_ref']['mean_longest']:.1f})")

tic("lag-1 vs lag-2")
lag = {}
lag["A(n),A(n+1)"] = dict(r=r([AN[n] for n in range(12, 24)], [AN[n + 1] for n in range(12, 24)]), n=12)
lag["A(n),A(n+2)"] = dict(r=r([AN[n] for n in range(12, 23)], [AN[n + 2] for n in range(12, 23)]), n=11)
lag["A(n),A(n+2) incl prov SC25"] = dict(r=r([AN[n] for n in range(12, 24)], [AN[n + 2] for n in range(12, 24)]), n=12)
lag["X1(-1..4)(n),A(n+2)"] = dict(r=r(x1v(range(12, 23), -1, 4), [AN[n + 2] for n in range(12, 23)]), n=11)
lag["Apf(n),Apf(n+1)"] = dict(r=r([APF[n] for n in range(15, 24)], [APF[n + 1] for n in range(15, 24)]), n=9)
lag["Apf(n),Apf(n+2)"] = dict(r=r([APF[n] for n in range(15, 23)], [APF[n + 2] for n in range(15, 23)]), n=8)
lag["Apf(n),A(n+2)"] = dict(r=r([APF[n] for n in range(15, 24)], [AN[n + 2] for n in range(15, 24)]), n=9)
lag["Apf(n),A(n+1)"] = dict(r=r([APF[n] for n in range(15, 24)], [AN[n + 1] for n in range(15, 24)]), n=9)
# partial view without machinery: does A(n-1) add to A(n) for A(n+1)?  residual of A(n+1) on A(n), against A(n-1)
nn = list(range(13, 24)); y = np.array([AN[n + 1] for n in nn]); x = np.array([AN[n] for n in nn]); z = np.array([AN[n - 1] for n in nn])
res_y = y - np.polyval(np.polyfit(x, y, 1), x)
lag["resid A(n+1)|A(n) vs A(n-1)"] = dict(r=r(res_y, z), n=11)
out["lags"] = lag
for k, v in lag.items():
    print(f"{k:32s} r={v['r']:+.2f} n={v['n']}")

tic("flip and near-flip cycles: low-latitude large regions vs the rest")
LL = band_mask(d, 500, np.inf, 0, 20)
fl = {}
for n in range(12, 26):
    fl[n] = dict(A_all=AN[n], A_lowlat_large=cycA(n, LL), A_rest=cycA(n, ~LL),
                 share_area_lowlat_large=float(d.area[LL & (d.dec >= SCB[n]) & (d.dec < SCB.get(n + 1, 9999))].sum() /
                                               d.area[(d.dec >= SCB[n]) & (d.dec < SCB.get(n + 1, 9999))].sum()))
out["flip_table"] = fl
for n in (13, 14, 16, 17, 18, 20, 21, 23, 24):
    v = fl[n]; print(n, {k: round(x, 3) for k, x in v.items()})

tic("SFT: SC13/14 and SC20/21 with and without the low-latitude large regions")
M = Model()
gL = (np.abs(CAT.latc.values) <= 20) & (CAT.area.values >= 500)
tilts = [None] + [tilt_draw(CAT, 7000 + k) for k in range(1, 26)]
sft = {}
for lab, cat in (("full", CAT), ("no_lowlat_large", CAT[~gL].reset_index(drop=True)), ("only_lowlat_large", CAT[gL].reset_index(drop=True))):
    vals = []
    for t in tilts:
        tt = None if t is None else (t if lab == "full" else (t[~gL] if lab == "no_lowlat_large" else t[gL]))
        rr = M.run(cat, tt, ic=(lab == "full"))
        vals.append([[apf_of(rr["full"][n]), d_of(rr["full"][n]), apf_of(rr["own"][n]), d_of(rr["own"][n])] for n in CYC])
    vals = np.array(vals)                       # (members, cycles, 4)
    sft[lab] = dict(det=dict(zip(CYC, vals[0])), p16=dict(zip(CYC, np.percentile(vals[1:], 16, 0))),
                    p84=dict(zip(CYC, np.percentile(vals[1:], 84, 0))),
                    frac_neg_D=dict(zip(CYC, (vals[1:, :, 1] < 0).mean(0))))
out["sft_flip"] = sft
x10 = json.load(open(os.path.join(P.CR, "out", "x10_sc14.json")))
out["sc14_observed_mwo"] = {k: x10[k] for k in x10 if not isinstance(x10[k], (list, dict))}
for n in (13, 14, 20, 21):
    print(n, "full A_pf det %+.3f  D det %+.3e | no LL-large D %+.3e | only LL-large D %+.3e | frac(D<0) full %.2f" % (
        sft["full"]["det"][n][0], sft["full"]["det"][n][1], sft["no_lowlat_large"]["det"][n][1], sft["only_lowlat_large"]["det"][n][1],
        sft["full"]["frac_neg_D"][n]))
print("SC14 observed (MWO, x10):", out["sc14_observed_mwo"])

tic("DATA CHECK: DPD and Temmer hemispheric sunspot numbers")
DPD = os.path.join(r"C:\Users\Shivansh Mishra\Downloads\scaling polar field\To discuss in meeting\Final Datasets to be used\Tilt Datasets", "dpdtiltdata.txt")
rows = []
with open(DPD) as fh:
    for line in fh:
        if not line.startswith("g "):
            continue
        t = line.split()
        y, mo, dd = int(t[1]), int(t[2]), int(t[3]); wc = float(t[11]); lat = float(t[12])
        if wc <= 0 or lat == 0 or wc >= 99999:
            continue
        rows.append((y + (pd.Timestamp(y, mo, dd).dayofyear - 0.5) / (366 if pd.Timestamp(y, 1, 1).is_leap_year else 365), t[7], wc, lat))
dp = pd.DataFrame(rows, columns=["dec", "gid", "area", "lat"])
dp["north"] = dp.lat > 0; dp["rmax"] = dp.groupby("gid").area.transform("max")
im = dp.groupby("gid").area.idxmax(); dp["rlat"] = dp.gid.map(dp.loc[im].set_index("gid").lat.abs())
dp.to_pickle(os.path.join(OUT, "dpd_groups.pkl"))
def A_of(df, n, mask=None):
    m = (df.dec >= SCB[n]) & (df.dec < SCB[n + 1]) if mask is None else (df.dec >= SCB[n]) & (df.dec < SCB[n + 1]) & mask
    N = df.area[m & df.north].sum(); S = df.area[m & ~df.north].sum(); return float((N - S) / (N + S)), float(N + S)
dpd = {}
for n in (20, 21, 22, 23):
    a, tot = A_of(dp, n); aLL, _ = A_of(dp, n, (dp.rmax >= 500) & (dp.rlat <= 20)); aR, _ = A_of(dp, n, ~((dp.rmax >= 500) & (dp.rlat <= 20)))
    cover = float(max(0.0, min(dp.dec.max(), SCB[n + 1]) - max(dp.dec.min(), SCB[n])) / (SCB[n + 1] - SCB[n]))
    dpd[n] = dict(A_dpd=a, A_dpd_lowlat_large=aLL, A_dpd_rest=aR, A_rgo_soon=AN[n], A_rgo_lowlat_large=fl[n]["A_lowlat_large"],
                  coverage=cover)
    print(n, {k: round(v, 3) for k, v in dpd[n].items()})
# yearly N-S comparison over 1977-2013 (USAF/NOAA era) and 1974-1976 (RGO era)
yy = np.arange(1974, 2014)
def yearly(df_dec, df_area, df_north, years):
    return np.array([df_area[(df_dec >= y) & (df_dec < y + 1) & df_north].sum() - df_area[(df_dec >= y) & (df_dec < y + 1) & ~df_north].sum() for y in years])
dpd_y = yearly(dp.dec.values, dp.area.values, dp.north.values, yy)
rgo_y = yearly(D.dec.values, D.area.values, (D.hemi == "N").values, yy)
okyr = np.array([y not in (2014, 2015, 2016) for y in yy])
out["dpd_yearly_NminusS"] = dict(years=yy, dpd=dpd_y, rgo_soon=rgo_y, r=r(dpd_y, rgo_y), n=len(yy),
                                 r_1977_1986=r(dpd_y[(yy >= 1977) & (yy <= 1986)], rgo_y[(yy >= 1977) & (yy <= 1986)]))
out["dpd_cycles"] = dpd
print("yearly N-S, DPD vs RGO/SOON r =", round(out["dpd_yearly_NminusS"]["r"], 3), "n =", len(yy), "; 1977-86:", round(out["dpd_yearly_NminusS"]["r_1977_1986"], 3))
# Temmer et al. 2006 (Catalogue B, t06 columns 1945-1991)
EXT = os.path.join(P.REV, "external", "silso_hem")
CBN = ["date", "n_comb", "s_comb", "n_comb_sm", "s_comb_sm", "n_area", "s_area", "n_area_sm", "s_area_sm",
       "n_t06", "s_t06", "n_t06_sm", "s_t06_sm", "n_silso", "s_silso", "n_silso_sm", "s_silso_sm"]
cb = pd.read_csv(os.path.join(EXT, "Catalogue_B.txt"), sep=r"\s+", header=None, names=CBN)
cb["t"] = cb.date.str[:4].astype(int) + (cb.date.str[5:7].astype(int) - 0.5) / 12
cb = cb[(cb.n_t06 >= 0) & (cb.s_t06 >= 0)]
tem = {}
for n in (19, 20, 21):
    m = (cb.t >= SCB[n]) & (cb.t < SCB[n + 1]); a = cb[m]
    tem[n] = dict(A_SN_temmer=float((a.n_t06.sum() - a.s_t06.sum()) / (a.n_t06.sum() + a.s_t06.sum())),
                  months=int(m.sum()), A_area_rgo_soon=AN[n])
out["temmer_cycles"] = tem; print("Temmer:", tem)

tic("decay-time sweep (tau 2-17 yr)")
tau_list = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 13.0, 17.0]
sw = {}
obs = np.array([APF[n] for n in C1])
for tau in tau_list:
    Mt = M if tau == 5.0 else Model(tau=tau, base=M)
    rr = Mt.run(CAT, keep=True)
    a = np.array([apf_of(rr["full"][n]) for n in CYC]); Dd = np.array([d_of(rr["full"][n]) for n in CYC])
    old = []
    for n in CYC:
        sel, c = rr["reg"][n]; N, S = rr["full"][n]; dD = np.sign(N) * c[:, 0] - np.sign(S) * c[:, 1]
        o = CAT.cycle.values[sel] < n; old.append(float(dD[o].sum() / dD.sum()) if dD.sum() != 0 else np.nan)
    ens_r, ens_l1 = [], []
    for k in range(1, 26):
        re = Mt.run(CAT, tilts[k]); ae = np.array([apf_of(re["full"][n]) for n in CYC])
        ens_r.append(r(ae[3:], obs)); ens_l1.append(r(ae[:-1], ae[1:]))
    sw[tau] = dict(lag1_model_Apf=r(a[:-1], a[1:]), lag2_model_Apf=r(a[:-2], a[2:]), lag1_model_D=r(Dd[:-1], Dd[1:]),
                   old_cycle_share_of_D_median=float(np.nanmedian(old)), old_cycle_share_of_D_mean_abs=float(np.nanmean(np.abs(old))),
                   r_model_vs_obs_Apf=r(a[3:], obs), r_model_vs_obs_Apf_ens_median=float(np.median(ens_r)),
                   lag1_model_Apf_ens_median=float(np.median(ens_l1)), lag1_model_Apf_ens_16_84=np.percentile(ens_l1, [16, 84]).tolist(),
                   r_modelApf_n_vs_A_n1=r(a[:-1], [AN[n + 1] for n in CYC[:-1]]),
                   r_abs_modelApf_vs_abs_A_n1=r(np.abs(a[:-1]), [abs(AN[n + 1]) for n in CYC[:-1]]),
                   n_runs_model_Apf_sign=len(runs(list(np.sign(a)))), model_Apf=dict(zip(CYC, a)))
    v = sw[tau]
    print(f"tau {tau:4.1f}: lag1 A_pf {v['lag1_model_Apf']:+.2f} (ens {v['lag1_model_Apf_ens_median']:+.2f}) old share {v['old_cycle_share_of_D_median']:+.2f} "
          f"r obs {v['r_model_vs_obs_Apf']:+.2f} (ens {v['r_model_vs_obs_Apf_ens_median']:+.2f}) r(Apf_mod, A(n+1)) {v['r_modelApf_n_vs_A_n1']:+.2f} runs {v['n_runs_model_Apf_sign']}")
out["tau_sweep"] = sw
dump("q3_persistence.json", out)
tic("done")
