"""figs_q.py -- one figure per question (Q0-Q3), the ensemble figure, P6 and the exploratory figures."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from core import *

style()
q0 = load("q0_reconcile.json"); q1 = load("q1_sources.json"); q2 = load("q2_window.json"); q3 = load("q3_persistence.json")
ev = load("ens_verify.json"); ex = load("explore.json")
odd = lambda n: n % 2 == 1
CO = {True: "k", False: "0.55"}

def lab_pts(ax, xs, ys, ns, dx=0.0, dy=0.0, fs=5.5):
    for x_, y_, n in zip(xs, ys, ns):
        ax.annotate(str(n), (x_, y_), xytext=(2, 2), textcoords="offset points", fontsize=fs)

def fitline(ax, x, y, color="0.3"):
    b = np.polyfit(x, y, 1); xx = np.linspace(min(x), max(x), 10); ax.plot(xx, np.polyval(b, xx), "-", color=color, lw=0.7)

# ------------------------------------------------------------------------------------------ Q0 ----
tic("Q0")
fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.0))
xh = x1v(C1, 0, 1, True) / 1e4; x01 = x1v(C1, 0, 1) / 1e4; apf = np.array([APF[n] for n in C1]); S = np.array([S_SIGNED[n] for n in C1])
mk = [("o" if odd(n) else "s") for n in C1]
a_ = axs[0, 0]
for x_, y_, m_, n in zip(xh, apf, mk, C1):
    a_.plot(x_, y_, m_, color=WONG["N"] if y_ > 0 else WONG["S"], ms=4)
fitline(a_, xh, apf); lab_pts(a_, xh, apf, C1)
a_.axhline(0, color="0.7", lw=0.4); a_.axvline(0, color="0.7", lw=0.4)
a_.set_title(f"(a) parity-signed: r = {q0['C1_parity_signed']['r']:+.2f} (n = 10), signs 10/10", loc="left")
a_.set_xlabel("X1hat = (-1)^(n+1) X1, t_R+0->+1 (1e4)"); a_.set_ylabel("A_pf(n) = (|B_N|-|B_S|)/(|B_N|+|B_S|)")
a_ = axs[0, 1]
for x_, y_, m_, n in zip(x01, S, mk, C1):
    a_.plot(x_, y_, m_, color="0.3", ms=4)
fitline(a_, x01, S); lab_pts(a_, x01, S, C1)
a_.axhline(0, color="0.7", lw=0.4); a_.axvline(0, color="0.7", lw=0.4)
a_.set_title(f"(b) same relation, signed field: r = {q0['C1_signed_field_S']['r']:+.2f} (n = 10), signs 0/10", loc="left")
a_.text(0.03, 0.05, f"first differences: r = {q0['C1_signed_field_S_first_differences']['r']:+.2f} (n = 9)\n"
        "identity: X1hat*A_pf = -X1*S at all 10 minima", transform=a_.transAxes, fontsize=6)
a_.set_xlabel("X1, unsigned, t_R+0->+1 (1e4)"); a_.set_ylabel("S(n) = (B_N+B_S)/(|B_N|+|B_S|)")
a_ = axs[1, 0]
x14 = x1v(C2, -1, 4) / 1e4; an1 = np.array([AN[n + 1] for n in C2])
for x_, y_, n in zip(x14, an1, C2):
    a_.plot(x_, y_, "o", color=WONG["N"] if y_ > 0 else WONG["S"], ms=4)
fitline(a_, x14, an1); lab_pts(a_, x14, an1, C2)
a_.axhline(0, color="0.7", lw=0.4); a_.axvline(0, color="0.7", lw=0.4)
a_.set_title(f"(c) next cycle, unsigned: r = {q0['C2_unsigned']['r']:+.2f} (n = 12), signs 10/12", loc="left")
a_.set_xlabel("X1, unsigned, t_R-1->+4 (1e4)"); a_.set_ylabel("A(n+1), whole-cycle area")
a_ = axs[1, 1]
c9 = list(range(15, 24)); xa = np.abs([APF[n] for n in c9]); ya = np.abs([AN[n + 1] for n in c9])
a_.plot(xa, ya, "o", color="0.3", ms=4); fitline(a_, xa, ya); lab_pts(a_, xa, ya, c9)
a_.set_title(f"(d) magnitudes: r = {q0['mag_|Apf(n)|_vs_|A(n+1)|']['r']:+.2f} (n = 9)", loc="left")
a_.text(0.03, 0.88, f"|A(n)| vs |A(n+1)|: r = {q0['mag_|A(n)|_vs_|A(n+1)|']['r']:+.2f} (n = 12)", transform=a_.transAxes, fontsize=6)
a_.set_xlabel("|A_pf(n)|"); a_.set_ylabel("|A(n+1)|")
fig.text(0.01, -0.01, "Circles: odd cycles, squares: even cycles in (a, b).  Colour: sign of the target (blue N, orange S).", fontsize=6)
fig.tight_layout(); savefig(fig, "Q0_reconcile"); plt.close(fig)

# ------------------------------------------------------------------------------------------ Q1 ----
tic("Q1")
sr = q1["series"]; cyc = C1; xs = np.arange(len(cyc))
fig = plt.figure(figsize=(7.1, 5.4)); gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.45, wspace=0.42)
a_ = fig.add_subplot(gs[0, :])
g = lambda k: np.array([sr[k][str(n)] for n in cyc])
a_.fill_between(xs, g("all_p16"), g("all_p84"), color=WONG["giant"], alpha=0.25, lw=0, label="(d) all: area + tilt scatter + transport asym. (16-84 %)")
a_.fill_between(xs, g("tilt_only_p16"), g("tilt_only_p84"), color="0.5", alpha=0.25, lw=0, label="(b) tilt scatter only (16-84 %)")
for k_, v in q1["c_transport_only"].items():
    a_.plot(xs, [v["apf"][str(n)] for n in cyc], "-", color=WONG["purple"], lw=0.5, alpha=0.7, label="(c) transport only (8 configs)" if k_ == "eta-0.2" else None)
a_.plot(xs, g("area"), "-o", color=WONG["N"], ms=3, lw=1.1, label="(a) area only (real emergence, Joy tilts)")
a_.plot(xs, g("amount"), "--", color=WONG["N"], lw=0.8, label="(a1) amount only")
a_.plot(xs, g("placement"), ":", color=WONG["N"], lw=1.0, label="(a2) placement only (flux equalised)")
a_.errorbar(xs, g("obs"), yerr=[SIG[n] for n in cyc], fmt="D", color="k", ms=3.5, lw=0.7, capsize=1.5, label="observed A_pf (1 sigma)", zorder=5)
a_.axhline(0, color="0.6", lw=0.5); a_.set_xticks(xs); a_.set_xticklabels([f"SC{n}" for n in cyc]); a_.set_ylabel("A_pf at the minimum ending n")
a_.set_ylim(-0.22, 0.5); a_.legend(ncol=2, fontsize=5.8, loc="upper left"); a_.set_title("(a) modelled polar-field asymmetry with one hemispheric difference switched on at a time", loc="left")
exps = [("(a) area only", q1["a_area"]["r_obs"], None, q1["a_area"]["rms"]),
        ("(a1) amount only", q1["a1_amount_only"]["r_obs"], None, q1["a1_amount_only"]["rms"]),
        ("(a2) placement only", q1["a2_placement_only"]["r_obs"], None, q1["a2_placement_only"]["rms"]),
        ("(b) tilt only", q1["b_tilt_only"]["r_obs_median"], q1["b_tilt_only"]["r_obs_16_84"], q1["b_tilt_only"]["rms_median"]),
        ("(c) transport only", float(np.median([v["r_obs"] for v in q1["c_transport_only"].values()])),
         [min(v["r_obs"] for v in q1["c_transport_only"].values()), max(v["r_obs"] for v in q1["c_transport_only"].values())],
         float(np.median([v["rms"] for v in q1["c_transport_only"].values()]))),
        ("area + tilt scatter", q1["area_plus_tilt_scatter"]["r_obs_median"], q1["area_plus_tilt_scatter"]["r_obs_16_84"], q1["area_plus_tilt_scatter"]["rms_median"]),
        ("(d) all", q1["d_all"]["r_obs_median"], q1["d_all"]["r_obs_16_84"], q1["d_all"]["rms_median"])]
a_ = fig.add_subplot(gs[1, 0]); ys = np.arange(len(exps))[::-1]
for y_, (nm, rv, rng_, _) in zip(ys, exps):
    if rng_ is not None:
        a_.plot(rng_, [y_, y_], "-", color="0.4", lw=1.2)
    a_.plot(rv, y_, "o", color=WONG["N"], ms=4)
a_.axvline(0, color="0.6", lw=0.5); a_.set_yticks(ys); a_.set_yticklabels([e[0] for e in exps], fontsize=6)
a_.set_xlabel("r(model, observed A_pf), n = 10"); a_.set_xlim(-0.6, 1); a_.set_title("(b) agreement", loc="left")
a_ = fig.add_subplot(gs[1, 1])
for y_, (nm, _, _, rms) in zip(ys, exps):
    a_.barh(y_, rms, color=WONG["N"], height=0.6)
a_.axvline(np.sqrt(np.mean(np.array([APF[n] for n in cyc]) ** 2)), color="k", lw=0.7, ls="--")
a_.set_yticks(ys); a_.set_yticklabels([]); a_.set_xlabel("rms of modelled A_pf"); a_.set_title("(c) amplitude (dashed: observed)", loc="left")
a_ = fig.add_subplot(gs[1, 2]); al = q1["algebraic_hemi_eta"]; dd = sorted(map(float, al))
a_.plot(dd, [al[f"{d}" if f"{d}" in al else str(d)]["r_obs"] for d in dd], "-o", color=WONG["green"], ms=2.5)
a_.axhline(q1["X1hat"]["r"], color="k", lw=0.6, ls="--"); a_.text(-0.6, q1["X1hat"]["r"] + 0.02, f"X1hat {q1['X1hat']['r']:.2f}", fontsize=5.5)
a_.axhline(q1["X1hat_times_Tlin"]["r"], color=WONG["purple"], lw=0.6, ls=":"); a_.text(-0.6, q1["X1hat_times_Tlin"]["r"] - 0.09, f"X1hat*T_lin {q1['X1hat_times_Tlin']['r']:.2f}", fontsize=5.5, color=WONG["purple"])
a_.text(0.0, 0.47, "algebraic dipole,\neta_N != eta_S", ha="center", fontsize=5.5, color=WONG["green"])
a_.set_xlabel("(eta_N - eta_S)/(eta_N + eta_S)"); a_.set_ylabel("r with observed A_pf (n = 10)", fontsize=6.5); a_.set_ylim(-0.2, 1.0)
a_.set_title("(d) algebraic, eta per hemisphere", loc="left")
savefig(fig, "Q1_what_sets_Apf"); plt.close(fig)

# ------------------------------------------------------------------------------------------ Q3 flips ----
tic("Q3 flips")
show = [(14, "flip S->N"), (21, "flip N->S"), (17, "near flip"), (18, "near flip"), (24, "near flip")]
fig, axs = plt.subplots(1, 5, figsize=(7.1, 2.5), sharey=True)
ft = q3["flip_table"]
for a_, (n, tag) in zip(axs, show):
    m = (CAT.cycle.values == n); c = CAT[m]; LLm = (np.abs(c.latc) <= 20) & (c.area >= 500)
    t = c.year - SCB[n]
    a_.scatter(t[~LLm], c.latc[~LLm], s=1.2, c="0.75", lw=0, rasterized=True)
    a_.scatter(t[LLm], c.latc[LLm], s=c.area[LLm] / 40, c=np.where(c.latc[LLm] > 0, WONG["N"], WONG["S"]), alpha=0.7, lw=0.3, edgecolor="w")
    a_.axhline(0, color="0.5", lw=0.4); a_.axvline(TR[n] - SCB[n], color="k", lw=0.5, ls=":")
    v = ft[str(n)]
    a_.set_title(f"SC{n} ({tag})", loc="left", fontsize=7)
    a_.text(0.02, 0.02, f"all {v['A_all']:+.3f}\nlow-lat >=500 {v['A_lowlat_large']:+.3f}\nrest {v['A_rest']:+.3f}", transform=a_.transAxes, fontsize=5.3, va="bottom")
    a_.set_xlabel("yr since minimum", fontsize=6); a_.set_xlim(-1, 12.5); a_.set_ylim(-45, 45)
axs[0].set_ylabel("latitude (deg)")
fig.text(0.01, -0.09, "Coloured: regions with peak area >= 500 MSH within 20 deg of the equator (area ~ symbol size; blue N, orange S); grey: all other regions. Dotted line: t_R.", fontsize=6)
savefig(fig, "Q3a_flip_butterflies"); plt.close(fig)

# ------------------------------------------------------------------------------------------ Q3 SFT + data check ----
tic("Q3 sft and data check")
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.6), gridspec_kw=dict(width_ratios=[1.25, 1.25, 1]))
sf = q3["sft_flip"]; cc = CYC; xs = np.arange(len(cc))
a_ = axs[0]
full = np.array([sf["full"]["det"][str(n)][1] for n in cc]); nol = np.array([sf["no_lowlat_large"]["det"][str(n)][1] for n in cc]); oll = np.array([sf["only_lowlat_large"]["det"][str(n)][1] for n in cc])
a_.bar(xs - 0.2, oll, 0.38, color=WONG["giant"], label="low-lat >= 500 MSH only")
a_.bar(xs + 0.2, nol, 0.38, color=WONG["rest"], label="all other regions")
a_.plot(xs, full, "kD", ms=3, label="full model D")
a_.fill_between(xs, [sf["full"]["p16"][str(n)][1] for n in cc], [sf["full"]["p84"][str(n)][1] for n in cc], color="0.5", alpha=0.25, lw=0, label="tilt scatter 16-84 %")
a_.axhline(0, color="0.5", lw=0.5); a_.set_xticks(xs); a_.set_xticklabels([str(n) for n in cc], fontsize=6)
for n in (14, 21):
    a_.axvspan(cc.index(n) - 0.5, cc.index(n) + 0.5, color="0.9", zorder=0)
a_.set_ylabel("modelled |B_N|-|B_S| (G)"); a_.set_xlabel("cycle"); a_.legend(fontsize=5.3, loc="lower left")
a_.set_title("(a) SFT polar imbalance at the minimum ending n", loc="left")
a_ = axs[1]
dy = q3["dpd_yearly_NminusS"]; yy = np.array(dy["years"]); ok = ~np.isin(yy, [2014, 2015, 2016])
a_.plot(yy[ok] + 0.5, np.array(dy["rgo_soon"])[ok] / 1e5, "-", color="k", lw=0.9, label="RGO -> USAF/NOAA (RGO_SOON)")
a_.plot(yy[ok] + 0.5, np.array(dy["dpd"])[ok] / 1e5, "--", color=WONG["green"], lw=1.0, label="Debrecen DPD")
a_.axvline(1977, color=WONG["S"], lw=0.7, ls=":"); a_.text(1977.5, 0.03, "RGO -> USAF splice", transform=a_.get_xaxis_transform(), fontsize=5.5, color=WONG["S"])
for n in (21, 22, 23):
    a_.axvline(SCB[n], color="0.85", lw=0.5, zorder=0)
a_.axhline(0, color="0.5", lw=0.4); a_.set_xlabel("year"); a_.set_ylabel("yearly A_N - A_S (1e5 MSH day)")
a_.set_title(f"(b) yearly N-S: r = {dy['r']:.2f} (n = {dy['n']})", loc="left"); a_.legend(fontsize=5.3, loc="upper left", bbox_to_anchor=(0.12, 1.0))
a_ = axs[2]
dpd = q3["dpd_cycles"]; tem = q3["temmer_cycles"]
rows_ = [("SC20 RGO", AN[20], "k"), ("SC20 Temmer SN", tem["20"]["A_SN_temmer"], WONG["purple"]),
         ("SC21 RGO/USAF", AN[21], "k"), ("SC21 DPD", dpd["21"]["A_dpd"], WONG["green"]), ("SC21 Temmer SN", tem["21"]["A_SN_temmer"], WONG["purple"]),
         ("SC21 low-lat>=500 RGO/USAF", dpd["21"]["A_rgo_lowlat_large"], WONG["giant"]), ("SC21 low-lat>=500 DPD", dpd["21"]["A_dpd_lowlat_large"], WONG["giant"]),
         ("SC22 RGO/USAF", AN[22], "k"), ("SC22 DPD", dpd["22"]["A_dpd"], WONG["green"])]
yv = np.arange(len(rows_))[::-1]
for y_, (nm, v, c_) in zip(yv, rows_):
    a_.plot([0, v], [y_, y_], "-", color=c_, lw=1.0); a_.plot(v, y_, "o", color=c_, ms=3.5)
    a_.text(0.29, y_, f"{v:+.3f}", fontsize=5.5, va="center", ha="right")
a_.axvline(0, color="0.5", lw=0.5); a_.set_yticks(yv); a_.set_yticklabels([r_[0] for r_ in rows_], fontsize=5.5); a_.set_xlim(-0.2, 0.3)
a_.set_xlabel("whole-cycle asymmetry (N-S)/(N+S)"); a_.set_title("(c) the SC21 flip, three records", loc="left")
fig.tight_layout(); savefig(fig, "Q3b_sft_and_datacheck"); plt.close(fig)

# ------------------------------------------------------------------------------------------ Q3 persistence + tau ----
tic("Q3 persistence")
fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.5), gridspec_kw=dict(width_ratios=[1.4, 1]))
a_ = axs[0]; R = q3["runs"]
rowsP = [("all regions", "area_all_SC12-24", list(range(12, 25))), ("peak >= 500 MSH", "area_large>=500", list(range(12, 25))),
         ("low-lat >= 500 MSH", "lowlat_large>=500_|lat|<=20", list(range(12, 25))), ("100-499 MSH", "area_medium100-499", list(range(12, 25))),
         ("< 100 MSH", "area_small<100", list(range(12, 25))), ("A_pf (polar)", "Apf_SC15-24", C1)]
for i, (nm, key, cyc_) in enumerate(rowsP[::-1]):
    for n, s_ in zip(cyc_, R[key]["seq"]):
        a_.add_patch(plt.Rectangle((n - 0.45, i - 0.4), 0.9, 0.8, color=WONG["N"] if s_ > 0 else WONG["S"], lw=0))
        a_.text(n, i, "N" if s_ > 0 else "S", ha="center", va="center", fontsize=5, color="w")
    a_.text(25.7, i, f"{R[key]['n_runs']} runs (indep. {R[key]['independent_ref']['mean_runs']:.1f})", fontsize=5.5, va="center")
a_.set_xlim(11.4, 29.5); a_.set_ylim(-0.6, len(rowsP) - 0.4); a_.set_yticks(range(len(rowsP))); a_.set_yticklabels([r_[0] for r_ in rowsP[::-1]], fontsize=6)
a_.set_xticks(range(12, 25)); a_.set_xticklabels([str(n) for n in range(12, 25)], fontsize=6); a_.spines[:].set_visible(False)
a_.set_title("(a) dominant hemisphere per cycle", loc="left")
a_ = axs[1]; sw = q3["tau_sweep"]; taus = sorted(map(float, sw))
gk = lambda k: [sw[str(t) if str(t) in sw else f"{t}"][k] for t in taus]
lo_ = [sw[str(t)]["lag1_model_Apf_ens_16_84"][0] for t in taus]; hi_ = [sw[str(t)]["lag1_model_Apf_ens_16_84"][1] for t in taus]
a_.fill_between(taus, lo_, hi_, color=WONG["N"], alpha=0.2, lw=0)
a_.plot(taus, gk("lag1_model_Apf_ens_median"), "-o", color=WONG["N"], ms=2.5, label="lag-1 r of modelled A_pf (n = 12)")
a_.plot(taus, gk("old_cycle_share_of_D_median"), "-s", color="0.4", ms=2.5, label="share of D from older cycles")
a_.plot(taus, gk("r_model_vs_obs_Apf_ens_median"), "-^", color=WONG["green"], ms=2.5, label="r(model, obs A_pf) (n = 10)")
a_.axhline(q3["lags"]["A(n),A(n+1)"]["r"], color="k", lw=0.6, ls="--"); a_.text(17, q3["lags"]["A(n),A(n+1)"]["r"] + 0.03, "data: lag-1 r of area A", fontsize=5.5, ha="right")
a_.axhline(0, color="0.6", lw=0.4); a_.set_xlabel("surface decay time tau (yr)"); a_.set_ylim(-0.5, 0.68); a_.legend(fontsize=5.3, loc="lower center")
a_.set_title("(b) SFT memory vs decay time", loc="left")
fig.tight_layout(); savefig(fig, "Q3c_persistence_tau"); plt.close(fig)

# ------------------------------------------------------------------------------------------ ensemble figure ----
tic("ensemble figure")
ed = np.load(os.path.join(OUT, "ens_data.npz")); em = np.load(os.path.join(OUT, "ens_model_emergence.npz"))
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.6))
for a_, XE, YE, xo, yo, ttl, key in ((axs[0], ed["XHE"], None, x1v(C1, 0, 1, True), np.array([APF[n] for n in C1]), "(a) H1: X1hat vs A_pf", "H1_emergence"),
                                      (axs[1], ed["X1E"], ed["AE"], x1v(C2, -1, 4), np.array([AN[n + 1] for n in C2]), "(b) H2: X1 vs A(n+1)", "H2_emergence")):
    xx = np.linspace(xo.min() * 1.1, xo.max() * 1.1, 50) / 1e4; lines = []
    for k in range(len(XE)):
        ytar = yo if YE is None else YE[k]
        b = np.polyfit(XE[k] / 1e4, ytar, 1); lines.append(np.polyval(b, xx))
    lines = np.array(lines)
    if YE is None:                      # H1: the A_pf errors are the spread that matters; area noise is a hairline
        rg = np.random.default_rng(3); sg = np.array([SIG[n] for n in C1]); le = []
        for k in range(300):
            b = np.polyfit(xo / 1e4, yo + rg.normal(0, 1, len(yo)) * sg, 1); le.append(np.polyval(b, xx))
        a_.fill_between(xx, np.percentile(le, 16, 0), np.percentile(le, 84, 0), color="0.6", alpha=0.3, lw=0, label="fit, A_pf redrawn within 1 sigma (16-84 %)")
    a_.fill_between(xx, np.percentile(lines, 16, 0), np.percentile(lines, 84, 0), color=WONG["N"], alpha=0.5, lw=0, label="fit, area noise (16-84 %)")
    a_.legend(fontsize=5.2, loc="upper left")
    xlo, xhi = np.percentile(XE, 16, 0) / 1e4, np.percentile(XE, 84, 0) / 1e4
    yerr = np.array([SIG[n] for n in C1]) if YE is None else np.clip(np.vstack([yo - np.percentile(YE, 16, 0), np.percentile(YE, 84, 0) - yo]), 0, None)
    xerr = np.clip(np.vstack([xo / 1e4 - xlo, xhi - xo / 1e4]), 0, None)      # observed value can sit at the edge of the band
    a_.errorbar(xo / 1e4, yo, xerr=xerr, yerr=yerr, fmt="o", ms=3, color="k", lw=0.6, capsize=1)
    e = ev[key]
    if key == "H1_emergence":
        xa = load("ens_extra.json")
        sub = (f"area noise [{xa['H1_area_only']['p16_84'][0]:.2f}, {xa['H1_area_only']['p16_84'][1]:.2f}]; "
               f"A_pf errors [{xa['H1_apf_error_only']['p16_84'][0]:.2f}, {xa['H1_apf_error_only']['p16_84'][1]:.2f}]")
    else:
        sub = f"area noise 16-84 % [{e['p16_84'][0]:.2f}, {e['p16_84'][1]:.2f}]"
    a_.set_title(f"{ttl}: r = {e['r_obs']:+.2f} (n = {e['n']})\n{sub}", loc="left", fontsize=6.3)
    a_.axhline(0, color="0.7", lw=0.4); a_.axvline(0, color="0.7", lw=0.4)
axs[0].set_xlabel("X1hat, t_R+0->+1 (1e4)"); axs[0].set_ylabel("A_pf"); axs[1].set_xlabel("X1, t_R-1->+4 (1e4)"); axs[1].set_ylabel("A(n+1)")
a_ = axs[2]; A_ = em["A"]; cyc_m = list(em["cycles"]); xs = np.array(cyc_m)
a_.fill_between(xs, np.percentile(A_, 16, 0), np.percentile(A_, 84, 0), color=WONG["giant"], alpha=0.3, lw=0, label="emergence ensemble 16-84 %")
tsr = np.load(os.path.join(OUT, "ens_transport_series.npz"))
a_.fill_between(tsr["cycles"], tsr["apf"].min(0), tsr["apf"].max(0), color=WONG["purple"], alpha=0.25, lw=0, label="transport ensemble, 27 configs (range)")
a_.plot(xs, A.ensemble_apf("A1", cyc_m, member=0), "-", color=WONG["giant"], lw=1.0, label="model, deterministic")
a_.errorbar(C1, [APF[n] for n in C1], yerr=[SIG[n] for n in C1], fmt="kD", ms=3, lw=0.6, capsize=1, label="observed A_pf")
a_.axhline(0, color="0.6", lw=0.5); a_.axvspan(20.5, 21.5, color="0.9", zorder=0)
a_.set_xlabel("cycle n"); a_.set_ylabel("A_pf"); a_.legend(fontsize=5.3, loc="upper left")
ts = ev["transport_summary"]
a_.set_title(f"(c) H3: SC21 polar sign S in {ev['model_emergence']['frac_D21_negative']*100:.0f} % of members,\n"
             f"{ts['D21_negative_configs']}/{ts['n_configs']} transport configs", loc="left", fontsize=6.8)
fig.tight_layout(); savefig(fig, "ENS_headlines"); plt.close(fig)

# ------------------------------------------------------------------------------------------ P6 concentration ----
tic("P6")
cz = np.load(os.path.join(OUT, "explore_concentration.npz"))
fig, ax = plt.subplots(figsize=(4.2, 2.8))
for n in CYC:
    c = cz[f"c{n}"]; k = np.arange(1, len(c) + 1)
    hl = n in (14, 21)
    ax.plot(k, c, "-", color=(WONG["S"] if n == 21 else WONG["N"]) if hl else "0.7", lw=1.3 if hl else 0.7, zorder=3 if hl else 1)
    if hl or n in (19, 23):
        ax.text(k[-1] * 1.05, c[-1], f"SC{n}", fontsize=5.5, va="center")
ax.axhline(1, color="k", lw=0.6, ls="--"); ax.set_xscale("log"); ax.set_xlim(1, 400)
ax.set_xlabel("k largest contributors (ranked by contribution to the imbalance)"); ax.set_ylabel("cumulative share of |D|")
cs = ex["C_summary"]
ax.set_title(f"Concentration: |D| reached by k = {cs['k_min']}-{cs['k_max']} regions (median {cs['k_median']:.0f})", loc="left", fontsize=7)
savefig(fig, "P6_concentration"); plt.close(fig)

# ------------------------------------------------------------------------------------------ exploratory ----
tic("exploratory")
fig, axs = plt.subplots(1, 4, figsize=(7.1, 2.4), gridspec_kw=dict(width_ratios=[1.4, 1, 1, 1]))
cf = ex["A_counterfactual"]; xs = np.arange(len(CYC))
a_ = axs[0]
a_.bar(xs - 0.2, [cf[str(n)]["dD_removed_frac"] for n in CYC], 0.38, color=WONG["giant"], label="removed")
a_.bar(xs + 0.2, [cf[str(n)]["dD_moved_frac"] for n in CYC], 0.38, color=WONG["purple"], label="moved to 35 deg")
for i, n in enumerate(CYC):
    if cf[str(n)]["sign_flip_removed"] or cf[str(n)]["sign_flip_moved"]:
        a_.text(i, 0.02 + max(0, cf[str(n)]["dD_removed_frac"]), "flip", fontsize=5, ha="center")
a_.axhline(0, color="0.5", lw=0.5); a_.set_xticks(xs); a_.set_xticklabels([str(n) for n in CYC], fontsize=5.5)
a_.set_ylabel("change in D / |D|"); a_.set_xlabel("cycle"); a_.legend(fontsize=5.3, loc="lower left"); a_.set_title("(A) largest low-lat region", loc="left")
a_ = axs[1]; B = ex["B_XG"]; xg = np.array([B["XG"][str(n)] for n in C1]); apf = np.array([APF[n] for n in C1])
a_.scatter(xg / np.abs(xg).max(), apf, c=["k" if odd(n) else "0.6" for n in C1], s=10)
lab_pts(a_, xg / np.abs(xg).max(), apf, C1, fs=5)
a_.axhline(0, color="0.7", lw=0.4); a_.axvline(0, color="0.7", lw=0.4)
a_.set_xlabel("X_G (normalised)"); a_.set_ylabel("A_pf")
a_.set_title(f"(B) r = {B['r_XG_vs_obsApf']:+.2f}; X1hat {B['r_X1hat_vs_obsApf']:+.2f}", loc="left", fontsize=6.5)
a_ = axs[2]; Dm = ex["D_memory"]
names = [("fair coin", dict(mean_longest_run=q3["runs"]["area_all_SC12-24"]["independent_ref"]["mean_longest"])),
         ("white noise\n(obs. mean)", ex["D_sim_white"]), ("AR(1)", ex["D_sim_ar1"]), ("AR(2) fit", ex["D_sim_ar2"])]
a_.bar(range(4), [v["mean_longest_run"] for _, v in names], color=["0.85", "0.7", WONG["rest"], WONG["N"]])
a_.axhline(7, color="k", lw=0.7, ls="--"); a_.text(3.4, 7.1, "observed 7", fontsize=5.5, ha="right")
a_.set_xticks(range(4)); a_.set_xticklabels([n_ for n_, _ in names], fontsize=5); a_.set_ylabel("mean longest run (13 cycles)")
a_.set_title(f"(D) a = {Dm['a']:+.2f}, b = {Dm['b']:+.2f}", loc="left", fontsize=6.5)
a_ = axs[3]; E_ = ex["E_timing"]; lagv = np.array([E_["per_cycle"][str(n)]["lag_N_minus_S"] for n in range(12, 25)])
Av = np.array([AN[n] for n in range(12, 25)])
a_.scatter(lagv, Av, c=[WONG["N"] if v > 0 else WONG["S"] for v in Av], s=10); lab_pts(a_, lagv, Av, range(12, 25), fs=5)
a_.axhline(0, color="0.7", lw=0.4); a_.axvline(0, color="0.7", lw=0.4)
a_.set_xlabel("$t_N - t_S$ (yr, area centroid)"); a_.set_ylabel("A(n)")
a_.set_title(f"(E) r = {E_['r_lag_vs_A']:+.2f} (n = 13)\nvs A(n+1) {E_['r_lag_vs_A_next']:+.2f} (n = 12)", loc="left", fontsize=6.3)
fig.tight_layout(); savefig(fig, "EXP_A_B_D_E"); plt.close(fig)
tic("done")
