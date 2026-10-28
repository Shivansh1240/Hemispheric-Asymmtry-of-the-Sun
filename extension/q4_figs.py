"""q4_figs.py -- Q4 figures: randomised emergence (SFT), the key correlations without the largest regions, and the
self-excited 2-D dynamo with hemispheric coupling."""
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from core import *

style()
q4 = load("q4_stochastic.json"); Z = np.load(os.path.join(OUT, "q4_randomised.npz"))
obs = Z["obs"]; det = Z["det"]; xs = np.arange(len(C1))
LAB = {"T_tilt_only": "real catalogue, tilt scatter", "R1_size_shuffle": "R1: sizes re-dealt in each hemisphere",
       "R2_slot_shuffle": "R2: timing/latitude re-dealt across hemispheres"}
COL = {"T_tilt_only": WONG["N"], "R1_size_shuffle": WONG["green"], "R2_slot_shuffle": WONG["S"]}

# ------------------------------------------------------------------- Q4a randomised emergence ----
tic("Q4a")
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.6), gridspec_kw=dict(width_ratios=[1.7, 0.9, 1.1]))
a_ = axs[0]
for k, off in (("T_tilt_only", -0.22), ("R1_size_shuffle", 0.0), ("R2_slot_shuffle", 0.22)):
    E = Z[k]; lo, md, hi = np.percentile(E, [16, 50, 84], 0)
    a_.errorbar(xs + off, md, yerr=[md - lo, hi - md], fmt="o", ms=2.5, lw=0.8, color=COL[k], label=LAB[k])
a_.plot(xs, obs, "kD", ms=3.5, label="observed A_pf", zorder=5)
a_.axhline(0, color="0.6", lw=0.5); a_.set_xticks(xs); a_.set_xticklabels([str(n) for n in C1]); a_.set_xlabel("cycle")
a_.set_ylabel("A_pf at the minimum ending n"); a_.legend(fontsize=5.2, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2)
a_.set_title("(a) 100 members each, hemispheric totals kept (16-50-84 %)", loc="left", fontsize=7)
a_ = axs[1]
for i, k in enumerate(("T_tilt_only", "R1_size_shuffle", "R2_slot_shuffle")):
    rr = [r(e, obs) for e in Z[k]]
    a_.boxplot([rr], positions=[i], widths=0.5, showfliers=False, patch_artist=True,
               boxprops=dict(facecolor=COL[k], alpha=0.35, lw=0.6), medianprops=dict(color="k", lw=0.8), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6))
a_.axhline(0, color="0.6", lw=0.5); a_.set_xticks(range(3)); a_.set_xticklabels(["tilt", "R1", "R2"]); a_.set_ylabel("r(member, observed A_pf), n = 10")
a_.set_title("(b) agreement", loc="left", fontsize=7)
a_ = axs[2]
for i, k in enumerate(("T_tilt_only", "R1_size_shuffle", "R2_slot_shuffle")):
    v = q4[k]
    a_.barh(i, v["signal_fraction"], color=COL[k], height=0.55)
    a_.text(v["signal_fraction"] + 0.02, i, f"{v['signal_fraction']:.2f}; obs sign kept in {v['cycles_sign_kept_84pct_obs']}/10", va="center", fontsize=5.5)
a_.set_yticks(range(3)); a_.set_yticklabels(["tilt", "R1", "R2"]); a_.set_xlim(0, 1.25); a_.set_xlabel("predetermined fraction of variance")
a_.set_title("(c) set by the catalogue vs noise", loc="left", fontsize=7)
fig.tight_layout(); savefig(fig, "Q4a_randomised_emergence"); plt.close(fig)

# ------------------------------------------------------------------- Q4b without the largest regions ----
tic("Q4b")
e2 = q4["E2_data"]; ks = [0, 1, 2, 3, 5, 10, 20]
key = lambda k: "all" if k == 0 else f"without_top{k}"
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.6), gridspec_kw=dict(width_ratios=[1.3, 1.1, 0.9]))
a_ = axs[2]
nruns = [len(e2[key(k)]["runs"]) for k in ks]; a21 = [e2[key(k)]["A21"] for k in ks]; a17 = []
xk = np.arange(len(ks))
a_.plot(xk, nruns, "-o", ms=3, color="0.3")
a_.axhline(7.0, color="0.6", lw=0.6, ls=":"); a_.text(len(ks) - 1, 7.15, "independent cycles: 7.0", fontsize=5.3, ha="right", color="0.4")
a_.set_xticks(xk); a_.set_xticklabels([str(k) for k in ks]); a_.set_ylim(0, 9)
a_.set_xlabel("largest regions removed per cycle (k)"); a_.set_ylabel("runs of the dominant hemisphere")
for x_, v in zip(xk, a21):
    a_.text(x_, 0.6, f"{v:+.3f}".replace("0.", "."), fontsize=4.8, ha="center", color=WONG["S"] if v < 0 else WONG["N"])
a_.text(-0.4, 1.3, "A(21) below:", fontsize=5.3)
a_.set_title("(c) the 3-run record needs single regions", loc="left", fontsize=6.5)
a_ = axs[0]
for nm, lab, c in (("H1_x1hat_apf", "X1hat -> A_pf (n = 10)", WONG["N"]), ("H2_x1_next", "X1 -> A(n+1) (n = 12)", WONG["S"]),
                   ("decline_unsigned_apf", "decline X1 -> A_pf, unsigned (n = 10)", WONG["green"]), ("persistence_lag1", "A(n) -> A(n+1) (n = 12)", "0.4")):
    a_.plot(ks, [e2[key(k)][nm] for k in ks], "-o", ms=3, color=c, label=lab)
    if e2.get("only_top10"):
        a_.plot(10.6, e2["only_top10"][nm], "x", color=c, ms=5)
a_.set_xscale("symlog", linthresh=1); a_.set_xticks(ks); a_.set_xticklabels([str(k) for k in ks])
a_.set_xlabel("largest regions removed per cycle (k)"); a_.set_ylabel("Pearson r"); a_.set_ylim(-0.2, 1.0); a_.axhline(0, color="0.6", lw=0.5)
a_.legend(fontsize=5.3, loc="lower left"); a_.set_title("(a) data: with and without the largest regions (x: top 10 only)", loc="left", fontsize=7)
a_ = axs[1]; s2 = q4["E2_sft"]; ks2 = [0, 1, 3, 5, 10, 20]
a_.plot(ks2, [s2["all"]["r_obs"]] + [s2[f"without_top{k}"]["r_obs"] for k in ks2[1:]], "-o", ms=3, color=WONG["N"], label="SFT without top k vs observed A_pf")
a_.plot(ks2[1:], [s2[f"only_top{k}"]["r_obs"] for k in ks2[1:]], "-s", ms=3, color=WONG["giant"], label="SFT with ONLY the top k vs observed")
a_.plot(ks2[1:], [s2[f"without_top{k}"]["r_with_full_model"] for k in ks2[1:]], "--", color="0.5", lw=0.9, label="without top k vs full model")
a_.set_xscale("symlog", linthresh=1); a_.set_xticks(ks2); a_.set_xticklabels([str(k) for k in ks2]); a_.axhline(0, color="0.6", lw=0.5)
a_.set_xlabel("largest regions per cycle (k)"); a_.set_ylabel("Pearson r (n = 10)"); a_.set_ylim(-0.4, 1.05); a_.legend(fontsize=5.3, loc="lower left")
a_.set_title("(b) SFT", loc="left", fontsize=7)
fig.tight_layout(); savefig(fig, "Q4b_without_largest"); plt.close(fig)

# ------------------------------------------------------------------- Q4c dynamo ----
DY = os.path.join(OUT, "dynamo")
runs_ = sorted(glob.glob(os.path.join(DY, "prod_*_stats.json")))
if runs_:
    tic("Q4c")
    from q4_dynamo_analyse import per_cycle
    S = {os.path.basename(p)[5:-11]: json.load(open(p)) for p in runs_}
    base = "k1" if "k1" in S else sorted(S)[0]
    TICK = {"k03": r"$\kappa$=0.3", "k1": r"$\kappa$=1", "k3": r"$\kappa$=3", "k1s2": r"$\kappa$=1" + "\nseed 2",
            "rho05": r"$\rho$=0.5", "rho08": r"$\rho$=0.8", "k1lr": r"$\kappa$=1" + "\nlow res", "k10lr": r"$\kappa$=10" + "\nlow res"}
    res = dict(np.load(os.path.join(DY, f"prod_{base}.npz"), allow_pickle=True))
    fig = plt.figure(figsize=(7.1, 5.4)); gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 0.8, 1.3], hspace=0.6, wspace=0.3)
    a_ = fig.add_subplot(gs[0, :]); t = res["t"]; m = (t >= 300) & (t <= 420)
    BB = res["Bbar"][m].T; lat = res["lat"]; vm = np.percentile(np.abs(BB), 99)
    a_.pcolormesh(t[m], lat, BB, cmap="RdBu_r", vmin=-vm, vmax=vm, shading="auto", rasterized=True)
    a_.set_ylim(-60, 60); a_.set_ylabel("latitude"); a_.set_title(f"(a) self-excited 2-D dynamo ({base}): tachocline toroidal field", loc="left", fontsize=7)
    a_ = fig.add_subplot(gs[1, :], sharex=a_)
    a_.plot(t[m], res["EN"][m], color=WONG["N"], lw=0.7, label="N emergence"); a_.plot(t[m], res["ES"][m], color=WONG["S"], lw=0.7, label="S emergence")
    a_.set_ylabel("BL flux (arb.)"); a_.legend(fontsize=5.3, ncol=2, loc="upper left"); a_.set_xlabel("model time (yr)")
    obs_stats = dict(lag1_A=load("q3_persistence.json")["lags"]["A(n),A(n+1)"]["r"], lag2_A=load("q3_persistence.json")["lags"]["A(n),A(n+2)"]["r"],
                     longest_run_mean_per_13=7.0, r_absApf_absnextA=load("q0_reconcile.json")["mag_|Apf(n)|_vs_|A(n+1)|"]["r"],
                     lag1_Apf=load("q3_persistence.json")["lags"]["Apf(n),Apf(n+1)"]["r"])
    order = ("k03", "k1", "k3", "rho05", "rho08", "k1lr", "k10lr")      # seed 1 of each; seed spread is in Q4d
    names = [k for k in order if k in S]
    a_ = fig.add_subplot(gs[2, 0])
    for j, (stat, lab, c) in enumerate((("lag1_A", "lag-1 r of A", WONG["N"]), ("lag2_A", "lag-2 r of A", WONG["S"]),
                                        ("lag1_Apf", "lag-1 r of A_pf", WONG["green"]), ("r_absApf_absnextA", "r(|A_pf(n)|, |A(n+1)|)", WONG["purple"]))):
        a_.plot(range(len(names)), [S[k]["stats"][stat] for k in names], "-o", ms=3, color=c, label=lab)
        a_.axhline(obs_stats[stat], color=c, lw=0.6, ls=":")
    a_.axhline(0, color="0.6", lw=0.4); a_.set_xticks(range(len(names))); a_.set_xticklabels([TICK.get(k, k) for k in names], fontsize=5.5)
    a_.set_ylabel("r (dotted: Sun)"); a_.legend(fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.3), ncol=2); a_.set_title("(b) memory vs coupling", loc="left", fontsize=7)
    a_ = fig.add_subplot(gs[2, 1])
    a_.bar(range(len(names)), [S[k]["stats"]["longest_run_mean_per_13"] for k in names], color=WONG["N"], width=0.6)
    for i, k in enumerate(names):
        a_.text(i, 0.3, f"{S[k]['sunlike']['frac_sunlike']*100:.0f}%", ha="center", fontsize=5, color="w")
    a_.axhline(7, color="k", lw=0.7, ls="--"); a_.text(len(names) - 0.6, 7.1, "Sun: 7", fontsize=5.5, ha="right")
    a_.axhline(4.0, color="0.5", lw=0.7, ls=":"); a_.text(len(names) - 0.6, 4.1, "independent: 4.0", fontsize=5.5, ha="right")
    a_.set_xticks(range(len(names))); a_.set_xticklabels([TICK.get(k, k) for k in names], fontsize=5.5); a_.set_ylabel("mean longest run in 13 cycles")
    a_.set_title("(c) sustained dominance (white: % Sun-like cycles)", loc="left", fontsize=7)
    savefig(fig, "Q4c_dynamo_coupling"); plt.close(fig)
tic("done")

# ------------------------------------------------------------------- Q4d dynamo seed ensembles ----
ET = os.path.join(OUT, "dynamo", "ensemble_table.json")
if os.path.exists(ET):
    tic("Q4d")
    T_ = json.load(open(ET)); labs = list(T_)
    q3_ = load("q3_persistence.json"); q0_ = load("q0_reconcile.json")
    sun = dict(lag1_A=q3_["lags"]["A(n),A(n+1)"]["r"], longest_run_mean_per_13=7.0,
               r_absApf_absnextA=q0_["mag_|Apf(n)|_vs_|A(n+1)|"]["r"], r_Apf_next_A=q3_["lags"]["Apf(n),A(n+1)"]["r"])
    TL = {"kappa=1": r"$\kappa$=1, 5 seeds", "rho=0.8": r"$\rho$=0.8, 5 seeds", "tau_c=3 yr": r"$\tau_c$=3 yr, 3 seeds", "tau_c=5 yr": r"$\tau_c$=5 yr, 3 seeds",
          "v0=10 m/s": "$v_0$=10 m/s, 3 seeds"}
    labs = [l for l in labs if l in TL]                       # the multi-seed settings
    fig, axs = plt.subplots(1, 4, figsize=(7.1, 2.7))
    for a_, key, ttl in ((axs[0], "lag1_A", "lag-1 r of A"), (axs[1], "longest_run_mean_per_13", "mean longest run / 13 cycles"),
                         (axs[2], "r_absApf_absnextA", "r(|A_pf(n)|, |A(n+1)|)"), (axs[3], "r_Apf_next_A", "r(A_pf(n), A(n+1))")):
        for i, l in enumerate(labs):
            ps = T_[l]["per_seed"][key]
            a_.plot([i, i], [ps["min"], ps["max"]], "-", color="0.55", lw=2.2, solid_capstyle="round")
            a_.plot(i, ps["median"], "o", color=WONG["N"], ms=3.5, zorder=3)
            if key in T_[l]["pooled"]:
                a_.plot(i + 0.18, T_[l]["pooled"][key], "D", color=WONG["S"], ms=3, zorder=3)
        a_.axhline(sun[key], color="k", lw=0.7, ls="--")
        if key == "longest_run_mean_per_13":
            a_.axhline(4.0, color="0.5", lw=0.6, ls=":")
        else:
            a_.axhline(0, color="0.7", lw=0.4)
        a_.set_xticks(range(len(labs))); a_.set_xticklabels([TL.get(l, l) for l in labs], rotation=45, ha="right", fontsize=5.5)
        a_.set_title(ttl, loc="left", fontsize=6.5)
    axs[0].plot([], [], "o", color=WONG["N"], ms=3.5, label="seed median"); axs[0].plot([], [], "-", color="0.55", lw=2.2, label="seed range")
    axs[0].plot([], [], "D", color=WONG["S"], ms=3, label="pooled seeds"); axs[0].plot([], [], "--", color="k", lw=0.7, label="Sun")
    axs[0].legend(fontsize=5, loc="upper right", bbox_to_anchor=(1.0, 0.9))
    fig.tight_layout(); savefig(fig, "Q4d_dynamo_seed_ensembles"); plt.close(fig)
