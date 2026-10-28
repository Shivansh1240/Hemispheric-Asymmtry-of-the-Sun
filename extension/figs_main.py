"""figs_main.py -- the seven requested plot types (P1-P7) and the Q2 'whole table'.

P1 butterfly per hemisphere (marker size = peak area) with the polar-field trace and the N-S differences underneath
P2 r over window start x end (data, both targets) and the model's own r(D_window, D)
P3 Green's function: emergence epoch x latitude -> contribution to the N pole, the S pole and the imbalance
P4 cumulative N-S area through each cycle (all / low-latitude large / rest)
P5 return map A(n) -> A(n+1) and A(n) -> A(n+2) with arrows
P6 concentration curve (top-k regions' share of the modelled imbalance)  [data from explore.py]
P7 stacked bars: each pole's field by size class of the contributing regions, per cycle
Q2 size class x latitude band table (r at the paper windows and at the physics window)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from core import *
import lib_asym as L                       # importable once core/asymcore has set the project path

style()
q2 = load("q2_window.json")
PHYS = (q2["physics_window"]["s"], q2["physics_window"]["e"])
MINS = [SCB[n] for n in range(12, 27)]

def cycle_ticks(ax, top=True):
    for t in MINS:
        ax.axvline(t, color="0.85", lw=0.5, zorder=0)
    if top:
        for n in range(12, 26):
            ax.text(0.5 * (SCB[n] + SCB[n + 1]), 1.02, str(n), transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6, color="0.35")

# ---------------------------------------------------------------- P1 butterfly + polar trace + N-S ----
tic("P1")
fig, ax = plt.subplots(4, 1, figsize=(7.1, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[2.3, 1.1, 1.0, 1.0], hspace=0.12))
o = np.argsort(CAT.area.values)
c = CAT.iloc[o]; col = np.where(c.latc > 0, WONG["N"], WONG["S"])
ax[0].scatter(c.year, c.latc, s=np.clip(c.area.values, 5, None) / 60.0, c=col, alpha=0.45, lw=0, rasterized=True)
ax[0].axhline(0, color="0.5", lw=0.4); ax[0].set_ylim(-50, 50); ax[0].set_ylabel("latitude (deg)")
for a_ in (500, 1500, 3000):
    ax[0].scatter([], [], s=a_ / 60.0, c="0.4", lw=0, label=f"{a_} MSH")
ax[0].legend(loc="upper left", ncol=3, fontsize=6, handletextpad=0.2, borderaxespad=0.2)
ax[0].text(0.995, 0.96, "North", color=WONG["N"], transform=ax[0].transAxes, ha="right", va="top", fontsize=7)
ax[0].text(0.995, 0.04, "South", color=WONG["S"], transform=ax[0].transAxes, ha="right", va="bottom", fontsize=7)
cycle_ticks(ax[0])
pf = L._csv_pf(); w = L._wso_daily()
for h, cc in (("N", WONG["N"]), ("S", WONG["S"])):
    t, b = pf[(h, "MWO")]; m = t < 1976.4
    ax[1].plot(t[m], b[m], "-o", ms=1.5, lw=0.8, color=cc, label=f"{h} (MWO faculae, calibrated)" if h == "N" else None)
    ax[1].plot(w[:, 0], w[:, 1 if h == "N" else 2] / 100.0, "-", lw=0.6, color=cc, alpha=0.8)
ax[1].axhline(0, color="0.5", lw=0.4); ax[1].set_ylabel("polar field (G)")
ax[1].text(0.005, 0.95, "B_N (blue), B_S (orange): MWO faculae to 1976, WSO after", transform=ax[1].transAxes, fontsize=6, va="top")
cycle_ticks(ax[1], False)
yrs = ANN.yr.values; dA = (ANN.AN - ANN.AS).values
ax[2].bar(yrs + 0.5, dA / 1e5, width=0.9, color=np.where(dA > 0, WONG["N"], WONG["S"]), lw=0)
ax[2].axhline(0, color="0.5", lw=0.4); ax[2].set_ylabel("A_N - A_S\n(1e5 MSH day)")
cycle_ticks(ax[2], False)
# |B_N| - |B_S|: yearly means, MWO before 1976.4, WSO after
yy = np.arange(1906, 2025); dB = []
for y in yy:
    if y + 1 <= 1976.4:
        tN, bN = pf[("N", "MWO")]; tS, bS = pf[("S", "MWO")]
        a = bN[(tN >= y) & (tN < y + 1)]; b = bS[(tS >= y) & (tS < y + 1)]
        dB.append(abs(a.mean()) - abs(b.mean()) if len(a) and len(b) else np.nan)
    else:
        m = (w[:, 0] >= y) & (w[:, 0] < y + 1)
        dB.append((abs(w[m, 1].mean()) - abs(w[m, 2].mean())) / 100.0 if m.sum() > 5 else np.nan)
dB = np.array(dB)
ax[3].bar(yy + 0.5, dB, width=0.9, color=np.where(dB > 0, WONG["N"], WONG["S"]), lw=0)
ax[3].axhline(0, color="0.5", lw=0.4); ax[3].set_ylabel("|B_N| - |B_S| (G)"); ax[3].set_xlabel("year")
for n in C1:
    ax[3].plot(SCB[n + 1], APF[n] * 2.0, "kD", ms=2.5)
ax[3].text(0.005, 0.95, "diamonds: A_pf at each minimum (x2, same sign convention)", transform=ax[3].transAxes, fontsize=6, va="top")
cycle_ticks(ax[3], False)
ax[3].set_xlim(1874, 2027)
savefig(fig, "P1_butterfly_polar_NS"); plt.close(fig)

# ---------------------------------------------------------------- P2 window heatmaps ----
tic("P2")
S_ = list(range(-5, 8)); E_ = list(range(-5, 9))
def grid(key, src):
    Z = np.full((len(E_), len(S_)), np.nan)
    for k, v in src.items():
        s, e = map(int, k.split(",")); Z[E_.index(e), S_.index(s)] = v if key is None else v[key]
    return Z
Za = grid("r_apf_n10", q2["scan_window"]); Zb = grid("r_next_n12", q2["scan_window"]); Zc = grid(None, q2["physics_window_r"])
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.9), sharey=True)
for a_, Z, ttl in ((axs[0], Za, "(a) data: X1hat vs A_pf (n = 10)"), (axs[1], Zb, "(b) data: X1 vs A(n+1) (n = 12)"),
                   (axs[2], Zc, "(c) SFT: D_window vs D (n = 13)")):
    im = a_.imshow(Z, origin="lower", cmap="RdBu_r", vmin=-1, vmax=1, extent=(-5.5, 7.5, -5.5, 8.5), aspect="auto")
    for i, e in enumerate(E_):
        for j, s in enumerate(S_):
            if np.isfinite(Z[i, j]):
                a_.text(s, e, f"{Z[i, j]:.2f}".replace("0.", "."), ha="center", va="center", fontsize=3.6,
                        color="k" if abs(Z[i, j]) < 0.75 else "w")
    a_.plot(0, 1, "s", mfc="none", mec="k", ms=7, mew=0.9); a_.plot(-1, 4, "o", mfc="none", mec="k", ms=7, mew=0.9)
    a_.plot(*PHYS, "D", mfc="none", mec=WONG["green"], ms=7, mew=1.2)
    a_.set_title(ttl, loc="left"); a_.set_xlabel("window start s (yr after t_R)")
axs[0].set_ylabel("window end e (yr after t_R)")
cb = fig.colorbar(im, ax=axs, shrink=0.85, pad=0.01); cb.set_label("Pearson r")
fig.text(0.01, -0.04, "square: t_R+0->+1 (paper, channel 1); circle: t_R-1->+4 (paper, channel 2); green diamond: physics window "
         f"t_R+{PHYS[0]}->+{PHYS[1]} (shortest window carrying r >= 0.9 of the model's own imbalance)", fontsize=6)
savefig(fig, "P2_window_heatmaps"); plt.close(fig)

# ---------------------------------------------------------------- P3 Green's-function maps ----
tic("P3")
G3 = np.load(os.path.join(OUT, "q2_green_maps.npz"))
E, Lb = G3["E"], G3["L"]
fig, axs = plt.subplots(1, 3, figsize=(7.1, 2.6), sharey=True)
for a_, key, ttl in ((axs[0], "HN", "(a) builds the N pole"), (axs[1], "HS", "(b) builds the S pole"), (axs[2], "HD", "(c) builds the imbalance |B_N|-|B_S|")):
    Z = G3[key].T * 100; vm = np.nanmax(np.abs(Z))
    im = a_.pcolormesh(E, Lb, Z, cmap="RdBu_r", norm=TwoSlopeNorm(0, -vm, vm), shading="flat", rasterized=True)
    a_.axvline(0, color="k", lw=0.5, ls=":"); a_.axhline(0, color="0.5", lw=0.4)
    a_.axvspan(PHYS[0], PHYS[1] + 1, color=WONG["green"], alpha=0.12, lw=0)
    a_.set_title(ttl, loc="left"); a_.set_xlabel("emergence epoch - t_R (yr)"); a_.set_xlim(-8, 10)
    cb = fig.colorbar(im, ax=a_, shrink=0.8, pad=0.02); cb.ax.tick_params(labelsize=5); cb.set_label("% of total per cell", fontsize=6)
axs[0].set_ylabel("emergence latitude (deg)")
q = q2["imbalance_epoch_quantiles"]
axs[2].text(0.02, 0.03, f"imbalance built at t_R+{q['0.16']:.1f} to +{q['0.84']:.1f} yr (16-84 %)", transform=axs[2].transAxes, fontsize=5.5)
savefig(fig, "P3_green_epoch_latitude"); plt.close(fig)

# ---------------------------------------------------------------- P4 cumulative N-S ----
tic("P4")
LL = band_mask(D, 500, np.inf, 0, 20)
fig, axs = plt.subplots(2, 7, figsize=(7.1, 2.9), sharex=True, sharey=True)
north = (D.hemi == "N").values
for k, n in enumerate(range(12, 26)):
    a_ = axs.flat[k]; m = ((D.dec >= SCB[n]) & (D.dec < SCB[n + 1])).values
    t = D.dec.values[m] - SCB[n]; o = np.argsort(t); tot = D.area.values[m].sum()
    for msk, cc, lab in ((np.ones(m.sum(), bool), "k", "all"), (LL.values[m], WONG["giant"], "low-lat >= 500 MSH"), (~LL.values[m], WONG["rest"], "rest")):
        v = np.where(north[m], 1.0, -1.0) * D.area.values[m] * msk
        a_.plot(t[o], np.cumsum(v[o]) / tot, color=cc, lw=0.9 if cc == "k" else 0.8, label=lab)
    a_.axhline(0, color="0.6", lw=0.4)
    flag = " flip" if n in (14, 21) else (" near" if n in (17, 18, 24) else "")
    a_.set_title(f"SC{n}{flag}" + (" (prov.)" if n == 25 else ""), fontsize=6.5, loc="left")
for a_ in axs[1]:
    a_.set_xlabel("yr since min", fontsize=6)
axs[0, 0].set_ylabel("cum. (N-S)/total", fontsize=6); axs[1, 0].set_ylabel("cum. (N-S)/total", fontsize=6)
axs[0, 0].legend(fontsize=5, loc="lower left", handlelength=1.2)
savefig(fig, "P4_cumulative_NS"); plt.close(fig)

# ---------------------------------------------------------------- P5 return map ----
tic("P5")
fig, axs = plt.subplots(1, 2, figsize=(7.1, 3.3))
for a_, lagk, ttl in ((axs[0], 1, "(a) A(n) -> A(n+1)"), (axs[1], 2, "(b) A(n) -> A(n+2)")):
    ns = [n for n in range(12, 26 - lagk)]
    x = np.array([AN[n] for n in ns]); y = np.array([AN[n + lagk] for n in ns])
    a_.axhline(0, color="0.6", lw=0.5); a_.axvline(0, color="0.6", lw=0.5)
    a_.fill_between([0, 0.3], 0, 0.3, color=WONG["N"], alpha=0.06, lw=0); a_.fill_between([-0.3, 0], -0.3, 0, color=WONG["S"], alpha=0.06, lw=0)
    for i in range(len(ns) - 1):
        a_.annotate("", xy=(x[i + 1], y[i + 1]), xytext=(x[i], y[i]), arrowprops=dict(arrowstyle="->", lw=0.6, color="0.45", shrinkA=3, shrinkB=3))
    prov = np.array([n + lagk == 25 for n in ns])
    a_.scatter(x[~prov], y[~prov], s=22, c=np.where(y[~prov] > 0, WONG["N"], WONG["S"]), zorder=3, edgecolor="w", lw=0.5)
    a_.scatter(x[prov], y[prov], s=22, facecolor="none", edgecolor="k", zorder=3, lw=0.6)
    for n, xi, yi in zip(ns, x, y):
        a_.text(xi + 0.006, yi + 0.006, str(n), fontsize=6)
    a_.plot([-0.3, 0.3], [-0.3, 0.3], ":", color="0.6", lw=0.6)
    lim = 0.21; a_.set_xlim(-lim, lim); a_.set_ylim(-lim, lim); a_.set_aspect("equal")
    full = [n for n in ns if n + lagk <= 24]
    a_.set_title(f"{ttl}: r = {r([AN[n] for n in full], [AN[n + lagk] for n in full]):+.2f} (n = {len(full)})", loc="left")
    a_.set_xlabel("A(n)"); a_.set_ylabel(f"A(n+{lagk})")
fig.text(0.01, -0.02, "Arrows join successive cycles (label = n). Shaded quadrants: same sign as the earlier cycle (dominance kept). Open circle: provisional SC25.", fontsize=6)
savefig(fig, "P5_return_map"); plt.close(fig)

# ---------------------------------------------------------------- P7 stacked bars by size class ----
tic("P7")
M = Model(); res = M.run(CAT, keep=True)
CLS = [("<100", 0, 100, "#c6dbef"), ("100-499", 100, 500, "#6baed6"), ("500-999", 500, 1000, "#2171b5"), (">=1000", 1000, np.inf, "#08306b")]
fig, ax = plt.subplots(figsize=(7.1, 2.8))
xs = np.arange(len(CYC)); wdt = 0.38
for p, (pole, off) in enumerate((("N", -0.2), ("S", 0.2))):
    posb = np.zeros(len(CYC)); negb = np.zeros(len(CYC))
    parts = []
    for n in CYC:
        sel, c = res["reg"][n]; full = res["full"][n]; sgn = np.sign(full[p]); v = sgn * c[:, p]
        own = CAT.cycle.values[sel] == n; tot = abs(full[p])
        row = [v[own & (CAT.area.values[sel] >= lo) & (CAT.area.values[sel] < hi)].sum() / tot for _, lo, hi, _ in CLS]
        row.append(v[~own].sum() / tot); row.append(sgn * M.IC[n][p] / tot)
        parts.append(row)
    parts = np.array(parts)
    labels = [f"{c_[0]} MSH" for c_ in CLS] + ["older cycles", "initial field"]
    colors = [c_[3] for c_ in CLS] + ["#bdbdbd", "#f0f0f0"]
    for j in range(parts.shape[1]):
        v = parts[:, j]; bot = np.where(v >= 0, posb, negb)
        ax.bar(xs + off, v, wdt, bottom=bot, color=colors[j], edgecolor="w", lw=0.4, label=labels[j] if p == 0 else None)
        posb += np.where(v >= 0, v, 0); negb += np.where(v < 0, v, 0)
    for i in range(len(CYC)):
        ax.text(xs[i] + off, posb[i] + 0.015, pole, ha="center", va="bottom", fontsize=5, color=WONG[pole])
ax.axhline(0, color="k", lw=0.5); ax.axhline(1, color="0.6", lw=0.5, ls=":")
ax.set_xticks(xs); ax.set_xticklabels([f"SC{n}" for n in CYC]); ax.set_ylabel("share of the pole's field at minimum")
ax.legend(ncol=6, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.13))
savefig(fig, "P7_pole_flux_by_size"); plt.close(fig)
shares = {}
for p, pole in enumerate("NS"):
    for n in CYC:
        sel, c = res["reg"][n]; full = res["full"][n]; v = np.sign(full[p]) * c[:, p]; own = CAT.cycle.values[sel] == n
        shares[f"{pole}{n}"] = {lab: float(v[own & (CAT.area.values[sel] >= lo) & (CAT.area.values[sel] < hi)].sum() / abs(full[p])) for lab, lo, hi, _ in CLS}
dump("p7_pole_shares.json", shares)

# ---------------------------------------------------------------- Q2 whole size x latitude table ----
tic("Q2 table")
tab = q2["scan_size_lat"]; sizes = list(dict.fromkeys(t["size"] for t in tab)); lats = list(dict.fromkeys(t["lat"] for t in tab))
# cycles with a non-zero predictor in the window (degenerate cells are greyed)
nz = {}
for t in tab:
    lo, hi = {"all": (0, np.inf), "<100": (0, 100), "100-499": (100, 500), "500-999": (500, 1000), ">=1000": (1000, np.inf), ">=500": (500, np.inf)}[t["size"]]
    a0, a1 = {"all": (0, 90), "0-10": (0, 10), "10-20": (10, 20), "20-30": (20, 30), ">=30": (30, 90), "0-20": (0, 20)}[t["lat"]]
    an = ann_of(band_mask(D, lo, hi, a0, a1))
    for lab, (s, e), cyc in (("apf_w01", (0, 1), C1), ("next_wm14", (-1, 4), C2), ("apf_phys", PHYS, C1), ("next_phys", PHYS, C2)):
        nz[(t["size"], t["lat"], lab)] = int(np.sum(np.abs(x1v(cyc, s, e, False, an)) > 0))
fig, axs = plt.subplots(2, 2, figsize=(7.1, 5.2))
for a_, key, ttl, ncyc in ((axs[0, 0], "apf_w01", "(a) X1hat vs A_pf, t_R+0->+1 (n = 10)", 10), (axs[0, 1], "next_wm14", "(b) X1 vs A(n+1), t_R-1->+4 (n = 12)", 12),
                           (axs[1, 0], "apf_phys", f"(c) X1hat vs A_pf, physics window t_R+{PHYS[0]}->+{PHYS[1]}", 10),
                           (axs[1, 1], "next_phys", f"(d) X1 vs A(n+1), physics window t_R+{PHYS[0]}->+{PHYS[1]}", 12)):
    Z = np.full((len(sizes), len(lats)), np.nan); Nr = np.zeros_like(Z)
    for t in tab:
        i, j = sizes.index(t["size"]), lats.index(t["lat"])
        ok = nz[(t["size"], t["lat"], key)] >= ncyc - 2
        Z[i, j] = t[f"r_{key}"] if ok else np.nan; Nr[i, j] = t["n_regions"]
    a_.imshow(Z, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    for i in range(len(sizes)):
        for j in range(len(lats)):
            txt = f"{Z[i, j]:+.2f}" if np.isfinite(Z[i, j]) else "--"
            a_.text(j, i - 0.12, txt, ha="center", va="center", fontsize=6.5, color="w" if np.isfinite(Z[i, j]) and abs(Z[i, j]) > 0.75 else "k")
            a_.text(j, i + 0.25, f"{int(Nr[i, j])}", ha="center", va="center", fontsize=4.5, color="0.4" if not (np.isfinite(Z[i, j]) and abs(Z[i, j]) > 0.75) else "0.85")
    a_.set_xticks(range(len(lats))); a_.set_xticklabels([f"{l}" for l in lats], fontsize=6); a_.set_yticks(range(len(sizes))); a_.set_yticklabels(sizes, fontsize=6)
    a_.set_title(ttl, loc="left", fontsize=7); a_.set_xlabel("|latitude| at peak (deg)", fontsize=6)
    a_.spines[:].set_visible(False)
for a_ in axs[:, 0]:
    a_.set_ylabel("peak area (MSH)", fontsize=6)
fig.text(0.01, 0.0, "Every cell is shown. Small grey numbers: regions in the cell. '--': the predictor is zero in more than two cycles of the window (r undefined or degenerate).", fontsize=6)
fig.tight_layout()
savefig(fig, "Q2_size_latitude_table"); plt.close(fig)
dump("q2_table_nonzero.json", {f"{a}|{b}|{c}": v for (a, b, c), v in nz.items()})
tic("done")
