"""q2b_phase_sign.py -- the effective (parity) sign, modelled explicitly with the physics window.

The SFT says the polar imbalance at a minimum is built by flux emerging t_R+2.3 ... +5.1 yr (q2_window.py), and the
model is polarity-blind, so its prediction for A_pf is the UNSIGNED area imbalance of that epoch.  This script
compares, for every cycle 12-24:
    X1 at maximum (t_R+0 -> +1)   and   X1 in the decline (physics window t_R+3 -> +4, chosen from the model only)
and asks whether the dominant hemisphere switches between the two, odd vs even.  Then A_pf against the unsigned
decline-phase X1 (no parity factor).  Window sensitivity is reported for neighbouring windows.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from core import *

style()
q2 = load("q2_window.json"); PS, PE = q2["physics_window"]["s"], q2["physics_window"]["e"]
out = {}
cyc = list(range(12, 25)); apf = np.array([APF[n] for n in C1])
xm = {n: x1(n, 0, 1) for n in cyc}; xd = {n: x1(n, PS, PE) for n in cyc}
sw = {n: bool(np.sign(xm[n]) != np.sign(xd[n])) for n in cyc}
ev = [n for n in cyc if n % 2 == 0]; od = [n for n in cyc if n % 2]
out["switch_max_vs_decline"] = dict(per_cycle=sw, even=f"{sum(sw[n] for n in ev)}/{len(ev)}", odd=f"{sum(sw[n] for n in od)}/{len(od)}",
                                    sc12_14_only=dict(even=f"{sum(sw[n] for n in (12, 14))}/2", odd=f"{int(sw[13])}/1"))
xd10 = np.array([xd[n] for n in C1]); xh = x1v(C1, 0, 1, True)
odd = np.array([n % 2 == 1 for n in C1])
out["Apf_vs_unsigned_decline_X1"] = dict(r=r(xd10, apf), n=10, signs=signs(xd10, apf), r_odd=r(xd10[odd], apf[odd]), r_even=r(xd10[~odd], apf[~odd]),
                                         failing_cycles=[n for n, a, b in zip(C1, xd10, apf) if np.sign(a) != np.sign(b)],
                                         r_without_SC19=r(np.delete(xd10, 4), np.delete(apf, 4)))
out["Apf_vs_parity_signed_max_X1hat"] = dict(r=r(xh, apf), n=10, signs=signs(xh, apf), r_without_SC19=r(np.delete(xh, 4), np.delete(apf, 4)))
sens = {}
for (a, b, c, d) in [(0, 1, 3, 4), (0, 1, 3, 5), (0, 1, 4, 5), (0, 1, 2, 4), (0, 1, 2, 3), (-1, 1, 3, 5), (0, 2, 3, 5)]:
    s_ = {n: np.sign(x1(n, a, b)) != np.sign(x1(n, c, d)) for n in cyc}; xdd = x1v(C1, c, d)
    sens[f"max[{a},{b}] vs decline[{c},{d}]"] = dict(even=f"{sum(s_[n] for n in ev)}/7", odd=f"{sum(s_[n] for n in od)}/6",
                                                     r_unsigned_decline_vs_Apf=r(xdd, apf), signs=signs(xdd, apf))
out["window_sensitivity"] = sens
# link to exploratory E: the N-S timing lag (area centroid) against the model's imbalance and the decline-phase X1
ex = load("explore.json"); lag = np.array([ex["E_timing"]["per_cycle"][str(n)]["lag_N_minus_S"] for n in cyc])
mod = apf_vec(Model().run(CAT), cyc)
out["timing_lag_links"] = dict(r_lag_vs_model_Apf=r(lag, mod), n=len(cyc), r_lag_vs_X1_decline=r(lag, [xd[n] for n in cyc]),
                               r_lag_vs_A_next=r(lag[:-1], [AN[n + 1] for n in cyc[:-1]]), n_next=len(cyc) - 1,
                               signs_lag_vs_A_next=signs(lag[:-1], [AN[n + 1] for n in cyc[:-1]]))
for k, v in out.items():
    print(k, v)
dump("q2b_phase_sign.json", out)

fig, axs = plt.subplots(1, 2, figsize=(7.1, 2.7), gridspec_kw=dict(width_ratios=[1.6, 1]))
a_ = axs[0]; xs = np.arange(len(cyc))
mm = np.array([xm[n] for n in cyc]) / 1e4; dd = np.array([xd[n] for n in cyc]) / 1e4
a_.bar(xs - 0.2, mm, 0.38, color="0.65", label="X1 at maximum, t_R+0->+1")
a_.bar(xs + 0.2, dd, 0.38, color=[WONG["N"] if v > 0 else WONG["S"] for v in dd], label=f"X1 in decline, t_R+{PS}->+{PE} (physics window)")
for i, n in enumerate(cyc):
    if sw[n]:
        a_.text(i, max(mm[i], dd[i], 0) + 0.15, "switch", ha="center", fontsize=5.5)
    a_.text(i, -3.9, "even" if n % 2 == 0 else "odd", ha="center", fontsize=5.5, color="0.35")
a_.axhline(0, color="0.5", lw=0.5); a_.set_xticks(xs); a_.set_xticklabels([str(n) for n in cyc]); a_.set_ylim(-4.1, 4.1)
a_.set_ylabel("X1 (1e4)"); a_.set_xlabel("cycle"); a_.legend(fontsize=5.5, loc="upper left")
a_.set_title(f"(a) dominant hemisphere switches maximum -> decline: even {out['switch_max_vs_decline']['even']}, odd {out['switch_max_vs_decline']['odd']}", loc="left", fontsize=7)
a_ = axs[1]
for n, x_, y_ in zip(C1, xd10 / 1e4, apf):
    a_.plot(x_, y_, "o" if n % 2 else "s", color=WONG["N"] if y_ > 0 else WONG["S"], ms=4)
    a_.annotate(str(n), (x_, y_), xytext=(2, 2), textcoords="offset points", fontsize=5.5)
b = np.polyfit(xd10 / 1e4, apf, 1); xx = np.linspace(xd10.min() / 1e4, xd10.max() / 1e4, 5); a_.plot(xx, np.polyval(b, xx), "-", color="0.3", lw=0.7)
a_.axhline(0, color="0.7", lw=0.4); a_.axvline(0, color="0.7", lw=0.4)
v = out["Apf_vs_unsigned_decline_X1"]
a_.set_title(f"(b) unsigned, no parity factor: r = {v['r']:+.2f} (n = 10), signs {v['signs']}/10", loc="left", fontsize=7)
a_.set_xlabel(f"X1, t_R+{PS}->+{PE} (1e4)"); a_.set_ylabel("A_pf")
fig.tight_layout(); savefig(fig, "Q2b_phase_sign"); plt.close(fig)
