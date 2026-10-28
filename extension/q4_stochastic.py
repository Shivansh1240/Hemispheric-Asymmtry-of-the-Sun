"""q4_stochastic.py -- Q4 (part 1): stochastic or predetermined?

E1  Randomised emergence in the SFT, every member keeping each hemisphere's region set (sizes, hence total flux) in
    every cycle, plus a fresh tilt-scatter draw:
      R1 size shuffle    the hemisphere's region sizes are re-dealt among that hemisphere's own emergence slots
                         (time, latitude, Joy tilt of the slot): the butterfly of each hemisphere is kept, which
                         region is big is random.
      R2 slot shuffle    the (time, |lat|) slots of both hemispheres of the cycle are pooled and re-dealt: each
                         hemisphere keeps its sizes and total, but its timing and latitudes become random.
      T  tilt only       the real catalogue, tilt scatter only (reference).
    If A_pf were predetermined by the hemispheric totals, R2 would keep it.
E2  The key correlations with and without the largest regions (top k per cycle by peak area, k = 1, 2, 3, 5, 10),
    from the daily rows; and the SFT with the same regions removed.
Plain statistics only.
"""
import itertools
import numpy as np
from core import *

out = {}
obs = np.array([APF[n] for n in C1])
M = Model(); det = apf_vec(M.run(CAT), C1)
cy = CAT.cycle.values; north = CAT.latc.values > 0
groups = {(n, h): np.flatnonzero((cy == n) & (north == h)) for n in np.unique(cy) for h in (True, False)}

def r1(rng):
    c = CAT.copy(); f = c.flux.values.copy(); d = c.ddeg.values.copy(); a = c.area.values.copy()
    for idx in groups.values():
        p = rng.permutation(idx); f[idx], d[idx], a[idx] = f[p], d[p], a[p]
    c["flux"], c["ddeg"], c["area"] = f, d, a
    return c

def r2(rng):
    c = CAT.copy(); yr = c.year.values.copy(); lat = c.latc.values.copy(); al = c.alpha.values.copy()
    for n in np.unique(cy):
        idx = np.flatnonzero(cy == n)
        p = rng.permutation(idx)                       # slot of region idx[k] is p[k]
        yr[idx] = CAT.year.values[p]; al[idx] = CAT.alpha.values[p]
        lat[idx] = np.abs(CAT.latc.values[p]) * np.where(north[idx], 1.0, -1.0)
    c["year"], c["latc"], c["alpha"] = yr, lat, al
    return c.sort_values("year").reset_index(drop=True)

def tilt_for(c, rng):
    return c.alpha.values + rng.normal(0, 1, len(c)) * scatter_sigma(c.area.values)

ens = {}
for lab, maker in (("T_tilt_only", None), ("R1_size_shuffle", r1), ("R2_slot_shuffle", r2)):
    tic(lab); rng = np.random.default_rng({"T_tilt_only": 11, "R1_size_shuffle": 12, "R2_slot_shuffle": 13}[lab])
    E = []
    for k in range(100):
        c = CAT if maker is None else maker(rng)
        E.append(apf_vec(M.run(c, tilt_for(c, rng)), C1))
    E = np.array(E); ens[lab] = E
    rr = [r(e, obs) for e in E]; rd = [r(e, det) for e in E]
    sig_frac = float(np.var(E.mean(0)) / (np.var(E.mean(0)) + E.var(0).mean()))
    out[lab] = dict(r_obs_median=float(np.median(rr)), r_obs_16_84=np.percentile(rr, [16, 84]).tolist(),
                    r_ensmean_obs=r(E.mean(0), obs), r_det_median=float(np.median(rd)),
                    same_sign_as_obs=dict(zip(C1, (np.sign(E) == np.sign(obs)).mean(0))),
                    same_sign_as_det=dict(zip(C1, (np.sign(E) == np.sign(det)).mean(0))),
                    cycles_sign_kept_84pct_obs=int(((np.sign(E) == np.sign(obs)).mean(0) >= 0.84).sum()),
                    signal_fraction=sig_frac, rms_member=float(np.median(np.sqrt((E ** 2).mean(1)))),
                    p16=dict(zip(C1, np.percentile(E, 16, 0))), p84=dict(zip(C1, np.percentile(E, 84, 0))),
                    median=dict(zip(C1, np.median(E, 0))))
    v = out[lab]
    print(f"  r(member, obs) {v['r_obs_median']:+.2f} [{v['r_obs_16_84'][0]:+.2f}, {v['r_obs_16_84'][1]:+.2f}]  r(ens mean, obs) {v['r_ensmean_obs']:+.2f}"
          f"  r(member, real-catalogue model) {v['r_det_median']:+.2f}  signal fraction {v['signal_fraction']:.2f}  cycles with obs sign in >=84 %: {v['cycles_sign_kept_84pct_obs']}/10")
np.savez(os.path.join(OUT, "q4_randomised.npz"), det=det, obs=obs, **ens)
out["det_model"] = dict(r_obs=r(det, obs), n=10)

tic("E2: key correlations without the largest regions")
rank = D.groupby(["cyc", "gid"]).rmax.first().reset_index()
rank["rk"] = rank.groupby("cyc").rmax.rank(ascending=False, method="first")
an1 = np.array([AN[n + 1] for n in C2])
def keys(mask):
    ann = ann_of(mask)
    A_ = {n: cycA(n, mask) for n in range(12, 26)}
    seq = [np.sign(A_[n]) for n in range(12, 25)]
    return dict(H1_x1hat_apf=r(x1v(C1, 0, 1, True, ann), obs),
                H2_x1_next=r(x1v(C2, -1, 4, False, ann), [A_[n + 1] for n in C2]),
                H2_x1_next_fixed_target=r(x1v(C2, -1, 4, False, ann), an1),
                decline_unsigned_apf=r(x1v(C1, 3, 4, False, ann), obs),
                decline_signs=signs(x1v(C1, 3, 4, False, ann), obs),
                persistence_lag1=r([A_[n] for n in C2], [A_[n + 1] for n in C2]),
                runs=[len(list(g)) for _, g in itertools.groupby(seq)],
                A21=A_[21], A14=A_[14])
e2 = {"all": keys(None)}
for k in (1, 2, 3, 5, 10, 20):
    big = set(rank.gid[rank.rk <= k]); m = ~D.gid.isin(big)
    e2[f"without_top{k}"] = keys(m)
    e2[f"only_top{k}"] = keys(~m) if k >= 5 else None
    e2[f"without_top{k}"]["regions_removed"] = len(big)
    e2[f"without_top{k}"]["area_share_removed"] = float(D.area[~m].sum() / D.area.sum())
for k_, v in e2.items():
    if v:
        print(f"  {k_:14s} H1 {v['H1_x1hat_apf']:+.2f}  H2 {v['H2_x1_next']:+.2f} (fixed target {v['H2_x1_next_fixed_target']:+.2f})  decline {v['decline_unsigned_apf']:+.2f} ({v['decline_signs']}/10)"
              f"  lag1 {v['persistence_lag1']:+.2f}  runs {v['runs']}  A21 {v['A21']:+.3f} A14 {v['A14']:+.3f}"
              + (f"  removed {v['regions_removed']} regions, {v['area_share_removed']*100:.1f} % of area" if "regions_removed" in v else ""))
out["E2_data"] = e2

tic("E2: SFT without the largest regions")
crank = CAT.groupby("cycle").area.rank(ascending=False, method="first").values
sft_e2 = {"all": dict(r_obs=r(det, obs), signs=signs(det, obs))}
for k in (1, 3, 5, 10, 20):
    keep = crank > k
    a = apf_vec(M.run(CAT[keep].reset_index(drop=True)), C1)
    a_only = apf_vec(M.run(CAT[~keep].reset_index(drop=True), ic=False), C1)
    sft_e2[f"without_top{k}"] = dict(r_obs=r(a, obs), signs=signs(a, obs), r_with_full_model=r(a, det), D21_sign=int(np.sign(a[C1.index(21)])))
    sft_e2[f"only_top{k}"] = dict(r_obs=r(a_only, obs), signs=signs(a_only, obs), r_with_full_model=r(a_only, det))
    print(f"  SFT without top{k}: r obs {sft_e2[f'without_top{k}']['r_obs']:+.2f} (signs {sft_e2[f'without_top{k}']['signs']}/10), r with full model {sft_e2[f'without_top{k}']['r_with_full_model']:+.2f};"
          f"  only top{k}: r obs {sft_e2[f'only_top{k}']['r_obs']:+.2f}, r with full model {sft_e2[f'only_top{k}']['r_with_full_model']:+.2f}")
out["E2_sft"] = sft_e2
dump("q4_stochastic.json", out)
tic("done")
