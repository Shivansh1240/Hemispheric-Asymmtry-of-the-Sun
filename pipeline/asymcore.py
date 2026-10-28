"""asymcore.py -- the reusable core of the hemispheric-asymmetry analysis, in one import.

For new work:  `import asymcore as A`  and build on the functions below.  They call the verified project modules
(lib_asym, common, p4lib, x_engine, bn_sft); nothing that produced a published number is re-implemented here, so a
result built with asymcore uses exactly the same data, conventions and engine as the paper.

    python asymcore.py          self-test: data checksum, the two headline correlations, Table 1, one SFT readout
    python asymcore.py --sft    also builds the SFT Green's function and runs one deterministic readout (~1-2 min)

Sections
  1 conventions     SCB, TR, C1, C2, parity(n)
  2 data            daily(), sft_catalogue(), annual(), cycle_asymmetry(), size_mask()
  3 predictor       x1(), x1_vector()
  4 polar field     apf(), apf_table(), apf_sigma(), polar_signed(), apf_and_S()
  5 targets         targets()  -> A_pf, S, A^n, A^(n+1) per cycle, in one dict
  6 statistics      pearson, spearman, corr, fisher_ci, boot_ci, runs_test, ols, loo_pred, nested_select, skill,
                    controls, perm_max_r, nested_window_selection
  7 SFT             SFTModel (Green's-function readout), HemiSFT (N/S-dependent transport), ensemble(), ensemble_apf()
  8 plotting        theme(), COLORS, save()
  9 checks          headline(), selftest()

Conventions that must not change (they were each the source of a real error; see pipeline/README.md):
  * parity factor (-1)^(n+1), applied ONLY to polar-field targets; area targets are unsigned.
  * the factor is an empirical regularity (= -sgn B_N at every minimum), NOT a Hale correction of A_pf.
  * A_pf(n) = (|B_N|-|B_S|)/(|B_N|+|B_S|), +-1.6 yr about the minimum that ENDS cycle n;
    MWO faculae for SC15-19, WSO daily filtered field for SC20-24.
  * size classes are per region (peak area), which needs the unique region ids of RGO_SOON_corrected.txt.
  * rows at Lat = 0 are dropped from hemispheric sums.
"""
import os, sys
import numpy as np

# ------------------------------------------------------------------------------------------- 0  paths ----
# the 'Asymmetry analysis' folder: $ASYM_PROJECT, else this repository if it holds paper/, else the original location
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROJECT = os.environ.get("ASYM_PROJECT") or (_REPO if os.path.isdir(os.path.join(_REPO, "paper")) else
          r"C:\Users\Shivansh Mishra\Downloads\Quantum Entanglement and Computation - Online LaTeX Editor Overleaf_files"
          r"\Asymmetry analysis")
PAPER = os.path.join(PROJECT, "paper")
CR = os.path.join(PAPER, "corrected_reanalysis"); REV = os.path.join(PAPER, "revision"); P4DIR = os.path.join(REV, "p4")
for _p in (P4DIR, REV, CR, os.path.join(CR, "engine")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
if not os.path.isdir(CR):
    raise ImportError(f"project not found under {PROJECT}; set ASYM_PROJECT to the 'Asymmetry analysis' folder")

import lib_asym as L                                   # noqa: E402  data layer, predictor, A_pf, statistics
import common as C                                     # noqa: E402  polar reductions, epochs, calibrations, CIs

def _p4():                                             # p4lib imports the SFT engine; load it only when needed
    import p4lib
    return p4lib

def _x():
    import x_engine
    return x_engine

# ------------------------------------------------------------------------------------ 1  conventions ----
SCB, TR = L.SCB, L.TR                                  # cycle minima, polar-reversal epochs (t_R)
C1 = list(range(15, 25))                               # cycles with a measured A_pf (channel 1, n = 10)
C2 = list(range(12, 24))                               # complete next-cycle pairs (channel 2, n = 12)
RGO = L.RGO_NEW                                        # the single sunspot source
OUT_REANALYSIS, OUT_REVISION, OUT_P4 = L.OUT, C.OUT, os.path.join(P4DIR, "out")

def parity(n):
    """(-1)^(n+1): +1 for odd cycles.  Polar-field targets only."""
    return (-1) ** (n + 1)

# ------------------------------------------------------------------------------------------- 2  data ----
def daily(equator="drop"):
    """Every daily row: yr, dec, hemi, area (MSH), abslat, lat, gid (unique region id), cyc, rmax (peak area of the
    region), rlat (|lat| of the region at its peak).  equator='drop' is the corrected convention."""
    return L.daily(RGO, equator=equator)

def sft_catalogue():
    """One record per region at maximum flux, as injected into the SFT (year, latc, area, flux, alpha, ddeg, parity,
    cycle, giant, trel ...)."""
    return _x().load_catalogue()

def annual(d=None, mask=None):
    """Annual table yr, AN, AS, ml (mean |lat|), X1y = (AN-AS)/ml, from all rows or rows where mask is True."""
    d = daily() if d is None else d
    return L.ann_of(d if mask is None else d[mask])

def cycle_asymmetry(n, mask=None, d=None):
    """Whole-cycle (minimum to minimum) area asymmetry (N-S)/(N+S) of cycle n, optionally of a sub-population."""
    d = daily() if d is None else d
    return L.AA(d, n, mask)

def size_mask(d, lo=0.0, hi=np.inf, latmax=None, by="peak"):
    """Rows of regions with lo <= area < hi (by='peak': the region's peak area; by='row': the daily row area, legacy)
    and, if latmax is given, |lat| at peak <= latmax."""
    a = d.rmax if by == "peak" else d.area
    m = (a >= lo) & (a < hi)
    if latmax is not None:
        m &= d.rlat <= latmax
    return m

# -------------------------------------------------------------------------------------- 3  predictor ----
def x1(n, s=0, e=1, signed=False, ann=None, normalised=False):
    """X1 = sum over whole years y with y+0.5 in [t_R+s, t_R+e+1) of (A_N-A_S)_y / mean|lat|_y.
    signed=True multiplies by parity(n) (X1hat, polar-field targets).  normalised divides by sum (A_N+A_S)_y/ml_y."""
    if normalised:
        return _p4().xw(n, s, e, signed=signed, ann=ann, norm=True)
    ann = annual() if ann is None else ann
    return L.X(ann, n, s, e, signed)

def x1_vector(cycles, s=0, e=1, signed=False, ann=None, normalised=False):
    ann = annual() if ann is None else ann
    return np.array([x1(n, s, e, signed, ann, normalised) for n in cycles], float)

# ------------------------------------------------------------------------------------ 4  polar field ----
def apf(variant="adopted"):
    """A_pf per cycle (15-24).  variant: adopted | mwo_csv | wso_csv | wso_daily | paper."""
    return L.apf(variant)

def apf_table():
    return L.apf_table()

def apf_sigma(correlated=False):
    """1-sigma of A_pf from the published SEMs (independent: SEM of the mean; correlated: conservative bound)."""
    return _p4().sigma_apf(correlated)

def polar_signed(reduction="adopted"):
    """Signed (B_N, B_S) at each minimum ending cycle n.  reduction: adopted | homogeneous | wso_yearly."""
    return C.polar()[reduction]

def apf_and_S(B):
    """(A_pf, S) from a signed pair (B_N, B_S); S = (B_N+B_S)/(|B_N|+|B_S|), the quadrupole form."""
    return C.apf_S(B)

# ------------------------------------------------------------------------------------------ 5  targets ----
def targets():
    """dict with Apf, S (cycles 15-24), An = whole-cycle area asymmetry A^n (12-25), Anext = {n: A^(n+1)} (12-23),
    B (signed polar), Dobs, dA, ann, d, sig_ind, sig_cor -- the cached observational layer of Phase 4."""
    D = _p4().data()
    out = dict(D)
    out["Anext"] = {n: D["An"][n + 1] for n in C2}
    return out

# --------------------------------------------------------------------------------------- 6  statistics ----
pearson, runs_test = L.pearson, L.runs_test
corr, fisher_ci, boot_ci = C.corr, C.fisher_ci, C.boot_ci

def spearman(x, y):
    from scipy.stats import spearmanr
    return float(spearmanr(x, y)[0])

def ols(x, y): return _p4().ols(x, y)
def loo_pred(x, y): return _p4().loo_pred(x, y)
def loo_clim(y): return _p4().loo_clim(y)
def multi_loo_pred(Xm, y): return _p4().multi_loo_pred(Xm, y)
def nested_select(cands, y): return _p4().nested_select(cands, y)
def skill(pred, y): return _p4().skill(pred, y)
def controls(x, y, cycles, xn=None): return _p4().controls(x, y, cycles, xn)

def rows_r(X, y):
    """Pearson r of every row of X (k, n) with y (n,)."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    Xz = (X - X.mean(1, keepdims=True)) / X.std(1, keepdims=True); yz = (y - y.mean()) / y.std()
    return Xz @ yz / len(y)

def perm_max_r(X, y, n_perm=10000, seed=20261019):
    """Family-wise p of the in-sample maximum r over the candidate rows of X: permute y over cycles, recompute the max.
    Returns (r_max, argmax row, p)."""
    rng = np.random.default_rng(seed); r = rows_r(X, y); rmax = float(np.nanmax(r)); cnt = 0
    for _ in range(n_perm):
        cnt += np.nanmax(rows_r(X, rng.permutation(y))) >= rmax - 1e-12
    return rmax, int(np.nanargmax(r)), (cnt + 1) / (n_perm + 1)

def loo_mse_rows(X, y):
    """Closed-form LOO MSE of y ~ a + b x for every row of X (k, n).  (verbatim from revision/p4/ws_g.py)"""
    X = np.asarray(X, float); y = np.asarray(y, float)
    n = len(y); xm = X.mean(1, keepdims=True); xc = X - xm; Sxx = (xc ** 2).sum(1); ok = Sxx > 1e-300
    b = np.where(ok, (xc * (y - y.mean())).sum(1) / np.where(ok, Sxx, 1), 0.0)
    r = y - y.mean() - b[:, None] * xc; h = 1.0 / n + xc ** 2 / np.where(ok, Sxx, 1)[:, None]
    out = np.mean((r / (1 - h)) ** 2, axis=1); out[~ok] = np.inf; return out

def nested_window_selection(X, y):
    """Nested leave-one-cycle-out over candidate rows of X (k, n): in each fold choose the row with the lowest inner
    LOO MSE on the training cycles, refit, predict the held-out cycle.  Returns (predictions, chosen row per fold).
    This is the criterion now used for every nested value in the manuscript.  (verbatim from ws_g.py)"""
    X = np.asarray(X, float); y = np.asarray(y, float); n = len(y); p = np.empty(n); ch = []
    for i in range(n):
        tr = np.arange(n) != i; j = int(np.argmin(loo_mse_rows(X[:, tr], y[tr]))); ch.append(j)
        b, a = np.polyfit(X[j, tr], y[tr], 1); p[i] = a + b * X[j, i]
    return p, ch

# ---------------------------------------------------------------------------------------------- 7  SFT ----
class SFTModel:
    """1-D surface flux transport (bn_sft: eta 250 km^2/s, u0 15 m/s, lambda0 75 deg, tau 5 yr) driven by every observed
    region, read out as the observations are taken (+-1.6 yr about each minimum), through the cap Green's function.

        m = SFTModel(); r = m.readout(seed=None)          # seed None = deterministic Joy tilts
        r['full'][21] -> (N, S) cap readout incl. the initial field; r['own'][21] -> cycle 21's own regions only
        A.apf_of(r['full'][21])

    Building the Green's function takes ~1 minute the first time."""
    def __init__(self, wname="cap70"):
        X = _x(); P = _p4()
        self.X, self.P, self.wname = X, P, wname
        self.solver = X.S
        self.cat = X.load_catalogue()
        self.G, self.nk = P.green(self.solver)
        self.IC = P.ic(self.solver)

    def tilts(self, seed=None):
        """Tilt of every region: catalogue Joy tilts (seed None) or one draw of the size-dependent scatter."""
        return self.X.tilts(self.cat, seed)

    def readout(self, seed=None, tilt=None, keep_regions=False):
        tilt = self.tilts(seed) if tilt is None else tilt
        return self.P.run(self.solver, self.G, self.nk, self.IC, self.cat, tilt, keep_regions=keep_regions, wname=self.wname)

def HemiSFT(etaN=250.0, etaS=250.0, v0N=15.0, v0S=15.0, tau_yr=5.0):
    """The bn_sft solver with hemisphere-dependent diffusivity and flow amplitude (Phase 4, WS-A)."""
    return _p4().HemiSFT(etaN, etaS, v0N, v0S, tau_yr)

def apf_of(NS): return _p4().apf_of(NS)
def d_of(NS): return _p4().d_of(NS)

def ensemble():
    """The pre-registered 101-member ensemble (out/x1_ensemble.npz): full/own readouts (config, member, cycle, N/S),
    per-region contributions for analogs A1-A3, region metadata."""
    return np.load(os.path.join(L.OUT, "x1_ensemble.npz"), allow_pickle=True)

def ensemble_apf(analog="A1", cycles=C1, member=None):
    """Modelled A_pf per cycle from the ensemble: member 0 = deterministic, None = all members (shape members x cycles)."""
    Z = ensemble(); F = Z["full"][list(Z["cfgs"]).index(analog)][:, [list(Z["cycles"]).index(n) for n in cycles]]
    a = (np.abs(F[..., 0]) - np.abs(F[..., 1])) / (np.abs(F[..., 0]) + np.abs(F[..., 1]))
    return a if member is None else a[member]

# ----------------------------------------------------------------------------------------- 8  plotting ----
COLORS = dict(N="#0072B2", S="#D55E00", giant="#E69F00", rest="#56B4E9", total="#000000", green="#009E73",
              purple="#CC79A7", yellow="#F0E442")           # Wong colour-blind-safe palette, as in every figure

def theme():
    """plotnine theme of the paper figures: Times 8 pt, no minor grid, legend on top."""
    return _p4().theme()

def save(fig, name, w_mm=85, h_mm=65, outdir="."):
    """Save a plotnine figure as vector PDF (TrueType fonts) + 600-dpi PNG; 85 mm = one column, 180 mm = full width."""
    import matplotlib
    matplotlib.rcParams["pdf.fonttype"] = 42; matplotlib.rcParams["ps.fonttype"] = 42
    os.makedirs(outdir, exist_ok=True)
    for ext, kw in (("pdf", {}), ("png", dict(dpi=600))):
        fig.save(os.path.join(outdir, f"{name}.{ext}"), width=w_mm / 25.4, height=h_mm / 25.4, units="in", verbose=False, **kw)

# ------------------------------------------------------------------------------------------- 9  checks ----
def headline():
    """The two channels of the paper, recomputed from scratch."""
    T = targets(); ann = T["ann"]
    x_ab = x1_vector(C1, 0, 1, signed=True, ann=ann); y_ab = np.array([T["Apf"][n] for n in C1])
    x_ac = x1_vector(C2, -1, 4, signed=False, ann=ann); y_ac = np.array([T["Anext"][n] for n in C2])
    return dict(channel1=dict(window="t_R+0..+1, parity-signed", n=len(C1), r=corr(x_ab, y_ab),
                              signs=int((np.sign(x_ab) == np.sign(y_ab)).sum())),
                channel2=dict(window="t_R-1..+4, unsigned", n=len(C2), r=corr(x_ac, y_ac),
                              signs=int((np.sign(x_ac) == np.sign(y_ac)).sum())))

def selftest(with_sft=False):
    import pandas as pd
    ok = True
    def check(name, cond, detail=""):
        nonlocal ok
        ok &= bool(cond); print(f"  {'PASS' if cond else 'FAIL'}  {name}  {detail}")
    check("catalogue checksum", L._md5(RGO) == L.MD5[RGO])
    h = headline()
    check("channel 1 r = 0.927, 10/10 signs", abs(h["channel1"]["r"] - 0.927) < 5e-4 and h["channel1"]["signs"] == 10,
          f"(r = {h['channel1']['r']:.4f}, signs {h['channel1']['signs']}/10)")
    check("channel 2 r = 0.771, 10/12 signs", abs(h["channel2"]["r"] - 0.771) < 5e-4 and h["channel2"]["signs"] == 10,
          f"(r = {h['channel2']['r']:.4f}, signs {h['channel2']['signs']}/12)")
    d = daily(); t1 = pd.read_csv(os.path.join(L.OUT, "table1_region.csv")).set_index("cycle")
    v = cycle_asymmetry(12, size_mask(d, 500), d)
    check("Table 1, SC12 large (per region)", abs(v - t1.loc[12, "A_large"]) < 5e-5, f"({v:+.4f})")
    e = ensemble_apf("A1", [21], member=0)[0]
    check("ensemble: modelled A_pf(21), deterministic A1", abs(e - (-0.059)) < 5e-4, f"({e:+.4f})")
    if with_sft:
        m = SFTModel(); r = m.readout(None)
        a21 = apf_of(r["full"][21])
        check("SFT readout reproduces the ensemble member 0 for SC21", abs(a21 - e) < 1e-6, f"({a21:+.6f})")
    print("all checks passed" if ok else "SOME CHECKS FAILED")
    return ok

if __name__ == "__main__":
    sys.exit(0 if selftest("--sft" in sys.argv) else 1)
