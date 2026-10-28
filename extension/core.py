"""core.py -- shared layer for analysis_v3 (Q1-Q3, ensembles, exploratory track A-E).

Everything is built on the verified project modules through master_code/asymcore.py (same catalogue
RGO_SOON_corrected.txt, same A_pf reduction, same bn_sft engine and cap-70 Green's-function readout), so a number
here is comparable with the manuscript.  New code is limited to:
  * a cap-70-only Green's function that can be rescaled for the decay time (G_tau = G_5 exp(-k dt (1/tau-1/5)),
    exact because the decay is a scalar factor per step; declining-phase pipeline, dpp/green.py),
  * a readout that keeps the contribution of EVERY region to every minimum (not only the cycle's own regions),
  * catalogue operations (mirror, symmetrise, rescale, move a region).

Statistics are deliberately plain: Pearson r with n, sign counts, run lengths.
"""
import os, sys, json, time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out"); FIG = os.path.join(HERE, "figs")
os.makedirs(OUT, exist_ok=True); os.makedirs(FIG, exist_ok=True)
for _d in ("master_code", "pipeline"):                      # pipeline/ is master_code/ in the repository
    sys.path.insert(0, os.path.join(HERE, "..", _d))
import asymcore as A                                          # noqa: E402

P = A._p4(); X = A._x(); BN = P.BN; CM = A.C
from sft_asym import T0, scatter_sigma                        # noqa: E402

SCB, TR = A.SCB, A.TR
CYC = list(X.CYCLES)                 # 12..24: cycles with a model readout at the minimum that ends them
C1, C2 = A.C1, A.C2                  # 15..24 (measured A_pf), 12..23 (next-cycle pairs)
parity = A.parity
SEC = BN.SEC_PER_YEAR
WCAP = np.column_stack(X.WEIGHTS["cap70"])                    # (nlat, 2): area-weighted mean over |lat| >= 70
LAT = X.LAT

def r(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    return float(np.corrcoef(x[m], y[m])[0, 1])

def signs(x, y):
    return int(np.sum(np.sign(np.asarray(x)) == np.sign(np.asarray(y))))

def dump(name, obj):
    def conv(v):
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, (np.floating,)): return float(v)
        if hasattr(v, "tolist"): return v.tolist()
        return str(v)
    json.dump(obj, open(os.path.join(OUT, name), "w"), indent=1, default=conv)

def load(name):
    return json.load(open(os.path.join(OUT, name)))

# ------------------------------------------------------------------------------------------ observations ----
T = A.targets()
D = T["d"]                           # daily rows: yr dec hemi area abslat lat gid cyc rmax rlat
ANN = T["ann"]
APF = {n: float(T["Apf"][n]) for n in C1}
S_SIGNED = {n: float(T["S"][n]) for n in C1}
AN = {n: float(T["An"][n]) for n in range(12, 26)}
SIG = {n: float(T["sig_ind"][n]) for n in C1}

def x1(n, s=0, e=1, signed=False, ann=None):
    return float(A.x1(n, s, e, signed, ANN if ann is None else ann))

def x1v(cyc, s=0, e=1, signed=False, ann=None):
    return np.array([x1(n, s, e, signed, ann) for n in cyc])

def band_mask(d, lo=0.0, hi=np.inf, latlo=0.0, lathi=90.0):
    """regions with peak area in [lo, hi) and |lat| at peak in [latlo, lathi)."""
    return (d.rmax >= lo) & (d.rmax < hi) & (d.rlat >= latlo) & (d.rlat < lathi)

def ann_of(mask=None, d=None):
    d = D if d is None else d
    return A.annual(d, mask)

def cycA(n, mask=None, d=None):
    return float(A.cycle_asymmetry(n, mask, D if d is None else d))

# ---------------------------------------------------------------------------------------------- catalogue ----
CAT = X.load_catalogue()             # one row per region at maximum flux (year, cycle, latc, area, ddeg, flux, parity, ...)

def mirror(cat):
    c = cat.copy(); c["latc"] = -c.latc.values; return c

def scaled(cat, f):
    """flux x f and polarity separation x sqrt(f) (d ~ sqrt(A) in the catalogue)."""
    c = cat.copy(); f = np.asarray(f, float)
    c["flux"] = c.flux.values * f; c["ddeg"] = c.ddeg.values * np.sqrt(f); c["area"] = c.area.values * f
    return c

def symmetrised(cat):
    """every region replaced by two half-flux copies at +|lat| and -|lat|: identical hemispheres."""
    a = cat.copy(); b = cat.copy()
    a["latc"] = np.abs(cat.latc.values); b["latc"] = -np.abs(cat.latc.values)
    out = pd.concat([a, b], ignore_index=True)
    out["flux"] = out.flux.values * 0.5; out["area"] = out.area.values * 0.5
    out["src"] = np.r_[np.arange(len(cat)), np.arange(len(cat))]
    return out

# ------------------------------------------------------------------------------------------------- SFT ----
def solver(eta=250.0, v0=15.0, tau=5.0, etaN=None, etaS=None, v0N=None, v0S=None):
    if etaN is None and etaS is None and v0N is None and v0S is None:
        return BN.BNSFT(eta=eta * 1e10, v0=v0 * 100.0, tau_yr=tau)
    return P.HemiSFT(etaN if etaN is not None else eta, etaS if etaS is not None else eta,
                     v0N if v0N is not None else v0, v0S if v0S is not None else v0, tau)

def green(s, years=45.0):
    """G[j, (N,S), k]: cap-70 N/S response k steps after a unit source in latitude cell j."""
    nk = int(np.ceil(years * SEC / s.dt)) + 1
    B = np.eye(len(s.lat)); G = np.empty((len(s.lat), 2, nk))
    for k in range(nk):
        G[:, :, k] = B @ WCAP; B = s.step(B)
    return G

def retau(G, s, tau_new, tau_old):
    k = np.arange(G.shape[2]); dt = s.dt / SEC
    return G * np.exp(-k * dt * (1.0 / tau_new - 1.0 / tau_old))[None, None, :]

def cumG(G):
    Cc = np.zeros(G.shape[:2] + (G.shape[2] + 1,)); np.cumsum(G, axis=2, out=Cc[:, :, 1:]); return Cc

def ic_readout(s, h=1.6, tau_new=None, tau_old=None):
    """N/S cap readout of the evolved initial field (-sin lat at T0(12)), averaged over +-h about each minimum."""
    t0 = T0[12]; B = -np.sin(s.latr); ns = int(np.ceil((max(X.READ.values()) + h + 0.1 - t0) * SEC / s.dt)) + 1
    vals = np.empty((ns, 2))
    for k in range(ns):
        vals[k] = B @ WCAP; B = s.step(B)
    if tau_new is not None:
        vals = vals * np.exp(-np.arange(ns) * s.dt / SEC * (1.0 / tau_new - 1.0 / tau_old))[:, None]
    tk = t0 + np.arange(ns) * s.dt / SEC
    return {n: vals[np.abs(tk - X.READ[n]) <= h].mean(0) for n in CYC}

class Model:
    """One transport configuration: Green's function + initial-field readout.  m.run(cat, tilt) gives, for every
    cycle n = 12..24, the N/S cap readout at the minimum ending n (+-1.6 yr), its own-cycle part and (keep=True) the
    contribution of every region that reaches that readout."""
    def __init__(self, eta=250.0, v0=15.0, tau=5.0, etaN=None, etaS=None, v0N=None, v0S=None, G=None, base=None):
        self.s = solver(eta, v0, tau, etaN, etaS, v0N, v0S)
        self.par = dict(eta=eta, v0=v0, tau=tau, etaN=etaN, etaS=etaS, v0N=v0N, v0S=v0S)
        if base is not None:                      # same transport, different tau: exact rescaling of base's kernel
            self.G = retau(base.G, self.s, tau, base.par["tau"])
            self.IC = ic_readout(base.s, tau_new=tau, tau_old=base.par["tau"])
        else:
            self.G = green(self.s) if G is None else G
            self.IC = ic_readout(self.s)
        self.C = cumG(self.G); self.nk = self.G.shape[2]

    def run(self, cat, tilt=None, keep=False, ic=True):
        s = self.s; tilt = cat.alpha.values if tilt is None else tilt
        Ss = X.sources(cat, tilt).astype(np.float64)
        ks = np.floor((cat.year.values - T0[12]) * SEC / s.dt).astype(int); cyc = cat.cycle.values
        full, own, reg = {}, {}, {}
        for n in CYC:
            sel, lo, hi, nK = P.ridx(ks, s, self.nk, 1.6, X.READ[n])
            g = (self.C[:, :, hi + 1] - self.C[:, :, lo]) / nK            # (nlat, 2, nsel)
            c = np.einsum("ij,jki->ik", Ss[sel], g)                     # (nsel, 2)
            full[n] = c.sum(0) + (np.asarray(self.IC[n]) if ic else 0.0)
            o = cyc[sel] == n; own[n] = c[o].sum(0)
            if keep:
                reg[n] = (sel, c)
        return dict(full=full, own=own, reg=reg)

def apf_of(NS):
    return P.apf_of(NS)

def d_of(NS):
    return P.d_of(NS)

def apf_vec(res, cyc=C1, key="full"):
    return np.array([apf_of(res[key][n]) for n in cyc])

def d_vec(res, cyc=C1, key="full"):
    return np.array([d_of(res[key][n]) for n in cyc])

def tilt_draw(cat, seed):
    """catalogue Joy tilt + one draw of the size-dependent scatter sigma(A) (the project's 101-member ensemble)."""
    return X.tilts(cat, seed)

# ------------------------------------------------------------------------------------------------- plots ----
WONG = dict(N="#0072B2", S="#D55E00", giant="#E69F00", rest="#56B4E9", total="#000000", green="#009E73",
            purple="#CC79A7", grey="#8c8c8c")

def style():
    import matplotlib as mpl
    mpl.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "DejaVu Serif"], "font.size": 8,
                         "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
                         "ytick.labelsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.6,
                         "lines.linewidth": 1.2, "legend.frameon": False, "mathtext.fontset": "stix",
                         "axes.spines.top": False, "axes.spines.right": False})

import re as _re
_PRETTY = [(r"X1hat", r"$\\hat{X}_1$"), (r"\bX_G\b", r"$X_G$"), (r"\bX1\b", r"$X_1$"), (r"A_pf", r"$A_{\\rm pf}$"),
           (r"t_R", r"$t_R$"), (r"\bB_N\b", r"$B_N$"), (r"\bB_S\b", r"$B_S$"), (r"\bA_N\b", r"$A_N$"), (r"\bA_S\b", r"$A_S$"),
           (r"D_window", r"$D_w$"), (r"T_lin", r"$T_{\\rm lin}$"), (r"\beta_N\b", r"$\\eta_N$"), (r"\beta_S\b", r"$\\eta_S$"),
           (r"->", r"$\\rightarrow$"), (r">=", r"$\\geq$"), (r"\bsigma\b", r"$\\sigma$"), (r"\btau\b", r"$\\tau$")]

def pretty(s):
    """plain figure text -> mathtext for the symbols used in this project (A_pf, X1hat, t_R, >= ...)."""
    if not s or "$" in s:
        return s
    for a, b in _PRETTY:
        s = _re.sub(a, b, s)
    return s

def savefig(fig, name, dpi=300):
    import matplotlib.text as mtext
    from matplotlib.ticker import FixedFormatter, FuncFormatter
    for ax in fig.axes:                       # string tick labels are re-created at draw time: convert at the source
        for axis in (ax.xaxis, ax.yaxis):
            fmt = axis.get_major_formatter()
            if isinstance(fmt, (FixedFormatter, FuncFormatter)):     # mpl 3.9 stores set_*ticklabels lists as FuncFormatter
                axis.set_major_formatter(FuncFormatter(lambda x, pos, f=fmt: pretty(f(x, pos))))
    for t in fig.findobj(mtext.Text):
        t.set_text(pretty(t.get_text()))
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, name + ".png"), dpi=dpi, bbox_inches="tight")

def tic(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)
