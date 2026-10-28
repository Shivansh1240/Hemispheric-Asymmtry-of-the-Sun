"""sft_bhowmik.py -- the SFT analyses rerun with Bhowmik's SFT settings (v8).  New file; nothing existing is changed.

Bhowmik (2019, Sects. 3.1 and 6; the data-driven SFT of Bhowmik & Nandy 2018) differs from our SFT in two places:
  * no decay term: our tau = 5 yr is removed (tau = infinity);
  * tilt alpha = g T_n sqrt|lat| with g = 0.7 (inflows, Cameron et al. 2010; her Sect. 3.1.3) and the cycle-amplitude
    law of Jiang et al. (2011) Eq. 15, T_n = 1.73 - 0.0035 S_n, S_n the cycle maximum of the 12-month running mean group
    sunspot number (their Table 1, cycles 11-23).  S_24 and S_25 are not in that table: they are the SILSO V2 maxima of
    the same running mean scaled by S_23 / V2_23 = 119 / 182.2 (76.7 and 105.2).  She adds no tilt scatter in the
    data-driven run; our scatter ensembles keep the Jiang et al. (2014) scatter she uses in her Sect. 4.3, now around
    her mean tilt.
Transport (eta 250 km^2/s, v0 15 m/s, lambda0 75 deg) and flux (7e19 A Mx) are already hers.  Her 2019 hemispheric
law T = 1.71 - 0.12 F (F the peak monthly hemispheric flux) is calibrated on synthetic cycles about 20 times weaker
than the RGO record and gives T < 0 for every observed cycle; it is not used.

How: each analysis script runs unchanged except that, before it runs,
  * core.CAT.alpha is replaced by her tilt (every deterministic tilt, scatter draw around it, mirrored / symmetrised /
    re-dealt catalogue and T_n = alpha / sqrt|lat| inherits it; the catalogue tilt is kept as CAT.alpha_cat),
  * core.Model defaults to tau = infinity with a 150-yr Green's function (without decay old flux does not fade, so the
    45-yr kernel is too short; 1878-2021 fits in 150 yr); kernels are cached in out/v8_sft/G_cache,
  * outputs go to out/v8_sft/ (figures of the scripts to out/v8_sft/figs) instead of out/ and figs/,
  * a few literal decay settings are read as the base tau (listed in SUBS, applied to the executed copy only and
    recorded in out/v8_sft/substitutions.json),
  * Python runs with -O, so the scripts' asserts that pin the tau = 5 yr numbers are off.

    python -O sft_bhowmik.py run all            every script in ORDER
    python -O sft_bhowmik.py run q3_persistence closure_v5:n1 ...
"""
import os, sys, json, shutil, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import core                                                   # noqa: E402

OUT8 = os.path.join(HERE, "out", "v8_sft"); GC = os.path.join(OUT8, "G_cache"); FIG8 = os.path.join(OUT8, "figs")
for d in (OUT8, GC, FIG8, os.path.join(OUT8, "v5")):
    os.makedirs(d, exist_ok=True)
G_YEARS = 150.0
G_INFLOW = 0.7
SN = {11: 101, 12: 67, 13: 96, 14: 64, 15: 109, 16: 81, 17: 123, 18: 143, 19: 186, 20: 108, 21: 154, 22: 156, 23: 119,
      24: 76.7, 25: 105.2}                                   # Jiang et al. 2011 Table 1; 24, 25 from SILSO V2 (see above)
BASE_TAU = np.inf

def T_n(n):
    return 1.73 - 0.0035 * SN[int(n)]

def her_tilt(cat):
    return G_INFLOW * np.array([T_n(c) for c in cat.cycle.values]) * np.sqrt(np.abs(cat.latc.values)) * np.sign(cat.alpha.values)

# ------------------------------------------------------------------------------------------------- the patch ----
_Model0 = core.Model

class ModelB(_Model0):
    """core.Model with Bhowmik's default tau (none) and a 150-yr cached Green's function."""
    def __init__(self, eta=250.0, v0=15.0, tau=BASE_TAU, etaN=None, etaS=None, v0N=None, v0S=None, G=None, base=None):
        if G is None and base is None:
            key = f"eta{eta:g}_v0{v0:g}_tau{tau:g}_eN{etaN}_eS{etaS}_vN{v0N}_vS{v0S}"
            f = os.path.join(GC, key + ".npy")
            if os.path.exists(f):
                G = np.load(f)
            else:
                t0 = time.time()
                G = core.green(core.solver(eta, v0, tau, etaN, etaS, v0N, v0S), years=G_YEARS); np.save(f, G)
                print(f"[sft_bhowmik] Green's function {key}: {time.time() - t0:.0f} s", flush=True)
        super().__init__(eta, v0, tau, etaN, etaS, v0N, v0S, G=G, base=base)

def patch():
    if "alpha_cat" not in core.CAT:
        core.CAT["alpha_cat"] = core.CAT["alpha"].values.copy()
        core.CAT["alpha"] = her_tilt(core.CAT)
    core.Model = ModelB
    core.OUT = OUT8; core.FIG = FIG8

# ------------------------------------------------------------------------------- scripts and substitutions ----
SUBS = {
    "q1_sources": [("R_SUN, TAU, U0, LAM0, SIG_SPOT = 6.96e10, 5.0,", "R_SUN, TAU, U0, LAM0, SIG_SPOT = 6.96e10, _BASE_TAU,")],
    "q3_persistence": [("tau_list = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 13.0, 17.0]", "tau_list = [2.0, 3.0, 4.0, 5.0, 7.0, 10.0, 13.0, 17.0, _BASE_TAU]"),
                       ("tau == 5.0", "tau == _BASE_TAU")],
    "ens_verify": [("for tau in (3.0, 5.0, 10.0):", "for tau in (_BASE_TAU,):"), ("tau=5.0)", "tau=_BASE_TAU)"), ("tau == 5.0", "tau == _BASE_TAU")],
    "ens_transport_series": [("for tau in (3.0, 5.0, 10.0):", "for tau in (_BASE_TAU,):"), ("tau=5.0)", "tau=_BASE_TAU)"), ("tau == 5.0", "tau == _BASE_TAU")],
    "closure_v5": [("tau == 5.0", "tau == _BASE_TAU")],
    # the per-cycle model series of the figures: the base model (no decay) instead of tau = 5 yr; the key name is kept
    "plots_src/export_series": [('q3["tau_sweep"]["5.0"]', 'q3["tau_sweep"]["Infinity"]')],
}
ORDER = ["q3_persistence", "q1_sources", "q2_window", "q2b_phase_sign", "q4_stochastic", "ens_verify", "ens_extra",
         "ens_transport_series", "explore", "figs_main", "closure_v5:n1", "closure_v5:n1x", "closure_v5:n2", "closure_v5:n3",
         "closure_v5:n3b", "v0_checks", "plots_src/export_series"]

def inherit():
    """copy the top-level data products once, so data-only inputs exist; every SFT product is then overwritten."""
    mark = os.path.join(OUT8, "inherited.json")
    if os.path.exists(mark):
        return
    rec = {}
    for sub in ("", "v5"):
        src = os.path.join(HERE, "out", sub)
        for f in os.listdir(src):
            p = os.path.join(src, f)
            if os.path.isfile(p) and f.split(".")[-1] in ("json", "npz", "pkl"):
                shutil.copy2(p, os.path.join(OUT8, sub, f)); rec[os.path.join(sub, f)] = os.path.getmtime(p)
    json.dump(rec, open(mark, "w"), indent=1)

def run_one(item):
    name, _, func = item.partition(":")
    src = open(os.path.join(HERE, name + ".py"), encoding="utf-8").read()
    done = []
    for old, new in SUBS.get(name, []):
        k = src.count(old)
        if k == 0:
            raise RuntimeError(f"{name}: substitution target not found: {old!r}")
        src = src.replace(old, new); done.append(dict(old=old, new=new, count=k))
    log = os.path.join(OUT8, "substitutions.json")
    S = json.load(open(log)) if os.path.exists(log) else {}
    S[name] = done; json.dump(S, open(log, "w"), indent=1)
    g = dict(__name__="__main__", __file__=os.path.join(HERE, name + ".py"), _BASE_TAU=BASE_TAU)
    argv0 = sys.argv
    sys.argv = [name + ".py"] + ([func] if func else [])
    t0 = time.time()
    try:
        exec(compile(src, os.path.join(HERE, name + ".py"), "exec"), g)
    finally:
        sys.argv = argv0
    print(f"[sft_bhowmik] {item} done in {time.time() - t0:.0f} s", flush=True)

def main():
    if sys.flags.optimize < 1:
        raise SystemExit("run with python -O (the tau = 5 yr asserts must be off)")
    items = ORDER if sys.argv[2:] == ["all"] else sys.argv[2:]
    patch(); inherit()
    for it in items:
        run_one(it)

if __name__ == "__main__":
    {"run": main}[sys.argv[1]]()
