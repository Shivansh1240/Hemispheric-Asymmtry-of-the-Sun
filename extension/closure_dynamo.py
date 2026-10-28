"""closure_dynamo.py -- the dynamo runs for closure tests T3 and T4 (new file; q4_dynamo.py is imported, not changed).

T3  memory the dynamo needs: forcing coherence time tau_c = 1, 3, 5, 10, 30 yr, kappa = 1, rho = 0, 3 seeds each.
    sigma is re-tuned per tau_c so that the rms cycle asymmetry rms A is 0.10 (the Sun's, and that of the kappa = 1
    production setting at tau_c = 1).  One calibration step: two 600-yr runs at a first-guess sigma, then
    sigma_prod = sigma_cal * 0.10 / rms A_cal.  tau_c = 1 reuses the existing production seeds 1-3 (sigma 0.6).
T4  constant hemispheric bias: a deterministic difference delta = 1, 3, 10 % between the hemispheres in
      alpha0  (BL source amplitude: N x (1 + delta/2), S x (1 - delta/2)),
      v0      (meridional flow: the stream function times 1 + (delta/2) tanh(lat / 3 deg); the flow stays
               divergence-free because both face fluxes still come from one corner array of the stream function),
      eta     (turbulent diffusivity, poloidal and toroidal, radial and latitudinal: N x (1 + delta/2), S x (1 - delta/2);
               the radial Crank-Nicolson operator is factorised per hemisphere; latitudinal diffusion in flux form),
    with the production noise on (sigma 0.6, tau_c 1 yr, rho 0, kappa 1), 3 seeds each.  The unbiased reference is
    the existing production setting kappa = 1, seeds 1-3.

    python closure_dynamo.py check     delta = 0 must reproduce SXD exactly; delta > 0 must run
    python closure_dynamo.py launch    all T3 / T4 runs (process pool), then the analysis
    python closure_dynamo.py analyse   per-cycle statistics of every finished run -> out/closure/dynamo_T3_T4.json
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import q4_dynamo as QD                                          # noqa: E402
from q4_dynamo import SXD                                        # noqa: E402

OUTC = os.path.join(HERE, "out", "closure"); OUTD = os.path.join(OUTC, "dynamo"); os.makedirs(OUTD, exist_ok=True)
BASE = dict(alpha0=40.0, wlat=35.0, kappa=1.0, nr=65, nth=121, v0=1500.0)        # the production model
NOISE = dict(sigma=0.6, rho=0.0, tau_c=1.0)
TAUS = [1.0, 3.0, 5.0, 10.0, 30.0]
SIG_GUESS = {3.0: 0.25, 5.0: 0.18, 10.0: 0.12, 30.0: 0.10}       # first guesses; the calibration step sets the value used
SEEDS = [1, 2, 3]
REF_TAU1 = ["prod_k1", "prod_k1s2", "prod_k1_s3"]                 # existing kappa = 1 production runs, seeds 1-3
DELTAS = [0.01, 0.03, 0.10]
KINDS = ["alpha", "v0", "eta"]
TARGET_RMS_A = 0.10


class SXB(SXD):
    """SXD with a constant hemispheric difference in alpha0, v0 or eta (kind), N larger by delta."""

    def __init__(self, kind="none", delta=0.0, w_h=3.0, **kw):
        self.kind, self.delta, self.w_h = kind, float(delta), w_h
        super().__init__(**kw)

    def _hsmooth(self, lat):
        return 1.0 + 0.5 * self.delta * np.tanh(np.asarray(lat) / self.w_h)

    # ---- v0: hemispheric factor on the stream function, normalisation restored so N = v0 (1 + delta/2)
    def _stream(self, rr, tt):
        out = super()._stream(rr, tt)
        return out * self._hsmooth(90.0 - np.degrees(tt)) if self.kind == "v0" else out

    def _flow(self):
        super()._flow()
        if self.kind != "v0":
            return
        g = float(self._hsmooth(45.0))
        for a in ("Fr_f", "Fth_f", "vr_f", "vth_f"):
            setattr(self, a, getattr(self, a) * g)
        self.psi_scale *= g
        self.vr_c = 0.5 * (self.vr_f[1:] + self.vr_f[:-1]); self.vth_c = 0.5 * (self.vth_f[:, 1:] + self.vth_f[:, :-1])
        self.Fr = self.Fr_f[1:-1, 1:-1]; self.Fth = self.Fth_f[1:-1, 1:-1]
        self.vr_if = self.vr_f[1:-1, 1:-1]; self.vth_if = self.vth_f[1:-1, 1:-1]

    # ---- eta: column factor hc (N 1 + delta/2, S 1 - delta/2, equator 1), face factor hf
    def _profiles(self):
        super()._profiles()
        h = 0.5 * self.delta if self.kind == "eta" else 0.0
        self.hc = np.where(self.lat > 1e-9, 1.0 + h, np.where(self.lat < -1e-9, 1.0 - h, 1.0))
        self.hf = (0.5 * (self.hc[1:] + self.hc[:-1]))[None, :]

    def _timestep(self):
        super()._timestep()
        if self.kind == "eta":
            self.dt /= 1.0 + 0.5 * abs(self.delta)

    def _factorise(self):
        super()._factorise()
        if self.kind == "eta":
            groups = [np.flatnonzero(np.isclose(self.hc, h)) for h in np.unique(self.hc)]
            self.TAg = [(g, self._tri(self.eta_p * self.hc[g[0]], True)) for g in groups]
            self.TBg = [(g, self._tri(self.eta_t * self.hc[g[0]], False)) for g in groups]

    @staticmethod
    def _solve(T, d):
        out = np.empty_like(d)
        for g, M in T:
            out[:, g] = M.solve(d[:, g])
        return out

    def _lap_theta_h(self, F_, eta_r):
        sF = self.sin * F_
        G = eta_r[:, None] * self.facf * self.hf * (sF[:, 1:] - sF[:, :-1]) / (self.dth * np.sin(self.thh)[None, :])
        out = np.zeros_like(F_)
        out[:, 1:-1] = (G[:, 1:] - G[:, :-1]) / (self.dth * self.R2[:, 1:-1] ** 2)
        return out

    # ---- alpha: hemispheric factor on the BL source amplitude
    def source(self, B, aN, aS):
        if self.kind == "alpha":
            aN, aS = aN * (1.0 + 0.5 * self.delta), aS * (1.0 - 0.5 * self.delta)
        return super().source(B, aN, aS)

    def step(self, A, B, S=None):
        if self.kind != "eta":
            return super().step(A, B, S)
        dt, H = self.dt, self.hc[None, :]
        slope = self.M_top @ A[-1]
        U = self.s * A
        rhs = self._upwind_A(U, A, slope) + self._lap_theta_h(A, self.eta_p)
        if S is not None:
            rhs = rhs + S
        Dr = np.zeros_like(A)
        Dr[1:-1] = self.ETA_P[1:-1] * H * self._lap_r(A)[1:-1]
        Dr[-1] = self.eta_p[-1] * self.hc * self._lap_r_top(A, slope)
        star = A + dt * rhs + 0.5 * dt * Dr
        star[-1] += self.hc * (0.5 * dt * self.eta_p[-1] * (2 * self.dr * slope) / (self.dr ** 2) +
                               0.5 * dt * self.eta_p[-1] * (2 * self.dr * slope) / (2 * self.dr) * 2.0 / self.r[-1])
        An = self._solve(self.TAg, star)
        Br, Bth = self.bpol(A)
        shear = self.s * (Br * self.dom_dr + Bth * self.dom_dth / self.R2)
        rB = self.R2 * B
        dgrad = np.zeros_like(B)
        dgrad[1:-1] = (self.deta_t[1:-1, None] * H * (rB[2:] - rB[:-2]) / (2 * self.dr)) / self.R2[1:-1]
        rhsB = self._adv_B(B) + self._lap_theta_h(B, self.eta_t) + shear + dgrad
        DrB = np.zeros_like(B)
        DrB[1:-1] = self.ETA_T[1:-1] * H * self._lap_r(B)[1:-1]
        Bn = self._solve(self.TBg, B + dt * rhsB + 0.5 * dt * DrB)
        An[0] = 0.0; An[:, 0] = 0.0; An[:, -1] = 0.0
        Bn[0] = 0.0; Bn[-1] = 0.0; Bn[:, 0] = 0.0; Bn[:, -1] = 0.0
        return An, Bn


def run_one(name, years, seed, sigma, tau_c, kind="none", delta=0.0):
    path = os.path.join(OUTD, f"{name}.npz")
    if os.path.exists(path):
        return name, "exists"
    m = SXB(kind=kind, delta=delta, **BASE)
    res = m.run(years, seed=seed, sigma=sigma, rho=0.0, tau_c=tau_c)
    meta = dict(name=name, years=years, seed=seed, sigma=sigma, tau_c=tau_c, kind=kind, delta=delta, **BASE)
    np.savez_compressed(path, **{k: v for k, v in res.items() if k not in ("A", "B")}, meta=json.dumps(meta))
    return name, f"{res['wall']:.0f} s"


def load(name):
    for d in (OUTD, QD.OUT):
        p = os.path.join(d, f"{name}.npz")
        if os.path.exists(p):
            return dict(np.load(p, allow_pickle=True))
    raise FileNotFoundError(name)


def rms_A(names, spin=100.0):
    from q4_dynamo_analyse import per_cycle
    A = np.concatenate([[c["A"] for c in per_cycle(load(n), spin=spin)] for n in names])
    return float(np.sqrt(np.mean(A ** 2))), len(A)


def check():
    out = {}
    for kind in KINDS:
        a = SXD(**BASE).run(15.0, seed=7, sigma=0.6, tau_c=1.0)
        b = SXB(kind=kind, delta=0.0, **BASE).run(15.0, seed=7, sigma=0.6, tau_c=1.0)
        d0 = float(max(np.max(np.abs(a[k] - b[k])) / np.max(np.abs(a[k])) for k in ("PN", "PS", "EN", "ES")))
        c = SXB(kind=kind, delta=0.10, **BASE).run(15.0, seed=7, sigma=0.6, tau_c=1.0)
        out[kind] = dict(delta0_max_rel_diff=d0, delta10_finite=bool(np.isfinite(c["A"]).all()),
                         delta10_rel_change_PN=float(np.max(np.abs(np.interp(a["t"], c["t"], c["PN"]) - a["PN"])) / np.max(np.abs(a["PN"]))))
        print(kind, out[kind], flush=True)
    m = SXB(kind="v0", delta=0.10, **BASE); j = int(np.argmin(np.abs(m.lat - 45.0))); k = int(np.argmin(np.abs(m.lat + 45.0)))
    out["v0_surface_speed_45N_45S_m_s"] = [float(m.vth_f[-1, j] / 100), float(m.vth_f[-1, k] / 100)]
    print("v0 delta 0.10: surface v_theta at 45N / 45S (m/s)", out["v0_surface_speed_45N_45S_m_s"])
    json.dump(out, open(os.path.join(OUTC, "dynamo_checks.json"), "w"), indent=1)


def launch(workers=11):
    from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
    log = open(os.path.join(OUTC, "launch.log"), "a")
    def say(s):
        print(s, flush=True); log.write(time.strftime("%H:%M:%S ") + s + "\n"); log.flush()
    cal = {tc: [f"T3cal_tc{tc:g}_s{s}" for s in (101, 102)] for tc in SIG_GUESS}
    jobs = {}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for tc, names in cal.items():
            for nm, s in zip(names, (101, 102)):
                jobs[ex.submit(run_one, nm, 600.0, s, SIG_GUESS[tc], tc)] = nm
        for kind in KINDS:
            for dl in DELTAS:
                for s in SEEDS:
                    nm = f"T4_{kind}_d{int(round(dl * 100))}_s{s}"
                    jobs[ex.submit(run_one, nm, 1100.0, s, NOISE["sigma"], NOISE["tau_c"], kind, dl)] = nm
        say(f"submitted {len(jobs)} runs")
        cal_names = {n for v in cal.values() for n in v}; prod_sent = False
        pending = set(jobs)
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for f in done:
                say(f"finished {jobs[f]}: {f.result()[1]}")
            if not prod_sent and all(not any(jobs[f] == n for f in pending) for n in cal_names):
                sig = {}
                for tc, names in cal.items():
                    r, n = rms_A(names)
                    sig[tc] = dict(sigma_cal=SIG_GUESS[tc], rms_A_cal=r, n_cycles_cal=n, sigma_prod=SIG_GUESS[tc] * TARGET_RMS_A / r)
                    say(f"calibration tau_c {tc:g}: sigma {SIG_GUESS[tc]} -> rms A {r:.3f} (n = {n}) -> sigma {sig[tc]['sigma_prod']:.3f}")
                json.dump({f"{k:g}": v for k, v in sig.items()}, open(os.path.join(OUTC, "T3_calibration.json"), "w"), indent=1)
                for tc in SIG_GUESS:
                    for s in SEEDS:
                        nm = f"T3_tc{tc:g}_s{s}"
                        f = ex.submit(run_one, nm, 1100.0, s, sig[tc]["sigma_prod"], tc); jobs[f] = nm; pending.add(f)
                prod_sent = True
                say("T3 production submitted")
    say("all runs finished")


if __name__ == "__main__":
    {"check": check, "launch": launch}[sys.argv[1]]()
