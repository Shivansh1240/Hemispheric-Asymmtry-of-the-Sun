"""q4_dynamo.py -- Q4 (part 2): the declining-phase 2-D flux-transport dynamo, closed into a self-excited
Babcock-Leighton loop on the full sphere, with stochastic forcing per hemisphere and a hemispheric-coupling knob.

The transport operator is the verified one (sc24_final/ftd_core.py, imported read-only; its eight checks pass:
q4_ftd_verify.py).  Added here, and nothing else:

  * BL source (non-local, Dikpati & Charbonneau 1999 / Chatterjee et al. 2004 form):
        S(r, th, t) = alpha_h(t) * cos(th) sin(th) * W(lat) * F(r) * Bbar(th) / (1 + (Bbar/B0)^2)
    Bbar = toroidal field averaged over the tachocline layer 0.70-0.72 R at the same colatitude, F(r) the
    near-surface deposition envelope of ftd_core.deposit (r_dep 0.95 R), W confines emergence to |lat| < 45 deg.
  * stochastic forcing per hemisphere: alpha_h = alpha0 (1 + sigma xi_h), xi_h ~ N(0,1) renewed every tau_c,
    xi_S = rho xi_N + sqrt(1-rho^2) xi'  (rho = 0: independent hemispheres).
  * hemispheric coupling kappa: the LATITUDINAL diffusivity (both A and B) is multiplied by
        1 + (kappa - 1) exp(-(lat/lat_c)^2),  lat_c = 10 deg,
    i.e. the cross-equatorial diffusive path is made weaker (kappa < 1) or stronger (kappa > 1).  The radial
    operator is untouched, so the verified implicit radial solver still applies.  kappa = 1 is the published model.

Readouts every `rec` years: the BL flux deposited per hemisphere (the model's "emergence"), the area-weighted surface
B_r over |lat| >= 70 deg per hemisphere (the polar field), and the tachocline Bbar(lat) for butterflies.

    python q4_dynamo.py tune                         short deterministic runs for the period and growth
    python q4_dynamo.py run NAME years seed kappa sigma rho [alpha0 v0 nr nth]
"""
import os, sys, json, time
import numpy as np
from scipy.special import erf

HERE = os.path.dirname(os.path.abspath(__file__))
FTD_DIR = os.environ.get("FTD_DIR", r"C:\Users\Shivansh Mishra\OneDrive\Desktop\declining phase paper\sc24_final")
sys.path.insert(0, FTD_DIR)
import ftd_core as F                                       # noqa: E402

OUT = os.path.join(HERE, "out", "dynamo"); os.makedirs(OUT, exist_ok=True)
YR = F.SEC_PER_YEAR


class SXD(F.FTD):
    def __init__(self, alpha0=25.0, B0=1.0e4, kappa=1.0, lat_c=10.0, wlat=45.0, **kw):
        self.kappa, self.lat_c, self.wlat = kappa, lat_c, wlat
        super().__init__(**kw)
        self.alpha0, self.B0 = alpha0, B0
        r = self.r / F.RSUN
        self.Fr = 0.5 * (1.0 + erf((r - self.p["r_dep"]) / self.p["d_dep"]))           # deposition envelope
        lat = self.lat
        self.Wth = np.cos(self.th) * np.sin(self.th) * 0.5 * (1.0 - erf((np.abs(lat) - self.wlat) / 5.0))
        self.tach = (r >= 0.70) & (r <= 0.72)
        self.north = lat > 0; self.south = lat < 0
        w = np.sin(self.th)
        self.capN = np.where(lat >= 70.0, w, 0.0); self.capN /= self.capN.sum()
        self.capS = np.where(lat <= -70.0, w, 0.0); self.capS /= self.capS.sum()
        self.wsin = np.sin(self.th) * self.dth

    # ---- hemispheric coupling: latitudinal diffusivity with an equatorial factor, flux form at theta faces
    def _profiles(self):
        super()._profiles()
        latf = 90.0 - np.degrees(0.5 * (self.th[1:] + self.th[:-1])) if hasattr(self, "th") else None
        fac = 1.0 + (self.kappa - 1.0) * np.exp(-(latf / self.lat_c) ** 2)
        self.facf = fac[None, :]                                   # at theta faces j+1/2

    def _timestep(self):
        super()._timestep()
        if self.kappa > 1.0:                                       # explicit theta diffusion with the larger eta
            self.dt /= self.kappa

    def _lap_theta_eta(self, F_, eta_r):
        sF = self.sin * F_
        G = eta_r[:, None] * self.facf * (sF[:, 1:] - sF[:, :-1]) / (self.dth * np.sin(self.thh)[None, :])
        out = np.zeros_like(F_)
        out[:, 1:-1] = (G[:, 1:] - G[:, :-1]) / (self.dth * self.R2[:, 1:-1] ** 2)
        return out

    def step(self, A, B, S=None):
        """identical to FTD.step except that the latitudinal diffusion uses the coupling factor."""
        if self.kappa == 1.0:
            return super().step(A, B, S)
        dt = self.dt
        slope = self.M_top @ A[-1]
        U = self.s * A
        rhs = self._upwind_A(U, A, slope) + self._lap_theta_eta(A, self.eta_p)
        if S is not None:
            rhs = rhs + S
        Dr = np.zeros_like(A)
        Dr[1:-1] = self.ETA_P[1:-1] * self._lap_r(A)[1:-1]
        Dr[-1] = self.eta_p[-1] * self._lap_r_top(A, slope)
        star = A + dt * rhs + 0.5 * dt * Dr
        star[-1] += 0.5 * dt * self.eta_p[-1] * (2 * self.dr * slope) / (self.dr ** 2) + \
            0.5 * dt * self.eta_p[-1] * (2 * self.dr * slope) / (2 * self.dr) * 2.0 / self.r[-1]
        An = self.TA.solve(star)
        Br, Bth = self.bpol(A)
        shear = self.s * (Br * self.dom_dr + Bth * self.dom_dth / self.R2)
        rB = self.R2 * B
        dgrad = np.zeros_like(B)
        dgrad[1:-1] = (self.deta_t[1:-1, None] * (rB[2:] - rB[:-2]) / (2 * self.dr)) / self.R2[1:-1]
        rhsB = self._adv_B(B) + self._lap_theta_eta(B, self.eta_t) + shear + dgrad
        DrB = np.zeros_like(B)
        DrB[1:-1] = self.ETA_T[1:-1] * self._lap_r(B)[1:-1]
        Bn = self.TB.solve(B + dt * rhsB + 0.5 * dt * DrB)
        An[0] = 0.0; An[:, 0] = 0.0; An[:, -1] = 0.0
        Bn[0] = 0.0; Bn[-1] = 0.0; Bn[:, 0] = 0.0; Bn[:, -1] = 0.0
        return An, Bn

    # ---- the Babcock-Leighton source
    def source(self, B, aN, aS):
        Bbar = B[self.tach].mean(0)
        q = Bbar / (1.0 + (Bbar / self.B0) ** 2)
        alpha = np.where(self.north, aN, np.where(self.south, aS, 0.0))
        return (alpha * self.Wth * q)[None, :] * self.Fr[:, None], Bbar, alpha * self.Wth * q

    def run(self, years, seed=0, sigma=0.0, rho=0.0, tau_c=1.0, rec=0.05, A0=None, B0f=None, seed_asym=0.0, sym=False, swap=False):
        rng = np.random.default_rng(seed)
        if A0 is None:
            A = F.RSUN * np.sin(self.T2) * (self.R2 / F.RSUN) ** 4      # dipole seed, B_r ~ 1 G at the surface
            A *= (1.0 + seed_asym * np.cos(self.T2))                    # optional small quadrupolar seed
            B = np.zeros_like(A)
        else:
            A, B = A0.copy(), B0f.copy()
        nsteps = int(years * YR / self.dt); krec = max(1, int(rec * YR / self.dt)); kren = max(1, int(tau_c * YR / self.dt))
        xiN = xiS = 0.0
        t, EN, ES, PN, PS, BB, XI = [], [], [], [], [], [], []
        t0 = time.time()
        for k in range(nsteps):
            if sigma > 0 and k % kren == 0:
                xiN = rng.normal(); xiS = rho * xiN + np.sqrt(max(0.0, 1 - rho ** 2)) * rng.normal()
                if sym:
                    xiS = xiN
                if swap:                                            # mirror check: the hemispheres trade forcing
                    xiN, xiS = xiS, xiN
            aN = self.alpha0 * (1.0 + sigma * xiN); aS = self.alpha0 * (1.0 + sigma * xiS)
            S, Bbar, s1d = self.source(B, aN, aS)
            A, B = self.step(A, B, S)
            if k % krec == 0:
                Br, _ = self.bpol(A)
                t.append(k * self.dt / YR)
                EN.append(float(np.sum(np.abs(s1d[self.north]) * self.wsin[self.north])))
                ES.append(float(np.sum(np.abs(s1d[self.south]) * self.wsin[self.south])))
                PN.append(float(Br[-1] @ self.capN)); PS.append(float(Br[-1] @ self.capS))
                BB.append(Bbar[::4].astype(np.float32)); XI.append((xiN, xiS))
            if not np.isfinite(A).all():
                raise FloatingPointError("dynamo diverged at step %d" % k)
        return dict(t=np.array(t), EN=np.array(EN), ES=np.array(ES), PN=np.array(PN), PS=np.array(PS),
                    Bbar=np.array(BB), lat=self.lat[::4], xi=np.array(XI), A=A, B=B, wall=time.time() - t0, dt=self.dt)


def cycles(res, smooth_yr=1.0, min_rev_sep=None):
    """Cycle minima anchored to the polar reversals: the axial-dipole proxy P_N - P_S (1-yr running mean) changes
    sign once per cycle, near maximum; the minimum of the smoothed total emergence between two successive reversals
    is the cycle minimum.  (Minima of the emergence alone are ambiguous once the forcing is noisy.)
    When the hemispheres drift out of phase the proxy can cross zero without a reversal, so crossings closer than
    0.6 of the run's own N-polar reversal period are merged (min_rev_sep None = that rule)."""
    t, E = res["t"], res["EN"] + res["ES"]
    dtr = t[1] - t[0]; w = max(1, int(smooth_yr / dtr)); k = np.ones(w) / w
    Es = np.convolve(E, k, mode="same")
    Dp = np.convolve(res["PN"] - res["PS"], k, mode="same")
    if min_rev_sep is None:
        zN = np.flatnonzero(np.diff(np.sign(np.convolve(res["PN"], k, mode="same"))) != 0)
        min_rev_sep = 0.6 * (t[-1] - t[0]) / max(len(zN), 1) if len(zN) > 2 else 4.0
    rev = []
    for i in np.flatnonzero(np.diff(np.sign(Dp)) != 0):
        if w <= i < len(Dp) - w and (not rev or t[i] - t[rev[-1]] >= min_rev_sep):
            rev.append(i)
    mins = [a + int(np.argmin(Es[a:b])) for a, b in zip(rev[:-1], rev[1:])]
    return np.array(mins), Es


def tune():
    out = {}
    for v0 in (1500.0, 2000.0):
        for a0 in (20.0, 40.0):
            m = SXD(alpha0=a0, nr=65, nth=121, v0=v0)
            res = m.run(80.0, rec=0.1)
            mins, Es = cycles(res)
            per = float(np.median(np.diff(res["t"][mins]))) if len(mins) > 2 else float("nan")
            E = res["EN"] + res["ES"]
            amp = [float(E[(res["t"] > a) & (res["t"] < a + 20)].max()) for a in (20, 40, 60)]
            asym = float(np.max(np.abs(res["EN"] - res["ES"])) / np.max(E))
            out[f"v{v0:.0f}_a{a0:.0f}"] = dict(period=per, amp_20_40_60=amp, max_asym=asym, wall=res["wall"], dt_days=m.dt / 86400)
            print(f"v0 {v0/100:.0f} m/s alpha0 {a0}: period {per:.1f} yr, E max over 20-yr blocks {['%.2e' % x for x in amp]}, "
                  f"max |N-S|/E {asym:.1e}, wall {res['wall']:.0f} s, dt {m.dt/86400:.2f} d", flush=True)
    json.dump(out, open(os.path.join(OUT, "tune.json"), "w"), indent=1)


if __name__ == "__main__":
    if sys.argv[1] == "tune":
        tune()
    elif sys.argv[1] == "run":
        import argparse
        ap = argparse.ArgumentParser()
        ap.add_argument("cmd"); ap.add_argument("name"); ap.add_argument("--years", type=float, default=1100.0)
        ap.add_argument("--seed", type=int, default=1); ap.add_argument("--kappa", type=float, default=1.0)
        ap.add_argument("--sigma", type=float, default=0.5); ap.add_argument("--rho", type=float, default=0.0)
        ap.add_argument("--tau_c", type=float, default=1.0); ap.add_argument("--alpha0", type=float, default=30.0)
        ap.add_argument("--v0", type=float, default=1500.0); ap.add_argument("--nr", type=int, default=65)
        ap.add_argument("--nth", type=int, default=121); ap.add_argument("--swap", action="store_true")
        ap.add_argument("--wlat", type=float, default=45.0)
        a = ap.parse_args()
        m = SXD(alpha0=a.alpha0, kappa=a.kappa, nr=a.nr, nth=a.nth, v0=a.v0, wlat=a.wlat)
        res = m.run(a.years, seed=a.seed, sigma=a.sigma, rho=a.rho, tau_c=a.tau_c, swap=a.swap)
        keep = {k: v for k, v in res.items() if k not in ("A", "B")}
        np.savez_compressed(os.path.join(OUT, f"{a.name}.npz"), **keep, meta=json.dumps(vars(a)))
        print(f"{a.name}: done in {res['wall']:.0f} s, dt {res['dt']/86400:.2f} d", flush=True)
