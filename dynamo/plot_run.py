"""plot_run.py -- butterfly diagrams and a north-south check for one Surya run.

    python plot_run.py [RUNDIR] [--out figures]

Reads rad.dat (surface B_r), butbot.dat (toroidal field B_phi at r = 0.706 R), series_b and ber.dat from RUNDIR (default
run/), writes butterfly.png / butterfly.pdf and prints the hemispheric asymmetry of both fields,
A(t) = (sum_N |B| - sum_S |B|) / (sum_N |B| + sum_S |B|), the same form as A_pf in the observational analysis.
"""
import argparse, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 256                                         # grid intervals in surya_2.f (nmax = 257)
YR = 1.0e8 / (3600 * 24 * 365)                  # code time unit (1e8 s) in years, 'fac' in surya_2.f


def load_map(path):
    """t (code units) and a (time, colatitude) array; Surya writes j = 2..n for every output time."""
    d = np.loadtxt(path)
    nj = N - 1
    nt = len(d) // nj
    d = d[:nt * nj].reshape(nt, nj, 3)
    return d[:, 0, 0], d[0, :, 1], d[:, :, 2]


def ns_asymmetry(f):
    """(N - S)/(N + S) of sum |f| at each time; column k and column nj-1-k are mirror points about the equator."""
    nj = f.shape[1]
    s, n = np.abs(f[:, :nj // 2]).sum(1), np.abs(f[:, nj - nj // 2:]).sum(1)     # q runs from pi (south) to 0
    return (n - s) / np.where(n + s > 0, n + s, 1)


def parity_residual(f, sign):
    """max |f(theta) - sign * f(pi - theta)| / max |f|: 0 for a field of pure parity."""
    return np.abs(f - sign * f[:, ::-1]).max() / np.abs(f).max()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir", nargs="?", default="run")
    ap.add_argument("--out", default="figures")
    a = ap.parse_args()

    tb, q, bphi = load_map(os.path.join(a.rundir, "butbot.dat"))
    tr, _, br = load_map(os.path.join(a.rundir, "rad.dat"))
    lat = 90 - np.degrees(q)

    print(f"{len(tb)} snapshots, t = {tb[0] * YR:.2f} .. {tb[-1] * YR:.2f} yr")
    for name, f, sgn in (("B_phi(0.706 R)", bphi, -1), ("B_r(surface)", br, -1)):    # both odd about the equator
        A = ns_asymmetry(f)
        print(f"  {name:15s} max |A_NS| = {np.abs(A).max():.2e}   dipolar-parity residual = {parity_residual(f, sgn):.2e}")
    ser = os.path.join(a.rundir, "series_b")
    if os.path.exists(ser):
        s = np.loadtxt(ser)
        print(f"  |B_phi| at r = 0.70 R, lat 14.8 N: mean {s[:, 1].mean():.3f}, max {s[:, 1].max():.3f} (code units)")
    ber = os.path.join(a.rundir, "ber.dat")
    if os.path.exists(ber):
        with open(ber) as fh:                   # written with a Fortran D exponent
            e = np.atleast_2d(np.loadtxt([l.replace("D", "E") for l in fh]))
        print(f"  buoyant eruptions: {len(e)} ({(np.cos(e[:, 1]) > 0).sum()} north, {(np.cos(e[:, 1]) < 0).sum()} south)")

    os.makedirs(a.out, exist_ok=True)
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 5.4), sharex=True, constrained_layout=True)
    for axi, t, f, title in ((ax[0], tb, bphi, r"toroidal field $B_\phi$ at $r = 0.706\,R_\odot$"),
                             (ax[1], tr, br, r"surface radial field $B_r$")):
        v = np.percentile(np.abs(f), 99)
        m = axi.pcolormesh(t * YR, lat, f.T, cmap="RdBu_r", vmin=-v, vmax=v, shading="auto", rasterized=True)
        axi.set_ylabel("latitude (deg)"); axi.set_yticks([-90, -45, 0, 45, 90]); axi.set_title(title, fontsize=9, loc="left")
        fig.colorbar(m, ax=axi, pad=0.01, label="code units")
    ax[1].set_xlabel("time (yr)")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(a.out, f"butterfly.{ext}"), dpi=200)
    print(f"wrote {a.out}/butterfly.png, .pdf")


if __name__ == "__main__":
    main()
