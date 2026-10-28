"""q0_candidates.py -- which simple relation gives r ~ -0.98 with a sign reversal?

The brief refers to an "earlier C1 result (r ~ -0.98, sign reversal)".  No file in the project stores such a
number, so this script computes every plausible candidate (predictor or polar asymmetry against the next-cycle
area asymmetry, signed and unsigned, full sample and parity subsets) and prints them with n.
"""
import os, sys, json
import numpy as np
for _d in ("master_code", "pipeline"):                      # pipeline/ is master_code/ in the repository
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", _d))
import asymcore as A

T = A.targets(); ann = T["ann"]
r = lambda x, y: float(np.corrcoef(x, y)[0, 1])
out = []

def add(label, cyc, x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    out.append((label, len(cyc), r(x, y), int((np.sign(x) == np.sign(y)).sum()), f"{cyc[0]}-{cyc[-1]}"))

C2 = A.C2; C9 = list(range(15, 24))                       # pairs with a measured A_pf: 15->16 ... 23->24
An1 = lambda cyc: [T["An"][n + 1] for n in cyc]
odd = lambda cyc: [n for n in cyc if n % 2]
even = lambda cyc: [n for n in cyc if n % 2 == 0]

for (s, e) in [(0, 1), (-1, 4)]:
    for signed in (False, True):
        for lab, cyc in (("all", C2), ("odd", odd(C2)), ("even", even(C2))):
            add(f"X1[{s:+d},{e:+d}] {'signed' if signed else 'unsigned'} vs A(n+1) [{lab}]", cyc,
                A.x1_vector(cyc, s, e, signed, ann), An1(cyc))

for lab, cyc in (("all", C9), ("odd", odd(C9)), ("even", even(C9))):
    add(f"A_pf(n) vs A(n+1) [{lab}]", cyc, [T["Apf"][n] for n in cyc], An1(cyc))
    add(f"parity*A_pf(n) vs A(n+1) [{lab}]", cyc, [A.parity(n) * T["Apf"][n] for n in cyc], An1(cyc))
    add(f"S(n) vs A(n+1) [{lab}]", cyc, [T["S"][n] for n in cyc], An1(cyc))
    add(f"A_pf(n) vs A(n) same cycle [{lab}]", cyc, [T["Apf"][n] for n in cyc], [T["An"][n] for n in cyc])
    add(f"parity*A_pf(n) vs A(n) [{lab}]", cyc, [A.parity(n) * T["Apf"][n] for n in cyc], [T["An"][n] for n in cyc])

C10 = A.C1
for lab, cyc in (("all", C10), ("odd", odd(C10)), ("even", even(C10))):
    add(f"X1[+0,+1] unsigned vs A_pf(n) [{lab}]", cyc, A.x1_vector(cyc, 0, 1, False, ann), [T["Apf"][n] for n in cyc])
    add(f"X1[+0,+1] signed vs A_pf(n) [{lab}]", cyc, A.x1_vector(cyc, 0, 1, True, ann), [T["Apf"][n] for n in cyc])
    add(f"X1[+0,+1] unsigned vs S(n) [{lab}]", cyc, A.x1_vector(cyc, 0, 1, False, ann), [T["S"][n] for n in cyc])

for lab, cyc in (("all", C2), ("odd", odd(C2)), ("even", even(C2))):
    add(f"A(n) vs A(n+1) persistence [{lab}]", cyc, [T["An"][n] for n in cyc], An1(cyc))

out.sort(key=lambda t: t[2])
print(f"{'relation':58s} {'n':>3s} {'r':>7s} signs  cycles")
for lab, n, rv, sg, cy in out:
    print(f"{lab:58s} {n:3d} {rv:+7.3f} {sg:2d}/{n:<2d} {cy}")
json.dump([dict(relation=a, n=b, r=c, signs=d, cycles=e) for a, b, c, d, e in out],
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "q0_candidates.json"), "w"), indent=1)
