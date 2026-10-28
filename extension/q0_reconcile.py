"""q0_reconcile.py -- C1 in its signed-field form vs the unsigned next-cycle result.

C1 written against the SIGNED polar field S = (B_N+B_S)/(|B_N|+|B_S|) needs no parity factor and comes out
negative (every sign reversed).  Because (-1)^(n+1) = -sgn(B_N) at every minimum, X1hat*A_pf = -X1*S identically,
so r(X1, S) and r(X1hat, A_pf) are the same relation.  This script prints both, the first-difference form, and the
pieces that connect C1 to the next-cycle result: magnitudes and signs.
"""
import numpy as np
from core import *

out = {}
x01 = x1v(C1, 0, 1); xh = x1v(C1, 0, 1, True)
apf = np.array([APF[n] for n in C1]); S = np.array([S_SIGNED[n] for n in C1])
out["C1_parity_signed"] = dict(r=r(xh, apf), n=10, signs=signs(xh, apf))
out["C1_signed_field_S"] = dict(r=r(x01, S), n=10, signs=signs(x01, S))
out["C1_signed_field_S_first_differences"] = dict(r=r(np.diff(x01), np.diff(S)), n=9)
out["C1_parity_signed_first_differences"] = dict(r=r(np.diff(xh), np.diff(apf)), n=9)
out["C1_unsigned_no_parity"] = dict(r=r(x01, apf), n=10, signs=signs(x01, apf))
odd = np.array([n % 2 == 1 for n in C1])
out["C1_unsigned_odd"] = dict(r=r(x01[odd], apf[odd]), n=5)
out["C1_unsigned_even"] = dict(r=r(x01[~odd], apf[~odd]), n=5)
out["identity_parity_eq_minus_sgnBN"] = all(parity(n) == -np.sign(T["B"][n][0]) for n in C1)

x14 = x1v(C2, -1, 4); an1 = np.array([AN[n + 1] for n in C2]); an = np.array([AN[n] for n in C2])
out["C2_unsigned"] = dict(r=r(x14, an1), n=12, signs=signs(x14, an1))
out["C2_window_0_1"] = dict(r=r(x1v(C2, 0, 1), an1), n=12)
out["persistence_A(n)_A(n+1)"] = dict(r=r(an, an1), n=12, signs=signs(an, an1))

# the chain X1 -> A_pf -> A(n+1): magnitudes
c9 = list(range(15, 24))
apf9 = np.array([APF[n] for n in c9]); an19 = np.array([AN[n + 1] for n in c9])
out["mag_|Apf(n)|_vs_|A(n+1)|"] = dict(r=r(np.abs(apf9), np.abs(an19)), n=9)
out["mag_|X1(n)|_vs_|A(n+1)|"] = dict(r=r(np.abs(x14), np.abs(an1)), n=12)
out["mag_|X1hat01(n)|_vs_|A(n+1)|"] = dict(r=r(np.abs(x1v(C2, 0, 1)), np.abs(an1)), n=12)
out["mag_|A(n)|_vs_|A(n+1)|"] = dict(r=r(np.abs(an), np.abs(an1)), n=12)
out["signed_Apf(n)_vs_A(n+1)"] = dict(r=r(apf9, an19), n=9, signs=signs(apf9, an19))
out["parity*Apf(n)_vs_A(n+1)"] = dict(r=r([parity(n) * APF[n] for n in c9], an19), n=9,
                                       signs=signs([parity(n) * APF[n] for n in c9], an19))
out["S(n)_vs_A(n+1)"] = dict(r=r([S_SIGNED[n] for n in c9], an19), n=9, signs=signs([S_SIGNED[n] for n in c9], an19))
for k, v in out.items():
    print(f"{k:45s} {v}")
dump("q0_reconcile.json", out)
