"""ens_transport_series.py -- per-cycle modelled A_pf (deterministic Joy tilts) for the 27 transport configurations of
ens_verify.py, saved for the band in the ensemble figure.  Same models, same readout."""
import numpy as np
from core import *

rows, labels = [], []
for eta in (150.0, 250.0, 400.0):
    for v0 in (10.0, 15.0, 20.0):
        Mb = Model(eta=eta, v0=v0, tau=5.0)
        for tau in (3.0, 5.0, 10.0):
            Mt = Mb if tau == 5.0 else Model(eta=eta, v0=v0, tau=tau, base=Mb)
            rows.append(apf_vec(Mt.run(CAT), CYC)); labels.append(f"eta{eta:.0f}_v{v0:.0f}_tau{tau:.0f}")
        tic(f"eta {eta} v0 {v0} done")
np.savez(os.path.join(OUT, "ens_transport_series.npz"), apf=np.array(rows), labels=np.array(labels), cycles=np.array(CYC))
