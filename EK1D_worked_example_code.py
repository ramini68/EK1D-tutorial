"""
EK1D_worked_example_code.py
Condensed runner for the one-dimensional Cu(II) electrokinetic reactive-transport
worked example. The full, documented model lives in
EK1D_tutorial.py; this script just runs it, writes the figure, and saves the
snapshot arrays (ek1d_out.npz) used to typeset Figure 7. numpy + matplotlib only.
"""
import numpy as np
import EK1D_tutorial as ek

res = ek.run(ek.P)                      # transport + reaction, with STEP-10 diagnostics
ek.plot(res, "EK1D_results.png")        # Figure 7 panels (pH, Cu, current)

# save snapshots for external plotting / reproducibility
st = res["snaps_t"]
np.savez("ek1d_out.npz",
         x=res["x"], snaps_t=np.array(st),
         **{f"pH_{k}":   res["snaps"][k]["pH"]     for k in st},
         **{f"Cuaq_{k}": res["snaps"][k]["Cu_aq"]  for k in st},
         **{f"Cus_{k}":  res["snaps"][k]["Cu_s"]   for k in st},
         **{f"Cuppt_{k}":res["snaps"][k]["Cu_ppt"] for k in st},
         t_hist=res["t_hist"], I_hist=res["j_hist"])
print("saved ek1d_out.npz")
