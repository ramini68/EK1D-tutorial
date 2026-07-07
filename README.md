# EK1D-tutorial

A compact, self-contained one-dimensional **electrokinetic reactive transport (EKRT)**
model that accompanies the review

> Mahyapour, R. and Sadat-Noori, M. *Electrokinetic Reactive Transport Modeling in
> Porous Media: A Unified Approach for Coupled Transport, Hydrogeochemistry, and
> Model Evaluation.* (submitted).

The script reproduces the worked example (copper transport and pH-front evolution)
and the corresponding figure, and demonstrates the twelve-step workflow described
in Section 4 of the paper.

## What it does

For five aqueous species (H+, OH-, Cu2+, Na+, NO3-) the code solves, on a
finite-volume grid:

* the **Nernst-Planck** flux (diffusion + electromigration + electroosmosis);
* the electric potential from **electroneutral current continuity including the
  diffusion (junction) potential**, so the discrete ionic current is
  divergence-free and charge is conserved;
* **Faradaic electrode boundaries** (anode: 2H2O -> O2 + 4H+ + 4e-;
  cathode: 2H2O + 2e- -> H2 + 2OH-), with the electrode flux set equal to the
  ionic current at the boundary face;
* operator-split, charge-conserving geochemistry: water self-ionisation,
  **proton-coupled Cu sorption** (Cu2+ <-> 2H+), and Cu(OH)2 precipitation solved
  as a joint Cu-OH equilibrium.

Numerics: first-order upwind electromigration with adaptive CFL sub-stepping
(no velocity cap, no concentration clipping); the electric field is re-solved each
sub-step. Dependencies are **numpy + matplotlib only**.

## Files

| File | Purpose |
|------|---------|
| `EK1D_tutorial.py` | Full, documented model. Runs a verification benchmark, the Cu(II) scenario, prints STEP-10 diagnostics, and writes `EK1D_results.png`. |
| `EK1D_worked_example_code.py` | Condensed runner: imports the tutorial, produces the figure and saves snapshot arrays (`ek1d_out.npz`). |
| `EK1D_results.png` | Figure produced by the model (pH fronts, Cu redistribution, current). |
| `requirements.txt`, `LICENSE`, `.gitignore` | Environment, license (MIT), housekeeping. |

## Install and run

```bash
python -m venv venv && source venv/bin/activate   # optional
pip install -r requirements.txt
python EK1D_tutorial.py
```

Expected console output (verification and conservation diagnostics):

```
STEP 10  benchmark (pure electromigration of a tracer, transport only)
  v_analytic = 3.127e-06 m/s   v_numeric = 3.127e-06 m/s   rel.err = 0.00%
STEP 10  diagnostics
  Cu  mass conservation  : rel. error ~1e-15
  Na / NO3 mass (closed) : rel. error ~1e-15
  charge imbalance |q|   : max ~3 mol/m^3  (localised at the cathode reaction zone)
  min concentration seen : >= 0  (no clipping used)
```

All model parameters are collected in the class `P` at the top of
`EK1D_tutorial.py` (Table 5 of the paper) and can be edited directly.

## Notes / scope

The example is **educational and reproducible**, not a calibrated site model. It is
intended to illustrate model construction, verification (analytical benchmark,
Faraday's law, mass and charge balance), and the interpretation of coupled pH-front
and metal-redistribution behaviour. The small residual charge imbalance is confined
to the cathode reaction zone and stems from sequential reaction operator splitting,
a limitation discussed in the paper.

## License

MIT (see `LICENSE`).
