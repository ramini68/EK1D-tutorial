# EK1D Tutorial

This repository contains the one-dimensional electrokinetic reactive transport model accompanying the manuscript:

**Electrokinetic Reactive Transport Modeling in Porous Media: A Tutorial and Review of Coupled Flow, Charge Transport, and Hydrogeochemistry**

## Purpose

The model provides a transparent pedagogical example of coupled electrokinetic and reactive transport in a saturated porous medium.

It simulates:

* Nernst-Planck transport of H+, OH-, Cu2+, Na+, and NO3-
* concentration-dependent electrical conductivity
* electroosmotic and hydraulic advection
* current-dependent Faradaic electrode fluxes
* water self-ionization
* pH-dependent copper sorption
* Cu(OH)2 precipitation

The model is intended for teaching and methodological demonstration. It is not calibrated to a specific laboratory or field system.

## Numerical method

The implementation uses:

* one-dimensional finite volumes
* explicit transport integration
* first-order upwinding for migration and advection
* sequential non-iterative transport-reaction coupling
* an ohmic current-continuity formulation under local electroneutrality

## Requirements

* Python 3
* NumPy
* Matplotlib

Install the required packages with:

```bash
python -m pip install -r requirements.txt
```

## Running the model

Run:

```bash
python EK1D_tutorial.py
```

The script performs:

1. An electromigration benchmark
2. A 1.5-day copper-transport simulation
3. Mass- and charge-balance diagnostics
4. Generation of `EK1D_results.png`

## Model assumptions

* saturated and isothermal porous medium
* constant porosity
* local electroneutrality
* no gas-phase transport
* simplified equilibrium copper sorption
* equilibrium Cu(OH)2 precipitation
* pedagogical rather than calibrated parameters

## Repository status

This repository is under development. Numerical verification and charge-balance diagnostics should be reviewed before the model is used for predictive applications.

## License

MIT License
