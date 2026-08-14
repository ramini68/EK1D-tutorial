# EK1D Tutorial

**A compact, reproducible one-dimensional electrokinetic reactive transport model for copper redistribution and pH-front evolution**

EK1D is an educational electrokinetic reactive transport (EKRT) model developed to accompany the critical review:

> Mahyapour, R., and Sadat-Noori, M. *Electrokinetic reactive transport modeling for contaminant fate and remediation in low-permeability subsurface media: A critical review and evaluation framework.*

The repository illustrates selected formulation, implementation, and verification elements of the review's evidence-informed **twelve-stage EKRT workflow**. It is designed for transparent inspection and modification rather than as a calibrated laboratory, field, or engineering-design model.

## Purpose

The default simulation demonstrates how coupled electrical, transport, and geochemical processes can produce:

* acidic and alkaline fronts generated at the electrodes;
* electromigration of dissolved Cu(II);
* electroosmotic transport of all aqueous species;
* pH-dependent redistribution between aqueous and sorbed copper;
* solubility-controlled formation of Cu(OH)2(s);
* evolving electrical current under constant applied voltage; and
* contaminant re-immobilization before extraction from the domain.

The example is intended to support mechanistic interpretation and reproducible teaching. Its outputs should not be interpreted as quantitative predictions for a specific soil, sediment, aquifer, or remediation system.

## Model formulation

### Mobile aqueous species

The model transports five aqueous species:

* `H+`
* `OH-`
* `Cu2+`
* `Na+`
* `NO3-`

`Na+` and `NO3-` provide the background electrolyte and contribute to conductivity and charge balance.

### Nernst-Planck transport

For each aqueous species `i`, the one-dimensional balance is

```math
\theta\frac{\partial c_i}{\partial t}
= -\frac{\partial J_i}{\partial x}+\theta R_i,
```

with flux

```math
J_i =
-D_{e,i}\frac{\partial c_i}{\partial x}
+\left(
\frac{z_iF}{RT}D_{e,i}E+k_{eo}E+q_h
\right)c_i,
```

where

```math
D_{e,i}=D_{0,i}\frac{\theta}{\tau},
\qquad
E=-\frac{\partial\phi}{\partial x}.
```

The implemented transport model uses the ideal dilute-solution form. Activity-coefficient gradients, aqueous complexation, mechanical dispersion, and species-specific nonideal mobility corrections are not included.

### Electrical closure

The electric potential is obtained from locally electroneutral, nonzero current continuity while retaining the diffusion or junction-current contribution:

```math
\frac{\partial}{\partial x}
\left(\sigma\frac{\partial\phi}{\partial x}\right)
=
\frac{\partial i_{\mathrm{diff}}}{\partial x},
```

with

```math
i_{\mathrm{diff}}
=-F\sum_i z_iD_{e,i}\frac{\partial c_i}{\partial x},
```

and

```math
\sigma
=\frac{F^2}{RT}\sum_i z_i^2D_{e,i}c_i.
```

The potential is fixed at the two ends of the domain. The discrete potential and species-flux operators are constructed consistently so that the ionic through-current is spatially uniform during each transport sub-step. Local electroneutrality residuals introduced by the subsequent operator-split reaction update are monitored and reported rather than assumed to be zero.

### Electrode reactions

Water electrolysis is represented through current-linked Faradaic aqueous boundary fluxes:

```text
Anode:   2 H2O -> O2 + 4 H+ + 4 e-
Cathode: 2 H2O + 2 e- -> H2 + 2 OH-
```

The default implementation assumes unit Faradaic efficiency for these reactions. Oxygen and hydrogen gases are not transported, and electrode overpotential, Butler-Volmer kinetics, competing electrode reactions, and finite reservoir mass balances are not represented.

### Geochemistry

The local reaction update includes:

1. water self-ionization;
2. an empirical pH-dependent, proton-coupled Cu sorption edge; and
3. solubility-controlled Cu(OH)2 precipitation.

The sorption and precipitation stoichiometries are charge-balanced. Cu(OH)2 is added when the aqueous state is supersaturated, using a bounded bisection calculation constrained by `Ksp`. In the current implementation, accumulated precipitate is not redissolved during later undersaturation; this is therefore a pedagogical precipitation treatment rather than a general reversible mineral-equilibrium solver.

## Numerical implementation

The model uses:

* a one-dimensional finite-volume grid;
* central differencing for molecular diffusion;
* first-order upwinding for electromigration and electroosmotic advection;
* adaptive CFL-controlled transport sub-stepping;
* recalculation of the electric field during every transport sub-step;
* sequential non-iterative transport-reaction splitting;
* nonnegative transport updates without routine concentration clipping; and
* NumPy and Matplotlib only.

The internal `STEP` labels in the source code are navigation labels. The repository illustrates selected elements of the manuscript's twelve-stage workflow; it is not a one-to-one implementation of every lifecycle stage, including calibration, independent validation, and uncertainty quantification.

## Default case

| Parameter                          |           Default value |
| ---------------------------------- | ----------------------: |
| Column length                      |                `0.10 m` |
| Finite-volume cells                |                    `50` |
| Volumetric water content, `theta`  |                  `0.45` |
| Geometric tortuosity factor, `tau` |                   `1.6` |
| Temperature                        |              `298.15 K` |
| Applied voltage                    |                  `20 V` |
| Nominal initial voltage gradient   |             `200 V m-1` |
| Electroosmotic coefficient, `k_eo` |     `1.2e-9 m2 V-1 s-1` |
| Prescribed Darcy flux, `q_h`       |               `0 m s-1` |
| Initial pH                         |                   `5.0` |
| Initial total Cu                   | `5 mol m-3` (`0.005 M`) |
| Background NaNO3                   | `10 mol m-3` (`0.01 M`) |
| Cu sorption-edge midpoint          |                `pH 5.5` |
| Cu(OH)2 solubility product         |   `2.2e-20 (mol L-1)^3` |
| Macro time step                    |                  `15 s` |
| Target CFL number                  |                   `0.4` |
| Simulated duration                 |                  `36 h` |

All model parameters are collected in class `P` near the top of `EK1D_tutorial.py`.

## Repository contents

| File                          | Purpose                                                                                                                   |
| ----------------------------- | ------------------------------------------------------------------------------------------------------------------------- |
| `EK1D_tutorial.py`            | Full documented tutorial model, isolated electromigration benchmark, Cu(II) scenario, diagnostics, and figure generation. |
| `EK1D_worked_example_code.py` | Condensed runner that executes the model, writes the figure, and saves snapshot arrays.                                   |
| `EK1D_results.png`            | Reference figure showing pH evolution, copper redistribution, and current.                                                |
| `requirements.txt`            | Minimum Python package requirements.                                                                                      |
| `CITATION.cff`                | Software and manuscript citation metadata.                                                                                |
| `LICENSE`                     | MIT License.                                                                                                              |
| `README.md`                   | Repository description, usage instructions, assumptions, and limitations.                                                 |

Running `EK1D_worked_example_code.py` also creates `ek1d_out.npz`, containing the spatial grid, saved pH and copper snapshots, and current history for independent plotting or analysis.

## Installation

A current Python 3 environment is required.

```bash
git clone https://github.com/ramini68/EK1D-tutorial.git
cd EK1D-tutorial
python -m venv .venv
```

Activate the environment on macOS or Linux:

```bash
source .venv/bin/activate
```

Activate it on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

The current requirements are:

```text
numpy>=1.20
matplotlib>=3.3
```

## Running the model

Run the full tutorial, benchmark, diagnostics, and figure generation:

```bash
python EK1D_tutorial.py
```

Run the condensed worked-example script and save the numerical snapshots:

```bash
python EK1D_worked_example_code.py
```

Generated files:

```text
EK1D_results.png
ek1d_out.npz    # generated by the condensed runner
```

## Verification and diagnostics

The full tutorial reports four types of numerical evidence.

### 1. Isolated electromigration benchmark

A nonreactive Cu(II) tracer is transported in a prescribed uniform electric field. The numerical centroid velocity is compared with

```math
v=\frac{zFD_eE}{RT}.
```

This test verifies the isolated migration-diffusion discretization. It does not, by itself, verify the complete coupled EKRT model.

### 2. Component-inventory checks

The script reports relative inventory residuals for:

* total Cu across aqueous, sorbed, and precipitated pools;
* Na in the closed species balance; and
* NO3 in the closed species balance.

### 3. Electrical consistency diagnostics

The script reports:

* maximum, root-mean-square, and median local charge residual;
* charge residual normalized by local ionic charge concentration; and
* the spatially uniform ionic current produced by the transport-step potential solve.

### 4. Positivity diagnostic

The minimum aqueous concentration encountered during the run is reported. The transport algorithm does not use routine post-update concentration clipping.

A representative run of the current code produces output close to:

```text
STEP 10  benchmark (pure electromigration of a tracer, transport only)
  v_analytic = 3.127e-06 m/s   v_numeric = 3.127e-06 m/s   rel.err = 0.00%

STEP 11  scenario simulation (Cu(II) worked example)
STEP 10  diagnostics
  Cu  mass conservation  : rel. error 5.55e-15
  Na  mass (closed)      : rel. error 3.21e-15
  NO3 mass (closed)      : rel. error 2.58e-15
  charge imbalance |q|   : max 3.30e+00  rms 5.15e-01  median 1.31e-07 mol/m^3
  charge imbalance / IS  : max 2.59e-02  mean 1.35e-03
                            (localised at cathode reaction zone)
  min concentration seen : 5.15e-15 mol/m^3 (>=0 -> no clipping used)
```

Small floating-point differences may occur across Python, NumPy, and operating-system versions.

## Interpretation of the default output

The simulation predicts inward-moving acidic and alkaline fronts. Acidification releases part of the initially sorbed copper near the anode. Dissolved Cu(II) then migrates toward the cathode, where alkaline conditions promote renewed sorption and Cu(OH)2 formation. The current evolves because ionic concentrations and calculated conductivity change with time.

The main mechanistic result is that **mobilization is not equivalent to extraction**. Copper can be mobilized from one part of the domain and subsequently re-immobilized elsewhere before reaching a collection boundary.

## Assumptions and limitations

The repository is intentionally compact. The following limitations should be considered when interpreting or extending it:

* one-dimensional, saturated, isothermal, single-phase porous medium;
* constant porosity, tortuosity, electroosmotic coefficient, and temperature;
* ideal dilute-solution transport without activity corrections or aqueous complexes;
* local bulk electroneutrality; no Poisson-resolved space charge or electrical double layer;
* fixed applied voltage and concentration-dependent Nernst-Einstein conductivity;
* no surface conductivity or solid-phase electrical conduction;
* no mechanical dispersion;
* no explicit hydraulic-pressure solution or electroosmotic back-pressure;
* `q_h = 0` denotes no imposed Darcy flux, while electroosmotic advection is retained constitutively;
* empirical pH-dependent Cu sorption rather than a calibrated surface-complexation model;
* precipitation accumulates under supersaturation but accumulated solid does not redissolve;
* unit Faradaic efficiency, no electrode kinetics, no competing electrode reactions, and no gas transport;
* no finite electrode-reservoir mass balance;
* sequential operator splitting introduces a localized electroneutrality residual near the cathode reaction zone;
* the included analytical benchmark isolates the transport discretization rather than the complete coupled formulation;
* no bundled mesh-refinement, macro-time-step refinement, parameter-calibration, independent-validation, or uncertainty-analysis study; and
* no calibration to a specific experiment, material, or field site.

The model should therefore be treated as an inspectable tutorial and starting point. Predictive or design use requires application-specific chemistry, boundary conditions, hydraulic treatment, numerical verification, calibration, independent validation, and uncertainty assessment.

## Reproducibility

The simulation is deterministic for a fixed software environment and parameter set. For reproducible use:

1. record the repository release tag or commit hash;
2. retain the Python and dependency versions;
3. preserve any changes to class `P` and the reaction functions;
4. archive the generated figure and `ek1d_out.npz`; and
5. report all code modifications when publishing derived results.

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). Cite both the software version used and the accompanying manuscript.

Before creating the archival software release, update `CITATION.cff` so that its title, version, release date, and preferred article citation match the final manuscript and tagged repository version.

## License

This project is distributed under the MIT License. See [`LICENSE`](LICENSE).

## Contact

**Ramin Mahyapour**
College of Science and Engineering, James Cook University
Email: [ramin.mahyapour@my.jcu.edu.au](mailto:ramin.mahyapour@my.jcu.edu.au)
ORCID: 0000-0001-5025-2132
