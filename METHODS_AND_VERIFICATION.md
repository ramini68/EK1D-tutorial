# EK1D: formulation and numerical verification

This five-species copper example accompanies the EKRT review. Model settings correspond to manuscript Table 6 and selected numerical checks to Supplementary Table S3. The selected Figure 7 dataset is `data/rk2_N400_cfl01`; measured checks and refinement comparisons are recorded in `verification/publication_assessment.json`. Parameters are illustrative, without experimental calibration or validation.

The `data/`, `verification/` and reference-figure files mentioned below are supplied in `EK1D_Numerical_Supplement.zip`. Extract that archive into this code directory to use the retained results and additional tests.

## Governing formulation

The domain is a homogeneous, rigid, saturated column containing H⁺, OH⁻, Cu²⁺, Na⁺ and NO₃⁻. Aqueous concentrations are expressed per pore-water volume; sorbed and precipitated Cu use the same water-volume-equivalent basis. Transport satisfies

$$\theta\frac{\partial c_i}{\partial t}=-\frac{\partial J_i}{\partial x}+\theta R_i,\qquad J_i=-D_i^*\frac{\partial c_i}{\partial x}+\frac{z_iF}{RT}D_i^*E c_i,\qquad D_i^*=\frac{\theta}{\tau}D_i^0.$$

Here $E=-\partial\phi/\partial x$, and positive flux is directed toward increasing $x$. The bulk-area effective diffusivity already contains $\theta$; the finite-volume flux divergence is divided by $\theta$. Phase inventories are $M_k=\theta\int_0^L C_k\,dx$, in mol m⁻² of bulk cross-sectional area.

Impermeable hydraulic endpoints and water continuity require zero net water flux: $q=0$. The hydraulic contribution balances electroosmosis, $q_h=-k_{eo}E$, in a generalized Darcy representation (Tamagnini et al., 2010). Hydraulic permeability is unspecified, so pressure magnitude is not calculated. The retained electroosmotic coefficient does not drive net advection in this closed-flow limit.

Electrical closure uses spatially constant $j=F\sum_i z_iJ_i$ and the imposed electrolyte potential difference (Roy et al., 2023). Centered diffusion and MC-limited migration share face states with the current calculation. A monotone piecewise-linear voltage relation determines current, with a bracketed fallback at rounded zero-field breakpoints. Diffusion current and either local field direction are retained. Electrical endpoint half-cells extend adjacent-cell conductivity over half a cell without a separate diffusive-current contribution. The applied 20 V is an electrolyte potential difference, excluding electrode overpotentials.

For positive current, $J_H(0)=j/F$ and $J_{OH}(L)=-j/F$ represent unit-efficiency Faradaic sources. Cu, Na and nitrate have blocking boundaries. Electrode-source identities switch if current becomes negative. Gas, heat, double-layer storage and reaction-induced changes in water volume are omitted.

Accumulated charge is $Q_A(t)=\int_0^t j(s)\,ds$. It is calculated by method-consistent quadrature over accepted stages. Under positive current, each electrode generates $Q_A/F$ mol m⁻² of H or OH; this is the generated dose, not the acid/base remaining after reaction or the amount of recovered Cu.

## Reaction reconstruction and time integration

Each reaction update conserves $T_{Cu}=c_{Cu}+C_s+C_p$ and $A=H-OH+2c_{Cu}$. Sorption adopts an idealized two-proton exchange: sorbing one Cu releases two aqueous protons, and desorption reverses the exchange. This closure does not resolve surface sites or finite sorption capacity. Water self-ionization and reversible Cu(OH)₂ precipitation preserve the same invariants.

With all concentrations in the following equilibrium equations expressed in mol L⁻¹ and unit activity coefficients,

$$OH=\frac{K_w}{H},\qquad K_s(H)=\left(\frac{10^{-pH_{50}}}{H}\right)^n,\qquad c_{Cu}(H)=\min\left[\frac{T_{Cu}}{1+K_s(H)},\frac{K_{sp}H^2}{K_w^2}\right],$$

$$H-\frac{K_w}{H}+2c_{Cu}(H)-A=0,\qquad C_s=K_s(H)c_{Cu}(H),\qquad C_p=T_{Cu}-c_{Cu}(H)-C_s.$$

A bracketed scalar solve enforces water, sorption and mineral equilibrium jointly, including dissolution of existing precipitate. Failed states are not repaired by concentration clipping or spectator-ion adjustment. Aqueous Cu represents free Cu²⁺ only; hydroxo and other complexes, redox chemistry and host-mineral buffering are omitted. The 2.2 mol m⁻³ sorbed plateau in precipitate-bearing cells follows from $n=2$ and the selected equilibrium constants, not a measured adsorption capacity.

SSPRK2 advances conserved components, with equilibrium reconstruction $\mathcal P$ and conservative transport operator $\mathcal L$:

$$U^{(1)}=\mathcal P\left(U^n+\Delta t\,\mathcal L(U^n)\right),\qquad U^{n+1}=\mathcal P\left[\frac12U^n+\frac12\left(U^{(1)}+\Delta t\,\mathcal L(U^{(1)})\right)\right].$$

The state mixture includes both immobile Cu pools. Electrical coefficients and outgoing-flux positivity bounds are recomputed at each stage. Trials violating the second-stage bound are discarded and retried with a smaller step. SSP theory motivates the construction (Gottlieb et al., 2008); second-order accuracy through moving phase boundaries is not assumed.

## Parameters and initial/boundary conditions

Initial concentrations are spatially uniform. Concentration values below use mol m⁻³ of water or water-equivalent volume unless indicated otherwise.

| Quantity | Setting |
|---|---|
| Length; water content; tortuosity | 0.10 m; $\theta=0.45$; $\tau=1.6$ |
| Species charges, H/OH/Cu/Na/NO₃ order | +1, −1, +2, +1, −1 |
| Free-water diffusivities, same order | (9.31, 5.27, 0.714, 1.33, 1.90) × 10⁻⁹ m² s⁻¹ |
| Temperature; constants | 298.15 K; $F=96485$ C mol⁻¹; $R=8.314$ J mol⁻¹ K⁻¹ |
| Electrical endpoints | $\phi(0)=20$ V; $\phi(0.10)=0$ V; nominal mean field 200 V m⁻¹ |
| Initial pH; H; OH | 5.0; 0.010000; 0.000001 |
| Initial aqueous/sorbed/precipitated Cu | 4.545454545 / 0.4545454545 / 0 |
| Total initial Cu | 5 mol m⁻³; column inventory 0.225 mol m⁻² |
| Initial Na; nitrate | 10; 19.100908091, satisfying initial charge balance |
| Equilibrium constants | $K_w=10^{-14}$ M²; $K_{sp}=2.2\times10^{-20}$ M³ |
| Sorption | $pH_{50}=5.5$; $n=2$; dimensionless $K_s=(10^{-pH_{50}}/H)^n$ with both numerator and denominator in M |
| Hydraulic condition | $q=0$; $q_h=-k_{eo}E$; $k_{eo}=1.2\times10^{-9}$ m² V⁻¹ s⁻¹ |
| Solute endpoints | Cu, Na and nitrate blocked; Faradaic H/OH sources as defined above |
| Selected spatial discretization | 400 cells; $\Delta x=0.25$ mm; centered diffusion, MC-limited migration |
| Selected time integration | SSPRK2; CFL 0.1; maximum step 3.75 s; stage-wise equilibrium and positivity checks |
| Duration; saved samples | 36 h; 60-s instantaneous output; charge quadrature over all accepted stages |

## Verification and acceptance criteria

The packaged driver includes documentation corrections and a diagnostic-plot fix for a fifth saved profile. Its transport, reaction, electrical and time-integration code is unchanged. `SOURCE_UPDATE.json` in the numerical supplement records the exact changes; historical dataset hashes continue to identify their original source files.

Retained verification covers chemistry, transport, electrical closure, Faradaic sources, native/Python agreement, driver integration/export and invalid-input handling. Native/Python agreement checks implementation consistency, not physical validity. An independent DOP853 integration on a smooth unsaturated fixed-grid problem gives SSPRK2 temporal orders 2.00194, 2.00097 and 2.0005. This result does not establish the same order at precipitation fronts.

Tutorial-specific acceptance criteria are:

| Check | Limit |
|---|---|
| Relative conserved-inventory error | 10⁻⁸ |
| Normalized local charge, $\lvert\sum z_ic_i\rvert/\sum\lvert z_ic_i\rvert$ | 10⁻⁵ |
| Relative equilibrium-law residual | 10⁻⁹ |
| Current-continuity / voltage residual | 10⁻⁹ A m⁻² / 10⁻⁹ V |
| Final current difference, spatial / temporal | 5% / 0.1% |
| Integrated current difference, spatial / temporal | 2% / 0.1% |
| Final phase-inventory difference / initial Cu, spatial | 1% |
| Phase-history difference / initial Cu, spatial / temporal | 5% / 0.1% |
| Charge-history difference / refined final charge, spatial / temporal | 2% / 0.1% |
| Final pH = 7 front displacement, spatial | 0.25 mm |
| Sampled-current history relative L1 difference, temporal | 1% |

There is no spatial acceptance criterion for the complete instantaneous-current history or raw precipitate peak. All runs must complete 36 h with finite states, nonnegative transported concentrations, reconciled phase inventories and consistent integrated current/Faradaic input. Positive current is checked for this positive-voltage case.

For the selected run, maximum relative Cu inventory error is $3.03046\times10^{-11}$, normalized local charge residual $1.37762\times10^{-6}$, equilibrium residual $8.37961\times10^{-16}$, current-continuity residual $3.90799\times10^{-14}$ A m⁻² and voltage residual $3.12639\times10^{-13}$ V. Na and nitrate inventory errors are below $4.8\times10^{-12}$.

## Five-run refinement comparison

| Cells | CFL | Final current, A m⁻² | Charge, kC m⁻² | Final precipitate, mol m⁻² | Runtime, min |
|---|---|---|---|---|---|
| 100 | 0.2 | 0.400402 | 130.019 | 0.141195 | 0.699348 |
| 200 | 0.2 | 0.395202 | 128.736 | 0.140678 | 6.19406 |
| 200 | 0.1 | 0.395201 | 128.736 | 0.140678 | 13.2274 |
| 400 | 0.2 | 0.385250 | 127.768 | 0.140212 | 44.2971 |
| 400 | 0.1 | 0.385249 | 127.768 | 0.140212 | 84.1671 |

Each run spans 36 simulated hours. CFL 0.2 and 0.1 use maximum steps of 7.5 and 3.75 s, respectively. Spatial comparisons hold time controls constant; temporal comparisons halve both controls at fixed grid. Runtimes describe the execution environment, not a portable benchmark.

Between 200 and 400 cells, final current changes by 2.5832%, charge by 0.75807%, phase histories by at most 0.72658% of initial Cu, and the final pH = 7 front by 0.0912676 mm. Halving time controls at 400 cells changes final current by 0.00013049%, charge by 0.0000011546%, and phase histories by 0.00086529% of initial Cu. These satisfy the declared output-specific criteria.

Phase-history differences use the maximum over phases and saved times, normalized by initial Cu. Charge-history differences use the refined final charge. Sampled-current L1 differences integrate absolute differences on common 60-s samples; charge comparisons use accepted-stage quadrature. The pH front is linearly interpolated between adjacent cell centers.

## Figures and interpretation limits

`Figure7_Revised.png` shows pH at 0, 6, 24 and 36 h, Cu phase inventories and cumulative charge. At 36 h, aqueous, sorbed and precipitated fractions are 20.963%, 16.720% and 62.317%. `FigureS1_GridSensitivity.png` compares grids; `FigureS2_TimeSensitivity.png` compares time controls. Raw states and current histories remain in the corresponding datasets.

Raw precipitate peaks change from 164.808 to 271.303 mol m⁻³ between 200 and 400 cells. Spatial current-history L1 differences are 6.3294% for 100–200 cells and 7.3965% for 200–400 cells, without monotonic reduction. Neither peak height, band width nor narrow current features are certified as mesh independent. Phase inventories and cumulative charge do not establish pointwise convergence.

The selected explicit calculation required 28,138,187 accepted steps and 19,973,072 retried trials; the mean accepted step was 0.00460584 s. Localized fields impose severe positivity restrictions. Unit activities, restricted speciation and absent buffering further limit physical interpretation. Blocking Cu boundaries mean zero collected Cu recovery. The results demonstrate coupled redistribution and phase transfer, not experimentally validated treatment performance.

## References

- Roy, T., Andrej, J., and Beck, V. A. (2023). A scalable DG solver for the electroneutral Nernst–Planck equations. *Journal of Computational Physics*, 475, 111859. https://doi.org/10.1016/j.jcp.2022.111859
- Gottlieb, S., Ketcheson, D. I., and Shu, C. W. (2008). *High Order Strong Stability Preserving Time Discretizations*. Author manuscript, 17 April 2008. https://www.davidketcheson.info/assets/papers/sspreview.pdf
- Tamagnini, C., Jommi, C., and Cattaneo, F. (2010). A model for coupled electro-hydro-mechanical processes in fine grained soils accounting for gas generation and transport. *Anais da Academia Brasileira de Ciências*, 82, 169–193. https://doi.org/10.1590/S0001-37652010000100014
