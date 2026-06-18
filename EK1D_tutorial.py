"""
================================================================================
EK1D_tutorial.py
One-dimensional electrokinetic (EK) reactive transport model — Cu(II) worked example
================================================================================
Accompanies: "Electrokinetic Reactive Transport Modeling: A Tutorial and Survey
for Coupled Flow, Ion Transport, Geochemical Reactions, and Subsurface
Applications" (Section 4.7 / Figure 9 results).

The script is deliberately compact and self-contained (numpy + matplotlib only)
and is organised to follow the TWELVE-STEP modelling workflow of Figure 4:

    1  Problem definition        7  Coupling strategy
    2  Data collection           8  Solver setup
    3  Geometry and mesh          9  Calibration (hook)
    4  Governing equations       10  Validation (benchmark + balances)
    5  Reaction network          11  Scenario simulation
    6  Boundary/initial cond.    12  Sensitivity / uncertainty (hook)

--------------------------------------------------------------------------------
PHYSICS SOLVED
--------------------------------------------------------------------------------
For each mobile aqueous species i in {H+, OH-, Cu2+, Na+, NO3-}, the 1-D
Nernst–Planck reactive transport equation (continuum scale) is

    d(c_i)/dt = -(1/theta) dJ_i/dx + R_i                              (mass conservation)

with the molar flux (per unit bulk cross-sectional area)

    J_i = -Dstar_i dc_i/dx              <- diffusion / dispersion
          + w_i c_i                     <- advection: migration + electroosmosis
    w_i = z_i (F/RT) Dstar_i E + k_eo E + q_h          (effective advective velocity)
    E   = -dphi/dx                                     (electric field)
    Dstar_i = D0_i * (theta / tau)                     (effective diffusion coeff.)

The electric potential phi follows from current continuity under the local
electroneutrality (ohmic) approximation:

    d/dx( sigma dphi/dx ) = 0,   sigma(c) = (F^2/RT) * sum_i z_i^2 Dstar_i c_i
    phi(0) = V (anode),  phi(L) = 0 (cathode)

Electrode water-electrolysis reactions enter as FARADAIC boundary fluxes that
scale with the current density j = mean(-sigma dphi/dx):

    anode  : 2 H2O -> O2 + 4 H+ + 4 e-     -> H+  flux into the domain  = j/F
    cathode: 2 H2O + 2 e- -> H2 + 2 OH-    -> OH- flux into the domain  = j/F

Geochemistry (operator split; local equilibrium each step):
    (a) water self-ionisation   [H+][OH-] = Kw   (proton condition H+-OH- conserved)
    (b) pH-dependent Cu sorption edge (aqueous <-> immobile sorbed)
    (c) Cu(OH)2 precipitation when [Cu2+][OH-]^2 > Ksp

CONVENTIONS / SIMPLIFYING ASSUMPTIONS (state these in any write-up):
  * porosity theta is constant (no clogging feedback on theta in this minimal version);
  * the medium is locally electroneutral, so phi is ohmic (no explicit Poisson);
  * isothermal, saturated, single-phase; gas (O2/H2) is not transported;
  * concentrations are stored in mol/m^3; equilibrium constants use mol/L (note the
    1000x conversion via molL()).

Run:  python EK1D_tutorial.py
================================================================================
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------------------------------------------------ constants
F   = 96485.0      # Faraday constant            [C/mol]
Rg  = 8.314        # universal gas constant      [J/mol/K]
T   = 298.15       # absolute temperature        [K]
FRT = F / (Rg * T) # F/RT                         [1/V]
molL = lambda c: c / 1000.0     # mol/m^3 -> mol/L


# ============================================================================
# STEP 1-2 :  PROBLEM DEFINITION  &  DATA COLLECTION   (cf. Table 6)
# ============================================================================
class P:
    """All model parameters in one place (Table 6 of the manuscript)."""
    # --- geometry & medium ---
    L     = 0.10                 # column length                 [m]
    N     = 50                   # number of finite-volume cells [-]
    theta = 0.45                 # porosity (vol. water content) [-]
    tau   = 1.6                  # tortuosity factor             [-]
    # --- operation ---
    V     = 20.0                 # applied voltage (constant)    [V]   (=> 200 V/m)
    keo   = 1.2e-9               # electroosmotic permeability   [m^2/V/s]
    qh    = 0.0                  # Darcy (hydraulic) flux        [m/s] (closed column)
    # --- species: H+, OH-, Cu2+, Na+, NO3- ---
    names = ['H', 'OH', 'Cu', 'Na', 'NO3']
    z     = np.array([+1, -1, +2, +1, -1], dtype=float)          # charge numbers
    D0    = np.array([9.31e-9, 5.27e-9, 0.714e-9, 1.33e-9, 1.90e-9])  # free-soln D [m^2/s]
    # --- geochemistry ---
    Kw    = 1.0e-14              # water ion product             [(mol/L)^2]
    Ksp   = 2.2e-20             # Cu(OH)2 solubility product    [(mol/L)^3]
    pH50  = 5.5                  # Cu sorption-edge midpoint
    nedge = 2.0                  # Cu sorption-edge sharpness
    # --- initial state ---
    pH0     = 5.0                # initial pore-water pH
    Cu_tot0 = 5.0                # initial TOTAL Cu (aq+sorbed)  [mol/m^3] (0.005 M)
    Na0     = 10.0               # background NaNO3              [mol/m^3] (0.01 M)
    # --- time stepping ---
    dt    = 15.0                 # time step                     [s]
    tend  = 1.5 * 24 * 3600.0    # simulated duration            [s] (1.5 days)


# Derived grid / effective coefficients (STEP 3 + part of STEP 4)
def grid(p=P):
    dx = p.L / p.N
    x  = (np.arange(p.N) + 0.5) * dx          # cell centres [m]
    Dstar = p.D0 * (p.theta / p.tau)          # effective diffusion coeff. [m^2/s]
    return dx, x, Dstar


# ============================================================================
# STEP 5 :  REACTION NETWORK
# ============================================================================
def sorb_frac(pH, p=P):
    """Fraction of total Cu that is SORBED (immobile) as a function of pH.
    Sigmoidal 'sorption edge': ~0 at low pH (acid mobilises Cu), ->1 at high pH."""
    e = np.clip(-p.nedge * (pH - p.pH50), -30, 30)
    return 1.0 / (1.0 + 10.0 ** e)

def react(c, Cu_s, Cu_ppt, p=P):
    """Operator-split local-equilibrium geochemistry applied AFTER transport.
    c        : (5,N) aqueous concentrations [mol/m^3]  (H,OH,Cu,Na,NO3)
    Cu_s     : (N,)  immobile sorbed Cu     [mol/m^3]
    Cu_ppt   : (N,)  precipitated Cu        [mol/m^3]
    Returns updated (c, Cu_s, Cu_ppt, pH)."""
    # (a) water self-ionisation: conserve d = [H+]-[OH-], enforce [H+][OH-]=Kw
    h = molL(c[0]); o = molL(c[1]); d = h - o
    h = (d + np.sqrt(d * d + 4.0 * p.Kw)) / 2.0
    o = np.maximum(h - d, 1e-14)
    c[0] = h * 1000.0
    c[1] = o * 1000.0
    pH = -np.log10(np.maximum(h, 1e-14))

    # (b) Cu sorption edge: repartition the TOTAL Cu pool (aqueous + sorbed)
    pool = c[2] + Cu_s
    fa   = 1.0 - sorb_frac(pH, p)            # aqueous fraction
    c[2] = pool * fa
    Cu_s = pool - c[2]

    # (c) Cu(OH)2 precipitation where supersaturated: bring [Cu2+] to saturation
    cu = molL(c[2]); o = molL(c[1])
    over = cu * o * o > p.Ksp
    if np.any(over):
        cu_eq = p.Ksp / np.maximum(o * o, 1e-30)        # equilibrium [Cu2+] [mol/L]
        ppt   = np.maximum(cu - cu_eq, 0.0) * over       # amount to precipitate
        c[2]  = np.maximum(c[2] - ppt * 1000.0, 0.0)
        Cu_ppt = Cu_ppt + ppt * 1000.0
    return c, Cu_s, Cu_ppt, pH


# ============================================================================
# STEP 4 :  GOVERNING EQUATIONS  —  electric potential (ohmic / electroneutral)
# ============================================================================
def solve_phi(c, dx, Dstar, p=P):
    """Solve  d/dx(sigma dphi/dx)=0  with phi(0)=V, phi(L)=0 by the Thomas algorithm.
    Returns (phi[N], sigma[N]).  sigma is the concentration-dependent conductivity."""
    sigma = (F ** 2 / (Rg * T)) * np.sum((p.z ** 2)[:, None] * Dstar[:, None]
                                         * np.maximum(c, 0.0), axis=0) + 1e-8
    sf = 0.5 * (sigma[:-1] + sigma[1:])          # face conductivities (N-1)
    N = p.N
    a = np.zeros(N); b = np.zeros(N); cc = np.zeros(N); d = np.zeros(N)
    for k in range(N):
        if k > 0:     a[k]  = sf[k-1] / dx; b[k] -= sf[k-1] / dx
        if k < N - 1: cc[k] = sf[k]   / dx; b[k] -= sf[k]   / dx
    # Dirichlet half-cell boundary conditions (phi=V at x=0, phi=0 at x=L)
    gL = 2 * sigma[0]  / dx; gR = 2 * sigma[-1] / dx
    b[0]  -= gL; d[0]  -= gL * p.V
    b[-1] -= gR
    # Thomas (tridiagonal) solve
    cp = np.zeros(N); dp = np.zeros(N)
    cp[0] = cc[0] / b[0]; dp[0] = d[0] / b[0]
    for k in range(1, N):
        m = b[k] - a[k] * cp[k-1]
        cp[k] = cc[k] / m if k < N - 1 else 0.0
        dp[k] = (d[k] - a[k] * dp[k-1]) / m
    phi = np.zeros(N); phi[-1] = dp[-1]
    for k in range(N - 2, -1, -1):
        phi[k] = dp[k] - cp[k] * phi[k+1]
    return phi, sigma


# ============================================================================
# STEP 4 :  GOVERNING EQUATIONS  —  Nernst–Planck transport (one explicit step)
# ============================================================================
def transport_step(c, dx, Dstar, p=P):
    """Advance the aqueous concentrations by ONE explicit finite-volume time step.
    Returns (c, j) with j the current density [A/m^2].
    Scheme: central diffusion + first-order UPWIND for migration/electroosmosis;
            current-tied Faradaic fluxes at the electrode boundary faces."""
    phi, sigma = solve_phi(c, dx, Dstar, p)
    Ef = -(phi[1:] - phi[:-1]) / dx                  # electric field at interior faces
    sf = 0.5 * (sigma[:-1] + sigma[1:])
    j  = float(np.mean(sf * Ef))                     # current density [A/m^2]
    vcap = 0.5 * dx / p.dt                            # CFL velocity cap for stability

    J = np.zeros((5, p.N - 1))
    for i in range(5):
        grad = (c[i, 1:] - c[i, :-1]) / dx
        # effective advective velocity = electromigration + electroosmosis + Darcy
        w = p.z[i] * FRT * Dstar[i] * Ef + p.keo * Ef + p.qh
        w = np.clip(w, -vcap, vcap)                  # cap to keep CFL < 1
        c_up = np.where(w >= 0.0, c[i, :-1], c[i, 1:])   # upwind value
        J[i] = -Dstar[i] * grad + w * c_up

    # flux divergence -> dc/dt (interior); closed boundaries except Faradaic fluxes
    dc = np.zeros((5, p.N))
    dc[:, 1:-1] = -(J[:, 1:] - J[:, :-1]) / dx
    dc[:, 0]    = -(J[:, 0]) / dx
    dc[:, -1]   =  (J[:, -1]) / dx
    ffar = j / F                                     # Faradaic molar flux [mol/m^2/s]
    dc[0, 0]  += ffar / dx                            # H+  produced at the anode
    dc[1, -1] += ffar / dx                            # OH- produced at the cathode

    c = np.maximum(c + p.dt * dc, 0.0)
    return c, j


# ============================================================================
# STEP 6 :  INITIAL & BOUNDARY CONDITIONS
# ============================================================================
def initialise(p=P):
    """Build the initial state, enforcing electroneutrality for the NO3- balance."""
    c = np.zeros((5, p.N))
    c[0] = 10.0 ** (-p.pH0) * 1000.0          # H+  [mol/m^3]  (pH0)
    c[1] = p.Kw / 10.0 ** (-p.pH0) * 1000.0    # OH- [mol/m^3]
    c[3] = p.Na0                               # Na+
    # partition initial Cu by the sorption edge at pH0
    pH = -np.log10(molL(c[0]))
    c[2] = p.Cu_tot0 * (1.0 - sorb_frac(pH, p))
    Cu_s = p.Cu_tot0 - c[2]
    Cu_ppt = np.zeros(p.N)
    # NO3- closes electroneutrality:  sum z_i c_i = 0
    c[4] = c[0] + 2 * c[2] + c[3] - c[1]
    return c, Cu_s, Cu_ppt


# ============================================================================
# STEP 10 :  VALIDATION HELPERS  —  mass/charge balance + analytical benchmark
# ============================================================================
def balances(c, Cu_s, Cu_ppt, dx, p=P):
    """Return (total Cu mass [mol/m^2], max |charge imbalance| [mol/m^3])."""
    cu_mass = np.sum((c[2] + Cu_s + Cu_ppt)) * p.theta * dx
    charge  = np.max(np.abs(np.sum(p.z[:, None] * c, axis=0)))
    return cu_mass, charge

def benchmark_electromigration(p=P):
    """STEP 10 benchmark: pure electromigration of a non-reacting tracer in a uniform
    field has the analytical front velocity  v = z F Dstar E /(RT).
    Run a short transport-only test and compare numerical vs analytical displacement."""
    import copy
    q = copy.deepcopy(p); q.keo = 0.0; q.qh = 0.0; q.N = 200; q.dt = 5.0
    dx, x, Dstar = grid(q)
    # uniform field via a uniform-conductivity background electrolyte (Na+, NO3-)
    c = np.zeros((5, q.N)); c[3] = 10.0; c[4] = 10.0
    # a small Cu2+ pulse near the anode end
    c[2, 5:10] = 1.0; c[4, 5:10] += 2.0          # keep electroneutral
    E = q.V / q.L                                 # ~uniform field [V/m]
    v_analytic = q.z[2] * FRT * Dstar[2] * E      # analytical migration velocity [m/s]
    t = 0.0; tmax = 0.3 * q.L / v_analytic        # travel ~30% of the column
    x0 = np.sum(x * c[2]) / np.sum(c[2])           # initial centroid
    while t < tmax:
        c, _ = transport_step(c, dx, Dstar, q)
        # suppress reactions in the benchmark (transport only)
        t += q.dt
    x1 = np.sum(x * c[2]) / np.sum(c[2])           # final centroid
    v_numeric = (x1 - x0) / t
    return v_analytic, v_numeric, abs(v_numeric - v_analytic) / v_analytic


# ============================================================================
# STEP 7-8-11 :  COUPLING (operator splitting) + SOLVER + SCENARIO RUN
# ============================================================================
def run(p=P, snaps_h=(0, 6, 24, 36), verbose=True):
    dx, x, Dstar = grid(p)
    c, Cu_s, Cu_ppt = initialise(p)
    snaps_t = [int(h * 3600) for h in snaps_h]
    snaps = {}
    t_hist, j_hist = [], []
    cu0, _ = balances(c, Cu_s, Cu_ppt, dx, p)

    def record(tt):
        snaps[tt] = dict(pH=(-np.log10(np.maximum(molL(c[0]), 1e-14))).copy(),
                         Cu_aq=molL(c[2]).copy(), Cu_s=molL(Cu_s).copy(),
                         Cu_ppt=molL(Cu_ppt).copy())
    record(0); nxt = 1
    nsteps = int(p.tend / p.dt)
    for n in range(1, nsteps + 1):
        t = n * p.dt
        # ---- STEP 7: operator splitting = TRANSPORT then REACT ----
        c, j = transport_step(c, dx, Dstar, p)
        c, Cu_s, Cu_ppt, pH = react(c, Cu_s, Cu_ppt, p)
        if n % 40 == 0:
            t_hist.append(t / 3600.0); j_hist.append(j)
        if nxt < len(snaps_t) and t >= snaps_t[nxt]:
            record(snaps_t[nxt]); nxt += 1

    cuF, charge = balances(c, Cu_s, Cu_ppt, dx, p)
    if verbose:
        print("STEP 10  diagnostics")
        print(f"  Cu mass conservation : init {cu0:.4e}  final {cuF:.4e}  "
              f"(rel. error {abs(cuF-cu0)/cu0:.2e})")
        print(f"  max |charge imbalance|: {charge:.3e} mol/m^3")
    return dict(x=x, snaps=snaps, snaps_t=snaps_t,
                t_hist=np.array(t_hist), j_hist=np.array(j_hist))


# ============================================================================
# STEP 12 :  SENSITIVITY (hook) + OUTPUT PLOTTING
# ============================================================================
def plot(res, fname="EK1D_results.png"):
    x = res['x'] * 100.0                          # cm
    st = res['snaps_t']
    lab = {st[0]: '0 h', st[1]: '6 h', st[2]: '24 h', st[3]: '36 h'}
    col = {st[0]: '#999999', st[1]: '#4a6fa5', st[2]: '#c08a2e', st[3]: '#b04a4a'}
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.7))
    for k in st:
        ax[0].plot(x, res['snaps'][k]['pH'], color=col[k], lw=2, label=lab[k])
    ax[0].set_xlabel('Distance from anode (cm)'); ax[0].set_ylabel('pH')
    ax[0].set_title('(a) pH evolution', loc='left', fontweight='bold'); ax[0].legend(frameon=False)
    k = st[-1]
    ax[1].plot(x, res['snaps'][st[0]]['Cu_s']*1000, '--', color='#aaa', lw=1.5, label='sorbed (initial)')
    ax[1].plot(x, res['snaps'][k]['Cu_s']*1000,  color='#4a8a5f', lw=2, label='sorbed')
    ax[1].plot(x, res['snaps'][k]['Cu_aq']*1000, color='#c08a2e', lw=2, label='aqueous')
    ax[1].plot(x, res['snaps'][k]['Cu_ppt']*1000,color='#b04a4a', lw=2, label='precipitated')
    ax[1].set_xlabel('Distance from anode (cm)'); ax[1].set_ylabel('Cu (mol m$^{-3}$)')
    ax[1].set_title('(b) Copper redistribution', loc='left', fontweight='bold'); ax[1].legend(frameon=False)
    ax[2].plot(res['t_hist'], res['j_hist'], color='#2b4a6f', lw=2)
    ax[2].set_xlabel('Time (h)'); ax[2].set_ylabel('Current density (A m$^{-2}$)')
    ax[2].set_title('(c) Current (constant voltage)', loc='left', fontweight='bold')
    for a in ax: a.spines['top'].set_visible(False); a.spines['right'].set_visible(False)
    plt.tight_layout(); plt.savefig(fname, dpi=200, bbox_inches='tight')
    print(f"  saved {fname}")


# ============================================================================
# MAIN
# ============================================================================
if __name__ == "__main__":
    print("STEP 10  benchmark (pure electromigration of a tracer)")
    va, vn, err = benchmark_electromigration()
    print(f"  v_analytic = {va:.3e} m/s   v_numeric = {vn:.3e} m/s   rel.err = {err:.2%}")

    print("\nSTEP 11  scenario simulation (Cu(II) worked example)")
    res = run()
    plot(res)
    print("\nDone.")
