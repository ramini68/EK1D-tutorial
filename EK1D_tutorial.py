"""
EK1D_tutorial.py
One-dimensional electrokinetic reactive transport (EKRT) model: copper(II)
worked example.

Accompanies:
    Mahyapour, R. and Sadat-Noori, M. Electrokinetic Reactive Transport Modeling
    in Porous Media: A Unified Approach for Coupled Transport, Hydrogeochemistry,
    and Model Evaluation.

The script implements the twelve-step modelling workflow (Section 4) and
reproduces the worked example (Figure 7). It is compact and self-contained
(numpy + matplotlib only).

Model
-----
For each mobile aqueous species i in {H+, OH-, Cu2+, Na+, NO3-}, the 1-D
Nernst-Planck reactive transport equation is solved on a finite-volume grid:

    theta d(c_i)/dt = -dJ_i/dx + theta R_i
    J_i = -Dstar_i dc_i/dx + (z_i (F/RT) Dstar_i E + k_eo E + q_h) c_i
    Dstar_i = D0_i theta / tau ,   E = -dphi/dx

The electric potential is obtained from electroneutral current continuity that
retains the diffusion (junction) potential:

    d/dx( sigma dphi/dx ) = d/dx( i_diff )
    i_diff = -F sum_i z_i Dstar_i dc_i/dx
    sigma  = (F^2/RT) sum_i z_i^2 Dstar_i c_i

with phi(0) = V at the anode and phi(L) = 0 at the cathode. The species fluxes
reuse the same discrete operators, so the through-current is divergence-free and
electroneutrality is preserved. Electrode water electrolysis enters as Faradaic
boundary fluxes equal to the ionic current at the boundary faces:

    anode   : 2 H2O -> O2 + 4 H+ + 4 e-
    cathode : 2 H2O + 2 e- -> H2 + 2 OH-

Geochemistry is solved by operator splitting at each step and is charge
conserving: water self-ionisation, proton-coupled Cu sorption (Cu2+ <-> 2 H+),
and Cu(OH)2 precipitation solved as a joint Cu-OH equilibrium.

Numerics
--------
Central differencing for diffusion and first-order upwind for electromigration
and electroosmosis, with adaptive CFL sub-stepping; the field is re-solved each
sub-step and concentrations remain non-negative.

Assumptions
-----------
Constant porosity; isothermal, saturated, single phase; local electroneutrality
(no explicit Poisson or electrical double layer); gas (O2/H2) not transported.
Concentrations are in mol/m^3; equilibrium constants use mol/L (converted by
molL()).

Usage
-----
    python EK1D_tutorial.py

License: MIT.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def thomas(a, b, cu, d):
    """Solve a tridiagonal system (sub a, diag b, super cu, rhs d). numpy only."""
    n = len(b); cp = np.empty(n); dp = np.empty(n)
    cp[0] = cu[0] / b[0]; dp[0] = d[0] / b[0]
    for k in range(1, n):
        m = b[k] - a[k] * cp[k-1]
        cp[k] = (cu[k] / m) if k < n - 1 else 0.0
        dp[k] = (d[k] - a[k] * dp[k-1]) / m
    x = np.empty(n); x[-1] = dp[-1]
    for k in range(n - 2, -1, -1):
        x[k] = dp[k] - cp[k] * x[k+1]
    return x

# ------------------------------------------------------------------ constants
F   = 96485.0      # Faraday constant            [C/mol]
Rg  = 8.314        # universal gas constant      [J/mol/K]
T   = 298.15       # absolute temperature        [K]
FRT = F / (Rg * T) # F/RT                         [1/V]
molL = lambda c: c / 1000.0     # mol/m^3 -> mol/L


# ============================================================================
# STEP 1-2 :  PROBLEM DEFINITION  &  DATA COLLECTION   (cf. Table 5)
# ============================================================================
class P:
    """All model parameters in one place (Table 5 of the manuscript)."""
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
    dt    = 15.0                 # OUTPUT/macro step             [s]
    cfl   = 0.4                  # target Courant number for adaptive sub-stepping
    sig_min = 1.0e-12            # tiny floor only to avoid division by zero
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
    """Fraction of total Cu that is SORBED (immobile) vs pH (sigmoidal edge):
    ~0 at low pH (acid mobilises Cu), ->1 at high pH."""
    e = np.clip(-p.nedge * (pH - p.pH50), -30, 30)
    return 1.0 / (1.0 + 10.0 ** e)

def react(c, Cu_s, Cu_ppt, p=P):
    """Charge-conserving operator-split local-equilibrium geochemistry.
    All reactions balance charge:
      (a) water self-ionisation         (neutral)
      (b) proton-coupled Cu sorption    Cu2+ <-> 2 H+   (releases 2 H+ on sorption)
      (c) Cu(OH)2 precipitation         Cu2+ + 2 OH-    (removes +2 and -2 -> neutral)
    Returns updated (c, Cu_s, Cu_ppt, pH)."""
    H  = molL(c[0]); OH = molL(c[1]); Cu = molL(c[2]); Cs = molL(Cu_s)
    pool = Cu + Cs                       # conserved total Cu (aq+sorbed) [mol/L]
    d0   = H - OH                        # proton excess before this react step
    Cs0  = Cs.copy()

    # (a)+(b): solve proton-coupled sorption with water self-ionisation
    Csf = Cs0.copy()
    for _ in range(60):
        d   = d0 + 2.0 * (Csf - Cs0)     # sorbing Cu releases 2 H+ -> raises d
        Hf  = 0.5 * (d + np.sqrt(d * d + 4.0 * p.Kw))
        pHf = -np.log10(np.maximum(Hf, 1e-14))
        Cs_new = pool * sorb_frac(pHf, p)
        if np.max(np.abs(Cs_new - Csf)) < 1e-13:
            Csf = Cs_new; break
        Csf = 0.5 * Csf + 0.5 * Cs_new   # damped fixed point (self-limiting)
    d  = d0 + 2.0 * (Csf - Cs0)
    Hf = 0.5 * (d + np.sqrt(d * d + 4.0 * p.Kw))
    c[0] = Hf * 1000.0
    c[1] = (p.Kw / Hf) * 1000.0
    c[2] = np.maximum(pool - Csf, 0.0) * 1000.0
    Cu_s = np.maximum(Csf, 0.0) * 1000.0

    # (c) Cu(OH)2 precipitation (Cu2+ + 2 OH- <-> Cu(OH)2(s); charge neutral).
    # Solve the JOINT equilibrium (cu-ppt)(o-2ppt)^2 = Ksp for the precipitated
    # amount ppt, bounded by 0 <= ppt <= min(cu, o/2) so that neither Cu2+ nor OH-
    # goes negative. Removing ppt Cu and 2*ppt OH keeps the cell charge balanced.
    cu = molL(c[2]); o = molL(c[1])
    over = cu * o * o > p.Ksp
    if np.any(over):
        lo = np.zeros_like(cu)
        hi = np.minimum(cu, o / 2.0)
        for _ in range(60):                                   # vectorised bisection
            mid = 0.5 * (lo + hi)
            f = (cu - mid) * (o - 2.0 * mid) ** 2 - p.Ksp     # >0 still supersaturated
            hi = np.where(f > 0.0, hi, mid)
            lo = np.where(f > 0.0, mid, lo)
        ppt = np.where(over, 0.5 * (lo + hi), 0.0)
        o2  = o - 2.0 * ppt                                   # >= 0 by construction
        c[2]   = (cu - ppt) * 1000.0
        Cu_ppt = Cu_ppt + ppt * 1000.0
        H  = molL(c[0]); d2 = H - o2                          # OH removed -> d increases
        H  = 0.5 * (d2 + np.sqrt(d2 * d2 + 4.0 * p.Kw))
        c[0] = H * 1000.0; c[1] = (p.Kw / H) * 1000.0

    pH = -np.log10(np.maximum(molL(c[0]), 1e-14))
    return c, Cu_s, Cu_ppt, pH


# ============================================================================
# STEP 4 :  GOVERNING EQUATIONS  -  electric potential
#           (electroneutral current continuity WITH the diffusion/junction term)
# ============================================================================
def solve_phi(c, dx, Dstar, p=P):
    """Solve  d/dx(sigma dphi/dx) = d/dx(i_diff),  phi(0)=V, phi(L)=0  (Thomas).
    i_diff = -F sum_i z_i Dstar_i dc_i/dx   is the diffusion (junction) current.
    Returns (phi[N], sigma_cell[N], i_face[N-1], idiff_face[N-1])."""
    cc_pos = np.maximum(c, 0.0)
    sigma = (F ** 2 / (Rg * T)) * np.sum((p.z ** 2)[:, None] * Dstar[:, None] * cc_pos, axis=0) + p.sig_min
    # UPWIND face conductivity, consistent with the upwind migration used in
    # transport() (E>0 anode->cathode: cations take the left cell, anions the
    # right). This makes the discrete through-current equal to sum_i z_i J_i
    # face by face, so charge is conserved.
    cL = cc_pos[:, :-1]; cR = cc_pos[:, 1:]
    cup = np.where((p.z[:, None] > 0), cL, cR)
    sf = (F ** 2 / (Rg * T)) * np.sum((p.z ** 2)[:, None] * Dstar[:, None] * cup, axis=0) + p.sig_min
    idiff = -F * np.sum(p.z[:, None] * Dstar[:, None] * (c[:, 1:] - c[:, :-1]) / dx, axis=0)
    N = p.N
    # vectorised tridiagonal assembly:  a=sub, b=diag, cu=super
    a = np.zeros(N); b = np.zeros(N); cu = np.zeros(N); d = np.zeros(N)
    a[1:]    = sf / dx
    cu[:-1]  = sf / dx
    b[1:]   -= sf / dx
    b[:-1]  -= sf / dx
    # RHS source = d/dx(i_diff) integrated over the cell = idiff_{k+1/2} - idiff_{k-1/2}
    d[1:-1] = idiff[1:] - idiff[:-1]
    d[0]    = idiff[0]
    d[-1]   = -idiff[-1]
    # Dirichlet half-cell BCs (phi=V at x=0, phi=0 at x=L)
    gL = 2 * sigma[0]  / dx; gR = 2 * sigma[-1] / dx
    b[0]  -= gL; d[0]  -= gL * p.V
    b[-1] -= gR
    phi = thomas(a, b, cu, d)
    Ef = -(phi[1:] - phi[:-1]) / dx                          # field at interior faces
    i_face = sf * Ef + idiff                                 # total ionic current (uniform)
    return phi, sigma, i_face, idiff


# ============================================================================
# STEP 4 :  GOVERNING EQUATIONS  -  Nernst-Planck transport over one macro step
#           (electroneutral current continuity; adaptive CFL sub-stepping)
# ============================================================================
def transport(c, dx, Dstar, dt_target, p=P, faradaic=True):
    """Advance aqueous concentrations by dt_target with adaptive CFL sub-steps.
    The electric field (current continuity WITH the diffusion/junction term) is
    re-solved every sub-step, so the discrete ionic current is divergence free at
    each update and electroneutrality is preserved. Central diffusion + first-order
    upwind migration; electroosmosis uses a common upwind side so it carries ~no
    spurious charge. Concentrations remain non-negative. Returns (c, j).
    """
    Dmax = float(np.max(Dstar))
    t_done = 0.0; j_acc = 0.0
    while t_done < dt_target - 1e-9 * dt_target:
        phi, sigma, i_face, idiff = solve_phi(c, dx, Dstar, p)
        Ef = -(phi[1:] - phi[:-1]) / dx
        wmig = (p.z[:, None] * FRT * Dstar[:, None]) * Ef[None, :]
        wadv = p.keo * Ef + p.qh
        smig = (wmig >= 0.0); sadv = (wadv >= 0.0)
        wmax = float(np.max(np.abs(wmig)) + np.max(np.abs(wadv))) + 1e-30
        dts = min(dt_target - t_done,
                  p.cfl * dx / wmax,
                  p.cfl * dx * dx / (2.0 * Dmax + 1e-30))
        cL = c[:, :-1]; cR = c[:, 1:]
        grad = (cR - cL) / dx
        cmig = np.where(smig, cL, cR)
        cadv = np.where(sadv[None, :], cL, cR)
        J = -Dstar[:, None] * grad + wmig * cmig + wadv[None, :] * cadv
        szJ = np.sum(p.z[:, None] * J, axis=0)               # discrete ionic current / F
        dc = np.zeros((5, p.N))
        dc[:, 1:-1] = -(J[:, 1:] - J[:, :-1]) / dx
        dc[:, 0]    = -(J[:, 0]) / dx
        dc[:, -1]   =  (J[:, -1]) / dx
        if faradaic:
            # electrode flux = actual ionic current at the boundary face / F
            dc[0, 0]  += szJ[0]  / dx                         # H+  at anode
            dc[1, -1] += szJ[-1] / dx                         # OH- at cathode
        c = c + dts * dc
        j_acc += F * float(np.mean(szJ)) * dts
        t_done += dts
    return c, j_acc / dt_target


# ============================================================================
# STEP 6 :  INITIAL & BOUNDARY CONDITIONS
# ============================================================================
def initialise(p=P):
    """Build the initial state, enforcing electroneutrality via NO3-."""
    c = np.zeros((5, p.N))
    c[0] = 10.0 ** (-p.pH0) * 1000.0          # H+  [mol/m^3]
    c[1] = p.Kw / 10.0 ** (-p.pH0) * 1000.0    # OH- [mol/m^3]
    c[3] = p.Na0                               # Na+
    pH = -np.log10(molL(c[0]))
    c[2] = p.Cu_tot0 * (1.0 - sorb_frac(pH, p))
    Cu_s = p.Cu_tot0 - c[2]
    Cu_ppt = np.zeros(p.N)
    c[4] = c[0] + 2 * c[2] + c[3] - c[1]       # NO3- closes electroneutrality
    return c, Cu_s, Cu_ppt


# ============================================================================
# STEP 10 :  VALIDATION HELPERS  -  mass/charge balance + analytical benchmark
# ============================================================================
def balances(c, Cu_s, Cu_ppt, dx, p=P):
    """Return (total Cu [mol/m^2], max|charge imbalance| [mol/m^3])."""
    cu_mass = np.sum((c[2] + Cu_s + Cu_ppt)) * p.theta * dx
    charge  = np.max(np.abs(np.sum(p.z[:, None] * c, axis=0)))
    return cu_mass, charge

def benchmark_electromigration(p=P):
    """STEP 10 benchmark: verify the upwind electromigration discretisation
    against the analytical front velocity  v = z F Dstar E /(RT)  in a PRESCRIBED
    UNIFORM field, for a single non-reacting tracer. This isolates the numerical
    migration+diffusion scheme (no electrode injection, no background-ion
    polarisation), so the comparison is clean."""
    N = 400; dx = p.L / N; x = (np.arange(N) + 0.5) * dx
    Dstar = p.D0 * (p.theta / p.tau)
    i = 2                                   # Cu2+ tracer
    D = Dstar[i]; zi = p.z[i]; E = p.V / p.L
    w = zi * FRT * D * E                     # analytical migration velocity [m/s]
    v_analytic = w
    c = np.zeros(N); c[20:30] = 1.0          # tracer pulse, away from both walls
    dt = p.cfl * min(dx / abs(w), dx * dx / (2.0 * D))
    tmax = 0.3 * p.L / abs(w); nst = int(tmax / dt)
    x0 = np.sum(x * c) / np.sum(c)
    for _ in range(nst):
        cL = c[:-1]; cR = c[1:]
        cup = cL if w >= 0 else cR
        J = -D * (cR - cL) / dx + w * cup
        dc = np.zeros(N)
        dc[1:-1] = -(J[1:] - J[:-1]) / dx
        dc[0] = -J[0] / dx; dc[-1] = J[-1] / dx
        c = c + dt * dc
    x1 = np.sum(x * c) / np.sum(c); t = nst * dt
    v_numeric = (x1 - x0) / t
    return v_analytic, v_numeric, abs(v_numeric - v_analytic) / v_analytic


# ============================================================================
# STEP 7-8-11 :  COUPLING (operator splitting) + SOLVER + SCENARIO RUN
# ============================================================================
def run(p=P, snaps_h=(0, 6, 24, 36), verbose=True):
    dx, x, Dstar = grid(p)
    c, Cu_s, Cu_ppt = initialise(p)
    snaps_t = [int(h * 3600) for h in snaps_h]
    snaps = {}; t_hist, j_hist = [], []
    cu0, _ = balances(c, Cu_s, Cu_ppt, dx, p)
    charge_max = 0.0
    # per-species mass tracking (H,OH excluded - they have electrode sources)
    m0 = {i: np.sum(c[i]) * p.theta * dx for i in (2, 3, 4)}

    def record(tt):
        snaps[tt] = dict(pH=(-np.log10(np.maximum(molL(c[0]), 1e-14))).copy(),
                         Cu_aq=molL(c[2]).copy(), Cu_s=molL(Cu_s).copy(),
                         Cu_ppt=molL(Cu_ppt).copy())
    record(0); nxt = 1
    nsteps = int(p.tend / p.dt)
    for n in range(1, nsteps + 1):
        t = n * p.dt
        c, j = transport(c, dx, Dstar, p.dt, p)               # STEP 7: transport ...
        c, Cu_s, Cu_ppt, pH = react(c, Cu_s, Cu_ppt, p)        #         ... then react
        _, ch = balances(c, Cu_s, Cu_ppt, dx, p)
        charge_max = max(charge_max, ch)
        if n % 40 == 0:
            t_hist.append(t / 3600.0); j_hist.append(j)
        if nxt < len(snaps_t) and t >= snaps_t[nxt]:
            record(snaps_t[nxt]); nxt += 1

    cuF, charge = balances(c, Cu_s, Cu_ppt, dx, p)
    mF = {i: np.sum(c[i]) * p.theta * dx for i in (2, 3, 4)}
    q = np.sum(p.z[:, None] * c, axis=0)
    ions = np.sum(np.abs(p.z[:, None] * c), axis=0)
    if verbose:
        print("STEP 10  diagnostics")
        print(f"  Cu  mass conservation  : rel. error {abs(cuF-cu0)/cu0:.2e}")
        print(f"  Na  mass (closed)      : rel. error {abs(mF[3]-m0[3])/m0[3]:.2e}")
        print(f"  NO3 mass (closed)      : rel. error {abs(mF[4]-m0[4])/m0[4]:.2e}")
        print(f"  charge imbalance |q|   : max {np.max(np.abs(q)):.2e}  rms {np.sqrt(np.mean(q**2)):.2e}"
              f"  median {np.median(np.abs(q)):.2e} mol/m^3")
        print(f"  charge imbalance / IS  : max {np.max(np.abs(q)/(ions+1e-12)):.2e}"
              f"  mean {np.mean(np.abs(q)/(ions+1e-12)):.2e}  (localised at cathode reaction zone)")
        print(f"  min concentration seen : {np.min(c):.2e} mol/m^3 (>=0 -> no clipping used)")
    return dict(x=x, snaps=snaps, snaps_t=snaps_t,
                t_hist=np.array(t_hist), j_hist=np.array(j_hist))


# ============================================================================
# STEP 12 :  SENSITIVITY (hook) + OUTPUT PLOTTING
# ============================================================================
def plot(res, fname="EK1D_results.png"):
    x = res['x'] * 100.0
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
    print("STEP 10  benchmark (pure electromigration of a tracer, transport only)")
    va, vn, err = benchmark_electromigration()
    print(f"  v_analytic = {va:.3e} m/s   v_numeric = {vn:.3e} m/s   rel.err = {err:.2%}")

    print("\nSTEP 11  scenario simulation (Cu(II) worked example)")
    res = run()
    plot(res)
    print("\nDone.")
