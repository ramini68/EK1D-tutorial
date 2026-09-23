"""Joint water--copper sorption--Cu(OH)2 equilibrium for the EK1D tutorial.

Inputs and returned concentrations use mol/m^3 pore water.  Equilibrium
constants use mol/L as in the original tutorial.  The empirical sorption law
is retained; it is an idealised proton-exchange closure, not a surface
complexation model.  No spectator-ion correction or concentration clipping is
used to enforce charge conservation.
"""
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ChemistryParameters:
    Kw: float = 1.0e-14
    Ksp: float = 2.2e-20
    pH50: float = 5.5
    nedge: float = 2.0


@dataclass(frozen=True)
class EquilibriumResult:
    H: np.ndarray
    OH: np.ndarray
    Cu: np.ndarray
    Cu_s: np.ndarray
    Cu_ppt: np.ndarray
    pH: np.ndarray
    iterations: int
    max_charge_residual: float
    max_copper_residual: float


def _water_h(d, kw):
    """Positive root of H - Kw/H = d, avoiding cancellation for d < 0."""
    r = np.hypot(d, 2.0 * np.sqrt(kw))
    # Only evaluate the chosen branch: the unused negative-d numerator would
    # lose significance, whereas evaluating both via np.where can overflow.
    out = np.empty_like(d)
    positive = d >= 0.0
    out[positive] = 0.5 * d[positive] + 0.5 * r[positive]
    out[~positive] = kw / (0.5 * r[~positive] - 0.5 * d[~positive])
    return out


def equilibrate(total_cu, proton_cu_charge, params=ChemistryParameters(), initial_H=None):
    """Return the unique joint equilibrium at fixed Cu and aqueous charge.

    ``total_cu = Cuaq + Cs + Cp`` and
    ``proton_cu_charge = H - OH + 2*Cuaq`` are in mol/m^3.  Na and NO3 are
    passive in this local calculation, so conserving the latter also
    conserves full aqueous charge when their concentrations are unchanged.

    The scalar unknown y = ln(H / (mol/L)) solves the strictly increasing
    equation H - Kw/H + 2*Cuaq(H) - A = 0, with
      Cuaq(H) = min(T / [1 + (10^(-pH50)/H)^nedge], Ksp*H^2/Kw^2).
    This simultaneously permits precipitation AND dissolution.  A bounded
    Newton method in log-H prevents negative concentrations and falls back
    to bisection when a Newton step leaves its root bracket.

    ChemistryParameters or an object exposing Kw, Ksp, pH50 and nedge may be
    supplied. ``initial_H`` optionally provides a previous/transported aqueous
    H concentration in mol/m^3 to warm-start the same bounded solve. Invalid
    hints use the cold-start bracket midpoint; hints never alter invariants
    or tolerances. Failed convergence raises RuntimeError; it is never ignored.
    """
    total, charge = np.broadcast_arrays(
        np.asarray(total_cu, dtype=float),
        np.asarray(proton_cu_charge, dtype=float),
    )
    if np.any(~np.isfinite(total)) or np.any(~np.isfinite(charge)):
        raise ValueError("Equilibrium totals must be finite.")
    if np.any(total < 0.0):
        raise ValueError("Total copper must be nonnegative.")
    if not (params.Kw > 0.0 and params.Ksp > 0.0 and params.nedge > 0.0):
        raise ValueError("Kw, Ksp and sorption exponent must be positive.")
    if not np.all(np.isfinite([params.Kw, params.Ksp, params.pH50, params.nedge])):
        raise ValueError("Equilibrium parameters must be finite.")

    shape = total.shape
    t = total.ravel() / 1000.0
    a = charge.ravel() / 1000.0
    kw, n = float(params.Kw), float(params.nedge)
    log_kw = np.log(kw)
    log_ksp = np.log(float(params.Ksp))
    log_h50 = -float(params.pH50) * np.log(10.0)
    log_t = np.full_like(t, -np.inf)
    np.log(t, out=log_t, where=t > 0.0)

    # Since 0 <= Cuaq <= T, these stable water roots rigorously bracket H.
    lower = np.log(_water_h(a - 2.0 * t, kw))
    upper = np.log(_water_h(a, kw))
    y = 0.5 * (lower + upper)
    if initial_H is not None:
        hint = np.broadcast_to(np.asarray(initial_H, dtype=float), shape).ravel()
        valid_hint = np.isfinite(hint) & (hint > 0.0)
        log_hint = y.copy()
        # Log before unit conversion avoids underflow for subnormal hints.
        np.log(hint, out=log_hint, where=valid_hint)
        log_hint[valid_hint] -= np.log(1000.0)
        y = np.where(valid_hint, np.minimum(upper, np.maximum(lower, log_hint)), y)
    previous_step = upper - lower
    eps = np.finfo(float).eps

    def evaluate(log_h):
        h = np.exp(log_h)
        oh = np.exp(log_kw - log_h)
        u = n * (log_h - log_h50)
        # Stable complementary logistic fractions, without an exponent cap.
        ex = np.exp(-np.abs(u))
        den = 1.0 + ex
        aq_fraction = np.where(u >= 0.0, 1.0 / den, ex / den)
        sorb_fraction = np.where(u >= 0.0, ex / den, 1.0 / den)
        aq_unsat = t * aq_fraction
        sorb_unsat = t * sorb_fraction
        log_aq_sat = log_ksp + 2.0 * log_h - 2.0 * log_kw
        # Compare in logarithms; only exponentiate a saturated concentration
        # where the dissolved+sorbed saturated pool is smaller than T.
        log_aq_unsat = log_t - np.logaddexp(0.0, -u)
        saturated = (t > 0.0) & (log_aq_sat < log_aq_unsat)
        aq = aq_unsat.copy()
        sorbed = sorb_unsat.copy()
        aq[saturated] = np.exp(log_aq_sat[saturated])
        sorbed[saturated] = np.exp(log_aq_sat[saturated] - u[saturated])
        precipitated = np.zeros_like(t)
        precipitated[saturated] = (
            t[saturated] - (aq[saturated] + sorbed[saturated])
        )
        # At an exactly coincident phase boundary, floating-point exp/log
        # evaluation can leave a negative remainder of order roundoff.  Use
        # the mathematically equivalent unsaturated branch for that cell.
        at_boundary = saturated & (precipitated < 0.0)
        if np.any(at_boundary):
            if np.any(precipitated[at_boundary] < -64.0 * eps * t[at_boundary]):
                raise RuntimeError("Inconsistent phase selection.")
            aq[at_boundary] = aq_unsat[at_boundary]
            sorbed[at_boundary] = sorb_unsat[at_boundary]
            precipitated[at_boundary] = 0.0
            saturated[at_boundary] = False
        residual = (h - oh) + 2.0 * aq - a
        derivative = h + oh + np.where(
            saturated, 4.0 * aq, 2.0 * n * aq * sorb_fraction
        )
        scale = h + oh + 2.0 * aq + np.abs(a)
        return h, oh, aq, sorbed, precipitated, residual, derivative, scale

    for iteration in range(1, 81):
        values = evaluate(y)
        h, oh, aq, sorbed, precipitated, residual, derivative, scale = values
        # The absolute term corresponds to 1e-21 mol/m^3.  The relative
        # allowance accommodates cancellation in a physically large charge
        # pool; actual residuals remain available to the caller.
        converged = np.abs(residual) <= 64.0 * eps * scale + 1.0e-24
        if np.all(converged):
            break
        lower = np.where((residual < 0.0) & ~converged, y, lower)
        upper = np.where((residual > 0.0) & ~converged, y, upper)
        proposal = y - residual / derivative
        outside = (~np.isfinite(proposal)) | (proposal <= lower) | (proposal >= upper)
        # A noncontracting Newton cycle can straddle a phase boundary.
        # Require its step to contract, otherwise bisect the bracket.
        slow = np.abs(proposal - y) > 0.5 * previous_step
        proposal = np.where(outside | slow, 0.5 * (lower + upper), proposal)
        previous_step = np.where(converged, previous_step, np.abs(proposal-y))
        y = np.where(converged, y, proposal)
    else:
        raise RuntimeError(
            "Joint equilibrium did not converge in 80 iterations; "
            f"max |charge residual| = {np.max(np.abs(residual))*1000.0:.6g} mol/m^3; "
            f"worst T,A,y,lo,hi: {t[np.argmax(abs(residual))]}, {a[np.argmax(abs(residual))]}, "
            f"{y[np.argmax(abs(residual))]}, {lower[np.argmax(abs(residual))]}, {upper[np.argmax(abs(residual))]}"
        )

    copper_residual = aq + sorbed + precipitated - t
    if np.any(np.abs(copper_residual) > 64.0 * eps * t + 1.0e-18):
        raise RuntimeError("Joint equilibrium did not conserve total copper.")
    def concentration(value):
        return (1000.0 * value).reshape(shape)
    return EquilibriumResult(
        H=concentration(h), OH=concentration(oh), Cu=concentration(aq),
        Cu_s=concentration(sorbed), Cu_ppt=concentration(precipitated),
        pH=(-y / np.log(10.0)).reshape(shape), iterations=iteration,
        max_charge_residual=float(np.max(np.abs(residual), initial=0.0) * 1000.0),
        max_copper_residual=float(np.max(np.abs(copper_residual), initial=0.0) * 1000.0),
    )


def react(c, Cu_s, Cu_ppt, p=ChemistryParameters()):
    """Original tutorial-compatible reaction wrapper; preserves spectator ions."""
    aqueous = np.array(c, dtype=float, copy=True)
    if aqueous.shape[0] != 5:
        raise ValueError("Expected species order H, OH, Cu, Na, NO3.")
    if np.any(~np.isfinite(aqueous)) or np.any(aqueous < 0.0):
        raise ValueError("Aqueous concentrations must be finite and nonnegative.")
    sorbed = np.asarray(Cu_s, dtype=float)
    ppt = np.asarray(Cu_ppt, dtype=float)
    if np.any(~np.isfinite(sorbed)) or np.any(~np.isfinite(ppt)):
        raise ValueError("Solid concentrations must be finite.")
    if np.any(sorbed < 0.0) or np.any(ppt < 0.0):
        raise ValueError("Solid concentrations must be nonnegative.")
    total = aqueous[2] + sorbed + ppt
    invariant = aqueous[0] - aqueous[1] + 2.0 * aqueous[2]
    eq = equilibrate(total, invariant, p, initial_H=aqueous[0])
    aqueous[0], aqueous[1], aqueous[2] = eq.H, eq.OH, eq.Cu
    return aqueous, eq.Cu_s, eq.Cu_ppt, eq.pH
