"""Optional C++17 acceleration of the unchanged corrected EK1D equations.

Build uses a standard C++ compiler and ctypes; NumPy is the sole Python
dependency. No fast-math optimizations are enabled. The Python implementation
remains an independently callable reference. Each call advances one interval
and returns copies of the new states plus diagnostics over every accepted step.
"""
from pathlib import Path
import ctypes
import hashlib
import os
import shutil
import subprocess
import tempfile
import numpy as np

_FUNCTION = None


def _load():
    global _FUNCTION
    if _FUNCTION is not None:
        return _FUNCTION
    source = Path(__file__).with_name("native_kernel.cpp")
    compiler = os.environ.get("CXX") or shutil.which("c++") or shutil.which("g++")
    if not compiler:
        raise RuntimeError("Native acceleration requires a C++17 compiler; use the Python backend instead.")
    tag = hashlib.sha256(source.read_bytes()).hexdigest()[:16]
    cache = Path(tempfile.gettempdir()) / "ek1d_native"
    cache.mkdir(parents=True, exist_ok=True)
    shared = cache / f"kernel_{tag}.so"
    if not shared.exists():
        temporary = cache / f"kernel_{tag}_{os.getpid()}.so"
        command = [compiler, "-std=c++17", "-O3", "-ffp-contract=off", "-shared", "-fPIC", str(source), "-o", str(temporary)]
        compiled = subprocess.run(command, text=True, capture_output=True)
        if compiled.returncode:
            raise RuntimeError(f"Native compiler failed:\n{compiled.stderr}")
        temporary.replace(shared)
    library = ctypes.CDLL(str(shared))
    function = library.advance
    pointer = ctypes.POINTER(ctypes.c_double)
    function.argtypes = [ctypes.c_int, pointer, pointer, pointer, ctypes.c_double,
                         pointer, pointer, ctypes.c_char_p, ctypes.c_int]
    function.restype = ctypes.c_int
    function._library_reference = library
    _FUNCTION = function
    return function


def advance(c, cs, cp, p, duration, initial_inventories=None):
    """Return (c_new, cs_new, cp_new, diagnostics) after ``duration`` seconds.

    ``initial_inventories`` is the global initial Cu/Na/NO3 areal inventory
    (mol/m2) used for cumulative conservation checks. If omitted, the incoming
    chunk inventories are the reference. Inputs are never modified.
    """
    aqueous = np.array(c, dtype=np.float64, order="C", copy=True)
    sorbed = np.array(cs, dtype=np.float64, order="C", copy=True)
    precipitated = np.array(cp, dtype=np.float64, order="C", copy=True)
    if aqueous.shape != (5, p.N) or sorbed.shape != (p.N,) or precipitated.shape != (p.N,):
        raise ValueError("Inconsistent native state array dimensions.")
    if not np.array_equal(p.z, np.array([1., -1., 2., 1., -1.])):
        raise ValueError("Native species charges are fixed to H/OH/Cu/Na/NO3.")
    diffusivities = np.asarray(p.D0, dtype=float)
    if diffusivities.shape != (5,) or not np.isfinite(diffusivities).all() or np.any(diffusivities <= 0):
        raise ValueError("Expected five finite positive molecular diffusivities.")
    if initial_inventories is None:
        initial_inventories = p.theta*p.L/p.N*np.array([
            np.sum(aqueous[2]+sorbed+precipitated), aqueous[3].sum(), aqueous[4].sum()])
    initial = np.asarray(initial_inventories, dtype=float)
    if initial.shape != (3,):
        raise ValueError("Expected three initial inventories: Cu, Na, NO3.")
    params = np.array([p.L, p.theta, p.tau, p.V, p.Kw, p.Ksp, p.pH50,
                       p.nedge, p.dt, p.cfl, *diffusivities, *initial, getattr(p,"spatial_order",1),
                       float(getattr(p,"algebraic_chemistry",True)),getattr(p,"time_order",1)], dtype=np.float64)
    metrics = np.zeros(25, dtype=np.float64)
    error = ctypes.create_string_buffer(1024)
    pointer = ctypes.POINTER(ctypes.c_double)
    ptr = lambda a: a.ctypes.data_as(pointer)
    result = _load()(p.N, ptr(aqueous), ptr(sorbed), ptr(precipitated), float(duration),
                     ptr(params), ptr(metrics), error, len(error))
    if result:
        raise RuntimeError(error.value.decode("utf-8", errors="replace"))
    diagnostics = dict(
        accepted_steps=int(metrics[0]), step_min_s=float(metrics[1]),
        step_max_s=float(metrics[2]), minimum_aqueous_mol_m3=float(metrics[3]),
        max_charge_mol_m3=float(metrics[4]), max_normalized_charge=float(metrics[5]),
        max_current_residual_A_m2=float(metrics[6]), max_voltage_residual_V=float(metrics[7]),
        max_relative_inventory_errors=dict(zip(("Cu", "Na", "NO3"), metrics[8:11].tolist())),
        max_chemistry_residuals=dict(zip(("water", "sorption", "precipitation"), metrics[11:14].tolist())),
        faradaic_moles_each_electrode_per_m2=float(metrics[14]), advanced_duration_s=float(metrics[15]),
        max_local_chemistry_charge_residual_mol_m3=float(metrics[16]),
        max_local_chemistry_Cu_residual_mol_m3=float(metrics[17]),
        max_chemistry_iterations=int(metrics[18]), algebraic_chemistry_calls=int(metrics[19]),
        generic_chemistry_calls=int(metrics[20]), rejected_steps=int(metrics[21]),
        # Retain the historical alias for callers outside the tutorial runner.
        rejected_RK2_attempts=int(metrics[21]),
        current_stage_min_A_m2=float(metrics[22]),
        current_stage_max_A_m2=float(metrics[23]),
        integrated_signed_current_C_m2=float(metrics[24]),
        evaluated_state_diagnostics_include_rejected_trials=True,
        time_order=int(getattr(p,"time_order",1)),backend="C++17 exact-equation acceleration")
    return aqueous, sorbed, precipitated, diagnostics


def equilibrate_native(total, charge, p, initial_H=None, algebraic=True):
    """Isolated native equilibrium, returning H/OH/Cu/Cs/Cp and diagnostics.

    Set ``algebraic=False`` to retain the original bounded log-H solver.
    The optional n=2 algebraic evaluation uses the same invariant equations;
    failures fall back to that generic log-H solver.
    """
    total, charge = np.broadcast_arrays(np.asarray(total,dtype=float),np.asarray(charge,dtype=float))
    shape=total.shape
    tt=np.ascontiguousarray(total.ravel());aa=np.ascontiguousarray(charge.ravel())
    hh=np.full(tt.size,np.nan) if initial_H is None else np.ascontiguousarray(np.broadcast_to(initial_H,shape).ravel(),dtype=float)
    pp=np.zeros(20,dtype=float);pp[4:8]=[p.Kw,p.Ksp,p.pH50,p.nedge];pp[19]=float(algebraic)
    values=np.empty((5,tt.size),dtype=float);metrics=np.zeros(21,dtype=float)
    error=ctypes.create_string_buffer(1024);pointer=ctypes.POINTER(ctypes.c_double)
    function=_load()._library_reference.equilibrate_native
    function.argtypes=[ctypes.c_int,pointer,pointer,pointer,pointer,pointer,pointer,ctypes.c_char_p,ctypes.c_int]
    function.restype=ctypes.c_int
    ptr=lambda a:a.ctypes.data_as(pointer)
    result=function(tt.size,ptr(tt),ptr(aa),ptr(hh),ptr(pp),ptr(values),ptr(metrics),error,len(error))
    if result:
        raise RuntimeError(error.value.decode("utf-8",errors="replace"))
    diagnostics=dict(max_charge_residual_mol_m3=float(metrics[16]),max_Cu_residual_mol_m3=float(metrics[17]),
                     max_iterations=int(metrics[18]),algebraic_calls=int(metrics[19]),generic_calls=int(metrics[20]))
    return values.reshape((5,)+shape),diagnostics
