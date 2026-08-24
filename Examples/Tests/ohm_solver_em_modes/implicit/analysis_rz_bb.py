#!/usr/bin/env python3
#
# --- Analysis for the pc_block_banded arm of the RZ theta-implicit hybrid
# --- EM-modes deck. The preconditioner only changes the Krylov path, so the
# --- delivered physics must agree with the unpreconditioned parent test
# --- (which this test declares as a CTest dependency) to within the
# --- accumulated linear-solve tolerance, while the iteration counts must
# --- show the preconditioner actually working. The operator-correctness
# --- gates themselves (extraction self-check, LU/Apply round-trips, FD-JVP
# --- comparison) run inside the solver at every rebuild (precond.bb_verify)
# --- and abort the run on failure, so a completed run has already passed
# --- them.

import numpy as np
from openpmd_viewer import OpenPMDTimeSeries

BASE = "../test_rz_ohm_solver_em_modes_implicit_picmi"

# Iteration-count regression bounds, calibrated 2026-08-18 on a 2-rank CPU
# run of this deck: PC arm Newton max 6 / GMRES mean 15.6 max 19 per step,
# parent Newton max 6 / GMRES mean 22.6 max 27 per step. At dt x1 the
# whistler is not stiff, so the margin between the arms is modest and these
# ceilings are a coarse tripwire only; the dt x8 stiff arm
# (test_rz_ohm_solver_em_modes_implicit_bb_stiff_picmi) carries the sharp
# preconditioner-quality bound (11 vs 44 GMRES/step).
NEWTON_MAX_CEIL = 8
GMRES_MEAN_CEIL = 20.0
GMRES_MAX_CEIL = 25

# Physics agreement ceilings (relative L2 over all components at the last
# common output step). A right preconditioner only changes the Krylov path,
# so with both arms converged to the same Newton tolerance the fields agree
# far below these values (measured 2026-08-18: E 4.5e-6, B 1.3e-10 at step
# 30, the B denominator carrying the static B0); the ceilings guard the
# silent-solution-corruption class (seam write-back, ghost propagation),
# which shows up at O(1e-2) or worse, with headroom for platform and
# noise-realization differences.
E_REL_CEIL = 1e-3
B_REL_CEIL = 1e-7


def load_newton_diag(path):
    data = np.loadtxt(path)
    data = np.atleast_2d(data)
    return {
        "newton_iters": data[:, 2],
        "gmres_iters": data[:, 6],
    }


pc = load_newton_diag("diags/newton_diag.txt")
base = load_newton_diag(f"{BASE}/diags/newton_diag.txt")

newton_max = pc["newton_iters"].max()
gmres_mean = pc["gmres_iters"].mean()
gmres_max = pc["gmres_iters"].max()

print("pc_block_banded arm:")
print(
    f"  Newton iters/step: max {newton_max:.0f}, mean {pc['newton_iters'].mean():.2f}"
)
print(f"  GMRES iters/step:  max {gmres_max:.0f}, mean {gmres_mean:.2f}")
print("unpreconditioned parent:")
print(
    f"  Newton iters/step: max {base['newton_iters'].max():.0f}, mean {base['newton_iters'].mean():.2f}"
)
print(
    f"  GMRES iters/step:  max {base['gmres_iters'].max():.0f}, mean {base['gmres_iters'].mean():.2f}"
)

assert newton_max <= NEWTON_MAX_CEIL, (
    f"Newton iterations regressed: max {newton_max} > ceiling {NEWTON_MAX_CEIL}"
)
assert gmres_mean <= GMRES_MEAN_CEIL, (
    f"mean GMRES/step regressed: {gmres_mean:.2f} > ceiling {GMRES_MEAN_CEIL}"
)
assert gmres_max <= GMRES_MAX_CEIL, (
    f"max GMRES/step regressed: {gmres_max} > ceiling {GMRES_MAX_CEIL}"
)
# strict win on the per-step mean (the arms may run different step counts,
# so totals are not comparable)
assert gmres_mean < base["gmres_iters"].mean(), (
    "preconditioned arm did not beat the unpreconditioned parent on mean "
    f"GMRES/step ({gmres_mean:.2f} vs {base['gmres_iters'].mean():.2f})"
)

# --- physics agreement at the last common output step
ts_pc = OpenPMDTimeSeries("diags/field_diags", check_all_files=True)
ts_base = OpenPMDTimeSeries(f"{BASE}/diags/field_diags", check_all_files=True)
common = sorted(set(ts_pc.iterations) & set(ts_base.iterations))
assert len(common) > 1, "arms share no output iterations beyond the initial one"
it = common[-1]

for name, ceil in (("E", E_REL_CEIL), ("B", B_REL_CEIL)):
    num = 0.0
    den = 0.0
    for coord in ("x", "y", "z"):
        f_pc, _ = ts_pc.get_field(name, coord, iteration=it)
        f_base, _ = ts_base.get_field(name, coord, iteration=it)
        num += np.sum((f_pc - f_base) ** 2)
        den += np.sum(f_base**2)
    rel = np.sqrt(num / den)
    print(f"{name}: relative L2 difference vs parent at step {it} = {rel:.3e}")
    assert rel <= ceil, (
        f"{name} disagrees with the unpreconditioned parent: {rel:.3e} > {ceil:.1e}"
    )

print("pc_block_banded CI gates passed")
