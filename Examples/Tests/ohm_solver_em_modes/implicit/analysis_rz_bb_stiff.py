#!/usr/bin/env python3
#
# --- Analysis for the whistler-stiff (dt x8) pc_block_banded arm of the RZ
# --- theta-implicit hybrid EM-modes deck. At this time step the grid-scale
# --- whistler CFL makes the unpreconditioned linear solve expensive
# --- (~44-60 GMRES iterations/step measured on this deck), while the
# --- block-banded direct preconditioner holds the count flat, so the
# --- iteration ceiling below is a sharp regression tripwire for the
# --- preconditioner itself. Convergence is enforced in the run
# --- (newton.require_convergence) and the operator-correctness gates
# --- (extraction self-check, LU/Apply round-trips, FD-JVP comparison)
# --- abort the run on failure, so this script only checks the counts.

import numpy as np

data = np.atleast_2d(np.loadtxt("diags/newton_diag.txt"))
newton_iters = data[:, 2]
gmres_iters = data[:, 6]

newton_max = newton_iters.max()
gmres_mean = gmres_iters.mean()
gmres_max = gmres_iters.max()

print("pc_block_banded stiff arm (dt x8):")
print(f"  Newton iters/step: max {newton_max:.0f}, mean {newton_iters.mean():.2f}")
print(f"  GMRES iters/step:  max {gmres_max:.0f}, mean {gmres_mean:.2f}")

# Calibrated 2026-08-18 (2-rank CPU, OMP_NUM_THREADS=1): PC-on measured a
# flat 11 GMRES / 3 Newton per step; the PC-off control on the same deck
# (line search on) measured a flat 44 GMRES / 3 Newton per step. The mean
# ceiling sits at ~2x the preconditioned count and well below half the
# unpreconditioned count.
NEWTON_MAX_CEIL = 8
GMRES_MEAN_CEIL = 20.0
GMRES_MAX_CEIL = 30

assert newton_max <= NEWTON_MAX_CEIL, (
    f"Newton iterations regressed: max {newton_max} > ceiling {NEWTON_MAX_CEIL}"
)
assert gmres_mean <= GMRES_MEAN_CEIL, (
    f"mean GMRES/step regressed: {gmres_mean:.2f} > ceiling {GMRES_MEAN_CEIL}"
)
assert gmres_max <= GMRES_MAX_CEIL, (
    f"max GMRES/step regressed: {gmres_max} > ceiling {GMRES_MAX_CEIL}"
)

print("pc_block_banded stiff-arm CI gates passed")
