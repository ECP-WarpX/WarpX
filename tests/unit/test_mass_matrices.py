# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Remi Lehe
# License: BSD-3-Clause-LBNL

"""Unit tests for the mass matrices of the implicit solvers.

The mass matrices ``S`` are the linear response of the deposited current
density to the electric field, ``dJ = S dE``, which the implicit solvers use
in place of pushing and depositing the particles at every linear iteration.
They are deposited in ``Source/Particles/Deposition/MassMatricesDeposition.H``
and applied in ``ImplicitSolver::ApplyMassMatrices``. The tests below check
them against the thing they stand in for: the current that WarpX itself
deposits after pushing the particles in that same electric field.
"""

import numpy as np
import pytest

import pywarpx
from pywarpx import picmi

constants = picmi.constants

# RZ deposits with an inverse volume scaling and rotates the mass matrices
# into cylindrical components; make_sim does not build that geometry yet.
pytestmark = pytest.mark.skipif(
    pywarpx.libwarpx.geometry_dim not in ("1d", "2d", "3d"),
    reason="the mass matrices comparison is set up for Cartesian geometries",
)


def _theta_implicit_with_mass_matrices():
    """Evolve scheme that allocates and can deposit the full mass matrices."""
    newton = picmi.NewtonNonlinearSolver(
        linear_solver=picmi.GMRESLinearSolver(),
        use_mass_matrices_jacobian=True,
    )
    return picmi.ThetaImplicitEMEvolveScheme(nonlinear_solver=newton, theta=0.5)


def _alloc_like(sim, name, template, n_grow_extra=0):
    """Register a zeroed vector field with the layout of vector field ``template``.

    Same box arrays, staggering and guard cells (plus ``n_grow_extra``) for
    each of the three components, so that the new field can stand in for
    ``template`` wherever the C++ side expects that staggering.
    """
    fields = sim.fields
    for direction in ("x", "y", "z"):
        mf = fields.get(template, direction, 0)
        fields.alloc_init(
            name,
            direction,
            0,
            mf.box_array(),
            mf.dm(),
            mf.n_comp,
            mf.n_grow_vect + n_grow_extra,
            0.0,
            False,
            False,
        )


def _fill_periodic_random(mf, n_cell, rng, amplitude):
    """Fill ``mf``, guard cells included, with a random periodic field.

    One random value is drawn per cell of the domain and every point of the
    MultiFab, whether valid or guard, nodal or cell-centered, reads the value
    of the cell it wraps to. A nodal point on the upper boundary therefore
    equals its periodic image on the lower boundary, and the guard cells hold
    what a periodic fill would put there, without any communication.
    """
    table = rng.uniform(-amplitude, amplitude, size=n_cell)
    # imesh puts cell-centered points at half-integers; floor gives the cell
    wrapped = [
        np.floor(mf.imesh(idir, include_ghosts=True)).astype(int) % n_cell[idir]
        for idir in range(len(n_cell))
    ]
    mf[()] = table[np.ix_(*wrapped)]


@pytest.mark.parametrize("particle_shape", ["linear", "quadratic", "cubic"])
def test_mass_matrices_match_push_and_deposit(
    make_sim, uniform_particles, rtol, particle_shape
):
    """``S dE`` must equal the current deposited after a push in ``dE``.

    The particles start at rest, which makes this an exact identity rather
    than a linearization: with ``u = 0`` the Lorentz factor is exactly one on
    both sides, and the current does not depend on the positions to first
    order, so the response of the deposited current to the electric field is
    the Boris rotation matrix that the mass matrices are built from.

    Both sides use WarpX's own kernels: ``push_p`` gathers ``dE`` and ``B``
    at the particles and rotates their momentum, ``deposit_current`` puts the
    result back on the grid. What the test pins down is that the mass
    matrices, deposited and applied as banded stencils, reproduce that
    gather-push-deposit chain: kernel, shape factors, staggering of every
    (J, E) pair, guard cell exchange between boxes and periodic wrapping.

    The mass matrices give the response of the time-centered current
    ``(u^n + u^{n+1}) / 2``, whereas the push from rest leaves ``u^{n+1}`` on
    the particles, hence the factor 1/2 on the reference.
    """
    # The mass matrices use plain shape factors for the field gather (see
    # MassMatricesDeposition.H), so the pusher must too. WarpX already turns
    # the Galerkin (energy-conserving) gather off for the direct deposition
    # with an electromagnetic solver; this pins it down independently of that.
    pywarpx.interpolation.galerkin_scheme = 0

    n_axes = {"1d": 1, "2d": 2, "3d": 3}[pywarpx.libwarpx.geometry_dim]
    # 8 cells and 4 cells per box: two boxes per axis, so that contributions
    # crossing a box boundary and the periodic boundary are both exercised
    sim = make_sim(
        n_cell=[8] * n_axes,
        max_grid_size=4,
        particle_shape=particle_shape,
        current_deposition_algo="direct",
        evolve_scheme=_theta_implicit_with_mass_matrices(),
    )
    warpx = sim.extension.warpx
    fields = sim.fields
    dt = warpx.getdt(0)
    n_cell = list(warpx.Geom(0).domain.size)

    # the particles start at rest, see the docstring
    electrons, _ = uniform_particles(sim)

    # A uniform magnetic field with all three components, strong enough that
    # the normalized gyration ``b = q dt B / (2 m)`` is of order one: the
    # off-diagonal blocks of the mass matrices and their ``1 / (1 + b^2)``
    # denominator only matter for a magnetized push.
    b_unit = 2.0 * constants.m_e / (constants.q_e * dt)
    for direction, b in zip(("x", "y", "z"), (0.6, -0.8, 1.1)):
        fields.get("Bfield_fp", direction, 0).set_val(b * b_unit)

    # Deposit the mass matrices the way the Darwin solver does: the deposit
    # leaves the contributions of particles near a box edge in the guard
    # cells, so sum those into the valid cells before mirroring the symmetric
    # half of the diagonal blocks. The valid cells are then complete.
    solver = warpx.implicit_solver()
    warpx.deposit_mass_matrices()
    warpx.sync_mass_matrices()
    solver.finish_mass_matrices()

    # ApplyMassMatrices reads ``dE`` as far as the band of each (J, E) pair
    # reaches, and silently truncates the band at the guard cells of ``dE``.
    # Along a direction where J is nodal and E is cell-centered (or the other
    # way round) the band is one component wider, so it reaches nox + 1 cells:
    # one more than the guard cells of J and, for the quadratic shape, also one
    # more than ``Efield_fp`` has. Give ``dE`` enough guard cells for the full
    # band, so that this test checks the mass matrices themselves and not the
    # guard cells of ``Efield_fp``.
    n_grow_j = fields.get("current_fp", "x", 0).n_grow_vect
    n_grow_e = fields.get("Efield_fp", "x", 0).n_grow_vect
    n_grow_extra = max(
        0, max(n_grow_j[idir] + 1 - n_grow_e[idir] for idir in range(n_axes))
    )
    _alloc_like(sim, "dE", "Efield_fp", n_grow_extra=n_grow_extra)
    _alloc_like(sim, "dJ", "current_fp")
    _alloc_like(sim, "J_ref", "current_fp")

    # The amplitude keeps the push non-relativistic: q dE dt / m is a fraction
    # of a meter per second, so gamma is one to far better than the tolerance.
    rng = np.random.default_rng(seed=42)
    for direction in ("x", "y", "z"):
        _fill_periodic_random(fields.get("dE", direction, 0), n_cell, rng, 1.0)

    # mass matrices: dJ = S dE
    solver.apply_mass_matrices("dJ", "dE", zero_out_first=True)

    # reference: push from rest in (dE, B), then deposit the current
    sim.particles.push_p(
        0,
        dt,
        *(fields.get("dE", direction, 0) for direction in ("x", "y", "z")),
        *(fields.get("Bfield_fp", direction, 0) for direction in ("x", "y", "z")),
    )
    electrons.deposit_current("J_ref", 0, dt, 0.0)

    for direction in ("x", "y", "z"):
        dj_mass_matrices = fields.get("dJ", direction, 0)[...]
        dj_reference = 0.5 * fields.get("J_ref", direction, 0)[...]

        # a magnetized push from rest in a random field drives all components
        scale = np.max(np.abs(dj_reference))
        assert scale > 0.0

        # atol scaled to the current that is actually flowing, so that the
        # entries that happen to be small are held to the same accuracy as
        # the large ones and the assertion is not vacuous
        error = np.max(np.abs(dj_mass_matrices - dj_reference)) / scale
        assert np.allclose(
            dj_mass_matrices, dj_reference, rtol=rtol, atol=rtol * scale
        ), f"J{direction}: max |dJ_mm - dJ_ref| / max |dJ_ref| = {error}"
