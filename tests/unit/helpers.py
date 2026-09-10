# Copyright 2026 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Axel Huebl
# License: BSD-3-Clause-LBNL

"""Setup helpers for the WarpX unit tests.

The test modules import from here, which relies on pytest's default ``prepend``
import mode putting ``tests/unit`` on ``sys.path``. That is what
``add_warpx_pytest()`` runs; ``--import-mode=importlib`` would not resolve it.
"""

import numpy as np

import pywarpx
from pywarpx import picmi

# PICMI grid per Cartesian geometry, and how many axes it takes. RZ is
# deliberately absent: CylindricalGrid needs n_azimuthal_modes and a different
# set of boundary conditions, so it wants its own helper once a test needs it.
GRID_CLASS = {
    "1d": picmi.Cartesian1DGrid,
    "2d": picmi.Cartesian2DGrid,
    "3d": picmi.Cartesian3DGrid,
}
N_AXES = {"1d": 1, "2d": 2, "3d": 3}

# spreads particles over the domain without aligning them with the cells
GOLDEN_RATIO = 0.6180339887498949


def make_sim(
    n_cell=None,
    lower_bound=None,
    upper_bound=None,
    max_grid_size=8,
    particle_shape="quadratic",
    current_deposition_algo=None,
    dt=None,
):
    """Build and initialize a minimal simulation of this process' dimensionality.

    Only one simulation can live at a time. ``warpx_lifecycle`` in
    ``conftest.py`` tears it down after each test, which is what makes
    parametrizing over e.g. the current deposition algorithm possible.
    """
    dims = pywarpx.libwarpx.geometry_dim
    if dims not in GRID_CLASS:
        raise NotImplementedError(
            f"make_sim does not build a {dims} geometry yet. Tests that rely on "
            "it should skip on a non-Cartesian geometry."
        )

    # A second simulation in one test would reach amrex_init() with AMReX
    # already initialized, which aborts the whole pytest process. Fail here
    # instead, with a traceback that points at the test.
    assert not pywarpx.libwarpx.initialized, (
        "a simulation is already running; only one simulation can live at "
        "a time, and warpx_lifecycle finalizes it after the test"
    )

    n_axes = N_AXES[dims]
    n_cell = [16] * n_axes if n_cell is None else list(n_cell)
    lower_bound = [-1.0e-3] * n_axes if lower_bound is None else list(lower_bound)
    upper_bound = [1.0e-3] * n_axes if upper_bound is None else list(upper_bound)

    # dt=None lets WarpX pick the CFL-limited time step, which keeps the
    # per-step particle displacement a sizeable fraction of a cell
    grid = GRID_CLASS[dims](
        number_of_cells=n_cell,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
        lower_boundary_conditions=["periodic"] * n_axes,
        upper_boundary_conditions=["periodic"] * n_axes,
        lower_boundary_conditions_particles=["periodic"] * n_axes,
        upper_boundary_conditions_particles=["periodic"] * n_axes,
        # more than one box, so that the guard cell exchange is exercised
        warpx_max_grid_size=max_grid_size,
    )
    solver = picmi.ElectromagneticSolver(grid=grid, method="Yee", cfl=0.9)
    electrons = picmi.Species(particle_type="electron", name="electrons")

    sim = picmi.Simulation(
        solver=solver,
        time_step_size=dt,
        max_steps=1,
        verbose=0,
        particle_shape=particle_shape,
        warpx_current_deposition_algo=current_deposition_algo,
    )
    # no particles are injected: the tests add them explicitly
    sim.add_species(electrons, layout=None)

    sim.initialize_inputs()

    # AMReX runtime parameters, mirroring the ones ImpactX and pyAMReX use
    # in their pytest suites
    pywarpx.warpx.get_bucket("tiny_profiler").enabled = 0
    #   throw exceptions instead of writing Backtrace files, so a debugger
    #   can be attached
    pywarpx.amrex.throw_exception = 1
    pywarpx.amrex.signal_handling = 0
    #   abort GPU runs on out-of-memory instead of swapping to host RAM
    pywarpx.amrex.abort_on_out_of_gpu_memory = 1
    #   allocate GPU memory on demand instead of pre-allocating 3/4th, so
    #   that tests can share a GPU
    pywarpx.amrex.the_arena_init_size = 0

    sim.initialize_warpx()

    return sim


def uniform_particles(sim, n_per_dim=4, weight=1.0e6, ux=0.0, uy=0.0, uz=0.0):
    """Add a lattice of macro particles to the ``electrons`` species.

    ``n_per_dim`` positions per grid axis, so this yields ``n_per_dim``
    particles in 1D and ``n_per_dim**3`` in 3D. Returns the particle container
    together with the arrays it was built from, which the tests need to form
    the expected values.
    """
    geom = sim.extension.warpx.Geom(0)
    lo = np.array(geom.ProbLo())
    hi = np.array(geom.ProbHi())

    # Fractional positions in (0, 1). A regular lattice would be commensurate
    # with the grid and put every particle at the same offset inside its cell,
    # so use a golden ratio sequence instead: it spreads the particles evenly
    # over the domain while giving each of them a different sub-cell offset,
    # which keeps the shape factors away from their symmetric special case.
    frac = ((np.arange(n_per_dim) + 1) * GOLDEN_RATIO) % 1.0
    axes = np.meshgrid(*(frac,) * len(lo), indexing="ij")
    coords = [(lo[i] + axis * (hi[i] - lo[i])).ravel() for i, axis in enumerate(axes)]

    n_part = coords[0].size
    zeros = np.zeros(n_part)
    if len(lo) == 1:
        x, y, z = zeros, zeros, coords[0]
    elif len(lo) == 2:
        x, y, z = coords[0], zeros, coords[1]
    else:
        x, y, z = coords

    w = np.full(n_part, weight)
    uxp = np.full(n_part, ux)
    uyp = np.full(n_part, uy)
    uzp = np.full(n_part, uz)

    electrons = sim.particles.get("electrons")
    electrons.add_particles(
        x=x, y=y, z=z, ux=uxp, uy=uyp, uz=uzp, w=w, unique_particles=False
    )

    return electrons, dict(x=x, y=y, z=z, ux=uxp, uy=uyp, uz=uzp, w=w)
