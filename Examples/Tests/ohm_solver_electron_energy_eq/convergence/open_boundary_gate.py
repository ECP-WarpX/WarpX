#!/usr/bin/env python3
"""Acceptance gate for the open (outflow) domain-boundary treatment.

2D hybrid deck, x periodic, z-lo = PMC (the symmetry side), z-hi = Open
(the outflow side, zero-gradient continuation). Pokes a z-LINEAR ramp
into Bx at step 0 (ghost slots garbage-marked so a surviving stale
ghost cannot fake a pass), steps once, then reads all six E/B
components WITH ghosts and asserts on the z-hi (Open) side:

  (V1) every component's z-hi ghosts are the exact staggering-aware
       even mirror of the valid data: nodal-in-z components mirror
       across the boundary NODE, cell-centered-in-z components across
       the boundary FACE (parity inferred per component from the
       array extent, so the gate is staggering-agnostic);
  (V2) for nodal-in-z components the centered end-node z-derivative
       vanishes identically (the Neumann/continuation statement);
  (V3) valid data is untouched by the fill (pure-ghost property,
       implied by V1 using post-step valid data as the reference).

The z-lo (PMC) side is deliberately not asserted here: PMC applies its
own per-component image (upstream-tested); this gate only verifies
that the Open fill acts on the Open side.

Usage: python3 open_boundary_gate.py
"""

import sys

import numpy as np

import pywarpx
from pywarpx import callbacks, fields, picmi

constants = picmi.constants

NX, NZ = 16, 32
LX, LZ = 0.5, 1.0
n0 = 1.0e18
Te0_eV = 10.0
B0 = 0.05  # T

grid = picmi.Cartesian2DGrid(
    number_of_cells=[NX, NZ],
    lower_bound=[0.0, 0.0],
    upper_bound=[LX, LZ],
    lower_boundary_conditions=["periodic", "dirichlet"],
    upper_boundary_conditions=["periodic", "dirichlet"],
    lower_boundary_conditions_particles=["periodic", "reflecting"],
    upper_boundary_conditions_particles=["periodic", "absorbing"],
    warpx_max_grid_size=32,
)

solver = picmi.HybridPICSolver(
    grid=grid,
    gamma=5.0 / 3.0,
    Te=Te0_eV,
    n0=n0,
    n_floor=0.01 * n0,
    plasma_resistivity=1.0e-6,
    substeps=4,
)

ions = picmi.Species(
    particle_type="H",
    name="ions",
    charge_state=1,
    mass=1.0e6 * constants.m_p,
    initial_distribution=picmi.AnalyticDistribution(
        density_expression=f"{n0}",
        momentum_expressions=["0.0", "0.0", "0.0"],
    ),
)

sim = picmi.Simulation(
    solver=solver,
    time_step_size=1.0e-9,
    max_steps=1,
    verbose=0,
    warpx_serialize_initial_conditions=True,
    warpx_current_deposition_algo="direct",
)
sim.add_species(
    ions, layout=picmi.PseudoRandomLayout(grid=grid, n_macroparticles_per_cell=4)
)
sim.add_applied_field(
    picmi.AnalyticInitialField(
        Bx_expression="0.0", By_expression="0.0", Bz_expression=f"{B0}"
    )
)

sim.initialize_inputs()
# the boundary bucket is written by initialize_inputs from the picmi
# grid, so the field-BC override must come AFTER it (and before
# initialize_warpx): z-lo = pmc (symmetry side), z-hi = open (outflow
# continuation side)
pywarpx.boundary.field_lo = ["periodic", "pmc"]
pywarpx.boundary.field_hi = ["periodic", "open"]
sim.initialize_warpx()

# read the _fp state directly (the boundary fill acts on it); the
# wrapper's [()] view includes the ghost cells
Bx_wrap = fields.MultiFabWrapper(mf_name="Bfield_fp", idir=0, level=0)
By_wrap = fields.MultiFabWrapper(mf_name="Bfield_fp", idir=1, level=0)
Bz_wrap = fields.MultiFabWrapper(mf_name="Bfield_fp", idir=2, level=0)
Ex_wrap = fields.MultiFabWrapper(mf_name="Efield_fp", idir=0, level=0)
Ey_wrap = fields.MultiFabWrapper(mf_name="Efield_fp", idir=1, level=0)
Ez_wrap = fields.MultiFabWrapper(mf_name="Efield_fp", idir=2, level=0)


def poke():
    from pywarpx import libwarpx

    if libwarpx.libwarpx_so.get_instance().getistep(0) != 0:
        return
    # z-linear ramp on Bx (valid region only computed from the extent),
    # ghost slots garbage-marked
    full = Bx_wrap[()]
    nz_valid = NZ + 1 if (full.shape[1] - NZ) % 2 == 1 else NZ
    nx_valid = NX + 1 if (full.shape[0] - NX) % 2 == 1 else NX
    ngx = (full.shape[0] - nx_valid) // 2
    ngz = (full.shape[1] - nz_valid) // 2
    z = np.linspace(0.0, LZ, nz_valid)
    body = np.tile(0.2 * B0 * (1.0 + z / LZ)[None, :], (nx_valid, 1))
    arr = np.full_like(full, -1.0e7)
    arr[ngx : ngx + nx_valid, ngz : ngz + nz_valid] = body
    Bx_wrap[()] = arr


callbacks.installparticleinjection(poke)

sim.step(1)

fails = []


def check(name, wrap):
    a = wrap[()]
    nz_valid = NZ + 1 if (a.shape[1] - NZ) % 2 == 1 else NZ
    nx_valid = NX + 1 if (a.shape[0] - NX) % 2 == 1 else NX
    ngx = (a.shape[0] - nx_valid) // 2
    ngz = (a.shape[1] - nz_valid) // 2
    # exclude the DUPLICATED x-node of the periodic seam for nodal-x
    # fields: that column is owned by the nodal-sync/comms machinery
    # (owner-image override at communication precision), not by the
    # local ghost fill, so bitwise exactness is not guaranteed there
    nodal_x = nx_valid == NX + 1
    xs = slice(ngx, ngx + nx_valid - (1 if nodal_x else 0))
    nodal_z = nz_valid == NZ + 1
    if nodal_z:
        ihb = ngz + NZ  # boundary node
        for k in range(1, ngz + 1):
            err = np.abs(a[xs, ihb + k] - a[xs, ihb - k]).max()
            if err != 0.0:
                fails.append(f"V1 {name} (nodal-z) ghost row {k}: {err:.3e}")
        gd = np.abs(a[xs, ihb + 1] - a[xs, ihb - 1]).max()
        if gd != 0.0:
            fails.append(f"V2 {name} end-node centered d/dz: {gd:.3e}")
    else:
        ihc = ngz + NZ - 1  # last valid cell; boundary face at its high edge
        for k in range(1, ngz + 1):
            err = np.abs(a[xs, ihc + k] - a[xs, ihc - k + 1]).max()
            if err != 0.0:
                fails.append(f"V1 {name} (cc-z) ghost row {k}: {err:.3e}")


for nm, w in (
    ("Bx", Bx_wrap),
    ("By", By_wrap),
    ("Bz", Bz_wrap),
    ("Ex", Ex_wrap),
    ("Ey", Ey_wrap),
    ("Ez", Ez_wrap),
):
    check(nm, w)

if fails:
    print("[open_boundary_gate] FAIL:")
    for f in fails:
        print("   ", f)
    sys.exit(1)
print(
    "[open_boundary_gate] PASS: all six E/B components even-mirror exact "
    "on the z-hi Open side (staggering-aware), end-node centered "
    "derivatives zero on nodal-z components"
)
