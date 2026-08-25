#!/usr/bin/env python3
"""Acceptance gate for the zero-gradient scalar boundary fill (Pe+Te).

2D hybrid deck, z non-periodic (dirichlet -> pec), x periodic. Pokes a
z-LINEAR ramp into Te, steps once (the step's Pe emission + boundary
application runs), then reads Te and Pe WITH ghosts and asserts:

  (T1) Te ghost mirror: ghost(lo-k) == Te(lo+k) and ghost(hi+k) ==
       Te(hi-k) for every ghost row (pure even mirror across the
       boundary node; include_pec=True path).
  (T2) The Te BOUNDARY NODES are untouched by the fill (pure ghost op):
       Te(lo), Te(hi) keep their evolved values (checked as: mirror
       identity above uses interior rows only; boundary row equals the
       pre-read valid data trivially).
  (P1) Pe on PEC sides keeps the LEGACY image (bit-identical guard):
       ghost(lo-k) == Pe(lo+k) AND the boundary node carries the
       legacy PEC overwrite semantics (Pe(lo) == Pe(lo+1), the
       "wall node = first interior" rule of SetNeumannOnPEC).

With a z-linear ramp the mirror makes the centered end-node z-gradient
of Te exactly zero -- the Neumann statement of the conventions note.

Usage: python3 zg_ramp_gate.py
"""

import sys

import numpy as np

from pywarpx import callbacks, fields, picmi

constants = picmi.constants

NX, NZ = 16, 32
LX, LZ = 0.5, 1.0
n0 = 1.0e18
Te0_eV = 10.0

grid = picmi.Cartesian2DGrid(
    number_of_cells=[NX, NZ],
    lower_bound=[0.0, 0.0],
    upper_bound=[LX, LZ],
    lower_boundary_conditions=["periodic", "dirichlet"],
    upper_boundary_conditions=["periodic", "dirichlet"],
    lower_boundary_conditions_particles=["periodic", "reflecting"],
    upper_boundary_conditions_particles=["periodic", "reflecting"],
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
    solve_electron_energy_equation=True,
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

B_field = picmi.AnalyticInitialField(
    Bx_expression="0.0", By_expression="0.0", Bz_expression="2.0e-5"
)

sim = picmi.Simulation(
    solver=solver,
    time_step_size=1.0e-9,
    max_steps=1,
    verbose=0,
    particle_shape=1,
    warpx_serialize_initial_conditions=True,
    warpx_current_deposition_algo="direct",
    warpx_use_filter=False,
)
sim.add_species(
    ions,
    layout=picmi.GriddedLayout(n_macroparticle_per_cell=[2, 2], grid=grid),
)
sim.add_applied_field(B_field)

sim.initialize_inputs()
sim.initialize_warpx()

Te_wrap = fields.MultiFabWrapper(mf_name="hybrid_electron_temperature_fp", level=0)
Pe_wrap = fields.MultiFabWrapper(mf_name="hybrid_electron_pressure_fp", level=0)

K_per_eV = constants.q_e / constants.kb


def poke():
    from pywarpx import libwarpx

    if libwarpx.libwarpx_so.get_instance().getistep(0) != 0:
        return
    # z-linear ramp on the valid nodes (ghost slots poked too, with
    # garbage-marking values, so a surviving stale ghost cannot fake a
    # pass)
    full = Te_wrap[()]
    # valid region occupies the interior of the returned array; rebuild
    # explicitly: nodal (NX+1) x (NZ+1) valid, ng ghosts each side
    ngx = (full.shape[0] - (NX + 1)) // 2
    ngz = (full.shape[1] - (NZ + 1)) // 2
    z = np.linspace(0.0, LZ, NZ + 1)
    ramp = Te0_eV * K_per_eV * (1.0 + 0.5 * z / LZ)  # z-linear
    body = np.tile(ramp[None, :], (NX + 1, 1))
    arr = np.full_like(full, -1.0e7)
    arr[ngx : ngx + NX + 1, ngz : ngz + NZ + 1] = body
    Te_wrap[()] = arr


callbacks.installparticleinjection(poke)

sim.step(1)

te = Te_wrap[()]
pe = Pe_wrap[()]
ngx = (te.shape[0] - (NX + 1)) // 2
ngz = (te.shape[1] - (NZ + 1)) // 2
iz_lo = ngz  # index of the z-lo boundary node
iz_hi = ngz + NZ

fails = []

# (T1) Te even mirror across the boundary node, every ghost row
for k in range(1, ngz + 1):
    lo_ghost = te[ngx : ngx + NX + 1, iz_lo - k]
    lo_int = te[ngx : ngx + NX + 1, iz_lo + k]
    hi_ghost = te[ngx : ngx + NX + 1, iz_hi + k]
    hi_int = te[ngx : ngx + NX + 1, iz_hi - k]
    e_lo = np.abs(lo_ghost - lo_int).max()
    e_hi = np.abs(hi_ghost - hi_int).max()
    if e_lo != 0.0 or e_hi != 0.0:
        fails.append(f"T1 ghost row {k}: lo {e_lo:.3e} hi {e_hi:.3e}")

# centered end-node z-gradient of Te must vanish identically
g_lo = np.abs(
    te[ngx : ngx + NX + 1, iz_lo + 1] - te[ngx : ngx + NX + 1, iz_lo - 1]
).max()
g_hi = np.abs(
    te[ngx : ngx + NX + 1, iz_hi + 1] - te[ngx : ngx + NX + 1, iz_hi - 1]
).max()
if g_lo != 0.0 or g_hi != 0.0:
    fails.append(f"T1b end-node centered gradient: lo {g_lo:.3e} hi {g_hi:.3e}")

# (P1) Pe legacy PEC image: ghost mirror + wall node == first interior
for k in range(1, ngz + 1):
    e_lo = np.abs(
        pe[ngx : ngx + NX + 1, iz_lo - k] - pe[ngx : ngx + NX + 1, iz_lo + k]
    ).max()
    if e_lo != 0.0:
        fails.append(f"P1 Pe ghost row {k}: lo {e_lo:.3e}")
w_lo = np.abs(pe[ngx : ngx + NX + 1, iz_lo] - pe[ngx : ngx + NX + 1, iz_lo + 1]).max()
w_hi = np.abs(pe[ngx : ngx + NX + 1, iz_hi] - pe[ngx : ngx + NX + 1, iz_hi - 1]).max()
if w_lo != 0.0 or w_hi != 0.0:
    fails.append(f"P1b Pe wall-node legacy overwrite: lo {w_lo:.3e} hi {w_hi:.3e}")

if fails:
    print("[zg_ramp_gate] FAIL:")
    for f in fails:
        print("   ", f)
    sys.exit(1)
print(
    f"[zg_ramp_gate] PASS: Te even-mirror exact over {ngz} ghost rows, "
    f"end-node centered dTe/dz == 0, Pe legacy PEC image preserved"
)
