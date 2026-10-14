#!/usr/bin/env python3
"""Exercise material evolution and current sources on static refined grids."""

import argparse
import sys
from pathlib import Path

import numpy as np

from pywarpx import amr, libwarpx, picmi

parser = argparse.ArgumentParser()
parser.add_argument("--dims", choices=("2", "3", "RZ"), required=True)
parser.add_argument("--source", choices=("port", "antenna", "material"), required=True)
parser.add_argument(
    "--placement", choices=("inside", "outside", "crossing"), default="inside"
)
parser.add_argument(
    "--method", choices=("backwardeuler", "laxwendroff"), default="backwardeuler"
)
parser.add_argument("--restart")
args, remaining = parser.parse_known_args()
sys.argv = sys.argv[:1] + remaining

ndim = 3 if args.dims == "3" else 2
is_rz = args.dims == "RZ"
lower = [0.0, -0.5] if is_rz else [-0.5] * ndim
upper = [0.5, 0.5] if is_rz else [0.5] * ndim
refined_lower = [0.0625, -0.375] if is_rz else [-0.375] * ndim
refined_upper = [0.4375, 0.375] if is_rz else [0.375] * ndim
if args.placement == "outside":
    refined_upper[0] = 0.125 if is_rz else -0.25
elif args.placement == "crossing" or args.source == "antenna":
    refined_upper[0] = 0.25 if is_rz else 0.0

grid_type = {
    "2": picmi.Cartesian2DGrid,
    "3": picmi.Cartesian3DGrid,
    "RZ": picmi.CylindricalGrid,
}[args.dims]
grid = grid_type(
    number_of_cells=[32] * ndim,
    lower_bound=lower,
    upper_bound=upper,
    lower_boundary_conditions=["none", "periodic"] if is_rz else ["periodic"] * ndim,
    upper_boundary_conditions=["dirichlet", "periodic"]
    if is_rz
    else ["periodic"] * ndim,
    lower_boundary_conditions_particles=["none", "periodic"]
    if is_rz
    else ["periodic"] * ndim,
    upper_boundary_conditions_particles=["reflecting", "periodic"]
    if is_rz
    else ["periodic"] * ndim,
    refined_regions=[[1, refined_lower, refined_upper]],
    warpx_max_grid_size=16,
    warpx_blocking_factor=4,
)
amr.n_error_buf = 0

dt = 1.0e-12
simulation = picmi.Simulation(
    solver=picmi.ElectromagneticSolver(
        grid=grid, method="Yee", divE_cleaning=args.source == "antenna"
    ),
    max_steps=2,
    time_step_size=dt,
    particle_shape="linear",
    warpx_current_deposition_algo="esirkepov",
    warpx_use_filter=0,
    warpx_serialize_initial_conditions=1,
    warpx_amr_restart=args.restart,
)
if not is_rz:
    # Parser fields distinguish fine/coarse coordinates and exercise non-vacuum mu.
    for name, expression in (
        ("epsilon", f"{picmi.constants.ep0}*(2+x)"),
        ("sigma", f"{picmi.constants.ep0}*1.e10*(2+x)"),
        ("mu", f"{picmi.constants.mu0}*(2+0.5*z)"),
    ):
        simulation.add_macroscopic_property(
            picmi.MacroscopicProperty(
                name=name,
                implicit_function=expression,
                method=args.method if name == "sigma" else None,
            )
        )

waveform = str(Path(__file__).resolve().parent / "current_profile.txt")
if args.source == "port":
    transverse_lo = 0.1875 if is_rz else -0.125
    transverse_hi = 0.3125 if is_rz else 0.125
    simulation.add_current_controlled_port(
        picmi.CurrentControlledPort(
            direction=2,
            terminal_0_lower_bound=[transverse_lo, -0.125, -0.125],
            terminal_0_upper_bound=[transverse_hi, 0.125, -0.125],
            terminal_1_lower_bound=[transverse_lo, -0.125, 0.125],
            terminal_1_upper_bound=[transverse_hi, 0.125, 0.125],
            file=waveform,
        )
    )
elif args.source == "antenna":
    simulation.add_prescribed_current_injection(
        picmi.PrescribedCurrentInjection(
            file=waveform,
            drives=[
                picmi.PrescribedCurrentDrive(
                    lower_bound=[0.125 if is_rz else -0.25, -0.25, -0.25],
                    upper_bound=[0.375 if is_rz else 0.25, 0.25, 0.25],
                    direction=2,
                    area=np.pi * (0.375**2 - 0.125**2) if is_rz else 0.25,
                )
            ],
        )
    )
else:
    assert not is_rz
    simulation.add_applied_field(
        picmi.AnalyticInitialField(
            Ex_expression="0", Ey_expression="1.0", Ez_expression="0"
        )
    )

diagnostic_fields = ["E", "B", "J", "divB", "rho"]
if not is_rz:
    diagnostic_fields += ["sigma", "epsilon", "mu"]
simulation.add_diagnostic(
    picmi.FieldDiagnostic(
        name="diag1",
        grid=grid,
        period=1,
        data_list=diagnostic_fields,
        write_dir="diags",
        warpx_file_prefix="diag1",
    )
)
simulation.add_diagnostic(
    picmi.Checkpoint(name="chk", period=1, write_dir="diags", warpx_file_prefix="chk")
)
simulation.step(1 if args.restart else 2)

if args.source == "antenna":
    # Discrete continuity is a stronger oracle than normalization alone. The
    # nodal endpoint charge must match the divergence of the staggered Jz.
    for level in (0, 1):
        rho = np.asarray(simulation.fields.get("rho_fp", level=level)[...])
        jz = np.asarray(simulation.fields.get("current_fp", dir="z", level=level)[...])
        rho = np.squeeze(rho)
        jz = np.squeeze(jz)
        dz = (upper[-1] - lower[-1]) / (32 * 2**level)
        interior = (slice(1, -1),) * ndim
        rate = (rho[..., 1] - rho[..., 0])[interior] / dt
        divergence = np.diff(jz, axis=ndim - 1) / dz
        divergence = divergence[(slice(1, -1),) * (ndim - 1) + (slice(None),)]
        scale = max(np.max(np.abs(rate)), np.max(np.abs(divergence)))
        assert scale > 0
        error = np.max(np.abs(rate + divergence)) / scale
        print(f"level {level} discrete continuity error: {error}")
        assert error < 2.0e-10

if args.source == "port":
    status = np.asarray(libwarpx.warpx.get_current_controlled_port_status())
    print(f"terminal status: {status}")
    np.testing.assert_allclose(status, [2000.0] * 3, rtol=2.0e-12, atol=0)

    # Measure native B independently of the status API on every represented patch.
    for level, patch in ((0, "fp"), (1, "fp"), (1, "cp"), (1, "aux")):
        if level == 1 and args.placement == "outside":
            continue
        dx = np.array(upper) - lower
        dx /= 32 * (2 if level == 1 and patch != "cp" else 1)

        reference = simulation.fields.get(
            f"Bfield_{patch}",
            dir="x" if args.dims == "3" else "y" if not is_rz else "theta",
            level=level,
        )
        origin = np.array(reference.box_array().minimal_box().small_end)

        def field(component):
            return np.squeeze(
                np.asarray(
                    simulation.fields.get(
                        f"Bfield_{patch}", dir=component, level=level
                    )[...]
                )
            )

        def index(coordinate, axis, offset=0):
            return round((coordinate - lower[axis]) / dx[axis] - offset) - origin[axis]

        planes = [index(z, ndim - 1) for z in (-0.125, 0, 0.125)]
        if args.dims == "3":
            bx, by = field("x"), field("y")
            # Yee curl contours surround the vector-potential top hat.
            lo, hi = index(-0.125, 0), index(0.125, 0)
            z = lower[2] + (np.array(planes) + origin[2]) * dx[2]
            mu = picmi.constants.mu0 * (2 + 0.5 * (z + 0.5 * dx[2]))
            currents = (
                (
                    np.sum(bx[lo : hi + 1, lo - 1][:, planes], axis=0)
                    + np.sum(by[hi, lo : hi + 1][:, planes], axis=0)
                    - np.sum(bx[lo : hi + 1, hi][:, planes], axis=0)
                    - np.sum(by[lo - 1, lo : hi + 1][:, planes], axis=0)
                )
                * dx[0]
                / mu
            )
        elif is_rz:
            bt = field("theta")
            lo, hi = index(transverse_lo, 0), index(transverse_hi, 0)
            currents = (
                2
                * np.pi
                * (transverse_hi * bt[hi, planes] - transverse_lo * bt[lo, planes])
                / picmi.constants.mu0
            )
        else:
            by = field("y")
            lo, hi = index(transverse_lo, 0), index(transverse_hi, 0)
            z = lower[1] + (np.array(planes) + origin[1]) * dx[1]
            mu = picmi.constants.mu0 * (2 + 0.5 * (z + 0.5 * dx[1]))
            currents = (by[hi, planes] - by[lo, planes]) / mu
        print(f"independent level {level} {patch} currents: {currents}")
        np.testing.assert_allclose(currents, 2000.0, rtol=2.0e-10, atol=0)
