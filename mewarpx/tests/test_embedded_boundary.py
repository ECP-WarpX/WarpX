"""Test the mewarpx wrapper for embedded boundaries. This test is the same
as the test in warpx/Examples/Tests/ElectrostaticSphereEB/inputs_3d but in 2d.
"""
import pytest
import numpy as np
import os

from mewarpx import util as mwxutil

def test_embedded_cylinder():
    name = "Embedded_cylinder_solve"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import testing_util, assemblies
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(92160881)

    # Specific numbers match older run for consistency
    DT = 1e-6
    DIAG_STEPS = 1
    D_CA = 1  # m
    NX = 64
    NZ = 64
    run = diode_setup.DiodeRun_V1(
        dim=dim,
        #V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        NX=NX,
        NZ=NZ,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=1,
        DIAG_STEPS=DIAG_STEPS,
        FIELD_DIAG_DATA_LIST=['phi'],
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_electron_plasma=False,
        init_solver=False,
        init_scraper=False,
        init_injectors=False,
        init_field_diag=True,
        init_simcontrol=True,
        init_simulation=False,
        init_warpx=False
    )

    # Install the embedded boundary
    cylinder = assemblies.Cylinder(
        center_x=0.0, center_z=0.5, radius=0.3,
        V=-2.0, T=300, WF=4.7, name="Cylinder"
    )

    # Initialize solver
    run.init_solver()

    # Initialize the simulation
    run.init_simulation()
    run.init_warpx()

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Check phi results against reference data                            #
    #######################################################################

    phi = mwxrun.get_gathered_phi_grid()[0]
    ref_phi = np.load(os.path.join(
        testing_util.test_dir, 'embedded_boundary', 'embedded_cylinder_phi.npy'
    ))
    assert np.allclose(phi, ref_phi)
