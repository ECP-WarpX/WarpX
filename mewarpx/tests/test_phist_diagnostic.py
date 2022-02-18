"""Test the particle histogram diagnostic used to track positions where
particles are scraped.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from mewarpx import assemblies, emission
from mewarpx.diags_store.particle_histogram_diagnostic import ZPlanePHistDiag
from mewarpx.mwxrun import mwxrun
from mewarpx.setups_store import diode_setup
from mewarpx.utils_store import testing_util


def test_particle_hist_diag():
    name = "particle_hist_diag"
    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(2216083)

    # Specific numbers match older run for consistency
    D_CA = 20e-6  # m
    run = diode_setup.DiodeRun_V1(
        GEOM_STR='XZ',
        V_ANODE_CATHODE=10,
        D_CA=D_CA,
        NX=64,
        NZ=64,
        DT=1e-12,
        TOTAL_TIMESTEPS=500,
        DIAG_STEPS=500,
        FIELD_DIAG_DATA_LIST=['phi'],
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_simcontrol=True,
        init_simulation=False
    )

    # Install the embedded boundary
    cylinder = assemblies.InfCylinderY(
        center_x=0.0, center_z=5e-6, radius=1.5e-6,
        V=-2.0, T=300, WF=4.7, name="Cylinder"
    )
    run.electrons.save_particles_at_eb = 1

    # Initialize solver
    run.init_solver()
    run.init_conductors()

    # Add particle histogram diagnostic for the anode
    phist = ZPlanePHistDiag(run.DIAG_STEPS, run.anode)

    # Initialize the simulation
    run.init_runinfo()
    run.init_fluxdiag()
    run.init_simulation()
    run.init_warpx()

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Check histogram results against reference                           #
    #######################################################################

    ref_hist = np.load(os.path.join(
        testing_util.test_dir, 'histograms', 'Anode_electrons.npy'
    ))
    hist = np.load('diags/histograms/Anode_electrons_0000000500.npy')
    assert np.allclose(hist, ref_hist)
