"""Test utility to create voltage pulses."""
import os

import numpy as np
from pywarpx import picmi

from mewarpx.utils_store import util as mwxutil


def test_utils_pulsing_expr():

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from mewarpx.utils_store import pulsing

    pulse_expr = pulsing.linear_pulse_function(
        V_off=-1.0, V_on=7.5, pulse_period=20e-9, pulse_length=10e-9,
        t_rise=5e-9, t_fall=5e-9, plot=False
    )
    print(pulse_expr)
    assert pulse_expr == "((fmod((t-5e-09),2e-08)>=0 and fmod((t-5e-09),2e-08)<5e-09)*(-1.0+1700000000.0*fmod((t-5e-09),2e-08))+(fmod((t-5e-09),2e-08)>=5e-09 and fmod((t-5e-09),2e-08)<1.5000000000000002e-08)*7.5+(fmod((t-5e-09),2e-08)>=1.5000000000000002e-08 and fmod((t-5e-09),2e-08)<2e-08)*(7.5-1700000000.0*(fmod((t-5e-09),2e-08)-1.5000000000000002e-08))+(if((fmod((t-5e-09),2e-08)>=0 and fmod((t-5e-09),2e-08)<5e-09) or (fmod((t-5e-09),2e-08)>=5e-09 and fmod((t-5e-09),2e-08)<1.5000000000000002e-08) or (fmod((t-5e-09),2e-08)>=1.5000000000000002e-08 and fmod((t-5e-09),2e-08)<2e-08), 0, 1))*-1.0)"

def test_utils_pulsing_sim():
    name = "pulsing_utility_test"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import assemblies, mespecies
    from mewarpx.mwxrun import mwxrun
    from mewarpx.poisson_pseudo_1d import PoissonSolverPseudo1D
    from mewarpx.utils_store import pulsing, testing_util

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(9216024)

    nx = 8
    nz = 128
    D_CA = 1e-4 # m

    xmin = 0.0
    zmin = 0.0

    xmax = D_CA / nz * nx
    zmax = D_CA

    max_steps = 50

    DT = 1e-10

    #####################################
    # grid, solver and timesteps
    #####################################

    mwxrun.init_grid([xmin, zmin], [xmax, zmax], [nx, nz])

    solver = PoissonSolverPseudo1D(grid=mwxrun.grid)

    mwxrun.simulation.solver = solver
    mwxrun.init_timestep(DT=DT)
    mwxrun.simulation.max_steps = max_steps

    pulse_expr = pulsing.linear_pulse_function(
        V_off=-1.0, V_on=7.5, pulse_period=20e-9,
        pulse_length=2e-9, t_rise=1e-9, t_fall=1e-9,
        wait_time=0.5e-9, plot=True
    )
    anode = assemblies.Anode(D_CA, pulse_expr, 100, 3.0)

    electrons = mespecies.Species(
        particle_type='electron',
        name='electrons',
    )

    # mwxrun.simulation.write_input_file()
    # exit()

    #################################
    # simulation setup
    ################################

    mwxrun.init_run()

    # check that plot was saved successfully
    assert os.path.isfile('diags/circuit/pulse_schematic.png')

    ################################
    # Simulation run
    ################################

    times = np.arange(max_steps)*DT
    Vs = np.zeros(max_steps)

    # Run the main loop
    for ii in range(max_steps):
        Vs[ii] = anode.getvoltage()
        mwxrun.simulation.step(1)

    print(repr(Vs))
    ref_Vs = np.array(
        [-1.  , -1.  , -1.  , -1.  , -1.  , -1.  , -1.  , -0.15,  0.7 ,
        1.55,  2.4 ,  3.25,  4.1 ,  4.95,  5.8 ,  6.65,  7.5 ,  7.5 ,
        7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,
        7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,  7.5 ,
        7.5 ,  6.65,  5.8 ,  4.95,  4.1 ,  3.25,  2.4 ,  1.55,  0.7 ,
       -0.15, -1.  , -1.  , -1.  , -1.  ]
    )
    assert np.allclose(Vs, ref_Vs, atol=0)
