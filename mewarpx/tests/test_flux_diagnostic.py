import collections
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas
import pytest

from mewarpx.utils_store import util as mwxutil


def test_timeseries():
    name = "timeseries"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx.utils_store import testing_util
    from mewarpx.diags_store import timeseries

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed, Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(47239315)

    dt = 1e-12
    data1 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    data2 = np.array([10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
    data3 = np.array([20, 21, 22, 23, 24, 25, 26, 27, 28, 29])
    datadict = {"a": data1, "b": data2, "c": data3}

    ts = timeseries.Timeseries(0, 10, dt, datadict)

    assert np.all(np.isclose(ts.get_timeseries_by_key("a", False), data1))
    assert np.mean(data2) == ts.get_averagevalue_by_key("b")

    resampled_ts = ts.resample(dt * 2)

    assert np.all(np.isclose(
        np.array([[0.0e-12, 20],[2.0e-12, 22], [4.0e-12, 24], [6.0e-12, 25], [8.0e-12, 27]]),
        resampled_ts.get_timeseries_by_key("c", True)
    ))


def test_injector_flux_diagnostic():
    name = "injectorfluxDiagnostic"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx.utils_store import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    import dill

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed, Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(47239475)

    TOTAL_TIME = 1e-10 # s
    DIAG_INTERVAL = 1e-10 # s
    DT = 1e-12 # s

    P_INERT = 1 # torr
    T_INERT = 300 # K

    D_CA = 5e-4 # m
    VOLTAGE = 25 # V
    CATHODE_TEMP = 1100 + 273.15 # K
    CATHODE_PHI = 2.1 # work function in eV
    NX = 8
    NZ = 128

    DIRECT_SOLVER = True

    max_steps = int(TOTAL_TIME / DT)
    diag_steps = int(DIAG_INTERVAL / DT)

    run = diode_setup.DiodeRun_V1(
        dim=dim,
        CATHODE_TEMP=CATHODE_TEMP,
        CATHODE_PHI=CATHODE_PHI,
        V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        P_INERT=P_INERT,
        T_INERT=T_INERT,
        NPPC=50,
        NX=NX,
        NZ=NZ,
        DIRECT_SOLVER=DIRECT_SOLVER,
        PERIOD=D_CA * NX / NZ,
        DT=DT,
        TOTAL_TIMESTEPS=max_steps,
        DIAG_STEPS=diag_steps,
        DIAG_INTERVAL=DIAG_INTERVAL,
        CHECK_CHARGE_CONSERVATION=False
    )

    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=True,
        init_scraper=False,
        init_runinfo=True,
        init_fluxdiag=True,
        init_warpx=True
    )

    mwxrun.simulation.step(max_steps)

    # check that current in CSV is correct
    filename = "diags/fluxes/thermionic_injector_electrons_injected.csv"

    assert os.path.isfile(filename), "Could not find output CSV."

    df = pandas.read_csv(filename)

    q = df["q"]

    cathode_area = run.emitter.area
    J = mwxutil.J_RD(run.injector.T, run.injector.WF, run.injector.A)
    current = J * cathode_area * -1

    assert np.allclose(q / DT, current, rtol=0.01, atol=0)

    # check that pickle file current is correct
    with open("diags/fluxes/fluxdata_{:010d}.dpkl".format(100), "rb") as pfile:
        flux_diags = dill.load(pfile)

    # leave off first element (time step 0) where there is 0 current
    J_pickle = (flux_diags['fullhist_dict'][('inject', 'cathode', 0)]
                .get_timeseries_by_key("J", False)[1:])

    # convert to cm^2
    J /= -1.0e4

    assert np.allclose(J_pickle, J, rtol=0.01, atol=0)

    # check that plotfile exists
    assert os.path.isfile("diags/fluxes/flux_plots_{:010d}.png".format(100))