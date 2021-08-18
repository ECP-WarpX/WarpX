import os
import numpy as np
import pandas
import re

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
    J_pickle = (flux_diags['fullhist_dict'][('inject', 'cathode', 'electrons')]
                .get_timeseries_by_key("J", False)[1:])

    # convert to cm^2
    J /= -1.0e4

    assert np.allclose(J_pickle, J, rtol=0.01, atol=0)

    # check that plotfile exists
    assert os.path.isfile("diags/fluxes/flux_plots_{:010d}.png".format(100))


def test_flux_diag_accuracy(capsys):
    name = "FluxDiagRun"
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from mewarpx.utils_store import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun
    from mewarpx.diags_store import flux_diagnostic

    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. For initial
    # statistics, instead use a plain random seed (commented out here).
    np.random.seed(18672033)

    # these values come from the default values in the mewarp diode setup
    TOTAL_TIME = 3e-10 # s
    DIAG_INTERVAL = 1e-10 # s
    DT = 1e-12 # s not from mewarp diode setup
    CATHODE_TEMP = 1550 # K
    ANODE_TEMP = 773 # K

    P_INERT = 4 # torr
    T_INERT = .5 * (CATHODE_TEMP + ANODE_TEMP) # K
    NPPC = 2
    D_CA = 1.0e-4

    CATHODE_PHI = 2.4 # work function in eV
    NX = 8 # not from mewarp diode setup
    NZ = 128 # not from mewarp diode setup

    DIRECT_SOLVER = True

    max_steps = int(TOTAL_TIME / DT)
    diag_steps = int(DIAG_INTERVAL / DT)

    # Run with relatively large bias to have stronger signals.
    run = diode_setup.DiodeRun_V1(
        dim=2,
        CATHODE_TEMP=CATHODE_TEMP,
        CATHODE_PHI=CATHODE_PHI,
        V_ANODE_CATHODE=4.,
        D_CA=D_CA,
        P_INERT=P_INERT,
        T_INERT=T_INERT,
        NPPC=NPPC,
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

    run.setup_run(
        init_runinfo=True,
        init_fluxdiag=True,
        init_warpx=True
    )

    mwxrun.simulation.step(max_steps)

    captured = capsys.readouterr()
    outtext = captured[0]
    print(outtext)

    # Check expected print output
    # Match a pattern to allow for numbers to change
    pattern = (
        r"""THIS DIAGNOSTIC PERIOD:
Cathode Electrons Current Emitted: -4\.5\d* A/cm\^2
Cathode Electrons Current Collected: 4\.[0123]\d* A/cm\^2
Cathode Electrons Current Net: -0\.[234]\d* A/cm\^2
Anode Electrons Current Collected: 0\.3\d* A/cm\^2
Total Current Emitted: -4\.5\d* A/cm\^2
Total Current Collected: 4\.[456]\d* A/cm\^2
Total Current Net: -?0\.[01]\d* A/cm\^2
Total Power Net: -[01]\.[019]\d* W/cm\^2 """
    ).replace("\n", " ")

    match = re.search(pattern, outtext.replace("\n", " "))

    assert match is not None

    print("Diagnostic output match:")
    print(match.group(0))

    #Check diagnostic files are present
    filelist = [
        'diags/fluxes/Anode_scraped.csv',
        'diags/fluxes/Cathode_scraped.csv',
        'diags/fluxes/flux_plots_0000000100.png',
        'diags/fluxes/flux_plots_0000000200.png',
        'diags/fluxes/flux_plots_0000000300.png',
        'diags/fluxes/fluxdata_0000000100.dpkl',
        'diags/fluxes/fluxdata_0000000200.dpkl',
        'diags/fluxes/fluxdata_0000000300.dpkl',
        'diags/fluxes/thermionic_injector_electrons_injected.csv',
    ]

    for filename in filelist:
        assert os.path.isfile(filename)

    # Check that Qs plus powers sum to about 0 for most recent period
    Q_emit = run.fluxdiag.ts_dict[
            ('inject', 'cathode', 'electrons')].get_averagevalue_by_key('dQ')
    Q_abs_cathode = run.fluxdiag.ts_dict[
            ('scrape', 'cathode', 'electrons')].get_averagevalue_by_key('dQ')
    Q_abs_anode = run.fluxdiag.ts_dict[
            ('scrape', 'anode', 'electrons')].get_averagevalue_by_key('dQ')
    P_anode = run.fluxdiag.ts_dict[
            ('scrape', 'anode', 'electrons')].get_averagevalue_by_key('P')
    conservation = Q_emit + Q_abs_cathode + Q_abs_anode - P_anode
    print(f"Q_emit: {Q_emit} W/cm^2")
    print(f"Q_abs_cathode: {Q_abs_cathode} W/cm^2")
    print(f"Q_abs_anode: {Q_abs_anode} W/cm^2")
    print(f"P_anode: {P_anode} W/cm^2")
    print(f"Conservation: {conservation} W/cm^2")
    assert abs(conservation) < 0.4

    # Gather results to check.  The index call here ensures there's one row to
    # assign into in the new DataFrame.

    df = pandas.DataFrame(index=list(range(1)))

    for fluxtype in ['J', 'P', 'dQ', 'n']:
        df['inject_cathode_' + fluxtype] = run.fluxdiag.fullhist_dict[
            ('inject', 'cathode', 'electrons')].get_averagevalue_by_key(fluxtype)
        df['scrape_cathode_' + fluxtype] = run.fluxdiag.fullhist_dict[
            ('scrape', 'cathode', 'electrons')].get_averagevalue_by_key(fluxtype)
        df['scrape_anode_' + fluxtype] = run.fluxdiag.fullhist_dict[
            ('scrape', 'anode', 'electrons')].get_averagevalue_by_key(fluxtype)

    assert testing_util.test_df_vs_ref(name, df)

    # Check that loaded results are identical
    fluxdiag_2 = flux_diagnostic.FluxDiagFromFile()

    for key in run.fluxdiag.fullhist_dict:
        for fluxtype in ['J', 'P', 'dQ', 'n']:
            assert np.isclose(
                run.fluxdiag.fullhist_dict[
                    key].get_averagevalue_by_key(fluxtype),
                fluxdiag_2.fullhist_dict[
                    key].get_averagevalue_by_key(fluxtype)
            )

            assert np.isclose(
                run.fluxdiag.ts_dict[
                    key].get_averagevalue_by_key(fluxtype),
                fluxdiag_2.ts_dict[
                    key].get_averagevalue_by_key(fluxtype)
            )

    reftext = fluxdiag_2.print_fluxes(fluxdiag_2.ts_dict)
    print(reftext)

    pattern = (
        r"""Cathode Electrons Current Emitted: -4\.5\d* A/cm\^2
Cathode Electrons Current Collected: 4\.[0123]\d* A/cm\^2
Cathode Electrons Current Net: -0\.[234]\d* A/cm\^2
Anode Electrons Current Collected: 0\.3\d* A/cm\^2
Total Current Emitted: -4\.5\d* A/cm\^2
Total Current Collected: 4\.[456]\d* A/cm\^2
Total Current Net: -?0\.[01]\d* A/cm\^2
Total Power Net: -[01]\.[019]\d* W/cm\^2 """
    ).replace("\n", " ")

    match = re.search(pattern, reftext.replace("\n", " "))

    assert match is not None

    print("Postprocess output match:")
    print(match.group(0))
