"""Test functionality in mewarpx.emission.py"""
import collections
import os
from matplotlib.pyplot import plot
import numpy as np
import pandas
import pytest

from mewarpx.utils_store import util as mwxutil


def test_thermionic_emission():
    name = "thermionicEmission"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx.utils_store import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    import mewarpx.utils_store.mwxconstants as constants

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed, Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(47239475)

    TOTAL_TIME = 1e-10 # s
    DIAG_INTERVAL = 5e-10 # s
    DT = 0.5e-12 # s

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
        DIAG_INTERVAL=DIAG_INTERVAL
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=True,
        init_scraper=False,
        init_warpx=True
    )

    mwxrun.simulation.step(max_steps)

    net_rho_grid = mwxrun.get_gathered_rho_grid(include_ghosts=False)[:, :, 0]
    ref_path = os.path.join(testing_util.test_dir,
                            "thermionic_emission",
                            "thermionic_emission.npy")

    # slice the reference data to just get the non ghost cells
    ref_rho_grid = np.load(ref_path)[2:11, 2:131]

    assert np.allclose(net_rho_grid, ref_rho_grid)


def test_thermionic_emission_with_Schottky():
    name = "thermionicEmissionSchottky"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx.utils_store import testing_util
    from mewarpx.setups_store import diode_setup
    from mewarpx.mwxrun import mwxrun

    import mewarpx.utils_store.mwxconstants as constants

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed, Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(47239475)

    TOTAL_TIME = 1e-10 # s
    DIAG_INTERVAL = 5e-10 # s
    DT = 0.5e-12 # s

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
        USE_SCHOTTKY=True,
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
        DIAG_INTERVAL=DIAG_INTERVAL
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=True,
        init_scraper=False,
        init_warpx=True
    )

    mwxrun.simulation.step(max_steps)

    net_rho_grid = mwxrun.get_gathered_rho_grid(include_ghosts=False)[:, :, 0]
    # np.save('thermionic_emission_Schottky.npy', net_rho_grid)
    ref_path = os.path.join(testing_util.test_dir,
                            "thermionic_emission",
                            "thermionic_emission_Schottky.npy")
    ref_rho_grid = np.load(ref_path)

    assert np.allclose(net_rho_grid, ref_rho_grid)


def test_circle_emitter():
    name = "circleEmitter"
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from pywarpx import picmi
    from mewarpx import assemblies, emission, mepicmi
    from mewarpx.utils_store import testing_util

    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(671237741)

    #####################################
    # embedded boundary, grid, and solver
    #####################################

    T_cylinder = 2173.15 # K
    cylinder = assemblies.Cylinder(
        center_x=0.2, center_z=0.4, radius=0.1, V=0, T=T_cylinder,
        WF=1.2, name='circle'
    )

    mwxrun.init_grid(0, 0.6, 0, 1, 20, 20)
    solver = picmi.ElectrostaticSolver(
        grid=mwxrun.grid, method='Multigrid', required_precision=1e-6
    )

    #################################
    # physics components
    ################################

    electrons = mepicmi.Species(particle_type='electron', name='electrons')

    #################################
    # simulation setup
    ################################

    mwxrun.simulation.solver = solver
    mwxrun.init_run()

    ######################################
    # Add ME emission
    #####################################

    emitter = emission.ArbitraryEmitter2D(
        conductor=cylinder, T=T_cylinder, res_fac=10
    )

    res_dict = emitter.get_newparticles(
        10000, 1, electrons.sq, electrons.sm, randomdt=False, velhalfstep=False
    )
    df = pandas.DataFrame(index=list(range(1)))

    # Compare main results, leave out E_total, since variation is too high
    for label in ['vx', 'vy', 'vz', 'x', 'y', 'z']:
        df[label + '_min'] = np.min(res_dict[label])
        df[label + '_max'] = np.max(res_dict[label])
        df[label + '_mean'] = np.mean(res_dict[label])
        df[label + '_std'] = np.std(res_dict[label])

    assert testing_util.test_df_vs_ref(
        testname="circle_emitter", df=df, margin=0.3
    )

    emitter.plot_contours()
    plotfile =  f"{cylinder.name}_contour_plot.png"

    assert os.path.exists(plotfile)


def test_rectangle_emitter():
    name = "rectangleEmitter"
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from pywarpx import picmi
    from mewarpx import assemblies, emission, mepicmi
    from mewarpx.utils_store import testing_util

    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(274737523)

    #####################################
    # embedded boundary, grid, and solver
    #####################################

    T_rectangle = 2173.15 # K
    # cylinder = assemblies.Cylinder(
    #     center_x=0.2, center_z=0.4, radius=0.1, V=0, T=T_box,
    #     WF=1.2, name='circle'
    # )

    rectangle = assemblies.Rectangle(
        center_x=0.2, center_z=0.4, length_x=0.2, length_z =0.4, V=0,
        T=T_rectangle, WF=1.2, name="rectangle")

    mwxrun.init_grid(0, 0.6, 0, 1, 20, 20)
    solver = picmi.ElectrostaticSolver(
        grid=mwxrun.grid, method='Multigrid', required_precision=1e-6
    )

    #################################
    # physics components
    ################################

    electrons = mepicmi.Species(particle_type='electron', name='electrons')

    #################################
    # simulation setup
    ################################

    mwxrun.simulation.solver = solver
    mwxrun.init_run()

    ######################################
    # Add ME emission
    #####################################

    emitter = emission.ArbitraryEmitter2D(
        conductor=rectangle, T=T_rectangle, res_fac=10
    )

    res_dict = emitter.get_newparticles(
        10000, 1, electrons.sq, electrons.sm, randomdt=False, velhalfstep=False
    )
    df = pandas.DataFrame(index=list(range(1)))

    # Compare main results, leave out E_total, since variation is too high
    for label in ['vx', 'vy', 'vz', 'x', 'y', 'z']:
        df[label + '_min'] = np.min(res_dict[label])
        df[label + '_max'] = np.max(res_dict[label])
        df[label + '_mean'] = np.mean(res_dict[label])
        df[label + '_std'] = np.std(res_dict[label])

    assert testing_util.test_df_vs_ref(
        testname="rectangle_emitter", df=df, margin=0.3
    )

    emitter.plot_contours()
    plotfile =  f"{rectangle.name}_contour_plot.png"

    assert os.path.exists(plotfile)


def test_plasma_injector():
    name = "plasmainjector"
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from pywarpx import picmi, _libwarpx
    from mewarpx import assemblies, emission
    from mewarpx.utils_store import testing_util
    from mewarpx.setups_store import diode_setup

    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(11874889)

    run = diode_setup.DiodeRun_V1(
        dim=2, PERIOD=0.8e-06*256, D_CA=0.8e-06*256, DT=1e-12,
        TOTAL_TIMESTEPS=3
    )

    run.setup_run(init_inert_gas=True, init_scraper=False, init_injectors=False)

    surfaceemitter = emission.ZPlaneEmitter(
        conductor=run.cathode, T=run.CATHODE_TEMP, ymin=0.0, ymax=0.0,
        transverse_fac=1.0,
    )
    volemitter = emission.ZSinDistributionVolumeEmitter(
        T=1550, zmin=0, zmax=run.D_CA,
    )
    volemitter_uniform = emission.UniformDistributionVolumeEmitter(
        T=1550, zmin=0, zmax=run.D_CA,
    )

    # Test under-specified error
    with pytest.raises(ValueError, match='Invalid plasma_density'):
        emission.PlasmaInjector(
            emitter=volemitter, species1=run.electrons,
            species2=run.ions, npart=40000
        )

    # Test inadequate pressure/temp error
    with pytest.raises(ValueError, match='Must specify positive'):
        emission.PlasmaInjector(
            emitter=volemitter, species1=run.electrons,
            species2=run.ions, npart=40000,
            ionization_frac=0.1, P_neutral=3,
        )
    # Test invalid ionization_frac with a surface emitter
    with pytest.raises(RuntimeError, match='Cannot use ionization_frac'):
        emission.PlasmaInjector(
            emitter=surfaceemitter, species1=run.electrons,
            species2=run.ions, npart=40000,
            ionization_frac=0.1, P_neutral=3, T_neutral=800,
        )

    # Test overspecified input
    with pytest.raises(
        ValueError, match='Specify ionization_frac or plasma_density'
    ):
        emission.PlasmaInjector(
            emitter=volemitter, species1=run.electrons,
            species2=run.ions, npart=5000,
            plasma_density=1e19,
            ionization_frac=0.1, P_neutral=3, T_neutral=800,
        )

    # Three real injectors
    # Volume emission with spec'd plasma_density
    emission.PlasmaInjector(
        emitter=volemitter, species1=run.electrons, species2=run.ions,
        npart=2000, plasma_density=1e19,
    )
    # Volume emission with spec'd ionization_frac
    emission.PlasmaInjector(
        emitter=volemitter_uniform, species1=run.electrons, species2=run.ions,
        npart=5000, ionization_frac=0.001, P_neutral=3, T_neutral=800,
        injectoffset=2, name='ionization_frac_based',
    )
    # Surface emission
    with pytest.warns(UserWarning,
                      match="Using a surface emitter with the PlasmaInjector "
                            "has not been tested for accuracy."):
        emission.PlasmaInjector(
            emitter=surfaceemitter, species1=run.electrons, species2=run.ions,
            npart=10000, plasma_density=1e16, injectoffset=1, injectfreq=1,
            name='surface_based',
        )

    run.init_warpx()
    mwxrun.simulation.step(3)

    npart_dict = mwxrun.get_npart_species_dict()

    # Gather results to check.  The index call here ensures there's one row to
    # assign into in the new DataFrame.
    df = pandas.DataFrame(index=list(range(1)))

    # grab particle data
    for species in [run.electrons, run.ions]:
        sname = species.name
        res_dict = collections.OrderedDict()
        res_dict['x'] = np.concatenate(_libwarpx.get_particle_x(sname), axis=0)
        res_dict['y'] = np.concatenate(_libwarpx.get_particle_y(sname), axis=0)
        res_dict['z'] = np.concatenate(_libwarpx.get_particle_z(sname), axis=0)
        res_dict['ux'] = np.concatenate(_libwarpx.get_particle_ux(sname), axis=0)
        res_dict['uy'] = np.concatenate(_libwarpx.get_particle_uy(sname), axis=0)
        res_dict['uz'] = np.concatenate(_libwarpx.get_particle_uz(sname), axis=0)
        res_dict['w'] = np.concatenate(_libwarpx.get_particle_weight(sname), axis=0)

        label_base = f'{species.name}_'
        df[label_base + 'npart'] = npart_dict[sname]

        # Compare main results
        for key, val in res_dict.items():
            df[label_base + key + '_min'] = np.min(val)
            df[label_base + key + '_max'] = np.max(val)
            df[label_base + key + '_mean'] = np.mean(val)
            df[label_base + key + '_std'] = np.std(val)

    assert testing_util.test_df_vs_ref(testname=name, df=df, margin=0.3)


def test_plasma_injector_fixedT2():
    name = "plasmainjector_fixedT2"
    mwxutil.init_libwarpx(ndim=2, rz=False)
    from pywarpx import _libwarpx
    from mewarpx import emission
    from mewarpx.utils_store import testing_util
    from mewarpx.setups_store import diode_setup

    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(11874889)

    run = diode_setup.DiodeRun_V1(
        dim=2, PERIOD=0.8e-06*256, D_CA=0.8e-06*256, DT=1e-12,
        TOTAL_TIMESTEPS=3
    )

    run.setup_run(init_inert_gas=True, init_scraper=False, init_injectors=False)

    volemitter = emission.ZSinDistributionVolumeEmitter(
        T=1550, zmin=0, zmax=run.D_CA,
    )

    # Volume emission with spec'd plasma_density
    emission.PlasmaInjector(
        emitter=volemitter, species1=run.electrons, species2=run.ions,
        npart=2000, T_2=500, plasma_density=1e19,
    )

    run.init_warpx()
    mwxrun.simulation.step(3)

    npart_dict = mwxrun.get_npart_species_dict()

    # Gather results to check.  The index call here ensures there's one row to
    # assign into in the new DataFrame.
    df = pandas.DataFrame(index=list(range(1)))

    # grab particle data
    for species in [run.electrons, run.ions]:
        sname = species.name
        res_dict = collections.OrderedDict()
        res_dict['x'] = np.concatenate(_libwarpx.get_particle_x(sname), axis=0)
        res_dict['y'] = np.concatenate(_libwarpx.get_particle_y(sname), axis=0)
        res_dict['z'] = np.concatenate(_libwarpx.get_particle_z(sname), axis=0)
        res_dict['ux'] = np.concatenate(_libwarpx.get_particle_ux(sname), axis=0)
        res_dict['uy'] = np.concatenate(_libwarpx.get_particle_uy(sname), axis=0)
        res_dict['uz'] = np.concatenate(_libwarpx.get_particle_uz(sname), axis=0)
        res_dict['w'] = np.concatenate(_libwarpx.get_particle_weight(sname), axis=0)

        label_base = f'{species.name}_'
        df[label_base + 'npart'] = npart_dict[sname]

        # Compare main results
        for key, val in res_dict.items():
            df[label_base + key + '_min'] = np.min(val)
            df[label_base + key + '_max'] = np.max(val)
            df[label_base + key + '_mean'] = np.mean(val)
            df[label_base + key + '_std'] = np.std(val)

    assert testing_util.test_df_vs_ref(testname=name, df=df, margin=0.3)
