"""Test functionality in mewarpx.emission.py"""
import collections
import os

import numpy as np
import pandas
import pytest
from pywarpx import picmi

from mewarpx import assemblies, diags, emission, mespecies
from mewarpx.mwxrun import mwxrun
from mewarpx.setups_store import diode_setup
from mewarpx.utils_store import testing_util


def test_thermionic_emission():
    name = "thermionicEmission"
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

    D_CA = 5e-4 # m
    VOLTAGE = 25 # V
    CATHODE_TEMP = 1100 + 273.15 # K
    CATHODE_PHI = 2.1 # work function in eV

    DIRECT_SOLVER = True

    max_steps = int(TOTAL_TIME / DT)
    diag_steps = int(DIAG_INTERVAL / DT)

    run = diode_setup.DiodeRun_V1(
        GEOM_STR='XZ',
        CATHODE_TEMP=CATHODE_TEMP,
        CATHODE_PHI=CATHODE_PHI,
        V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        NPPC=50,
        NX=8,
        NZ=128,
        DIRECT_SOLVER=DIRECT_SOLVER,
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
    # np.save('thermionic_emission.npy', net_rho_grid)
    ref_path = os.path.join(testing_util.test_dir,
                            "emission",
                            "thermionic_emission.npy")
    ref_rho_grid = np.load(ref_path)

    assert np.allclose(net_rho_grid, ref_rho_grid, rtol=1e-4)


def test_thermionic_emission_with_Schottky():
    name = "thermionicEmissionSchottky"
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

    D_CA = 5e-4 # m
    VOLTAGE = 25 # V
    CATHODE_TEMP = 1100 + 273.15 # K
    CATHODE_PHI = 2.1 # work function in eV

    DIRECT_SOLVER = True

    max_steps = int(TOTAL_TIME / DT)
    diag_steps = int(DIAG_INTERVAL / DT)

    run = diode_setup.DiodeRun_V1(
        GEOM_STR='XZ',
        CATHODE_TEMP=CATHODE_TEMP,
        CATHODE_PHI=CATHODE_PHI,
        USE_SCHOTTKY=True,
        V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        NPPC=50,
        NX=8,
        NZ=128,
        DIRECT_SOLVER=DIRECT_SOLVER,
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
                            "emission",
                            "thermionic_emission_Schottky.npy")
    ref_rho_grid = np.load(ref_path)

    assert np.allclose(net_rho_grid, ref_rho_grid, rtol=5e-4)


def test_vacuum_thermionic_diode():
    name = "vacuum_thermionic_diode"
    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed.
    np.random.seed(12178517)

    import sys
    sys.path.append(testing_util.example_dir)

    from thermionic_diode import PlanarVacuumTEC

    run = PlanarVacuumTEC(
        V_ANODE_CATHODE=12.5, TOTAL_TIMESTEPS=1000, SAVE=False, DIAG_STEPS=500,
        USE_SCHOTTKY=True,
    )
    run.NPPC = 2
    run.setup_run()
    run.run_sim()

    key = ('scrape', 'anode', 'electrons')
    J_diode = run.run.fluxdiag.ts_dict[key].get_averagevalue_by_key('J')
    assert np.isclose(J_diode, 2.3287, atol=1e-4)


def test_thermionic_emission_disc_rz():
    name = "thermionicEmissionDiscRZ"
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

    D_CA = 5e-4 # m
    VOLTAGE = 25 # V
    CATHODE_TEMP = 1100 + 273.15 # K
    CATHODE_PHI = 2.1 # work function in eV

    DIRECT_SOLVER = False

    max_steps = int(TOTAL_TIME / DT)
    diag_steps = int(DIAG_INTERVAL / DT)

    run = diode_setup.DiodeRun_V1(
        GEOM_STR='RZ',
        CATHODE_TEMP=CATHODE_TEMP,
        CATHODE_PHI=CATHODE_PHI,
        V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        NPPC=50,
        NX=8,
        NZ=128,
        DIRECT_SOLVER=DIRECT_SOLVER,
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
    # np.save("thermionic_emission_disc_rz.npy", net_rho_grid)
    ref_path = os.path.join(testing_util.test_dir,
                            "emission",
                            "thermionic_emission_disc_rz.npy")
    ref_rho_grid = np.load(ref_path)

    assert np.allclose(net_rho_grid, ref_rho_grid)


def test_circle_emitter():
    name = "circleEmitter"
    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(671237741)

    #####################################
    # grid, solver and embedded boundary
    #####################################

    mwxrun.init_grid(
        lower_bound=[0, 0], upper_bound=[0.6, 1], number_of_cells=[20, 20],
        min_tiles=16
    )
    solver = picmi.ElectrostaticSolver(
        grid=mwxrun.grid, method='Multigrid', required_precision=1e-6
    )

    T_cylinder = 2173.15 # K
    cylinder = assemblies.InfCylinderY(
        center_x=0.2, center_z=0.4, radius=0.1, V=0, T=T_cylinder,
        WF=1.2, name='circle'
    )

    #################################
    # physics components
    ################################

    electrons = mespecies.Species(particle_type='electron', name='electrons')

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


def test_xplane_emitter():
    name = "xplaneEmitter"
    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(277812124)

    #####################################
    # grid, solver and embedded boundary
    #####################################

    mwxrun.init_grid(
        lower_bound=[0, 0], upper_bound=[0.6, 1.0], number_of_cells=[20, 20],
        min_tiles=16
    )
    solver = picmi.ElectrostaticSolver(
        grid=mwxrun.grid, method='Multigrid', required_precision=1e-6
    )

    T_rectangle = 2173.15 # K
    rectangle = assemblies.Rectangle(
        center_x=0.15, center_z=0.4, length_x=0.2, length_z=0.4, V=0,
        T=T_rectangle, WF=1.2, name="rectangle")

    #################################
    # physics components
    ################################

    electrons = mespecies.Species(particle_type='electron', name='electrons')

    #################################
    # simulation setup
    ################################

    mwxrun.simulation.solver = solver
    mwxrun.init_run()

    ######################################
    # Add ME emission
    #####################################

    emitter = emission.XPlaneEmitter(conductor=rectangle)

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
        testname="xplane_emitter", df=df, margin=0.3
    )


def test_zcylinder_emitter():
    name = "zcylinderEmitter"
    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(27781552)

    #####################################
    # grid, solver and embedded boundary
    #####################################

    mwxrun.init_grid(
        lower_bound=[0, 0], upper_bound=[0.6, 1.0], number_of_cells=[20, 20],
        min_tiles=16, use_rz=True
    )
    solver = picmi.ElectrostaticSolver(
        grid=mwxrun.grid, method='Multigrid', required_precision=1e-6
    )

    T_cylinder = 2173.15 # K
    cylinder = assemblies.CylinderZ(
        r_outer=0.7, r_inner=0.6, V=0, T=T_cylinder, WF=1.2, name="cylinder")

    #################################
    # physics components
    ################################

    electrons = mespecies.Species(particle_type='electron', name='electrons')

    #################################
    # simulation setup
    ################################

    mwxrun.simulation.solver = solver
    mwxrun.init_run()

    ######################################
    # Add ME emission
    #####################################

    emitter = emission.ZCylinderEmitter(
        conductor=cylinder, zmin=0.1, zmax=0.9, rdir=-1
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
        testname="zcylinder_emitter", df=df, margin=0.3
    )


def test_rectangle_emitter():
    name = "rectangleEmitter"
    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(274737523)

    #####################################
    # grid, solver and embedded boundary
    #####################################

    mwxrun.init_grid(
        lower_bound=[0, 0], upper_bound=[0.6, 1], number_of_cells=[20, 20],
        min_tiles=16
    )
    solver = picmi.ElectrostaticSolver(
        grid=mwxrun.grid, method='Multigrid', required_precision=1e-6
    )

    T_rectangle = 2173.15 # K
    rectangle = assemblies.Rectangle(
        center_x=0.2, center_z=0.4, length_x=0.2, length_z=0.4, V=0,
        T=T_rectangle, WF=1.2, name="rectangle")

    #################################
    # physics components
    ################################

    electrons = mespecies.Species(particle_type='electron', name='electrons')

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
    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(11874889)

    run = diode_setup.DiodeRun_V1(
        GEOM_STR='XZ', PERIOD=0.8e-06*256, D_CA=0.8e-06*256, DT=1e-12,
        TOTAL_TIMESTEPS=3
    )

    run.setup_run(init_inert_gas=True, init_scraper=False, init_injectors=False)

    surfaceemitter = emission.ZPlaneEmitter(
        conductor=run.cathode, T=run.CATHODE_TEMP, ymin=0.0, ymax=0.0,
        transverse_fac=1.0,
    )
    volemitter_uniform = emission.UniformDistributionVolumeEmitter(
        T=1550, zmin=0, zmax=run.D_CA,
    )
    volemitter = emission.ZSinDistributionVolumeEmitter(
        T=1550, zmin=0, zmax=run.D_CA,
    )
    volemitter_gauss = emission.XGaussZSinDistributionVolumeEmitter(
        T=1550, zmin=0, zmax=run.D_CA, x_sigma=100e-6
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

    # Four real injectors
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
    # Guassian emitter
    emission.PlasmaInjector(
        emitter=volemitter_gauss, species1=run.electrons, species2=run.ions,
        npart=2000, plasma_density=1e19,
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
        res_dict['x'] = np.concatenate(mwxrun.sim_ext.get_particle_x(sname), axis=0)
        res_dict['y'] = np.concatenate(mwxrun.sim_ext.get_particle_y(sname), axis=0)
        res_dict['z'] = np.concatenate(mwxrun.sim_ext.get_particle_z(sname), axis=0)
        res_dict['ux'] = np.concatenate(mwxrun.sim_ext.get_particle_ux(sname), axis=0)
        res_dict['uy'] = np.concatenate(mwxrun.sim_ext.get_particle_uy(sname), axis=0)
        res_dict['uz'] = np.concatenate(mwxrun.sim_ext.get_particle_uz(sname), axis=0)
        res_dict['w'] = np.concatenate(mwxrun.sim_ext.get_particle_weight(sname), axis=0)

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
    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    np.random.seed(11874889)

    run = diode_setup.DiodeRun_V1(
        GEOM_STR='XZ', PERIOD=0.8e-06*256, D_CA=0.8e-06*256, DT=1e-12,
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
        res_dict['x'] = np.concatenate(mwxrun.sim_ext.get_particle_x(sname), axis=0)
        res_dict['y'] = np.concatenate(mwxrun.sim_ext.get_particle_y(sname), axis=0)
        res_dict['z'] = np.concatenate(mwxrun.sim_ext.get_particle_z(sname), axis=0)
        res_dict['ux'] = np.concatenate(mwxrun.sim_ext.get_particle_ux(sname), axis=0)
        res_dict['uy'] = np.concatenate(mwxrun.sim_ext.get_particle_uy(sname), axis=0)
        res_dict['uz'] = np.concatenate(mwxrun.sim_ext.get_particle_uz(sname), axis=0)
        res_dict['w'] = np.concatenate(mwxrun.sim_ext.get_particle_weight(sname), axis=0)

        label_base = f'{species.name}_'
        df[label_base + 'npart'] = npart_dict[sname]

        # Compare main results
        for key, val in res_dict.items():
            df[label_base + key + '_min'] = np.min(val)
            df[label_base + key + '_max'] = np.max(val)
            df[label_base + key + '_mean'] = np.mean(val)
            df[label_base + key + '_std'] = np.std(val)

    assert testing_util.test_df_vs_ref(testname=name, df=df, margin=0.3)


def test_arbitrary_distribution_emitter():
    name = "arbitrary_distribution_emitter"
    # Include a random run number to allow parallel runs to not collide.  Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed(11874889)

    NX = NZ = 512
    NPPC = 10

    run = diode_setup.DiodeRun_V1(
        GEOM_STR='XZ', PERIOD=0.8e-06 * NX, D_CA=0.8e-06 * NZ, DT=1e-12,
        TOTAL_TIMESTEPS=1, DIAG_STEPS=1
    )

    run.setup_run(init_inert_gas=True, init_scraper=False, init_injectors=False)

    z = np.arange(0, NZ) + 0.5
    x = np.arange(0, NX) + 0.5

    zz, xx = np.meshgrid(z, x)

    plasma_density = np.pi/2 * 1e13 * (
        np.sin(np.pi / NX * xx) * np.sin(np.pi / NZ * zz))

    emitter = emission.ArbitraryDistributionVolumeEmitter(
        d_grid=plasma_density,
        T=1550, zmin=0, zmax=run.D_CA,
    )

    emission.PlasmaInjector(
        emitter=emitter, species1=run.electrons, species2=run.ions,
        npart=NX*NZ*NPPC
    )

    diags.FieldDiagnostic(
        diag_steps=1, style='roelof', save_pdf=False,
        plot=True,
    )

    diags.TextDiag(
        diag_steps=1, preset_string='perfdebug'
    )

    run.init_warpx()
    mwxrun.simulation.step(2)

    electron_density = np.load(os.path.join(
        os.curdir, "diags", "fields",
        "electrons_particle_density_0000000002.npy"
    ))
    electron_xavg = np.mean(electron_density, axis=0)
    electron_zavg = np.mean(electron_density, axis=1)

    ion_density = np.load(os.path.join(
        os.curdir, "diags", "fields",
        "ar_ions_particle_density_0000000002.npy"
    ))
    ion_xavg = np.mean(ion_density, axis=0)
    ion_zavg = np.mean(ion_density, axis=1)

    sin_xavg = np.mean(plasma_density, axis=0)
    sin_zavg = np.mean(plasma_density, axis=1)

    assert np.allclose(electron_xavg[10:502], sin_xavg[10:502], rtol=0.06)
    assert np.allclose(electron_zavg[10:502], sin_zavg[10:502], rtol=0.06)
    assert np.allclose(ion_xavg[10:502], sin_xavg[10:502], rtol=0.06)
    assert np.allclose(ion_zavg[10:502], sin_zavg[10:502], rtol=0.06)
