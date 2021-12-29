"""Test the mewarpx wrapper for embedded boundaries. This test is the same
as the test in warpx/Examples/Tests/ElectrostaticSphereEB/inputs_3d but in 2d.
"""
import os

import matplotlib.pyplot as plt
import numpy as np

from mewarpx.utils_store import util as mwxutil


def test_embedded_cylinder():
    name = "Embedded_cylinder_solve"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import assemblies
    from mewarpx.mwxrun import mwxrun
    from mewarpx.setups_store import diode_setup
    from mewarpx.utils_store import testing_util

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(92160881)

    # Specific numbers match older run for consistency
    D_CA = 1  # m
    run = diode_setup.DiodeRun_V1(
        dim=dim,
        #V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        NX=64,
        NZ=64,
        DT=1e-6,
        TOTAL_TIMESTEPS=1,
        DIAG_STEPS=1,
        FIELD_DIAG_DATA_LIST=['phi'],
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_scraper=False,
        init_electrons=False,
        init_solver=False,
        init_injectors=False,
        init_field_diag=True,
        init_simcontrol=True,
        init_simulation=False
    )

    # Install the embedded boundary
    cylinder = assemblies.InfCylinderY(
        center_x=0.0, center_z=0.5, radius=0.1,
        V=-2.0, T=300, WF=4.7, name="Cylinder"
    )

    # Initialize solver
    run.init_solver()
    run.init_conductors()

    # Initialize the simulation
    run.init_simulation()
    run.init_warpx()

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Check phi results against reference data                            #
    #######################################################################

    phi = mwxrun.get_gathered_phi_grid(include_ghosts=False)
    # np.save('embedded_cylinder_phi.npy', phi)
    ref_phi = np.load(os.path.join(
        testing_util.test_dir, 'embedded_boundary', 'embedded_cylinder_phi.npy'
    ))
    assert np.allclose(phi, ref_phi, rtol=0.001)


def test_embedded_rectangle():
    name = "Embedded_rectangle_solve"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import assemblies
    from mewarpx.mwxrun import mwxrun
    from mewarpx.setups_store import diode_setup
    from mewarpx.utils_store import testing_util

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(92160881)

    # Specific numbers match older run for consistency
    D_CA = 1  # m
    run = diode_setup.DiodeRun_V1(
        D_CA=D_CA,
        NX=64,
        NZ=64,
        DT=1e-6,
        TOTAL_TIMESTEPS=1,
        DIAG_STEPS=1,
        FIELD_DIAG_DATA_LIST=['phi'],
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_scraper=False,
        init_electrons=False,
        init_solver=False,
        init_injectors=False,
        init_field_diag=True,
        init_simcontrol=True,
        init_simulation=False
    )

    # Install the embedded boundary
    cylinder = assemblies.Rectangle(
        center_x=0.0, center_z=0.5, length_x=0.3, length_z=0.3,
        V=-2.0, T=300, WF=4.7, name="Box"
    )

    # Initialize solver
    run.init_solver()
    run.init_conductors()

    # Initialize the simulation
    run.init_simulation()
    run.init_warpx()

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Check phi results against reference data                            #
    #######################################################################

    phi = mwxrun.get_gathered_phi_grid(include_ghosts=False)
    # np.save('embedded_rectangle_phi.npy', phi)
    ref_phi = np.load(os.path.join(
        testing_util.test_dir, 'embedded_boundary', 'embedded_rectangle_phi.npy'
    ))

    assert np.allclose(phi, ref_phi, rtol=0.001)


def test_two_embedded_cylinders():
    name = "two_embedded_cylinders_solve"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import assemblies
    from mewarpx.mwxrun import mwxrun
    from mewarpx.setups_store import diode_setup
    from mewarpx.utils_store import testing_util

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(42578825)

    # Specific numbers match older run for consistency
    D_CA = 1  # m
    run = diode_setup.DiodeRun_V1(
        #V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        NX = 64,
        NZ = 64,
        DT=1e-6,
        TOTAL_TIMESTEPS=1,
        DIAG_STEPS=1,
        FIELD_DIAG_DATA_LIST=['phi'],
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_scraper=False,
        init_electrons=False,
        init_solver=False,
        init_injectors=False,
        init_field_diag=True,
        init_simcontrol=True,
        init_simulation=False
    )

    # Install the embedded boundaries
    cylinder1 = assemblies.InfCylinderY(
        center_x=-0.25, center_z=0.5, radius=0.1,
        V=-2.0, T=300, WF=4.7, name="Cylinder1"
    )
    cylinder2 = assemblies.InfCylinderY(
        center_x=0.25, center_z=0.5, radius=0.1,
        V=5.0, T=300, WF=4.7, name="Cylinder2"
    )

    # Initialize solver
    run.init_solver()
    run.init_conductors()

    # Initialize the simulation
    run.init_simulation()
    run.init_warpx()

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Check phi results against reference data                            #
    #######################################################################

    phi = mwxrun.get_gathered_phi_grid(include_ghosts=False)
    # np.save('two_embedded_cylinders_phi.npy', phi)
    ref_phi = np.load(os.path.join(
        testing_util.test_dir, 'embedded_boundary',
        'two_embedded_cylinders_phi.npy'
    ))
    assert np.allclose(phi, ref_phi, rtol=0.001)


def test_two_embedded_cylinders_scraping():
    name = "two_embedded_cylinders_scraping"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import assemblies, emission
    from mewarpx.mwxrun import mwxrun
    from mewarpx.setups_store import diode_setup
    from mewarpx.utils_store import testing_util

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(42147820)

    # Specific numbers match older run for consistency
    D_CA = 0.025  # m
    run = diode_setup.DiodeRun_V1(
        #V_ANODE_CATHODE=VOLTAGE,
        D_CA=D_CA,
        NX=64,
        NZ=64,
        DT=1e-10,
        TOTAL_TIMESTEPS=15,
        DIAG_STEPS=15,
        FIELD_DIAG_DATA_LIST=['phi'],
        # FIELD_DIAG_PLOT=True,
        INERT_GAS_TYPE='positron'
    )
    # Only the functions we change from defaults are listed here
    run.setup_run(
        init_conductors=False,
        init_scraper=False,
        init_electrons=True,
        init_inert_gas=True,
        init_solver=False,
        init_injectors=False,
        init_field_diag=True,
        init_simcontrol=True,
        init_simulation=False
    )

    # Install the embedded boundaries
    cylinder1 = assemblies.InfCylinderY(
        center_x=-0.25*D_CA, center_z=0.5*D_CA, radius=0.1*D_CA,
        V=-0.5, T=300, WF=4.7, name="Cylinder1"
    )
    cylinder2 = assemblies.InfCylinderY(
        center_x=0.25*D_CA, center_z=0.5*D_CA, radius=0.1*D_CA,
        V=0.2, T=300, WF=4.7, name="Cylinder2"
    )

    # Initialize solver
    run.init_solver()
    run.init_conductors()
    run.electrons.save_particles_at_eb = 1
    run.ions.save_particles_at_eb = 1

    # Inject particles in the simulation
    volemitter = emission.ZSinDistributionVolumeEmitter(
        T=3000, zmin=0, zmax=run.D_CA,
    )
    emission.PlasmaInjector(
        emitter=volemitter, species1=run.electrons, species2=run.ions,
        npart=4000, plasma_density=1e14,
    )

    # Initialize the simulation
    run.init_simulation()
    run.init_warpx()

    # Run the main WARP loop
    while run.control.check_criteria():
        mwxrun.simulation.step()

    #######################################################################
    # Check flux on each cylinder                                         #
    #######################################################################

    cylinder1.init_scrapedparticles(cylinder1.fields)
    cylinder1.record_scrapedparticles()
    cyl1_scraped = cylinder1.get_scrapedparticles()

    cylinder2.init_scrapedparticles(cylinder2.fields)
    cylinder2.record_scrapedparticles()
    cyl2_scraped = cylinder2.get_scrapedparticles()

    assert np.allclose(cyl1_scraped['n'], np.array([3, 0]))
    assert np.allclose(cyl2_scraped['n'], np.array([1, 2]))

    #######################################################################
    # Check rho results against reference data                            #
    #######################################################################

    rho = mwxrun.get_gathered_rho_grid(include_ghosts=False)[:,:,0]
    # np.save('two_embedded_cylinders_rho.npy', rho)
    ref_rho = np.load(os.path.join(
        testing_util.test_dir, 'embedded_boundary',
        'two_embedded_cylinders_rho.npy'
    ))
    assert np.allclose(rho, ref_rho)


def test_infinite_cylinder_z():
    name = "infinite_cylinder_z_scraping"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=True)

    # Initialize and import only when we know dimension
    from mewarpx.utils_store import testing_util

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed.
    np.random.seed(42147820)

    import sys
    sys.path.append(testing_util.example_dir)

    from thermionic_diode_rz_cylinder import CylinderVacuumTEC

    run = CylinderVacuumTEC(
        V_ANODE_CATHODE=5.0, TOTAL_TIMESTEPS=1000, SAVE=True, USE_EB=True,
        DIAG_STEPS=500
    )
    run.NR = 256

    run.setup_run()
    run.run_sim()

    key = ('scrape', 'anode', 'electrons')
    J_diode = run.fluxdiag.ts_dict[key].get_averagevalue_by_key('J')
    assert np.isclose(J_diode, 1.5154551964454426)
