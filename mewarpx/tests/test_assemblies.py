"""Test the mewarpx wrapper for embedded boundaries. This test is the same
as the test in warpx/Examples/Tests/ElectrostaticSphereEB/inputs_3d but in 2d.
"""
import numpy as np
import os
import matplotlib.pyplot as plt

from mewarpx.utils_store import util as mwxutil


def test_embedded_cylinder():
    name = "Embedded_cylinder_solve"
    dim = 2

    # Initialize and import only when we know dimension
    mwxutil.init_libwarpx(ndim=dim, rz=False)
    from mewarpx import assemblies
    from mewarpx.setups_store import diode_setup
    from mewarpx.utils_store import testing_util
    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(92160881)

    # Specific numbers match older run for consistency
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
    cylinder = assemblies.Cylinder(
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
    from mewarpx.setups_store import diode_setup
    from mewarpx.utils_store import testing_util
    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(92160881)

    # Specific numbers match older run for consistency
    D_CA = 1  # m
    NX = 64
    NZ = 64
    run = diode_setup.DiodeRun_V1(
        dim=dim,
        D_CA=D_CA,
        NX=NX,
        NZ=NZ,
        # This gives equal spacing in x & z
        PERIOD=D_CA * NX / NZ,
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
    from mewarpx.setups_store import diode_setup
    from mewarpx.utils_store import testing_util
    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(42578825)

    # Specific numbers match older run for consistency
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
    cylinder1 = assemblies.Cylinder(
        center_x=-0.25, center_z=0.5, radius=0.1,
        V=-2.0, T=300, WF=4.7, name="Cylinder1"
    )
    cylinder2 = assemblies.Cylinder(
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
    from mewarpx.setups_store import diode_setup
    from mewarpx.utils_store import testing_util
    from mewarpx.mwxrun import mwxrun

    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(42147825)

    # Specific numbers match older run for consistency
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
        DT=1e-9,
        TOTAL_TIMESTEPS=2,
        DIAG_STEPS=2,
        FIELD_DIAG_DATA_LIST=['phi'],
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
    cylinder1 = assemblies.Cylinder(
        center_x=-0.25, center_z=0.5, radius=0.1,
        V=-2.0, T=300, WF=4.7, name="Cylinder1"
    )
    cylinder2 = assemblies.Cylinder(
        center_x=0.25, center_z=0.5, radius=0.1,
        V=5.0, T=300, WF=4.7, name="Cylinder2"
    )

    # Initialize solver
    run.init_solver()
    run.init_conductors()
    run.electrons.save_particles_at_eb = 1
    run.ions.save_particles_at_eb = 1

    # Inject particles in the simulation
    volemitter = emission.ZSinDistributionVolumeEmitter(
        T=1550, zmin=0, zmax=run.D_CA,
    )
    emission.PlasmaInjector(
        emitter=volemitter, species1=run.electrons, species2=run.ions,
        npart=2000, plasma_density=1e19,
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

    cylinder2.init_scrapedparticles(cylinder1.fields)
    cylinder2.record_scrapedparticles()
    cyl2_scraped = cylinder2.get_scrapedparticles()

    assert np.all(cyl1_scraped['n'] == 48)
    assert np.all(cyl2_scraped['n'] == 42)

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
