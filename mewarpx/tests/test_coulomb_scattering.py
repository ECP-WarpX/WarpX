"""Tests for functionality in coulomb_scattering.py"""
import numpy as np
import pytest
from pywarpx import callbacks, picmi
import yt

from mewarpx import coulomb_scattering, diags, emission, mespecies
from mewarpx.mwxrun import mwxrun
from mewarpx.poisson_solvers import DummyPoissonSolver
from mewarpx.utils_store import mwxconstants, testing_util
from mewarpx.utils_store import util as mwxutil


def test_coulomb_scattering():

    name = "coulomb_scattering_langevin"
    # Include a random run number to allow parallel runs to not collide. Using
    # python randint prevents collisions due to numpy rseed below
    testing_util.initialize_testingdir(name)

    # Initialize each run with consistent, randomly-chosen, rseed. Use a random
    # seed instead for initial dataframe generation.
    # np.random.seed()
    np.random.seed(92160000)


    class BeamEmitter(emission.UniformDistributionVolumeEmitter):

        """Child class of UniformDistributionVolumeEmitter used to overwrite the
        get_velocities function when sampling velocities for the electrons since
        we want to inject them with a delta function."""

        def __init__(self, T, v_beam):

            self.v_beam = v_beam

            super(BeamEmitter, self).__init__(T=T)

        def _get_xv_coords(self, npart, m, rseed):
            """Get velocities and call specialized function for position."""
            x_coords = self._get_x_coords(npart)

            if m == 1:
                # ion velocities can be Maxwellian
                v_coords = mwxutil.get_velocities(
                    npart, self.T, m=m, emission_type='random', rseed=None)
            else:
                # electron velocities should be a delta function, we need to
                # have small non-zero values for x and y otherwise the Coulomb
                # scattering algorithm fails due to divide by zero errors
                v_coords = np.zeros((3, npart))
                v_coords[0] = np.random.normal(0.0, 1.0, npart)
                v_coords[1] = np.random.normal(0.0, 1.0, npart)
                v_coords[2] = self.v_beam

            return (
                x_coords[:, 0], x_coords[:, 1], x_coords[:, 2],
                v_coords[0], v_coords[1], v_coords[2]
            )


    class VarianceTracker(diags.WarpXDiagnostic):

        def __init__(self, total_steps, diag_steps):
            """Diagnostic to record the electron velocity variance at every
            diagnostic step."""
            self.i = 0
            self.results_array = np.zeros((total_steps//diag_steps, 2))

            super(VarianceTracker, self).__init__(diag_steps)

            callbacks.installafterstep(self._record_variance)

        def _record_variance(self):

            if not self.check_timestep():
                return

            ux = np.concatenate(mwxrun.sim_ext.get_particle_ux('electrons'))

            self.results_array[self.i, 0] = mwxrun.get_t()
            self.results_array[self.i, 1] = np.std(ux)
            self.i += 1


    #######################################################################
    # Simulation parameters                                               #
    #######################################################################

    ZMAX = 250e-6
    NZ = 256
    NX = 8

    DT = 1e-13

    MAX_STEPS = 150
    DIAG_STEPS = 15

    SEED_DENSITY = 5e17
    SEED_COUNT = 20000
    LOGLAMBDA = 7.5

    VBEAM = 1e6

    # Derived quantities

    PERIOD = ZMAX / NZ * NX

    #######################################################################
    # Set geometry, boundary conditions and timestep                      #
    #######################################################################

    mwxrun.init_grid(
        lower_bound=[-PERIOD/2.0, 0.], upper_bound=[PERIOD/2.0, ZMAX],
        number_of_cells=[NX, NZ],
        bc_fields_z_min='periodic',
        bc_fields_z_max='periodic',
        bc_particles_z_min='periodic',
        bc_particles_z_max='periodic'
    )
    mwxrun.init_timestep(DT=DT)
    mwxrun.simulation.max_steps = MAX_STEPS

    #######################################################################
    # Dummy Poisson solver to effectively turn off the field solve        #
    #######################################################################

    solver = DummyPoissonSolver(grid=mwxrun.grid)
    mwxrun.simulation.solver = solver

    #######################################################################
    # Particle types setup                                                #
    #######################################################################

    electrons = mespecies.Species(
        particle_type='electron', name='electrons'
    )
    # artificially increase ion mass to match Lorentz gas better
    ions = mespecies.Species(
        particle_type='Xe', name='xe_ions', charge='q_e', mass=1.0
    )

    #######################################################################
    # Seed simulation with quasi-neutral plasma                           #
    #######################################################################

    volume_emitter = BeamEmitter(1.0, VBEAM)
    emission.PlasmaInjector(
        volume_emitter, electrons, ions, npart=SEED_COUNT*2,
        plasma_density=SEED_DENSITY
    )

    #######################################################################
    # Coulomb scattering                                                  #
    #######################################################################

    langevin = coulomb_scattering.LangevinElectronIonScattering(
        electron_species=electrons, ion_species=ions,
        log_lambda=LOGLAMBDA, subcycling_steps=1
    )

    #######################################################################
    # Diagnostics                                                         #
    #######################################################################

    diags.TextDiag(MAX_STEPS // 5, preset_string='perfdebug')

    variance_tracker = VarianceTracker(MAX_STEPS, DIAG_STEPS)

    #######################################################################
    # Run simulation                                                      #
    #######################################################################

    mwxrun.init_run()
    mwxrun.simulation.step()

    #######################################################################
    # Compare results to theory                                           #
    #######################################################################

    D = (
        SEED_DENSITY * mwxconstants.e**4 * LOGLAMBDA
        / (4 * np.pi * mwxconstants.epsilon_0**2 * mwxconstants.m_e**2 * VBEAM)
    )
    times = variance_tracker.results_array[:, 0]
    calculated_values = variance_tracker.results_array[:, 1] / VBEAM
    expected_values = np.sqrt(D*times) / VBEAM

    print("Calculated ratio: ", calculated_values )
    print("Expected ratio: ", expected_values)
    assert np.allclose(calculated_values, expected_values, rtol=0.01)
    '''
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    times_full = np.linspace(0, times[-1]*1.1, 100)
    ax.plot(times_full*1e12, np.sqrt(D*times_full) / VBEAM, 'k--')
    ax.plot(times*1e12, variance_tracker.results_array[:, 1] / VBEAM, 'o')
    ax.set_ylabel(r'Standard deviation of $v_x$ / $v_0$')
    ax.set_xlabel('Time (ps)')
    ax.grid()
    plt.savefig('benchmark_plot.png', dpi=300)
    plt.show()
    '''
