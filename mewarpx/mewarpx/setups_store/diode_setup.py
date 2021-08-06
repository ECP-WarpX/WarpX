"""Framework that sets up a 1D or 2D run based on input parameters, with
defaults available.
"""
from pywarpx import picmi

from mewarpx.mwxrun import mwxrun
from mewarpx.sim_control import SimControl
from mewarpx import mcc_wrapper, poisson_pseudo_1d, emission, assemblies, mepicmi
from mewarpx.diags_store import field_diagnostic, particle_diagnostic


class DiodeRun_V1(object):

    """A combination of settings and initialization functions to standardize
    unit tests (hopefully).

    Explicitly versioned with V1 so that if API is changed this can be updated,
    but if "Best practices" change in a way that would change results, a new
    version can be created while existing unit tests use this.
    """

    # ### ELECTRODES ###
    # Cathode temperature in K
    CATHODE_TEMP = 1550
    # Richardson constant in A/m^2/K^2
    CATHODE_A = 120.17e4
    # Work function in eV.
    CATHODE_PHI = 2.4
    # Anode temperature in K
    ANODE_TEMP = 773
    # Work Function of Anode.
    ANODE_PHI = 1.4
    # Anode vacuum bias relative to cathode
    V_ANODE_CATHODE = 2.0
    # Instead of constant voltage, use an expression. In this case, V_CATHODE
    # is set to 0.
    V_ANODE_EXPRESSION = None

    # ### VACUUM ###
    # Cathode anode distance
    D_CA = 1e-3

    # Gas properties
    # String, eg 'He'
    INERT_GAS_TYPE = 'Ar'
    # Specify N_INERT or P_INERT
    N_INERT = None  # m^-3
    P_INERT = None  # in Torr
    # Neutral temperature. NOTE if either cathode or anode temperatures are
    # changed, T_INERT needs to be changed afterwards too.
    T_INERT = 0.5 * (CATHODE_TEMP + ANODE_TEMP)  # rough approximation

    # ### PHYSICS MODEL SETTINGS ###
    # Beam velocity spread in transverse vs longitudinal
    TRANSVERSE_FAC = 1.1

    # Reflection physics
    # [[[TODO once in WarpX]]]

    # Plasma seed density and electron temperature
    PLASMA_DENSITY = None  # m^-3
    T_ELEC = None  # K
    SEED_NPPC = None # number of particles per cell in the seed plasma

    # Particle Species
    # A list of particle species created
    SPECIES = None

    # ### INJECTION SETTINGS ###
    # Number of particles to inject TOTAL per timestep
    NPARTPERSTEP = 10
    # 2D only: Normal thermionic injection - number of particles to inject per
    # cell per timestep
    NPPC = 2
    # 2D only: Normal thermionic injection - specify PERIOD instead to ignore
    # NPARTPERSTEP and force this period
    PERIOD = None

    # [[[TODO ONCE IN WARPX]]]
    # # Whether to use Schottky injection
    USE_SCHOTTKY = False
    # # Noninteracting run - only inject trace particles and do field solve once
    NONINTERAC = False
    # # Number of particles injected each noninteracting wave
    # NONINTERAC_WAVENUM = 50000

    # ### RUN SETTINGS ###
    # [[[TODO once in warpx]]]
    # # P_CUTOFF is the value at which the simulation is terminated at a
    # # diagnostic interval, if power output is below the cutoff
    # P_CUTOFF = -1e6
    # # J_TOLERANCE is the percentage change allowed in currents between
    # # diag_steps results for which the run will terminate rather than continue
    # # (ie the acceptable percentage deviation for convergence). Set to None to
    # # disallow.
    # J_TOLERANCE = None
    # # J_TOTAL_TOLERANCE disallows termination due to J convergence if there's
    # # still a net current (eg emitted currents != collected currents) above a
    # # given threshold. Presently defaults to 1e-3 if not specified here.
    # J_TOTAL_TOLERANCE = None
    # # CFL_FACTOR determines the time step at which 99% of particles should
    # # satisfy the CFL condition. 1.0 would be typical; usually we use 0.4 for a
    # # safety margin
    # CFL_FACTOR = 0.4
    # RES_LENGTH is the cell size in meters. It always reflects the z spacing,
    # so is overwritten if NZ is given. If NX and/or NY are specified, the cell
    # spacing in those dimensions is not guaranteed to match RES_LENGTH.
    RES_LENGTH = 0.8e-06
    # MAX_GRID_SIZE_FACTOR is the factor by which nz (the number of cells
    # in the z-direction) is reduced to obtain the max grid size parameter
    # for Cartesian grids.
    MAX_GRID_SIZE_FACTOR = 4
    # OFFSET pads the cathode and anode by extending the simulation limits.
    # However, if the boundaries are serving as the conductors this should be
    # 0.
    OFFSET = 0.0
    # # Check Debye Length allows disabling this check to run with large
    # # resolutions
    # CHECK_DEBYE_LENGTH = True
    # DT allows the auto-calculated DT to be overridden. CFL_FACTOR is ignored
    # in this case.
    DT = None
    # # CHECK_CHARGE_CONSERVATION is a heuristic only for detecting (usually)
    # # violations of CFL. However, it can also be triggered during short
    # # diagnostic periods that cover transient behavior, and should thus often be
    # # disabled in short unit tests.
    # CHECK_CHARGE_CONSERVATION = True
    # Get diagnostic information every DIAG_STEPS setps.
    DIAG_STEPS = None
    # The data list for field diagnostics
    FIELD_DIAG_DATA_LIST = None
    # Whether or not to plot parameters on each diagnostic step.
    FIELD_DIAG_PLOT_ON_DIAG_STEPS = False
    # If FIELD_DIAG_PLOT_ON_DIAG_STEPS is True, what parameters to plot. Current
    # accepted values are 'rho' and 'phi'
    FIELD_DIAG_PLOT_DATA_LIST = None
    # Whether or not to plot the parameters in FIELD_DIAG_DATA_LIST after
    # the simulation is complete using yt files generated during the run at
    # diagnostic steps.
    FIELD_DIAG_PLOT_AFTER_RUN = False
    # The data list for particle diagnostics
    PARTICLE_DIAG_DATA_LIST = None
    # Species to be plotted from the particle diagnostics
    PARTICLE_PLOT_SPECIES = None
    # Data to be plotted for each selected species from the particle diagnostics
    PARTICLE_DIAG_PLOT_DATA_LIST = None
    # Whether or not to load yt particle diagnostic data and save plots in
    # post-processing
    PARTICLE_DIAG_PLOT_AFTER_RUN = False

    # The total timesteps of the simulation
    TOTAL_TIMESTEPS = None
    # number of cells in the x direction
    NX = None
    # number of cells in the y direction
    NY = None
    # number of cells in the z direction
    NZ = None
    # should the direct solver be used?
    DIRECT_SOLVER = False

    def __init__(self, dim=1, rz=False, **kwargs):
        for kw in list(kwargs.keys()):
            setattr(self, kw, kwargs[kw])

        self.dim = dim
        if self.dim not in [1, 2, 3]:
            raise ValueError(f"Unavailable dimension {self.dim}")
        self.rz = rz
        if self.dim != 2 and self.rz:
            raise ValueError("dim=2 required for RZ")

        if self.dim > 1:
            if self.PERIOD is None:
                raise ValueError("Derived period not yet supported.")
                # self.PERIOD = (
                #     np.ceil(float(self.NPARTPERSTEP)/self.NPPC)
                #     * self.RES_LENGTH
                # )

        if self.V_ANODE_EXPRESSION is None:
            # V_CATHODE is set so that its Fermi level is at V=0.
            self.V_CATHODE = -self.CATHODE_PHI
            self.V_ANODE = self.V_CATHODE + self.V_ANODE_CATHODE
            self.V_ANODE_EXPRESSION = self.V_ANODE
        else:
            self.V_CATHODE = 0.0
            self.V_ANODE = self.V_ANODE_EXPRESSION

        print('Setting up simulation with')
        print(f'  dt = {self.DT:.3e} s')
        print(f'  Total time = {self.TOTAL_TIMESTEPS * self.DT:.3e} '
              f's ({self.TOTAL_TIMESTEPS} timesteps)')
        print('  Diag time = {self.DIAG_STEPS * self.DT:.3e} s '
              f'({self.DIAG_STEPS} timesteps)')

    def setup_run(
        self,
        init_base=True,
        init_electrons=True,
        init_solver=True,
        init_conductors=True,
        init_scraper=True,
        init_injectors=True,
        init_neutral_plasma=False,
        init_mcc=False,
        init_reflection=False,
        init_inert_gas=False,
        init_merging=False,
        init_traceparticles=False,
        init_runinfo=False,
        init_fluxdiag=False,
        init_field_diag=False,
        init_particle_diag=False,
        init_simcontrol=False,
        init_simulation=True,
        init_warpx=False
    ):
        """Perform each part of setup if requested.

        Ones marked True are generally necessary for a standard run.
        init_simulation and init_warpx are necessary, but must come after any
        other user customization, so by default we keep them False.
        """
        if init_base:
            self.init_base()
        if init_electrons:
            self.init_electrons()
        if init_inert_gas:
            self.init_inert_gas()
        if init_solver:
            self.init_solver()
        if init_conductors:
            self.init_conductors()
        if init_neutral_plasma:
            self.init_neutral_plasma()
        if init_mcc:
            self.init_MCC()
        if init_scraper:
            self.init_scraper()
        if init_reflection:
            self.init_reflection()
        if init_merging:
            self.init_merging()
        if init_traceparticles:
            self.init_traceparticles()
        if init_runinfo:
            self.init_runinfo()
        if init_fluxdiag:
            self.init_fluxdiag()
        if init_simcontrol:
            self.init_simcontrol()
        if init_simulation:
            self.init_simulation()
        if init_field_diag:
            self.init_field_diag()
        if init_particle_diag:
            self.init_particle_diag()
        if init_warpx:
            self.init_warpx()
        if init_injectors:
            self.init_injectors()

    def init_base(self):
        print('### Init Diode Base Setup ###', flush=True)
        # Set grid boundaries
        if self.dim == 1:
            # Translational symmetry in x & y, 1D simulation in z. Note warp
            # immediately sets xmin and ymin to 0 when
            # warp.w3d.initdecompositionw3d() is called, so we work with that
            # offset here, though nothing should depend on it.
            xmin = 0
            xmax = 1.0
        else:
            xmin = -self.PERIOD/2.0
            xmax = self.PERIOD/2.0

        if self.dim < 3 and not self.rz:
            ymin = 0
            ymax = 1.0
        else:
            ymin = -self.PERIOD/2.0
            ymax = self.PERIOD/2.0

        zmin = -self.OFFSET
        zmax = self.D_CA + self.OFFSET

        # Grid parameters - set grid counts

        # If NZ was passed in recalculate RES_LENGTH
        # otherwise calculate NX, NY, and NZ
        if self.NZ is not None:
            self.RES_LENGTH = (zmax-zmin)/self.NZ
        else:
            self.NZ = int(round((zmax - zmin)/self.RES_LENGTH))

        if self.NX is None:
            if self.dim == 1:
                self.NX = 0
            else:
                self.NX = int(round(self.PERIOD/self.RES_LENGTH))
        if self.NY is None:
            if self.dim < 3 and not self.rz:
                self.NY = 0
            else:
                self.NY = int(round(self.PERIOD/self.RES_LENGTH))

        print(
            f"Creating grid with NX={self.NX}, NY={self.NY}, NZ={self.NZ} and "
            f"x, y, z limits of [[{xmin:.4g}, {xmax:.4g}], [{ymin:.4g}, "
            f"{ymax:.4g}], [{zmin:.4g}, {zmax:.4g}]]"
        )

        # create the grid
        if self.dim == 1:
            raise NotImplementedError(
                "1D grid is not yet implemented in mewarpx")
        elif self.dim == 2:
            mwxrun.init_grid(
                xmin, xmax, zmin, zmax, self.NX, self.NZ,
                max_grid_size=self.MAX_GRID_SIZE_FACTOR
            )
        elif self.dim == 3:
            raise NotImplementedError("3d grid not yet implemented.")
            '''
            self.grid = picmi.Cartesian3DGrid(
                number_of_cells=[self.NX, self.NY, self.NZ],
                lower_bound=[xmin, ymin, zmin],
                upper_bound=[xmax, ymax, zmax],
                lower_boundary_conditions=['periodic', 'periodic', 'dirichlet'],
                upper_boundary_conditions=['periodic', 'periodic', 'dirichlet'],
                warpx_potential_lo_z=self.V_CATHODE,
                warpx_potential_hi_z=self.V_ANODE_EXPRESSION,
                lower_boundary_conditions_particles=[
                    'periodic', 'periodic', 'absorbing'],
                upper_boundary_conditions_particles=[
                    'periodic', 'periodic', 'absorbing'],
                moving_window_velocity=None,
                warpx_max_grid_size=self.NZ//self.MAX_GRID_SIZE_FACTOR
            )
            '''

        #######################################################################
        # Run setup calculations and print diagnostic info                    #
        #######################################################################

        # [[[TODO once implemented in WarpX]]] Most of this can be commented
        # out. For now diag_steps and total_timesteps may need to be class
        # input variables as other things above are.

        # Must be run after warp boundaries (xmin etc.) and cell counts (nx
        # etc.) are set.
        # Full decomp always used for 1D parallel runs
        # Particle decomposition - each processor keeps particles across full
        # domain
        # No field decomposition - each processor solves field itself

        # self.setupinfo = gentools.SetupCalcs(
        #     cathode_temp=self.CATHODE_TEMP,
        #     cathode_phi=self.CATHODE_PHI,
        #     cathode_A=self.CATHODE_A,
        #     V_anode=self.V_ANODE_CATHODE,
        #     V_grid=self.V_ANODE_CATHODE,
        #     D_CG=self.D_CA,
        #     CFL_factor=self.CFL_FACTOR,
        #     transverse_fac=self.TRANSVERSE_FAC,
        #     check_debye_length=self.CHECK_DEBYE_LENGTH,
        #     dt=self.DT,
        #     part_full_decomp=self.PART_FULL_DECOMP,
        #     no_fs_decomp=self.NO_FS_DECOMP,
        #     use_B=False,
        # )

    def init_electrons(self):
        print('### Init Electrons ###', flush=True)
        self.electrons = mepicmi.Species(
            particle_type='electron', name='electrons'
        )

    def init_inert_gas(self):
        print('### Init Inert Gas Ions ###', flush=True)
        if self.INERT_GAS_TYPE not in ['He', 'Ar', 'Xe']:
            raise NotImplementedError(
                f"Inert gas is not yet implemented in mewarpx with "
                f"INERT_GAS_TYPE = {self.INERT_GAS_TYPE}"
            )

        self.ions = mepicmi.Species(
            particle_type=self.INERT_GAS_TYPE,
            name=f'{self.INERT_GAS_TYPE.lower()}_ions',
            charge='q_e'
        )

    def init_solver(self):
        print('### Init Diode Solver Setup ###', flush=True)
        if self.dim == 1:
            raise NotImplementedError(
                "1D solving is not yet implemented in mewarpx")
            # self.solver = poisson1d.PoissonSolver1D()
        elif self.dim in [2, 3]:
            if self.DIRECT_SOLVER:
                self.solver = poisson_pseudo_1d.PoissonSolverPseudo1D(
                    grid=mwxrun.grid
                )
            else:
                self.solver = picmi.ElectrostaticSolver(
                    grid=mwxrun.grid,
                    method='Multigrid',
                    required_precision=1e-6,
                    maximum_iterations=10000
                )
        self.solver.self_fields_verbosity = 2 if self.NONINTERAC else 0

    def init_conductors(self):
        print('### Init Diode Conductors Setup ###', flush=True)
        self.cathode = assemblies.Cathode(
            V=self.V_CATHODE, T=self.CATHODE_TEMP, WF=self.CATHODE_PHI
        )

        self.anode = assemblies.Anode(
            z=self.D_CA, V=self.V_ANODE, T=self.ANODE_TEMP, WF=self.ANODE_PHI
        )

        self.surface_list = [self.cathode, self.anode]

        # Set solver boundary potentials
        mwxrun.grid.potential_zmin = self.cathode.V
        mwxrun.grid.potential_zmax = self.anode.V

    def init_scraper(self):
        raise NotImplementedError(
            "Diode scraper is not yet implemented in mewarpx")
        print('### Init Diode Scraper Setup ###', flush=True)
        # profile_decorator = None
        # if self.PROFILE_DIAG:
        #     profile_decorator = util.get_runtime_profile_decorator()

        # if self.dim == 1 and not self.FORCE_2D_SCRAPER:
        #     self.scraper = particlescraper1d.ParticleScraper1D(
        #         cathode=self.cathode, anode=self.anode_plane,
        #         profile_decorator=profile_decorator
        #     )

        # else:
        #     # Small Aura required to make sure 1 cell (isnear) boundary is
        #     # computed properly for conductors that terminate on a grid line
        #     dmin = min(self.setupinfo.dx, self.setupinfo.dz)/100.
        #     # ltunnel to False because it requires a griddistance method that
        #     # ZSrfrv doesn't have
        #     self.scraper = warp.ParticleScraper(
        #         self.surface_list,
        #         lcollectlpdata=True,
        #         ltunnel=False,
        #         aura=dmin,
        #         nxscale=self.NSCALE,
        #         nyscale=self.NSCALE,
        #         profile_decorator=profile_decorator
        #     )

    def init_injectors(self):
        print('### Init Diode Injectors Setup ###', flush=True)
        if self.rz:
            self.emitter = emission.ZDiscEmitter(
                conductor=self.cathode, T=self.CATHODE_TEMP,
                transverse_fac=self.TRANSVERSE_FAC,
                use_Schottky=self.USE_SCHOTTKY,
            )
        elif self.dim < 3:
            self.emitter = emission.ZPlaneEmitter(
                conductor=self.cathode,
                T=self.CATHODE_TEMP,
                ymin=0.0, ymax=0.0,
                transverse_fac=self.TRANSVERSE_FAC,
                use_Schottky=self.USE_SCHOTTKY,
            )
        else:
            self.emitter = emission.ZPlaneEmitter(
                conductor=self.cathode,
                T=self.CATHODE_TEMP,
                transverse_fac=self.TRANSVERSE_FAC,
                use_Schottky=self.USE_SCHOTTKY,
            )

        if self.NONINTERAC:
            # runtools.set_noninteracting(solvefreq=self.DIAG_STEPS)
            # weight = (gentools.J_RD(
            #     self.CATHODE_TEMP, self.CATHODE_PHI, self.CATHODE_A)
            #     * self.emitter.area
            #     * warp.top.dt * self.DIAG_STEPS
            #     / e / self.NONINTERAC_WAVENUM)
            # self.injector = emission.FixedNumberInjector(
            #     self.emitter, 'beam', npart=self.NONINTERAC_WAVENUM,
            #     injectfreq=self.DIAG_STEPS,
            #     weight=weight
            # )
            raise NotImplementedError
        else:
            if self.dim == 1:
                npart = self.NPARTPERSTEP
            else:
                npart = self.NPPC

            self.injector = emission.ThermionicInjector(
                self.emitter, self.electrons, npart, self.CATHODE_TEMP,
                self.CATHODE_PHI, self.CATHODE_A
            )

    def init_neutral_plasma(self):
        print('### Init Neutral Seed Plasma Setup ###', flush=True)

        self.vol_emitter = emission.UniformDistributionVolumeEmitter(
            T=self.T_ELEC
        )
        self.plasma_injector = emission.PlasmaInjector(
            emitter=self.vol_emitter, species1=self.electrons,
            species2=self.ions, npart=2 * self.SEED_NPPC * self.NX * self.NZ,
            T_2=self.T_INERT, plasma_density=self.PLASMA_DENSITY
        )

    def init_MCC(self):
        print('### Init MCC Setup ###', flush=True)
        if not hasattr(self, "exclude_collisions"):
            self.exclude_collisions = None

        mcc_wrapper.MCC(
            electron_species=self.electrons,
            ion_species=self.ions,
            T_INERT=self.T_INERT,
            P_INERT=self.P_INERT,
            N_INERT=self.N_INERT,
            exclude_collisions=self.exclude_collisions)

    def init_reflection(self):
        raise NotImplementedError(
            "Diode reflection is not yet implemented in mewarpx")
        print('### Init Diode Reflection Setup ###', flush=True)

    def init_merging(self):
        raise NotImplementedError(
            "Diode merging is not yet implemented in mewarpx")
        print('### Init Diode Merging ###', flush=True)

    def init_traceparticles(self):
        raise NotImplementedError(
            "Diode TraceParticles is not yet implemented in mewarpx")
        print('### Init Diode TraceParticles ###', flush=True)

    def init_runinfo(self):
        raise NotImplementedError(
            "Diode Runinfo is not yet implemented in mewarpx")
        print('### Init Diode Runinfo Setup ###', flush=True)

    def init_fluxdiag(self):
        raise NotImplementedError(
            "Diode FluxDiag is not yet implemented in mewarpx")
        print('### Init Diode FluxDiag ###', flush=True)

    def init_field_diag(self):
        print('### Init Diode FieldDiag ###', flush=True)

        self.field_diag = field_diagnostic.FieldDiagnostic(
            name='fields',
            grid=mwxrun.grid,
            diag_steps=self.DIAG_STEPS,
            diag_data_list=self.FIELD_DIAG_DATA_LIST,
            write_dir='diags/fields',
            plot_on_diag_step=self.FIELD_DIAG_PLOT_ON_DIAG_STEPS,
            plot_data_list=self.FIELD_DIAG_PLOT_DATA_LIST,
            post_processing=self.FIELD_DIAG_PLOT_AFTER_RUN
        )

    def init_particle_diag(self):
        print('### Init Diode ParticleDiag', flush=True)

        self.particle_diag = particle_diagnostic.ParticleDiagnostic(
            name='particles',
            diag_steps=self.DIAG_STEPS,
            write_dir='diags/particles',
            data_list=self.PARTICLE_DIAG_DATA_LIST,
            plot_species=self.PARTICLE_PLOT_SPECIES,
            plot_data_list=self.PARTICLE_DIAG_PLOT_DATA_LIST,
            post_processing=self.PARTICLE_DIAG_PLOT_AFTER_RUN
        )

    def init_simcontrol(self):
        print('### Init Diode SimControl ###', flush=True)
        self.control = SimControl(max_steps=self.TOTAL_TIMESTEPS)

    def init_simulation(self):
        print('### Init Simulation Setup ###', flush=True)
        mwxrun.simulation.solver = self.solver
        mwxrun.simulation.time_step_size = self.DT
        mwxrun.simulation.max_steps = self.TOTAL_TIMESTEPS

    def init_warpx(self):
        print('### Init Simulation Run ###', flush=True)
        mwxrun.init_run()
