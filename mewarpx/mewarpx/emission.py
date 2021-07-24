"""
Module for various types of particle emission in WarpX.
"""
# import collections
import logging
import warnings

import numba
import numpy as np

import skimage.measure
from pywarpx import callbacks, _libwarpx, picmi

import mewarpx.util as mwxutil
from mewarpx.mwxrun import mwxrun
import mewarpx.mwxconstants as constants

# Get module-level logger
logger = logging.getLogger(__name__)


class Injector(object):

    """Base class for injection.
    All injectors must include an emitter object, and should also include a
    'name' field for diagnostics.
    """

    emitter = None

    # This is overridden if a diagnostic is installed to record injected
    # current.
    injector_diag = None

    # fields is used by the diags.FluxInjectorDiag to know what to write to
    # the CSV file. It can be overridden by child classes, but is not currently
    # adjustable by the user.

    # IF CHANGING THIS, CHANGE IN self.record_injectedparticles() AS WELL.
    fields = ['t', 'step', 'jsid', 'V_e', 'n', 'q', 'E_total']

    # @staticmethod
    # def setup_warp():
    #     """Stuff that needs to be set before injectors are used."""
    #     # Update warp derived quantities if needed.
    #     warp.derivqty()
    #     # Record E_total
    #     if 'E_total' not in warp.Species._addedpids:
    #         warp.Species.addpid('E_total')

    @staticmethod
    def compute_npart(npart_total, unique_particles):
        """Compute number of particles to insert at a given timestep.

        This function translates between total particle number and this
        processor's particle numbers. If particles are designated "unique",
        none are discarded by WarpX so we have logic here to give the processor
        the right number of particles, with additional logic to load-balance
        the remainder. If unique_particles is False, WarpX essentially does the
        particle discarding, so each processor should inject the whole number
        of particles to start.

        Arguments:
            npart_total (int): Integer number of total particles to insert this
                timestep.
            unique_particles (bool): If True, WarpX keeps all particles sent to
                it. If False, it only keeps a processor's fraction of total
                particles.

        Returns:
            npart (int): Integer number of total particles for this processor
                to insert this timestep.
        """
        if not unique_particles:
            return npart_total

        npart = npart_total // mwxrun.n_procs

        # Early-numbered processors add one additional particle if needed.
        # Particles get re-distributed between processors after injection, so
        # this shouldn't load-imbalance anything.
        if mwxrun.me < (npart_total % mwxrun.n_procs):
            npart += 1

        return npart

    def getvoltage_e(self):
        """Return the electrical voltage of the injector. Defaults to returning
        0, unless an emitter is associated with this injector (it should be) in
        which case return the emitter's electrical voltage.
        Child classes can override this if needed.
        """
        if self.emitter is not None:
            return self.emitter.getvoltage_e()

        return 0.

    def init_injectedparticles(self, fieldlist):
        """Set up the injected particles array. Call before
        append_injectedparticles.

        Arguments:
            fieldlist (list): List of string titles for the fields. Order is
                important; it must match the order for future particle appends
                that are made.
        """
        raise NotImplementedError
        # self._injectedparticles_fields = fieldlist
        # self._injectedparticles_data = warp.AppendableArray(
        #     typecode='d', unitshape=[len(fieldlist)])

    def record_injectedparticles(self, species, E_total=None, w=None,
                                 n=None):
        """Handles transforming raw particle information to the information
        used to record particles as a function of time. Also handles parallel
        sum and appending to the data array the current amount of injection.

        Note:
            Assumes the fixed form of fields given in Injector().  Doesn't
            check since this is called many times.
            Since a parallelsum is performed, call this with only the species
            argument if no particles are being added by this processor.

        Arguments:
            species: [TODO - adapt to WarpX]
            E_total (np.ndarray): Array of length npart with E_total values.
            w (np.ndarray): Array of length ``npart`` with variable weights.
            n (int): Number of macroparticles, _only_ needed if overriding the
                length of E_total. This is useful mostly in the case that
                E_total is already summed over particles, in which case a
                single number can be passed for it rather than an array.
        """
        raise NotImplementedError
        # data = np.zeros(7)

        # data[0] = warp.top.it*warp.top.dt
        # data[1] = warp.top.it
        # data[2] = species.jslist[0]
        # data[3] = self.getvoltage_e()

        # if E_total is not None:
        #     if n is not None:
        #         data[4] = n
        #     else:
        #         data[4] = len(E_total)

        #     if w is not None:
        #         data[5] = species.sq * species.sw * np.sum(w)
        #     else:
        #         data[5] = species.sq * species.sw * data[4]

        #     data[6] = np.sum(E_total)

        # self.append_injectedparticles(data)

    def append_injectedparticles(self, data):
        """Append one or more lines of injected particles data.

        Arguments:
            data (np.ndarray): Array of shape (m) or (n, m) where m is the
                number of fields and n is the number of rows of data to append.
        """
        self._injectedparticles_data.append(data)

    def get_injectedparticles(self, clear=False):
        """Retrieve a copy of injectedparticles data.

        Arguments:
            clear (bool): If True, clear the particle data rows entered (field
                names are still initialized as before). Default False.

        Returns:
            injectedparticles_dict (collections.OrderedDict): Keys are the
                originally passed field strings for lost particles. Values are
                an (n)-shape numpy array for each field.
        """
        raise NotImplementedError


class FixedNumberInjector(Injector):

    """Inject n particles every t timesteps."""

    def __init__(self, emitter, species, npart,
                 injectfreq=None, injectoffset=1,
                 weight=0., rseed=None,
                 name=None, unique_particles=True):
        """Sets up user-specified injection with fixed timestep and weights.

        Arguments:
            emitter (:class:`mewarpx.emission.Emitter`): Emitter object that
                will specify positions and velocities of particles to inject.
            species (picmi.Species): Premade species to inject particles of.
            npart (int): Number of particles to inject total
            injectfreq (int): Number of steps to wait for next injection.
                Default infinity.
            injectoffset (int): First timestep to inject. Default 1 (the
                first possible timestep in WarpX).
            weight (float): Macroparticle weight to be introduced.
            rseed (int): If specified, all injection should be repeatable using
                this rseed. At present each set of injected particles will have
                the same initial position and velocities as the previous set.
            name (str): Injector name for diagnostics. Constructed from
                speciesname if not given.
            unique_particles (bool): Whether WarpX will keep all particles
                given it from every processor (True) or keep only a fraction of
                particles based on processor count (False).
        """
        # Save class parameters
        self.emitter = emitter
        self.species = species
        self.npart_total = npart
        self.injectfreq = injectfreq
        if self.injectfreq is None:
            self.injectfreq = np.inf
        self.injectoffset = injectoffset
        self.weight = weight
        self.rseed = rseed
        self.name = name
        if self.name is None:
            self.name = "fixed_injector_" + self.species.name
        self.unique_particles = unique_particles

        print(f"Fixed injection of {self.npart_total} particles, "
              f"weight {self.weight}, every {self.injectfreq}"
              f"timesteps.")
        callbacks.installparticleinjection(self.inject_particles)

        # add E_total PID to this species
        self.species.add_pid("E_total")

    def inject_particles(self):
        """Perform the actual injection!"""
        effective_it = mwxrun.get_it() - self.injectoffset
        if effective_it >= 0 and effective_it % self.injectfreq == 0:

            # Adjust npart for processor number if needed
            npart = self.compute_npart(
                npart_total=self.npart_total,
                unique_particles=self.unique_particles
            )

            # TODO randomdt and velhalfstep are False simply because they're
            # not supported at present
            particles_dict = self.emitter.get_newparticles(
                npart=npart, w=self.weight,
                q=self.species.sq, m=self.species.sm,
                rseed=self.rseed,
                randomdt=False, velhalfstep=False
            )

            print(f"Inject {len(particles_dict['x'])} particles")

            # Note some parts of WarpX call the variables ux and some parts vx,
            # and they're referred to as momenta. But I don't see anywhere
            # they're actually used as momenta including the particle mass -
            # the actual update is in Source/Particles/Pusher/UpdatePosition.H
            self.species.add_particles(
                x=particles_dict['x'],
                y=particles_dict['y'],
                z=particles_dict['z'],
                ux=particles_dict['vx'],
                uy=particles_dict['vy'],
                uz=particles_dict['vz'],
                w=particles_dict['w'],
                unique_particles=self.unique_particles,
                E_total=particles_dict['E_total']
            )

            if self.injector_diag is not None:
                self.record_injectedparticles(
                    self.species,
                    particles_dict['E_total'],
                    particles_dict.get('w', None)
                )


class ThermionicInjector(Injector):

    """Performs standard every-timestep injection from a thermionic cathode."""

    def __init__(self, emitter, species, npart_per_cellstep, T,
                 WF, A=constants.A0*1e4,
                 allow_poisson=False, wfac=1.0,
                 name=None, profile_decorator=None,
                 unique_particles=True):
        """Sets up user-specified injection for warpX.

        Arguments:
            emitter (:class:`mewarpx.emission.Emitter`): Emitter object that
                will specify positions and velocities of particles to inject.
            species (mepicmi.Species): A premade species. Note only electrons
                will actually give physically meaningful weight calculations.
            npart_per_cellstep (int): Number of macroparticles to inject per
                cell on the cathode surface per timestep
            T (float): Cathode temperature (K)
            WF (float): Cathode work function (eV)
            A (float): Coefficient of emission in Amp/m^2/K^2. Default is
                the theoretical max, approximately 1.2e6.
            allow_poisson (bool): If True and < npart_per_cellstep electrons
                would be injected per cell, inject whole electrons with a
                Poisson distribution. If False, inject fractions of electrons.
                Default False.
            wfac (float): Constant factor applied to variable particle
                weights, which changes the actual injection weight from the
                physically calculated quantity. Currently used only for
                testing, or for e.g.  artificially lowering weight of trace
                particles.
            name (str or None): Injector name for diagnostics. Constructed from
                speciesname if not given.
            profile_decorator (decorator): A decorator used to profile the
                injection methods and related functions.
            unique_particles (bool): Whether WarpX will keep all particles
                given it from every processor (True) or keep only a fraction of
                particles based on processor count (False). Default True.
        """
        # Save class parameters
        self.emitter = emitter
        self.species = species
        self.T = T
        self.WF = WF
        self.A = A
        self.wfac = wfac

        if profile_decorator is not None:
            self.inject_particles = profile_decorator(self.inject_particles)
            self.record_injectedparticles = (
                profile_decorator(self.record_injectedparticles)
            )

        self.name = name
        if self.name is None:
            self.name = "thermionic_injector_" + self.species.name
        self.unique_particles = unique_particles

        area = self.emitter.area
        dt = mwxrun.get_dt()
        if (area is None) or (area <= 0.0) or (dt <= 0.0):
            raise ValueError(f"area {area} or dt {dt}"
                             f" is invalid for injection.")

        # Determine weight and injection numbers
        electrons_per_step = (mwxutil.J_RD(self.T, self.WF, self.A)
                              * area * dt / picmi.constants.q_e)
        print(f"Setting up thermionic paticle injection. Area {area:.3g} m^2, "
              f"dt {dt:.3e} s, J {mwxutil.J_RD(self.T, self.WF, self.A):.3g} "
              "A/m^2.")
        print(f"Emission current corresponds to injection of "
              f"{electrons_per_step:.2e} electrons per timestep")
        max_injections = int(round(npart_per_cellstep *
                                   self.emitter.cell_count))

        # If it was requested to inject more particles than we have electrons,
        # we instead inject electrons with a poisson distribution if allowed.
        if electrons_per_step < max_injections and allow_poisson:
            self.ptcl_per_step = electrons_per_step
            self.weight = self.wfac
            self.poisson = True
            print("Using stochastic injection of electrons with "
                  "Poisson sampling")
        else:
            self.ptcl_per_step = max_injections
            self.weight = self.wfac * electrons_per_step / self.ptcl_per_step
            self.poisson = False
            print(f"Using deterministic injection of {self.ptcl_per_step} "
                  f"particles per step, each with weight {self.weight}")
        callbacks.installparticleinjection(self.inject_particles)

        # add E_total PID to this species
        self.species.add_pid("E_total")

    def inject_particles(self):
        """Perform the actual injection!"""
        if self.poisson:
            num_injections = np.random.poisson(self.ptcl_per_step)
        else:
            num_injections = self.ptcl_per_step

        # Adjust npart for processor number if needed
        npart = self.compute_npart(
            npart_total=num_injections,
            unique_particles=self.unique_particles
        )

        # TODO randomdt and velhalfstep are False simply because they're
        # not supported at present
        particles_dict = self.emitter.get_newparticles(
            npart=npart, w=self.weight, q=self.species.sq, m=self.species.sm,
            randomdt=False, velhalfstep=False
        )

        # Note some parts of WarpX call the variables ux and some parts vx,
        # and they're referred to as momenta. But I don't see anywhere
        # they're actually used as momenta including the particle mass -
        # the actual update is in Source/Particles/Pusher/UpdatePosition.H
        self.species.add_particles(
            x=particles_dict['x'],
            y=particles_dict['y'],
            z=particles_dict['z'],
            ux=particles_dict['vx'],
            uy=particles_dict['vy'],
            uz=particles_dict['vz'],
            w=particles_dict['w'],
            unique_particles=self.unique_particles,
            E_total=particles_dict['E_total']
        )

        if self.injector_diag is not None:
            self.record_injectedparticles(
                self.species,
                particles_dict['E_total'],
                particles_dict.get('w', None)
            )


class BaseEmitter(object):

    """Parent class of both Emitter (which handles injection from a surface or
    other area) and VolumeEmitter (which handles injection throughout a
    volume).
    All BaseEmitter objects are expected to contain:
        - ``get_newparticles()`` returns coordinates, velocities, and KE in a
          dict - implemented here
        - ``_get_xv_coords()`` implements the subclass-specific particle
          injection logic
        - ``getvoltage()`` calculates the potential energy for particle
          energies.
        - ``getvoltage_e()`` calculates the potential energy for particle
          energies including the work function.
        - ``geoms`` is a property containing a list of simulation geometries
          supported by the Emitter, as strings
    """
    # Stores a list of functions that are used to adjust variable particle
    # weights.
    _wfnlist = None
    # Overridden only for surface injectors that use Schottky corrections
    use_Schottky = False
    # Needs to be overridden to specify acceptable geometries
    geoms = []

    def __init__(self):
        """Check geometry and any other universal initialization.
        """
        self.solver_geom = self.check_geom()

        # # Use to get E and phi as needed.
        # self.particle_helper = ParticleValHelper()

    def check_geom(self):
        """Return the current solver geometry, or throw an error if it is
        unsupported by the Emitter.
        """
        geom = mwxrun.geom_str

        if geom not in self.geoms:
            raise ValueError(
                f"{geom} geometry not supported by this Emitter")

        return geom

    def getvoltage(self):
        """This should return the potential energy at the injection site for
        fully accurate energetics.
        """
        raise NotImplementedError

    def getvoltage_e(self):
        """This should return the potential energy, including work function,
        at the injection site for fully accurate energetics.
        """
        raise NotImplementedError

    @staticmethod
    def _gen_particle_dict(x, y, z, vx, vy, vz, w, **kwargs):
        """Change standard arrays into format expected by an injector.
        The transfer to an injector uses a dict so that optional
        arguments can be passed, or additional arguments added.
        Arguments:
            x (np.ndarray): n-shape position array
            y (np.ndarray): n-shape position array
            z (np.ndarray): n-shape position array
            vx (np.ndarray): n-shape velocity array
            vy (np.ndarray): n-shape velocity array
            vz (np.ndarray): n-shape velocity array
            w (float or np.ndarray): Particle weight, either constant or
                per-particle.
            kwargs (np.ndarray): These are simply copied into the dictionary
        """
        particle_dict = {
            'x': x, 'y': y, 'z': z,
            'vx': vx, 'vy': vy, 'vz': vz,
            'w': np.ones_like(x) * w
        }
        particle_dict.update(kwargs)

        return particle_dict

    def _get_E_total(self, vx, vy, vz, q, m, w):
        """Calculate initial particle energies.

        Note:
            The conductor voltage V of the conductor the particle is ejected from
            must also be set for this object.

        Arguments:
            vx (np.ndarray): n-length array of velocity x-components
            vy (np.ndarray): n-length array of velocity y-components
            vz (np.ndarray): n-length array of velocity z-components
            q (float): Charge of the particles, usually species.sq.
            m (float): Mass of the particles, usually species.sm.
            w (np.ndarray): Variable particle weight, n-shape
        """
        V = self.getvoltage()

        E_total = w*(0.5*m*(vx**2 + vy**2 + vz**2) + q*V)

        return E_total

    def get_newparticles(self, npart, w, q, m, rseed=None,
                         randomdt=True, velhalfstep=True):
        """Return dict with coordinates, velocities, and KE

        Note:
            This function SHOULD (but doesn't in WarpX yet)  handle the random
            timestep advancement and the negative half-step velocity push. They
            can be turned off if desired.  No leapfrogging is done in the
            initial random advancement, which could be a (hopefully very minor)
            source of error.

        Arguments:
            npart (int): Total number of particles to inject
            w (float): Weight of the particles
            q (float): Charge of the particles, usually species.sq.
            m (float): Mass of the particles, usually species.sm. Equals
                electron mass if not otherwise specified.
            rseed (int): Random seed, if specified, can be used to provide
                reproducible results. Typically used for test / not production
                runs.
            randomdt (bool): If True, move each particle ahead a random delta t
                in [0, dt), advancing both position and velocity together.
                Default True.
            velhalfstep (bool): If True, push the velocities a negative
                half-step using the E-field. Aligns position and velocities
                correctly for the leapfrog algorithm.

        Returns:
            particle_dict (dict): Contains lists, each with length equal to the
            number of particles:

                - ``x``, ``y``, and ``z`` contain initial positions
                - ``vx``, ``vy``, and ``vz`` contain initial velocities
                - ``E_total`` contains initial energy of each particle, kinetic
                  & potential.
                - ``w`` contains particle weights.
        """
        if rseed is not None:
            nprstate = np.random.get_state()
            np.random.seed(rseed)
            rseedxv = np.random.randint(1000000000)
            rseedt = np.random.randint(1000000000)
        else:
            rseedxv = None
            rseedt = None

        x, y, z, vx, vy, vz = self._get_xv_coords(
            npart=npart, m=m, rseed=rseedxv
        )
        particle_dict = self._gen_particle_dict(
            x=x, y=y, z=z, vx=vx, vy=vy, vz=vz, w=w
        )

        if self._wfnlist is not None:
            for wfn in self._wfnlist:
                particle_dict['w'] = wfn(particle_dict)

        if self.use_Schottky:
            raise RuntimeError("Error: Schottky enhancement not currently implemented!")

        # if (self.use_Schottky and abs(q + e)/e < 1e-6
        #         and abs(m - m_e)/m_e < 1e-6):
        #     particle_dict['w'] = self.apply_Schottky_weights(particle_dict)

        particle_dict['E_total'] = self._get_E_total(
            vx=particle_dict['vx'],
            vy=particle_dict['vy'],
            vz=particle_dict['vz'],
            q=q, m=m, w=particle_dict['w']
        )

        # After E_total has been computed, we advance particles as needed.
        if randomdt:
            self.particle_helper.advance_random_deltat(
                particle_dict['x'], particle_dict['y'], particle_dict['z'],
                particle_dict['vx'], particle_dict['vy'], particle_dict['vz'],
                q=q, m=m, rseed=rseedt
            )

        if velhalfstep:
            self.particle_helper.push_v_minus_halfstep(
                particle_dict['x'], particle_dict['y'], particle_dict['z'],
                particle_dict['vx'], particle_dict['vy'], particle_dict['vz'],
                q=q, m=m
            )

        if rseed is not None:
            np.random.set_state(nprstate)

        return particle_dict

    def _update_params(self):
        """Update local parameters if needed based on WarpX settings.
        By default does nothing, but subclasses can implement it to update
        parameters before new particle coordinates are generated.
        """
        pass

    def _get_xv_coords(self, npart, m, rseed):
        """Per-subclass implementation of generating new particle data.
        See :func:`mewarpx.emission.BaseEmitter.get_newparticles` for details on
        arguments.

        Returns:
            x, y, z, vx, vy, vz (np.array): Each must be a 1D numpy array.
        """
        raise NotImplementedError(
            "BaseEmitter subclasses must implement _get_xv_coords")

    def add_wfn(self, wfn):
        """Add a variable weight function to the emitter.

        Arguments:
            wfn (function): This must take in a particle dictionary with
                positions, velocities, and existing weights, and return a new
                array of particle weights.
        """
        if self._wfnlist is None:
            self._wfnlist = []

        self._wfnlist.append(wfn)


class Emitter(BaseEmitter):

    """Parent class for emission from a surface.

    All Emitter objects are expected to contain:
        - ``area`` is a property containing the area in m^2
        - ``cell_count`` is a property containing the number of mesh cells
          spanned by the Emitter
        - ``geoms`` is a property containing a list of simulation geometries
          supported by the Emitter
        - ``_get_xv_coords()`` implements the subclass-specific particle
          injection logic
        - ``get_normals()`` returns the normals for a set of particle
          coordinates.
    """
    area = None
    cell_count = None
    geoms = []

    def __init__(self, T, conductor=None, use_Schottky=True,
                 emission_type='thermionic'):
        """Default initialization for all Emitter objects.

        Arguments:
            T (float): Emitter temperature in Kelvin. Determines particle
                velocity distribution.
            conductor (assemblies.Assembly): Conductor the emitter is attached
                to, used for recording initial voltages and energies. If None,
                V_e is set to 0. Since there's no current use case for this, a
                warning is printed.
            use_Schottky (bool): Flag specifying whether or not to augment the
                Emitter's emission current via field-dependent particle
                weights.  Defaults to True.
            emission_type (str): Distribution function type used to sample
                velocities of the emitted particles. Must be defined in
                :func:`mewarpx.util.get_velocities`. Defaults to 'thermionic'.
        """
        super(Emitter, self).__init__()
        self.T = T

        self.conductor = conductor
        if self.conductor is not None:
            if self.conductor.WF <= 0.:
                raise ValueError("Conductor WF must be set for emitters.")
        else:
            warnings.warn("No conductor set for emitter. Power will not be "
                          "correct.")

        self.use_Schottky = use_Schottky
        self.emission_type = emission_type

        if self.use_Schottky:
            raise NotImplementedError("Schottky enhancement not implemented"
                                      " in WarpX yet.")

    def getvoltage(self):
        if self.conductor is None:
            return 0.

        return self.conductor.getvoltage()

    def getvoltage_e(self):
        """Electrical voltage includes WF, eg the Fermi level voltage."""
        if self.conductor is None:
            return 0.

        return self.conductor.getvoltage() + self.conductor.WF

    def apply_Schottky_weights(self, particle_dict):
        """Variable weight function for field-enhanced Schottky emission.

        Notes:
            This function requires the "T" attribute of the Emitter to be set
            for specifying the emitter temperature in Kelvin, and also requires
            the Emitter subclass to implement a get_normals() function.

        Arguments:
            particle_dict (dict): Particle dictionary returned by
                get_newparticles() during standard thermionic injection. Any
                existing variable weights in the dictionary should be
                normalized to match the total zero-field saturation current of
                the Emitter without any Schottky enhancement.

        Returns:
            new_weights (np.ndarray): 1D array of updated particle weights.
        """
        raise NotImplementedError

    def get_normals(self, x, y, z):
        """Calculate local surface normal at specified coordinates.

        Arguments:
            x (np.ndarray): x-coordinates of emitted particles (in meters).
            y (np.ndarray): y-coordinates of emitted particles (in meters).
            z (np.ndarray): z-coordinates of emitted particles (in meters).

        Returns:
            normals (np.ndarray): nx3 array containing the outward surface
                normal vector at each particle location.
        """
        raise NotImplementedError('Normal calculations must be implemented by'
                                  + ' Emitter sub-classes.')


class ZPlaneEmitter(Emitter):
    """This is the standard injection for a planar cathode."""

    geoms = ['Z', 'XZ', 'XYZ']

    def __init__(self, conductor, T=None, xmin=None, xmax=None,
                 ymin=None, ymax=None, transverse_fac=1.0, **kwargs):
        """Initialize an emitter for a planar cathode.

        Arguments:
            conductor (:class:`mewarpx.assemblies.Assembly`): Conductor object, used to obtain work
                function. Can later grab other variables from this conductor.
            T (float): Temperature in Kelvin for the emitter; determines
                velocities.
            z (float): Position of the emitter along the z axis. If None, grab
                zcent from the conductor.
            zsign (float): -1 for the cathode, +1 for the anode, or None to
                grab automatically from the conductor.
            xmin (float): Minimum position of the rectangular emitter along x.
                Default mwxrun.xmin.
            xmax (float): Maximum position of the rectangular emitter along x.
                Default mwxrun.xmax.
            ymin (float): Minimum position of the rectangular emitter along y.
                Default mwxrun.ymin.
            ymax (float): Maximum position of the rectangular emitter along y.
                Default mwxrun.ymax.
            transverse_fac (float): Scale the transverse energy distribution by
                this factor. Default 1. See
                :func:`mewarpx.util.get_velocities` for details.
            kwargs (dict): Any other keyword arguments supported by the parent
                Emitter constructor (such as "use_Schottky" or
                "emission_type").
        """
        # Default initialization
        if T is None:
            T = conductor.T
        super(ZPlaneEmitter, self).__init__(T=T, conductor=conductor, **kwargs)

        self.z = conductor.z
        self.zsign = conductor.zsign
        self.transverse_fac = transverse_fac

        # Determine bounds
        # Will be 4 element array [xmin, xmax, ymin, ymax]
        self.bounds = []
        for coord, default in [(xmin, mwxrun.xmin),
                               (xmax, mwxrun.xmax),
                               (ymin, mwxrun.ymin),
                               (ymax, mwxrun.ymax)]:
            self.bounds.append(coord if coord is not None else default)

        # Compute area
        x_range = self.bounds[1] - self.bounds[0]
        y_range = self.bounds[3] - self.bounds[2]
        if self.solver_geom == 'Z':
            print("x/y span is 1m for purposes of charge injection")
            x_range = 1.
            y_range = 1.
        if self.solver_geom == 'XZ':
            print("y span is 1m for purposes of charge injection")
            y_range = 1.
        self.area = x_range * y_range

        # Compute cell count
        if self.solver_geom == 'Z':
            self.cell_count = 1
        elif self.solver_geom == 'XZ':
            self.cell_count = self.area / mwxrun.dx
        else:
            self.cell_count = self.area / (mwxrun.dx * mwxrun.dy)

    def _get_xv_coords(self, npart, m, rseed):
        """Get particle coordinates given particle number.
        See :func:`mewarpx.emission.BaseEmitter.get_newparticles` for details.
        """
        if rseed is not None:
            nprstate = np.random.get_state()
            np.random.seed(rseed)
            rseedv = np.random.randint(1000000000)
            rseedx = np.random.randint(1000000000)
        else:
            rseedv = None
            rseedx = None

        vx, vy, vz = mwxutil.get_velocities(
            npart, self.T, m=m, transverse_fac=self.transverse_fac,
            emission_type=self.emission_type, rseed=rseedv)
        x, y, z = mwxutil.get_positions(
            npart, xmin=self.bounds[0], xmax=self.bounds[1],
            ymin=self.bounds[2], ymax=self.bounds[3], z=self.z,
            rseed=rseedx)

        # Flip z velocities for anode emission. This appears to be faster than
        # an if statement for 10000 or fewer particles.
        vz = -self.zsign * vz

        if rseed is not None:
            np.random.set_state(nprstate)

        return x, y, z, vx, vy, vz

    def get_normals(self, x, y, z):
        """Calculate local surface normal at specified coordinates.

        Arguments:
            x (np.ndarray): x-coordinates of emitted particles (in meters).
            y (np.ndarray): y-coordinates of emitted particles (in meters).
            z (np.ndarray): z-coordinates of emitted particles (in meters).

        Returns:
            normals (np.ndarray): nx3 array containing the outward surface
                normal vector at each particle location.
        """
        normals = np.zeros((len(x), 3))
        normals[:, 2] = -self.zsign
        return normals


class ArbitraryEmitter2D(Emitter):

    """ ArbitraryEmitter2D class takes in a conductor, calculates an approximate
    surface that encloses the conductor and then sets up the appropriate
    emitting surfaces, given a number of particles to emit.
    """
    geoms = ['XZ']

    def __init__(self, conductor, T, res_fac=5., transverse_fac=1.0, **kwargs):
        """Construct the emitter based on conductor object and temperature.

        Arguments:
            conductor (mewarpx.assemblies object): Conductor to emit from.
            T (float): Temperature in Kelvin.
            res_fac (float): Level of resolution beyond the grid resolution to
                use for calculating shape contours.
            transverse_fac (float): Scale the transverse energy distribution by
                this factor. Default 1. See
                :func:`mewarpx.util.get_velocities` for details.
            kwargs (dict): Any other keyword arguments supported by the parent
                Emitter constructor (such as "use_Schottky" or
                "emission_type").
        """
        # Default Emitter initialization.
        super(ArbitraryEmitter2D, self).__init__(T=T, conductor=conductor, **kwargs)
        # Save input parameters
        self.res_fac = res_fac
        self.transverse_fac = transverse_fac

        # Generate grid enclosed in bounding box
        self.dx = mwxrun.dx/res_fac
        self.dy = 1.
        self.dz = mwxrun.dz/res_fac

        self.dA = np.sqrt(self.dx*self.dz)

        # A small delta is added to the maxima here; this ensures the last point
        # is included. Without it, floating point errors determine whether or
        # not the last point is included.
        self.xvec = np.arange(
            mwxrun.xmin, mwxrun.xmax + self.dx/1000., self.dx)
        self.yvec = [0.]
        self.zvec = np.arange(
            mwxrun.zmin, mwxrun.zmax + self.dz/1000., self.dz)

        [X, Y, Z] = np.squeeze(np.meshgrid(self.xvec, self.yvec, self.zvec,
                                           indexing='xy'))
        oshape = X.shape
        X = X.flatten()
        Y = Y.flatten()
        Z = Z.flatten()

        inside = np.reshape(
            self.conductor.isinside(X, Y, Z, aura=self.dA/5.),
            oshape)

        self.contours = np.squeeze(skimage.measure.find_contours(
            inside, 0.5))
        self.contours[:, 0] = np.interp(self.contours[:, 0],
                                        np.arange(self.xvec.size),
                                        self.xvec)
        self.contours[:, 1] = np.interp(self.contours[:, 1],
                                        np.arange(self.zvec.size),
                                        self.zvec)

        self.centers = np.array(
            [(self.contours[1:, 0] + self.contours[:-1, 0])/2.,
             (self.contours[1:, 1] + self.contours[:-1, 1])/2.]).T
        self.dvec = np.array(
            [self.contours[1:, 0] - self.contours[:-1, 0],
             self.contours[1:, 1] - self.contours[:-1, 1]]).T

        # Calculate the distance of each segment & sum to calculate the area
        self.distances = np.sqrt(self.dvec[:, 0]**2 + self.dvec[:, 1]**2)
        self.area = sum(self.distances)
        self.cell_count = self.area / min(mwxrun.dx, mwxrun.dz)
        self.CDF = np.cumsum(self.distances)/self.area

        # Calculate Normal Vector by taking cross product with y-hat
        ndvec = self.dvec/np.tile(self.distances, (2, 1)).T
        marching_normal = np.zeros(self.dvec.shape)
        marching_normal[:, 0] = -ndvec[:, 1]
        marching_normal[:, 1] = ndvec[:, 0]

        # Check to make sure normal plus center is outside of conductor
        partdist = self.dA*float(self.res_fac)/2.

        pos = self.centers + marching_normal*partdist
        px = pos[:, 0]
        py = np.zeros_like(px)
        pz = pos[:, 1]

        nhat = self.conductor.calculatenormal(px, py, pz)
        self.normal = nhat[[0, 2], :].T

    def _get_xv_coords(self, npart, m, rseed):
        """Get particle coordinates given particle number.
        See :func:`mewarpx.emitter.get_newparticles` for details.
        """
        if rseed is not None:
            nprstate = np.random.get_state()
            np.random.seed(rseed)
            # rseedv is passed to get velocities. The basic rseed here is used
            # for positions, below.
            rseedv = np.random.randint(1000000000)
        else:
            rseedv = None

        # Draw Random Numbers to determine which face to emit from
        self.contour_idx = np.searchsorted(self.CDF, np.random.rand(npart))

        vels = np.column_stack(mwxutil.get_velocities(num_samples=npart, T=self.T, m=m,
                                                        rseed=rseedv,
                                                        transverse_fac=self.transverse_fac,
                                                        emission_type=self.emission_type))

        # Rotate velocities based on angle of normal
        newvels = self.convert_vel_zhat_nhat(vels, self.normal[self.contour_idx])
        vx = np.asarray(newvels[:, 0], order="C")
        vy = np.asarray(newvels[:, 1], order="C")
        vz = np.asarray(newvels[:, 2], order="C")

        # Now get positions
        pos1 = self.contours[self.contour_idx, :]
        positions = (pos1 +
                     (np.tile(np.random.rand(npart), (2, 1)).T
                      * self.dvec[self.contour_idx, :]))

        x = np.asarray(positions[:, 0], order="C")
        y = np.asarray(0., order="C")
        z = np.asarray(positions[:, 1], order="C")

        if rseed is not None:
            np.random.set_state(nprstate)

        return x, y, z, vx, vy, vz

    @staticmethod
    # Synthetic tests showed 18 ms to 660us change from using np.dot +
    # numba compilation. Without these changes, this function was taking 2-4% of
    # some run times so the improvement is warranted.
    @numba.jit(nopython=True)
    def convert_vel_zhat_nhat(vels, nhat):
        """Create a rotation matrix for Zhat to Nhat"""
        Zhat = np.array([0., 1.])

        newvels = np.zeros(vels.shape)

        for ii in range(vels.shape[0]):
            Cvec = Zhat - nhat[ii, :]
            Cvec2 = np.dot(Cvec, Cvec)

            theta = np.arccos(1. - Cvec2/2.)

            # Check to see if normal is pointing toward -xhat
            # Resolves angle ambiguity in law of cosines
            if nhat[ii, 0] < 0.:
                theta = -theta

            # Rotate in XZ plane, keeping Y the same
            R = np.array([[np.cos(theta), 0., np.sin(theta)],
                          [0., 1., 0.],
                          [-np.sin(theta), 0., np.cos(theta)]])

            newvels[ii, :] = np.dot(R, vels[ii, :])

        return newvels

    def get_normals(self, x, y, z):
        """Calculate local surface normal at specified coordinates.
        Arguments:
            x (np.ndarray): x-coordinates of emitted particles (in meters).
            y (np.ndarray): y-coordinates of emitted particles (in meters).
            z (np.ndarray): z-coordinates of emitted particles (in meters).

        Returns:
            normals (np.ndarray): nx3 array containing the outward surface
                normal vector at each particle location.
        """
        # Since we've already pre-computed all the normals and already picked
        # the right ones during the call to _get_xv_coords(), we can ignore the
        # coordinate arguments here entirely and use the recently saved
        # "contour_idx" values for indexing the pre-tabulated normals. To
        # prevent this from being abused, we'll first check that the length of
        # the coordinate lists matches that of the contour_idx list.
        if len(x) != len(self.contour_idx):
            raise ValueError('Length of particle coordinate list does not match'
                             + ' the most recent number of emitted particles!')

        normals = np.zeros((len(x), 3))
        normals[:, 0] = self.normal[self.contour_idx, 0]
        normals[:, 2] = self.normal[self.contour_idx, 1]
        return normals
