"""
Module for various types of particle emission in WarpX.
"""
# import collections
import logging
import warnings

import numba
import numpy as np
import collections

import skimage.measure
import matplotlib.colors as colors
import matplotlib.pyplot as plt

from pywarpx import callbacks, _libwarpx, picmi

import mewarpx.utils_store.util as mwxutil
from mewarpx.mwxrun import mwxrun
import mewarpx.utils_store.mwxconstants as constants
from mewarpx.utils_store import appendablearray, parallel_util
from mewarpx.mespecies import Species

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
    fields = ['t', 'step', 'species_id', 'V_e', 'n', 'q', 'E_total']

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
        self._injectedparticles_fields = fieldlist
        self._injectedparticles_data = appendablearray.AppendableArray(
            typecode='d', unitshape=[len(fieldlist)])

    def record_injectedparticles(self, species, w, E_total=None,
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
            species (:class:`mewarpx.mespecies.Species`): Species of particle
            w (np.ndarray or float): Array of length npart with particle weights
            E_total (np.ndarray or float): Array of length npart with E_total
                values.
            n (int): Number of macroparticles, _only_ needed if overriding the
                length of E_total. This is useful mostly in the case that
                E_total is already summed over particles, in which case a
                single number can be passed for it rather than an array.
        """
        if n is not None and np.size(w) != 1:
            raise RuntimeError("Cannot pass array for w and specify n")

        if n is None and np.size(w) == 1:
            raise RuntimeError("Cannot pass single value for w and not specify n")

        data = np.zeros(7)

        # time for current step
        data[0] = mwxrun.get_it() * mwxrun.get_dt()
        # current step
        data[1] = mwxrun.get_it()
        # species ID
        data[2] = species.species_number
        # voltage of emitter
        data[3] = self.getvoltage_e()
        # number of macroparticles
        data[4] = n if np.size(w) == 1 else np.size(w)
        # total charge emitted
        data[5] = species.sq * np.sum(w)

        if E_total is not None:
            data[6] = np.sum(E_total)

        self.append_injectedparticles(data)

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
        lpdata = self._injectedparticles_data.data()

        # Sum all except t/step/species_id/V_e from all processors
        lpdata[:,4:] = parallel_util.parallelsum(np.array(lpdata[:,4:]))

        lpdict = collections.OrderedDict(
            [(fieldname, np.array(lpdata[:, ii], copy=True))
             for ii, fieldname in enumerate(self._injectedparticles_fields)])

        if clear:
            self._injectedparticles_data.cleardata()

        return lpdict


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

        logger.info(
            f"Fixed injection of {self.npart_total} particles, "
            f"weight {self.weight}, every {self.injectfreq}"
            f"timesteps."
        )
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

            logger.info(f"Inject {len(particles_dict['x'])} particles")

            # Note some parts of WarpX call the variables ux and some parts vx,
            # and they're referred to as momenta. But I don't see anywhere
            # they're actually used as momenta including the particle mass -
            # the actual update is in Source/Particles/Pusher/UpdatePosition.H
            _libwarpx.add_particles(
                self.species.name,
                x=particles_dict['x'],
                y=particles_dict['y'],
                z=particles_dict['z'],
                ux=particles_dict['vx'],
                uy=particles_dict['vy'],
                uz=particles_dict['vz'],
                w=particles_dict['w'],
                E_total=particles_dict['E_total'],
                unique_particles=self.unique_particles
            )

            if self.injector_diag is not None:
                self.record_injectedparticles(
                    species=self.species,
                    w=particles_dict['w'],
                    E_total=particles_dict['E_total'],
                )


class ThermionicInjector(Injector):

    """Performs standard every-timestep injection from a thermionic cathode."""

    def __init__(self, emitter, species, npart_per_cellstep, T=None,
                 WF=None, A=constants.A0*1e4, use_Schottky=True,
                 allow_poisson=False, wfac=1.0,
                 name=None, profile_decorator=None,
                 unique_particles=True):
        """Sets up user-specified injection for warpX.

        Arguments:
            emitter (:class:`mewarpx.emission.Emitter`): Emitter object that
                will specify positions and velocities of particles to inject.
            species (mewarpx.mespecies.Species): A premade species. Note only
                electrons will actually give physically meaningful weight
                calculations.
            npart_per_cellstep (int): Number of macroparticles to inject per
                cell on the cathode surface per timestep
            T (float): Cathode temperature (K). Uses emitter T if not specified.
            WF (float): Cathode work function (eV). Uses WF of the conductor
                associated with the emitter if not specified.
            A (float): Coefficient of emission in Amp/m^2/K^2. Default is
                the theoretical max, approximately 1.2e6.
            use_Schottky (bool): Flag specifying whether or not to augment the
                emission current via field-dependent particle weights.
                Defaults to True.
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
        # sanity check species
        if species.particle_type != 'electron':
            raise AttributeError(
                "Thermionic emission is only applicable with electrons as the "
                f"injection species, but species type {species.particle_type} "
                "was given."
            )

        # Save class parameters
        self.emitter = emitter
        self.species = species
        self.T = T
        self.WF = WF
        self.A = A
        self.use_Schottky = use_Schottky
        self.wfac = wfac

        # Get values from the emitter and its conductor if not specified
        if self.T is None:
            self.T = self.emitter.T
        if self.WF is None:
            self.WF = self.emitter.conductor.WF

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
        logger.info(
            f"Setting up thermionic paticle injection. Area {area:.3g} m^2, "
            f"dt {dt:.3e} s, J {mwxutil.J_RD(self.T, self.WF, self.A):.3g} "
            "A/m^2."
        )
        logger.info(
            "Emission current corresponds to injection of "
            f"{electrons_per_step:.2e} electrons per timestep"
        )
        max_injections = int(round(npart_per_cellstep *
                                   self.emitter.cell_count))

        # If it was requested to inject more particles than we have electrons,
        # we instead inject electrons with a poisson distribution if allowed.
        if electrons_per_step < max_injections and allow_poisson:
            self.ptcl_per_step = electrons_per_step
            self.weight = self.wfac
            self.poisson = True
            logger.info(
                "Using stochastic injection of electrons with "
                "Poisson sampling"
            )
        else:
            self.ptcl_per_step = max_injections
            self.weight = self.wfac * electrons_per_step / self.ptcl_per_step
            self.poisson = False
            logger.info(
                f"Using deterministic injection of {self.ptcl_per_step} "
                f"particles per step, each with weight {self.weight}"
            )

        # create new species that will be used to properly distribute new
        # particles and retrieve the electric field at their injection sites in
        # order to calculate Schottky enhancement
        if self.use_Schottky:
            self.injection_species = Species(
                particle_type='electron', name=self.species.name+'_injection'
            )
        else:
            self.injection_species = self.species

        callbacks.installparticleinjection(self.inject_particles)

        # add E_total PID to this species
        self.species.add_pid("E_total")
        self.injection_species.add_pid("E_total")

        if self.use_Schottky:
            # add PIDs to hold the normal vector
            # TODO work out a better way to handle these PIDs since this is not
            # a great use of memory
            self.species.add_pid("norm_x")
            self.species.add_pid("norm_y")
            self.species.add_pid("norm_z")
            self.injection_species.add_pid("norm_x")
            self.injection_species.add_pid("norm_y")
            self.injection_species.add_pid("norm_z")

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

        extra_pids = {}
        extra_pids['E_total'] = particles_dict['E_total']
        extra_pids['w'] = particles_dict['w']
        if self.use_Schottky:
            # Determine the local surface normal for each particle
            normal_vectors = self.emitter.get_normals(
                particles_dict['x'], particles_dict['y'], particles_dict['z']
            )
            extra_pids['norm_x'] = normal_vectors[:, 0]
            extra_pids['norm_y'] = normal_vectors[:, 1]
            extra_pids['norm_z'] = normal_vectors[:, 2]

        # Note some parts of WarpX call the variables ux and some parts vx,
        # and they're referred to as momenta. But I don't see anywhere
        # they're actually used as momenta including the particle mass -
        # the actual update is in Source/Particles/Pusher/UpdatePosition.H
        _libwarpx.add_particles(
            self.injection_species.name,
            x=particles_dict['x'],
            y=particles_dict['y'],
            z=particles_dict['z'],
            ux=particles_dict['vx'],
            uy=particles_dict['vy'],
            uz=particles_dict['vz'],
            unique_particles=self.unique_particles,
            **extra_pids
        )

        if self.use_Schottky:
            # Up-weight the particles by the local Schottky factor, calculated
            # as exp[sqrt(e / 4*pi*eps0) / (kT) * sqrt(max(-E, 0))]
            pre_fac = (
                np.sqrt(constants.e / (4.0 * np.pi * constants.epsilon_0))
                / (constants.kb_eV * self.emitter.T)
            )
            mwxrun.calc_Schottky_weight(
                self.injection_species.name, pre_fac
            )

            # get the total injected weight and energy
            total_weight = 0.
            total_energy = 0.
            npart = 0
            weight_arrays = _libwarpx.get_particle_arrays(
                self.injection_species.name, 'w', 0
            )
            ux_arrays = _libwarpx.get_particle_arrays(
                self.injection_species.name, 'ux', 0
            )
            uy_arrays = _libwarpx.get_particle_arrays(
                self.injection_species.name, 'uy', 0
            )
            uz_arrays = _libwarpx.get_particle_arrays(
                self.injection_species.name, 'uz', 0
            )
            for ii, w in enumerate(weight_arrays):
                npart += len(w)
                total_weight += np.sum(w)
                total_energy += np.sum(self.emitter._get_E_total(
                    ux_arrays[ii], uy_arrays[ii], uz_arrays[ii],
                    constants.e, constants.m_e, w
                ))

            # Move particles from temporary container to "real" container
            mwxrun.move_particles_between_species(
                self.injection_species.name, self.species.name
            )
        else:
            total_weight = np.sum(particles_dict['w'])
            total_energy = np.sum(particles_dict['E_total'])

        if self.injector_diag is not None:
            self.record_injectedparticles(
                species=self.species,
                w=total_weight,
                E_total=total_energy,
                n=npart
            )


class PlasmaInjector(Injector):

    """Inject particles at simulation start, or at regular timesteps, to
    seed a plasma. Can use any emitter object. The defining feature is that the
    2nd species positions and weights are copied from the first species, so the
    spatial distribution is always identical to start. Velocities are
    independent, however.
    """

    def __init__(self, emitter, species1, species2, npart, T_2=None,
                 plasma_density=None, ionization_frac=None,
                 P_neutral=None, T_neutral=None,
                 injectfreq=None, injectoffset=1,
                 rseed=None, name=None, unique_particles=True
                 ):
        """Initialize injection of a plasma with two species and given emitter.

        Arguments:
            emitter (:class:`mewarpx.emission.BaseEmitter`): BaseEmitter object
                that will specify positions and velocities of particles to
                inject.
            species1 (:class:`mewarpx.mespecies.Species`): First species, eg
                electron
            species2 (:class:`mewarpx.mespecies.Species`): Second species, eg ion
            npart (int): Number of macroparticles to inject total among all
                processors and species.
            T_2 (float): If specified, species2 will be injected at this
                temperature.
            plasma_density (float): Ion number density to inject. If using
                volumetric emitter, in m^(-3), if using surface emitter, in
                m^(-2)
            ionization_frac (float): Instead of plasma_density, use a specific
                ionization fraction of the neutral gas. Volumetric emitter
                only.
            P_neutral (float): If using ionization_frac only, the neutral gas
                density (*Torr*).
            T_neutral (float): If using ionization_frac only, the neutral gas
                temperature (K).
            injectfreq (int): Number of steps to wait for next injection.
                Default infinity.
            injectoffset (int): First timestep to inject. Default 1 (the first
                possible timestep in WarpX).
            rseed (int): If specified, all injection should be repeatable using
                this rseed. At present each set of injected particles will have
                the same initial position and velocities as the previous set.
            name (str or None): Injector name for diagnostics. Constructed from
                species names if not given.
            unique_particles (bool): Whether WarpX will keep all particles
                given it from every processor (True) or keep only a fraction of
                particles based on processor count (False). Default True.
        """
        # Save class parameters
        self.emitter = emitter
        self.npart_per_species = npart // 2
        self.species1 = species1
        self.species2 = species2
        self.T_2 = T_2
        if injectfreq is None:
            injectfreq = np.inf
        self.injectfreq = injectfreq
        self.injectoffset = injectoffset
        self.rseed = rseed

        self._calc_plasma_density(
            plasma_density=plasma_density,
            ionization_frac=ionization_frac,
            P_neutral=P_neutral,
            T_neutral=T_neutral,
        )

        self.name = name
        if self.name is None:
            self.name = (
                f"plasma_injector_{self.species1.name}_{self.species2.name}"
            )
        self.unique_particles = unique_particles

        logger.info(
            f"Plasma injection {self.name}: "
            f"{self.npart_per_species} particles each of {self.species1.name} "
            f"and {self.species2.name}, every {self.injectfreq} timesteps,"
        )

        # Surface emission
        if isinstance(self.emitter, Emitter):
            self.weight = (
                self.emitter.area * self.plasma_density / self.npart_per_species
            )
            warnings.warn(
                "Using a surface emitter with the PlasmaInjector has not been "
                "tested for accuracy."
            )
            logger.info(
                f"  full weight {self.weight:.4g}, surface density "
                f"{self.plasma_density:.4g} m^-2, area "
                f"{self.emitter.area:.4g} m^2."
            )
        # Volume emission
        else:
            self.weight = (
                self.emitter.volume * self.plasma_density
                / self.npart_per_species
            )
            logger.info(
                f"  full weight {self.weight:.4g}, volume density "
                f"{self.plasma_density:.4g} m^-3, volume "
                f"{self.emitter.volume:.4g} m^3."
            )
            debye_length = mwxutil.plasma_Debye_length(
                self.emitter.T, self.plasma_density)
            logger.info(
                f"  Corresponding plasma Debye length is {debye_length:.3e} m."
            )

        callbacks.installparticleinjection(self.inject_particles)

        # add E_total PID to the species involved
        self.species1.add_pid("E_total")
        self.species2.add_pid("E_total")

    def _calc_plasma_density(self, plasma_density, ionization_frac, P_neutral,
                             T_neutral):
        """Helper function to separate out part of initialization."""
        self.plasma_density = plasma_density
        if ionization_frac is not None:
            if self.plasma_density is not None:
                raise ValueError(
                    "Specify ionization_frac or plasma_density, not both.")
            if (
                (P_neutral is None) or (P_neutral <= 0) or
                (T_neutral is None) or (T_neutral <= 0)
            ):
                raise ValueError("Must specify positive neutral pressure and "
                                 "temperature to use ionization_frac.")
            if isinstance(self.emitter, Emitter):
                raise RuntimeError("Cannot use ionization_frac with a surface"
                                   " (area-based) Emitter.")
            n_neutral = mwxutil.ideal_gas_density(P_neutral, T_neutral)
            self.plasma_density = n_neutral * ionization_frac

        if (self.plasma_density is None) or (self.plasma_density <= 0):
            raise ValueError("Invalid plasma_density {}".format(
                self.plasma_density))

    def inject_particles(self):
        """Inject particles, same position & weight for each."""
        effective_it = mwxrun.get_it() - self.injectoffset
        if effective_it >= 0 and effective_it % self.injectfreq == 0:
            # Adjust npart for processor number if needed
            npart = self.compute_npart(
                npart_total=self.npart_per_species,
                unique_particles=self.unique_particles
            )

            # TODO randomdt and velhalfstep are False simply because they're
            # not supported at present
            particles1_dict = self.emitter.get_newparticles(
                npart=npart, w=self.weight, q=self.species1.sq,
                m=self.species1.sm, rseed=self.rseed,
                randomdt=False, velhalfstep=False
            )

            # if requested get particles for species2 at the specified
            # temperature
            if self.T_2 is not None:
                T_temp, self.emitter.T = self.emitter.T, self.T_2

            # TODO randomdt and velhalfstep are False simply because they're
            # not supported at present
            particles2_dict = self.emitter.get_newparticles(
                npart=npart, w=self.weight, q=self.species2.sq,
                m=self.species2.sm, rseed=self.rseed,
                randomdt=False, velhalfstep=False
            )
            if self.T_2 is not None:
                self.emitter.T = T_temp

            for key in ['x', 'y', 'z', 'w']:
                particles2_dict[key] = particles1_dict[key]

            logger.info(
                f"Inject {len(particles1_dict['x'])} particles each of "
                f"{self.species1.name} and {self.species2.name}."
            )
            _libwarpx.add_particles(
                self.species1.name,
                x=particles1_dict['x'],
                y=particles1_dict['y'],
                z=particles1_dict['z'],
                ux=particles1_dict['vx'],
                uy=particles1_dict['vy'],
                uz=particles1_dict['vz'],
                w=particles1_dict['w'],
                E_total=particles1_dict['E_total'],
                unique_particles=self.unique_particles
            )
            _libwarpx.add_particles(
                self.species2.name,
                x=particles2_dict['x'],
                y=particles2_dict['y'],
                z=particles2_dict['z'],
                ux=particles2_dict['vx'],
                uy=particles2_dict['vy'],
                uz=particles2_dict['vz'],
                w=particles2_dict['w'],
                E_total=particles2_dict['E_total'],
                unique_particles=self.unique_particles
            )

            if self.injector_diag is not None:
                self.record_injectedparticles(
                    species=self.species1,
                    w=particles1_dict['w'],
                    E_total=particles1_dict['E_total'],

                )
                self.record_injectedparticles(
                    species=self.species2,
                    w=particles2_dict['w'],
                    E_total=particles2_dict['E_total'],
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
            The conductor voltage V of the conductor the particle is ejected
            from must also be set for this object.

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

    def __init__(self, T, conductor=None, emission_type='thermionic'):
        """Default initialization for all Emitter objects.

        Arguments:
            T (float): Emitter temperature in Kelvin. Determines particle
                velocity distribution. If None, the temperature of the
                conductor will be used if one is specified.
            conductor (assemblies.Assembly): Conductor the emitter is attached
                to, used for recording initial voltages and energies. If None,
                V_e is set to 0. Since there's no current use case for this, a
                warning is printed.
            emission_type (str): Distribution function type used to sample
                velocities of the emitted particles. Must be defined in
                :func:`mewarpx.utils_store.util.get_velocities`. Defaults to
                'thermionic'.
        """
        super(Emitter, self).__init__()
        self.T = T
        if self.T is None and conductor is not None:
            self.T = conductor.T
        if self.T is None:
            raise ValueError(
                "No value for T given to the Emitter. An Emitter T must be "
                "specified directly, or on a conductor passed to the Emitter."
            )

        self.conductor = conductor
        if self.conductor is not None:
            if self.conductor.WF <= 0.:
                raise ValueError("Conductor WF must be set for emitters.")
        else:
            warnings.warn("No conductor set for emitter. Power will not be "
                          "correct.")

        self.emission_type = emission_type

    def getvoltage(self):
        if self.conductor is None:
            return 0.

        return self.conductor.getvoltage()

    def getvoltage_e(self):
        """Electrical voltage includes WF, eg the Fermi level voltage."""
        if self.conductor is None:
            return 0.

        return self.conductor.getvoltage() + self.conductor.WF

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
        raise NotImplementedError("Normal calculations must be implemented by "
                                  "Emitter sub-classes.")


class ZPlaneEmitter(Emitter):
    """This is the standard injection for a planar cathode."""

    geoms = ['Z', 'XZ', 'XYZ']

    def __init__(self, conductor, T=None, xmin=None, xmax=None,
                 ymin=None, ymax=None, transverse_fac=1.0, **kwargs):
        """Initialize an emitter for a planar cathode.

        Arguments:
            conductor (:class:`mewarpx.assemblies.Assembly`): Conductor object,
                used to obtain work function. Can later grab other variables
                from this conductor.
            T (float): Temperature in Kelvin for the emitter; determines
                velocities. If not specified the temperature of the conductor
                will be used.
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
                :func:`mewarpx.utils_store.util.get_velocities` for details.
            kwargs (dict): Any other keyword arguments supported by the parent
                Emitter constructor (such as "emission_type").
        """
        # Default initialization
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
            logger.info("x/y span is 1m for purposes of charge injection")
            x_range = 1.
            y_range = 1.
        if self.solver_geom == 'XZ':
            logger.info("y span is 1m for purposes of charge injection")
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

    def __init__(self, conductor, T=None, res_fac=5., transverse_fac=1.0,
                 **kwargs):
        """Construct the emitter based on conductor object and temperature.

        Arguments:
            conductor (mewarpx.assemblies object): Conductor to emit from.
            T (float): Temperature in Kelvin. If not specified the temperature
                of the conductor will be used.
            res_fac (float): Level of resolution beyond the grid resolution to
                use for calculating shape contours.
            transverse_fac (float): Scale the transverse energy distribution by
                this factor. Default 1. See
                :func:`mewarpx.utils_store.util.get_velocities` for details.
            kwargs (dict): Any other keyword arguments supported by the parent
                Emitter constructor (such as "emission_type").
        """
        # Default initialization
        super(ArbitraryEmitter2D, self).__init__(
            T=T, conductor=conductor, **kwargs
        )
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

        # level of 0.17 was chosen to keep the original ratio of 0.5:3 from warp
        # compared to 0.17:1 now in warpx
        # increasing the level causes particles to be injected inside cylinder
        self.contours = np.squeeze(skimage.measure.find_contours(
            inside, 0.17))

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
        partdist = self.dA * float(self.res_fac) / 2.

        pos = self.centers + marching_normal * partdist
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

        vels = np.column_stack(mwxutil.get_velocities(
            num_samples=npart, T=self.T, m=m,
            rseed=rseedv,
            transverse_fac=self.transverse_fac,
            emission_type=self.emission_type
        ))

        # Rotate velocities based on angle of normal
        newvels = self.convert_vel_zhat_nhat(
            vels, self.normal[self.contour_idx])
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

    def plot_contours(self):
        """Plots the contours generated for the assembly object and the
        assembly object. The object is plotted in yellow, and the contours
        are plotted in blue. The plot is saved in contours.png"""

        # calculate which tiles are inside of assembly object
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

        contours = np.array(skimage.measure.find_contours(inside, 0.17))

        # plot assembly object first
        assembly_cmap = colors.LinearSegmentedColormap.from_list('my_cmap',['white','#66c2a5'],256)
        fig, ax = plt.subplots()

        ax.imshow(inside, cmap=assembly_cmap, origin="lower")

        # plot contours
        for contour in contours:
            ax.plot(contour[:, 1], contour[:, 0], linewidth=2, color="#fc8d62")

        # set title and labels
        ax.set_title(f"{self.conductor.name} contour plot")

        x_range = [self.res_fac * mwxrun.zmin / mwxrun.dz, self.res_fac * mwxrun.zmax / mwxrun.dz]
        y_range = [self.res_fac * mwxrun.xmin / mwxrun.dx, self.res_fac * mwxrun.xmax / mwxrun.dx]

        x_step = mwxrun.dz / (self.res_fac)
        y_step = mwxrun.dx / (self.res_fac)

        minor_xticks = np.linspace(x_range[0], x_range[1], mwxrun.nz)
        minor_yticks = np.linspace(y_range[0], y_range[1], mwxrun.nx)

        major_xticks = np.linspace(x_range[0], x_range[1], 5)
        major_yticks = np.linspace(y_range[0], y_range[1], 5)

        ax.set_xlabel("Z (m)")
        ax.set_ylabel("X (m)")

        ax.set_xticks(major_xticks)
        ax.set_xticks(minor_xticks, minor=True)
        ax.set_xticklabels(np.round(major_xticks * x_step, 8), rotation=45)
        ax.set_yticks(major_yticks)
        ax.set_yticks(minor_yticks, minor=True)
        ax.set_yticklabels(np.round(major_yticks * y_step, 8))


        ax.grid(visible=True, which="minor")
        ax.set_aspect(mwxrun.dx/mwxrun.dz, adjustable='box')
        fig.tight_layout()
        fig.savefig(f"{self.conductor.name}_contour_plot.png")

class VolumeEmitter(BaseEmitter):

    """Parent class for volumetric particle injection coordinates.

        - ``volume`` gives the spatial volume in m^3
        - ``_get_xv_coords()`` implements the subclass-specific particle
          injection logic
    """

    volume = 0

    def __init__(self, T):
        """Default initialization for all Emitter objects.

        Arguments:
            T (float): Emitter temperature in Kelvin. Determines particle
                velocity distribution.
        """
        super(VolumeEmitter, self).__init__()
        self.T = T

    def getvoltage(self):
        """Ideally this is probably the local potential, but default to 0."""
        return 0.

    def getvoltage_e(self):
        """Ideally this is probably the local potential, but default to 0."""
        return self.getvoltage()


class CartesianVolumeEmitter(VolumeEmitter):

    """VolumeEmitter with bounds specified by x/y/z rectangular prism.

    A template ``_get_xv_coords()`` is also provided but can be overridden if
    needed. It calls ``_get_x_coords()``.
    """

    # RZ support should probably take into account the R wall, so we don't list
    # it as supported.
    geoms = ['Z', 'XZ', 'XYZ']

    def __init__(self, T, xmin=None, xmax=None, ymin=None, ymax=None,
                 zmin=None, zmax=None):
        """Initialize emitter boundaries."""
        super(CartesianVolumeEmitter, self).__init__(T=T)

        self.bounds = np.zeros((3, 2))
        for ii, (lim, defaultlim) in enumerate(
            zip([xmin, xmax, ymin, ymax, zmin, zmax],
                [mwxrun.xmin, mwxrun.xmax, mwxrun.ymin,
                 mwxrun.ymax, mwxrun.zmin, mwxrun.zmax])
        ):
            if lim is None:
                lim = defaultlim
            self.bounds[ii // 2, ii % 2] = lim

        self.volume = np.prod(self.bounds[:, 1] - self.bounds[:, 0])

        # Note the negation here will catch nans, checking <= 0 won't.
        if not (self.volume > 0):
            raise RuntimeError("Invalid warpX geometry limits.")

    def _get_xv_coords(self, npart, m, rseed):
        """Get velocities and call specialized function for position."""
        if rseed is not None:
            nprstate = np.random.get_state()
            np.random.seed(rseed)
            # rseedv is passed to get velocities. The basic rseed here is used
            # for positions, below.
            rseedv = np.random.randint(1000000000)
        else:
            rseedv = None

        x_coords = self._get_x_coords(npart)
        v_coords = mwxutil.get_velocities(
            npart, self.T, m=m, emission_type='random', rseed=rseedv)

        if rseed is not None:
            np.random.set_state(nprstate)

        return (
            x_coords[:, 0], x_coords[:, 1], x_coords[:, 2],
            v_coords[0], v_coords[1], v_coords[2]
        )


class UniformDistributionVolumeEmitter(CartesianVolumeEmitter):

    """Inject particles uniformly throughout the simulation at a specified
    temperature.
    """

    def _get_x_coords(self, npart):
        """Get coordinates uniformly distributed in space.

        rseed, if used, is handled by the parent function.
        """
        xyz_pos = [
            np.random.uniform(self.bounds[ii, 0], self.bounds[ii, 1],
                              npart)
            for ii in range(3)
        ]
        return np.array(xyz_pos).T


class ZSinDistributionVolumeEmitter(CartesianVolumeEmitter):

    """Vary density in z as a half-period sin wave."""

    def _get_x_coords(self, npart):
        """Get coordinates with sin distribution.

        rseed, if used, is handled by the parent function.
        """
        xpos = np.random.uniform(self.bounds[0, 0], self.bounds[0, 1], npart)
        ypos = np.random.uniform(self.bounds[1, 0], self.bounds[1, 1], npart)
        z_random_draw = np.random.random(npart)
        zpos = (
            np.arccos(1 - 2.0*z_random_draw) / np.pi
            * (self.bounds[2, 1] - self.bounds[2, 0])
            + self.bounds[2, 0]
        )
        return np.array([xpos, ypos, zpos]).T
