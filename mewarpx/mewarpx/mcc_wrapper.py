import glob, os
from pywarpx import picmi
from mewarpx import util as mwxutil

# For use later to sync with diode test template and access sim object in mwxrun
from mewarpx.mwxrun import mwxrun


class MCC():

    """Wrapper used to initialize Monte Carlo collision parameters"""

    def __init__(self, electron_species, ion_species, T_INERT,
                 P_INERT=None, N_INERT=None, scraper=None, **kwargs):
        """Initialize MCC parameters.

        Arguments:
            electron_species (picmi.Species): Species that will be producing the
                ions via impact ionization. This will normally be electrons.
            ion_species (picmi.Species): Ion species generated from ionization
                events. Charge state should be specified during Species
                construction. Also used to obtain the neutral mass.
            T_INERT (float): Temperature for injected ions in
                Kelvin.
            P_INERT (float): Pressure of the neutral "target" for
                impact ionization, in Torr. Assumed to be such that the density
                is much larger than both the electron and ion densities, so that
                the neutral dynamics can be ignored. Cannot be specified if
                N_INERT is specified.
            N_INERT (float): Neutral gas density in m^-3. Cannot be specified
                if P_INERT is specified.
            scraper (pywarpx.ParticleScraper): The particle scraper is
                instructed to save pid's for number of MCC events.
            **kwargs that can be included:
            exclude_collisions (list): A list of collision types to exclude.
        """
        self.electron_species = electron_species
        self.ion_species = ion_species
        self.T_INERT = T_INERT
        self.N_INERT = N_INERT
        self.P_INERT = P_INERT

        self.exclude_collisions = kwargs.get("exclude_collisions", None)
        if self.exclude_collisions is None:
            self.exclude_collisions = []

        if self.N_INERT is not None:
            # N and P cannot both be specified
            if self.P_INERT is not None:
                raise ValueError("Must specify N_INERT or P_INERT, not both")
            # if N is not None and P is None, everything is all good
        # N and P cannot both be unspecified
        elif self.P_INERT is None:
            raise ValueError("Must specify one of N_INERT or P_INERT")
        # set N using ideal gas law if only P is specified
        else:
            # convert from cm^-3 to m^-3
            self.N_INERT = (
                mwxutil.ideal_gas_density(self.P_INERT, self.T_INERT) * 1e6
            )

        self.scraper = scraper

        # Use environment variable if possible, otherwise look one
        # directory up from warpx
        path_name = os.environ.get(
            "MCC_CROSS_SECTIONS_DIR", os.path.join(
                mwxutil.mewarpx_dir, "../../../warpx-data/MCC_cross_sections"
            )
        )
        path_name = os.path.join(path_name, self.ion_species.particle_type)
        # include all collision processes that match species
        file_paths = glob.glob(os.path.join(path_name, "*.dat"))

        elec_collision_types = {
            "electron_scattering": "elastic",
            "excitation_1": "excitation1",
            "excitation_2": "excitation2",
            "ionization": "ionization",
        }
        ion_collision_types = {
            "ion_scattering": "elastic",
            "ion_back_scatter": "back",
            "charge_exchange": "charge_exchange"
        }
        required_energy = {
            "He": {
                "excitation_1": 19.82,
                "excitation_2": 20.61,
                "ionization": 24.55
            },
            "Ar": {
                "excitation_1": 11.5,
                "ionization": 15.7596112
            },
            "Xe": {
                "excitation_1": 8.315,
                "ionization": 12.1298431
            }
        }

        # build scattering process dictionaries
        elec_scattering_processes = {}
        ion_scattering_processes = {}

        for path in file_paths:
            file_name = os.path.basename(path)
            coll_key = file_name.split('.dat')[0]

            # if electron process
            if coll_key in elec_collision_types:
                coll_type = elec_collision_types[coll_key]
                # exclude collision type if specified
                if coll_type in self.exclude_collisions:
                    continue
                scatter_dict = {"cross_section": path}
                # add energy if needed
                ion = self.ion_species.particle_type
                if coll_key in required_energy[ion]:
                    scatter_dict["energy"] = required_energy[ion][coll_key]
                # specify species for ionization
                if coll_key == "ionization":
                    scatter_dict["species"] = self.ion_species
                elec_scattering_processes[coll_type] = scatter_dict

            # if ion process
            elif coll_key in ion_collision_types:
                coll_type = ion_collision_types[coll_key]
                # exclude collision type if specified
                if coll_type in self.exclude_collisions:
                    continue
                scatter_dict = {"cross_section": path}
                ion_scattering_processes[coll_type] = scatter_dict

            else:
                raise ValueError(
                    f"{path}: filename not recognized as an MCC cross-section "
                    "file. Please move outside this folder or end with "
                    "something other than .dat if it is not a cross-section "
                    "file."
                )

        # raise an error if no scattering processes exist
        if (not elec_scattering_processes) and (not ion_scattering_processes):
            raise ValueError(
                "No scattering processes for electron or ion species."
            )

        mwxrun.simulation.collisions = []

        if elec_scattering_processes:
            self.electron_mcc = picmi.MCCCollisions(
                name='coll_elec',
                species=self.electron_species,
                background_density=self.N_INERT,
                background_temperature=self.T_INERT,
                background_mass=self.ion_species.mass,
                scattering_processes=elec_scattering_processes
            )
            mwxrun.simulation.collisions.append(self.electron_mcc)

        if ion_scattering_processes:
            self.ion_mcc = picmi.MCCCollisions(
                name='coll_ion',
                species=self.ion_species,
                background_density=self.N_INERT,
                background_temperature=self.T_INERT,
                scattering_processes=ion_scattering_processes
            )
            mwxrun.simulation.collisions.append(self.ion_mcc)
