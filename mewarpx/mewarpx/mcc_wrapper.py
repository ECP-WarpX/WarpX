import glob, os
from pywarpx import picmi
from mewarpx import util as mwxutil

# For use later to sync with diode test template and access sim object in mwxrun
from mewarpx import mwxrun


class MCC():
    """Wrapper used to initialize MCC parameters"""

    def __init__(self, electron_species, ion_species, P_INERT,
                 T_INERT, scraper=None, **kwargs):
        """Initialize MCC"""
        self.electron_species = electron_species
        self.ion_species = ion_species
        self.P_INERT = P_INERT
        self.T_INERT = T_INERT
        self.N_INERT = mwxutil.ideal_gas_density(self.P_INERT, self.T_INERT)
        self.scraper = scraper

        # include all collision processes that match species
        path_name = os.path.join("../../../warpx-data/MCC_cross_sections", self.ion_species.particle_type)

        file_paths = glob.glob(path_name + "/*.dat")

        elec_collision_types = {
                            "electron_scattering.dat": "elastic",
                            "excitation_1.dat": "excitation1",
                            "excitation_2.dat": "excitation2",
                            "ionization.dat": "ionization",
                           }
        ion_collision_types = {
                            "ion_scattering.dat": "elastic",
                            "ion_back_scatter.dat": "back",
                            "charge_exchange.dat": "charge_exchange"
                            }
        requires_energy = {
                            "Ar": {
                                "excitation_1.dat": 11.5,
                                "ionization.dat": 15.7596112
                                },
                            "He": {
                                "excitation_1.dat": 19.82,
                                "excitation_2.dat": 20.61,
                                "ionization.dat": 24.55
                                },
                            "Xe": {
                                "excitation_1.dat": 8.315,
                                "ionization.dat": 12.1298431
                                }
                        }

        # build scattering process dictionaries
        elec_scattering_processes = {}
        ion_scattering_processes = {}

        for path in range(len(file_paths)):
            file_name = os.path.absename(path)
            # if electron process
            if file_name in elec_collision_types:
                scatter_dict = {"cross_section": path}
                # add energy if needed
                if file_name in requires_energy[self.ion_species.name]:
                    scatter_dict["energy"] = requires_energy[self.ion_species.name][file_name]
                # specify species for ionization
                if file_name == "ionization.dat":
                    scatter_dict["species"] = self.ion_species
                elec_scattering_processes[elec_collision_types[file_name]] = scatter_dict
            # if ion process
            elif file_name in ion_collision_types:
                scatter_dict = {"cross_section": path}
                ion_scattering_processes[ion_collision_types[file_name]] = scatter_dict

        self.mcc_electrons = picmi.MCCCollisions(
            name='coll_elec',
            species=self.electron_species,
            background_density=self.N_INERT,
            background_temperature=self.T_INERT,
            background_mass=self.ion_species.mass,
            scattering_processes=elec_scattering_processes
        )

        self.mcc_ions = picmi.MCCCollisions(
            name='coll_ion',
            species=self.ion_species,
            background_density=self.N_INERT,
            background_temperature=self.T_INERT,
            scattering_processes=ion_scattering_processes
        )
