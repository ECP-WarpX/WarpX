"""Diagnostic code that wraps the picmi.ParticleDiagnostics class"""

import os
import glob
import logging
import warnings

import matplotlib.pyplot as plt
import yt

from pywarpx import callbacks, picmi

from mewarpx.mwxrun import mwxrun
from mewarpx.diags_store.diag_base import WarpXDiagnostic

logger = logging.getLogger(__name__)


class ParticleDiagnostic(WarpXDiagnostic):
    """Output particle diagnostics every diagnostic period and produce
    plots from the diagnostics

    Contains:
        add_particle_diag (function): Function to add the particle diagnostic to
            the mwxrun.simulation

        post_processing: plot particle diagnostic data, not yet implemented
    """

    def __init__(self, diag_steps, name=None, species=None,
                 data_list=None, write_dir='.',
                 post_processing=False, plot_data_list=None,
                 plot_species=None, **kwargs):

        """Initializes the picmi.ParticleDiagnostic and adds the diagnostic to
        the simulation

        Arguments:
            diag_steps (int): Number of steps between each diagnostic output

            name (str): name of the diag output folder, defaults to
                ``particle_diag``

            species (:class:`mewarpx.mepicmi.Species`): species in the
                simulation, if None then uses all particles in the simulation

            data_list (list str): list of attributes to be outputted by the
                diagnostic, default uses ``["position", "momentum", "weighting"]``

            write_dir (str): where to write particle diagnostic data

            post_process (bool): generate plots for each diagnostic data
                directory produced

            plot_data_list (list str): list of data to be plotted for each
                diagnostic step ("particle_position_x", "particle_position_y",
                "particle_position_z", "particle_momentum_x",
                "particle_momentum_y", "particle_momentum_z")

            plot_species (list str): list of species names to be plotted,
                defaults to all species if not specified. Name variable in
                :class:`mewarpx.mepicmi.Species` must be set for each species
                in the simulation

        """
        self.name = name
        self.species = species
        self.data_list = data_list
        self.write_dir = write_dir

        self.post_processing = post_processing
        self.plot_data_list = plot_data_list
        self.plot_species = plot_species

        if self.name is None:
            self.name = 'particle_diag'
        if self.species is None:
            self.species = mwxrun.simulation.species
        if self.data_list is None:
            self.data_list = ['position', 'momentum', 'weighting']
        if self.plot_species is None:
            self.plot_species = []
            for specimen in self.species:
                self.plot_species.append(specimen.name)

        super(ParticleDiagnostic, self).__init__(
            diag_steps=diag_steps, **kwargs)
        self.add_particle_diag()

        if self.post_processing:
            if self.plot_data_list is None:
                warnings.warn('Warning: post_process is True,'
                              'but plot_data_list is None!')
            else:
                callbacks.installafterstep(self.check_for_end_of_sim)

    def add_particle_diag(self):
        particle_diag = picmi.ParticleDiagnostic(
            name=self.name,
            period=self.diag_steps,
            species=self.species,
            data_list=self.data_list,
            write_dir=self.write_dir
        )

        mwxrun.simulation.add_diagnostic(particle_diag)

    def check_for_end_of_sim(self):
        if mwxrun.get_it() == mwxrun.simulation.max_steps:
            self.do_post_processing()

    def do_post_processing(self):
        if mwxrun.me == 0:
            data_dirs = glob.glob(os.path.join(
                self.write_dir, self.name + '*'))

            if len(data_dirs) == 0:
                raise RuntimeError(f'No diagnostic data '
                                   f'found in {self.write_dir}')

            for datafolder in data_dirs:
                print('Reading', datafolder, '\n')
                step = datafolder.split("_")[-1]
                yt_data = yt.load(datafolder)
                grid_data = yt_data.covering_grid(
                    level=0,
                    left_edge=yt_data.domain_left_edge,
                    dims=yt_data.domain_dimensions
                )
                for species_name in self.plot_species:
                    for param in self.plot_data_list:
                        field = (species_name, param)
                        if field not in yt_data.field_list:
                            warnings.warn(f'{field} '
                                           'not found in yt_data field list')
                        else:
                            raw_data = grid_data[(species_name, param)]
                            np_data = raw_data.to_ndarray()
                            fig, ax = plt.subplots(figsize=(14, 14))
                            ax.hist(np_data, bins=500, histtype='stepfilled',
                                    alpha=0.85)
                            ax.set_xlabel(param)
                            ax.set_ylabel("Counts")
                            ax.set_title(species_name)
                            out_path = os.path.join(
                                self.write_dir,
                                f'{species_name}_{param}_{step}.jpg',
                            )
                            plt.savefig(out_path)
