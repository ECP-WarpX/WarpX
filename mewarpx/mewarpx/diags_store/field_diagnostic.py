"""Class for installing a field diagnostic with optional plotting"""

from mewarpx.mwxrun import mwxrun
from mewarpx.diags_store.diag_base import WarpXDiagnostic
from mewarpx import plotting

from pywarpx import callbacks, picmi
from mewarpx import plotting

import numpy as np
import matplotlib.pyplot as plt

import os
import yt
import glob
import warnings

class FieldDiagnostic(WarpXDiagnostic):
    def __init__(self, diag_steps, diag_data_list, grid, name, write_dir,
                 plot_on_diag_step=False, plot_data_list=None, plot_species_list=None,
                 post_processing=False,
                 **kwargs):
        """
        This class is a wrapper for both creating a picmi FieldDiagnostic,
        and optionally plotting parameters each diagnostic step.

        Arguments:
            diag_steps (int): Run the diagnostic with this period.
                Also plot on this period if enabled.
            diag_data_list (list (str)): A list of criteria to collect
                and print diagnostic information for. Accepted values can
                be found in the picmi standard for field diagnostics.
            grid (picmi grid): The grid to be passed into the picmi field diagnostic.
            name (str): The name of the diagnostic to be passed into the
                picmi field diagnostic.
            write_dir (str): The directory where diagnostic data will be written
                from the field diagnostic.
            plot_on_diag_step (bool): Whether or not to plot data on diagnostic steps.
            plot_data_list (list (str)): What parameters to plot if plotting on diagnostic steps.
            plot_species_list (list (str)): Plot rho for these species if plotting on diagnostic steps.
            post_processing (bool): Whether or not to plot data after simulation ends from any
                yt files generated during the run.
            kwargs: For a list of valid keyword arguments see diag_base.WarpXDiagnostic
        """
        self.diag_steps = diag_steps
        self.diag_data_list = diag_data_list
        self.grid = grid
        self.name = name
        self.write_dir = write_dir
        self.plot_on_diag_step = plot_on_diag_step
        self.plot_data_list = plot_data_list
        self.post_processing = post_processing

        if plot_species_list is not None:
            self.plot_species_list = plot_species_list
        else:
            self.plot_species_list = [species.name for species in mwxrun.simulation.species]

        super(FieldDiagnostic, self).__init__(diag_steps, **kwargs)

        self.add_field_diag()

        if self.plot_on_diag_step:
            if self.plot_data_list is None:
                raise ValueError(
                    "plot_on_diag_step was True but plot_data_list is None!")
            callbacks.installafterstep(self.plot_parameters)

        if self.post_processing:
            callbacks.installafterstep(self.check_for_end_of_sim)

    def add_field_diag(self):
        diagnostic = picmi.FieldDiagnostic(
            grid=self.grid,
            period=self.diag_steps,
            data_list=self.diag_data_list,
            name=self.name,
            write_dir=self.write_dir
        )

        mwxrun.simulation.add_diagnostic(diagnostic)

    def plot_parameters(self):
        if self.check_timestep():
            current_diag_num = int(mwxrun.get_it() / self.diag_steps)
            plot_name = "after_diag_step"
            if current_diag_num is not None:
                plot_name += "_" + str(current_diag_num)

            for param in self.plot_data_list:
                if param == 'rho':
                    data = mwxrun.get_gathered_rho_grid(include_ghosts=False)
                    if mwxrun.me == 0:
                        data = np.array(data[0])

                        fig, ax = plt.subplots(1, 1, figsize=(14, 14))
                        plotting.ArrayPlot(
                            array=data[:, :, 0],
                            template='rho', xaxis='z', yaxis='x', ax=ax,
                            plot_name="rho_" + plot_name
                        )

                    for species in self.plot_species_list:
                        data = mwxrun.get_gathered_rho_grid(include_ghosts=False, species_name=species)
                        if mwxrun.me == 0:
                            data = np.array(data[0])

                            fig, ax = plt.subplots(1, 1, figsize=(14, 14))
                            plotting.ArrayPlot(
                                array=data[:, :, 0],
                                template='rho', xaxis='z', yaxis='x', ax=ax,
                                plot_name="rho_" + species + "_" + plot_name
                            )

                elif param == 'phi':
                    data = mwxrun.get_gathered_phi_grid(include_ghosts=False)
                    if mwxrun.me == 0:
                        data = np.array(data[0])

                        fig, ax = plt.subplots(1, 1, figsize=(14, 14))
                        plotting.ArrayPlot(
                            array=data,
                            template='phi', xaxis='z', yaxis='x', ax=ax,
                            plot_name="phi_" + plot_name
                        )
                else:
                    warnings.warn("Attempted to plot a parameter for which plotting has not been implemented: " + param)

    def check_for_end_of_sim(self):
        if mwxrun.get_it() == mwxrun.simulation.max_steps:
            self.do_post_processing()

    def do_post_processing(self):
        if mwxrun.me == 0:
            constants = picmi.constants
            data_dirs = glob.glob(os.path.join(
                self.write_dir, self.name + '*'
            ))

            if len(data_dirs) == 0:
                raise RuntimeError("No data files found.")

            for datafolder in data_dirs:
                if "old" in datafolder:
                    continue

                print('Reading ', datafolder, '\n')
                ds = yt.load( datafolder )

                grid_data = ds.covering_grid(
                    level=0, left_edge=ds.domain_left_edge, dims=ds.domain_dimensions
                )

                for parameter in self.diag_data_list:
                    plot_name = os.path.join(self.write_dir, parameter + "_" + datafolder.replace(self.write_dir + "/" + self.name, ""))
                    data = np.mean(
                        grid_data[parameter].to_ndarray()[:, :, 0], axis=0
                    ) / constants.q_e

                    fig, ax = plt.subplots(1, 1, figsize=(14, 14))
                    template = "rho" if "rho" in parameter else "phi"
                    plotting.ArrayPlot(
                            array=data,
                            template=template, xaxis='z', yaxis='x', ax=ax,
                            plot_name=plot_name
                        )
