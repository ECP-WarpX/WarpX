"""Class for installing a field diagnostic with optional plotting"""

from mewarpx.mwxrun import mwxrun
from mewarpx.diags_store.diag_base import WarpXDiagnostic
from mewarpx.utils_store import mwxconstants, plotting

from pywarpx import callbacks, picmi

import numpy as np
import matplotlib.pyplot as plt

import os
import yt
import glob
import logging

logger = logging.getLogger(__name__)


class FieldDiagnostic(WarpXDiagnostic):

    FIELD_DIAG_DIR = "fields"

    def __init__(self, diag_steps, process_phi=True, process_E=False,
                 process_rho=True, species_list=None, plot=True,
                 barrier_slices=None, max_dim=16.0, min_dim=0.0, dpi=300,
                 install_field_diagnostic=False, post_processing=False,
                 **kwargs):
        """
        This class handles diagnostics for field quantities (output and
        plotting) typically of interest in Modern Electron simulations.
        Optionally a picmi FieldDiagnostic can also be installed.

        Arguments:
            diag_steps (int): Run the diagnostic with this period.
                Also plot on this period if enabled.
            process_phi (bool): If True, output phi. Default True.
            process_E (bool): If True, output E. Default False.
            process_rho (bool): If True, output rho. Default True.
            species_list (list): Optional list of picmi.Species objects for
                which output density. If not specified all species will be
                tabulated.
            plot (bool): If True, also generate plots. Default True.
            barrier_slices (list): If provided, also plot potential energy
                slices as a function of z for each value of x/r in the list.
                Units are in meters.
            max_dim (float): Maximum figure dimension in inches for field plots.
            min_dim (float): Minimum figure dimension in inches for field plots.
            dpi (int): Resolution to use for saved plots
            install_field_diagnostic (bool): If true, install a picmi
                FieldDiagnostic. All parameters for the diagnostic should be
                passed as keyword arguments.
            post_processing (bool): Whether or not to plot data after simulation
                ends from any yt files generated during the run.
            kwargs: For a list of valid keyword arguments see
                diag_base.WarpXDiagnostic
        """
        self.diag_steps = diag_steps
        self.write_dir = os.path.join(self.DIAG_DIR, self.FIELD_DIAG_DIR)
        # field quantities to process
        self.process_phi = process_phi
        self.process_E = process_E
        self.process_rho = process_rho
        self.species_list = species_list
        if self.species_list is None:
            self.species_list = mwxrun.simulation.species
        self.barrier_slices = barrier_slices

        self.plot = plot
        self.max_dim = max_dim
        self.min_dim = min_dim
        self.dpi = dpi
        self.a_ax = 'z'
        self.o_ax = 'x'

        super(FieldDiagnostic, self).__init__(
            diag_steps, kwargs.pop('diag_step_offset', 0),
            kwargs.pop('extended_interval_level', None),
            kwargs.pop('manual_timesteps', None)
        )

        self.install_field_diagnostic = install_field_diagnostic
        self.post_processing = post_processing
        self.kwargs = kwargs

        callbacks.installafterstep(self.fields_diag)

        if self.install_field_diagnostic:
            self.add_field_diag()

    def add_field_diag(self):
        """Add a picmi FieldDiagnostic to the simulation."""
        mwxrun.simulation.add_diagnostic(
            picmi.FieldDiagnostic(
                grid=self.kwargs['grid'],
                period=self.kwargs.pop('period', self.diag_steps),
                data_list=self.kwargs['data_list'],
                name=self.kwargs['name'],
                write_dir=self.kwargs.pop('write_dir', self.write_dir)
            )
        )

    def fields_diag(self):
        """Function to process (get, plot and save) field quantities. This
        function is called on every step, but only executes if check_timestep()
        evaluated to True.
        """
        if (self.post_processing
            and (mwxrun.get_it() == mwxrun.simulation.max_steps)
        ):
            self.do_post_processing()

        if not self.check_timestep():
            return

        logger.info("Analyzing fields...")
        self.it = mwxrun.get_it()

        if self.process_phi:
            data = mwxrun.get_gathered_phi_grid(include_ghosts=False)
            self.process_field(
                data=data,
                titlestr='Electrostatic potential',
                plottype='phi', draw_image=True, default_ticks=True,
                draw_contourlines=False)
            # Optionally generate barrier index plot if requested
            if self.plot and (self.barrier_slices is not None):
                self.plot_barrier_slices(data, self.barrier_slices)
        if self.process_E:
            raise NotImplementedError("E-field processing not yet implemented.")
            self.process_field(
                data=None,
                titlestr='Electric field strength',
                plottype='E', draw_image=True, default_ticks=True,
                draw_contourlines=False)
        if self.process_rho:
            # assume that rho_fp still holds the net charge density
            data = (
                mwxrun.get_gathered_rho_grid(include_ghosts=False)
                * 1e-6
            )[:,:,0]
            self.process_field(
                data=data,
                titlestr='Net charge density',
                plottype='rho', draw_image=True, default_ticks=True,
                draw_contourlines=False)

            # deposit the charge density for each species
            for species in self.species_list:
                data = (
                    mwxrun.get_gathered_rho_grid(
                        species_name=species.name, include_ghosts=False
                    ) / species.sq * 1e-6
                )[:,:,0]
                self.process_field(
                    data=data,
                    titlestr=f'{species.name} particle density',
                    plottype='n', draw_image=True, default_ticks=True,
                    draw_contourlines=False)

        logger.info("Finished analyzing fields")

    def process_field(self, data, titlestr, plottype=None,
                      **kwargs):
        """Save given field to file, and optionally plot it as well.

        Other kwargs are passed on to plotting.

        Arguments:
            data (numpy.ndarray): Array to output
            titlestr (string): String to use as first part of the file name.
                Also used for the title in plotting.
            plottype: Choose 'phi', 'E', or 'rho' to set correct labels
        """

        if mwxrun.me == 0:
            try:
                fileprefix = self.get_fileprefix(titlestr)
                # Save full field
                np.save(fileprefix + '.npy', data)
            except Exception as exc:
                logger.warning(f"Saving {titlestr} failed with error {exc}")

            try:
                # Plot if desired, and array is not all-0.
                if self.plot and not np.all(data == 0.):
                    self.plot_field(data=data, plottype=plottype,
                                    titlestr=titlestr, **kwargs)
            except Exception as exc:
                logger.warning(f"Plotting {titlestr} failed with error {exc}")

    def plot_barrier_slices(self, phi_array, slices, **kwargs):
        """Plot 1D potential energy slices from the phi array in order to
            visualize changes in the barrier index during the simulation.

        Other kwargs are passed on to plotting.

        Arguments:
            phi_array (numpy.ndarray): Full electrostatic potential array
            slices (float or list of floats): Positions (in m) along x/r for
                plotting barrier index lineouts.
        """
        if mwxrun.me == 0:
            try:
                self.plot_field(data=phi_array, plottype='barrier',
                                titlestr="Barrier index", plot1d=True,
                                points=slices, xaxis='z', yaxis='x', **kwargs)
            except Exception as exc:
                logger.warning(
                    f"Plotting barrier index failed with error {exc}"
                )

    def plot_field(self, data, plottype, titlestr, plot1d=False, **kwargs):
        """Plot given field and save to file as both pdf and png.

        Other kwargs are passed on to plotting.

        Arguments:
            data (numpy.ndarray): Full array to plot
            plottype: Choose 'phi', 'E', 'barrier' or 'rho' to set correct
                labels
            titlestr (string): Title for plot and filename.
            plot1d (bool): Used to determine figure size and plotting call.
        """
        # kwargs specified by user in initialization overwrite local kwargs
        kwargs.update(self.kwargs)
        if plot1d:
            fig, ax = plt.subplots()
        else:
            figsize = plotting.get_figsize_from_warpx(
                max_dim=self.max_dim, min_dim=self.min_dim
            )
            fig, ax = plt.subplots(1, figsize=figsize)

        fileprefix = self.get_fileprefix(titlestr)
        titleline2 = kwargs.pop('titleline2', f'Step {self.it:d}')

        xaxis = kwargs.pop('xaxis', self.a_ax)
        yaxis = kwargs.pop('yaxis', self.o_ax)

        plotting.ArrayPlot(array=data, template=plottype,
                           titlestr=titlestr, titleline2=titleline2,
                           plot1d=plot1d, xaxis=xaxis, yaxis=yaxis,
                           **kwargs)

        fig.set_tight_layout(True)
        if self.kwargs.get('save_pdf', True):
            fig.savefig(fileprefix + '.pdf', dpi=self.dpi)
        fig.savefig(fileprefix + '.png', dpi=self.dpi)
        plt.close()

    def get_fileprefix(self, title):
        """Return filepath except for the filetype.

        Arguments:
            title (str): String for title; whitespace will become underscores
        """
        return os.path.join(
            self.write_dir, '_'.join(title.split() + [f"{self.it:010d}"])
        )

    def do_post_processing(self):
        if mwxrun.me == 0:
            # we need to overwrite self.plot and the function
            # self.get_fileprefix in order to properly do the post-process
            # plotting
            self.plot = True
            self.get_fileprefix = (
                lambda title: os.path.join(self.write_dir, title)
            )

            data_dirs = glob.glob(os.path.join(
                self.write_dir, self.kwargs['name'] + '*'
            ))

            if len(data_dirs) == 0:
                raise RuntimeError("No data files found.")

            for datafolder in data_dirs:
                if "old" in datafolder:
                    continue

                logger.info(f"Reading {datafolder}\n")
                ds = yt.load( datafolder )

                timestep = datafolder[-5:]
                grid_data = ds.covering_grid(
                    level=0, left_edge=ds.domain_left_edge,
                    dims=ds.domain_dimensions
                )

                for parameter in self.kwargs['data_list']:
                    plot_name = os.path.join(
                        self.write_dir, parameter + "_" +
                        datafolder.replace(
                            os.path.join(self.write_dir, self.kwargs['name']),
                            ""
                        )
                    )
                    data = np.array(grid_data[parameter].to_ndarray()[:, :, 0])

                    if parameter == "phi":
                        plottype = "phi"
                    else:
                        plottype = "rho"
                        data = data / mwxconstants.e * 1e-6

                    self.process_field(
                        data=data,
                        titlestr=parameter+"_"+timestep,
                        plottype=plottype, draw_image=True, default_ticks=True,
                        draw_contourlines=False
                    )
