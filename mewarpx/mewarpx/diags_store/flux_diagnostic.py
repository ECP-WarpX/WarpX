import collections
import copy
import warnings
import os

import dill
import matplotlib.pyplot as plt
import numba
import numpy as np
import pandas
import glob

from mewarpx.mwxrun import mwxrun
import mewarpx.utils_store.util as mwxutil
from mewarpx.diags_store import diag_base, timeseries

from pywarpx import callbacks, _libwarpx

class ParticleCSVDiag(diag_base.WarpXDiagnostic):

    """Handles writing out, and accumulating real-time, a per-step record of
    charge injected or deposited. Provides a base for specific implementations.
    At present species_name is a required column and must be str.
    """

    # Columns not to accumulate internally
    columns_no_accumulate = ['t', 'step', 'species_id']

    def __init__(self, diag_steps, write_dir=None, save_name=None,
                 **kwargs):
        """Basic initialization.
        If ``write_dir`` or ``save_name`` are ``None``, no data is saved to file
        and only accumulated data is kept.

        Arguments:
            diag_steps (int): Number of steps between each output
            write_dir (string): Directory to write CSV files to.
            save_name(string): Full filename within the directory to write CSV
                files to.
            kwargs: See :class:`mewarpx.diags_store.diag_base.WarpXDiagnostic`
                for more timing options.
        """
        self.write_dir = write_dir
        self.save_name = save_name
        self.save = (
            (self.write_dir is not None) and (self.save_name is not None)
        )

        if self.save:
            # Initialize variables
            if mwxrun.me == 0:
                mwxutil.mkdir_p(self.write_dir)
            self.save_path = os.path.join(self.write_dir, self.save_name)

        self.accumulate_dict = {}

        super(ParticleCSVDiag, self).__init__(
            diag_steps=diag_steps, **kwargs)

    def _charge_accum_diag(self):
        """Present only for backwards compatibility with runs that installed
        this diagnostic.
        """
        pass

    def charge_accum_diag(self):
        """Generate flux dataframe; write to CSV if requested.
        Called by the FluxDiag owning this at appropriate timesteps.
        """
        df = self._collect_dataframe()
        self._accumulate_columns(df)
        if self.save:
            self._write_dataframe(df)

        return df

    def _collect_dataframe(self):
        """Gather the dataframe of per-step information and return it."""
        raise NotImplementedError("Must implement in child object.")

    def _accumulate_columns(self, df, accumulate_dict=None):
        """Track summed quantities only
        accumulate_dict not None will update that dict rather than the
        built-in dict.
        """
        if accumulate_dict is None:
            accumulate_dict = self.accumulate_dict

        species_list = np.unique(df['species_id'])

        for column in df.columns:
            if column not in self.columns_no_accumulate:
                for species_id in species_list:
                    accumulate_dict[(column, species_id)] = (
                        accumulate_dict.setdefault((column, species_id), 0.)
                        + np.sum(df[df['species_id'] == species_id][column])
                    )

    def get_species_list(self):
        """Get sorted list of all species names recorded by this diagnostic."""
        species_list = sorted(list(set(
            [x for _, x in list(self.accumulate_dict.keys())]
        )))
        return species_list

    def _write_dataframe(self, df):
        """Write the actual dataframe to file, append if possible."""
        if mwxrun.me == 0:
            # Write a header only if it's a new file
            if os.path.exists(self.save_path):
                df.to_csv(self.save_path, mode='a', header=False, index=False)
            else:
                df.to_csv(self.save_path, mode='w', header=True, index=False)

    def get_updated_accumulators(self):
        """Return up-to-date accumulated values in a separate dictionary,
        without clearing and/or writing the values.
        """
        acdict = self.accumulate_dict.copy()
        df = self._collect_dataframe(clear=False)
        self._accumulate_columns(df, accumulate_dict=acdict)
        return acdict


class InjectorFluxDiag(ParticleCSVDiag):

    """Handles writing out charge absorption by injectors."""

    def __init__(self, diag_steps, injector, write_dir, **kwargs):
        """Generate and install function to write out charge accumulation.

        Arguments:
            diag_steps (int): Number of steps between each output
            injector (emission.Injector): The particle injector stores
                intermediate values.
            write_dir (string): Directory to write CSV files to.
            kwargs: See :class:`mewarpx.diags_store.diag_base.WarpXDiagnostic`
                for more timing options.
        """
        # Initialize variables
        self.injector = injector

        # Register with the injector, make sure we don't have two diagnostics
        # attached, and initialize with fields.
        if self.injector.injector_diag is not None:
            raise RuntimeError("Cannot attach two InjectorFluxDiag objects "
                               "to the same injector.")
        self.injector.injector_diag = self
        self.injector.init_injectedparticles(self.injector.fields)

        # Setup write-out capability
        save_name = self.injector.name + '_injected.csv'

        super(InjectorFluxDiag, self).__init__(
            diag_steps=diag_steps,
            write_dir=write_dir,
            save_name=save_name,
            **kwargs
        )

    def _collect_dataframe(self, clear=True):
        partdict = self.injector.get_injectedparticles(clear=clear)
        # Convert species_name & step to int, because they come out as floats.
        partdict['species_id'] = partdict['species_id'].astype(int)
        partdict['step'] = partdict['step'].astype(int)
        df = pandas.DataFrame(partdict, columns=list(partdict.keys()))

        return df


class SurfaceFluxDiag(ParticleCSVDiag):

    """Handles writing out for a single surface.

    This class is used by FluxDiag; rarely used directly.
    """

    def __init__(self, diag_steps, surface, write_dir, **kwargs):
        """Initialize surface-specific features.

        Arguments:
            diag_steps (int): Number of steps between each output
            surface (mewarpx.Assembly object): The assembly object in which
                particles are scraped.
            write_dir (string): Directory to write CSV files to.
            kwargs: See :class:`mewarpx.diags_store.diag_base.WarpXDiagnostic`
                for more timing options.
        """
        self.surface = surface

        # Register with the surface, make sure we don't have two diagnostics
        # attached, and initialize with fields.
        if self.surface.scraper_diag is not None:
            raise RuntimeError("Cannot attach two SurfaceFluxDiag objects "
                               "to the same surface.")
        self.surface.scraper_diag = self
        self.surface.init_scrapedparticles(self.surface.fields)

        save_name = self.surface.name + '_scraped.csv'

        super(SurfaceFluxDiag, self).__init__(
            diag_steps=diag_steps,
            write_dir=write_dir,
            save_name=save_name,
            **kwargs
        )

        # install callback that will have the assembly object check the
        # particle buffer and record the scraped particles data
        # TODO change this to only process the particle buffer on a diagnostic
        # step instead of after every step
        callbacks.installbeforestep(_libwarpx.libwarpx.warpx_clearParticleBoundaryBuffer)
        callbacks.installbeforeEsolve(self.surface.record_scrapedparticles)

    def _collect_dataframe(self, clear=True):
        partdict = self.surface.get_scrapedparticles(clear=clear)
        # Convert species_id & step to int, because they come out as floats.
        partdict['species_id'] = partdict['species_id'].astype(int)
        partdict['step'] = partdict['step'].astype(int)
        df = pandas.DataFrame(partdict, columns=list(partdict.keys()))

        return df


class FluxDiagBase(diag_base.WarpXDiagnostic):

    """Handle generic printing & plotting of charge injection and absorption.
    This Base class should work with both postprocessing and during-run
    analyses.
    """

    default_printed_qtys = collections.OrderedDict([
        ('J', {'by_component': 'all',
               'total': 'all',
               'description': 'Current',
               'units': 'A/cm^2'}),
        ('P', {'by_component': None,
               'total': 'net',
               'description': 'Power',
               'units': 'W/cm^2'}),
        ('dQ', {'by_component': None,
                'total': None,
                'description': 'Heat transfer',
                'units': 'W/cm^2'}),
        ('n', {'by_component': None,
               'total': None,
               'description': 'Macroparticle rate',
               'units': 'particles per step'})
    ])

    def __init__(self, diag_steps, runinfo,
                 write_dir='diags/fluxes',
                 overwrite=True,
                 sig_figs=6,
                 printed_qtys=None,
                 fullhist_dict=None,
                 ts_dict=None,
                 **kwargs):
        """Generate and install function to write out fluxes.

        Arguments:
            diag_steps (int): Number of steps between each output
            runinfo (:class:`metools.runinfo.RunInfo`): RunInfo object is used
                to get the species, injectors, surfaces, and system area.
            write_dir (string): Directory to write CSV files and plots to.
            overwrite (bool): If True the dill pickled save file will overwrite
                the previous diagnostic period's saved file.
            sig_figs (int): Number of significant figures in text output.
                Default 6.
            printed_qtys (dict): Override individual values of
                default_printed_qtys; same input format but keys can be omitted
                to use defaults.
            fullhist_dict (dict): Dictionary of timeseries for the full run
            ts_dict (dict): Dictionary of timeseries for the last 8 steps
            kwargs: See :class:`mewarpx.diags_store.diag_base.WarpXDiagnostic`
                for more timing options.
        """
        # Save input variables
        self.write_dir = write_dir
        self.overwrite = overwrite
        self.format_str = '{:.%dg}' % sig_figs
        self.printed_qtys = copy.deepcopy(self.default_printed_qtys)
        if printed_qtys is not None:
            mwxutil.recursive_update(self.printed_qtys, printed_qtys)

        self.component_list = runinfo.component_list
        self.species_list = [spec.name for spec in mwxrun.simulation.species]

        self.fullhist_dict = fullhist_dict
        if self.fullhist_dict is None:
            self.fullhist_dict = collections.OrderedDict()

        self.ts_dict = ts_dict
        if self.ts_dict is None:
            self.ts_dict = collections.OrderedDict()

        super(FluxDiagBase, self).__init__(
            diag_steps=diag_steps,
            **kwargs
        )

    def refresh_species(self, runinfo):
        """When species are updated, this can be used to re-synchronize an
        existing FluxDiag with the new information. Will not update surfaces or
        injectors, however.

        Note:
            In the future, it may be best to replace this with something that
            more fully resets all references to runinfo, and/or just always
            explicitly reference runinfo rather than getting information from
            it at init. However, for now this is the quickest way to proceed.

        Arguments:
            runinfo (:class:`metools.runinfo.RunInfo`): An updated RunInfo
                object is solely used to get new species name information.
        """
        self.species_list = [spec.name for spec in mwxrun.simulation.species]

    def save(self):
        """Save only the critical variables to a pickle file. Since this can be
        done often, a specific function is used to minimize file size.
        """
        if self.overwrite:
            filename = 'fluxdata.dpkl'
        else:
            filename = 'fluxdata_{:010d}.dpkl'.format(mwxrun.get_it())

        filepath = os.path.join(self.write_dir, filename)
        dict_to_save = {
            'write_dir': self.write_dir,
            'format_str': self.format_str,
            'printed_qtys': self.printed_qtys,
            'fullhist_dict': self.fullhist_dict,
            'ts_dict': self.ts_dict,
        }
        with open(filepath, 'wb') as pfile:
            dill.dump(dict_to_save, pfile)

    def print_fluxes(self, tsdict):
        """Print the float results of calculations for the given dictionary of
        timeseries.

        Arguments:
            tsdict (dict of (str, str, str)): First string is either 'inject' or
                'scrape'; second string is the group name (eg cathode, accgrid);
                third value is species_name for the species. Note any cropping etc.
                should be done before this dict is passed in!

        Returns:
            resultstr (str): A string that can be printed for the performance.
        """
        emit_full = timeseries.concat_crop_timeseries([
            val for key, val in tsdict.items() if "inject" in key
        ])
        collect_full = timeseries.concat_crop_timeseries([
            val for key, val in tsdict.items() if "scrape" in key
        ])

        resultstr = ""
        for qty, flags in self.printed_qtys.items():
            if flags['by_component']:
                for component in self.component_list:
                    componentstr = self.runinfo.component_labels.get(
                        component, component
                    )

                    for species_name in self.species_list:
                        prefixstr = " ".join(
                            [componentstr, species_name.title(), flags['description']]
                        )
                        emit_ts = tsdict.get(('inject', component, species_name), None)
                        collect_ts = tsdict.get(('scrape', component, species_name), None)
                        suffixstr = flags['units']
                        net_only = (flags['by_component'] == 'net')
                        resultstr += self.print_fluxset(
                            prefixstr=prefixstr, key=qty, suffixstr=suffixstr,
                            emit_ts=emit_ts, collect_ts=collect_ts,
                            net_only=net_only
                        )

            if flags['total']:
                prefixstr = " ".join(["Total", flags['description']])
                suffixstr = flags['units']
                net_only = (flags['total'] == 'net')
                resultstr += self.print_fluxset(
                    prefixstr=prefixstr, key=qty, suffixstr=suffixstr,
                    emit_ts=emit_full, collect_ts=collect_full,
                    net_only=net_only
                )

        return resultstr

    def print_fluxset(self, prefixstr, key, suffixstr, emit_ts=None,
                      collect_ts=None, net_only=False):
        """Print emitted and collected flux for specific timeseries and specific
        flux quantity.

        Arguments:
            prefixstr (str): String that goes before numerical value. Should
                include a label for the key.
            key (str): Key for the quantity in the timeseries to print.
            suffixstr (str): String after the numerical value, usually units.
            emit_ts (:class:`mewarpx.diags_store.timeseries.Timeseries`):
                Emitted flux timeseries.
            collect_ts (:class:`mewarpx.diags_store.timeseries.Timeseries`):
                Collected flux timeseries.
            net_only (bool): Only print net value, not emitted and collected.

        Returns:
            resultstr (str): A string that can be printed for the performance.
        """
        print_emission = (emit_ts is not None) and (not net_only)
        print_collection = (collect_ts is not None) and (not net_only)
        print_net = (
            # Only print net and at least one ts is not None
            (net_only and ((emit_ts is not None) or (collect_ts is not None)))
            # Print net because have both emission and collection
            or (print_emission and print_collection)
        )

        emit_value = 0.
        if emit_ts is not None:
            emit_value = emit_ts.get_averagevalue_by_key(key)

        collect_value = 0.
        if collect_ts is not None:
            collect_value = collect_ts.get_averagevalue_by_key(key)

        resultstr = ""
        if print_emission:
            resultstr += " ".join([
                prefixstr,
                "Emitted:",
                self.format_str.format(emit_value),
                suffixstr])
            resultstr += "\n"

        if print_collection:
            resultstr += " ".join([
                prefixstr,
                "Collected:",
                self.format_str.format(collect_value),
                suffixstr])
            resultstr += "\n"

        if print_net:
            resultstr += " ".join([
                prefixstr,
                "Net:",
                self.format_str.format(emit_value + collect_value),
                suffixstr])
            resultstr += "\n"

        return resultstr

    def plot_fluxes(self, ts_dict, save=False):
        """
        Arguments:
            save (bool): If True, save and close figure. write_dir must be
                defined. If False, leave figure open and return it.
        """
        fig, axlist = plt.subplots(2, 2, figsize=(14, 8.5))

        # List of axes properties
        axlist = [x for y in axlist for x in y]
        qty_list = [
            {'key': 'J',
             'ylabel': r'Current (A/$\mathrm{cm}^2$)',
             'title': 'Simulation currents into system'},
            {'key': 'P',
             'ylabel': r'Power (W/$\mathrm{cm}^2$)',
             'title': 'Component power production'},
            {'key': 'dQ',
             'ylabel': r'Heat (W/$\mathrm{cm}^2$)',
             'title': 'Component heat transfer into system'},
            {'key': 'n',
             'ylabel': 'Macroparticles',
             'title': 'Simulation particles handled'},
        ]

        try:
            label_list = [
                '_'.join([
                    x[0].title(),
                    x[1].title(),
                    x[2].title()
                ])
                for x in list(ts_dict.keys())
            ]
        except KeyError:
            print(
                "NOTE: KeyErrors in flux plotting can be caused by creating "
                "species after RunInfo is saved (eg by initiating a "
                "TraceParticleInjector after RunInfo). All species should be "
                "initiated before RunInfo."
            )
            raise

        if mwxrun.get_it() * mwxrun.get_dt() < 0.5e-6:
            xlabel = 'Time (ns)'
        else:
            xlabel = r'Time ($\mu$s)'

        for ax, qtydict in zip(axlist, qty_list):
            typekey = qtydict['key']
            timeseries.TimeseriesPlot(
                array_list=[
                    (
                        label,
                        ts_dict[x].get_timeseries_by_key(typekey)
                    )
                    for (label, x) in
                    zip(label_list, list(ts_dict.keys()))
                ],
                ax=ax,
                xlabel=xlabel,
                ylabel=qtydict['ylabel'],
                title=qtydict['title'],
                titlesize=18,
                labelsize=16,
                alpha=0.7,
                legend=False
            )

        fig.legend(ax.get_lines(), label_list, 'lower center', fontsize=16,
                   frameon=True, ncol=3)
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.13 + int((len(label_list)-1) / 3) * 0.04,
                            hspace=0.33)

        if save:
            fig.savefig(
                os.path.join(
                    self.write_dir,
                    'flux_plots_{:010d}.png'.format(mwxrun.get_it())
                ), dpi=300
            )
            plt.close(fig)

        return fig

    def get_net_flux_timeseries(self, electrode_name, flux_type='J'):
        """Sum the emitted and absorbed flux across all species for a given
        electrode, and return a timeseries array containing the net flux from
        the electrode surface into the simulation volume.

        Arguments:
            electrode_name (str): Surface name such as 'cathode' or 'anode'.
                Must be present in the RunInfo injector_dict and surface_dict in
                order to be a valid key.
            flux_type (str): Either 'J', 'P', 'dQ', or 'n'. Defaults to 'J'.

        Returns:
            flux_timeseries (np.ndarray): n x 2 array containing the time in
                seconds in first column and the net flux from the electrode in
                second column.
        """
        flux_data = timeseries.concat_crop_timeseries([
            val for key, val in self.fullhist_dict.items()
            if electrode_name in key
        ])

        if flux_data is None:
            raise ValueError(
                'No flux data found for "%s" electrode!' % electrode_name
            )

        return flux_data.get_timeseries_by_key(flux_type)


class FluxDiagnostic(FluxDiagBase):

    """Handles writing out charge injection and absorption."""

    def __init__(self, diag_steps, runinfo,
                 write_dir='diags/fluxes', overwrite=True,
                 history_maxlen=5000, sig_figs=6,
                 printed_qtys=None,
                 check_charge_conservation=True,
                 print_per_diagnostic=True, print_total=False,
                 plot=True, profile_decorator=None, **kwargs):
        """Generate and install function to write out fluxes.

        Arguments:
            diag_steps (int): Number of steps between each output
            runinfo (:class:`metools.runinfo.RunInfo`): RunInfo object is used
                to get the species, injectors, surfaces, and system area.
            write_dir (string): Directory to write CSV files and plots to.
            overwrite (bool): If True the dill pickled save file will overwrite
                the previous diagnostic period's saved file.
            history_maxlen (int): Maximum length of full history to keep. If
                this is exceeded, history is resampled to a 2x lower frequency.
                Default 5000.
            sig_figs (int): Number of significant figures in text output.
                Default 6.
            printed_qtys (dict): Override individual values of
                default_printed_qtys; same input format but keys can be omitted
                to use defaults.
            check_charge_conservation (bool): Whether to check if charge is
                conserved in simulation.
            print_per_diagnostic (bool): Whether to print current results for
                the latest diagnostic period.
            print_total (bool): Whether to print total history of fluxes after a
                diagnostic period.
            plot (bool): Whether to save a plot of fluxes after each diagnostic
                period.
            profile_decorator (decorator): A decorator used to profile the
                timeseries update methods and related functions.
            kwargs: See :class:`metools.diags_store.diag_base.WarpDiagnostic`
                for more timing options.
        """
        # Save input variables
        self.runinfo = runinfo
        self.history_maxlen = history_maxlen
        self.history_dt = mwxrun.get_dt()
        self.check_charge_conservation = check_charge_conservation
        self.print_per_diagnostic = print_per_diagnostic
        self.print_total = print_total
        self.plot = plot

        if profile_decorator is not None:
            self.update_ts_dict = profile_decorator(self.update_ts_dict)
            self.update_fullhist_dict = (
                profile_decorator(self.update_fullhist_dict)
            )
            self.print_fluxes = profile_decorator(self.print_fluxes)
            self.plot_fluxes = profile_decorator(self.plot_fluxes)

        self.injector_dict = runinfo.injector_dict
        self.surface_dict = runinfo.surface_dict

        self.diags_dict = collections.OrderedDict()

        for key, val in self.injector_dict.items():
            self.injector_dict[key] = mwxutil.return_iterable(val)
            self.diags_dict[('inject', key)] = [
                InjectorFluxDiag(diag_steps=diag_steps,
                                 injector=injector,
                                 write_dir=write_dir,
                                 **kwargs)
                for injector in self.injector_dict[key]
            ]

        for key, val in self.surface_dict.items():
            self.surface_dict[key] = mwxutil.return_iterable(val)
            self.diags_dict[('scrape', key)] = [
                SurfaceFluxDiag(diag_steps=diag_steps,
                                surface=surface,
                                write_dir=write_dir,
                                **kwargs)
                for surface in self.surface_dict[key]
            ]

        # Initialize other variables
        # Note that at present, last_run_step is only updated for processor 0.
        self.last_run_step = 0

        super(FluxDiagnostic, self).__init__(
            diag_steps=diag_steps,
            runinfo=runinfo,
            write_dir=write_dir,
            overwrite=overwrite,
            sig_figs=sig_figs,
            printed_qtys=printed_qtys,
            **kwargs
        )

        callbacks.installafterstep(self._flux_ana)

    def _flux_ana(self):
        """Perform the calculation and processing of current data from the
        dataframes.
        The Timeseries object is used to store multiple keys for one collection
        of injectors and/or surfaces and one species. To facilitate that, we go
        by diagnostic object, create a timeseries by species, and concatenate
        timeseries by group of objects.
        We also keep a from-the-beginning tally of similar timeseries, so we can
        append the present timeseries to them, and then resample those if
        needed.
        """
        if self.check_timestep():
            self.update_ts_dict()

            if mwxrun.me == 0:
                self.update_fullhist_dict()

                if self.check_charge_conservation:
                    self._check_charge_conservation()

                if self.print_per_diagnostic:
                    print("THIS DIAGNOSTIC PERIOD:")
                    print(self.print_fluxes(self.ts_dict))
                if self.print_total:
                    print("TOTAL HISTORY:")
                    print(self.print_fluxes(self.fullhist_dict))

                if self.plot:
                    self.plot_fluxes(self.fullhist_dict, save=True)

                self.save()

            self.last_run_step = mwxrun.get_it()

    def update_ts_dict(self):
        """Run early in flux analysis to get this diagnostic period's
        timeseries.
        """
        self.ts_dict = collections.OrderedDict()

        for (keytype, key), diaglist in self.diags_dict.items():
            for diagobj in diaglist:
                df = diagobj.charge_accum_diag()
                # only need to hold a copy of the timeseries on root
                if mwxrun.me == 0:
                    species_list = diagobj.get_species_list()
                    for sp in species_list:
                        subdf = df[df['species_id'] == sp]
                        ts = FluxCalcDataframe(
                            df=subdf,
                            area=self.runinfo.area,
                            step_begin=self.last_run_step + 1,
                            step_end=mwxrun.get_it()+ 1
                        )

                        sp_name = mwxrun.simulation.species[sp].name
                        if (keytype, key, sp_name) in self.ts_dict:
                            self.ts_dict[(keytype, key, sp)] = (
                                timeseries.concat_crop_timeseries(
                                    [self.ts_dict[(keytype, key, sp)], ts]
                                )
                            )
                        else:
                            self.ts_dict[(keytype, key, sp_name)] = ts

    def update_fullhist_dict(self):
        """Once current diagnostic period is updated, update full history and
        resample if needed.
        """
        maxlen = 0

        # Build up strided, full-history timeseries.
        # Note fullkey refers to tuples of (keytype, key, species_name)
        for fullkey in self.ts_dict:
            if fullkey in self.fullhist_dict:
                self.fullhist_dict[fullkey] = (
                    timeseries.concat_crop_timeseries(
                        [self.fullhist_dict[fullkey],
                            self.ts_dict[fullkey]],
                        dt=self.history_dt
                    )
                )
            else:
                # We force all full histories to start at timestep 0. Otherwise
                # if eg anode absorption starts after first diagnostic period,
                # we'll have different times in denominators of different
                # values.
                self.fullhist_dict[fullkey] = timeseries.concat_crop_timeseries(
                    [self.ts_dict[fullkey]], step_begin=0, dt=self.history_dt)

            maxlen = max(
                maxlen,
                self.fullhist_dict[fullkey].step_end
                - self.fullhist_dict[fullkey].step_begin
            )

        # Resize full history if we're above the maximum length - this
        # resamples, including smoothing, to half the size. But only if the new
        # dt of the fullhist_dict will be less than a diagnostic interval.
        if maxlen > self.history_maxlen:
            if self.history_dt * 4. < mwxrun.get_dt() * self.diag_steps:
                self.history_dt *= 2.
                for val in list(self.fullhist_dict.values()):
                    val.resample(new_dt=self.history_dt, inplace=True)
            else:
                self.history_maxlen = maxlen

    def _check_charge_conservation(self):
        """Function to check net current flow into simulation during the last
        diagnostic period.
        """
        full_ts = timeseries.concat_crop_timeseries([
            val for key, val in self.ts_dict.items()])
        net_current = full_ts.get_averagevalue_by_key('J')
        if abs(net_current) > 0.5:
            if abs(net_current) > 1e3:
                raise RuntimeError(
                    ('Step %d: Net current exceeds 1000 A/cm^2, which is almost'
                     ' definitely an error.') % mwxrun.get_it()
                )
            warnings.warn(
                ('Step %d: Net current (%.3f A/cm^2) exceeds 0.5 A/cm^2, which '
                 'likely indicates a violation of the CFL condition.')
                % (mwxrun.get_it(), net_current)
            )


class FluxCalcDataframe(timeseries.Timeseries):

    """A FluxCalc created from a pandas Dataframe.
    Dataframes passed to this class should only consist of a single species; as
    a result multiple lines at the same timestep are not supported. Desired
    concatenation or summation across timeseries should be done on the
    processed timeseries, not on the raw dataframes.
    **INSERT TEXT ABOUT EXTRA COLUMNS HERE -- ESPECIALLY IMPLICATIONS FOR WHICH
    QUANTITIES MAKE SENSE BECAUSE THE NORMALIZATION BY WEIGHT MAKES SENSE FOR
    THEM**
    All calculations are done with reference to ground (V_e=0).
    """

    def __init__(self, df, area, step_begin=None, step_end=None, dt=None):
        """Initialize the FluxCalc object, transforming df into the appropriate
        array.

        Arguments:
            df (pandas.DataFrame): A single dataframe with columns:
            - t: Time in sec
            - step: Step number, integer, with constant dt necessary.
            - n: Number of macroparticles
            - q: Injected/deposited charge
            - E_part: Total KE + PE of the particles, relative to phi = 0,
                assuming potential energy = q*phi. For scraped data, this has
                already been multiplied by -1 in the scraper to denote energy
                leaving the system.
            - V_e: Fermi level of the emitting or collecting surface.
            - Optional extra columns
            area (float): Area to consider for calculating current densities,
                usually Lx*Ly from the simulation, in m^2. Values are
                transformed to /cm^2 internally.
            step_begin (int): Beginning timestep to consider. If not provided,
                take min step in df. Must be <= min step in df.
            step_end (int): Last timestep +1 to consider. If not provided,
                take max step in df +1. Must be > max step in df.
            dt (float): Timestep increment. If not provided, use warp.top.dt.
                Must provide in post-processing.
        """
        self.area = area
        self.fac = 1.0 / (1e4 * self.area)

        if step_begin is None:
            step_begin = df['step'].min()
        if step_end is None:
            step_end = df['step'].max() + 1

        if dt is None:
            dt = mwxrun.get_dt()

        special_keys = ['t', 'species_id', 'step', 'n', 'q', 'E_total', 'V_e']

        # we need to get the sum of collected weight for each species in order
        # to properly average extra arrays by the total number of particles
        spid_array = np.array((df['species_id']))

        sq_array = np.array([spec.sq for spec in mwxrun.simulation.species])

        q_array = np.array(df['q'])
        weight_sum = np.zeros_like(spid_array)
        idx = np.where(q_array != 0.0)
        weight_sum[idx] = np.abs(
            q_array[idx]
            / sq_array[spid_array[idx]]
        )
        idx = np.where(weight_sum == 0.0)
        weight_sum[idx] = np.array(df['n'])[idx]

        extra_keys = ['weight_sum']
        extra_arrays = [weight_sum]
        for key in df.keys():
            if key not in special_keys:
                extra_keys.append(key)
                extra_arrays.append(np.array(df[key]))
        num_extra_arrays=len(extra_arrays)

        # this is necessary or else numba doesn't play nice
        if num_extra_arrays == 0:
            extra_arrays = np.zeros((1, 1))
        else:
            extra_arrays = np.array(extra_arrays)

        timeseries_array = gen_timeseries(
            step_begin=step_begin, step_end=step_end,
            dt=dt,
            step_array=np.array(df['step'], dtype=int),
            n_array=np.array(df['n']),
            q_array=self.fac*q_array,
            E_array=self.fac*np.array(df['E_total']),
            V_e_array=np.array(df['V_e']),
            extra_arrays=extra_arrays,
            num_extra_arrays=num_extra_arrays
        )

        array_dict={
            'n': timeseries_array[:, 0],
            'J': timeseries_array[:, 1],
            'dQ': timeseries_array[:, 2],
            'P': timeseries_array[:, 3],
        }
        for ii, key in enumerate(extra_keys):
            array_dict[key] = timeseries_array[:,4+ii]

        super(FluxCalcDataframe, self).__init__(
            step_begin=step_begin, step_end=step_end, dt=dt,
            array_dict=array_dict
        )

class FluxDiagFromFile(FluxDiagBase):

    """Load a flux diag from its minimal saved file.
    This is meant for postprocessing ONLY; it will not initialize things
    properly for during-run operations.
    """

    def __init__(self, basedir='diags', fluxdatafile=None,
                 fluxdatafileformat='fluxes/fluxdata*', fs=None):
        """Load from fluxdata_XXXXXXXXXX.dpkl and runinfo.dpkl files.

        Arguments:
            basedir (str): Base directory of the diagnostic files. Has to
                contain a runinfo dill-pickled file.
            fluxdatafile (str): Filename dill-pickled with FluxDiagBase.save().
            fluxdatafileformat (str): If a specific fluxdatafile is not
                specified, this naming format will be used to search for the
                latest fluxdata file in the base directory.
            fs (s3fs filesystem): Optional S3 filesystem to load data directly
                from a S3 bucket.
        """
        if fs is None:
            self.open_command = open
            self.exists_command = os.path.exists
            self.glob_command = glob.glob
        else:
            self.open_command = fs.open
            self.exists_command = fs.exists
            self.glob_command = fs.glob

        self.fluxdatafileformat = fluxdatafileformat

        runinfofile = self._get_runinfofile(basedir)
        if fluxdatafile is None:
            fluxdatafile = self._get_fluxdatafile(runinfofile, fs)

        with self.open_command(runinfofile, 'rb') as pfile:
            self.runinfo = dill.load(pfile)

        self.component_list = self.runinfo.component_list
        self.species_list = [
            species.name for species in self.runinfo.species_list
        ]

        with self.open_command(fluxdatafile, 'rb') as pfile:
            dict_to_load = dill.load(pfile)

        self.__dict__.update(dict_to_load)

    def _get_runinfofile(self, basedir):
        """Function to check if a runinfo.dpkl file exists in the base
        directory."""
        runinfofile = os.path.join(basedir, 'runinfo.dpkl')
        if self.exists_command(runinfofile):
            return runinfofile
        raise IOError(f'No runinfo.dpkl found in base directory {basedir}')

    def _get_fluxdatafile(self, runinfofile, fs):
        """Function to look for latest fluxdata file from the base directory of
        the runinfo.dpkl file.
        """
        direc_base = runinfofile[:runinfofile.rfind('/')+1]
        flux_files = sorted(
            self.glob_command(os.path.join(direc_base, self.fluxdatafileformat))
        )
        if len(flux_files) == 0:
            raise IOError('No fluxdata files found in base directory: %s'
                          % (direc_base))
        return flux_files[-1]

# ### Numba functions for FluxCalcDataframe ###
@numba.jit(nopython=True)
def gen_timeseries(step_begin, step_end, dt, step_array, n_array,
                   q_array, E_array, V_e_array, extra_arrays, num_extra_arrays):
    """Transform dataframe entries, with potentially sparse rows (eg some
    timesteps have no associated row), into an array capturing all timesteps
    between step_begin and step_end. Multiple rows for the same timestep are
    not supported.

    Note:
        Signs are somewhat confusing here! For particles that don't have
        scattering or ionization, E_array is essentially KE_emit + q*V_emit. dQ
        for emission is then KE_emit + q*V_emit - q*V_e, but V_emit = V_e - WF,
        so dQ = KE_emit - q*WF. Since q is negative for electrons, we're adding
        |q|*WF potential energy into simulation, which does correctly represent
        surface cooling due to electron evaporation.
        dQ for absorption -- the arguments passed in, E_array and q_array, are
        already both negated to represent absorption vs injection. With
        standard physical signs for q, Q, KE, etc., then,
        dQ = -(KE_emit + q*V_cathode) + q*V_e_anode =
        -(KE_emit + q(V_cathode - (V_anode + WF))), where V_cathode and V_anode
        are vacuum biases. The front minus sign indicates absorption of this
        heat onto the surface/out of the system.  q*(V_cathode - V_anode)
        captures kinetic energy from accelerating the particle across the gap.
        -q*WF is the heat deposited from the electron falling back to the Fermi
        level. So it should all work, and tests confirm energy conservation --
        see test_diags_fluxdiag.py.

    Arguments:
        step_begin (int): First step to consider
        step_end (int): n + 1th step to consider (eg python indexing)
        dt (float): Timestep in seconds, necessary for normalization
        step_array (np.ndarray): m-length array of the step number for each row.
            All entries in step_array must be integers >= step_begin and <
            step_end. m is the number of rows in the original dataframe.
        n_array (np.ndarray): m-length array of the macroparticle count for each
            row.
        q_array (np.ndarray): m-length array of the charge (in C) for each row.
        E_array (np.ndarray): m-length array of the KE + PE for each row. For
            scraped data, this has already been multiplied by -1 in the scraper
            to denote energy leaving the system.
        V_e_array (np.ndarray): m-length array of the Fermi level voltage for
            each row.
        extra_arrays (np.ndarrays): m-length arrays that are averaged into a
            timeseries.
        num_extra_arrays (int): Number of extra arrays included.

    Returns:
        timeseries_array (np.ndarray): (step_end - step_begin)x(4
        + len(extra_arrays)) array with n, J, dQ, P at each timestep. Extra
        arrays are also made into a timeseries.
        Q=(E_part - q*V_e)/dt, and P=-(q*V_e)/dt (eg electrical power). Signs
        mean we get heat INTO system, and "standard" power production. (Note q
        is signed like Q: positive into system / negative out of system. That
        means electrons leaving the system give positive q.)
    """
    timeseries_array = np.zeros((step_end - step_begin, 4 + num_extra_arrays))
    for idx in range(len(step_array)):
        i_step = step_array[idx] - step_begin
        timeseries_array[i_step, 0] = n_array[idx]
        timeseries_array[i_step, 1] = q_array[idx]/dt
        timeseries_array[i_step, 2] = (
            E_array[idx] - q_array[idx]*V_e_array[idx]
        ) / dt
        timeseries_array[i_step, 3] = -q_array[idx]*V_e_array[idx]/dt

        for ii in range(0, num_extra_arrays):
            timeseries_array[i_step, 4+ii] = extra_arrays[ii][idx]

    return timeseries_array
