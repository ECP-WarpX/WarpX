import logging
import os

import matplotlib.pyplot as plt
import numpy as np
from pywarpx import callbacks

from mewarpx import assemblies
from mewarpx.diags_store import diag_base
from mewarpx.mwxrun import mwxrun
from mewarpx.utils_store import parallel_util

logger = logging.getLogger(__name__)


class BaseParticleHistDiag(diag_base.WarpXDiagnostic):

    """Base class for diagnostics that record a multi-dimensional histogram of
    particle properties.
    Child classes should define ``name`` (diagnostic name), ``linres`` (the
    resolution in each dimension) and ``domain`` (the boundaries of the binned
    region).
    """

    PHIST_DIAG_DIR = "histograms"

    def __init__(self, diag_steps, species_list=None, include_overflow=True,
        **kwargs):
        """Initialize diagnostic.

        Arguments:
            diag_steps (int): Number of steps between each output.
            species_list (list of species): Species for which histograms should
                be recorded. If None, all simulation species will be included.
            include_overflow (bool): If True overflow bins to -inf and +inf
                will be added. This is useful for quantities for which definite
                ranges are not known (such as velocity).
        """
        self.write_dir = os.path.join(self.DIAG_DIR, self.PHIST_DIAG_DIR)
        self.species_list = species_list
        if self.species_list is None:
            self.species_list = mwxrun.simulation.species
        self.include_overflow = include_overflow

        super(BaseParticleHistDiag, self).__init__(
            diag_steps=diag_steps, **kwargs)

        # Check that weight is not in the quantities list. This is needed to
        # get the "sample" dimensions correct since weight is not a
        # histogram dimension.
        if 'w' in self.quantities:
            self.quantities.remove('w')

        # install the scraper diag function
        self.assembly.add_diag_fn(
            self._process_scraped_particles, self.quantities + ['w']
        )

        # Create bins and data array after the simulation has been initialized
        # since the bin specifications are written to file
        callbacks.installafterinit(self.initialize)

        # counter to properly normalize data
        self.accumulated_steps = 0

    def initialize(self):
        """Set up the histogram bins and data array."""
        self.setup_bins()
        self.setup_array()

    def setup_bins(self):
        """Calculate the bin edges."""
        # sanity check
        if len(self.domain) != len(self.linres):
            raise AttributeError(
                "Number of domain ranges should match length of linres array."
            )

        self.bins = []
        for ii in range(len(self.domain)):
            if self.linres[ii] <= 2:
                self.bins.append([-np.inf, np.inf])
            else:
                if self.include_overflow:
                    bins = (
                        [-np.inf]
                        + list(np.linspace(self.domain[ii][0],
                               self.domain[ii][1], self.linres[ii]-1))
                        + [np.inf]
                    )
                else:
                    bins = list(
                        np.linspace(self.domain[ii][0], self.domain[ii][1],
                         self.linres[ii]+1)
                    )
                self.bins.append(bins)
        self.bins = np.array(self.bins, dtype=object)

        # save bin details to file
        if mwxrun.me == 0:
            fname = os.path.join(
                self.write_dir, f"{self.assembly.name}_histogram_bins.npy"
            )
            np.save(fname, self.bins)

    def setup_array(self):
        """Calculate the necessary array size; create the array itself.
        """
        shape = ((len(self.species_list),) + tuple(self.linres))
        self.Harray = np.zeros(shape)

    def get_sample(self, scraped_particle_dict, species_id):
        """Function to get appropriate "sample" for binning.

        Arguments:
            scraped_particle_dict (dict): Dictionary containing scraped
                particle properties.
            species_id (int): ID of the species for which a sample should be
                collected.
        """
        idxs = np.where(
            scraped_particle_dict['species_id'] == species_id
        )
        sample = np.zeros((np.size(idxs), len(self.quantities)))
        for ii, key in enumerate(self.quantities):
            sample[:,ii] = scraped_particle_dict[key][idxs]

        return sample, scraped_particle_dict['w'][idxs]

    def _process_scraped_particles(self, scraped_particle_dict):
        """Function to process the scraped particle data.

        Arguments:
            scraped_particle_dict (dict): Dictionary containing scraped
                particle data.
        """
        # metools.runtools.PHistDiag() implements very general histogram
        # functionality. Take a look at that class before implementing new
        # functionality since that will likely avoid the need to duplicate work.

        for ii, species in enumerate(self.species_list):
            sample, weights = self.get_sample(
                scraped_particle_dict, species.species_number
            )
            H, _ = np.histogramdd(
                sample=sample, bins=self.bins, weights=weights
            )
            self.Harray[ii,...] += H
        self.accumulated_steps += 1

        # if this is a diagnostic period save the data
        if self.check_timestep():
            self.save_and_reset_histogram()

    def save_and_reset_histogram(self):
        """Save and reset the histogram."""
        # sum the histograms from all the processors
        self.Harray = parallel_util.parallelsum(self.Harray)

        if mwxrun.me == 0:
            self.Harray[:] /= (self.accumulated_steps * mwxrun.dt)
            for ii, species in enumerate(self.species_list):
                fileprefix = self.get_fileprefix(species.name)
                np.save(fileprefix + '.npy', self.Harray[ii,...])

            if self.plot:
                self.plot_histograms()

        self.Harray[:] = 0.0
        self.accumulated_steps = 0

    def get_fileprefix(self, species):
        """Return filepath except for the filetype.

        Arguments:
            title (str): String for title; whitespace will become underscores
        """
        return os.path.join(
            self.write_dir, f"{self.name}_{species}_{mwxrun.get_it():010d}"
        )


class ZPlanePHistDiag(BaseParticleHistDiag):

    """Particle histogram diagnostic to track the location particles are
    scraped on a ZPlane assembly.
    """

    def __init__(self, diag_steps, assembly, species_list=None, linres_x=30,
                 linres_y=30, plot=True, **kwargs):
        """Initialize diagnostic.

        Arguments:
            diag_steps (int): Number of steps between each output.
            assembly (mewarpx.Assembly object): The assembly object in which
                particles are scraped.
            species_list (list of species): Species for which histograms should
                be recorded. If None, all simulation species will be included.
            linres_x (int): Number of bins to create along the x-dimension.
                Default 30.
            linres_y (int): Number of bins to create along the y-dimension.
                Default 30 unless the simulation is in 2D in which case only
                1 y-bin will be used regardless of this input parameter.
            plot (bool): If True plot the histograms at each diagnostic step.
            kwargs: See :class:`mewarpx.diags_store.diag_base.WarpXDiagnostic`
                for more timing options.
        """
        self.assembly = assembly
        self.plot = plot

        self.name = self.assembly.name

        # sanity check
        if not isinstance(self.assembly, assemblies.ZPlane):
            raise AttributeError(
                "ZPlanePHistDiag should only be used with a ZPlane assembly."
            )

        # set sensible linres based on simulation geometry - can override
        # user given arguments if they don't make sense
        if mwxrun.dim == 1:
            logger.warn(
                "ZPlane particle scraping histograms in 1d are trivial so this "
                "diagnostic will not be installed."
            )
            return
        elif mwxrun.dim == 2:
            linres_y = 1

        # set up domain and number of bins
        self.domain = [(mwxrun.xmin, mwxrun.xmax), (mwxrun.ymin, mwxrun.ymax)]
        self.linres = [linres_x, linres_y]

        # for now we hard code this diagnostic to only save weight but this
        # can be extended in the future
        self.quantities = ['x', 'y']

        super(ZPlanePHistDiag, self).__init__(
            diag_steps=diag_steps, species_list=species_list,
            include_overflow=False, **kwargs
        )

    def plot_histograms(self):
        """Function to plot histogram of scraped particle positions."""
        for ii, species in enumerate(self.species_list):
            # skip plotting if no data is present
            if np.all(self.Harray[ii,...] == 0):
                continue
            plt.title(f"{species.name} scraped on {self.assembly.name}")

            if mwxrun.dim == 2:
                plt.xlabel("x position (mm)")
                plt.ylabel("Current density (A/cm$^2$)")

                dx = (self.bins[0][1] - self.bins[0][0])
                data = (
                    self.Harray[ii,:,0] * abs(species.sq)
                    / (dx * (mwxrun.ymax - mwxrun.ymin))
                )
                plt.bar(
                    np.array((self.bins[0])[:-1] + dx/2)*1e3, data*1e-4,
                    width=dx*1e3, alpha=0.8
                )
            else:
                # 2d plotting functionality can be found in
                # metools.runtools.PHistDiag() if needed in the future.
                logger.warn("2D histogram plotting not yet implemented.")
                return

            fileprefix = self.get_fileprefix(species.name)
            plt.savefig(fileprefix + '.png', dpi=300)
