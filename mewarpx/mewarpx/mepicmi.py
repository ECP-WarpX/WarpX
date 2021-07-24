"""
Child classes of WarpX Picmi objects to support ME utility functions and
extended functionality without going through their standards.
"""
import logging
import numpy as np

from pywarpx import picmi, _libwarpx
from mewarpx.mwxrun import mwxrun

# Get module-level logger
logger = logging.getLogger(__name__)


class Species(picmi.Species):

    def init(self, kw):
        """PICMI relies on C++ WarpX to translate charge & mass strings to
        floats. To get around that we have our own variables sq/sm (species
        charge/mass) that are always floats.

        Automatically adds the species to the simulation
        """
        self.grid = kw.pop("warpx_grid", None)
        self.n_macroparticle_per_cell = kw.pop(
            "warpx_n_macroparticle_per_cell", [0, 0])

        super(Species, self).init(kw)

        if isinstance(self.charge, str):
            if self.charge == 'q_e':
                self.sq = picmi.constants.q_e
            elif self.charge == '-q_e':
                self.sq = -picmi.constants.q_e
            else:
                raise ValueError("Unrecognized charge {}".format(self.charge))
        else:
            self.sq = self.charge

        if isinstance(self.mass, str):
            if self.mass == 'm_e':
                self.sm = picmi.constants.m_e
            elif self.mass == 'm_p':
                self.sm = picmi.constants.m_p
            else:
                raise ValueError("Unrecognized mass {}".format(self.mass))
        else:
            self.sm = self.mass

        mwxrun.simulation.add_species(
            self,
            layout=picmi.GriddedLayout(
                n_macroparticle_per_cell=self.n_macroparticle_per_cell,
                grid=self.grid
            )
        )

    def init_pid_dict(self):
        """Function to build a dictionary of all the extra particle attributes
        (and weight). This has to happen after warpx has been initialized
        so that the particle container already exists.
        """
        self.pid_dict = {'w':0}
        self.nattribs = _libwarpx.get_nattr()
        # npids are used to inform the dimension of the attr array (used in
        # add_particles), so we need to subtract the non-extra pid arrays
        # (which include the weight) and then add 1 since the weight is treated
        # as an extra pid by AddNParticles()
        self.npids = _libwarpx.get_nattr_species(self.name) - self.nattribs + 1

    def add_pid(self, pid_name):
        """Wrapper to add a new PID to the particle data arrays at runtime.

        Arguments:
            pid_name (str): Name of the new PID.
        """
        # check if PID already exists
        if pid_name in self.pid_dict:
            return

        _libwarpx.add_real_comp(self.name, pid_name)
        # the attr array index (used by add_particles) is the component
        # index - the PIdx::nattribs (the non-extra pid arrays) + 1 (since
        # in AddNParticles the weight is treated as an extra array)
        self.pid_dict[pid_name] = _libwarpx.get_particle_comp_index(
            self.name, pid_name
        ) - (self.nattribs + 3) + 1
        self.npids += 1

    def get_array_from_pid(self, pid_name, level=0):
        """Wrapper to grab particle data for a specific PID.

        Arguments:
            pid_name (str): Name of the PID to grab data for.
            level (int): Level for which to grab the data.

        Returns:
            A list of numpy arrays. The list has one element for every tile
            with the numpy array holding the particle data for the requested
            PID.
        """
        return _libwarpx.get_particle_arrays(self.name, pid_name, level)

    def add_particles(self, x, y, z, ux, uy, uz, w, unique_particles=True,
                      **kwargs):
        """Wrapper for libwarpx function `add_particles`. It is specifically
        a function of the species class so that it can process and properly
        order the data for each pid.
        """

        pid_array = np.zeros((np.size(x), self.npids))
        pid_array[:,self.pid_dict['w']] = w
        for key, val in kwargs.items():
            pid_array[:,self.pid_dict[key]] = val

        _libwarpx.add_particles(
            self.name, x=x, y=y, z=z, ux=ux, uy=uy, uz=uz, attr=pid_array,
            unique_particles=unique_particles
        )
