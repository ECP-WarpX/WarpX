"""
Child classes of WarpX Picmi objects to support ME utility functions and
extended functionality without going through their standards.
"""
import logging
import numpy as np

from pywarpx import picmi, _libwarpx, callbacks
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
                n_macroparticle_per_cell=[0, 0],
                grid=mwxrun.grid
            )
        )
        self.pids_initialized = False
        # Only keys are used; ensures pids are both unique and ordered.
        self.waiting_extra_pids = {}
        # add a callback to initialize the extra PIDs after sim init
        callbacks.installafterinit(self.init_pid_dict)

    def init_pid_dict(self):
        """Function to build a list of all the extra particle attributes
        (and weight). This has to happen after warpx has been initialized
        so that the particle container already exists.
        """
        self.pids_initialized = True

        self.pid_list = ['w']
        self.nattribs = _libwarpx.get_nattr()
        # npids are used to inform the dimension of the attr array (used in
        # add_particles), so we need to subtract the non-extra pid arrays
        # (which include the weight) and then add 1 since the weight is treated
        # as an extra pid by AddNParticles()
        self.npids = _libwarpx.get_nattr_species(self.name) - self.nattribs + 1

        for pid in self.waiting_extra_pids.keys():
            self.add_pid(pid)
        del self.waiting_extra_pids

    def add_pid(self, pid_name):
        """Wrapper to add a new PID to the particle data arrays at runtime. The
        benefit of this wrapper is that we first check if the PID already
        exists so there is no risk in calling it multiple times for the same
        PID.

        Arguments:
            pid_name (str): Name of the new PID.
        """
        if not self.pids_initialized:
            self.waiting_extra_pids[pid_name] = None
            return

        # check if PID already exists
        if pid_name in self.pid_list:
            return

        _libwarpx.add_real_comp(self.name, pid_name)
        self.pid_list.append(pid_name)
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
