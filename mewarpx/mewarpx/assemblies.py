"""
Assembly implementations.
"""

import collections
import numpy as np

from mewarpx.mwxrun import mwxrun
from mewarpx.utils_store.appendablearray import AppendableArray
from mewarpx.utils_store import parallel_util
from pywarpx import picmi, _libwarpx


class Assembly(object):

    """An assembly represents any shape in the simulation; usually a conductor.

    While V, T, and WF are required, specific use cases may allow None to be
    used for these fields.
    """

    # This is overridden if a diagnostic is installed to record scraped
    # current.
    scraper_diag = None

    # fields is used by the diags.FluxInjectorDiag to know what to write to
    # the CSV file. It can be overridden by child classes, but is not currently
    # adjustable by the user.

    # IF CHANGING THIS, CHANGE IN self.record_scrapedparticles() AS WELL.
    fields = ['t', 'step', 'species_id', 'V_e', 'n', 'q', 'E_total']

    def __init__(self, V, T, WF, name):
        """Basic initialization.

        Arguments:
            V (float): Voltage (V)
            T (float): Temperature (K)
            WF (float): Work function (eV)
            name (str): Assembly name
        """
        self.V = V
        self.T = T
        self.WF = WF
        self.name = name

    def getvoltage(self):
        """Allows for time-dependent implementations to override this."""
        return mwxrun.eval_expression_t(self.V)

    def getvoltage_e(self):
        """Allows for time-dependent implementations to override this."""
        return self.getvoltage() + self.WF

    def init_scrapedparticles(self, fieldlist):
        """Set up the scraped particles array. Call before
        append_scrapedparticles.

        Arguments:
            fieldlist (list): List of string titles for the fields. Order is
                important; it must match the order for future particle appends
                that are made.
        """
        self._scrapedparticles_fields = fieldlist
        self._scrapedparticles_data = AppendableArray(
            typecode='d', unitshape=[len(fieldlist)])

    def record_scrapedparticles(self):
        """Handles transforming raw particle information from the WarpX
        scraped particle buffer to the information used to record particles as
        a function of time.

        Note:
            Assumes the fixed form of fields given in Assembly().  Doesn't
            check since this is called many times.

        Note:
            The total charge scraped and energy of the particles scraped
            is multiplied by -1 since these quantities are leaving the system.
        """

        # skip conductors that don't have a label to get scraped particles with
        if not hasattr(self, 'scraper_label'):
            return

        # loop over species and the scraped particle data from the buffer
        for species in mwxrun.simulation.species:
            data = np.zeros(7)
            empty = True
            idx_list = []

            # get the number of particles in the buffer - this is primarily
            # to avoid trying to access the buffer if it is not defined
            # which causes a segfault
            buffer_count = _libwarpx.get_particle_boundary_buffer_size(
                species.name, self.scraper_label
            )
            if buffer_count > 0:
                # get the timesteps at which particles were scraped
                comp_steps = _libwarpx.get_particle_boundary_buffer(
                    species.name, self.scraper_label, "step_scraped", mwxrun.lev
                )
                # get the particles that were scraped in this timestep
                for arr in comp_steps:
                    idx_list.append(np.where(arr == mwxrun.get_it())[0])
                    if len(idx_list[-1]) != 0:
                        empty = False

            data[0] = mwxrun.get_t()
            data[1] = mwxrun.get_it()
            data[2] = species.species_number
            data[3] = self.getvoltage_e()
            if not empty:
                data[4] = sum(np.size(idx) for idx in idx_list)
                w_arrays = _libwarpx.get_particle_boundary_buffer(
                    species.name, self.scraper_label, "w", mwxrun.lev
                )
                data[5] = -species.sq * sum(np.sum(w_arrays[i][idx_list[i]])
                     for i in range(len(idx_list))
                )
                E_arrays = _libwarpx.get_particle_boundary_buffer(
                    species.name, self.scraper_label, "E_total", mwxrun.lev
                )
                data[6] = -sum(np.sum(E_arrays[i][idx_list[i]])
                     for i in range(len(idx_list))
                )

            self.append_scrapedparticles(data)

    def append_scrapedparticles(self, data):
        """Append one or more lines of scraped particles data.

        Arguments:
            data (np.ndarray): Array of shape (m) or (n, m) where m is the
                number of fields and n is the number of rows of data to append.
        """
        self._scrapedparticles_data.append(data)

    def get_scrapedparticles(self, clear=False):
        """Retrieve a copy of scrapedparticles data.

        Arguments:
            clear (bool): If True, clear the particle data rows entered (field
                names are still initialized as before). Default False.

        Returns:
            scrapedparticles_dict (collections.OrderedDict): Keys are the
                originally passed field strings for lost particles. Values are
                an (n)-shape numpy array for each field.
        """
        lpdata = self._scrapedparticles_data.data()

        # Sum all except t/step/jsid/V_e from all processors
        lpdata[:,4:] = parallel_util.parallelsum(np.array(lpdata[:,4:]))

        lpdict = collections.OrderedDict(
            [(fieldname, np.array(lpdata[:, ii], copy=True))
             for ii, fieldname in enumerate(self._scrapedparticles_fields)])

        if clear:
            self._scrapedparticles_data.cleardata()

        return lpdict


class ZPlane(Assembly):

    """A semi-infinite plane."""

    def __init__(self, z, zsign, V, T, WF, name):
        """Basic initialization.

        Arguments:
            z (float): The edge of the semi-infinite plane (m)
            zsign (int): =1 to extend from z to +inf, or =-1 to extend from
                -inf to z.
            V (float): Voltage (V)
            T (float): Temperature (K)
            WF (float): Work function (eV)
            name (str): Assembly name
        """
        super(ZPlane, self).__init__(V=V, T=T, WF=WF, name=name)

        self.z = z

        self.zsign = int(round(zsign))
        if self.zsign not in [-1, 1]:
            raise ValueError(f"self.zsign = {self.zsign} is not either -1 or 1.")


class Cathode(ZPlane):
    """A basic wrapper to define a semi-infinite plane for the cathode."""

    def __init__(self, V, T, WF):
        self.scraper_label = 'z_lo'
        super(Cathode, self).__init__(
            z=0, zsign=-1, V=V, T=T, WF=WF, name="Cathode"
        )


class Anode(ZPlane):
    """A basic wrapper to define a semi-infinite plane for the anode."""

    def __init__(self, z, V, T, WF):
        self.scraper_label = 'z_hi'
        super(Anode, self).__init__(
            z=z, zsign=1, V=V, T=T, WF=WF, name="Anode"
        )


class Cylinder(Assembly):
    """An infinitely long Cylinder pointing in the y-direction."""

    def __init__(self, center_x, center_z, radius, V, T, WF, name,
                 install_in_simulation=True):
        """Basic initialization.

        Arguments:
            center_x (float): The x-coordinates of the center of the cylinder.
                Coordinates are in (m)
            center_z (float): The z-coordinates of the center of the cylinder.
                Coordinates are in (m)
            radius (float): The radius of the cylinder (m)
            V (float): Voltage (V)
            T (float): Temperature (K)
            WF (float): Work function (eV)
            name (str): Assembly name
            install_in_simulation (bool): If True and the Assembly is an
                embedded boundary it will be included in the WarpX simulation
        """
        super(Cylinder, self).__init__(V=V, T=T, WF=WF, name=name)
        self.center_x = center_x
        self.center_z = center_z
        self.radius = radius
        self.implicit_function = (
            f"-((x-{self.center_x})**2+(z-{self.center_z})**2-{self.radius}**2)"
        )

        self.scraper_label = 'eb'

        if install_in_simulation:
            self._install_in_simulation()

    def isinside(self, X, Y, Z, aura=0):
        """
        Calculates which grid tiles are within the cylinder.

        Arguments:
            X (np.ndarray): array of x coordinates of flattened grid.
            Y (np.ndarray): array of y coordinates of flattened grid.
            Z (np.ndarray): array of z coordinates of flattened grid.
            aura (float): extra space around the conductor that is considered inside. Useful
                for small, thin conductors that don't overlap any grid points. In
                units of meters.

        Returns:
            result (np.ndarray): array of flattened grid where all tiles
                inside the cylinder are 1, and other tiles are 0.
        """
        dist_to_center = (X - self.center_x)**2 + (Z - self.center_z)**2
        boundary = (self.radius + aura)**2
        result = np.where(dist_to_center <= boundary, 1, 0)

        return result

    def calculatenormal(self, px, py, pz):
        """
        Calculates Normal of particle inside/outside of conductor to nearest
        surface of conductor

        Arguments:
            px (np.ndarray): x-coordinate of particle in meters.
            py (np.ndarray): y-coordinate of particle in meters.
            pz (np.ndarray): z-coordinate of particle in meters.

        Returns:
            nhat (np.ndarray): Normal of particles of conductor to nearest
                surface of conductor.
        """
        # distance from center of cylinder
        dist = np.sqrt((px - self.center_x)**2 + (pz - self.center_z)**2)
        nhat = np.zeros([3, len(px)])
        nhat[0, :] = (px - self.center_x) / dist
        nhat[2, :] = (pz - self.center_z) / dist
        return nhat

    def _install_in_simulation(self):
        """Function to pass this EB object to the WarpX simulation."""

        if mwxrun.simulation.embedded_boundary is not None:
            raise RuntimeError('Currently only 1 EB is supported.')

        mwxrun.simulation.embedded_boundary = picmi.EmbeddedBoundary(
            implicit_function=self.implicit_function,
            potential=self.V
        )
