"""
Assembly implementations.
"""

from mewarpx.mwxrun import mwxrun
from pywarpx import picmi
import numpy as np

class Assembly(object):

    """An assembly represents any shape in the simulation; usually a conductor.

    While V, T, and WF are required, specific use cases may allow None to be
    used for these fields.
    """

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
        return self.V


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
        super(Cathode, self).__init__(
            z=0, zsign=-1, V=V, T=T, WF=WF, name="Cathode"
        )


class Anode(ZPlane):
    """A basic wrapper to define a semi-infinite plane for the anode."""

    def __init__(self, z, V, T, WF):
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

        if install_in_simulation:
            self._install_in_simulation()

    def isinside(self, X, Y, Z, aura):
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
