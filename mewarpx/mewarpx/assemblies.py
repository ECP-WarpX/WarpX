"""
Assembly implementations.
"""
import collections
import logging

import numpy as np
from pywarpx import callbacks, picmi

from mewarpx.mwxrun import mwxrun
from mewarpx.utils_store.appendablearray import AppendableArray

# Get module-level logger
logger = logging.getLogger(__name__)

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

    def __init__(self, V, T, WF, name, read_scraped_particles):
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
        self.read_scraped_particles = read_scraped_particles

        # list of diagnostic functions to call after collecting scraped
        # particles
        self.scraper_diag_fnlist = []
        # the step scraped will be recorded by default
        self.scraped_particle_attribs_list = ["step_scraped"]
        # if needed will be an appendable array to store scraped particle
        # properties for easy processing with diagnostic functions
        self.scraped_particle_array = None

        # currently the beforeEsolve callback is the most immediate after
        # scraping callback
        if self.read_scraped_particles:
            callbacks.installbeforeEsolve(self._run_scraper_diag_steps)

    def check_geom(self):
        """Throw an error if the simulation geometry is unsupported by the
        Assembly.
        """
        if mwxrun.geom_str not in self.geoms:
            raise ValueError(
                f"{geom} geometry not supported by this Assembly")

    def getvoltage(self):
        """Allows for time-dependent implementations to override this."""
        return mwxrun.eval_expression_t(self.V)

    def getvoltage_e(self):
        """Allows for time-dependent implementations to override this."""
        return self.getvoltage() + self.WF

    def _install_in_simulation(self):
        """Function to handle installation of embedded boundaries in the
        simulation.
        Care is taken to have the possibility to install multiple embedded
        conductors (EBs). This is done by making the overall implicit function
        that describe the EBs use the maximum of the individual EB's implicit
        functions - this works since for each individual implicit function -1
        means outside the boundary, 0 means on the boundary edge and 1 means
        inside the boundary.
        Furthermore, the potential for the embedded boundary is also constructed
        with support for multiple embedded boundaries by adding a check for
        which embedded boundary is considered when assigning the potential.
        """
        if mwxrun.simulation.embedded_boundary is not None:

            # wrap the current EB potential in an if statement with the
            # new implicit function
            old_str = mwxrun.simulation.embedded_boundary.potential.strip('"')
            mwxrun.simulation.embedded_boundary.potential = (
                f'"if({self.implicit_function}>0,{self.V},{old_str})"'
            )

            # add a nested max() statement with the current implicit function
            # to add the new embedded boundary
            old_str = mwxrun.simulation.embedded_boundary.implicit_function.strip('"')
            mwxrun.simulation.embedded_boundary.implicit_function = (
                f'"max({old_str},{self.implicit_function})"'
            )
        else:
            mwxrun.simulation.embedded_boundary = picmi.EmbeddedBoundary(
                implicit_function=f'"{self.implicit_function}"',
                potential=f'"{self.V}"'
            )

    def add_diag_fn(self, diag_fn, attrib_list):
        """Add a diagnostic function to the diagnostic function list.

        Arguments:
            diag_fn (function): Function to be run after
                _read_scraped_particles() that processes the scraped particle
                dictionary for a specific diagnostic purpose.
            attrib_list: (list of strings): Particle attributes that should
                be included in the scraped particle dictionary for the given
                diagnostic function.
        """
        if not callable(diag_fn):
            raise ValueError("Must provide a callable diagnostic function.")

        self.scraper_diag_fnlist.append(diag_fn)

        new_attribs = False
        for attrib in attrib_list:
            if attrib not in self.scraped_particle_attribs_list:
                self.scraped_particle_attribs_list.append(attrib)
                new_attribs = True

        # re-initialize the scraped_particle_array so that it has space for
        # the new particle attributes
        if new_attribs:
            if (self.scraped_particle_array is not None
                and self.scraped_particle_array._datalen != 0
            ):
                raise RuntimeError(
                    f"The scraped particle array of {self.name} is not empty; "
                    "cannot add new particle attribute to track."
                )
            # note the +1 in unitshape is to save species ID
            self.scraped_particle_array = AppendableArray(
                typecode='d',
                unitshape=[len(self.scraped_particle_attribs_list)+1]
            )

    def _run_scraper_diag_steps(self):
        """Perform the steps involved in diagnostics of scraped particles."""
        # return early if no particle attributes will be recorded
        if self.scraped_particle_array is None:
            return

        if not hasattr(self, 'scraper_label'):
            raise AttributeError(
                f"Cannot record particles scraped on {self.name} since it "
                "does not have a scraper label."
            )

        # accumulate the scraped particle data
        self._read_scraped_particles()
        # construct dictionary of scraped particle data
        lpdict = self._get_scraped_particles_dict()
        # run all diagnostic functions
        for func in self.scraper_diag_fnlist:
            func(lpdict)

    def _read_scraped_particles(self):
        """Function to read the scraped particle buffer and populate the
        scraped_particle_array."""
        # loop over species and get the scraped particle data from the buffer
        for species in mwxrun.simulation.species:
            idx_list = []

            # get the number of particles in the buffer - this is primarily
            # to avoid trying to access the buffer if it is not defined
            # which causes a segfault
            buffer_count = mwxrun.sim_ext.get_particle_boundary_buffer_size(
                species.name, self.scraper_label
            )
            # logger.info(f"{self.name} scraped {buffer_count} {species.name}")

            # if there are no particles in the buffer continue to next species
            if buffer_count == 0:
                continue

            # get the timesteps at which particles were scraped
            comp_steps = mwxrun.sim_ext.get_particle_boundary_buffer(
                species.name, self.scraper_label, "step_scraped", mwxrun.lev
            )
            # get the particles that were scraped in this timestep
            for arr in comp_steps:
                idx_list.append(np.where(arr == mwxrun.get_it())[0])

            raw_particle_data = {}

            # get the particle structs from the boundary buffer
            structs = mwxrun.sim_ext.get_particle_boundary_buffer_structs(
                species.name, self.scraper_label, mwxrun.lev
            )
            if mwxrun.geom_str == 'XZ' or mwxrun.geom_str == 'RZ':
                raw_particle_data['x'] = [struct['x'] for struct in structs]
                raw_particle_data['y'] = [
                    np.zeros(len(struct['y'])) for struct in structs]
                raw_particle_data['z'] = [struct['y'] for struct in structs]
            elif mwxrun.geom_str == 'XYZ':
                raw_particle_data['x'] = [struct['x'] for struct in structs]
                raw_particle_data['y'] = [struct['y'] for struct in structs]
                raw_particle_data['z'] = [struct['z'] for struct in structs]
            else:
                raise NotImplementedError(
                    f"Scraping not implemented for {mwxrun.geom_str}."
                )

            # sort the particles appropriately if this is an eb
            if self.scraper_label == 'eb':
                temp_idx_list = []

                for i in range(len(idx_list)):
                    is_inside = self.isinside(
                        raw_particle_data['x'][i][idx_list[i]],
                        raw_particle_data['y'][i][idx_list[i]],
                        raw_particle_data['z'][i][idx_list[i]]
                    )
                    temp_idx_list.append(idx_list[i][np.where(is_inside)])

                    # set the scraped timestep to -1 for particles in this
                    # EB so that they are not considered again
                    comp_steps[i][idx_list[i][np.where(is_inside)]] = -1

                idx_list = temp_idx_list

            for ii, attrib in enumerate(self.scraped_particle_attribs_list):
                # get the desired particle property
                if attrib not in ['x', 'y', 'z']:
                    raw_particle_data[attrib] = (
                        mwxrun.sim_ext.get_particle_boundary_buffer(
                            species.name, self.scraper_label, attrib, mwxrun.lev
                        )
                    )

            # loop over all tiles and append the particle data from that tile
            # to scraped_particle_data
            for ii, idxs in enumerate(idx_list):
                if len(idxs) == 0:
                    continue
                data = np.zeros(
                    (len(idxs), len(self.scraped_particle_attribs_list)+1)
                )
                data[:,0] = species.species_number
                for jj in range(0, len(self.scraped_particle_attribs_list)):
                    data[:,jj+1] = raw_particle_data[
                        self.scraped_particle_attribs_list[jj]][ii][idxs]
                self.scraped_particle_array.append(data)

    def _get_scraped_particles_dict(self):
        """Retrieve a dictionary containing the scraped particle data.

        Returns:
            scrapedparticles_dict (collections.OrderedDict): Keys are the
                originally passed field strings for lost particles. Values are
                an (n)-shape numpy array for each field.
        """
        lpdata = self.scraped_particle_array.data()

        lpdict = collections.OrderedDict(
            [(name, np.array(lpdata[:, ii], copy=True)) for ii, name in
            enumerate(["species_id"]+self.scraped_particle_attribs_list)]
        )

        self.scraped_particle_array.cleardata()
        return lpdict


class ZPlane(Assembly):

    """A semi-infinite plane."""
    geoms = ['RZ', 'X', 'XZ', 'XYZ']

    def __init__(self, z, zsign, V, T, WF, name, read_scraped_particles=True):
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
        super(ZPlane, self).__init__(
            V=V, T=T, WF=WF, name=name,
            read_scraped_particles=read_scraped_particles
        )

        self.z = z

        self.zsign = int(round(zsign))
        if self.zsign not in [-1, 1]:
            raise ValueError(f"self.zsign = {self.zsign} is not -1 or 1.")


class Cathode(ZPlane):
    """A basic wrapper to define a semi-infinite plane for the cathode."""

    def __init__(self, V, T, WF, read_scraped_particles=True):
        super(Cathode, self).__init__(
            z=0, zsign=-1, V=V, T=T, WF=WF, name="Cathode",
            read_scraped_particles=read_scraped_particles
        )
        # set solver boundary potential
        mwxrun.grid.potential_zmin = self.V
        # set label to correctly grab scraped particles from the buffer
        self.scraper_label = 'z_lo'


class Anode(ZPlane):
    """A basic wrapper to define a semi-infinite plane for the anode."""

    def __init__(self, z, V, T, WF, read_scraped_particles=True):
        super(Anode, self).__init__(
            z=z, zsign=1, V=V, T=T, WF=WF, name="Anode",
            read_scraped_particles=read_scraped_particles
        )
        # set solver boundary potential
        mwxrun.grid.potential_zmax = self.V
        # set label to correctly grab scraped particles from the buffer
        self.scraper_label = 'z_hi'


class InfCylinderY(Assembly):
    """An infinitely long Cylinder pointing in the y-direction."""
    geoms = ['XZ', 'XYZ']

    def __init__(self, center_x, center_z, radius, V, T, WF, name,
                 install_in_simulation=True, read_scraped_particles=True):
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
        super(InfCylinderY, self).__init__(
            V=V, T=T, WF=WF, name=name,
            read_scraped_particles=read_scraped_particles
        )
        self.center_x = center_x
        self.center_z = center_z
        self.radius = radius
        # negative means outside of object
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
            aura (float): extra space around the conductor that is considered
                inside. Useful for small, thin conductors that don't overlap any
                grid points. In units of meters.

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


class CylinderZ(Assembly):
    """A Cylinder pointing in the z-direction with center along the x and y
    origin line (to support RZ)."""
    geoms = ['RZ']

    def __init__(self, V, T, WF, name, r_outer, r_inner=0.0, zmin=None,
                 zmax=None, read_scraped_particles=True):
        """Basic initialization.

        Arguments:
            V (float): Voltage (V)
            T (float): Temperature (K)
            WF (float): Work function (eV)
            name (str): Assembly name
            r_outer (float): The outer radius of the cylinder (m)
            r_inner (float): The inner radius of the cylinder (m)
            zmin (float): Lower z limit of the cylinder (m). Defaults to
                mwxrun.zmin.
            zmax (float): Upper z limit of the cylinder (m). Defaults to
                mwxrun.zmax.
        """
        super(CylinderZ, self).__init__(
            V=V, T=T, WF=WF, name=name,
            read_scraped_particles=read_scraped_particles
        )
        self.r_inner = r_inner
        self.r_outer = r_outer
        self.r_center = (self.r_inner + self.r_outer) / 2.0

        self.zmin = zmin
        if self.zmin is None:
            self.zmin = mwxrun.zmin
        self.zmax = zmax
        if self.zmax is None:
            self.zmax = mwxrun.zmax

        # check if z-limits span the full simulation domain
        # if np.close(self.zmax, mwxrun.zmax) and np.close(self.zmin, mwxrun.zmin):
        infinite_cyl = np.allclose(
            [self.zmax - mwxrun.zmax, self.zmin - mwxrun.zmin], 0.)

        self.center_x = 0.0
        self.center_y = 0.0

        # sanity check
        if self.r_outer == 0:
            raise AttributeError('Cannot have a cylinder with 0 outer radius.')

        # check if the cylinder forms the simulation boundary
        if np.isclose(self.r_inner, mwxrun.rmax):
            if not infinite_cyl:
                raise AttributeError(
                    "CylinderZ on the simulation boundary must span the full "
                    "z domain."
                )
            # set solver boundary potential
            mwxrun.grid.potential_xmax = self.V
            # set label to correctly grab scraped particles from the buffer
            self.scraper_label = 'x_hi'
        elif np.isclose(self.r_outer, mwxrun.rmin):
            if not infinite_cyl:
                raise AttributeError(
                    "CylinderZ on the simulation boundary must span the full "
                    "z domain."
                )
            # set solver boundary potential
            mwxrun.grid.potential_xmin = self.V
            # set label to correctly grab scraped particles from the buffer
            self.scraper_label = 'x_lo'
        else:
            # handle the EB case
            # a negative value of the implicit function means outside of object
            if infinite_cyl:
                # no z dependence since the cylinder is infinitely long
                self.implicit_function = (
                    f"-(abs(x-{self.r_center})-{self.r_outer - self.r_center})"
                )
            else:
                self.implicit_function = (
                    f"if(z>{self.zmin} and z<{self.zmax}, "
                    f"-(abs(x-{self.r_center})-{self.r_outer - self.r_center}),"
                    "-1.0)"
                )
            # set label to correctly grab scraped particles from the buffer
            self.scraper_label = 'eb'
            # install the EB in the simulation
            self._install_in_simulation()

    def isinside(self, X, Y, Z, aura=0):
        """
        Determines whether the given coordinates are inside the assembly.

        Arguments:
            X (np.ndarray): array of x coordinates.
            Y (np.ndarray): array of y coordinates.
            Z (np.ndarray): array of z coordinates.
            aura (float): extra space around the conductor that is considered
                inside. Useful for small, thin conductors that don't overlap any
                grid points. In units of meters.

        Returns:
            result (np.ndarray): array of values corresponding to the input
                coordinates where all points inside the assembly are 1, and
                others are 0.
        """
        dr = self.r_outer - self.r_inner + 2.0*aura
        dz = self.zmax - self.zmin + 2.0*aura

        # checks if (r_inner - aura <= R <= r_outer + aura) AND
        # (zmin - aura <= Z <= z_max + aura), with R the distance to the
        # cylinder axis
        result = np.where(
            (abs(np.sqrt((X - self.center_x)**2 + (Y - self.center_y)**2)
            - (self.r_inner + self.r_outer) / 2.0) <= dr / 2.0)
            & (abs(Z - (self.zmin + self.zmax) / 2.0) <= dz / 2.0),
            1, 0
        )
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
        dist = np.sqrt((px - self.center_x)**2 + (py - self.center_y)**2)
        nhat = np.zeros([3, len(px)])
        nhat[0, :] = (px - self.center_x) / dist
        nhat[2, :] = (pz - self.center_z) / dist

        # Flip the normal vector for points inside the cylinder inner wall,
        # or inside the cylinder assembly but closer to the outer wall than the inner.
        idx = np.where(
            np.logical_or(
                (dist < self.r_inner),
                ((dist > self.r_center) & (dist < self.router))
            )
        )
        nhat[:,idx] *= -1

        return nhat


class Rectangle(Assembly):
    """A rectangular prism infinite in y."""
    geoms = ['XZ', 'XYZ']

    def __init__(self, center_x, center_z, length_x, length_z, V, T, WF, name,
                 install_in_simulation=True, read_scraped_particles=True):
        """Basic initialization.

        Arguments:
            center_x (float): The x-coordinates of the center of the rectangle.
                Coordinates are in (m)
            center_z (float): The z-coordinates of the center of the rectangle.
                Coordinates are in (m)
            length_x (float): The length the rectangle in the x direction (m)
            length_z (float): The length the rectangle in the z direction (m)
            V (float): Voltage (V)
            T (float): Temperature (K)
            WF (float): Work function (eV)
            name (str): Assembly name
            install_in_simulation (bool): If True and the Assembly is an
                embedded boundary it will be included in the WarpX simulation
        """
        super(Rectangle, self).__init__(
            V=V, T=T, WF=WF, name=name,
            read_scraped_particles=read_scraped_particles
        )
        self.center_x = float(center_x)
        self.center_z = float(center_z)
        self.length_x = float(length_x)
        self.length_z = float(length_z)
        self.xmin = float(center_x - length_x/2)
        self.xmax = float(center_x + length_x/2)
        self.zmin = float(center_z - length_z/2)
        self.zmax = float(center_z + length_z/2)

        self.scaled_h = (
            2.0 * max(self.length_x, self.length_z)
            / min(self.length_x, self.length_z)
        )

        # the regions change depending on whether x or z is longer, so we need
        # to adjust the preset normals for these regions depending on the x and
        # z lengths
        if self.length_x <= self.length_z:
            self.region_normals = np.array([
                (0, -1), (1, 0), (0, 1), (-1, 0)
            ])
        else:
            self.region_normals = np.array([
                (-1, 0), (0, -1), (1, 0), (0, 1)
            ])

        self.transformation_matrix = self._get_transformation_matrix()

        # Negative means outside of object
        self.implicit_function = (
            f"-max(max(x-({self.xmax}),({self.xmin})-x),"
            f"max(z-({self.zmax}),({self.zmin})-z))"
        )
        self.scraper_label = 'eb'

        if install_in_simulation:
            self._install_in_simulation()

    def calculatenormal(self, px, py, pz):
        """
        Calculates Normal of particle inside/outside of conductor to nearest
        surface of conductor. Nearest surface of conductor is determined by the
        region of the transformed rectangle that the transformed particle
        belongs to.

        The 4 regions in order from 0 - 3 are:
            0: bottom portion of transformed rectangle and below transformed
               rectangle (-y)
            1: right portion of transformed rectangle and to the right of
               transformed rectangle (+x)
            2: top portion of transformed rectangle and above transformed
               rectangle (+y)
            3: left portion of transformed rectangle and to the left of
               transformed rectangle (-x)

        Arguments:
            px (np.ndarray): x-coordinate of particle in meters.
            py (np.ndarray): y-coordinate of particle in meters.
            pz (np.ndarray): z-coordinate of particle in meters.

        Returns:
            nhat (np.ndarray): Normal of particles of conductor to nearest
                surface of conductor.
        """
        arr_length = px.shape[0]
        points = np.zeros((arr_length, 2))
        points[:, 0] = px
        points[:, 1] = pz
        transformed_points = self._transform_coordinates(points)
        idx = self._sort_particles(transformed_points)
        normals = self.region_normals[idx]
        nhat = np.zeros((3, (normals.shape[0])))
        nhat[0, :] = normals[:, 0]
        nhat[2, :] = normals[:, 1]
        return nhat

    def _get_transformation_matrix(self):
        """Tranform1 is used to move the rectangle to the appropriate location
        in the plane and scale it to the appropriate size. The rectangle is
        placed such that the longest side is along the z-axis and the bottom
        right corner is at (1, -1). This placement makes it simple to pick out
        positions that are in regions 0 or 1 as well as positions in region 4
        or 5.
        """
        # first transpose the box to the origin
        transpose1 = np.identity(3)
        transpose1[0, 2] = (-self.center_x + self.length_x / 2.0)
        transpose1[1, 2] = (-self.center_z + self.length_z / 2.0)

        # next scale the axes to make zone 0 equilateral with side length 1
        l = min(self.length_x, self.length_z)
        scale = np.identity(3)
        scale[0, 0] = 2.0 / l
        scale[1, 1] = 2.0 / l

        # next transpose again so the center is at the origin
        transpose2 = np.identity(3)
        transpose2[0, 2] = -1
        transpose2[1, 2] = -1

        # rotate by 90 deg if the length_z is smallest
        rotate = np.identity(3)
        if self.length_z < self.length_x:
            rotate[0,0] = np.cos(np.pi/2)
            rotate[0,1] = -np.sin(np.pi/2)
            rotate[1,0] = np.sin(np.pi/2)
            rotate[1,1] = np.cos(np.pi/2)

        transformation_matrix = np.dot(
            np.dot(rotate, transpose2),
            np.dot(scale, transpose1)
        )
        return transformation_matrix

    def _transform_coordinates(self, points_in):
        """Transforms the coordinates of the particles using the
        transformation matrix of the rectangle"""
        points = np.ones((len(points_in), 3))
        points[:, 0] = points_in[:, 0]
        points[:, 1] = points_in[:, 1]

        points = points.dot(self.transformation_matrix.T)

        return points[:, 0:2]

    def _sort_particles(self, points):
        """Finds the regions of the rectangle that each point belongs to"""
        # The checks below isolate the region in which the given position lies.
        # c0 and c2 are True iff a particle is in region 0 or 2, respectively
        # not_c1 indicates a particle may be in regions 0, 2, or 3.
        c0 = np.array(points[:, 1] < -abs(points[:, 0]), dtype=int)
        c2 = (
            np.array(points[:, 1] - self.scaled_h + 2.0
            > abs(points[:, 0]), dtype=int)
        )
        not_c1 = np.array(points[:, 0] <= 0, dtype=int)
        # the booleans below form bit values to build up the index of the
        # region in which a given position lies.
        # bit2 is True if in region 2 or 3.
        bit2 = c2 | (not_c1 & np.logical_not(c0))
        # bit1 is True if in region 1 or 3; equivalently, if it's not in 0 or 2.
        bit1 = np.logical_not(c0 | c2)
        idx = 2*bit2 + bit1

        return idx

    def isinside(self, X, Y, Z, aura=0):
        """
        Calculates which grid tiles are within the rectangle.

        Arguments:
            X (np.ndarray): array of x coordinates of flattened grid.
            Y (np.ndarray): array of y coordinates of flattened grid.
            Z (np.ndarray): array of z coordinates of flattened grid.
            aura (float): extra space around the conductor that is considered
                inside. Useful for small, thin conductors that don't overlap any
                grid points. In units of meters.

        Returns:
            result (np.ndarray): array of flattened grid where all tiles
                inside the rectangle are 1, and other tiles are 0.
        """

        X_in_bounds = np.where(
            np.maximum(X-self.xmax-aura, self.xmin-aura-X) <= 0, 1, 0)
        Z_in_bounds = np.where(
            np.maximum(Z-self.zmax-aura, self.zmin-aura-Z) <= 0, 1, 0)
        result = np.where(X_in_bounds + Z_in_bounds > 1, 1, 0)

        return result
