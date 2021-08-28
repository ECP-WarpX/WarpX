"""
This facilitates a central place to grab lots of simulation information for our
use.

SETUP order:
    - Before importing other mewarpx modules::

        from mewarpx import util as mwxutil
        mwxutil.init_libwarpx(ndim=ndim, rz=rz)

    - Import other mewarpx modules
    - Use ``from mewarpx.mwxrun import mwxrun`` to import the class holding
      all the simulation information, defined here.
    - Set up PICMI things, referring to mwxrun's ``picmi.Simulation`` object,
      ``mwxrun.simulation``. Add PICMI species to the PICMI Simulation.

    - Call this class's init_run::

        mwxrun.mwxrun.init_run()

    - Initialize any other mewarpx objects
    - Perform run with ``sim.step()``
"""
import ctypes
import numpy as np
import logging

from pywarpx import _libwarpx, picmi
from mewarpx.utils_store import mwxconstants as constants
from mewarpx.utils_store import parallel_util
from mewarpx.utils_store.util import compute_step

logger = logging.getLogger(__name__)


class MEWarpXRun(object):

    """Hold base information and wrappers for the simulation. This should be a
    singleton object.

    This also does the initialization steps, using init_run to enforce an
    initialization order.
    """

    def __init__(self):
        self.initialized = False
        self.simulation = picmi.Simulation(verbose=0)

    def init_grid(self, xmin, xmax, zmin, zmax, nx, nz, max_grid_size=4):
        """Function to set up the simulation grid."""
        self.grid = picmi.Cartesian2DGrid(
            number_of_cells=[nx, nz],
            lower_bound=[xmin, zmin],
            upper_bound=[xmax, zmax],
            lower_boundary_conditions=['periodic', 'dirichlet'],
            upper_boundary_conditions=['periodic', 'dirichlet'],
            lower_boundary_conditions_particles=['periodic', 'absorbing'],
            upper_boundary_conditions_particles=['periodic', 'absorbing'],
            warpx_max_grid_size=nz//max_grid_size
        )

        self._set_geom_str()
        self._set_grid_params()

    def init_timestep(self, V_anode, DT=None, CFL_factor=None, V_grid=5):
        """Calculate timestep size based on grid data and CFL parameter

         Arguments:
            V_anode (float): Vacuum bias of anode relative to cathode, V.
            DT (float): The dt of each step of the simulation, if not given it
                will be calculated.
            CFL_factor (float): Multiplier to determine the actual timestep
                given the CFL ratio. eg. ``dt = CFL_factor * CFL_ratio``.
            V_grid (float): The vaccum bias of highest-voltage grid relative to
                the cathode V. If not defined a value of 5 is used, which
                is a safe value for finding vmax in the case that the
                electron velocities are only set by their temperature
        """
        if self.grid is None:
            raise ValueError("init_grid must be called before init_timestep!")

        if DT is not None:
            self.simulation.time_step_size = DT
            return DT

        if CFL_factor is None:
            raise ValueError("No CFL_factor passed into init_timestep when DT is not defined!")

        vmax = np.sqrt(
               2*np.e/constants.m_e * max(0.0, abs(V_anode), abs(V_grid))
        )

        if self.geom_str == 'XZ' or self.geom_str == 'RZ':
            dt_local = min(self.dx, self.dz) / vmax * CFL_factor
        if self.geom_str == 'XYZ':
            dt_local = min(self.dx, self.dy, self.dz) / vmax * CFL_factor

        dt = parallel_util.mpiallreduce(data=dt_local, opstring="MIN")
        self.simulation.time_step_size = dt

        return dt


    def init_run(self):
        if self.initialized:
            raise RuntimeError(
                "Attempted to initialize the mwxrun class multiple times.")

        self.simulation.initialize_inputs()
        self.simulation.initialize_warpx()

        self.me = _libwarpx.libwarpx.warpx_getMyProc()
        self.n_procs = _libwarpx.libwarpx.warpx_getNProcs()

        # loop over all species to build their pid_dict's
        for species in self.simulation.species:
            species.init_pid_dict()

        # A level is needed for many things like level number. For now I'm
        # statically setting the default level here. I'm not sure of pitfalls
        # or how to handle it more generally yet.
        self.lev = _libwarpx.libwarpx.warpx_finestLevel()

        self.initialized = True

    def _set_geom_str(self):
        """Set the geom_str variable corresponding to the geometry used.

        Currently supports XZ, RZ, and XYZ.
        """
        # Note these can also be obtained from pywarpx.geometry.coord_sys (0
        # for Cartesian, 1 for RZ) and len(pywarpx.geometry.prob_lo) for num
        # dimensions. Those are set by this solver grid in any case.

        if isinstance(self.grid, picmi.Cartesian2DGrid):
            self.geom_str = 'XZ'
        elif isinstance(self.grid, picmi.Cartesian3DGrid):
            self.geom_str = 'XYZ'
        elif isinstance(self.grid, picmi.CylindricalGrid):
            self.geom_str = 'RZ'
        else:
            raise ValueError("Unrecognized type of pywarpx.picmi Grid.")

    def _set_grid_params(self):
        """Set xmin, xmax, ymin, ymax, zmin, zmax, rmin, rmax.
        Also set nx, ny, nz and dx, dy, dz.

        I'm not currently sure the best way to store and access these, so this
        is subject to change, but gives a basic interface for now.

        Note PICMI uses x/y for 2D, while WarpX uses x/z for 2D (yay!). We
        stick with WarpX, because emission in +z makes more sense for
        compatibility with RZ.
        """
        # Note similar information is stored by pywarpx.geometry.prob_lo,
        # .prob_hi, and .coord_sys, and by pywarpx.amr.n_cell. Those are
        # set by this solver grid in any case, and the present format matches
        # our existing Python framework better.

        if isinstance(self.grid, picmi.Cartesian2DGrid):
            self.xmin = self.grid.xmin
            self.xmax = self.grid.xmax
            self.ymin = 0.0
            self.ymax = 1.0
            self.zmin = self.grid.ymin
            self.zmax = self.grid.ymax
            self.rmin = None
            self.rmax = None

            self.nx = self.grid.nx
            self.ny = 0
            self.nz = self.grid.ny
            self.nr = None

            self.dx = (self.xmax - self.xmin) / self.nx
            self.dy = None
            self.dz = (self.zmax - self.zmin) / self.nz
            self.dr = None

        elif isinstance(self.grid, picmi.Cartesian3DGrid):
            self.xmin = self.grid.xmin
            self.xmax = self.grid.xmax
            self.ymin = self.grid.ymin
            self.ymax = self.grid.ymax
            self.zmin = self.grid.zmin
            self.zmax = self.grid.zmax
            self.rmin = None
            self.rmax = None

            self.nx = self.grid.nx
            self.ny = self.grid.ny
            self.nz = self.grid.ny
            self.nr = None

            self.dx = (self.xmax - self.xmin) / self.nx
            self.dy = (self.ymax - self.ymin) / self.ny
            self.dz = (self.zmax - self.zmin) / self.nz
            self.dr = None

        elif isinstance(self.grid, picmi.CylindricalGrid):
            self.xmin = self.ymin = self.rmin = self.grid.rmin
            self.xmax = self.ymax = self.rmax = self.grid.rmax
            self.zmin = self.grid.zmin
            self.zmax = self.grid.zmax

            self.nx = self.ny = self.nr = self.grid.nr
            self.nz = self.grid.nz

            self.dx = self.dy = self.dr = (self.rmax - self.rmin) / self.nr
            self.dz = (self.zmax - self.zmin) / self.nz

        else:
            raise ValueError("Unrecognized type of pywarpx.picmi Grid.")

    def step(self, sim_control, interval=None):
        sim_control.max_steps = self.simulation.max_steps
        self.step_interval = compute_step(self.simulation, interval)
        while sim_control.check_criteria():
            self.simulation.step(self.step_interval)


    def get_domain_area(self):
        """Return float of simulation domain area in X & Y directions or R depending
        on geometry. Used to get the surface area over which current is emitted or
        absorbed."""
        pos_lims = [mwxrun.xmin, mwxrun.xmax,
                    mwxrun.zmin, mwxrun.zmax]

        if mwxrun.geom_str == "XZ":
            return (pos_lims[1] - pos_lims[0])
        elif mwxrun.geom_str == "XYZ":
            return ((pos_lims[1] - pos_lims[0]) * (pos_lims[3] - pos_lims[2]))
        else:
            # TODO: implement RZ when needed
            raise NotImplementedError("get_domain_area not implemented for RZ geometry")

    def get_it(self):
        """Return the current integer iteration number."""
        return _libwarpx.libwarpx.warpx_getistep(self.lev)

    def get_dt(self):
        """Return the timestep."""
        return _libwarpx.libwarpx.warpx_getdt(self.lev)

    def get_t(self):
        """Return the simulation time."""
        return (self.get_it() - 1.0) * self.get_dt()

    def get_npart(self):
        """Get total number of particles in simulation, across all processors.
        """
        npart = 0
        for spec in self.simulation.species:
            npart += _libwarpx.get_particle_count(spec.name)

        return npart

    def get_npart_species_dict(self):
        """Get total number of particles in simulation per species, across all
        processors.
        """
        npart_dict = {}
        for spec in self.simulation.species:
            if spec.name is None:
                raise ValueError("Unnamed species are not supported.")
            npart_dict[spec.name] = _libwarpx.get_particle_count(spec.name)

        return npart_dict

    def get_gathered_rho_grid(self, species_name=None, include_ghosts=True):
        """Get the full rho on the grid on the root processor.

        Arguments:
            species_name (str or None): If specified the charge density for the
                specific species will be returned (deposited on the grid). If
                None, the current state of rho_fp will be returned.
            include_ghosts (bool): Whether or not to include ghost cells.

        Returns:
            A list with only 1 element - a numpy array with rho on the full
            domain. In place of the numpy array, a reference to an unpopulated
            multifab object is returned on processors other than root.

        """

        if species_name is not None:
            _libwarpx.libwarpx.warpx_depositRhoSpecies(
                ctypes.c_char_p(species_name.encode('utf-8'))
            )
        return _libwarpx._get_mesh_field_list(
            _libwarpx.libwarpx.warpx_getGatheredChargeDensityFP,
            self.lev, None, include_ghosts
        )

    def get_gathered_phi_grid(self, include_ghosts=True):
        """Get the full phi on the grid on the root processor.

        Returns:
            A list with only 1 element - a numpy array with phi on the full
            domain. In place of the numpy array, a reference to an unpopulated
            multifab object is returned on processors other than root.
        """
        return _libwarpx._get_mesh_field_list(
            _libwarpx.libwarpx.warpx_getGatheredPhiFP,
            self.lev, None, include_ghosts
        )

    def set_phi_grid(self, phi_data):
        """Sets phi segments on the grid to input phi data"""
        # only proc 0 has the gathered phi grid so only it should set
        # the phi grid
        if self.me == 0:
            # get phi multifab from warpx
            phi_ptr = _libwarpx._get_mesh_field_list(
                _libwarpx.libwarpx.warpx_getPointerFullPhiFP,
                self.lev, None, True
            )
            try:
                phi_ptr[0][:] = phi_data
            except ValueError as e:
                if 'could not broadcast input array from shape' in str(e):
                    logger.error("Phi data must be the same shape as the phi multifab")
                raise
        _libwarpx.libwarpx.warpx_setPhiGridFP(self.lev)

    def eval_expression_t(self, expr):
        """Function to evaluate an expression that depends on time, at the
        current simulation time using the WarpX parser.

        Arguments:
            expr (str or float): Expression to evaluate.

        Returns:
            (float) Value of the expression at the current simulation time.
        """
        if isinstance(expr, str):
            return _libwarpx.libwarpx.eval_expression_t(
                ctypes.c_char_p(expr.encode('utf-8')), self.lev
            )
        else:
            return expr

mwxrun = MEWarpXRun()
