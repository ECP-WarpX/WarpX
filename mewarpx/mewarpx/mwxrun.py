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
from pywarpx import _libwarpx, picmi


class MEWarpXRun(object):

    """Hold base information and wrappers for the simulation. This should be a
    singleton object.

    This also does the initialization steps, using init_run to enforce an
    initialization order.
    """

    def __init__(self):
        self.initialized = False
        self.simulation = picmi.Simulation(verbose=0)

    def init_run(self):
        if self.initialized:
            raise RuntimeError(
                "Attempted to initialize the mwxrun class multiple times.")
        self.initialized = True

        self.simulation.initialize_inputs()
        self.simulation.initialize_warpx()

        self.me = _libwarpx.libwarpx.warpx_getMyProc()
        self.n_procs = _libwarpx.libwarpx.warpx_getNProcs()

        self._set_geom_str()
        self._set_grid_params()

    def _set_geom_str(self):
        """Set the geom_str variable corresponding to the geometry used.

        Currently supports XZ, RZ, and XYZ.
        """
        # Note these can also be obtained from pywarpx.geometry.coord_sys (0
        # for Cartesian, 1 for RZ) and len(pywarpx.geometry.prob_lo) for num
        # dimensions. Those are set by this solver grid in any case.
        grid = self.simulation.solver.grid

        if isinstance(grid, picmi.Cartesian2DGrid):
            self.geom_str = 'XZ'
        elif isinstance(grid, picmi.Cartesian3DGrid):
            self.geom_str = 'XYZ'
        elif isinstance(grid, picmi.CylindricalGrid):
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
        grid = self.simulation.solver.grid

        if isinstance(grid, picmi.Cartesian2DGrid):
            self.xmin = grid.xmin
            self.xmax = grid.xmax
            self.ymin = 0.0
            self.ymax = 1.0
            self.zmin = grid.ymin
            self.zmax = grid.ymax
            self.rmin = None
            self.rmax = None

            self.nx = grid.nx
            self.ny = 0
            self.nz = grid.ny
            self.nr = None

            self.dx = (self.xmax - self.xmin) / self.nx
            self.dy = None
            self.dz = (self.zmax - self.zmin) / self.nz
            self.dr = None

        elif isinstance(grid, picmi.Cartesian3DGrid):
            self.xmin = grid.xmin
            self.xmax = grid.xmax
            self.ymin = grid.ymin
            self.ymax = grid.ymax
            self.zmin = grid.zmin
            self.zmax = grid.zmax
            self.rmin = None
            self.rmax = None

            self.nx = grid.nx
            self.ny = grid.ny
            self.nz = grid.ny
            self.nr = None

            self.dx = (self.xmax - self.xmin) / self.nx
            self.dy = (self.ymax - self.ymin) / self.ny
            self.dz = (self.zmax - self.zmin) / self.nz
            self.dr = None

        elif isinstance(grid, picmi.CylindricalGrid):
            self.xmin = self.ymin = self.rmin = grid.rmin
            self.xmax = self.ymax = self.rmax = grid.rmax
            self.zmin = grid.zmin
            self.zmax = grid.zmax

            self.nx = self.ny = self.nr = grid.nr
            self.nz = grid.nz

            self.dx = self.dy = self.dr = (self.rmax - self.rmin) / self.nr
            self.dz = (self.zmax - self.zmin) / self.nz

        else:
            raise ValueError("Unrecognized type of pywarpx.picmi Grid.")

        # A level is needed for many things like level number. For now I'm
        # statically setting the default level here. I'm not sure of pitfalls
        # or how to handle it more generally yet.
        self.lev = _libwarpx.libwarpx.warpx_finestLevel()

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
            npart += _libwarpx.libwarpx.warpx_getNumParticles(
                spec.species_number)

        return npart

    def get_npart_species_dict(self):
        """Get total number of particles in simulation per species, across all
        processors.
        """
        npart_dict = {}
        for spec in self.simulation.species:
            if spec.name is None:
                raise ValueError("Unnamed species are not supported.")
            npart_dict[spec.name] = _libwarpx.libwarpx.warpx_getNumParticles(
                spec.species_number)

        return npart_dict

    def get_rho_grid(self):
        """Get rho segments on the grid for each tile of each processor.

        Returns:
            A list of numpy arrays, the list has an array for every tile and
            each array has dimensions given by the number of cells in that tile.
        """
        return _libwarpx.get_mesh_charge_density_fp(self.lev)

    def get_gathered_rho_grid(self):
        """Get the full rho on the grid on the root processor.

        Returns:
            A list with only 1 element - a numpy array with rho on the full
            domain. In place of the numpy array, a reference to an unpopulated
            multifab object is returned on processors other than root.

        """
        return _libwarpx.get_gathered_charge_density_fp(self.lev)

    def get_phi_grid(self):
        """Get phi segments on the grid for each tile of each processor.

        Returns:
            A list of numpy arrays, the list has an array for every tile and
            each array has dimensions given by the number of cells in that tile.
        """
        return _libwarpx.get_mesh_phi_fp(self.lev)

    def get_gathered_phi_grid(self):
        """Get the full phi on the grid on the root processor.

        Returns:
            A list with only 1 element - a numpy array with phi on the full
            domain. In place of the numpy array, a reference to an unpopulated
            multifab object is returned on processors other than root.
        """
        return _libwarpx.get_gathered_phi_fp(self.lev)

    def set_phi_grid(self, phi_data):
        """Sets phi segments on the grid to input phi data"""
        # only proc 0 has the gathered phi grid so only it should set
        # the phi grid
        if self.me == 0:
            # get phi multifab from warpx
            phi_ptr = _libwarpx.get_pointer_full_phi_fp(self.lev)
            try:
                phi_ptr[0][:] = phi_data
            except ValueError as e:
                if 'could not broadcast input array from shape' in str(e):
                    print("Phi data must be the same shape as the phi multifab")
                raise
        _libwarpx.set_phi_grid_fp(self.lev)


mwxrun = MEWarpXRun()
