# Copyright 2017-2022 The WarpX Community
#
# This file is part of WarpX. It defines the wrapper functions that directly
# call the underlying compiled routines through pybind11.
#
# NOTE: We will reduce the libwarpx.py level of abstraction eventually!
# Please add new functionality directly to pybind11-bound modules
# and call them via sim.extension.libwarpx_so. ... and sim.extension.warpx.
# ... from user code.
#
# Authors: Axel Huebl, Andrew Myers, David Grote, Remi Lehe, Weiqun Zhang
#
# License: BSD-3-Clause-LBNL

import atexit
import os
import sys

import numpy as np

from .Geometry import geometry


class LibWarpX():
    """This class manages the WarpX classes, part of the Python module from the compiled C++ code.
    It will only load the library when it is referenced, and this can only be done after
    the geometry is defined so that the version of the library that is needed can be determined.
    Once loaded, all the settings of function call interfaces are set up.
    """

    def __init__(self):
        # Track whether amrex and warpx have been initialized
        self.initialized = False
        atexit.register(self.finalize)

        # set once libwarpx_so is loaded
        self.__version__ = None

    def __getattr__(self, attribute):
        if attribute == 'libwarpx_so':
            # If the 'libwarpx_so' is referenced, load it.
            # Once loaded, it gets added to the dictionary so this code won't be called again.
            self.load_library()
            return self.__dict__[attribute]
        else:
            # For any other attribute, call the built-in routine - this should always
            # return an AttributeException.
            return self.__getattribute__(attribute)

    def _get_package_root(self):
        '''
        Get the path to the installation location (where libwarpx.so would be installed).
        '''
        cur = os.path.abspath(__file__)
        while True:
            name = os.path.basename(cur)
            if name == 'pywarpx':
                return cur
            elif not name:
                return ''
            cur = os.path.dirname(cur)

    def load_library(self):

        if 'libwarpx_so' in self.__dict__:
            raise RuntimeError(
                "Invalid attempt to load the pybind11 bindings library multiple times. "
                "Note that multiple AMReX/WarpX geometries cannot be loaded yet into the same Python process. "
                "Please write separate scripts for each geometry."
            )

        # --- Use geometry to determine whether to import the 1D, 2D, 3D or RZ version.
        # --- The geometry must be setup before the lib warpx shared object can be loaded.
        try:
            _prob_lo = geometry.prob_lo
            _dims = geometry.dims
        except AttributeError:
            raise Exception('The shared object could not be loaded. The geometry must be setup before the WarpX pybind11 module can be accessesd. The geometry determines which version of the shared object to load.')

        if _dims == 'RZ':
            self.geometry_dim = 'rz'
        elif (_dims == '1' or _dims == '2' or _dims == '3'):
            self.geometry_dim = '%dd'%len(_prob_lo)
        else:
            raise Exception('Undefined geometry %d'%_dims)

        try:
            if self.geometry_dim == "1d":
                import amrex.space1d as amr
                self.amr = amr
                from . import warpx_pybind_1d as cxx_1d
                self.libwarpx_so = cxx_1d
                self.dim = 1
            elif self.geometry_dim == "2d":
                import amrex.space2d as amr
                self.amr = amr
                from . import warpx_pybind_2d as cxx_2d
                self.libwarpx_so = cxx_2d
                self.dim = 2
            elif self.geometry_dim == "rz":
                import amrex.space2d as amr
                self.amr = amr
                from . import warpx_pybind_rz as cxx_rz
                self.libwarpx_so = cxx_rz
                self.dim = 2
            elif self.geometry_dim == "3d":
                import amrex.space3d as amr
                self.amr = amr
                from . import warpx_pybind_3d as cxx_3d
                self.libwarpx_so = cxx_3d
                self.dim = 3
        except ImportError:
            raise Exception(f"Dimensionality '{self.geometry_dim}' was not compiled in this Python install. Please recompile with -DWarpX_DIMS={_dims}")

        self.__version__ = self.libwarpx_so.__version__

    def getNProcs(self):
        '''

        Get the number of processors

        '''
        return self.libwarpx_so.getNProcs()

    def getMyProc(self):
        '''

        Get the number of the processor

        '''
        return self.libwarpx_so.getMyProc()

    def get_nattr(self):
        '''

        Get the number of extra particle attributes.

        '''
        # --- The -3 is because the comps include the velocites
        return self.libwarpx_so.warpx_nComps() - 3

    def get_nattr_species(self, species_name):
        '''
        Get the number of real attributes for the given species.

        Parameters
        ----------

        species_name: str
            Name of the species
        '''
        warpx = self.libwarpx_so.get_instance()
        mpc = warpx.multi_particle_container()
        pc = mpc.get_particle_container_from_name(species_name)
        return pc.num_real_comps()

    def amrex_init(self, argv, mpi_comm=None):
        if mpi_comm is None: # or MPI is None:
            self.libwarpx_so.amrex_init(argv)
        else:
            raise Exception('mpi_comm argument not yet supported')

    def initialize(self, argv=None, mpi_comm=None):
        '''

        Initialize WarpX and AMReX. Must be called before doing anything else.

        '''
        if argv is None:
            argv = sys.argv
        self.amrex_init(argv, mpi_comm)
        self.warpx = self.libwarpx_so.get_instance()
        self.warpx.initialize_data()
        self.libwarpx_so.execute_python_callback("afterinit")
        self.libwarpx_so.execute_python_callback("particleloader")

        #self.libwarpx_so.warpx_init()

        self.initialized = True

    def finalize(self, finalize_mpi=1):
        '''

        Call finalize for WarpX and AMReX. Registered to run at program exit.

        '''
        # TODO: simplify, part of pyAMReX already
        if self.initialized:
            del self.warpx
            # The call to warpx_finalize causes a crash - don't know why
            #self.libwarpx_so.warpx_finalize()
            self.libwarpx_so.amrex_finalize()

            from pywarpx import callbacks
            callbacks.clear_all()

    def getistep(self, level=0):
        '''
        Get the current time step number for the specified level

        Parameter
        ---------

        level : int
            The refinement level to reference
        '''

        return self.warpx.getistep(level)

    def gett_new(self, level=0):
        '''

        Get the next time for the specified level.

        Parameters
        ----------

        level        : int
            The refinement level to reference
        '''

        return self.warpx.gett_new(level)

    def evolve(self, num_steps=-1):
        '''
        Evolve the simulation for num_steps steps. If num_steps=-1,
        the simulation will be run until the end as specified in the
        inputs file.

        Parameters
        ----------

        num_steps: int
            The number of steps to take
        '''

        self.warpx.evolve(num_steps)

    def getProbLo(self, direction, level=0):
        '''
        Get the values of the lower domain boundary.

        Parameters
        ----------

        direction    : int
            Direction of interest
        '''

        assert 0 <= direction < self.dim, 'Inappropriate direction specified'
        return self.warpx.Geom(level).ProbLo(direction)

    def getProbHi(self, direction, level=0):
        '''
        Get the values of the upper domain boundary.

        Parameters
        ----------

        direction    : int
            Direction of interest
        '''

        assert 0 <= direction < self.dim, 'Inappropriate direction specified'
        return self.warpx.Geom(level).ProbHi(direction)

    def getCellSize(self, direction, level=0):
        '''
        Get the cell size in the given direction and on the given level.

        Parameters
        ----------

        direction    : int
            Direction of interest

        level        : int
            The refinement level to reference
        '''

        assert 0 <= direction < 3, 'Inappropriate direction specified'
        assert 0 <= level and level <= self.libwarpx_so.warpx_finestLevel(), 'Inappropriate level specified'
        return self.libwarpx_so.warpx_getCellSize(direction, level)

    #def get_sigma(self, direction):
    #    '''
    #
    #    Return the 'sigma' PML coefficients for the electric field
    #    in a given direction.
    #
    #    '''
    #
    #    size = ctypes.c_int(0)
    #    data = self.libwarpx_so.warpx_getPMLSigma(direction, ctypes.byref(size))
    #    arr = np.ctypeslib.as_array(data, (size.value,))
    #    arr.setflags(write=1)
    #    return arr
    #
    #
    #def get_sigma_star(self, direction):
    #    '''
    #
    #    Return the 'sigma*' PML coefficients for the magnetic field
    #    in the given direction.
    #
    #    '''
    #
    #    size = ctypes.c_int(0)
    #    data = self.libwarpx_so.warpx_getPMLSigmaStar(direction, ctypes.byref(size))
    #    arr = np.ctypeslib.as_array(data, (size.value,))
    #    arr.setflags(write=1)
    #    return arr
    #
    #
    #def compute_pml_factors(self, lev, dt):
    #    '''
    #
    #    This recomputes the PML coefficients for a given level, using the
    #    time step dt. This needs to be called after modifying the coefficients
    #    from Python.
    #
    #    '''
    #
    #    self.libwarpx_so.warpx_ComputePMLFactors(lev, dt)



    def depositChargeDensity(self, species_name, level, clear_rho=True, sync_rho=True):
        '''
        Deposit the specified species' charge density in rho_fp in order to
        access that data via pywarpx.fields.RhoFPWrapper().

        Parameters
        ----------

        species_name   : str
            The species name that will be deposited.

        level          : int
            Which AMR level to retrieve scraped particle data from.

        clear_rho      : bool
            If True, zero out rho_fp before deposition.

        sync_rho       : bool
            If True, perform MPI exchange and properly set boundary cells for rho_fp.
        '''
        rho_fp = self.warpx.multifab(f'rho_fp[level={level}]')

        if rho_fp is None:
            raise RuntimeError("Multifab `rho_fp` is not allocated.")
            # ablastr::warn_manager::WMRecordWarning(
            #     "WarpXWrappers", "rho_fp is not allocated",
            #     ablastr::warn_manager::WarnPriority::low
            # );
            # return

        if clear_rho:
            rho_fp.set_val(0.0)

        # deposit the charge density from the desired species
        mypc = self.warpx.multi_particle_container()
        myspc = mypc.get_particle_container_from_name(species_name)
        myspc.deposit_charge(rho_fp, level)

        if self.geometry_dim == 'rz':
            self.warpx.apply_inverse_volume_scaling_to_charge_density(rho_fp, level)

        if sync_rho:
            self.warpx.sync_rho()


    def set_potential_EB(self, potential):
        """
        Set the expression string for the embedded boundary potential

        Parameters
        ----------

        potential : str
            The expression string
        """
        self.warpx.set_potential_on_eb(potential)


libwarpx = LibWarpX()
