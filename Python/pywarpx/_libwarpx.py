# Copyright 2017-2019 Andrew Myers, David Grote, Remi Lehe
# Weiqun Zhang
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import atexit
import ctypes
from ctypes.util import find_library as _find_library
# --- This defines the wrapper functions that directly call the underlying compiled routines
import os
import platform
import sys

import numpy as np
from numpy.ctypeslib import ndpointer as _ndpointer

from .Geometry import geometry

try:
    # --- If mpi4py is going to be used, this needs to be imported
    # --- before libwarpx is loaded, because mpi4py calls MPI_Init
    from mpi4py import MPI

    # --- Change MPI Comm type depending on MPICH (int) or OpenMPI (void*)
    if MPI._sizeof(MPI.Comm) == ctypes.sizeof(ctypes.c_int):
        _MPI_Comm_type = ctypes.c_int
    else:
        _MPI_Comm_type = ctypes.c_void_p
except ImportError:
    MPI = None
    _MPI_Comm_type = ctypes.c_void_p

if platform.system() == 'Windows':
    _path_libc = _find_library('msvcrt')
else:
    _path_libc = _find_library('c')
_libc = ctypes.CDLL(_path_libc)
_LP_c_int = ctypes.POINTER(ctypes.c_int)
_LP_c_char = ctypes.c_char_p


class LibWarpX():

    """This class manages the warpx shared object, the library from the compiled C++ code.
    It will only load the library when it is referenced, and this can only be done after
    the geometry is defined so that the version of the library that is needed can be determined.
    Once loaded, all of the setting of function call interfaces is setup.
    """

    def __init__(self):
        # Track whether amrex and warpx have been initialized
        self.initialized = False
        atexit.register(self.finalize)

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
                "Invalid attempt to load libwarpx_so library multiple times."
            )

        # --- Use geometry to determine whether to import the 1D, 2D, 3D or RZ version.
        # --- The geometry must be setup before the lib warpx shared object can be loaded.
        try:
            _prob_lo = geometry.prob_lo
            _dims = geometry.dims
        except AttributeError:
            raise Exception('The shared object could not be loaded. The geometry must be setup before the WarpX shared object can be accessesd. The geometry determines which version of the shared object to load.')

        if _dims == 'RZ':
            self.geometry_dim = 'rz'
        elif (_dims == '1' or _dims == '2' or _dims == '3'):
            self.geometry_dim = '%dd'%len(_prob_lo)
        else:
            raise Exception('Undefined geometry %d'%_dims)

        # this is a plain C/C++ shared library, not a Python module
        if os.name == 'nt':
            mod_ext = 'dll'
        else:
            mod_ext = 'so'
        libname = f'libwarpx.{self.geometry_dim}.{mod_ext}'

        try:
            self.libwarpx_so = ctypes.CDLL(os.path.join(self._get_package_root(), libname))
        except OSError as e:
            value = e.args[0]
            if f'{libname}: cannot open shared object file: No such file or directory' in value:
                raise Exception(f'"{libname}" was not installed. Installation instructions can be found here https://warpx.readthedocs.io/en/latest/install/users.html') from e
            else:
                print("Failed to load the libwarpx shared object library")
                raise

        # WarpX can be compiled using either double or float
        self.libwarpx_so.warpx_Real_size.restype = ctypes.c_int
        self.libwarpx_so.warpx_ParticleReal_size.restype = ctypes.c_int

        _Real_size = self.libwarpx_so.warpx_Real_size()
        _ParticleReal_size = self.libwarpx_so.warpx_ParticleReal_size()

        if _Real_size == 8:
            c_real = ctypes.c_double
            self._numpy_real_dtype = 'f8'
        else:
            c_real = ctypes.c_float
            self._numpy_real_dtype = 'f4'

        if _ParticleReal_size == 8:
            c_particlereal = ctypes.c_double
            self._numpy_particlereal_dtype = 'f8'
        else:
            c_particlereal = ctypes.c_float
            self._numpy_particlereal_dtype = 'f4'

        self.dim = self.libwarpx_so.warpx_SpaceDim()

        # our particle data type, depends on _ParticleReal_size
        _p_struct = (
            [(d, self._numpy_particlereal_dtype) for d in 'xyz'[:self.dim]]
            + [('idcpu', 'u8')]
        )
        self._p_dtype = np.dtype(_p_struct, align=True)

        _numpy_to_ctypes = {}
        _numpy_to_ctypes[self._numpy_particlereal_dtype] = c_particlereal
        _numpy_to_ctypes['u8'] = ctypes.c_uint64

        class Particle(ctypes.Structure):
            _fields_ = [(v[0], _numpy_to_ctypes[v[1]]) for v in _p_struct]

        # some useful typenames
        _LP_particle_p = ctypes.POINTER(ctypes.POINTER(Particle))
        _LP_LP_c_int = ctypes.POINTER(_LP_c_int)
        #_LP_c_void_p = ctypes.POINTER(ctypes.c_void_p)
        _LP_c_real = ctypes.POINTER(c_real)
        _LP_LP_c_real = ctypes.POINTER(_LP_c_real)
        _LP_c_particlereal = ctypes.POINTER(c_particlereal)
        _LP_LP_c_particlereal = ctypes.POINTER(_LP_c_particlereal)
        _LP_LP_c_char = ctypes.POINTER(_LP_c_char)

        # set the arg and return types of the wrapped functions
        self.libwarpx_so.amrex_init.argtypes = (ctypes.c_int, _LP_LP_c_char)
        self.libwarpx_so.amrex_init_with_inited_mpi.argtypes = (ctypes.c_int, _LP_LP_c_char, _MPI_Comm_type)
        self.libwarpx_so.warpx_getParticleStructs.restype = _LP_particle_p
        self.libwarpx_so.warpx_getParticleArrays.restype = _LP_LP_c_particlereal
        self.libwarpx_so.warpx_getParticleCompIndex.restype = ctypes.c_int
        self.libwarpx_so.warpx_getEfield.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getEfieldLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getEfieldCP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getEfieldCPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getEfieldFP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getEfieldFPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getEfieldCP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getEfieldCPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getEfieldFP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getEfieldFPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getBfield.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getBfieldLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getBfieldCP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getBfieldCPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getBfieldFP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getBfieldFPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getBfieldCP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getBfieldCPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getBfieldFP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getBfieldFPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getCurrentDensity.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getCurrentDensityLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getCurrentDensityCP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getCurrentDensityCPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getCurrentDensityFP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getCurrentDensityFPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getCurrentDensityCP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getCurrentDensityCPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getCurrentDensityFP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getCurrentDensityFPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getChargeDensityCP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getChargeDensityCPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getChargeDensityFP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getChargeDensityFPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getPhiFP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getPhiFPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getVectorPotentialFP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getVectorPotentialFPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getFfieldCP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getFfieldCPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getFfieldFP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getFfieldFPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getFfieldCP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getFfieldCPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getFfieldFP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getFfieldFPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getGfieldCP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getGfieldCPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getGfieldFP.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getGfieldFPLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getGfieldCP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getGfieldCPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getGfieldFP_PML.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getGfieldFPLoVects_PML.restype = _LP_c_int
        self.libwarpx_so.warpx_getEdgeLengths.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getEdgeLengthsLoVects.restype = _LP_c_int
        self.libwarpx_so.warpx_getFaceAreas.restype = _LP_LP_c_real
        self.libwarpx_so.warpx_getFaceAreasLoVects.restype = _LP_c_int

        self.libwarpx_so.warpx_sumParticleCharge.restype = c_real
        self.libwarpx_so.warpx_getParticleBoundaryBufferSize.restype = ctypes.c_int
        self.libwarpx_so.warpx_getParticleBoundaryBufferStructs.restype = _LP_LP_c_particlereal
        self.libwarpx_so.warpx_getParticleBoundaryBuffer.restype = _LP_LP_c_particlereal
        self.libwarpx_so.warpx_getParticleBoundaryBufferScrapedSteps.restype = _LP_LP_c_int

        self.libwarpx_so.warpx_getEx_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getEy_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getEz_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getBx_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getBy_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getBz_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getJx_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getJy_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getJz_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getAx_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getAy_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getAz_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getRho_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getPhi_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getF_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getG_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getF_pml_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_getG_pml_nodal_flag.restype = _LP_c_int

        self.libwarpx_so.warpx_get_edge_lengths_x_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_get_edge_lengths_y_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_get_edge_lengths_z_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_get_face_areas_x_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_get_face_areas_y_nodal_flag.restype = _LP_c_int
        self.libwarpx_so.warpx_get_face_areas_z_nodal_flag.restype = _LP_c_int

        #self.libwarpx_so.warpx_getPMLSigma.restype = _LP_c_real
        #self.libwarpx_so.warpx_getPMLSigmaStar.restype = _LP_c_real
        #self.libwarpx_so.warpx_ComputePMLFactors.argtypes = (ctypes.c_int, c_real)
        self.libwarpx_so.warpx_addNParticles.argtypes = (ctypes.c_char_p, ctypes.c_int,
                                                      _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                                      _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                                      _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                                      _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                                      _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                                      _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                                      ctypes.c_int,
                                                      _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                                      ctypes.c_int,
                                                      _ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                                      ctypes.c_int)
        self.libwarpx_so.warpx_convert_id_to_long.argtypes = (_ndpointer(ctypes.c_int64, flags="C_CONTIGUOUS"),
                                                              _ndpointer(Particle, flags="C_CONTIGUOUS"),
                                                              ctypes.c_int)
        self.libwarpx_so.warpx_convert_cpu_to_int.argtypes = (_ndpointer(ctypes.c_int32, flags="C_CONTIGUOUS"),
                                                              _ndpointer(Particle, flags="C_CONTIGUOUS"),
                                                              ctypes.c_int)
        self.libwarpx_so.warpx_getProbLo.restype = c_real
        self.libwarpx_so.warpx_getProbHi.restype = c_real
        self.libwarpx_so.warpx_getCellSize.restype = c_real
        self.libwarpx_so.warpx_getistep.restype = ctypes.c_int
        self.libwarpx_so.warpx_gett_new.restype = c_real
        self.libwarpx_so.warpx_getdt.restype = c_real
        self.libwarpx_so.warpx_maxStep.restype = ctypes.c_int
        self.libwarpx_so.warpx_stopTime.restype = c_real
        self.libwarpx_so.warpx_finestLevel.restype = ctypes.c_int
        self.libwarpx_so.warpx_getMyProc.restype = ctypes.c_int
        self.libwarpx_so.warpx_getNProcs.restype = ctypes.c_int

        self.libwarpx_so.warpx_EvolveE.argtypes = [c_real]
        self.libwarpx_so.warpx_EvolveB.argtypes = [c_real]
        self.libwarpx_so.warpx_FillBoundaryE.argtypes = []
        self.libwarpx_so.warpx_FillBoundaryB.argtypes = []
        self.libwarpx_so.warpx_UpdateAuxilaryData.argtypes = []
        self.libwarpx_so.warpx_SyncCurrent.argtypes = []
        self.libwarpx_so.warpx_PushParticlesandDepose.argtypes = [c_real]
        self.libwarpx_so.warpx_getProbLo.argtypes = [ctypes.c_int]
        self.libwarpx_so.warpx_getProbHi.argtypes = [ctypes.c_int]
        self.libwarpx_so.warpx_getCellSize.argtypes = [ctypes.c_int, ctypes.c_int]
        self.libwarpx_so.warpx_getistep.argtypes = [ctypes.c_int]
        self.libwarpx_so.warpx_setistep.argtypes = [ctypes.c_int, ctypes.c_int]
        self.libwarpx_so.warpx_gett_new.argtypes = [ctypes.c_int]
        self.libwarpx_so.warpx_sett_new.argtypes = [ctypes.c_int, c_real]
        self.libwarpx_so.warpx_getdt.argtypes = [ctypes.c_int]
        self.libwarpx_so.warpx_setPotentialEB.argtypes = [ctypes.c_char_p]

    def _get_boundary_number(self, boundary):
        '''

        Utility function to find the boundary number given a boundary name.

        Parameters
        ----------

        boundary       : str
            The boundary from which to get the scraped particle data. In the
            form x/y/z_hi/lo or eb.

        Returns
        -------
        int
            Integer index in the boundary scraper buffer for the given boundary.
        '''
        if self.geometry_dim == '3d':
            dimensions = {'x' : 0, 'y' : 1, 'z' : 2}
        elif self.geometry_dim == '2d' or self.geometry_dim == 'rz':
            dimensions = {'x' : 0, 'z' : 1}
        elif self.geometry_dim == '1d':
            dimensions = {'z' : 0}
        else:
            raise RuntimeError(f"Unknown simulation geometry: {self.geometry_dim}")

        if boundary != 'eb':
            boundary_parts = boundary.split("_")
            dim_num = dimensions[boundary_parts[0]]
            if boundary_parts[1] == 'lo':
                side = 0
            elif boundary_parts[1] == 'hi':
                side = 1
            else:
                raise RuntimeError(f'Unknown boundary specified: {boundary}')
            boundary_num = 2 * dim_num + side
        else:
            if self.geometry_dim == '3d':
                boundary_num = 6
            elif self.geometry_dim == '2d' or self.geometry_dim == 'rz':
                boundary_num = 4
            elif self.geometry_dim == '1d':
                boundary_num = 2

        return boundary_num

    @staticmethod
    def _array1d_from_pointer(pointer, dtype, size):
        '''

        Function for converting a ctypes pointer to a numpy array

        '''
        if not pointer:
            raise Exception(f'_array1d_from_pointer: pointer is a nullptr')
        if sys.version_info.major >= 3:
            # from where do I import these? this might only work for CPython...
            #PyBuf_READ  = 0x100
            PyBUF_WRITE = 0x200
            buffer_from_memory = ctypes.pythonapi.PyMemoryView_FromMemory
            buffer_from_memory.argtypes = (ctypes.c_void_p, ctypes.c_int, ctypes.c_int)
            buffer_from_memory.restype = ctypes.py_object
            buf = buffer_from_memory(pointer, dtype.itemsize*size, PyBUF_WRITE)
        else:
            buffer_from_memory = ctypes.pythonapi.PyBuffer_FromReadWriteMemory
            buffer_from_memory.restype = ctypes.py_object
            buf = buffer_from_memory(pointer, dtype.itemsize*size)
        return np.frombuffer(buf, dtype=dtype, count=size)

    def getNProcs(self):
        '''

        Get the number of processors

        '''
        return self.libwarpx_so.warpx_getNProcs()

    def getMyProc(self):
        '''

        Get the number of the processor

        '''
        return self.libwarpx_so.warpx_getMyProc()

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

        return self.libwarpx_so.warpx_nCompsSpecies(
            ctypes.c_char_p(species_name.encode('utf-8')))

    def amrex_init(self, argv, mpi_comm=None):
        # --- Construct the ctype list of strings to pass in
        argc = len(argv)

        # note: +1 since there is an extra char-string array element,
        #       that ANSII C requires to be a simple NULL entry
        #       https://stackoverflow.com/a/39096006/2719194
        argvC = (_LP_c_char * (argc+1))()
        for i, arg in enumerate(argv):
            enc_arg = arg.encode('utf-8')
            argvC[i] = _LP_c_char(enc_arg)
        argvC[argc] = _LP_c_char(b"\0")  # +1 element must be NULL

        if mpi_comm is None or MPI is None:
            self.libwarpx_so.amrex_init(argc, argvC)
        else:
            comm_ptr = MPI._addressof(mpi_comm)
            comm_val = _MPI_Comm_type.from_address(comm_ptr)
            self.libwarpx_so.amrex_init_with_inited_mpi(argc, argvC, comm_val)

    def initialize(self, argv=None, mpi_comm=None):
        '''

        Initialize WarpX and AMReX. Must be called before doing anything else.

        '''
        if argv is None:
            argv = sys.argv
        self.amrex_init(argv, mpi_comm)
        self.libwarpx_so.warpx_ConvertLabParamsToBoost()
        self.libwarpx_so.warpx_ReadBCParams()
        if self.geometry_dim == 'rz':
            self.libwarpx_so.warpx_CheckGriddingForRZSpectral()
        self.libwarpx_so.warpx_init()

        self.initialized = True

    def finalize(self, finalize_mpi=1):
        '''

        Call finalize for WarpX and AMReX. Registered to run at program exit.

        '''
        if self.initialized:
            self.libwarpx_so.warpx_finalize()
            self.libwarpx_so.amrex_finalize(finalize_mpi)

    def getistep(self, level=0):
        '''
        Get the current time step number for the specified level

        Parameter
        ---------

        level : int
            The refinement level to reference
        '''

        return self.libwarpx_so.warpx_getistep(level)

    def gett_new(self, level=0):
        '''

        Get the next time for the specified level.

        Parameters
        ----------

        level        : int
            The refinement level to reference
        '''

        return self.libwarpx_so.warpx_gett_new(level)

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

        self.libwarpx_so.warpx_evolve(num_steps);

    def getProbLo(self, direction):
        '''
        Get the values of the lower domain boundary.

        Parameters
        ----------

        direction    : int
            Direction of interest
        '''

        assert 0 <= direction < self.dim, 'Inappropriate direction specified'
        return self.libwarpx_so.warpx_getProbLo(direction)

    def getProbHi(self, direction):
        '''
        Get the values of the upper domain boundary.

        Parameters
        ----------

        direction    : int
            Direction of interest
        '''

        assert 0 <= direction < self.dim, 'Inappropriate direction specified'
        return self.libwarpx_so.warpx_getProbHi(direction)

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

    def add_particles(self, species_name, x=None, y=None, z=None, ux=None, uy=None,
                      uz=None, w=None, unique_particles=True, **kwargs):
        '''
        A function for adding particles to the WarpX simulation.

        Parameters
        ----------

        species_name     : str
            The type of species for which particles will be added

        x, y, z          : arrays or scalars
            The particle positions (default = 0.)

        ux, uy, uz       : arrays or scalars
            The particle momenta (default = 0.)

        w                : array or scalars
            Particle weights (default = 0.)

        unique_particles : bool
            Whether the particles are unique or duplicated on several processes
            (default = True)

        kwargs           : dict
            Containing an entry for all the extra particle attribute arrays. If
            an attribute is not given it will be set to 0.
        '''

        # --- Get length of arrays, set to one for scalars
        lenx = np.size(x)
        leny = np.size(y)
        lenz = np.size(z)
        lenux = np.size(ux)
        lenuy = np.size(uy)
        lenuz = np.size(uz)
        lenw = np.size(w)

        # --- Find the max length of the parameters supplied
        maxlen = 0
        if x is not None:
            maxlen = max(maxlen, lenx)
        if y is not None:
            maxlen = max(maxlen, leny)
        if z is not None:
            maxlen = max(maxlen, lenz)
        if ux is not None:
            maxlen = max(maxlen, lenux)
        if uy is not None:
            maxlen = max(maxlen, lenuy)
        if uz is not None:
            maxlen = max(maxlen, lenuz)
        if w is not None:
            maxlen = max(maxlen, lenw)

        # --- Make sure that the lengths of the input parameters are consistent
        assert x is None or lenx==maxlen or lenx==1, "Length of x doesn't match len of others"
        assert y is None or leny==maxlen or leny==1, "Length of y doesn't match len of others"
        assert z is None or lenz==maxlen or lenz==1, "Length of z doesn't match len of others"
        assert ux is None or lenux==maxlen or lenux==1, "Length of ux doesn't match len of others"
        assert uy is None or lenuy==maxlen or lenuy==1, "Length of uy doesn't match len of others"
        assert uz is None or lenuz==maxlen or lenuz==1, "Length of uz doesn't match len of others"
        assert w is None or lenw==maxlen or lenw==1, "Length of w doesn't match len of others"
        for key, val in kwargs.items():
            assert np.size(val)==1 or len(val)==maxlen, f"Length of {key} doesn't match len of others"

        # --- Broadcast scalars into appropriate length arrays
        # --- If the parameter was not supplied, use the default value
        if lenx == 1:
            x = np.full(maxlen, (x or 0.), self._numpy_particlereal_dtype)
        if leny == 1:
            y = np.full(maxlen, (y or 0.), self._numpy_particlereal_dtype)
        if lenz == 1:
            z = np.full(maxlen, (z or 0.), self._numpy_particlereal_dtype)
        if lenux == 1:
            ux = np.full(maxlen, (ux or 0.), self._numpy_particlereal_dtype)
        if lenuy == 1:
            uy = np.full(maxlen, (uy or 0.), self._numpy_particlereal_dtype)
        if lenuz == 1:
            uz = np.full(maxlen, (uz or 0.), self._numpy_particlereal_dtype)
        if lenw == 1:
            w = np.full(maxlen, (w or 0.), self._numpy_particlereal_dtype)
        for key, val in kwargs.items():
            if np.size(val) == 1:
                kwargs[key] = np.full(
                    maxlen, val, self._numpy_particlereal_dtype
                )

        # --- The -3 is because the comps include the velocites
        nattr = self.get_nattr_species(species_name) - 3
        attr = np.zeros((maxlen, nattr), self._numpy_particlereal_dtype)
        attr[:,0] = w

        for key, vals in kwargs.items():
            # --- The -3 is because components 1 to 3 are velocities
            attr[:,self.get_particle_comp_index(species_name, key)-3] = vals

        nattr_int = 0
        attr_int = np.empty([0], ctypes.c_int)

        # Iff x/y/z/ux/uy/uz are not numpy arrays of the correct dtype, new
        # array copies are made with the correct dtype
        x = x.astype(self._numpy_particlereal_dtype, copy=False)
        y = y.astype(self._numpy_particlereal_dtype, copy=False)
        z = z.astype(self._numpy_particlereal_dtype, copy=False)
        ux = ux.astype(self._numpy_particlereal_dtype, copy=False)
        uy = uy.astype(self._numpy_particlereal_dtype, copy=False)
        uz = uz.astype(self._numpy_particlereal_dtype, copy=False)

        self.libwarpx_so.warpx_addNParticles(
            ctypes.c_char_p(species_name.encode('utf-8')), x.size,
            x, y, z, ux, uy, uz, nattr, attr, nattr_int, attr_int, unique_particles
        )

    def get_particle_count(self, species_name, local=False):
        '''
        Get the number of particles of the specified species in the simulation.

        Parameters
        ----------

        species_name : str
            The species name that the number will be returned for

        local        : bool
            If True the particle count on this processor will be returned.
            Default False.

        Returns
        -------

        int
            An integer count of the number of particles
        '''

        return self.libwarpx_so.warpx_getNumParticles(
            ctypes.c_char_p(species_name.encode('utf-8')), local
        )

    def get_particle_structs(self, species_name, level):
        '''
        This returns a list of numpy arrays containing the particle struct data
        on each tile for this process. The particle data is represented as a structured
        numpy array and contains the particle 'x', 'y', 'z', and 'idcpu'.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

        species_name : str
            The species name that the data will be returned for

        level        : int
            The refinement level to reference

        Returns
        -------

        List of numpy arrays
            The requested particle struct data
        '''

        particles_per_tile = _LP_c_int()
        num_tiles = ctypes.c_int(0)
        data = self.libwarpx_so.warpx_getParticleStructs(
            ctypes.c_char_p(species_name.encode('utf-8')), level,
            ctypes.byref(num_tiles), ctypes.byref(particles_per_tile)
        )

        particle_data = []
        for i in range(num_tiles.value):
            if particles_per_tile[i] == 0:
                continue
            arr = self._array1d_from_pointer(data[i], self._p_dtype, particles_per_tile[i])
            particle_data.append(arr)

        _libc.free(particles_per_tile)
        _libc.free(data)
        return particle_data

    def get_particle_arrays(self, species_name, comp_name, level):
        '''
        This returns a list of numpy arrays containing the particle array data
        on each tile for this process.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

        species_name   : str
            The species name that the data will be returned for

        comp_name      : str
            The component of the array data that will be returned

        level        : int
            The refinement level to reference

        Returns
        -------

        List of numpy arrays
            The requested particle array data
        '''

        particles_per_tile = _LP_c_int()
        num_tiles = ctypes.c_int(0)
        data = self.libwarpx_so.warpx_getParticleArrays(
            ctypes.c_char_p(species_name.encode('utf-8')),
            ctypes.c_char_p(comp_name.encode('utf-8')),
            level, ctypes.byref(num_tiles), ctypes.byref(particles_per_tile)
        )

        particle_data = []
        for i in range(num_tiles.value):
            if particles_per_tile[i] == 0:
                continue
            if not data[i]:
                raise Exception(f'get_particle_arrays: data[i] for i={i} was not initialized')
            arr = np.ctypeslib.as_array(data[i], (particles_per_tile[i],))
            try:
                # This fails on some versions of numpy
                arr.setflags(write=1)
            except ValueError:
                pass
            particle_data.append(arr)

        _libc.free(particles_per_tile)
        _libc.free(data)
        return particle_data

    def get_particle_x(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle 'x'
        positions on each tile.

        '''
        structs = self.get_particle_structs(species_name, level)
        if self.geometry_dim == '3d' or self.geometry_dim == '2d':
            return [struct['x'] for struct in structs]
        elif self.geometry_dim == 'rz':
            return [struct['x']*np.cos(theta) for struct, theta in zip(structs, self.get_particle_theta(species_name))]
        elif self.geometry_dim == '1d':
            raise Exception('get_particle_x: There is no x coordinate with 1D Cartesian')

    def get_particle_y(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle 'y'
        positions on each tile.

        '''
        structs = self.get_particle_structs(species_name, level)
        if self.geometry_dim == '3d':
            return [struct['y'] for struct in structs]
        elif self.geometry_dim == 'rz':
            return [struct['x']*np.sin(theta) for struct, theta in zip(structs, self.get_particle_theta(species_name))]
        elif self.geometry_dim == '1d' or self.geometry_dim == '2d':
            raise Exception('get_particle_y: There is no y coordinate with 1D or 2D Cartesian')

    def get_particle_r(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle 'r'
        positions on each tile.

        '''
        structs = self.get_particle_structs(species_name, level)
        if self.geometry_dim == 'rz':
            return [struct['x'] for struct in structs]
        elif self.geometry_dim == '3d':
            return [np.sqrt(struct['x']**2 + struct['y']**2) for struct in structs]
        elif self.geometry_dim == '2d' or self.geometry_dim == '1d':
            raise Exception('get_particle_r: There is no r coordinate with 1D or 2D Cartesian')

    def get_particle_z(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle 'z'
        positions on each tile.

        '''
        structs = self.get_particle_structs(species_name, level)
        if self.geometry_dim == '3d':
            return [struct['z'] for struct in structs]
        elif self.geometry_dim == 'rz' or self.geometry_dim == '2d':
            return [struct['y'] for struct in structs]
        elif self.geometry_dim == '1d':
            return [struct['x'] for struct in structs]

    def get_particle_id(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle 'id'
        numbers on each tile.

        '''
        ids = []
        structs = self.get_particle_structs(species_name, level)
        for ptile_of_structs in structs:
            arr = np.empty(ptile_of_structs.shape, np.int64)
            self.libwarpx_so.warpx_convert_id_to_long(arr, ptile_of_structs, arr.size)
            ids.append(arr)
        return ids

    def get_particle_cpu(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle 'cpu'
        numbers on each tile.

        '''
        cpus = []
        structs = self.get_particle_structs(species_name, level)
        for ptile_of_structs in structs:
            arr = np.empty(ptile_of_structs.shape, np.int32)
            self.libwarpx_so.warpx_convert_cpu_to_int(arr, ptile_of_structs, arr.size)
            cpus.append(arr)
        return cpus

    def get_particle_weight(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle
        weight on each tile.

        '''

        return self.get_particle_arrays(species_name, 'w', level)

    def get_particle_ux(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle
        x momentum on each tile.

        '''

        return self.get_particle_arrays(species_name, 'ux', level)

    def get_particle_uy(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle
        y momentum on each tile.

        '''

        return self.get_particle_arrays(species_name, 'uy', level)

    def get_particle_uz(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle
        z momentum on each tile.

        '''

        return self.get_particle_arrays(species_name, 'uz', level)

    def get_particle_theta(self, species_name, level=0):
        '''

        Return a list of numpy arrays containing the particle
        theta on each tile.

        '''

        if self.geometry_dim == 'rz':
            return self.get_particle_arrays(species_name, 'theta', level)
        elif self.geometry_dim == '3d':
            structs = self.get_particle_structs(species_name, level)
            return [np.arctan2(struct['y'], struct['x']) for struct in structs]
        elif self.geometry_dim == '2d' or self.geometry_dim == '1d':
            raise Exception('get_particle_theta: There is no theta coordinate with 1D or 2D Cartesian')

    def get_particle_comp_index(self, species_name, pid_name):
        '''
        Get the component index for a given particle attribute. This is useful
        to get the corrent ordering of attributes when adding new particles using
        `add_particles()`.

        Parameters
        ----------

        species_name   : str
            The name of the species

        pid_name       : str
            Name of the component for which the index will be returned

        Returns
        -------

        int
            Integer corresponding to the index of the requested attribute
        '''

        return self.libwarpx_so.warpx_getParticleCompIndex(
            ctypes.c_char_p(species_name.encode('utf-8')),
            ctypes.c_char_p(pid_name.encode('utf-8'))
        )

    def add_real_comp(self, species_name, pid_name, comm=True):
        '''
        Add a real component to the particle data array.

        Parameters
        ----------

        species_name   : str
            The species name for which the new component will be added

        pid_name       : str
            Name that can be used to identify the new component

        comm           : bool
            Should the component be communicated
        '''

        self.libwarpx_so.warpx_addRealComp(
            ctypes.c_char_p(species_name.encode('utf-8')),
            ctypes.c_char_p(pid_name.encode('utf-8')), comm
        )

    def get_species_charge_sum(self, species_name, local=False):
        '''
        Returns the total charge in the simulation due to the given species.

        Parameters
        ----------

        species_name   : str
            The species name for which charge will be summed

        local          : bool
            If True return total charge per processor
        '''

        return self.libwarpx_so.warpx_sumParticleCharge(
            ctypes.c_char_p(species_name.encode('utf-8')), local
        )

    def get_particle_boundary_buffer_size(self, species_name, boundary, local=False):
        '''
        This returns the number of particles that have been scraped so far in the simulation
        from the specified boundary and of the specified species.

        Parameters
        ----------

        species_name   : str
            Return the number of scraped particles of this species

        boundary       : str
            The boundary from which to get the scraped particle data in the
            form x/y/z_hi/lo

        local          : bool
            Whether to only return the number of particles in the current
            processor's buffer
        '''

        return self.libwarpx_so.warpx_getParticleBoundaryBufferSize(
            ctypes.c_char_p(species_name.encode('utf-8')),
            self._get_boundary_number(boundary), local
        )

    def get_particle_boundary_buffer_structs(self, species_name, boundary, level):
        '''
        This returns a list of numpy arrays containing the particle struct data
        for a species that has been scraped by a specific simulation boundary. The
        particle data is represented as a structured numpy array and contains the
        particle 'x', 'y', 'z', and 'idcpu'.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

        species_name : str
            The species name that the data will be returned for

        boundary     : str
            The boundary from which to get the scraped particle data in the
            form x/y/z_hi/lo or eb.

        level        : int
            Which AMR level to retrieve scraped particle data from.
        '''

        particles_per_tile = _LP_c_int()
        num_tiles = ctypes.c_int(0)
        data = self.libwarpx_so.warpx_getParticleBoundaryBufferStructs(
                ctypes.c_char_p(species_name.encode('utf-8')),
                self._get_boundary_number(boundary), level,
                ctypes.byref(num_tiles), ctypes.byref(particles_per_tile)
        )

        particle_data = []
        for i in range(num_tiles.value):
            if particles_per_tile[i] == 0:
                continue
            arr = self._array1d_from_pointer(data[i], self._p_dtype, particles_per_tile[i])
            particle_data.append(arr)

        _libc.free(particles_per_tile)
        _libc.free(data)
        return particle_data

    def get_particle_boundary_buffer(self, species_name, boundary, comp_name, level):
        '''
        This returns a list of numpy arrays containing the particle array data
        for a species that has been scraped by a specific simulation boundary.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            species_name   : str
                The species name that the data will be returned for.

            boundary       : str
                The boundary from which to get the scraped particle data in the
                form x/y/z_hi/lo or eb.

            comp_name      : str
                The component of the array data that will be returned. If
                "step_scraped" the special attribute holding the timestep at
                which a particle was scraped will be returned.

            level          : int
                Which AMR level to retrieve scraped particle data from.
        '''

        particles_per_tile = _LP_c_int()
        num_tiles = ctypes.c_int(0)
        if comp_name == 'step_scraped':
            data = self.libwarpx_so.warpx_getParticleBoundaryBufferScrapedSteps(
                ctypes.c_char_p(species_name.encode('utf-8')),
                self._get_boundary_number(boundary), level,
                ctypes.byref(num_tiles), ctypes.byref(particles_per_tile)
            )
        else:
            data = self.libwarpx_so.warpx_getParticleBoundaryBuffer(
                ctypes.c_char_p(species_name.encode('utf-8')),
                self._get_boundary_number(boundary), level,
                ctypes.byref(num_tiles), ctypes.byref(particles_per_tile),
                ctypes.c_char_p(comp_name.encode('utf-8'))
            )

        particle_data = []
        for i in range(num_tiles.value):
            if particles_per_tile[i] == 0:
                continue
            if not data[i]:
                raise Exception(f'get_particle_arrays: data[i] for i={i} was not initialized')
            arr = np.ctypeslib.as_array(data[i], (particles_per_tile[i],))
            try:
                # This fails on some versions of numpy
                arr.setflags(write=1)
            except ValueError:
                pass
            particle_data.append(arr)

        _libc.free(particles_per_tile)
        _libc.free(data)
        return particle_data

    def clearParticleBoundaryBuffer(self):
        '''

        Clear the buffer that holds the particles lost at the boundaries.

        '''
        self.libwarpx_so.warpx_clearParticleBoundaryBuffer()

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

        if clear_rho:
            from . import fields
            fields.RhoFPWrapper(level, True)[...] = 0.0
        self.libwarpx_so.warpx_depositChargeDensity(
            ctypes.c_char_p(species_name.encode('utf-8')), level
        )
        if sync_rho:
            self.libwarpx_so.warpx_SyncRho()

    def set_potential_EB(self, potential):
        """
        Set the expression string for the embedded boundary potential

        Parameters
        ----------

        potential : str
            The expression string
        """
        self.libwarpx_so.warpx_setPotentialEB(ctypes.c_char_p(potential.encode('utf-8')))

    def _get_mesh_field_list(self, warpx_func, level, direction, include_ghosts):
        """
        Generic routine to fetch the list of field data arrays.
        """
        shapes = _LP_c_int()
        size = ctypes.c_int(0)
        ncomps = ctypes.c_int(0)
        ngrowvect = _LP_c_int()
        if direction is None:
            data = warpx_func(level,
                            ctypes.byref(size), ctypes.byref(ncomps),
                            ctypes.byref(ngrowvect), ctypes.byref(shapes))
        else:
            data = warpx_func(level, direction,
                            ctypes.byref(size), ctypes.byref(ncomps),
                            ctypes.byref(ngrowvect), ctypes.byref(shapes))
        if not data:
            raise Exception('object was not initialized')

        ngvect = [ngrowvect[i] for i in range(self.dim)]
        grid_data = []
        shapesize = self.dim
        if ncomps.value > 1:
            shapesize += 1
        for i in range(size.value):
            shape = tuple([shapes[shapesize*i + d] for d in range(shapesize)])
            # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
            if shape[::-1] == 0:
                continue
            if not data[i]:
                raise Exception(f'get_particle_arrays: data[i] for i={i} was not initialized')
            arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
            try:
                # This fails on some versions of numpy
                arr.setflags(write=1)
            except ValueError:
                pass
            if include_ghosts:
                grid_data.append(arr)
            else:
                grid_data.append(arr[tuple([slice(ngvect[d], -ngvect[d]) for d in range(self.dim)])])

        _libc.free(shapes)
        _libc.free(data)
        return grid_data

    def get_mesh_electric_field(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh electric field
        data on each grid for this process.

        This version is for the full "auxiliary" solution on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getEfield, level, direction, include_ghosts)

    def get_mesh_electric_field_cp(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh electric field
        data on each grid for this process. This version returns the field on
        the coarse patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getEfieldCP, level, direction, include_ghosts)

    def get_mesh_electric_field_fp(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh electric field
        data on each grid for this process. This version returns the field on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getEfieldFP, level, direction, include_ghosts)

    def get_mesh_electric_field_cp_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh electric field
        data on each grid for this process. This version returns the field on
        the coarse patch for the PML for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getEfieldCP_PML, level, direction, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_electric_field_fp_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh electric field
        data on each grid for this process. This version returns the field on
        the fine patch for the PML for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getEfieldFP_PML, level, direction, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_magnetic_field(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh magnetic field
        data on each grid for this process.

        This version is for the full "auxiliary" solution on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getBfield, level, direction, include_ghosts)

    def get_mesh_magnetic_field_cp(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh magnetic field
        data on each grid for this process. This version returns the field on
        the coarse patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getBfieldCP, level, direction, include_ghosts)

    def get_mesh_magnetic_field_fp(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh magnetic field
        data on each grid for this process. This version returns the field on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getBfieldFP, level, direction, include_ghosts)

    def get_mesh_vector_potential_fp(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh vector potential
        data on each grid for this process. This version returns the field on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getVectorPotentialFP, level, direction, include_ghosts)

    def get_mesh_magnetic_field_cp_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh magnetic field
        data on each grid for this process. This version returns the field on
        the coarse patch for the PML for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getBfieldCP_PML, level, direction, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_magnetic_field_fp_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh magnetic field
        data on each grid for this process. This version returns the field on
        the fine patch for the PML for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getBfieldFP_PML, level, direction, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_current_density(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh current density
        data on each grid for this process.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getCurrentDensity, level, direction, include_ghosts)

    def get_mesh_current_density_cp(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh current density
        data on each grid for this process. This version returns the density for
        the coarse patch on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getCurrentDensityCP, level, direction, include_ghosts)

    def get_mesh_current_density_fp(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh current density
        data on each grid for this process. This version returns the density on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getCurrentDensityFP, level, direction, include_ghosts)

    def get_mesh_current_density_cp_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh current density
        data on each grid for this process. This version returns the density for
        the coarse patch for the PML for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getCurrentDensityCP_PML, level, direction, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_current_density_fp_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh current density
        data on each grid for this process. This version returns the density on
        the fine patch for the PML for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getCurrentDensityFP_PML, level, direction, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_charge_density_cp(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh charge density
        data on each grid for this process. This version returns the density for
        the coarse patch on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getChargeDensityCP, level, None, include_ghosts)

    def get_mesh_charge_density_fp(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh charge density
        data on each grid for this process. This version returns the density on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getChargeDensityFP, level, None, include_ghosts)

    def get_mesh_phi_fp(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh electrostatic
        potential data on each grid for this process. This version returns the
        potential on the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''
        return self._get_mesh_field_list(self.libwarpx_so.warpx_getPhiFP, level, None, include_ghosts)

    def get_mesh_F_cp(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh F field
        data on each grid for this process. This version returns the F field for
        the coarse patch on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getFfieldCP, level, None, include_ghosts)

    def get_mesh_F_fp(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh F field
        data on each grid for this process. This version returns the F field on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getFfieldFP, level, None, include_ghosts)

    def get_mesh_F_fp_pml(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh F field
        data on each grid for this process. This version returns the F field on
        the fine patch for the PML on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''
        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getFfieldFP_PML, level, None, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_F_cp_pml(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh F field
        data on each grid for this process. This version returns the F field on
        the coarse patch for the PML on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''
        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getFfieldCP_PML, level, None, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_G_cp(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh G field
        data on each grid for this process. This version returns the G field for
        the coarse patch on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getGfieldCP, level, None, include_ghosts)

    def get_mesh_G_fp(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh G field
        data on each grid for this process. This version returns the G field on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getGfieldFP, level, None, include_ghosts)

    def get_mesh_G_cp_pml(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh G field
        data on each grid for this process. This version returns the G field on
        the coarse patch for the PML on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''
        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getGfieldCP_PML, level, None, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_G_fp_pml(self, level, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh G field
        data on each grid for this process. This version returns the G field on
        the fine patch for the PML on the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''
        try:
            return self._get_mesh_field_list(self.libwarpx_so.warpx_getGfieldFP_PML, level, None, include_ghosts)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_edge_lengths(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh edge lengths
        data on each grid for this process. This version returns the density on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getEdgeLengths, level, direction, include_ghosts)

    def get_mesh_face_areas(self, level, direction, include_ghosts=True):
        '''

        This returns a list of numpy arrays containing the mesh face areas
        data on each grid for this process. This version returns the density on
        the fine patch for the given level.

        The data for the numpy arrays are not copied, but share the underlying
        memory buffer with WarpX. The numpy arrays are fully writeable.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A List of numpy arrays.

        '''

        return self._get_mesh_field_list(self.libwarpx_so.warpx_getFaceAreas, level, direction, include_ghosts)

    def _get_mesh_array_lovects(self, level, direction, include_ghosts=True, getlovectsfunc=None):
        assert(0 <= level and level <= self.libwarpx_so.warpx_finestLevel())

        size = ctypes.c_int(0)
        ngrowvect = _LP_c_int()
        if direction is None:
            data = getlovectsfunc(level, ctypes.byref(size), ctypes.byref(ngrowvect))
        else:
            data = getlovectsfunc(level, direction, ctypes.byref(size), ctypes.byref(ngrowvect))

        if not data:
            raise Exception('object was not initialized')

        lovects_ref = np.ctypeslib.as_array(data, (size.value, self.dim))

        # --- Make a copy of the data to avoid memory problems
        # --- Also, take the transpose to give shape (dims, number of grids)
        lovects = lovects_ref.copy().T

        ng = []
        if include_ghosts:
            for d in range(self.dim):
                ng.append(ngrowvect[d])
        else:
            for d in range(self.dim):
                ng.append(0)
                lovects[d,:] += ngrowvect[d]

        del lovects_ref
        _libc.free(data)
        return lovects, ng

    def get_mesh_electric_field_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        This version is for the full "auxiliary" solution on the given level.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getEfieldLoVects)

    def get_mesh_electric_field_cp_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getEfieldCPLoVects)

    def get_mesh_electric_field_fp_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getEfieldFPLoVects)

    def get_mesh_electric_field_cp_lovects_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getEfieldCPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_electric_field_fp_lovects_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getEfieldFPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_magnetic_field_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        This version is for the full "auxiliary" solution on the given level.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getBfieldLoVects)

    def get_mesh_magnetic_field_cp_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getBfieldCPLoVects)

    def get_mesh_magnetic_field_fp_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getBfieldFPLoVects)

    def get_mesh_vector_potential_fp_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh vector potential field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getVectorPotentialFPLoVects)

    def get_mesh_magnetic_field_cp_lovects_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getBfieldCPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_magnetic_field_fp_lovects_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getBfieldFPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_current_density_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getCurrentDensityLoVects)

    def get_mesh_current_density_cp_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getCurrentDensityCPLoVects)

    def get_mesh_current_density_fp_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getCurrentDensityFPLoVects)

    def get_mesh_current_density_cp_lovects_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getCurrentDensityCPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_current_density_fp_lovects_pml(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getCurrentDensityFPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_charge_density_cp_lovects(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh electric field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getChargeDensityCPLoVects)

    def get_mesh_charge_density_fp_lovects(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh
        charge density data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getChargeDensityFPLoVects)

    def get_mesh_phi_fp_lovects(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh
        electrostatic potential data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getPhiFPLoVects)

    def get_mesh_F_cp_lovects(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh F field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getFfieldCPLoVects)

    def get_mesh_F_fp_lovects(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh F field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getFfieldFPLoVects)

    def get_mesh_F_cp_lovects_pml(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh F field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getFfieldCPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_F_fp_lovects_pml(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh F field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getFfieldFPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_G_cp_lovects(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh G field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getGfieldCPLoVects)

    def get_mesh_G_fp_lovects(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh G field
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getGfieldFPLoVects)

    def get_mesh_G_cp_lovects_pml(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh G field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getGfieldCPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_G_fp_lovects_pml(self, level, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh G field
        data on each PML grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        try:
            return self._get_mesh_array_lovects(level, None, include_ghosts, self.libwarpx_so.warpx_getGfieldFPLoVects_PML)
        except ValueError:
            raise Exception('PML not initialized')

    def get_mesh_edge_lengths_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh edge lengths
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getEdgeLengthsLoVects)

    def get_mesh_face_areas_lovects(self, level, direction, include_ghosts=True):
        '''

        This returns a list of the lo vectors of the arrays containing the mesh face areas
        data on each grid for this process.

        Parameters
        ----------

            level          : the AMR level to get the data for
            direction      : the component of the data you want
            include_ghosts : whether to include ghost zones or not

        Returns
        -------

            A 2d numpy array of the lo vector for each grid with the shape (dims, number of grids)

        '''
        return self._get_mesh_array_lovects(level, direction, include_ghosts, self.libwarpx_so.warpx_getFaceAreasLoVects)

    def _get_nodal_flag(self, getdatafunc):
        data = getdatafunc()
        if not data:
            raise Exception('object was not initialized')

        nodal_flag_ref = np.ctypeslib.as_array(data, (self.dim,))

        # --- Make a copy of the data to avoid memory problems
        nodal_flag = nodal_flag_ref.copy()

        del nodal_flag_ref
        _libc.free(data)
        return nodal_flag

    def get_Ex_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Ex along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getEx_nodal_flag)

    def get_Ey_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Ey along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getEy_nodal_flag)

    def get_Ez_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Ez along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getEz_nodal_flag)

    def get_Bx_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Bx along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getBx_nodal_flag)

    def get_By_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for By along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getBy_nodal_flag)

    def get_Bz_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Bz along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getBz_nodal_flag)

    def get_Jx_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Jx along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getJx_nodal_flag)

    def get_Jy_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Jy along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getJy_nodal_flag)

    def get_Jz_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Jz along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getJz_nodal_flag)

    def get_Ax_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Ax along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getAx_nodal_flag)

    def get_Ay_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Ay along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getAy_nodal_flag)

    def get_Az_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Az along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getAz_nodal_flag)

    def get_Rho_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Rho along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getRho_nodal_flag)

    def get_Phi_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for Phi along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getPhi_nodal_flag)

    def get_F_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for F along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getF_nodal_flag)

    def get_G_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for G along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getG_nodal_flag)

    def get_edge_lengths_x_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for the x edge lengths along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_get_edge_lengths_x_nodal_flag)

    def get_edge_lengths_y_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for the y edge lengths along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_get_edge_lengths_y_nodal_flag)

    def get_edge_lengths_z_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for the z edge lengths along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_get_edge_lengths_z_nodal_flag)

    def get_face_areas_x_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for the x face areas along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_get_face_areas_x_nodal_flag)

    def get_face_areas_y_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for the y face areas along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_get_face_areas_y_nodal_flag)

    def get_face_areas_z_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for the z face areas along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_get_face_areas_z_nodal_flag)

    def get_F_pml_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for F in the PML along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getF_pml_nodal_flag)

    def get_G_pml_nodal_flag(self):
        '''
        This returns a 1d array of the nodal flags for G in the PML along each direction. A 1 means node centered, and 0 cell centered.
        '''
        return self._get_nodal_flag(self.libwarpx_so.warpx_getG_pml_nodal_flag)


libwarpx = LibWarpX()
