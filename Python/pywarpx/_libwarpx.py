# Copyright 2017-2019 Andrew Myers, David Grote, Remi Lehe
# Weiqun Zhang
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# --- This defines the wrapper functions that directly call the underlying compiled routines
import os
import sys
import platform
import atexit
import ctypes
from ctypes.util import find_library as _find_library
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

# --- Is there a better way of handling constants?
clight = 2.99792458e+8 # m/s

def _get_package_root():
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

# --- Use geometry to determine whether to import the 2D or 3D version.
# --- This assumes that the input is setup before this module is imported,
# --- which should normally be the case.
# --- Default to 3D if geometry is not setup yet.
try:
    _prob_lo = geometry.prob_lo
    _coord_sys = geometry.coord_sys
except AttributeError:
    geometry_dim = '3d'
else:
    if _coord_sys == 0:
        geometry_dim = '%dd'%len(_prob_lo)
    elif _coord_sys == 1:
        geometry_dim = 'rz'
    else:
        raise Exception('Undefined coordinate system %d'%_coord_sys)
    del _prob_lo, _coord_sys

if platform.system() == 'Windows':
    path_libc = _find_library('msvcrt')
else:
    path_libc = _find_library('c')
_libc = ctypes.CDLL(path_libc)

# this is a plain C/C++ shared library, not a Python module
if os.name == 'nt':
    mod_ext = "dll"
else:
    mod_ext = "so"
libname = "libwarpx.{0}.{1}".format(geometry_dim, mod_ext)

try:
    libwarpx = ctypes.CDLL(os.path.join(_get_package_root(), libname))
except OSError as e:
    value = e.args[0]
    if f'{libname}: cannot open shared object file: No such file or directory' in value:
        raise Exception(f'"{libname}" was not installed. Installation instructions can be found here https://warpx.readthedocs.io/en/latest/install/users.html') from e
    else:
        print("Failed to load the libwarpx shared object library")
        raise

# WarpX can be compiled using either double or float
libwarpx.warpx_Real_size.restype = ctypes.c_int
libwarpx.warpx_ParticleReal_size.restype = ctypes.c_int

_Real_size = libwarpx.warpx_Real_size()
_ParticleReal_size = libwarpx.warpx_ParticleReal_size()

if _Real_size == 8:
    c_real = ctypes.c_double
    _numpy_real_dtype = 'f8'
else:
    c_real = ctypes.c_float
    _numpy_real_dtype = 'f4'

if _ParticleReal_size == 8:
    c_particlereal = ctypes.c_double
    _numpy_particlereal_dtype = 'f8'
else:
    c_particlereal = ctypes.c_float
    _numpy_particlereal_dtype = 'f4'

dim = libwarpx.warpx_SpaceDim()

# our particle data type, depends on _ParticleReal_size
_p_struct = [(d, _numpy_particlereal_dtype) for d in 'xyz'[:dim]] + [('id', 'i4'), ('cpu', 'i4')]
_p_dtype = np.dtype(_p_struct, align=True)

_numpy_to_ctypes = {}
_numpy_to_ctypes[_numpy_particlereal_dtype] = c_particlereal
_numpy_to_ctypes['i4'] = ctypes.c_int

class Particle(ctypes.Structure):
    _fields_ = [(v[0], _numpy_to_ctypes[v[1]]) for v in _p_struct]


# some useful typenames
_LP_particle_p = ctypes.POINTER(ctypes.POINTER(Particle))
_LP_c_int = ctypes.POINTER(ctypes.c_int)
_LP_LP_c_int = ctypes.POINTER(_LP_c_int)
_LP_c_void_p = ctypes.POINTER(ctypes.c_void_p)
_LP_c_real = ctypes.POINTER(c_real)
_LP_LP_c_real = ctypes.POINTER(_LP_c_real)
_LP_c_particlereal = ctypes.POINTER(c_particlereal)
_LP_LP_c_particlereal = ctypes.POINTER(_LP_c_particlereal)
_LP_c_char = ctypes.POINTER(ctypes.c_char)
_LP_LP_c_char = ctypes.POINTER(_LP_c_char)

# this is a function for converting a ctypes pointer to a numpy array
def _array1d_from_pointer(pointer, dtype, size):
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


# set the arg and return types of the wrapped functions
libwarpx.amrex_init.argtypes = (ctypes.c_int, _LP_LP_c_char)
libwarpx.amrex_init_with_inited_mpi.argtypes = (ctypes.c_int, _LP_LP_c_char, _MPI_Comm_type)
libwarpx.warpx_getParticleStructs.restype = _LP_particle_p
libwarpx.warpx_getParticleArrays.restype = _LP_LP_c_particlereal
libwarpx.warpx_getParticleCompIndex.restype = ctypes.c_int
libwarpx.warpx_getEfield.restype = _LP_LP_c_real
libwarpx.warpx_getEfieldLoVects.restype = _LP_c_int
libwarpx.warpx_getEfieldCP.restype = _LP_LP_c_real
libwarpx.warpx_getEfieldCPLoVects.restype = _LP_c_int
libwarpx.warpx_getEfieldFP.restype = _LP_LP_c_real
libwarpx.warpx_getEfieldFPLoVects.restype = _LP_c_int
libwarpx.warpx_getEfieldCP_PML.restype = _LP_LP_c_real
libwarpx.warpx_getEfieldCPLoVects_PML.restype = _LP_c_int
libwarpx.warpx_getEfieldFP_PML.restype = _LP_LP_c_real
libwarpx.warpx_getEfieldFPLoVects_PML.restype = _LP_c_int
libwarpx.warpx_getBfield.restype = _LP_LP_c_real
libwarpx.warpx_getBfieldLoVects.restype = _LP_c_int
libwarpx.warpx_getBfieldCP.restype = _LP_LP_c_real
libwarpx.warpx_getBfieldCPLoVects.restype = _LP_c_int
libwarpx.warpx_getBfieldFP.restype = _LP_LP_c_real
libwarpx.warpx_getBfieldFPLoVects.restype = _LP_c_int
libwarpx.warpx_getBfieldCP_PML.restype = _LP_LP_c_real
libwarpx.warpx_getBfieldCPLoVects_PML.restype = _LP_c_int
libwarpx.warpx_getBfieldFP_PML.restype = _LP_LP_c_real
libwarpx.warpx_getBfieldFPLoVects_PML.restype = _LP_c_int
libwarpx.warpx_getCurrentDensity.restype = _LP_LP_c_real
libwarpx.warpx_getCurrentDensityLoVects.restype = _LP_c_int
libwarpx.warpx_getCurrentDensityCP.restype = _LP_LP_c_real
libwarpx.warpx_getCurrentDensityCPLoVects.restype = _LP_c_int
libwarpx.warpx_getCurrentDensityFP.restype = _LP_LP_c_real
libwarpx.warpx_getCurrentDensityFPLoVects.restype = _LP_c_int
libwarpx.warpx_getCurrentDensityCP_PML.restype = _LP_LP_c_real
libwarpx.warpx_getCurrentDensityCPLoVects_PML.restype = _LP_c_int
libwarpx.warpx_getCurrentDensityFP_PML.restype = _LP_LP_c_real
libwarpx.warpx_getCurrentDensityFPLoVects_PML.restype = _LP_c_int
libwarpx.warpx_getChargeDensityCP.restype = _LP_LP_c_real
libwarpx.warpx_getChargeDensityCPLoVects.restype = _LP_c_int
libwarpx.warpx_getChargeDensityFP.restype = _LP_LP_c_real
libwarpx.warpx_getChargeDensityFPLoVects.restype = _LP_c_int
libwarpx.warpx_getPhiFP.restype = _LP_LP_c_real
libwarpx.warpx_getPhiFPLoVects.restype = _LP_c_int
libwarpx.warpx_getFfieldCP.restype = _LP_LP_c_real
libwarpx.warpx_getFfieldCPLoVects.restype = _LP_c_int
libwarpx.warpx_getFfieldFP.restype = _LP_LP_c_real
libwarpx.warpx_getFfieldFPLoVects.restype = _LP_c_int
libwarpx.warpx_getGfieldCP.restype = _LP_LP_c_real
libwarpx.warpx_getGfieldCPLoVects.restype = _LP_c_int
libwarpx.warpx_getGfieldFP.restype = _LP_LP_c_real
libwarpx.warpx_getGfieldFPLoVects.restype = _LP_c_int
libwarpx.warpx_getParticleBoundaryBufferSize.restype = ctypes.c_int
libwarpx.warpx_getParticleBoundaryBuffer.restype = _LP_LP_c_particlereal
libwarpx.warpx_getParticleBoundaryBufferScrapedSteps.restype = _LP_LP_c_int

libwarpx.warpx_getEx_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getEy_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getEz_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getBx_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getBy_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getBz_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getJx_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getJy_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getJz_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getRho_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getPhi_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getF_nodal_flag.restype = _LP_c_int
libwarpx.warpx_getG_nodal_flag.restype = _LP_c_int

#libwarpx.warpx_getPMLSigma.restype = _LP_c_real
#libwarpx.warpx_getPMLSigmaStar.restype = _LP_c_real
#libwarpx.warpx_ComputePMLFactors.argtypes = (ctypes.c_int, c_real)
libwarpx.warpx_addNParticles.argtypes = (ctypes.c_char_p, ctypes.c_int,
                                         _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                         _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                         _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                         _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                         _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                         _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                         ctypes.c_int,
                                         _ndpointer(c_particlereal, flags="C_CONTIGUOUS"),
                                         ctypes.c_int)

libwarpx.warpx_getProbLo.restype = c_real
libwarpx.warpx_getProbHi.restype = c_real
libwarpx.warpx_getCellSize.restype = c_real
libwarpx.warpx_getistep.restype = ctypes.c_int
libwarpx.warpx_gett_new.restype = c_real
libwarpx.warpx_getdt.restype = c_real
libwarpx.warpx_maxStep.restype = ctypes.c_int
libwarpx.warpx_stopTime.restype = c_real
libwarpx.warpx_finestLevel.restype = ctypes.c_int
libwarpx.warpx_getMyProc.restype = ctypes.c_int
libwarpx.warpx_getNProcs.restype = ctypes.c_int

libwarpx.warpx_EvolveE.argtypes = [c_real]
libwarpx.warpx_EvolveB.argtypes = [c_real]
libwarpx.warpx_FillBoundaryE.argtypes = []
libwarpx.warpx_FillBoundaryB.argtypes = []
libwarpx.warpx_UpdateAuxilaryData.argtypes = []
libwarpx.warpx_SyncCurrent.argtypes = []
libwarpx.warpx_PushParticlesandDepose.argtypes = [c_real]
libwarpx.warpx_getProbLo.argtypes = [ctypes.c_int]
libwarpx.warpx_getProbHi.argtypes = [ctypes.c_int]
libwarpx.warpx_getCellSize.argtypes = [ctypes.c_int, ctypes.c_int]
libwarpx.warpx_getistep.argtypes = [ctypes.c_int]
libwarpx.warpx_setistep.argtypes = [ctypes.c_int, ctypes.c_int]
libwarpx.warpx_gett_new.argtypes = [ctypes.c_int]
libwarpx.warpx_sett_new.argtypes = [ctypes.c_int, c_real]
libwarpx.warpx_getdt.argtypes = [ctypes.c_int]

def get_nattr():
    '''

    Get the number of extra attributes.

    '''
    # --- The -3 is because the comps include the velocites
    return libwarpx.warpx_nComps() - 3

def get_nattr_species(species_name):
    '''

    Get the number of real attributes for the given species.

    '''
    return libwarpx.warpx_nCompsSpecies(
        ctypes.c_char_p(species_name.encode('utf-8')))

def amrex_init(argv, mpi_comm=None):
    # --- Construct the ctype list of strings to pass in
    argc = len(argv)
    argvC = (_LP_c_char * (argc+1))()
    for i, arg in enumerate(argv):
        enc_arg = arg.encode('utf-8')
        argvC[i] = ctypes.create_string_buffer(enc_arg)

    if mpi_comm is None or MPI is None:
        libwarpx.amrex_init(argc, argvC)
    else:
        comm_ptr = MPI._addressof(mpi_comm)
        comm_val = _MPI_Comm_type.from_address(comm_ptr)
        libwarpx.amrex_init_with_inited_mpi(argc, argvC, comm_val)

def initialize(argv=None, mpi_comm=None):
    '''

    Initialize WarpX and AMReX. Must be called before
    doing anything else.

    '''
    if argv is None:
        argv = sys.argv
    amrex_init(argv, mpi_comm)
    libwarpx.warpx_ConvertLabParamsToBoost()
    libwarpx.warpx_ReadBCParams()
    if geometry_dim == 'rz':
        libwarpx.warpx_CheckGriddingForRZSpectral()
    libwarpx.warpx_init()


@atexit.register
def finalize(finalize_mpi=1):
    '''

    Call finalize for WarpX and AMReX. Must be called at
    the end of your script.

    '''
    libwarpx.warpx_finalize()
    libwarpx.amrex_finalize(finalize_mpi)


def evolve(num_steps=-1):
    '''

    Evolve the simulation for num_steps steps. If num_steps=-1,
    the simulation will be run until the end as specified in the
    inputs file.

    Parameters
    ----------

    num_steps: int, the number of steps to take

    '''

    libwarpx.warpx_evolve(num_steps);


def getProbLo(direction):
    assert 0 <= direction < dim, 'Inappropriate direction specified'
    return libwarpx.warpx_getProbLo(direction)


def getProbHi(direction):
    assert 0 <= direction < dim, 'Inappropriate direction specified'
    return libwarpx.warpx_getProbHi(direction)


def getCellSize(direction, level=0):
    assert 0 <= direction < 3, 'Inappropriate direction specified'
    assert 0 <= level and level <= libwarpx.warpx_finestLevel(), 'Inappropriate level specified'
    return libwarpx.warpx_getCellSize(direction, level)


#def get_sigma(direction):
#    '''
#
#    Return the 'sigma' PML coefficients for the electric field
#    in a given direction.
#
#    '''
#
#    size = ctypes.c_int(0)
#    data = libwarpx.warpx_getPMLSigma(direction, ctypes.byref(size))
#    arr = np.ctypeslib.as_array(data, (size.value,))
#    arr.setflags(write=1)
#    return arr
#
#
#def get_sigma_star(direction):
#    '''
#
#    Return the 'sigma*' PML coefficients for the magnetic field
#    in the given direction.
#
#    '''
#
#    size = ctypes.c_int(0)
#    data = libwarpx.warpx_getPMLSigmaStar(direction, ctypes.byref(size))
#    arr = np.ctypeslib.as_array(data, (size.value,))
#    arr.setflags(write=1)
#    return arr
#
#
#def compute_pml_factors(lev, dt):
#    '''
#
#    This recomputes the PML coefficients for a given level, using the
#    time step dt. This needs to be called after modifying the coefficients
#    from Python.
#
#    '''
#
#    libwarpx.warpx_ComputePMLFactors(lev, dt)

def add_particles(species_name, x=None, y=None, z=None, ux=None, uy=None, uz=None, w=None,
                  unique_particles=True, **kwargs):
    '''

    A function for adding particles to the WarpX simulation.

    Parameters
    ----------

    species_name     : the species to add the particle to
    x, y, z          : arrays or scalars of the particle positions (default = 0.)
    ux, uy, uz       : arrays or scalars of the particle momenta (default = 0.)
    w                : array or scalar of particle weights (default = 0.)
    unique_particles : whether the particles are unique or duplicated on
                       several processes. (default = True)
    kwargs           : dictionary containing an entry for all the extra particle
                       attribute arrays. If an attribute is not given it will be
                       set to 0.

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

    # --- If the length of the input is zero, then quietly return
    # --- This is not an error - it just means that no particles are to be injected.
    if maxlen == 0:
        return

    # --- Broadcast scalars into appropriate length arrays
    # --- If the parameter was not supplied, use the default value
    if lenx == 1:
        x = np.full(maxlen, (x or 0.), float)
    if leny == 1:
        y = np.full(maxlen, (y or 0.), float)
    if lenz == 1:
        z = np.full(maxlen, (z or 0.), float)
    if lenux == 1:
        ux = np.full(maxlen, (ux or 0.), float)
    if lenuy == 1:
        uy = np.full(maxlen, (uy or 0.), float)
    if lenuz == 1:
        uz = np.full(maxlen, (uz or 0.), float)
    if lenw == 1:
        w = np.full(maxlen, (w or 0.), float)
    for key, val in kwargs.items():
        if np.size(val) == 1:
            kwargs[key] = np.full(maxlen, val, float)

    # --- The -3 is because the comps include the velocites
    nattr = get_nattr_species(species_name) - 3
    attr = np.zeros((maxlen, nattr))
    attr[:,0] = w

    for key, vals in kwargs.items():
        # --- The -3 is because components 1 to 3 are velocities
        attr[:,get_particle_comp_index(species_name, key)-3] = vals

    libwarpx.warpx_addNParticles(
        ctypes.c_char_p(species_name.encode('utf-8')), x.size,
        x, y, z, ux, uy, uz, nattr, attr, unique_particles
    )


def get_particle_count(species_name):
    '''

    This returns the number of particles of the specified species in the
    simulation.

    Parameters
    ----------

        species_name : the species name that the number will be returned for

    Returns
    -------

        An integer count of the number of particles

    '''
    return libwarpx.warpx_getNumParticles(
        ctypes.c_char_p(species_name.encode('utf-8'))
    )


def get_particle_structs(species_name, level):
    '''

    This returns a list of numpy arrays containing the particle struct data
    on each tile for this process. The particle data is represented as a structured
    numpy array and contains the particle 'x', 'y', 'z', 'id', and 'cpu'.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        species_name : the species name that the data will be returned for

    Returns
    -------

        A List of numpy arrays.

    '''

    particles_per_tile = _LP_c_int()
    num_tiles = ctypes.c_int(0)
    data = libwarpx.warpx_getParticleStructs(
        ctypes.c_char_p(species_name.encode('utf-8')), level,
        ctypes.byref(num_tiles), ctypes.byref(particles_per_tile)
    )

    particle_data = []
    for i in range(num_tiles.value):
        arr = _array1d_from_pointer(data[i], _p_dtype, particles_per_tile[i])
        particle_data.append(arr)

    _libc.free(particles_per_tile)
    _libc.free(data)
    return particle_data


def get_particle_arrays(species_name, comp_name, level):
    '''

    This returns a list of numpy arrays containing the particle array data
    on each tile for this process.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        species_name   : the species name that the data will be returned for
        comp_name      : the component of the array data that will be returned.

    Returns
    -------

        A List of numpy arrays.

    '''

    particles_per_tile = _LP_c_int()
    num_tiles = ctypes.c_int(0)
    data = libwarpx.warpx_getParticleArrays(
        ctypes.c_char_p(species_name.encode('utf-8')),
        ctypes.c_char_p(comp_name.encode('utf-8')),
        level, ctypes.byref(num_tiles), ctypes.byref(particles_per_tile)
    )

    particle_data = []
    for i in range(num_tiles.value):
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


def get_particle_x(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle 'x'
    positions on each tile.

    '''
    structs = get_particle_structs(species_name, level)
    if geometry_dim == '3d' or geometry_dim == '2d':
        return [struct['x'] for struct in structs]
    elif geometry_dim == 'rz':
        return [struct['x']*np.cos(theta) for struct, theta in zip(structs, get_particle_theta(species_name))]


def get_particle_y(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle 'y'
    positions on each tile.

    '''
    structs = get_particle_structs(species_name, level)
    if geometry_dim == '3d' or geometry_dim == '2d':
        return [struct['y'] for struct in structs]
    elif geometry_dim == 'rz':
        return [struct['x']*np.sin(theta) for struct, theta in zip(structs, get_particle_theta(species_name))]


def get_particle_r(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle 'r'
    positions on each tile.

    '''
    structs = get_particle_structs(species_name, level)
    if geometry_dim == 'rz':
        return [struct['x'] for struct in structs]
    elif geometry_dim == '3d':
        return [np.sqrt(struct['x']**2 + struct['y']**2) for struct in structs]
    elif geometry_dim == '2d':
        raise Exception('get_particle_r: There is no r coordinate with 2D Cartesian')


def get_particle_z(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle 'z'
    positions on each tile.

    '''
    structs = get_particle_structs(species_name, level)
    if geometry_dim == '3d':
        return [struct['z'] for struct in structs]
    elif geometry_dim == 'rz' or geometry_dim == '2d':
        return [struct['y'] for struct in structs]


def get_particle_id(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle 'id'
    positions on each tile.

    '''
    structs = get_particle_structs(species_name, level)
    return [struct['id'] for struct in structs]


def get_particle_cpu(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle 'cpu'
    positions on each tile.

    '''
    structs = get_particle_structs(species_name, level)
    return [struct['cpu'] for struct in structs]


def get_particle_weight(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle
    weight on each tile.

    '''

    return get_particle_arrays(species_name, 'w', level)


def get_particle_ux(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle
    x momentum on each tile.

    '''

    return get_particle_arrays(species_name, 'ux', level)


def get_particle_uy(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle
    y momentum on each tile.

    '''

    return get_particle_arrays(species_name, 'uy', level)


def get_particle_uz(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle
    z momentum on each tile.

    '''

    return get_particle_arrays(species_name, 'uz', level)


def get_particle_theta(species_name, level=0):
    '''

    Return a list of numpy arrays containing the particle
    theta on each tile.

    '''

    if geometry_dim == 'rz':
        return get_particle_arrays(species_name, 'theta', level)
    elif geometry_dim == '3d':
        return [np.arctan2(struct['y'], struct['x']) for struct in structs]
    elif geometry_dim == '2d':
        raise Exception('get_particle_r: There is no theta coordinate with 2D Cartesian')


def get_particle_comp_index(species_name, pid_name):
    '''

    Get the component index for a given particle attribute. This is useful
    to get the corrent ordering of attributes when adding new particles using
    `add_particles()`.

    Parameters
    ----------

        species_name   : the species name that the data will be returned for
        pid_name       : string that is used to identify the new component

    Returns
    -------

        Integer corresponding to the index of the requested attribute

    '''
    return libwarpx.warpx_getParticleCompIndex(
        ctypes.c_char_p(species_name.encode('utf-8')),
        ctypes.c_char_p(pid_name.encode('utf-8'))
    )


def add_real_comp(species_name, pid_name, comm=True):
    '''

    Add a real component to the particle data array.

    Parameters
    ----------

        species_name   : the species name for which the new component will be added
        pid_name       : string that is used to identify the new component
        comm           : should the component be communicated

    '''
    libwarpx.warpx_addRealComp(
        ctypes.c_char_p(species_name.encode('utf-8')),
        ctypes.c_char_p(pid_name.encode('utf-8')), comm
    )


def _get_boundary_number(boundary):
    '''

    Utility function to find the boundary number given a boundary name.

    Parameters
    ----------

        boundary       : the boundary from which to get the scraped particle data.
                         In the form x/y/z_hi/lo or eb.

    Returns
    -------

    Integer index in the boundary scraper buffer for the given boundary.
    '''
    if geometry_dim == '3d':
        dimensions = {'x' : 0, 'y' : 1, 'z' : 2}
    elif geometry_dim == '2d':
        dimensions = {'x' : 0, 'z' : 1}
    else:
        raise NotImplementedError("RZ is not supported for particle scraping.")

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
        boundary_num = 4 if geometry_dim == '2d' else 6

    return boundary_num


def get_particle_boundary_buffer_size(species_name, boundary):
    '''

    This returns the number of particles that have been scraped so far in the simulation
    from the specified boundary and of the specified species.

    Parameters
    ----------

        species_name   : return the number of scraped particles of this species
        boundary       : the boundary from which to get the scraped particle data.
                         In the form x/y/z_hi/lo

    Returns
    -------

        The number of particles scraped so far from a boundary and of a species.

    '''
    return libwarpx.warpx_getParticleBoundaryBufferSize(
        ctypes.c_char_p(species_name.encode('utf-8')),
        _get_boundary_number(boundary)
    )


def get_particle_boundary_buffer(species_name, boundary, comp_name, level):
    '''

    This returns a list of numpy arrays containing the particle array data
    for a species that has been scraped by a specific simulation boundary.

    The data for the numpy arrays are not copied, but share the underlying
    memory buffer with WarpX. The numpy arrays are fully writeable.

    Parameters
    ----------

        species_name   : the species name that the data will be returned for.
        boundary       : the boundary from which to get the scraped particle data.
                         In the form x/y/z_hi/lo or eb.
        comp_name      : the component of the array data that will be returned.
                         If "step_scraped" the special attribute holding the
                         timestep at which a particle was scraped will be
                         returned.
        level          : Which AMR level to retrieve scraped particle data from.
    Returns
    -------

        A List of numpy arrays.

    '''
    particles_per_tile = _LP_c_int()
    num_tiles = ctypes.c_int(0)
    if comp_name == 'step_scraped':
        data = libwarpx.warpx_getParticleBoundaryBufferScrapedSteps(
            ctypes.c_char_p(species_name.encode('utf-8')),
            _get_boundary_number(boundary), level,
            ctypes.byref(num_tiles), ctypes.byref(particles_per_tile)
        )
    else:
        data = libwarpx.warpx_getParticleBoundaryBuffer(
            ctypes.c_char_p(species_name.encode('utf-8')),
            _get_boundary_number(boundary), level,
            ctypes.byref(num_tiles), ctypes.byref(particles_per_tile),
            ctypes.c_char_p(comp_name.encode('utf-8'))
        )

    particle_data = []
    for i in range(num_tiles.value):
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


def _get_mesh_field_list(warpx_func, level, direction, include_ghosts):
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
    ngvect = [ngrowvect[i] for i in range(dim)]
    grid_data = []
    shapesize = dim
    if ncomps.value > 1:
        shapesize += 1
    for i in range(size.value):
        shape = tuple([shapes[shapesize*i + d] for d in range(shapesize)])
        # --- The data is stored in Fortran order, hence shape is reversed and a transpose is taken.
        arr = np.ctypeslib.as_array(data[i], shape[::-1]).T
        try:
            # This fails on some versions of numpy
            arr.setflags(write=1)
        except ValueError:
            pass
        if include_ghosts:
            grid_data.append(arr)
        else:
            grid_data.append(arr[tuple([slice(ngvect[d], -ngvect[d]) for d in range(dim)])])

    _libc.free(shapes)
    _libc.free(data)
    return grid_data


def get_mesh_electric_field(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getEfield, level, direction, include_ghosts)


def get_mesh_electric_field_cp(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getEfieldCP, level, direction, include_ghosts)


def get_mesh_electric_field_fp(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getEfieldFP, level, direction, include_ghosts)


def get_mesh_electric_field_cp_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_field_list(libwarpx.warpx_getEfieldCP_PML, level, direction, include_ghosts)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_electric_field_fp_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_field_list(libwarpx.warpx_getEfieldFP_PML, level, direction, include_ghosts)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_magnetic_field(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getBfield, level, direction, include_ghosts)


def get_mesh_magnetic_field_cp(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getBfieldCP, level, direction, include_ghosts)


def get_mesh_magnetic_field_fp(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getBfieldFP, level, direction, include_ghosts)


def get_mesh_magnetic_field_cp_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_field_list(libwarpx.warpx_getBfieldCP_PML, level, direction, include_ghosts)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_magnetic_field_fp_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_field_list(libwarpx.warpx_getBfieldFP_PML, level, direction, include_ghosts)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_current_density(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getCurrentDensity, level, direction, include_ghosts)


def get_mesh_current_density_cp(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getCurrentDensityCP, level, direction, include_ghosts)


def get_mesh_current_density_fp(level, direction, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getCurrentDensityFP, level, direction, include_ghosts)


def get_mesh_current_density_cp_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_field_list(libwarpx.warpx_getCurrentDensityCP_PML, level, direction, include_ghosts)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_current_density_fp_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_field_list(libwarpx.warpx_getCurrentDensityFP_PML, level, direction, include_ghosts)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_charge_density_cp(level, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getChargeDensityCP, level, None, include_ghosts)


def get_mesh_charge_density_fp(level, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getChargeDensityFP, level, None, include_ghosts)


def get_mesh_phi_fp(level, include_ghosts=True):
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
    return _get_mesh_field_list(libwarpx.warpx_getPhiFP, level, None, include_ghosts)


def get_mesh_F_cp(level, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getFfieldCP, level, None, include_ghosts)


def get_mesh_F_fp(level, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getFfieldFP, level, None, include_ghosts)


def get_mesh_G_cp(level, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getGfieldCP, level, None, include_ghosts)


def get_mesh_G_fp(level, include_ghosts=True):
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

    return _get_mesh_field_list(libwarpx.warpx_getGfieldFP, level, None, include_ghosts)


def _get_mesh_array_lovects(level, direction, include_ghosts=True, getlovectsfunc=None):
    assert(0 <= level and level <= libwarpx.warpx_finestLevel())

    size = ctypes.c_int(0)
    ngrowvect = _LP_c_int()
    if direction is None:
        data = getlovectsfunc(level, ctypes.byref(size), ctypes.byref(ngrowvect))
    else:
        data = getlovectsfunc(level, direction, ctypes.byref(size), ctypes.byref(ngrowvect))

    lovects_ref = np.ctypeslib.as_array(data, (size.value, dim))

    # --- Make a copy of the data to avoid memory problems
    # --- Also, take the transpose to give shape (dims, number of grids)
    lovects = lovects_ref.copy().T

    ng = []
    if include_ghosts:
        for d in range(dim):
            ng.append(ngrowvect[d])
    else:
        for d in range(dim):
            ng.append(0)
            lovects[d,:] += ngrowvect[d]

    del lovects_ref
    _libc.free(data)
    return lovects, ng


def get_mesh_electric_field_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getEfieldLoVects)


def get_mesh_electric_field_cp_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getEfieldCPLoVects)


def get_mesh_electric_field_fp_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getEfieldFPLoVects)


def get_mesh_electric_field_cp_lovects_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getEfieldCPLoVects_PML)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_electric_field_fp_lovects_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getEfieldFPLoVects_PML)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_magnetic_field_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getBfieldLoVects)


def get_mesh_magnetic_field_cp_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getBfieldCPLoVects)


def get_mesh_magnetic_field_fp_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getBfieldFPLoVects)


def get_mesh_magnetic_field_cp_lovects_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getBfieldCPLoVects_PML)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_magnetic_field_fp_lovects_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getBfieldFPLoVects_PML)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_current_density_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getCurrentDensityLoVects)


def get_mesh_current_density_cp_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getCurrentDensityCPLoVects)

def get_mesh_current_density_fp_lovects(level, direction, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getCurrentDensityFPLoVects)


def get_mesh_current_density_cp_lovects_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getCurrentDensityCPLoVects_PML)
    except ValueError:
        raise Exception('PML not initialized')

def get_mesh_current_density_fp_lovects_pml(level, direction, include_ghosts=True):
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
        return _get_mesh_array_lovects(level, direction, include_ghosts, libwarpx.warpx_getCurrentDensityFPLoVects_PML)
    except ValueError:
        raise Exception('PML not initialized')


def get_mesh_charge_density_cp_lovects(level, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, None, include_ghosts, libwarpx.warpx_getChargeDensityCPLoVects)


def get_mesh_charge_density_fp_lovects(level, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, None, include_ghosts, libwarpx.warpx_getChargeDensityFPLoVects)


def get_mesh_phi_fp_lovects(level, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, None, include_ghosts, libwarpx.warpx_getPhiFPLoVects)


def get_mesh_F_cp_lovects(level, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, None, include_ghosts, libwarpx.warpx_getFfieldCPLoVects)


def get_mesh_F_fp_lovects(level, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, None, include_ghosts, libwarpx.warpx_getFfieldFPLoVects)


def get_mesh_G_cp_lovects(level, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, None, include_ghosts, libwarpx.warpx_getGfieldCPLoVects)


def get_mesh_G_fp_lovects(level, include_ghosts=True):
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
    return _get_mesh_array_lovects(level, None, include_ghosts, libwarpx.warpx_getGfieldFPLoVects)


def _get_nodal_flag(getdatafunc):
    data = getdatafunc()
    nodal_flag_ref = np.ctypeslib.as_array(data, (dim,))

    # --- Make a copy of the data to avoid memory problems
    nodal_flag = nodal_flag_ref.copy()

    del nodal_flag_ref
    _libc.free(data)
    return nodal_flag


def get_Ex_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Ex along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getEx_nodal_flag)


def get_Ey_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Ey along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getEy_nodal_flag)


def get_Ez_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Ez along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getEz_nodal_flag)


def get_Bx_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Bx along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getBx_nodal_flag)


def get_By_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for By along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getBy_nodal_flag)


def get_Bz_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Bz along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getBz_nodal_flag)


def get_Jx_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Jx along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getJx_nodal_flag)


def get_Jy_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Jy along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getJy_nodal_flag)


def get_Jz_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Jz along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getJz_nodal_flag)


def get_Rho_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Rho along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getRho_nodal_flag)

def get_Phi_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for Phi along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getPhi_nodal_flag)

def get_F_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for F along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getF_nodal_flag)

def get_G_nodal_flag():
    '''
    This returns a 1d array of the nodal flags for G along each direction. A 1 means node centered, and 0 cell centered.
    '''
    return _get_nodal_flag(libwarpx.warpx_getG_nodal_flag)
