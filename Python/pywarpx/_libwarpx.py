# Copyright 2017-2022 The WarpX Community
#
# This file is part of WarpX. It defines the wrapper functions that directly
# call the underlying compiled routines through pybind11.
#
# Authors: Axel Huebl, Andrew Myers, David Grote, Remi Lehe, Weiqun Zhang
#
# License: BSD-3-Clause-LBNL

import atexit
import ctypes
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
                "Invalid attempt to load libwarpx_so... library multiple times."
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
                from . import warpx_pybind_1d as cxx_1d
                self.libwarpx_so = cxx_1d
                self.dim = 1
            elif self.geometry_dim == "2d":
                from . import warpx_pybind_2d as cxx_2d
                self.libwarpx_so = cxx_2d
                self.dim = 2
            elif self.geometry_dim == "rz":
                from . import warpx_pybind_rz as cxx_rz
                self.libwarpx_so = cxx_rz
                self.dim = 2
            elif self.geometry_dim == "3d":
                from . import warpx_pybind_3d as cxx_3d
                self.libwarpx_so = cxx_3d
                self.dim = 3
        except ImportError:
            raise Exception(f"Dimensionality '{self.geometry_dim}' was not compiled in this Python install. Please recompile with -DWarpX_DIMS={_dims}")

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
        if mpi_comm is None or MPI is None:
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
        self.libwarpx_so.convert_lab_params_to_boost()
        self.libwarpx_so.read_BC_params()
        if self.geometry_dim == 'rz':
            self.libwarpx_so.check_gridding_for_RZ_spectral()
        self.warpx = self.libwarpx_so.WarpX()
        self.warpx.initialize_data()
        self.libwarpx_so.execute_python_callback("afterinit");
        self.libwarpx_so.execute_python_callback("particleloader");

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

        self.warpx.evolve(num_steps);

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

        # --- The number of built in attributes
        # --- The three velocities
        built_in_attrs = 3
        if self.geometry_dim == 'rz':
            # --- With RZ, there is also theta
            built_in_attrs += 1

        # --- The number of extra attributes (including the weight)
        nattr = self.get_nattr_species(species_name) - built_in_attrs
        attr = np.zeros((maxlen, nattr), self._numpy_particlereal_dtype)
        attr[:,0] = w

        # --- Note that the velocities are handled separately and not included in attr
        # --- (even though they are stored as attributes in the C++)
        for key, vals in kwargs.items():
            attr[:,self.get_particle_comp_index(species_name, key) - built_in_attrs] = vals

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
        wx = self.libwarpx_so.WarpX.get_instance()
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


libwarpx = LibWarpX()
