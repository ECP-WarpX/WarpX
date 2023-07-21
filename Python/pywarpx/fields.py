# Copyright 2017-2023 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

"""Provides wrappers around multiFABs

Available routines:

ExWrapper, EyWrapper, EzWrapper
BxWrapper, ByWrapper, BzWrapper
JxWrapper, JyWrapper, JzWrapper

ExFPWrapper, EyFPWrapper, EzFPWrapper
BxFPWrapper, ByFPWrapper, BzFPWrapper
JxFPWrapper, JyFPWrapper, JzFPWrapper

RhoFPWrapper, PhiFPWrapper
FFPWrapper, GFPWrapper
AxFPWrapper, AyFPWrapper, AzFPWrapper

ExCPWrapper, EyCPWrapper, EzCPWrapper
BxCPWrapper, ByCPWrapper, BzCPWrapper
JxCPWrapper, JyCPWrapper, JzCPWrapper

RhoCPWrapper
FCPWrapper, GCPWrapper

EdgeLengthsxWrapper, EdgeLengthsyWrapper, EdgeLengthszWrapper
FaceAreasxWrapper, FaceAreasyWrapper, FaceAreaszWrapper

ExFPPMLWrapper, EyFPPMLWrapper, EzFPPMLWrapper
BxFPPMLWrapper, ByFPPMLWrapper, BzFPPMLWrapper
JxFPPMLWrapper, JyFPPMLWrapper, JzFPPMLWrapper
JxFPAmpereWrapper, JyFPAmpereWrapper, JzFPAmpereWrapper
FFPPMLWrapper, GFPPMLWrapper

ExCPPMLWrapper, EyCPPMLWrapper, EzCPPMLWrapper
BxCPPMLWrapper, ByCPPMLWrapper, BzCPPMLWrapper
JxCPPMLWrapper, JyCPPMLWrapper, JzCPPMLWrapper
FCPPMLWrapper, GCPPMLWrapper
"""
import numpy as np

try:
    import cupy as cp
except ImportError:
    cp = None

try:
    from mpi4py import MPI as mpi
    comm_world = mpi.COMM_WORLD
    npes = comm_world.Get_size()
except ImportError:
    npes = 1

from ._libwarpx import libwarpx


class _MultiFABWrapper(object):
    """Wrapper around MultiFabs
    This provides a convenient way to query and set data in the MultiFabs.
    The indexing is based on global indices.

    Parameters
    ----------
     mf_name: string
         The name of the MultiFab to be accessed, the tag specified when the
         MultiFab is allocated

     level: int
         The refinement level

     include_ghosts: bool, default=False
         Whether to include the ghost cells.
         Note that when True, the first n-ghost negative indices will refer to the lower
         ghost cells.
    """
    def __init__(self, mf_name, level, include_ghosts=False):
        self.mf_name = mf_name
        self.level = level
        self.include_ghosts = include_ghosts

        self.dim = libwarpx.dim

        # The overlaps list is one along the axes where the grid boundaries overlap the neighboring grid,
        # which is the case with node centering.
        ix_type = self.mf.box_array().ix_type()
        self.overlaps = self._get_indices([int(ix_type.node_centered(i)) for i in range(self.dim)], 0)

    def __len__(self):
        "Returns the number of blocks"
        return self.mf.size

    def __iter__(self):
        "The iteration is over the MultiFab"
        return self.mf.__iter__()

    @property
    def mf(self):
        # Always fetch this anew in case the C++ MultiFab is recreated
        warpx = libwarpx.libwarpx_so.get_instance()
        # All MultiFab names have the level suffix
        return warpx.multifab(f'{self.mf_name}[level={self.level}]')

    @property
    def shape(self):
        """Returns the shape of the global array
        """
        min_box = self.mf.box_array().minimal_box()
        shape = list(min_box.size - min_box.small_end)
        if self.include_ghosts:
            nghosts = self.mf.n_grow_vect()
            shape = [shape[i] + 2*nghosts[i] for i in range(self.dim)]
        shape.append(self.mf.nComp)
        return tuple(shape)

    def mesh(self, direction):
        """Returns the mesh along the specified direction with the appropriate centering.

        Parameters
        ----------
        direction: string
            In 3d, one of 'x', 'y', or 'z'.
            In 2d, Cartesian, one of 'x', or 'z'.
            In RZ, one of 'r', or 'z'
            In Z, 'z'.
        """

        try:
            if libwarpx.geometry_dim == '3d':
                idir = ['x', 'y', 'z'].index(direction)
                celldir = idir
            elif libwarpx.geometry_dim == '2d':
                idir = ['x', 'z'].index(direction)
                celldir = 2*idir
            elif libwarpx.geometry_dim == 'rz':
                idir = ['r', 'z'].index(direction)
                celldir = 2*idir
            elif libwarpx.geometry_dim == '1d':
                idir = ['z'].index(direction)
                celldir = idir
        except ValueError:
            raise Exception('Inappropriate direction given')

        min_box = self.mf.box_array().minimal_box()
        ilo = min_box.small_end[idir]
        ihi = min_box.big_end[idir]

        if self.include_ghosts:
            # The ghost cells are added to the upper and lower end of the global domain.
            nghosts = self.mf.n_grow_vect()
            ilo = list(ilo)
            ihi = list(ihi)
            min_box = self.mf.box_array().minimal_box()
            imax = min_box.big_end
            for i in range(self.dim):
                if ilo[i] == 0:
                    ilo[i] -= nghosts[i]
                if ihi[i] == imax[i]:
                    ihi[i] += nghosts[i]

        # Cell size in the direction
        warpx = libwarpx.libwarpx_so.get_instance()
        dd = warpx.Geom(self.level).data().CellSize(idir)

        # The centering shift
        ix_type = self.mf.box_array().ix_type()
        if ix_type.node_centered(idir):
            # node centered
            shift = 0.
        else:
            # cell centered
            shift = 0.5*dd

        lo = warpx.Geom(self.level).ProbLo(idir)
        return lo + np.arange(ilo,ihi+1)*dd + shift

    def _get_indices(self, index, missing):
        """Expand the index list to length three.

        Parameters
        ----------
        index: sequence of length dims
            The indices for each dim

        missing:
            The value used to fill in the extra dimensions added
        """
        if self.dim == 1:
            return index[0], missing, missing
        elif self.dim == 2:
            return index[0], index[1], missing
        elif self.dim == 3:
            return index[0], index[1], index[2]

    def _get_min_indices(self):
        """Returns the minimum indices, expanded to length 3"""
        min_box = self.mf.box_array().minimal_box()
        imin = self._get_indices(min_box.small_end, 0)
        if self.include_ghosts:
            nghosts = self._get_n_ghosts()
            imin = [imin[i] - nghosts[i] for i in range(3)]
        return imin

    def _get_max_indices(self):
        """Returns the maximum indices, expanded to length 3.
        """
        min_box = self.mf.box_array().minimal_box()
        imax = self._get_indices(min_box.big_end, 0)
        if self.include_ghosts:
            nghosts = self._get_n_ghosts()
            imax = [imax[i] + nghosts[i] for i in range(3)]
        return imax

    def _get_n_ghosts(self):
        nghosts = list(self._get_indices(self.mf.n_grow_vect(), 0))
        # The components always has nghosts = 0
        nghosts.append(0)
        return nghosts

    def _fix_index(self, ii, imax, d):
        """Handle negative index, wrapping them as needed.
        This depends on whether ghost cells are included. When true, the indices are
        shifted by the number of ghost cells before being wrapped.
        """
        nghosts = self._get_n_ghosts()
        if self.include_ghosts:
            ii += nghosts[d]
        if ii < 0:
            ii += imax
        if self.include_ghosts:
            ii -= nghosts[d]
        return ii

    def _find_start_stop(self, ii, imin, imax, d):
        """Given the input index, calculate the start and stop range of the indices.

        Parameters
        ----------
        ii: None, slice, integer
            Input index, either None, a slice object, or an integer.
            Note that ii can be negative.

        imin: integer
            The global lowest lower bound in the specified direction.
            This can include the ghost cells.

        imax: integer
            The global highest upper bound in the specified direction.
            This can include the ghost cells.
            This should be the max index + 1.

        d: integer
            The dimension number, 0, 1, 2, or 3 (3 being the components)

        If ii is a slice, the start and stop values are used directly,
        unless they are None, then the lower or upper bound is used.
        An assertion checks if the indices are within the bounds.
        """
        if ii is None:
            iistart = imin
            iistop = imax
        elif isinstance(ii, slice):
            if ii.start is None:
                iistart = imin
            else:
                iistart = self._fix_index(ii.start, imax, d)
            if ii.stop is None:
                iistop = imax
            else:
                iistop = self._fix_index(ii.stop, imax, d)
        else:
            ii = self._fix_index(ii, imax, d)
            iistart = ii
            iistop = ii + 1
        assert imin <= iistart <= imax, Exception(f'Dimension {d+1} lower index is out of bounds')
        assert imin <= iistop <= imax, Exception(f'Dimension {d+1} upper index is out of bounds')
        return iistart, iistop

    def _get_intersect_slice(self, box, starts, stops, ic):
        """Return the slices where the block intersects with the global slice.
        If the block does not intersect, return None.
        This also shifts the block slices by the number of ghost cells in the
        MultiFab arrays since the arrays include the ghost cells.

        Parameters
        ----------
        box: Box instance
            The box defining the block

        starts: sequence
            The minimum indices of the global slice.
            These can be negative.

        stops: sequence
            The maximum indices of the global slice.
            These can be negative.

        ic: integer
            The indices of the components
            This can be negative.

        Returns
        -------
        block_slices:
            The slice of the intersection relative to the block

        global_slices:
            The slice of the intersection relative to the global array where the data from individual block will go
        """
        ilo = self._get_indices(box.small_end, 0)
        ilo_g = self._get_indices(box.small_end, 0)
        ihi_g = self._get_indices(box.big_end, 0)
        nghosts = self._get_indices(self.mf.n_grow_vect(), 0)
        if self.include_ghosts:
            # The starts and stops could include the ghost cells, so add the ghost cells
            # to the lo and hi that they are compared against (to get i1 and i2 below).
            # The ghost cells are only added to the lower and upper end of the global domain.
            ilo_g = list(ilo_g)
            ihi_g = list(ihi_g)
            min_box = self.mf.box_array().minimal_box()
            imax = self._get_indices(min_box.big_end, 0)
            for i in range(3):
                if ilo_g[i] == 0:
                    ilo_g[i] -= nghosts[i]
                if ihi_g[i] == imax[i]:
                    ihi_g[i] += nghosts[i]
        # Add 1 to the upper end to be consistent with the slicing notation
        ihi_g_p1 = [i + 1 for i in ihi_g]
        i1 = np.maximum(starts, ilo_g)
        i2 = np.minimum(stops, ihi_g_p1)

        if np.all(i1 < i2):

            block_slices = []
            global_slices = []
            for i in range(3):
                block_slices.append(slice(i1[i] - ilo[i] + nghosts[i], i2[i] - ilo[i] + nghosts[i]))
                global_slices.append(slice(i1[i] - starts[i], i2[i] - starts[i]))

            if ic is None:
                ic = slice(None)

            block_slices.append(ic)
            global_slices.append(ic)

            return tuple(block_slices), tuple(global_slices)
        else:
            return None, None

    def __getitem__(self, index):
        """Returns slice of the MultiFab using global indexing.
        The shape of the object returned depends on the number of ix, iy and iz specified, which
        can be from none to all three. Note that the values of ix, iy and iz are
        relative to the fortran indexing, meaning that 0 is the lower boundary
        of the whole domain, and in fortran ordering, i.e. [ix,iy,iz].
        This allows negative indexing, though with ghosts cells included, the first n-ghost negative
        indices will refer to the lower guard cells.

        Parameters
        ----------
        index: integer, or sequence of integers or slices, or Ellipsis
            Index of the slice to return
        """
        # Note that the index can have negative values (which wrap around) and has 1 added to the upper
        # limit using python style slicing
        if index == Ellipsis:
            index = self.dim*[slice(None)]
        elif isinstance(index, slice):
            # If only one slice passed in, it was not wrapped in a list
            index = [index]

        if len(index) < self.dim+1:
            # Add extra dims to index, including for the component.
            # These are the dims left out and assumed to extend over the full size of the dim
            index = list(index)
            while len(index) < self.dim+1:
                index.append(slice(None))
        elif len(index) > self.dim+1:
            raise Exception('Too many indices given')

        # Expand the indices to length 3
        ii = self._get_indices(index, None)
        ic = index[-1]

        # Global extent. These include the ghost cells when include_ghosts is True
        ixmin, iymin, izmin = self._get_min_indices()
        ixmax, iymax, izmax = self._get_max_indices()

        # Setup the size of the array to be returned
        ixstart, ixstop = self._find_start_stop(ii[0], ixmin, ixmax+1, 0)
        iystart, iystop = self._find_start_stop(ii[1], iymin, iymax+1, 1)
        izstart, izstop = self._find_start_stop(ii[2], izmin, izmax+1, 2)
        icstart, icstop = self._find_start_stop(ic, 0, self.mf.n_comp(), 3)

        # Gather the data to be included in a list to be sent to other processes
        starts = [ixstart, iystart, izstart]
        stops = [ixstop, iystop, izstop]
        datalist = []
        for mfi in self.mf:
            box = mfi.tilebox()
            block_slices, global_slices = self._get_intersect_slice(box, starts, stops, ic)
            if global_slices is not None:
                # Note that the array will always have 4 dimensions,
                # the three dimensions plus the components, even when
                # self.dim < 3. The transpose is taken since the thing
                # returned by self.mf.array(mfi) is in C ordering.
                device_arr4 = self.mf.array(mfi)
                if cp is not None:
                    device_arr = cp.array(device_arr4, copy=False).T
                    # Copy the data from the device to the host
                    slice_arr = cp.asnumpy(device_arr[block_slices])
                else:
                    device_arr = np.array(device_arr4, copy=False).T
                    slice_arr = device_arr[block_slices]
                datalist.append((global_slices, slice_arr))

        # Gather the data from all processors
        if npes == 1:
            all_datalist = [datalist]
        else:
            all_datalist = comm_world.allgather(datalist)

        # Create the array to be returned
        result_shape = (max(0, ixstop - ixstart),
                        max(0, iystop - iystart),
                        max(0, izstop - izstart),
                        max(0, icstop - icstart))

        # Now, copy the data into the result array
        result_global = None
        for datalist in all_datalist:
            for global_slices, f_arr in datalist:
                if result_global is None:
                    # Delay allocation to here so that the type can be obtained
                    result_global = np.zeros(result_shape, dtype=f_arr.dtype)
                result_global[global_slices] = f_arr

        if result_global is None:
            # Something went wrong with the index and no data was found. Return an empty array.
            result_global = np.zeros(0)

        # Remove dimensions of length 1, and if all dimensions
        # are removed, return a scalar (that's what the [()] does)
        return result_global.squeeze()[()]

    def __setitem__(self, index, value):
        """Sets slices of a decomposed array.
        The shape of the input object depends on the number of arguments specified, which can
        be from none to all three.
        This allows negative indexing, though with ghosts cells included, the first n-ghost negative
        indices will refer to the lower guard cells.

        Parameters
        ----------
        index: integer, or sequence of integers or slices, or Ellipsis
            The slice to set

        value: scalar or array
            Input value to assign to the specified slice of the MultiFab
        """
        # Note that the index can have negative values (which wrap around) and has 1 added to the upper
        # limit using python style slicing
        if index == Ellipsis:
            index = tuple(self.dim*[slice(None)])
        elif isinstance(index, slice):
            # If only one slice passed in, it was not wrapped in a list
            index = [index]

        if len(index) < self.dim+1:
            # Add extra dims to index, including for the component.
            # These are the dims left out and assumed to extend over the full size of the dim.
            index = list(index)
            while len(index) < self.dim+1:
                index.append(slice(None))
        elif len(index) > self.dim+1:
            raise Exception('Too many indices given')

        # Expand the indices to length 3
        ii = self._get_indices(index, None)
        ic = index[-1]

        # Global extent. These include the ghost cells when include_ghosts is True
        ixmin, iymin, izmin = self._get_min_indices()
        ixmax, iymax, izmax = self._get_max_indices()

        # Setup the size of the global array to be set
        ixstart, ixstop = self._find_start_stop(ii[0], ixmin, ixmax+1, 0)
        iystart, iystop = self._find_start_stop(ii[1], iymin, iymax+1, 1)
        izstart, izstop = self._find_start_stop(ii[2], izmin, izmax+1, 2)
        icstart, icstop = self._find_start_stop(ic, 0, self.mf.n_comp(), 3)

        if isinstance(value, np.ndarray):
            # Expand the shape of the input array to match the shape of the global array
            # (it needs to be 4-D).
            # This converts value to an array if needed, and the [...] grabs a view so
            # that the shape change below doesn't affect value.
            value3d = np.array(value)[...]
            global_shape = list(value3d.shape)
            # The shape of 1 is added for the extra dimensions and when index is an integer
            # (in which case the dimension was not in the input array).
            if not isinstance(ii[0], slice): global_shape[0:0] = [1]
            if not isinstance(ii[1], slice): global_shape[1:1] = [1]
            if not isinstance(ii[2], slice): global_shape[2:2] = [1]
            if not isinstance(ic   , slice) or len(global_shape) < 4: global_shape[3:3] = [1]
            value3d.shape = global_shape

        starts = [ixstart, iystart, izstart]
        stops = [ixstop, iystop, izstop]
        for mfi in self.mf:
            box = mfi.tilebox()
            block_slices, global_slices = self._get_intersect_slice(box, starts, stops, ic)
            if global_slices is not None:
                if cp is not None:
                    mf_arr = cp.array(self.mf.array(mfi), copy=False).T
                else:
                    mf_arr = np.array(self.mf.array(mfi), copy=False).T
                if isinstance(value, np.ndarray):
                    slice_value = value3d[global_slices]
                    if cp is not None:
                        # Copy data from host to device
                        slice_value = cp.asarray(value3d[global_slices])
                    mf_arr[block_slices] = slice_value
                else:
                    mf_arr[block_slices] = value

    def min(self, *args):
        return self.mf.min(*args)

    def max(self, *args):
        return self.mf.max(*args)

    def sum(self, *args):
        return self.mf.sum(*args)

    def min_index(self, *args):
        return self.mf.minIndex(*args)

    def max_index(self, *args):
        return self.mf.maxIndex(*args)

    def norm0(self, *args):
        return self.mf.norm0(*args)

def ExWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_aux[x]', level=level, include_ghosts=include_ghosts)

def EyWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_aux[y]', level=level, include_ghosts=include_ghosts)

def EzWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_aux[z]', level=level, include_ghosts=include_ghosts)

def BxWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_aux[x]', level=level, include_ghosts=include_ghosts)

def ByWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_aux[y]', level=level, include_ghosts=include_ghosts)

def BzWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_aux[z]', level=level, include_ghosts=include_ghosts)

def JxWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[x]', level=level, include_ghosts=include_ghosts)

def JyWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[y]', level=level, include_ghosts=include_ghosts)

def JzWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[z]', level=level, include_ghosts=include_ghosts)

def ExFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_fp[x]', level=level, include_ghosts=include_ghosts)

def EyFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_fp[y]', level=level, include_ghosts=include_ghosts)

def EzFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_fp[z]', level=level, include_ghosts=include_ghosts)

def BxFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_fp[x]', level=level, include_ghosts=include_ghosts)

def ByFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_fp[y]', level=level, include_ghosts=include_ghosts)

def BzFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_fp[z]', level=level, include_ghosts=include_ghosts)

def JxFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[x]', level=level, include_ghosts=include_ghosts)

def JyFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[y]', level=level, include_ghosts=include_ghosts)

def JzFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[z]', level=level, include_ghosts=include_ghosts)

def RhoFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'rho_fp', level=level, include_ghosts=include_ghosts)

def PhiFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'phi_fp', level=level, include_ghosts=include_ghosts)

def FFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'f_fp', level=level, include_ghosts=include_ghosts)

def GFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'g_fp', level=level, include_ghosts=include_ghosts)

def AxFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'vector_potential_fp_nodal[x]', level=level, include_ghosts=include_ghosts)

def AyFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'vector_potential_fp_nodal[y]', level=level, include_ghosts=include_ghosts)

def AzFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'vector_potential_fp_nodal[z]', level=level, include_ghosts=include_ghosts)

def ExCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_cp[x]', level=level, include_ghosts=include_ghosts)

def EyCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_cp[y]', level=level, include_ghosts=include_ghosts)

def EzCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_cp[z]', level=level, include_ghosts=include_ghosts)

def BxCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_cp[x]', level=level, include_ghosts=include_ghosts)

def ByCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_cp[y]', level=level, include_ghosts=include_ghosts)

def BzCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_cp[z]', level=level, include_ghosts=include_ghosts)

def JxCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_cp[x]', level=level, include_ghosts=include_ghosts)

def JyCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_cp[y]', level=level, include_ghosts=include_ghosts)

def JzCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_cp[z]', level=level, include_ghosts=include_ghosts)

def RhoCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'rho_cp', level=level, include_ghosts=include_ghosts)

def FCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'f_cp', level=level, include_ghosts=include_ghosts)

def GCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'g_cp', level=level, include_ghosts=include_ghosts)

def EdgeLengthsxWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_edge_lengths[x]', level=level, include_ghosts=include_ghosts)

def EdgeLengthsyWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_edge_lengths[y]', level=level, include_ghosts=include_ghosts)

def EdgeLengthszWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_edge_lengths[z]', level=level, include_ghosts=include_ghosts)

def FaceAreasxWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_face_areas[x]', level=level, include_ghosts=include_ghosts)

def FaceAreasyWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_face_areas[y]', level=level, include_ghosts=include_ghosts)

def FaceAreaszWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_face_areas[z]', level=level, include_ghosts=include_ghosts)

def JxFPAmpereWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp_ampere[x]', level=level, include_ghosts=include_ghosts)

def JyFPAmpereWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp_ampere[y]', level=level, include_ghosts=include_ghosts)

def JzFPAmpereWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp_ampere[z]', level=level, include_ghosts=include_ghosts)

def ExFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_E_fp[x]', level=level, include_ghosts=include_ghosts)

def EyFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_E_fp[y]', level=level, include_ghosts=include_ghosts)

def EzFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_E_fp[z]', level=level, include_ghosts=include_ghosts)

def BxFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_B_fp[x]', level=level, include_ghosts=include_ghosts)

def ByFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_B_fp[y]', level=level, include_ghosts=include_ghosts)

def BzFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_B_fp[z]', level=level, include_ghosts=include_ghosts)

def JxFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_j_fp[x]', level=level, include_ghosts=include_ghosts)

def JyFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_j_fp[y]', level=level, include_ghosts=include_ghosts)

def JzFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_j_fp[z]', level=level, include_ghosts=include_ghosts)

def FFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_f_fp[z]', level=level, include_ghosts=include_ghosts)

def GFPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_g_fp[z]', level=level, include_ghosts=include_ghosts)

def ExCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_E_cp[x]', level=level, include_ghosts=include_ghosts)

def EyCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_E_cp[y]', level=level, include_ghosts=include_ghosts)

def EzCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_E_cp[z]', level=level, include_ghosts=include_ghosts)

def BxCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_B_cp[x]', level=level, include_ghosts=include_ghosts)

def ByCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_B_cp[y]', level=level, include_ghosts=include_ghosts)

def BzCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_B_cp[z]', level=level, include_ghosts=include_ghosts)

def JxCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_j_cp[x]', level=level, include_ghosts=include_ghosts)

def JyCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_j_cp[y]', level=level, include_ghosts=include_ghosts)

def JzCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_j_cp[z]', level=level, include_ghosts=include_ghosts)

def FCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_f_cp[z]', level=level, include_ghosts=include_ghosts)

def GCPPMLWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name='pml_g_cp[z]', level=level, include_ghosts=include_ghosts)
