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
         Whether to include the ghost cells
    """
    def __init__(self, mf_name, level, include_ghosts=False):
        self.mf_name = mf_name
        self.level = level
        self.include_ghosts = include_ghosts

        warpx = libwarpx.libwarpx_so.get_instance()
        self.mf = warpx.multifab(self.mf_name)
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

    def _find_start_stop(self, ii, imin, imax, overlap, d):
        """Given the input index, calculate the start and stop range of the indices.

        Parameters
        ----------
        ii: slice or integer
            Input index, either a slice object or an integer

        imin: integer
            The global lowest index value in the specified direction

        imax: integer
            The global highest index value in the specified direction

        overlap: integer
            Amount neighboring domains overlap

        d: integer
            The dimension number, 1, 2, 3, or 4 (4 being the components)

        If ii is a slice, the start and stop values are used directly,
        unless they are None, then the lower or upper bound is used.
        An assertion checks if the indices are within the bounds.
        """
        if isinstance(ii, slice):
            if ii.start is None:
                iistart = imin
            else:
                iistart = ii.start
            if ii.stop is None:
                iistop = imax + overlap
            else:
                iistop = ii.stop
        else:
            iistart = ii
            iistop = ii + 1
        assert imin <= iistart <= imax + overlap, Exception(f'Dimension {d} lower index is out of bounds')
        assert imin <= iistop <= imax + overlap, Exception(f'Dimension {d} upper index is out of bounds')
        return iistart, iistop

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
            return missing, missing, index[0]
        elif self.dim == 2:
            return index[0], missing, index[1]
        elif self.dim == 3:
            return index[0], index[1], index[2]

    def _get_min_indices(self):
        "Returns the minimum indices, expanded to length 3"
        min_box = self.mf.box_array().minimal_box()
        return self._get_indices(min_box.small_end, 0)

    def _get_max_indices(self):
        "Returns the maximum indices, expanded to length 3"
        min_box = self.mf.box_array().minimal_box()
        return self._get_indices(min_box.big_end, 1)

    def _get_intersect_slice(self, box, starts, stops, ic):
        """Return the slices where the block intersects with the global slice.
        If the block does not intersect, return None.

        Parameters
        ----------
        box: Box instance
            The box defining the block

        starts: sequence
            The minimum indices of the global slice

        stops: sequence
            The maximum indices of the global slice

        ic: integer
            The indices of the components

        Returns
        -------
        block_slices:
            The slice of the intersection relative to the block

        global_slices:
            The slice of the intersection relative to the global array where the data from individual block will go
        """
        ilo = self._get_indices(box.small_end, 0)
        ihi = self._get_indices(box.big_end, 1)
        i1 = np.maximum(starts, ilo)
        i2 = np.minimum(stops, ihi)

        if np.all(i1 < i2):

            block_slices = []
            global_slices = []
            for i in range(3):
                block_slices.append(slice(i1[i] - ilo[i], i2[i] - ilo[i]))
                global_slices.append(slice(i1[i] - starts[i], i2[i] - starts[i]))

            block_slices.append(ic)

            return tuple(block_slices), tuple(global_slices)
        else:
            return None, None

    def __getitem__(self, index):
        """Returns slice of the MultiFab using global indexing.
        The shape of the object returned depends on the number of ix, iy and iz specified, which
        can be from none to all three. Note that the values of ix, iy and iz are
        relative to the fortran indexing, meaning that 0 is the lower boundary
        of the whole domain, and in fortran ordering, i.e. [ix,iy,iz].

        Parameters
        ----------
        index: integer, or sequence of integers or slices, or Ellipsis
            Index of the slice to return
        """
        if index == Ellipsis:
            index = self.dim*[slice(None)]

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

        ixmin, iymin, izmin = self._get_min_indices()
        ixmax, iymax, izmax = self._get_max_indices()

        # Setup the size of the array to be returned
        ixstart, ixstop = self._find_start_stop(ii[0], ixmin, ixmax, self.overlaps[0], 1)
        iystart, iystop = self._find_start_stop(ii[1], iymin, iymax, self.overlaps[1], 2)
        izstart, izstop = self._find_start_stop(ii[2], izmin, izmax, self.overlaps[2], 3)
        icstart, icstop = self._find_start_stop(ic, 0, self.mf.n_comp(), 0, 4)

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
                mf_arr = np.array(self.mf.array(mfi), copy=False).T
                datalist.append((global_slices, mf_arr[block_slices]))

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
        result_global = np.zeros(result_shape, dtype=all_datalist[0][0][1].dtype)

        # Now, copy the data into the result array
        for datalist in all_datalist:
            for global_slices, f_arr in datalist:
                result_global[global_slices] = f_arr

        # Remove dimensions of length 1, and if all dimensions
        # are removed, return a scalar (that's what the [()] does)
        return result_global.squeeze()[()]

    def __setitem__(self, index, value):
        """Sets slices of a decomposed array.
        The shape of the input object depends on the number of arguments specified, which can
        be from none to all three.

        Parameters
        ----------
        index: integer, or sequence of integers or slices, or Ellipsis
            The slice to set

        value: scalar or array
            Input value to assign to the specified slice of the MultiFab
        """
        if index == Ellipsis:
            index = tuple(self.dim*[slice(None)])

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

        ixmin, iymin, izmin = self._get_min_indices()
        ixmax, iymax, izmax = self._get_max_indices()

        # Setup the size of the global array to be set
        ixstart, ixstop = self._find_start_stop(ii[0], ixmin, ixmax, self.overlaps[0], 1)
        iystart, iystop = self._find_start_stop(ii[1], iymin, iymax, self.overlaps[1], 2)
        izstart, izstop = self._find_start_stop(ii[2], izmin, izmax, self.overlaps[2], 3)
        icstart, icstop = self._find_start_stop(ic, 0, self.mf.n_comp(), 0, 4)

        if isinstance(value, np.ndarray):
            # Expand the shape of the input array to match the shape of the global array
            # (it needs to be 4-D).
            # The shape of 1 is added for the extra dimensions and when index is an integer
            # (in which case the dimension was not in the input array).
            value3d = np.array(value, copy=False)
            global_shape = list(value3d.shape)
            if not isinstance(ii[0], slice): global_shape[0:0] = [1]
            if not isinstance(ii[1], slice): global_shape[1:1] = [1]
            if not isinstance(ii[2], slice): global_shape[2:2] = [1]
            if not isinstance(ic   , slice): global_shape[3:3] = [1]
            value3d.shape = global_shape

        starts = [ixstart, iystart, izstart]
        stops = [ixstop, iystop, izstop]
        for mfi in self.mf:
            box = mfi.tilebox()
            block_slices, global_slices = self._get_intersect_slice(box, starts, stops, ic)
            if global_slices is not None:
                mf_arr = np.array(self.mf.array(mfi), copy=False).T
                if isinstance(value, np.ndarray):
                    mf_arr[block_slices] = value3d[global_slices]
                else:
                    mf_arr[block_slices] = value


def ExWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_aux[x][l={level}]', level=level, include_ghosts=include_ghosts)

def EyWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_aux[y][l={level}]', level=level, include_ghosts=include_ghosts)

def EzWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_aux[z][l={level}]', level=level, include_ghosts=include_ghosts)

def BxWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_aux[x][l={level}]', level=level, include_ghosts=include_ghosts)

def ByWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_aux[y][l={level}]', level=level, include_ghosts=include_ghosts)

def BzWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_aux[z][l={level}]', level=level, include_ghosts=include_ghosts)

def JxWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[x][l={level}]', level=level, include_ghosts=include_ghosts)

def JyWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[y][l={level}]', level=level, include_ghosts=include_ghosts)

def JzWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[z][l={level}]', level=level, include_ghosts=include_ghosts)

def ExFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_fp[x][l={level}]', level=level, include_ghosts=include_ghosts)

def EyFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_fp[y][l={level}]', level=level, include_ghosts=include_ghosts)

def EzFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_fp[z][l={level}]', level=level, include_ghosts=include_ghosts)

def BxFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_fp[x][l={level}]', level=level, include_ghosts=include_ghosts)

def ByFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_fp[y][l={level}]', level=level, include_ghosts=include_ghosts)

def BzFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_fp[z][l={level}]', level=level, include_ghosts=include_ghosts)

def JxFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[x][l={level}]', level=level, include_ghosts=include_ghosts)

def JyFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[y][l={level}]', level=level, include_ghosts=include_ghosts)

def JzFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp[z][l={level}]', level=level, include_ghosts=include_ghosts)

def RhoFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'rho_fp[l={level}]', level=level, include_ghosts=include_ghosts)

def PhiFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'phi_fp[l={level}]', level=level, include_ghosts=include_ghosts)

def FFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'f_fp[l={level}]', level=level, include_ghosts=include_ghosts)

def GFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'g_fp[l={level}]', level=level, include_ghosts=include_ghosts)

def AxFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'vector_potential_fp_nodal[x][l={level}]', level=level, include_ghosts=include_ghosts)

def AyFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'vector_potential_fp_nodal[y][l={level}]', level=level, include_ghosts=include_ghosts)

def AzFPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'vector_potential_fp_nodal[z][l={level}]', level=level, include_ghosts=include_ghosts)

def ExCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_cp[x][l={level}]', level=level, include_ghosts=include_ghosts)

def EyCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_cp[y][l={level}]', level=level, include_ghosts=include_ghosts)

def EzCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Efield_cp[z][l={level}]', level=level, include_ghosts=include_ghosts)

def BxCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_cp[x][l={level}]', level=level, include_ghosts=include_ghosts)

def ByCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_cp[y][l={level}]', level=level, include_ghosts=include_ghosts)

def BzCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'Bfield_cp[z][l={level}]', level=level, include_ghosts=include_ghosts)

def JxCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_cp[x][l={level}]', level=level, include_ghosts=include_ghosts)

def JyCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_cp[y][l={level}]', level=level, include_ghosts=include_ghosts)

def JzCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_cp[z][l={level}]', level=level, include_ghosts=include_ghosts)

def RhoCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'rho_cp[l={level}]', level=level, include_ghosts=include_ghosts)

def FCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'f_cp[l={level}]', level=level, include_ghosts=include_ghosts)

def GCPWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'g_cp[l={level}]', level=level, include_ghosts=include_ghosts)

def EdgeLengthsxWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_edge_lengths[x][l={level}]', level=level, include_ghosts=include_ghosts)

def EdgeLengthsyWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_edge_lengths[y][l={level}]', level=level, include_ghosts=include_ghosts)

def EdgeLengthszWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_edge_lengths[z][l={level}]', level=level, include_ghosts=include_ghosts)

def FaceAreasxWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_face_areas[x][l={level}]', level=level, include_ghosts=include_ghosts)

def FaceAreasyWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_face_areas[y][l={level}]', level=level, include_ghosts=include_ghosts)

def FaceAreaszWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'm_face_areas[z][l={level}]', level=level, include_ghosts=include_ghosts)

def JxFPAmpereWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp_ampere[x][l={level}]', level=level, include_ghosts=include_ghosts)

def JyFPAmpereWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp_ampere[y][l={level}]', level=level, include_ghosts=include_ghosts)

def JzFPAmpereWrapper(level=0, include_ghosts=False):
    return _MultiFABWrapper(mf_name=f'current_fp_ampere[z][l={level}]', level=level, include_ghosts=include_ghosts)

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
