from mpi4py import MPI as mpi
import numpy as np


from mewarpx.mwxrun import mwxrun


comm_world = mpi.COMM_WORLD

def mpiallreduce(data=None, op=None, comm=None):
    if op is None:
        op = mpi.SUM
    if comm is None:
        comm = comm_world
    # --- "fast" version was removed because it produced bugs
    if isinstance(data, np.ndarray) and data.dtype is not np.dtype('object'):
        result = np.empty_like(data)
        comm.Allreduce(data, result, op=op)
    else:
        result = comm.allreduce(data, op=op)

    return result

def parallelsum(a, comm=None):
    if mwxrun.n_procs <= 1:
        return a
    return mpiallreduce(a, op=mpi.SUM, comm=comm)
