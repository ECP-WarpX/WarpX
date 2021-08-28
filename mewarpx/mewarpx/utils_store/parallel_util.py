from mpi4py import MPI as mpi
import numpy as np

comm_world = mpi.COMM_WORLD


def mpiallreduce(data=None, opstring="SUM", comm=None):
    if opstring is None or opstring == "SUM":
        op = mpi.SUM
    elif opstring == "MIN":
        op = mpi.MIN
    else:
        raise NotImplementedError("The opstring is unrecognized or has not been implemented yet.")

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
    from mewarpx.mwxrun import mwxrun
    if mwxrun.n_procs <= 1:
        return a
    return mpiallreduce(a, opstring="SUM", comm=comm)
