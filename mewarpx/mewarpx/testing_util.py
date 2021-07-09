"""
Utility functions for pytest in mewarpx. Make sure to import after
init_libwarpx is run!
"""
import os
import random

from mewarpx import util as mwxutil
from mewarpx.mwxrun import mwxrun
from mpi4py import MPI


def initialize_testingdir(name):
    """Change into an appropriate directory for testing. This means placing it
    in util.temp_dir, and attaching "_XXXXXX", a random 6-digit integer, to the
    end of the path.
    Arguments:
        name (str): The base string used in the new directory.
    Returns:
        origwd, newdir (str, str): Original working directory, and current
        working directory.
    """
    wd = None
    #TODO: change this to check mwxrun.me once things are merged
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        # Use standard python random here, in case numpy rseed is being
        # used that might fix this randint call.
        wd = os.path.join(
            mwxutil.mewarpx_dir,
            '../tests/temp_files',
            name + "_{:06d}".format(random.randint(0, 1000000))
        )

    if comm.Get_size() > 1:
        wd = comm.bcast(wd, root=0)

    if wd is None:
        raise ValueError("Working directory not properly set or broadcast.")

    origwd = change_to_warpxdir(wd)

    return origwd, os.getcwd()


def change_to_warpxdir(wd):
    """Handle logic of moving to correct directory. The directory is created if
    needed.
    Arguments:
        wd (str): Path of desired working directory.
    Returns:
        origwd (str): Path of original working directory
    """
    origwd = os.getcwd()
    #TODO: change this to check mwxrun.me once things are merged
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    if rank == 0:
        mwxutil.mkdir_p(wd)

    comm.Barrier()

    os.chdir(wd)
    print("Change to working directory {}".format(os.getcwd()))

    comm.Barrier()

    return origwd
