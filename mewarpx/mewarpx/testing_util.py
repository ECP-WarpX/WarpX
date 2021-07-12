"""
Utility functions for pytest in mewarpx. Make sure to import after
init_libwarpx is run!
"""
import os
import random
import pandas
import numpy as np

from mewarpx import util as mwxutil
from mewarpx.mwxrun import mwxrun
from mpi4py import MPI

test_dir = os.path.join(mwxutil.mewarpx_dir, "../tests/test_files")
temp_dir = os.path.join(mwxutil.mewarpx_dir, "../tests/temp_files")

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
            temp_dir,
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
    print(f"Change to working directory {os.getcwd()}")

    comm.Barrier()

    return origwd


def get_min_wmargin(array, margin=0.1):
    """Get a minimum value with a multiplicative margin of safety.

    Arguments:
        array (np.ndarray): Array to find minimum of
        margin (float): Multiplicative safety margin

    Returns:
        minwmargin (float): If min < 0, this is min*(1-margin). If min > 0, this
        is min*(1+margin).
    """
    amin = np.min(array)
    return amin*(1.0 - np.sign(amin)*margin)


def get_max_wmargin(array, margin=0.1):
    """Get a maximum value with a multiplicative margin of safety.

    Arguments:
        array (np.ndarray): Array to find maximum of
        margin (float): Multiplicative safety margin

    Returns:
        maxwmargin (float): If max < 0, this is max*(1-margin). If max > 0, this
        is max*(1+margin).
    """
    amax = np.max(array)
    return amax*(1.0 + np.sign(amax)*margin)


def test_df_vs_ref(testname, df, suffix=None, margin=0.1):
    """Handle loading and test a reference DataFrame against test results.

    This creates or appends to a ``<testname>.csv`` or
    ``<testname>_<suffix>.csv`` file in ``tests/temp_files/Result_Stats``. It
    then compares the results to a dataframe stored in
    ``tests/test_files/Result_Stats`` under the same filename. This means that
    running and failing the test many times will generate a csv file with all
    the information needed to check that stats don't change in the future. See
    e.g.  ``test_runsuite_vweights.py::test_runsuite_vweights2diag``.

    Arguments:
        testname (str): Name of the test
        df (pandas.DataFrame): DataFrame containing the test results to compare
            against a reference
        suffix (str): An append to test names, allowing multiple DataFrames to
            be referenced in a single test.
        margin (float): A multiplicative safety margin: The tested df values
        should be within `(1 +- margin) * [min or max](ref_val)`. Default 0.1.

    Returns:
        equal_within_margin (bool): True if equal within the margin. False if
        reference file does not exist or data frames do not agree, meaning that
        the test value does not fall within min and max of reference values
        within some error margin.
    """
    csvname = (
        '{}.csv'.format(testname) if suffix is None
        else '{}_{}.csv'.format(testname, suffix)
    )

    # Record for future use if needed
    mwxutil.mkdir_p(os.path.join(temp_dir, 'Result_Stats'))
    record_file = os.path.join(temp_dir, 'Result_Stats', csvname)
    if os.path.exists(record_file):
        df.to_csv(record_file, mode='a', header=False, index=False)
    else:
        df.to_csv(record_file, mode='w', header=True, index=False)

    # Load authoritative records
    record_file = os.path.join(test_dir, 'Result_Stats', csvname)

    if os.path.isfile(record_file):
        refdf = pandas.read_csv(record_file)
        return test_df_in_margin(testdf=df, refdf=refdf, margin=margin)

    else:
        print(f"Record file {record_file} does not exist. Reference data "
              f"frame not loaded.")
        return False


def test_df_in_margin(testdf, refdf, margin=0.1):
    """Test that a DataFrame is within min/max of another, within a margin.

    Arguments:
        testdf (pandas.DataFrame): DataFrame to check
        refdf (pandas.DataFrame): Reference DataFrame
        margin (float): A multiplicative safety margin: The tested df values
        should be within `(1 +- margin) * [min or max](ref_val)`. Default 0.1.

    Returns:
        equal_within_margin (bool): True if equal within the margin
    """
    equal = True
    for column in testdf.columns:
        if column not in refdf.columns:
            raise ValueError("No column {} in reference dataframe".format(
                column))
    for column in refdf.columns:
        testmin = testdf[column].min()
        testmax = testdf[column].max()
        if not testmin >= get_min_wmargin(refdf[column], margin=margin):
            equal = False
            print(
                f"Column {column} with test min {testmin} is not within a "
                f"{margin} margin of reference min {refdf[column].min()}"
            )
        if not testmax <= get_max_wmargin(refdf[column], margin=margin):
            equal = False
            print(
                f"Column {column} with test max {testmax} is not within a "
                f"{margin} margin of reference max {refdf[column].max()}"
            )
    return equal
