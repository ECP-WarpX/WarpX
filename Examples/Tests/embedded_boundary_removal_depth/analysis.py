#!/usr/bin/env python

"""
This analysis script checks for any spurious charge build-up at the embedded boundary, when particles are removed in 3D.
It averages the divergence of the electric field (divE) over the last 20 time steps and compares the results with a specified tolerance.
"""

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import tqdm
from openpmd_viewer import OpenPMDTimeSeries

# yt.funcs.mylog.setLevel(0)
sys.path.insert(1, "../../../../warpx/Regression/Checksum/")
import checksumAPI

# Open plotfile specified in command line
filename = sys.argv[1]
test_name = os.path.split(os.getcwd())[1]
checksumAPI.evaluate_checksum(test_name, filename, output_format="openpmd")
print(os.getcwd())


def get_avg_divE(ts, start_avg_iter, end_avg_iter, ar_size):
    avg_divE = np.zeros((ar_size, ar_size))
    for iteration in tqdm.tqdm(ts.iterations[start_avg_iter:end_avg_iter]):
        divE = ts.get_field(
            "divE", iteration=iteration, slice_across="y", plot=False, cmap="RdBu"
        )
        avg_divE += divE[0]
    return avg_divE / (end_avg_iter - start_avg_iter)


def plot(array, vmax=1e-9):
    x = np.linspace(-7, 7, 400)
    y = np.sqrt(7**2 - x**2)

    fig, ax = plt.subplots()
    ax.plot(x, y, "k--")
    ax.plot(x, -y, "k--")
    cax = ax.imshow(array, cmap="RdBu", extent=[-10, 10, -10, 10], origin="lower")
    cbar = fig.colorbar(cax, ax=ax)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title("Averaged divE")


ts = OpenPMDTimeSeries("./diags/diag1/")

ar_size = 32
start_avg_iter = 20
end_avg_iter = 50

divE_avg = get_avg_divE(ts, start_avg_iter, end_avg_iter, ar_size)
plot(divE_avg, vmax=1e-9)
plt.savefig("AverageddivE.png")

tolerance = 1e-9


def check_tolerance(array, tolerance):
    assert np.all(
        array <= tolerance
    ), f"Test did not pass: one or more elements exceed the tolerance of {tolerance}."
    print("All elements of are within the tolerance.")


check_tolerance(divE_avg, tolerance)
