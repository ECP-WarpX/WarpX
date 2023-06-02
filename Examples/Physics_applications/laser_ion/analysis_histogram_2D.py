#!/usr/bin/env python3

# This script displays a 2D histogram.

import argparse

import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

ts = OpenPMDTimeSeries("hist2D")

it = ts.iterations
data, info = ts.get_field(field = "data", iteration = 0, plot = True)
print('The available iterations of the simulation are:', it)
print('The axis of the histogram are (0: ordinate ; 1: abscissa):', info.axes)
print('The data shape is:', data.shape)

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument("iter", help="Iteration number of the simulation that is plotted. Enter a number from the list of iterations or 'All' if you want all plots.")
args = parser.parse_args()

if args.iter == 'All' :
    for i in it :
        plt.figure()
        ts.get_field(field = "data", iteration = i, plot = True)
        plt.savefig('Histogram_2D_iteration_' + str(i) + '.pdf')
else :
    i = int(args.iter)
    plt.figure()
    ts.get_field(field = "data", iteration = i, plot = True)
    plt.savefig('Histogram_2D_iteration_' + str(i) + '.pdf')
