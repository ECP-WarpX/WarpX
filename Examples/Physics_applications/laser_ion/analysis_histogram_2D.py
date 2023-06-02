#!/usr/bin/env python3

# This script displays a 2D histogram.

import argparse

import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

parser = argparse.ArgumentParser(description='Process a 2D histogram name and an integer.')
parser.add_argument("hist2D", help="Folder name of the reduced diagnostic.")
parser.add_argument("iter", help="Iteration number of the simulation that is plotted. Enter a number from the list of iterations or 'All' if you want all plots.")
args = parser.parse_args()

path = 'diags/reducedfiles/' + args.hist2D

ts = OpenPMDTimeSeries(path)

it = ts.iterations
data, info = ts.get_field(field = "data", iteration = 0, plot = True)
print('The available iterations of the simulation are:', it)
print('The axis of the histogram are (0: ordinate ; 1: abscissa):', info.axes)
print('The data shape is:', data.shape)

if args.iter == 'All' :
    for i in it :
        plt.figure()
        ts.get_field(field = "data", iteration = i, plot = True)
        plt.savefig('Histogram_2D_' + args.hist2D + '_iteration_' + str(i) + '.pdf')
else :
    i = int(args.iter)
    plt.figure()
    ts.get_field(field = "data", iteration = i, plot = True)
    plt.savefig('Histogram_2D_' + args.hist2D + '_iteration_' + str(i) + '.pdf')
