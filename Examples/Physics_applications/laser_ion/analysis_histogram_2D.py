#!/usr/bin/env python3

# This script displays a 2D histogram.

import argparse

import matplotlib.colors as colors
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

# Add the simulation time to the title once this information
# is available in the "info" FieldMetaInformation object.
if args.iter == 'All' :
    for i in it :
        plt.figure()
        data, info = ts.get_field(field = "data", iteration = i, plot = False)
        plt.imshow(data, aspect="auto", norm=colors.LogNorm(), extent=info.imshow_extent)
        plt.title(args.hist2D + " (iteration %d)" % i)
        plt.xlabel(info.axes[1] + ' (m)')
        plt.ylabel(info.axes[0] + ' (m.s-1)')
        plt.colorbar()
        plt.savefig('Histogram_2D_' + args.hist2D + '_iteration_' + str(i) + '.pdf')
else :
    i = int(args.iter)
    plt.figure()
    data, info = ts.get_field(field = "data", iteration = i, plot = False)
    plt.imshow(data, aspect="auto", norm=colors.LogNorm(), extent=info.imshow_extent)
    plt.title(args.hist2D + " (iteration %d)" % i)
    plt.xlabel(info.axes[1] + ' (m)')
    plt.ylabel(info.axes[0] + ' (m.s-1)')
    plt.colorbar()
    plt.savefig('Histogram_2D_' + args.hist2D + '_iteration_' + str(i) + '.pdf')
