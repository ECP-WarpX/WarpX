#!/usr/bin/env python3

# This script displays a 2D histogram.

import matplotlib.pyplot as plt
from openpmd_viewer import OpenPMDTimeSeries

ts = OpenPMDTimeSeries("hist2D")

it = ts.iterations
data, info = ts.get_field(field = "data", iteration = 0, plot = True)
print('The available iterations of the simulation are:', it)
print('The axis of the histogram are (0: ordinate ; 1: abscissa):', info.axes)
x, y = data.shape
print('The data shape is:', (x - 1, y - 1))

iter = input("Which iteration do you want to plot? [enters a number from the list above or 'All' if you want all plots]")
if iter == 'All' :
    for i in it :
        fig = plt.figure()
        data, info = ts.get_field(field = "data", iteration = i, plot = True)
        plt.savefig('Histogram_2D_iteration_' + str(i) + '.pdf')
else :
    i = int(iter)
    fig = plt.figure()
    data, info = ts.get_field(field = "data", iteration = i, plot = True)
    plt.savefig('Histogram_2D_iteration_' + str(i) + '.pdf')
