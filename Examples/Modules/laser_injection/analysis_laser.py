#! /usr/bin/env python

# Copyright 2019 Andrew Myers, Jean-Luc Vay, Maxence Thevenet
# Remi Lehe, Weiqun Zhang
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

import os
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]

# you can save an image to be displayed on the website
t = np.arange(0.0, 2.0, 0.01)
s = 1 + np.sin(2*np.pi*t)
plt.plot(t, s)
plt.savefig("laser_analysis.png")

# Reset benchmark?
reset = ( os.getenv('CHECKSUM_RESET', 'False').lower() in
          ['true', '1', 't', 'y', 'yes', 'on'] )

# Run checksum regression test or reset
test_name = fn[:-9] # Could also be os.path.split(os.getcwd())[1]
if reset:
    checksumAPI.reset_benchmark(test_name, fn)
else:
    checksumAPI.evaluate_checksum(test_name, fn)
