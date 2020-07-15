#! /usr/bin/env python

# Copyright 2020 Michael Rowan
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL


import sys
import yt
yt.funcs.mylog.setLevel(50)

sys.path.insert(1, '../../../../warpx/Regression/Checksum/')
import checksumAPI

# Name of the plotfile
fn = sys.argv[1]

test_name = fn[:-9] # Could also be os.path.split(os.getcwd())[1]

checksumAPI.evaluate_checksum(test_name, fn, rtol=1e-4, do_particles=False)
