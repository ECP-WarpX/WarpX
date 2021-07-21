#! /usr/bin/env python

# Copyright 2021 Modern Electron

# This script just checks that the PICMI file executed successfully.
# If it did there will be a plotfile for the final step.

import yt

plotfile = 'Python_particle_attr_access_plt00010'
ds = yt.load( plotfile )

assert True
