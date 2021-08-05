#! /usr/bin/env python

# Copyright 2021 Modern Electron
#
# License: BSD-3-Clause-LBNL

# This script just checks that the PICMI file executed successfully.
# If it did there will be a plotfile for the final step.

import yt

plotfile = 'Python_particle_attr_access_plt_000010'
ds = yt.load( plotfile )  # noqa

plotfile = 'Python_particle_attr_access_plt_100010'
ds = yt.load( plotfile )  # noqa

assert True
