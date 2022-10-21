# Copyright 2016-2022 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Andrew Myers, David Grote, Lorenzo Giacomel, Axel Huebl
# License: BSD-3-Clause-LBNL

from .Algo import algo
from .Amr import amr
from .Amrex import amrex
from .Boundary import boundary
from .Collisions import collisions
from .Constants import my_constants
from .Diagnostics import diagnostics, reduced_diagnostics
from .EB2 import eb2
from .Geometry import geometry
from .HybridPICModel import hybridpicmodel
from .Interpolation import interpolation
from .Lasers import lasers
from .PSATD import psatd
from .Particles import newspecies, particles
from .WarpX import warpx
from ._libwarpx import found_cxx, libwarpx

# This is a circular import and must happen after the import of libwarpx
from . import picmi # isort:skip

# TODO
#__version__ = cxx.__version__
#__doc__ = cxx.__doc__
#__license__ = cxx.__license__
#__author__ = cxx.__author__
