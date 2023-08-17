# Copyright 2016-2023 The WarpX Community
#
# This file is part of WarpX.
#
# Authors: Andrew Myers, David Grote, Lorenzo Giacomel, Axel Huebl
# License: BSD-3-Clause-LBNL

import os

# Python 3.8+ on Windows: DLL search paths for dependent
# shared libraries
# Refs.:
# - https://github.com/python/cpython/issues/80266
# - https://docs.python.org/3.8/library/os.html#os.add_dll_directory
if os.name == "nt":
    # add anything in the current directory
    pwd = __file__.rsplit(os.sep, 1)[0] + os.sep
    os.add_dll_directory(pwd)
    # add anything in PATH
    paths = os.environ.get("PATH", "")
    for p in paths.split(";"):
        if os.path.exists(p):
            os.add_dll_directory(p)

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
from ._libwarpx import libwarpx

# This is a circular import and must happen after the import of libwarpx
from . import picmi # isort:skip

# TODO
#__version__ = cxx.__version__
#__doc__ = cxx.__doc__
#__license__ = cxx.__license__
#__author__ = cxx.__author__
