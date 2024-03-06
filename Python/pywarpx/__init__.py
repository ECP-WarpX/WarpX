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
        p_abs = os.path.abspath(os.path.expanduser(os.path.expandvars(p)))
        if os.path.exists(p_abs):
            os.add_dll_directory(p_abs)

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
from .LoadThirdParty import load_cupy
from .PSATD import psatd
from .Particles import newspecies, particles
from .WarpX import warpx
from ._libwarpx import libwarpx

# This is a circular import and must happen after the import of libwarpx
from . import picmi  # isort:skip

# intentionally query the value - only set once sim dimension is known
def __getattr__(name):
    # https://stackoverflow.com/a/57263518/2719194
    if name == '__version__':
        return libwarpx.__version__
    raise AttributeError(f"module '{__name__}' has no attribute '{name}'")

# TODO
#__doc__ = cxx.__doc__
#__license__ = cxx.__license__
#__author__ = cxx.__author__
