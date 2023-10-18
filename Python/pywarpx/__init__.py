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

# `noqa` below is to prevent `pyflakes` from raising an
# `imported but not used` warning
from .Algo import algo  # noqa
from .Amr import amr  # noqa
from .Amrex import amrex  # noqa
from .Boundary import boundary  # noqa
from .Collisions import collisions  # noqa
from .Constants import my_constants  # noqa
from .Diagnostics import diagnostics, reduced_diagnostics  # noqa
from .EB2 import eb2  # noqa
from .Geometry import geometry  # noqa
from .HybridPICModel import hybridpicmodel  # noqa
from .Interpolation import interpolation  # noqa
from .Lasers import lasers  # noqa
from .LoadThirdParty import load_cupy  # noqa
from .PSATD import psatd  # noqa
from .Particles import newspecies, particles  # noqa
from .WarpX import warpx  # noqa
from ._libwarpx import libwarpx  # noqa

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
