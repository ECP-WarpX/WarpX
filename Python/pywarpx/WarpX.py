# Copyright 2016-2020 Andrew Myers, David Grote, Maxence Thevenet
# Remi Lehe
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

from . import Particles
from .Algo import algo
from .Amr import amr
from .Boundary import boundary
from .Bucket import Bucket
from .Collisions import collisions, collisions_list
from .Constants import my_constants
from .Diagnostics import diagnostics
from .Geometry import geometry
from .Interpolation import interpolation
from .Langmuirwave import langmuirwave
from .Lasers import lasers, lasers_list
from .PSATD import psatd
from .Particles import particles, particles_list
from ._libwarpx import libwarpx


class WarpX(Bucket):
    """
    A Python wrapper for the WarpX C++ class
    """

    def create_argv_list(self):
        argv = []
        argv += warpx.attrlist()
        argv += my_constants.attrlist()
        argv += amr.attrlist()
        argv += geometry.attrlist()
        argv += boundary.attrlist()
        argv += algo.attrlist()
        argv += langmuirwave.attrlist()
        argv += interpolation.attrlist()
        argv += psatd.attrlist()

        # --- Search through species_names and add any predefined particle objects in the list.
        particles_list_names = [p.instancename for p in particles_list]
        for pstring in particles.species_names:
            if pstring in particles_list_names:
                # --- The species is already included in particles_list
                continue
            elif hasattr(Particles, pstring):
                # --- Add the predefined species to particles_list
                particles_list.append(getattr(Particles, pstring))
                particles_list_names.append(pstring)
            else:
                raise Exception('Species %s listed in species_names not defined'%pstring)

        argv += particles.attrlist()
        for particle in particles_list:
            argv += particle.attrlist()

        argv += collisions.attrlist()
        for collision in collisions_list:
            argv += collision.attrlist()

        argv += lasers.attrlist()
        for laser in lasers_list:
            argv += laser.attrlist()

        diagnostics.diags_names = diagnostics._diagnostics_dict.keys()
        argv += diagnostics.attrlist()
        for diagnostic in diagnostics._diagnostics_dict.values():
            diagnostic.species = diagnostic._species_dict.keys()
            argv += diagnostic.attrlist()
            for species_diagnostic in diagnostic._species_dict.values():
                argv += species_diagnostic.attrlist()

        return argv

    def init(self, mpi_comm=None):
        argv = ['warpx'] + self.create_argv_list()
        libwarpx.initialize(argv, mpi_comm=mpi_comm)

    def evolve(self, nsteps=-1):
        libwarpx.evolve(nsteps)

    def finalize(self, finalize_mpi=1):
        libwarpx.finalize(finalize_mpi)

    def getProbLo(self, direction):
        return libwarpx.libwarpx_so.warpx_getProbLo(direction)

    def getProbHi(self, direction):
        return libwarpx.libwarpx_so.warpx_getProbHi(direction)

    def write_inputs(self, filename='inputs', **kw):
        argv = self.create_argv_list()

        for k, v in kw.items():
            argv.append(f'{k} = {v}')

        # Sort the argv list to make it more human readable
        argv.sort()

        with open(filename, 'w') as ff:

            for arg in argv:
                ff.write(f'{arg}\n')

warpx = WarpX('warpx')
