# Copyright 2017 David Grote
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# The WarpXPIC class is the beginnings of an implementation of a standard interface
# for running PIC codes. This is the run time equivalent to the PICMI standard.
# This standard would specify the API for calling the various pieces of typical
# time step loops, for example get_self_fields and put_Efields. Ideally, a user
# could write a loop using the standard and, importing one of compliant codes, be
# able to run a customized PIC simulation with that code.

from warp.run_modes.timestepper import PICAPI

from ._libwarpx import libwarpx


class WarpXPIC(PICAPI):

    def get_time(self):
        return libwarpx.libwarpx_so.warpx_gett_new(0)

    def set_time(self, time):
        for i in range(libwarpx.libwarpx_so.warpx_finestLevel()+1):
            libwarpx.libwarpx_so.warpx_sett_new(i, time)

    def get_step_size(self):
        libwarpx.libwarpx_so.warpx_ComputeDt()
        return libwarpx.libwarpx_so.warpx_getdt(0)

    def get_step_number(self):
        return libwarpx.libwarpx_so.warpx_getistep(0)

    def set_step_number(self, it):
        for i in range(libwarpx.libwarpx_so.warpx_finestLevel()+1):
            libwarpx.libwarpx_so.warpx_setistep(i, it)

    def push_positions(self, dt):
        libwarpx.libwarpx_so.warpx_PushX(0, dt)

    def push_velocities_withE(self, dt):
        libwarpx.libwarpx_so.warpx_EPushV(0, dt)

    def push_velocities_withB(self, dt):
        libwarpx.libwarpx_so.warpx_BPushV(0, dt)

    def get_self_fields(self):
        libwarpx.libwarpx_so.warpx_FieldGather(0)

    def calculate_source(self):
        libwarpx.libwarpx_so.warpx_CurrentDeposition(0)

    def push_Efields(self, dt):
        libwarpx.libwarpx_so.warpx_EvolveE(0, dt)
        libwarpx.libwarpx_so.warpx_FillBoundaryE(0, True)

    def push_Bfields(self, dt):
        libwarpx.libwarpx_so.warpx_EvolveB(0, dt)
        libwarpx.libwarpx_so.warpx_FillBoundaryB(0, True)

    def apply_particle_boundary_conditions(self):
        libwarpx.libwarpx_so.mypc_Redistribute() # Redistribute particles
        libwarpx.libwarpx_so.warpx_MoveWindow(self.istep,True) # !!! not the correct place yet
