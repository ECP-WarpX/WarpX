#! /usr/bin/env python

import yt

# This test shoots a beam of electrons at cubic embedded boundary geometry
# At time step 40, none of the particles have hit the boundary yet. At time
# step 60, all of them should have been absorbed by the boundary. In the
# absence of the cube, none of the particles would have had time to exit
# the problem domain yet.

# all particles are still there
ds40 = yt.load("particle_absorption_plt00040")
np40 = ds40.index.particle_headers['electrons'].num_particles
assert(np40 == 612)

# all particles have been removed
ds60 = yt.load("particle_absorption_plt00060")
np60 = ds60.index.particle_headers['electrons'].num_particles
assert(np60 == 0)
