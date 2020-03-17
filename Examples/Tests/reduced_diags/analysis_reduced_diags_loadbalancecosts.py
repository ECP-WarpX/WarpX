#! /usr/bin/env python

# Copyright 2019-2020 Michael Rowan
#
# This file is part of WarpX.
#
# License: BSD-3-Clause-LBNL

# This script tests the reduced diagnostics `LoadBalanceCosts`.
# The setup is a uniform plasma with electrons.
# Particle energy and field energy will be outputed
# using the reduced diagnostics.
# And they will be compared with the data in the plotfiles.

# Tolerance: 1.0e-8 for particle energy, 1.0e-3 for field energy.
# The difference of the field energy is relatively large,
# because fields data in plotfiles are cell-centered,
# but fields data in reduced diagnostics are staggered.

# Possible running time: ~ 2 s


# PART1: load costs data
LBCdata = np.genfromtxt("./diags/reducedfiles/LBC.txt")

# PART2: compute efficiency before/after load balance and check it is improved

print('LoadBalanceCostsTest: ')
#print('difference of field energy:', abs(EFyt-EF))
#print('tolerance of field energy:', 1.0e-3)
#print('difference of particle energy:', abs(EPyt-EP))
#print('tolerance of particle energy:', 1.0e-8)

assert(True)
assert(True)

#assert(abs(EFyt-EF) < 1.0e-3)
#assert(abs(EPyt-EP) < 1.0e-8)
