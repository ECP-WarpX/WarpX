#!/usr/bin/env python3

# 2023 TAE Technologies

import os
import sys

import numpy as np

sys.path.append('../../../../warpx/Regression/Checksum/')

import checksumAPI

# this will be the name of the plot file
fn = sys.argv[1]
test_name = os.path.split(os.getcwd())[1]

my_check = checksumAPI.evaluate_checksum(test_name, fn, do_particles=True)

ref_density = np.array([
       1.27953969e+14, 2.23553999e+14, 2.55384510e+14, 2.55663110e+14,
       2.55805760e+14, 2.55812087e+14, 2.55813911e+14, 2.55754104e+14,
       2.55929601e+14, 2.56085472e+14, 2.55932867e+14, 2.55828121e+14,
       2.55901711e+14, 2.55985074e+14, 2.56182697e+14, 2.56446847e+14,
       2.56483696e+14, 2.56301187e+14, 2.56245301e+14, 2.56797584e+14,
       2.57257907e+14, 2.57023627e+14, 2.56500876e+14, 2.56106851e+14,
       2.56283546e+14, 2.56723967e+14, 2.56960855e+14, 2.56825486e+14,
       2.56674669e+14, 2.56567191e+14, 2.56310927e+14, 2.56361171e+14,
       2.56692197e+14, 2.56743606e+14, 2.56653108e+14, 2.56883854e+14,
       2.56763228e+14, 2.56343726e+14, 2.56385489e+14, 2.56570110e+14,
       2.56538112e+14, 2.56472179e+14, 2.56322922e+14, 2.56195384e+14,
       2.56474576e+14, 2.56764233e+14, 2.56533016e+14, 2.56257170e+14,
       2.56362463e+14, 2.56363962e+14, 2.56311292e+14, 2.56678788e+14,
       2.57061138e+14, 2.56785892e+14, 2.56406603e+14, 2.56334908e+14,
       2.56120051e+14, 2.56003269e+14, 2.56132187e+14, 2.56329572e+14,
       2.56535713e+14, 2.56708950e+14, 2.56661860e+14, 2.56448986e+14,
       2.56386823e+14, 2.56233660e+14, 2.56137632e+14, 2.56206263e+14,
       2.56364996e+14, 2.56483536e+14, 2.56308741e+14, 2.56447231e+14,
       2.56896301e+14, 2.56691405e+14, 2.56170780e+14, 2.56122216e+14,
       2.56427399e+14, 2.56897558e+14, 2.56928868e+14, 2.56659033e+14,
       2.56749993e+14, 2.56952497e+14, 2.56798907e+14, 2.56377081e+14,
       2.56453057e+14, 2.56796632e+14, 2.56944576e+14, 2.57248469e+14,
       2.57279426e+14, 2.56849516e+14, 2.56601834e+14, 2.56850545e+14,
       2.56953072e+14, 2.56442586e+14, 2.56329006e+14, 2.56790661e+14,
       2.57083582e+14, 2.57075550e+14, 2.56719615e+14, 2.56220486e+14,
       2.56222323e+14, 2.56547365e+14, 2.56499423e+14, 2.56434041e+14,
       2.56378587e+14, 2.56249892e+14, 2.56380492e+14, 2.56504513e+14,
       2.56337631e+14, 2.56204891e+14, 2.56325116e+14, 2.56297798e+14,
       2.56112782e+14, 2.56054218e+14, 2.56320120e+14, 2.56580938e+14,
       2.56446800e+14, 2.56267011e+14, 2.56372853e+14, 2.56617592e+14,
       2.56630745e+14, 2.56615242e+14, 2.56625259e+14, 2.56561320e+14,
       2.56640072e+14, 2.56693273e+14, 2.56613237e+14, 2.24169847e+14,
       1.27683197e+14
])

density_data = np.load( 'ion_density_case_1.npy' )
print(repr(density_data))
assert np.allclose(density_data, ref_density, rtol=1e-3)
