#!/usr/bin/env python3

# Copyright 2022 Modern Electron, David Grote

import numpy as np

ref_density = np.array([
       1.29556694e+14, 2.24358818e+14, 2.55381745e+14, 2.55655005e+14,
       2.55796268e+14, 2.55819108e+14, 2.55819686e+14, 2.55751184e+14,
       2.55920806e+14, 2.56072344e+14, 2.55937266e+14, 2.55849080e+14,
       2.55918981e+14, 2.55980835e+14, 2.56054153e+14, 2.56074694e+14,
       2.56036953e+14, 2.56181153e+14, 2.56322618e+14, 2.56253541e+14,
       2.56196224e+14, 2.56353090e+14, 2.56256022e+14, 2.55928997e+14,
       2.56110988e+14, 2.56658917e+14, 2.56832584e+14, 2.56551871e+14,
       2.56491186e+14, 2.56469928e+14, 2.56418625e+14, 2.56541071e+14,
       2.56513773e+14, 2.56424507e+14, 2.56302757e+14, 2.56242392e+14,
       2.56270399e+14, 2.56178952e+14, 2.56071407e+14, 2.56141949e+14,
       2.56419808e+14, 2.56606936e+14, 2.56437774e+14, 2.56252443e+14,
       2.56309513e+14, 2.56383484e+14, 2.56265140e+14, 2.56167674e+14,
       2.56466922e+14, 2.56924871e+14, 2.56901781e+14, 2.56631494e+14,
       2.56643458e+14, 2.56523464e+14, 2.56378273e+14, 2.56571301e+14,
       2.56794308e+14, 2.56788543e+14, 2.56549712e+14, 2.56303156e+14,
       2.56210811e+14, 2.56418363e+14, 2.57314552e+14, 2.58471405e+14,
       2.58169740e+14, 2.56946418e+14, 2.56726550e+14, 2.56853119e+14,
       2.56613698e+14, 2.56509538e+14, 2.56692976e+14, 2.56705132e+14,
       2.56372135e+14, 2.56167561e+14, 2.56296953e+14, 2.56498746e+14,
       2.56523099e+14, 2.56404333e+14, 2.56227098e+14, 2.56399004e+14,
       2.56614905e+14, 2.56436650e+14, 2.56388608e+14, 2.56553683e+14,
       2.56637912e+14, 2.56407782e+14, 2.56104130e+14, 2.56082338e+14,
       2.56095272e+14, 2.56278448e+14, 2.56808134e+14, 2.57127896e+14,
       2.56858173e+14, 2.56326991e+14, 2.56296032e+14, 2.56563348e+14,
       2.56482274e+14, 2.56667483e+14, 2.57072448e+14, 2.56767529e+14,
       2.56433245e+14, 2.56586564e+14, 2.56636403e+14, 2.56765624e+14,
       2.56868122e+14, 2.56783435e+14, 2.56714527e+14, 2.56651030e+14,
       2.56528399e+14, 2.56227514e+14, 2.56163300e+14, 2.56408217e+14,
       2.56433124e+14, 2.56374737e+14, 2.56542023e+14, 2.56748800e+14,
       2.56715205e+14, 2.56298166e+14, 2.56042658e+14, 2.56292458e+14,
       2.56352283e+14, 2.56370559e+14, 2.56487462e+14, 2.56483655e+14,
       2.56741185e+14, 2.56665111e+14, 2.56523794e+14, 2.24741566e+14,
       1.28486948e+14
])

density_data = np.load( 'ion_density_case_1.npy' )
print(repr(density_data))
assert np.allclose(density_data, ref_density)
