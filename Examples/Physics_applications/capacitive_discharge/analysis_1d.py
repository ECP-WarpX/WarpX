#!/usr/bin/env python3

# Copyright 2022 Modern Electron, David Grote

import numpy as np

ref_density = np.array([
       1.27843016e+14, 2.23906986e+14, 2.55922987e+14, 2.55777674e+14,
       2.55641776e+14, 2.55595178e+14, 2.55699554e+14, 2.55760902e+14,
       2.55835838e+14, 2.56074685e+14, 2.56065817e+14, 2.56039872e+14,
       2.56107369e+14, 2.55981947e+14, 2.56040212e+14, 2.56275067e+14,
       2.56271145e+14, 2.56079609e+14, 2.55962779e+14, 2.55932729e+14,
       2.56251376e+14, 2.56486639e+14, 2.56313309e+14, 2.56322615e+14,
       2.56270260e+14, 2.56108530e+14, 2.56310429e+14, 2.56655195e+14,
       2.56564720e+14, 2.56314263e+14, 2.56263961e+14, 2.56270219e+14,
       2.56195273e+14, 2.56164765e+14, 2.56392754e+14, 2.56662516e+14,
       2.56857056e+14, 2.56855350e+14, 2.56789810e+14, 2.56491612e+14,
       2.56108113e+14, 2.56433846e+14, 2.56782798e+14, 2.56595759e+14,
       2.56502258e+14, 2.56563228e+14, 2.56473791e+14, 2.56199181e+14,
       2.56190646e+14, 2.56597957e+14, 2.56864480e+14, 2.56899355e+14,
       2.56752355e+14, 2.56478253e+14, 2.56314325e+14, 2.56147898e+14,
       2.56104949e+14, 2.56218409e+14, 2.56363454e+14, 2.56569672e+14,
       2.56490681e+14, 2.56417902e+14, 2.56669621e+14, 2.56836780e+14,
       2.56813501e+14, 2.56602412e+14, 2.56249585e+14, 2.56385946e+14,
       2.56771598e+14, 2.56499283e+14, 2.56223825e+14, 2.56398451e+14,
       2.56479461e+14, 2.56360871e+14, 2.56175899e+14, 2.56193241e+14,
       2.56482548e+14, 2.56520505e+14, 2.56451686e+14, 2.56663997e+14,
       2.56489183e+14, 2.56229048e+14, 2.56364939e+14, 2.56230450e+14,
       2.56163926e+14, 2.56412272e+14, 2.56385497e+14, 2.56392760e+14,
       2.56663003e+14, 2.56708819e+14, 2.56419244e+14, 2.56289577e+14,
       2.56424426e+14, 2.56665943e+14, 2.56828725e+14, 2.56510990e+14,
       2.56288146e+14, 2.56304280e+14, 2.56290620e+14, 2.56568300e+14,
       2.56684637e+14, 2.56375972e+14, 2.56093499e+14, 2.56019550e+14,
       2.56036310e+14, 2.56064238e+14, 2.56179083e+14, 2.56234436e+14,
       2.56312122e+14, 2.56731003e+14, 2.56818818e+14, 2.56297593e+14,
       2.56190175e+14, 2.56383310e+14, 2.56264712e+14, 2.56197850e+14,
       2.56224777e+14, 2.56356622e+14, 2.56646032e+14, 2.56861241e+14,
       2.56720467e+14, 2.56256466e+14, 2.56174261e+14, 2.56565621e+14,
       2.56792482e+14, 2.56868549e+14, 2.56940701e+14, 2.24201861e+14,
       1.27339197e+14
])

density_data = np.load( 'ion_density_case_1.npy' )
print(repr(density_data))
assert np.allclose(density_data, ref_density)
