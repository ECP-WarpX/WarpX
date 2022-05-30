#!/usr/bin/env python3

# Copyright 2021 Modern Electron

import numpy as np

ref_density = np.array([
    1.29561598e+14, 2.24364831e+14, 2.55382204e+14, 2.55655510e+14,
    2.55796087e+14, 2.55814532e+14, 2.55819653e+14, 2.55749600e+14,
    2.55955007e+14, 2.56238304e+14, 2.56142653e+14, 2.55919674e+14,
    2.55914356e+14, 2.56009469e+14, 2.56194906e+14, 2.56364589e+14,
    2.56435346e+14, 2.56408809e+14, 2.56163214e+14, 2.56087527e+14,
    2.56234636e+14, 2.56400067e+14, 2.56483388e+14, 2.56624061e+14,
    2.56870567e+14, 2.56753021e+14, 2.56496358e+14, 2.56606663e+14,
    2.56936370e+14, 2.56773379e+14, 2.56358559e+14, 2.56433958e+14,
    2.56636315e+14, 2.56536939e+14, 2.56276859e+14, 2.56216543e+14,
    2.56429319e+14, 2.56561273e+14, 2.56405387e+14, 2.56274938e+14,
    2.56394238e+14, 2.56608718e+14, 2.56716157e+14, 2.56761833e+14,
    2.56606540e+14, 2.56238417e+14, 2.55944200e+14, 2.55846937e+14,
    2.55975172e+14, 2.56175314e+14, 2.56339088e+14, 2.56460787e+14,
    2.56556831e+14, 2.56423240e+14, 2.56356133e+14, 2.56447669e+14,
    2.56298131e+14, 2.56253242e+14, 2.56368052e+14, 2.56299330e+14,
    2.56203592e+14, 2.56337291e+14, 2.56410342e+14, 2.56274238e+14,
    2.56243795e+14, 2.56336237e+14, 2.56610219e+14, 2.56983613e+14,
    2.57076018e+14, 2.56744188e+14, 2.56321139e+14, 2.56314582e+14,
    2.56554890e+14, 2.56469851e+14, 2.56285766e+14, 2.56332856e+14,
    2.56348966e+14, 2.56364854e+14, 2.56359671e+14, 2.56601013e+14,
    2.57250144e+14, 2.57573442e+14, 2.57243442e+14, 2.56730921e+14,
    2.56523922e+14, 2.56383456e+14, 2.56138443e+14, 2.56176785e+14,
    2.56486745e+14, 2.56649320e+14, 2.56531626e+14, 2.56385344e+14,
    2.56295084e+14, 2.56211438e+14, 2.56261221e+14, 2.56294633e+14,
    2.56408309e+14, 2.56886655e+14, 2.57027225e+14, 2.56542718e+14,
    2.56237741e+14, 2.56253270e+14, 2.56188038e+14, 2.56226339e+14,
    2.56620619e+14, 2.57042054e+14, 2.57242664e+14, 2.57326786e+14,
    2.57060063e+14, 2.56515658e+14, 2.56162572e+14, 2.56182514e+14,
    2.56457757e+14, 2.56527241e+14, 2.56311643e+14, 2.56194572e+14,
    2.56387490e+14, 2.56728856e+14, 2.56791230e+14, 2.56499033e+14,
    2.56080399e+14, 2.56168249e+14, 2.56573015e+14, 2.56584914e+14,
    2.56583250e+14, 2.56427216e+14, 2.56337983e+14, 2.24478445e+14,
    1.28237654e+14
])

density_data = np.load( 'ion_density_case_1.npy' )
print(repr(density_data))
assert np.allclose(density_data, ref_density)
