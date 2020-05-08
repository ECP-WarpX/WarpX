# Third Schwinger test with intermediate electric field such that average created pair per cell
# is 1. A Poisson distribution is used to obtain the weights of the particles.

import schwinger_test_analyzer

Ex_test3 = 0.
Ey_test3 = 1.0321239524474501e+17
Ez_test3 = 0.
Bx_test3 = 0.
By_test3 = 0.
Bz_test3 = 0.

schwinger_test_analyzer.do_analysis(Ex_test3, Ey_test3, Ez_test3, Bx_test3, By_test3, Bz_test3)