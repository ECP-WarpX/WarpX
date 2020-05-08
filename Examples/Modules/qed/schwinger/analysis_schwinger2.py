# Second Schwinger test with stronger EM field. Many pairs are created and a Gaussian
# distribution is used to get the weights of the particles. This is the most sensitive test
# because the relative std is extremely low.

import schwinger_test_analyzer

Ex_test2 = 1.e18
Ey_test2 = 0.
Ez_test2 = 0.
Bx_test2 = 1679288857.0516706
By_test2 = 525665014.1557486
Bz_test2 = 1836353079.9561853

schwinger_test_analyzer.do_analysis(Ex_test2, Ey_test2, Ez_test2, Bx_test2, By_test2, Bz_test2)