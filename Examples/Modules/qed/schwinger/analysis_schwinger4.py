# Fourth Schwinger test with extremely strong EM field but with E and B perpendicular and nearly
# equal so that the pair production rate is fairly low. A Gaussian distribution is used in this
# case.

import schwinger_test_analyzer

Ex_test4 = 0.
Ey_test4 = 0.
Ez_test4 = 2.5e+20
Bx_test4 = 0.
By_test4 = 833910154604.3563
Bz_test4 = 0.

schwinger_test_analyzer.do_analysis(Ex_test4, Ey_test4, Ez_test4, Bx_test4, By_test4, Bz_test4)