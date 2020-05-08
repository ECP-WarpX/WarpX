# First Schwinger test with "weak" EM field. No pair should be created.

import schwinger_test_analyzer

Ex_test1 = 1.e16
Ey_test1 = 0.
Ez_test1 = 0.
Bx_test1 = 16792888.570516706
By_test1 = 5256650.141557486
Bz_test1 = 18363530.799561853

schwinger_test_analyzer.do_analysis(Ex_test1, Ey_test1, Ez_test1, Bx_test1, By_test1, Bz_test1)