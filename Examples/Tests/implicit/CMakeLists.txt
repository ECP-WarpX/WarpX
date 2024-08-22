# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_semi_implicit_picard # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/implicit/inputs_test_1d_semi_implicit_picard # inputs
    Examples/Tests/implicit/analysis_1d.py # analysis
    diags/diag1000100 # output
)

add_warpx_test(
    test_1d_theta_implicit_picard # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/implicit/inputs_test_1d_theta_implicit_picard # inputs
    Examples/Tests/implicit/analysis_1d.py # analysis
    diags/diag1000100 # output
)

add_warpx_test(
    test_2d_theta_implicit_jfnk_vandb # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/implicit/inputs_test_2d_theta_implicit_jfnk_vandb # inputs
    Examples/Tests/implicit/analysis_vandb_jfnk_2d.py # analysis
    diags/diag1000020 # output
)

add_warpx_test(
    test_2d_theta_implicit_jfnk_vandb_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/implicit/inputs_test_2d_theta_implicit_jfnk_vandb_picmi.py # inputs
    Examples/Tests/implicit/analysis_vandb_jfnk_2d.py # analysis
    diags/diag1000020 # output
)
