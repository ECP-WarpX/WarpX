# Add tests (alphabetical order) ##############################################
#

# semi_implicit_picard_1d #####################################################
#
add_warpx_test(
    semi_implicit_picard_1d # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/implicit/inputs_1d_test_semi_implicit_picard_1d # inputs
    Examples/Tests/implicit/analysis_1d.py # analysis
    diags/diag1000100 # output
)

# theta_implicit_jfnk_vandb_2d ################################################
#
add_warpx_test(
    theta_implicit_jfnk_vandb_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/implicit/inputs_2d_test_theta_implicit_jfnk_vandb_2d # inputs
    Examples/Tests/implicit/analysis_vandb_jfnk_2d.py # analysis
    diags/diag1000020 # output
)

# theta_implicit_jfnk_vandb_2d_picmi ##########################################
#
add_warpx_test(
    theta_implicit_jfnk_vandb_2d_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/implicit/inputs_2d_test_theta_implicit_jfnk_vandb_2d_picmi.py # inputs
    Examples/Tests/implicit/analysis_vandb_jfnk_2d.py # analysis
    diags/diag1000020 # output
)

# theta_implicit_picard_1d ####################################################
#
add_warpx_test(
    theta_implicit_picard_1d # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/implicit/inputs_1d_test_theta_implicit_picard_1d # inputs
    Examples/Tests/implicit/analysis_1d.py # analysis
    diags/diag1000100 # output
)
