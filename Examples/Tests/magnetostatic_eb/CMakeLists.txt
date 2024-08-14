# Add tests (alphabetical order) ##############################################
#

# magnetostatic_eb_3d #########################################################
#
add_warpx_test(
    magnetostatic_eb_3d # name
    3 # dims
    1 # nprocs
    ON # eb
    Examples/Tests/magnetostatic_eb/inputs_3d_test_magnetostatic_eb_3d # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000001 # output
)

# magnetostatic_eb_3d_picmi ###################################################
#
add_warpx_test(
    magnetostatic_eb_3d_picmi # name
    3 # dims
    1 # nprocs
    ON # eb
    Examples/Tests/magnetostatic_eb/inputs_3d_test_magnetostatic_eb_3d_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000001 # output
)

# magnetostatic_eb_rz_picmi ###################################################
#
add_warpx_test(
    magnetostatic_eb_rz_picmi # name
    RZ # dims
    1  # nprocs
    ON # eb
    Examples/Tests/magnetostatic_eb/inputs_rz_test_magnetostatic_eb_rz_picmi.py # inputs
    Examples/Tests/magnetostatic_eb/analysis_rz.py # analysis
    diags/diag1000001 # output
)
