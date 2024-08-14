# Add tests (alphabetical order) ##############################################
#

# electrostatic_sphere_eb #####################################################
#
add_warpx_test(
    electrostatic_sphere_eb # name
    3 # dims
    2 # nprocs
    ON # eb
    Examples/Tests/electrostatic_sphere_eb/inputs_3d_test_electrostatic_sphere_eb # inputs
    Examples/Tests/electrostatic_sphere_eb/analysis.py # analysis
    diags/diag1000001 # output
)

# electrostatic_sphere_eb_mixed_bc ############################################
#
add_warpx_test(
    electrostatic_sphere_eb_mixed_bc # name
    3 # dims
    2 # nprocs
    ON # eb
    Examples/Tests/electrostatic_sphere_eb/inputs_3d_test_electrostatic_sphere_eb_mixed_bc # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000001 # output
)

# electrostatic_sphere_eb_rz ##################################################
#
add_warpx_test(
    electrostatic_sphere_eb_rz # name
    RZ # dims
    2  # nprocs
    ON # eb
    Examples/Tests/electrostatic_sphere_eb/inputs_rz_test_electrostatic_sphere_eb_rz # inputs
    Examples/Tests/electrostatic_sphere_eb/analysis_rz.py # analysis
    diags/diag1000001 # output
)

# electrostatic_sphere_eb_rz_mr ###############################################
#
add_warpx_test(
    electrostatic_sphere_eb_rz_mr # name
    RZ # dims
    2  # nprocs
    ON # eb
    Examples/Tests/electrostatic_sphere_eb/inputs_rz_test_electrostatic_sphere_eb_rz_mr # inputs
    Examples/Tests/electrostatic_sphere_eb/analysis_rz_mr.py # analysis
    diags/diag1/ # output
)
