# Add tests (alphabetical order) ##############################################
#

# laser_acceleration_1d #######################################################
#
add_warpx_test(
    laser_acceleration_1d # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_1d_test_laser_acceleration_1d # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000100 # output
)

# laser_acceleration_1d_fluid #################################################
#
add_warpx_test(
    laser_acceleration_1d_fluid # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_1d_test_laser_acceleration_1d_fluid # inputs
    Examples/Physics_applications/laser_acceleration/analysis_1d_fluid.py # analysis
    diags/diag1040000 # output
)

# laser_acceleration_1d_fluid_boosted #########################################
#
add_warpx_test(
    laser_acceleration_1d_fluid_boosted # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_1d_test_laser_acceleration_1d_fluid_boosted # inputs
    Examples/Physics_applications/laser_acceleration/analysis_1d_fluid_boosted.py # analysis
    diags/diag1000001 # output
)

# laser_acceleration_1d_picmi #################################################
#
add_warpx_test(
    laser_acceleration_1d_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_1d_test_laser_acceleration_1d_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000100 # output
)

# laser_acceleration_2d_boosted ###############################################
#
add_warpx_test(
    laser_acceleration_2d_boosted # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_2d_test_laser_acceleration_2d_boosted # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000002 # output
)

# laser_acceleration_2d_mr ####################################################
#
add_warpx_test(
    laser_acceleration_2d_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_2d_test_laser_acceleration_2d_mr # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000200 # output
)

# laser_acceleration_3d #######################################################
#
add_warpx_test(
    laser_acceleration_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_3d_test_laser_acceleration_3d # inputs
    Examples/analysis_default_openpmd_regression.py # analysis
    diags/diag1/ # output
)

# laser_acceleration_3d_single_precision_comms ################################
#
add_warpx_test(
    laser_acceleration_3d_single_precision_comms # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_3d_test_laser_acceleration_3d_single_precision_comms # inputs
    Examples/analysis_default_openpmd_regression.py # analysis
    diags/diag1/ # output
)

# laser_acceleration_rz #######################################################
#
add_warpx_test(
    laser_acceleration_rz # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_rz_test_laser_acceleration_rz # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

# laser_acceleration_rz_opmd ##################################################
#
add_warpx_test(
    laser_acceleration_rz_opmd # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_rz_test_laser_acceleration_rz_opmd # inputs
    Examples/Physics_applications/laser_acceleration/analysis_openpmd_rz.py # analysis
    diags/diag1/ # output
)

# refined_injection ###########################################################
#
add_warpx_test(
    refined_injection # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_2d_test_refined_injection # inputs
    Examples/Physics_applications/laser_acceleration/analysis_refined_injection.py # analysis
    diags/diag1000200 # output
)
