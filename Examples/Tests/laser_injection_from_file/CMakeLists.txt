# Add tests (alphabetical order) ##############################################
#

# laser_injection_from_binary_file ############################################
#
add_warpx_test(
    laser_injection_from_binary_file # name
    2 # dims
    1 # nprocs
    Examples/Tests/laser_injection_from_file/inputs_2d_test_laser_injection_from_binary_file # inputs
    Examples/Tests/laser_injection_from_file/analysis_2d_binary.py # analysis
    diags/diag1000250 # output
)

# laser_injection_from_lasy_file_1d ###########################################
#
add_warpx_test(
    laser_injection_from_lasy_file_1d # name
    1 # dims
    1 # nprocs
    Examples/Tests/laser_injection_from_file/inputs_1d_test_laser_injection_from_lasy_file_1d # inputs
    Examples/Tests/laser_injection_from_file/analysis_1d.py # analysis
    diags/diag1000250 # output
)

# laser_injection_from_lasy_file_1d_boost #####################################
#
add_warpx_test(
    laser_injection_from_lasy_file_1d_boost # name
    1 # dims
    1 # nprocs
    Examples/Tests/laser_injection_from_file/inputs_1d_test_laser_injection_from_lasy_file_1d_boost # inputs
    Examples/Tests/laser_injection_from_file/analysis_1d_boost.py # analysis
    diags/diag1000250 # output
)

# laser_injection_from_lasy_file_2d ###########################################
#
add_warpx_test(
    laser_injection_from_lasy_file_2d # name
    2 # dims
    1 # nprocs
    Examples/Tests/laser_injection_from_file/inputs_2d_test_laser_injection_from_lasy_file_2d # inputs
    Examples/Tests/laser_injection_from_file/analysis_2d.py # analysis
    diags/diag1000250 # output
)

# laser_injection_from_lasy_file_3d ###########################################
#
add_warpx_test(
    laser_injection_from_lasy_file_3d # name
    3 # dims
    1 # nprocs
    Examples/Tests/laser_injection_from_file/inputs_3d_test_laser_injection_from_lasy_file_3d # inputs
    Examples/Tests/laser_injection_from_file/analysis_3d.py # analysis
    diags/diag1000250 # output
)

# laser_injection_from_lasy_file_rz ###########################################
#
add_warpx_test(
    laser_injection_from_lasy_file_rz # name
    RZ # dims
    1  # nprocs
    Examples/Tests/laser_injection_from_file/inputs_rz_test_laser_injection_from_lasy_file_rz # inputs
    Examples/Tests/laser_injection_from_file/analysis_RZ.py # analysis
    diags/diag1000250 # output
)

# laser_injection_from_rz_lasy_file ###########################################
#
add_warpx_test(
    laser_injection_from_rz_lasy_file # name
    RZ # dims
    1  # nprocs
    Examples/Tests/laser_injection_from_file/inputs_rz_test_laser_injection_from_rz_lasy_file # inputs
    Examples/Tests/laser_injection_from_file/analysis_from_RZ_file.py # analysis
    diags/diag1000250 # output
)
