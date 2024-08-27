# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_laser_injection_from_binary_file_prepare # name
    2 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_2d_laser_injection_from_binary_file_prepare.py # inputs
    OFF # analysis
    OFF # output
    OFF # dependency
)

add_warpx_test(
    test_2d_laser_injection_from_binary_file # name
    2 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_2d_laser_injection_from_binary_file # inputs
    Examples/Tests/laser_injection_from_file/analysis_2d_binary.py # analysis
    diags/diag1000250 # output
    test_2d_laser_injection_from_binary_file_prepare # dependency
)

add_warpx_test(
    test_1d_laser_injection_from_lasy_file_prepare # name
    1 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_1d_laser_injection_from_lasy_file_prepare.py # inputs
    OFF # analysis
    OFF # output
    OFF # dependency
)

add_warpx_test(
    test_1d_laser_injection_from_lasy_file # name
    1 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_1d_laser_injection_from_lasy_file # inputs
    Examples/Tests/laser_injection_from_file/analysis_1d.py # analysis
    diags/diag1000251 # output
    test_1d_laser_injection_from_lasy_file_prepare # dependency
)

add_warpx_test(
    test_1d_laser_injection_from_lasy_file_boost_prepare # name
    1 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_1d_laser_injection_from_lasy_file_boost_prepare.py # inputs
    OFF # analysis
    OFF # output
    OFF # dependency
)

add_warpx_test(
    test_1d_laser_injection_from_lasy_file_boost # name
    1 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_1d_laser_injection_from_lasy_file_boost # inputs
    Examples/Tests/laser_injection_from_file/analysis_1d_boost.py # analysis
    diags/diag1000001 # output
    test_1d_laser_injection_from_lasy_file_boost_prepare # dependency
)

add_warpx_test(
    test_2d_laser_injection_from_lasy_file_prepare # name
    2 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_2d_laser_injection_from_lasy_file_prepare.py # inputs
    OFF # analysis
    OFF # output
    OFF # dependency
)

add_warpx_test(
    test_2d_laser_injection_from_lasy_file # name
    2 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_2d_laser_injection_from_lasy_file # inputs
    Examples/Tests/laser_injection_from_file/analysis_2d.py # analysis
    diags/diag1000251 # output
    test_2d_laser_injection_from_lasy_file_prepare # dependency
)

add_warpx_test(
    test_3d_laser_injection_from_lasy_file_prepare # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_3d_laser_injection_from_lasy_file_prepare.py # inputs
    OFF # analysis
    OFF # output
    OFF # dependency
)

add_warpx_test(
    test_3d_laser_injection_from_lasy_file # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_3d_laser_injection_from_lasy_file # inputs
    Examples/Tests/laser_injection_from_file/analysis_3d.py # analysis
    diags/diag1000251 # output
    test_3d_laser_injection_from_lasy_file_prepare # dependency
)

add_warpx_test(
    test_rz_laser_injection_from_lasy_file_prepare # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_rz_laser_injection_from_lasy_file_prepare.py # inputs
    OFF # analysis
    OFF # output
    OFF # dependency
)

add_warpx_test(
    test_rz_laser_injection_from_lasy_file # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_rz_laser_injection_from_lasy_file # inputs
    Examples/Tests/laser_injection_from_file/analysis_rz.py # analysis
    diags/diag1000252 # output
    test_rz_laser_injection_from_lasy_file_prepare # dependency
)

add_warpx_test(
    test_rz_laser_injection_from_RZ_lasy_file_prepare # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_rz_laser_injection_from_RZ_lasy_file_prepare.py # inputs
    OFF # analysis
    OFF # output
    OFF # dependency
)

add_warpx_test(
    test_rz_laser_injection_from_RZ_lasy_file # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/laser_injection_from_file/inputs_test_rz_laser_injection_from_RZ_lasy_file # inputs
    Examples/Tests/laser_injection_from_file/analysis_from_RZ_file.py # analysis
    diags/diag1000612 # output
    test_rz_laser_injection_from_RZ_lasy_file_prepare # dependency
)
