# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_laser_acceleration # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_1d_laser_acceleration # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000100 # output
    OFF # dependency
)

add_warpx_test(
    test_1d_laser_acceleration_fluid # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_1d_laser_acceleration_fluid # inputs
    Examples/Physics_applications/laser_acceleration/analysis_1d_fluid.py # analysis
    diags/diag1040000 # output
    OFF # dependency
)

add_warpx_test(
    test_1d_laser_acceleration_fluid_boosted # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_1d_laser_acceleration_fluid_boosted # inputs
    Examples/Physics_applications/laser_acceleration/analysis_1d_fluid_boosted.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)

add_warpx_test(
    test_1d_laser_acceleration_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_1d_laser_acceleration_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000100 # output
    OFF # dependency
)

add_warpx_test(
    test_2d_laser_acceleration_boosted # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_2d_laser_acceleration_boosted # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000002 # output
    OFF # dependency
)

add_warpx_test(
    test_2d_laser_acceleration_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_2d_laser_acceleration_mr # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000200 # output
    OFF # dependency
)

add_warpx_test(
    test_2d_laser_acceleration_mr_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_2d_laser_acceleration_mr_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000200 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_laser_acceleration # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_3d_laser_acceleration # inputs
    Examples/analysis_default_openpmd_regression.py # analysis
    diags/diag1/ # output
    OFF # dependency
)

add_warpx_test(
    test_3d_laser_acceleration_picmi # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_3d_laser_acceleration_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000100 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_laser_acceleration_single_precision_comms # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_3d_laser_acceleration_single_precision_comms # inputs
    Examples/analysis_default_openpmd_regression.py # analysis
    diags/diag1/ # output
    OFF # dependency
)

add_warpx_test(
    test_rz_laser_acceleration # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_rz_laser_acceleration # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
    OFF # dependency
)

add_warpx_test(
    test_rz_laser_acceleration_opmd # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_rz_laser_acceleration_opmd # inputs
    Examples/Physics_applications/laser_acceleration/analysis_openpmd_rz.py # analysis
    diags/diag1/ # output
    OFF # dependency
)

add_warpx_test(
    test_rz_laser_acceleration_picmi # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_rz_laser_acceleration_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
    OFF # dependency
)

add_warpx_test(
    test_2d_refined_injection # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/inputs_test_2d_refined_injection # inputs
    Examples/Physics_applications/laser_acceleration/analysis_refined_injection.py # analysis
    diags/diag1000200 # output
    OFF # dependency
)
