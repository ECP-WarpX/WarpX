# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_collision_z # name
    1 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/collision/inputs_test_1d_collision_z # inputs
    Examples/Tests/collision/analysis_collision_1d.py # analysis
    diags/diag1000600 # output
)

add_warpx_test(
    test_3d_collision_iso # name
    3 # dims
    1 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/collision/inputs_test_3d_collision_iso # inputs
    Examples/Tests/collision/analysis_collision_3d_isotropization.py # analysis
    diags/diag1000100 # output
)

add_warpx_test(
    test_rz_collision # name
    RZ # dims
    1  # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/collision/inputs_test_rz_collision # inputs
    Examples/Tests/collision/analysis_collision_rz.py # analysis
    diags/diag1000150 # output
)

add_warpx_test(
    test_3d_collision_xyz # name
    3 # dims
    1 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/collision/inputs_test_3d_collision_xyz # inputs
    Examples/Tests/collision/analysis_collision_3d.py # analysis
    diags/diag1000150 # output
)

add_warpx_test(
    test_2d_collision_xz # name
    2 # dims
    1 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/collision/inputs_test_2d_collision_xz # inputs
    Examples/Tests/collision/analysis_collision_2d.py # analysis
    diags/diag1000150 # output
)

add_warpx_test(
    test_2d_collision_xz_picmi # name
    2 # dims
    1 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/collision/inputs_test_2d_collision_xz_picmi.py # inputs
    Examples/Tests/collision/analysis_collision_2d.py # analysis
    diags/diag1000150 # output
)
