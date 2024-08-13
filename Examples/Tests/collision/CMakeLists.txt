# Add tests (alphabetical order) ##############################################
#

# collision_z #################################################################
#
add_warpx_test(
    collision_z # name
    1 # dims
    2 # nprocs
    Examples/Tests/collision/inputs_1d_test_collision_z # inputs
    Examples/Tests/collision/analysis_collision_1d.py # analysis
    diags/diag1000600 # output
)

# collision_iso ###############################################################
#
add_warpx_test(
    collision_iso # name
    3 # dims
    1 # nprocs
    Examples/Tests/collision/inputs_3d_test_collision_iso # inputs
    Examples/Tests/collision/analysis_collision_3d_isotropization.py # analysis
    diags/diag1000100 # output
)

# collision_rz ################################################################
#
add_warpx_test(
    collision_rz # name
    RZ # dims
    1  # nprocs
    Examples/Tests/collision/inputs_rz_test_collision_rz # inputs
    Examples/Tests/collision/analysis_collision_rz.py # analysis
    diags/diag1000150 # output
)

# collision_xyz ###############################################################
#
add_warpx_test(
    collision_xyz # name
    3 # dims
    1 # nprocs
    Examples/Tests/collision/inputs_3d_test_collision_xyz # inputs
    Examples/Tests/collision/analysis_collision_3d.py # analysis
    collision_xyz_plt000150 # output
)

# collision_xz ################################################################
#
add_warpx_test(
    collision_xz # name
    2 # dims
    1 # nprocs
    Examples/Tests/collision/inputs_2d_test_collision_xz # inputs
    Examples/Tests/collision/analysis_collision_2d.py # analysis
    collision_xz_plt000150 # output
)

# collision_xz_picmi ##########################################################
#
add_warpx_test(
    collision_xz_picmi # name
    2 # dims
    1 # nprocs
    Examples/Tests/collision/inputs_2d_test_collision_xz_picmi.py # inputs
    Examples/Tests/collision/analysis_collision_2d.py # analysis
    collision_xz_plt000150 # output
)
