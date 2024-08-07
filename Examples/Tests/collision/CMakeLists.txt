# Add tests (alphabetical order) ##############################################
#

# collisionZ ##################################################################
#
add_warpx_test(
    collisionZ # name
    1  # dims
    ON # mpi
    Examples/Tests/collision/inputs_1d_test_collisionZ # inputs
    Examples/Tests/collision/analysis_collision_1d.py # analysis
    diags/diag1000600 # output
)

# collisionISO ################################################################
#
add_warpx_test(
    collisionISO # name
    3   # dims
    OFF # mpi
    Examples/Tests/collision/inputs_3d_test_collisionISO # inputs
    Examples/Tests/collision/analysis_collision_3d_isotropization.py # analysis
    diags/diag1000100 # output
)

# collisionRZ #################################################################
#
add_warpx_test(
    collisionRZ # name
    RZ  # dims
    OFF # mpi
    Examples/Tests/collision/inputs_rz_test_collisionRZ # inputs
    Examples/Tests/collision/analysis_collision_rz.py # analysis
    diags/diag1000150 # output
)

# collisionXYZ ################################################################
#
add_warpx_test(
    collisionXYZ # name
    3   # dims
    OFF # mpi
    Examples/Tests/collision/inputs_3d_test_collisionXYZ # inputs
    Examples/Tests/collision/analysis_collision_3d.py # analysis
    collisionXYZ_plt000150 # output
)

# collisionXZ #################################################################
#
add_warpx_test(
    collisionXZ # name
    2   # dims
    OFF # mpi
    Examples/Tests/collision/inputs_2d_test_collisionXZ # inputs
    Examples/Tests/collision/analysis_collision_2d.py # analysis
    collisionXZ_plt000150 # output
)
