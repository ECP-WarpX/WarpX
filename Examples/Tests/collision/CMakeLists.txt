# Add tests (alphabetical order) ##############################################
#

# collisionISO ################################################################
#
add_warpx_test(
    collisionISO # name
    3   # dims
    OFF # mpi
    Examples/Tests/collision/inputs_3d_isotropization_base # inputs
    Examples/Tests/collision/analysis_collision_3d_isotropization.py # analysis
    diags/diag1000100 # output
)

# collisionRZ #################################################################
#
add_warpx_test(
    collisionRZ # name
    RZ  # dims
    OFF # mpi
    Examples/Tests/collision/inputs_rz_base # inputs
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
