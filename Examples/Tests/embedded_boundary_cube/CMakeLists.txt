# Add tests (alphabetical order) ##############################################
#

# embedded_boundary_cube_2d ###################################################
#
add_warpx_test(
    embedded_boundary_cube_2d # name
    2 # dims
    1 # nprocs
    ON # eb
    Examples/Tests/embedded_boundary_cube/inputs_2d_test_embedded_boundary_cube_2d # inputs
    Examples/Tests/embedded_boundary_cube/analysis_fields_2d.py # analysis
    diags/diag1000114 # output
)

# embedded_boundary_cube_3d ###################################################
#
add_warpx_test(
    embedded_boundary_cube_3d # name
    3 # dims
    1 # nprocs
    ON # eb
    Examples/Tests/embedded_boundary_cube/inputs_3d_test_embedded_boundary_cube_3d # inputs
    Examples/Tests/embedded_boundary_cube/analysis_fields.py # analysis
    diags/diag1000208 # output
)

# embedded_boundary_cube_3d_macroscopic #######################################
#
add_warpx_test(
    embedded_boundary_cube_3d_macroscopic # name
    3 # dims
    1 # nprocs
    ON # eb
    Examples/Tests/embedded_boundary_cube/inputs_3d_test_embedded_boundary_cube_3d_macroscopic # inputs
    Examples/Tests/embedded_boundary_cube/analysis_fields.py # analysis
    diags/diag1000208 # output
)
