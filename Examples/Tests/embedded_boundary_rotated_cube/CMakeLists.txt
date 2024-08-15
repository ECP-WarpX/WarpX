# Add tests (alphabetical order) ##############################################
#

# embedded_boundary_rotated_cube_2d ###########################################
#
if(WarpX_EB)
    add_warpx_test(
        embedded_boundary_rotated_cube_2d # name
        2 # dims
        1 # nprocs
        ON # eb
        Examples/Tests/embedded_boundary_rotated_cube/inputs_2d_test_embedded_boundary_rotated_cube_2d # inputs
        Examples/Tests/embedded_boundary_rotated_cube/analysis_fields_2d.py # analysis
        diags/diag1000068
    )
endif()

# embedded_boundary_rotated_cube_3d ###########################################
#
if(WarpX_EB)
    add_warpx_test(
        embedded_boundary_rotated_cube_3d # name
        3 # dims
        1 # nprocs
        ON # eb
        Examples/Tests/embedded_boundary_rotated_cube/inputs_3d_test_embedded_boundary_rotated_cube_3d # inputs
        Examples/Tests/embedded_boundary_rotated_cube/analysis_fields_3d.py # analysis
        diags/diag1000111
    )
endif()
