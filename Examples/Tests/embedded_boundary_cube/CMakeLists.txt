# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
    add_warpx_test(
        test_2d_embedded_boundary_cube # name
        2 # dims
        1 # nprocs
        ON # eb
        Examples/Tests/embedded_boundary_cube/inputs_test_2d_embedded_boundary_cube # inputs
        Examples/Tests/embedded_boundary_cube/analysis_fields_2d.py # analysis
        diags/diag1000114 # output
    )
endif()

if(WarpX_EB)
    add_warpx_test(
        test_3d_embedded_boundary_cube # name
        3 # dims
        1 # nprocs
        ON # eb
        Examples/Tests/embedded_boundary_cube/inputs_test_3d_embedded_boundary_cube # inputs
        Examples/Tests/embedded_boundary_cube/analysis_fields.py # analysis
        diags/diag1000208 # output
    )
endif()

if(WarpX_EB)
    add_warpx_test(
        test_3d_embedded_boundary_cube_macroscopic # name
        3 # dims
        1 # nprocs
        ON # eb
        Examples/Tests/embedded_boundary_cube/inputs_test_3d_embedded_boundary_cube_macroscopic # inputs
        Examples/Tests/embedded_boundary_cube/analysis_fields.py # analysis
        diags/diag1000208 # output
    )
endif()
