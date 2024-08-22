# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
    add_warpx_test(
        test_3d_point_of_contact_eb # name
        3 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/point_of_contact_eb/inputs_test_3d_point_of_contact_eb # inputs
        Examples/Tests/point_of_contact_eb/analysis.py # analysis
        diags/diag1/ # output
    )
endif()

if(WarpX_EB)
    add_warpx_test(
        test_rz_point_of_contact_eb # name
        RZ # dims
        2  # nprocs
        ON # eb
        Examples/Tests/point_of_contact_eb/inputs_test_rz_point_of_contact_eb # inputs
        Examples/Tests/point_of_contact_eb/analysis.py # analysis
        diags/diag1/ # output
    )
endif()
