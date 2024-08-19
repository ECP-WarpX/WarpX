# Add tests (alphabetical order) ##############################################
#

# point_of_contact_eb_3d ######################################################
#
if(WarpX_EB)
    add_warpx_test(
        point_of_contact_eb_3d # name
        3 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/point_of_contact_eb/inputs_3d_test_point_of_contact_eb_3d # inputs
        Examples/Tests/point_of_contact_eb/analysis.py # analysis
        diags/diag1/ # output
    )
endif()

# point_of_contact_eb_rz ######################################################
#
if(WarpX_EB)
    add_warpx_test(
        point_of_contact_eb_rz # name
        RZ # dims
        2  # nprocs
        ON # eb
        Examples/Tests/point_of_contact_eb/inputs_rz_test_point_of_contact_eb_rz # inputs
        Examples/Tests/point_of_contact_eb/analysis.py # analysis
        diags/diag1/ # output
    )
endif()
