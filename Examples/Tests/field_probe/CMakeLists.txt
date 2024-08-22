# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
    add_warpx_test(
        test_2d_field_probe # name
        2 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/field_probe/inputs_test_2d_field_probe # inputs
        Examples/Tests/field_probe/analysis.py # analysis
        diags/diag1000544 # output
    )
endif()
