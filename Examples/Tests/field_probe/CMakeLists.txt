# Add tests (alphabetical order) ##############################################
#

# field_probe #################################################################
#
if(WarpX_EB)
    add_warpx_test(
        field_probe # name
        2 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/field_probe/inputs_2d_test_field_probe # inputs
        Examples/Tests/field_probe/analysis_field_probe.py # analysis
        diags/diag1000544 # output
    )
endif()
