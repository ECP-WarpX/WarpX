# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
    add_warpx_test(
        test_2d_embedded_circle # name
        2 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/embedded_circle/inputs_test_2d_embedded_circle # inputs
        Examples/Tests/embedded_circle/analysis.py # analysis
        diags/diag1000011
    )
endif()
