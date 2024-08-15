# Add tests (alphabetical order) ##############################################
#

# embedded_circle #############################################################
#
if(WarpX_EB)
    add_warpx_test(
        embedded_circle # name
        2 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/embedded_circle/inputs_2d_test_embedded_circle # inputs
        Examples/Tests/embedded_circle/analysis.py # analysis
        diags/diag1000011
    )
endif()
