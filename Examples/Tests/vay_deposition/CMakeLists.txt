# Add tests (alphabetical order) ##############################################
#

# vay_deposition_2d ###########################################################
#
if(WarpX_FFT)
    add_warpx_test(
        vay_deposition_2d # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/vay_deposition/inputs_2d_test_vay_deposition_2d # inputs
        Examples/Tests/vay_deposition/analysis.py # analysis
        diags/diag1000050 # output
    )
endif()

# vay_deposition_3d ###########################################################
#
if(WarpX_FFT)
    add_warpx_test(
        vay_deposition_3d # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/vay_deposition/inputs_3d_test_vay_deposition_3d # inputs
        Examples/Tests/vay_deposition/analysis.py # analysis
        diags/diag1000025 # output
    )
endif()
