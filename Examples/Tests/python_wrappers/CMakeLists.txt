# Add tests (alphabetical order) ##############################################
#

# python_wrappers_picmi #######################################################
#
if(WarpX_FFT)
    add_warpx_test(
        python_wrappers_picmi # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/python_wrappers/inputs_2d_test_python_wrappers_picmi.py # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000100 # output
    )
endif()
