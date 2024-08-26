# Add tests (alphabetical order) ##############################################
#

if(WarpX_FFT)
    add_warpx_test(
        test_2d_python_wrappers_picmi # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/python_wrappers/inputs_test_2d_python_wrappers_picmi.py # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000100 # output
    OFF # dependency
    )
endif()
