# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_id_cpu_read_picmi # name
    2 # dims
    1 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/restart/inputs_test_2d_id_cpu_read_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

# FIXME
#add_warpx_test(
#    test_2d_runtime_components_picmi # name
#    2 # dims
#    1 # nprocs
#    OFF # eb
#    OFF # preparation script
#    Examples/Tests/restart/inputs_test_2d_runtime_components_picmi.py # inputs
#    Examples/analysis_default_regression.py # analysis
#    diags/diag1000010 # output
#)

add_warpx_test(
    test_3d_acceleration # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/restart/inputs_test_3d_acceleration # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

add_warpx_test(
    test_3d_acceleration_restart # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/restart/inputs_test_3d_acceleration_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

if(WarpX_FFT)
    add_warpx_test(
        test_3d_acceleration_psatd # name
        3 # dims
        2 # nprocs
        OFF # eb
        OFF # preparation script
        Examples/Tests/restart/inputs_test_3d_acceleration_psatd # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000010 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_acceleration_psatd_restart # name
        3 # dims
        2 # nprocs
        OFF # eb
        OFF # preparation script
        Examples/Tests/restart/inputs_test_3d_acceleration_psatd_restart # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000010 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_acceleration_psatd_time_avg # name
        3 # dims
        2 # nprocs
        OFF # eb
        OFF # preparation script
        Examples/Tests/restart/inputs_test_3d_acceleration_psatd_time_avg # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000010 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_acceleration_psatd_time_avg_restart # name
        3 # dims
        2 # nprocs
        OFF # eb
        OFF # preparation script
        Examples/Tests/restart/inputs_test_3d_acceleration_psatd_time_avg_restart # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000010 # output
    )
endif()
