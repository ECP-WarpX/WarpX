# Add tests (alphabetical order) ##############################################
#

# id_cpu_read_picmi ###########################################################
#
add_warpx_test(
    id_cpu_read_picmi # name
    2 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/restart/inputs_2d_test_id_cpu_read_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

# FIXME
# restart_runtime_components_picmi ############################################
#
#add_warpx_test(
#    restart_runtime_components_picmi # name
#    2 # dims
#    1 # nprocs
#    OFF # eb
#    Examples/Tests/restart/inputs_2d_test_restart_runtime_components_picmi.py # inputs
#    Examples/analysis_default_regression.py # analysis
#    diags/diag1000010 # output
#)

# acceleration ################################################################
#
add_warpx_test(
    acceleration # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/restart/inputs_3d_test_acceleration # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

# acceleration_restart ########################################################
#
add_warpx_test(
    acceleration_restart # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/restart/inputs_3d_test_acceleration_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
)

# acceleration_psatd ##########################################################
#
if(WarpX_FFT)
    add_warpx_test(
        acceleration_psatd # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/restart/inputs_3d_test_acceleration_psatd # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000010 # output
    )
endif()

# acceleration_psatd_restart ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        acceleration_psatd_restart # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/restart/inputs_3d_test_acceleration_psatd_restart # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000010 # output
    )
endif()

# acceleration_psatd_time_avg #################################################
#
if(WarpX_FFT)
    add_warpx_test(
        acceleration_psatd_time_avg # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/restart/inputs_3d_test_acceleration_psatd_time_avg # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000010 # output
    )
endif()

# acceleration_psatd_time_avg_restart #########################################
#
if(WarpX_FFT)
    add_warpx_test(
        acceleration_psatd_time_avg_restart # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/restart/inputs_3d_test_acceleration_psatd_time_avg_restart # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000010 # output
    )
endif()
