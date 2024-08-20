# Add tests (alphabetical order) ##############################################
#

# reduced_diags ###############################################################
#
add_warpx_test(
    reduced_diags # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/reduced_diags/inputs_3d_test_reduced_diags # inputs
    Examples/Tests/reduced_diags/analysis_reduced_diags.py # analysis
    diags/diag1000200 # output
)

# reduced_diags_load_balance_costs_heuristic ##################################
#
add_warpx_test(
    reduced_diags_load_balance_costs_heuristic # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/reduced_diags/inputs_3d_test_reduced_diags_load_balance_costs_heuristic # inputs
    Examples/Tests/reduced_diags/analysis_reduced_diags_load_balance_costs.py # analysis
    diags/diag1000003 # output
)

# reduced_diags_load_balance_costs_timers #####################################
#
add_warpx_test(
    reduced_diags_load_balance_costs_timers # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/reduced_diags/inputs_3d_test_reduced_diags_load_balance_costs_timers # inputs
    Examples/Tests/reduced_diags/analysis_reduced_diags_load_balance_costs.py # analysis
    diags/diag1000003 # output
)

# reduced_diags_load_balance_costs_timers_picmi ###############################
#
add_warpx_test(
    reduced_diags_load_balance_costs_timers_picmi # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/reduced_diags/inputs_3d_test_reduced_diags_load_balance_costs_timers_picmi.py # inputs
    Examples/Tests/reduced_diags/analysis_reduced_diags_load_balance_costs.py # analysis
    diags/diag1000003 # output
)

# reduced_diags_load_balance_costs_timers_psatd ###############################
#
if(WarpX_FFT)
    add_warpx_test(
        reduced_diags_load_balance_costs_timers_psatd # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/reduced_diags/inputs_3d_test_reduced_diags_load_balance_costs_timers_psatd # inputs
        Examples/Tests/reduced_diags/analysis_reduced_diags_load_balance_costs.py # analysis
        diags/diag1000003 # output
    )
endif()
