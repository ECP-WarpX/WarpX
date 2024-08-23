# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_background_mcc_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    "Examples/Physics_applications/capacitive_discharge/inputs_base_1d_picmi.py --test --pythonsolver" # inputs
    Examples/Physics_applications/capacitive_discharge/analysis_1d.py # analysis
    diags/diag1000050 # output
)

add_warpx_test(
    test_2d_background_mcc # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/capacitive_discharge/inputs_test_2d_background_mcc # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000050 # output
)

# FIXME: can we make this a single precision for now?
#add_warpx_test(
#    test_2d_background_mcc_dp_psp # name
#    2 # dims
#    2 # nprocs
#    OFF # eb
#    OFF # preparation script
#    Examples/Physics_applications/capacitive_discharge/inputs_test_2d_background_mcc_dp_psp # inputs
#    Examples/analysis_default_regression.py # analysis
#    diags/diag1000050 # output
#)

add_warpx_test(
    test_2d_background_mcc_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Physics_applications/capacitive_discharge/inputs_test_2d_background_mcc_picmi.py # inputs
    Examples/Physics_applications/capacitive_discharge/analysis_2d.py # analysis
    diags/diag1000050 # output
)

add_warpx_test(
    test_1d_dsmc_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    "Examples/Physics_applications/capacitive_discharge/inputs_base_1d_picmi.py --test -dsmc" # inputs
    Examples/Physics_applications/capacitive_discharge/analysis_dsmc.py # analysis
    diags/diag1000050 # output
)
