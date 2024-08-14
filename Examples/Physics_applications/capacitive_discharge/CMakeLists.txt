# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    background_mcc # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/capacitive_discharge/inputs_2d_test_background_mcc # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000050 # output
)

# FIXME
# background_mcc_dp_psp #######################################################
#
#add_warpx_test(
#    background_mcc_dp_psp # name
#    2 # dims
#    2 # nprocs
#    OFF # eb
#    Examples/Physics_applications/capacitive_discharge/inputs_2d_test_background_mcc_dp_psp # inputs
#    Examples/analysis_default_regression.py # analysis
#    diags/diag1000050 # output
#)

# background_mcc_picmi ########################################################
#
add_warpx_test(
    background_mcc_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/capacitive_discharge/inputs_2d_test_background_mcc_picmi.py # inputs
    Examples/Physics_applications/capacitive_discharge/analysis_2d.py # analysis
    diags/diag1000050 # output
)
