add_warpx_test(
    background_mcc # name
    2  # dims
    ON # mpi
    2  # mpi_np
    Examples/Physics_applications/capacitive_discharge/inputs_2d_base # inputs_base
    Examples/analysis_default_regression.py # analysis
    diags/diag1000050 # output
    Examples/Physics_applications/capacitive_discharge/inputs_2d_test_background_mcc # inputs_test
)
