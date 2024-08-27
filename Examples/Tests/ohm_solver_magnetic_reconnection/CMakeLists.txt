# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_ohm_solver_magnetic_reconnection_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    "Examples/Tests/ohm_solver_magnetic_reconnection/inputs_test_2d_ohm_solver_magnetic_reconnection_picmi.py --test" # inputs
    Examples/Tests/ohm_solver_magnetic_reconnection/analysis.py # analysis
    diags/diag1000020 # output
    OFF # dependency
)
