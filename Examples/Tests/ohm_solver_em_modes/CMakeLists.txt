# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_ohm_solver_em_modes_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    "Examples/Tests/ohm_solver_em_modes/inputs_test_1d_ohm_solver_em_modes_picmi.py --test --dim 1 --bdir z" # inputs
    Examples/Tests/ohm_solver_em_modes/analysis.py # analysis
    diags/field_diag000250 # output
)

add_warpx_test(
    test_rz_ohm_solver_em_modes_picmi # name
    RZ # dims
    2  # nprocs
    OFF # eb
    OFF # preparation script
    "Examples/Tests/ohm_solver_em_modes/inputs_test_rz_ohm_solver_em_modes_picmi.py --test" # inputs
    Examples/Tests/ohm_solver_em_modes/analysis_rz.py # analysis
    diags/diag1000100 # output
)
