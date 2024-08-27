# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_ohm_solver_landau_damping_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    "Examples/Tests/ohm_solver_ion_Landau_damping/inputs_test_2d_ohm_solver_landau_damping_picmi.py --test --dim 2 --temp_ratio 0.1" # inputs
    Examples/Tests/ohm_solver_ion_Landau_damping/analysis.py # analysis
    diags/diag1000100 # output
    OFF # dependency
)
