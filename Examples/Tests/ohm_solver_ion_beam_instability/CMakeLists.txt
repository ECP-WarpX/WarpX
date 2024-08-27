# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_ohm_solver_ion_beam_picmi # name
    1 # dims
    2 # nprocs
    OFF # eb
    "Examples/Tests/ohm_solver_ion_beam_instability/inputs_test_1d_ohm_solver_ion_beam_picmi.py --test --dim 1 --resonant" # inputs
    Examples/Tests/ohm_solver_ion_beam_instability/analysis.py # analysis
    diags/diag1002500 # output
    OFF # dependency
)
