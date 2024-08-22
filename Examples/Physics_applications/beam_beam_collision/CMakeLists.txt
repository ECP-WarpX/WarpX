# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_beam_beam_collision # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/beam_beam_collision/inputs_test_3d_beam_beam_collision # inputs
    Examples/analysis_default_openpmd_regression.py # analysis
    diags/diag1/ # output
)
