# Add tests (alphabetical order) ##############################################
#

# beam_beam_collision #########################################################
#
add_warpx_test(
    beam_beam_collision # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/beam-beam_collision/inputs_3d_test_beam_beam_collision # inputs
    Examples/analysis_default_openpmd_regression.py # analysis
    diags/diag1/ # output
)
