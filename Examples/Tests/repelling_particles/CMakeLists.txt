# Add tests (alphabetical order) ##############################################
#

# repelling_particles #########################################################
#
add_warpx_test(
    repelling_particles # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/repelling_particles/inputs_2d_test_repelling_particles # inputs
    Examples/Tests/repelling_particles/analysis.py # analysis
    diags/diag1000200 # output
)
