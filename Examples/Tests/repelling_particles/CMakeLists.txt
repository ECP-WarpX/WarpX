# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_repelling_particles # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/repelling_particles/inputs_test_2d_repelling_particles # inputs
    Examples/Tests/repelling_particles/analysis.py # analysis
    diags/diag1000200 # output
)
