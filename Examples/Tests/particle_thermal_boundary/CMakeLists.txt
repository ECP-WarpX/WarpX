# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_particle_thermal_boundary # name
    2 # dims
    2 # nprocs
    OFF # eb
    inputs_test_2d_particle_thermal_boundary # inputs
    analysis.py # analysis
    diags/diag1002000 # output
    OFF # dependency
)
