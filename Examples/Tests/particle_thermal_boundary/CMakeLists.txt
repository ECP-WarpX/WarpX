# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_particle_thermal_boundary # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/particle_thermal_boundary/inputs_test_2d_particle_thermal_boundary # inputs
    Examples/Tests/particle_thermal_boundary/analysis.py # analysis
    diags/diag1002000 # output
)
