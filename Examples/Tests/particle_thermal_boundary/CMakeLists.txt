# Add tests (alphabetical order) ##############################################
#

# particle_thermal_boundary ###################################################
#
add_warpx_test(
    particle_thermal_boundary # name
    2 # dims
    2 # nprocs
    Examples/Tests/particle_thermal_boundary/inputs_2d_test_particle_thermal_boundary # inputs
    Examples/Tests/particle_thermal_boundary/analysis_2d.py # analysis
    diags/diag1002000 # output
)
