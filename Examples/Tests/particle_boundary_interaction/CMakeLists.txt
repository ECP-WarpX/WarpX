# Add tests (alphabetical order) ##############################################
#

# particle_boundary_interaction ###############################################
#
add_warpx_test(
    particle_boundary_interaction # name
    RZ # dims
    2  # nprocs
    Examples/Tests/particle_boundary_interaction/inputs_rz_test_particle_boundary_interaction_picmi.py # inputs
    Examples/Tests/particle_boundary_interaction/analysis.py # analysis
    diags/diag1000023 # output
)
