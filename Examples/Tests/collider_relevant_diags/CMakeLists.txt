# Add tests (alphabetical order) ##############################################
#

# collider_diagnostics ########################################################
#
add_warpx_test(
    collider_diagnostics # name
    3  # dims
    ON # mpi
    2  # mpi_np
    Examples/Tests/collider_relevant_diags/inputs_3d_multiple_particles_test_collider_diagnostics # inputs
    Examples/Tests/collider_relevant_diags/analysis_multiple_particles.py # analysis
    diags/diag1000001 # output
)
