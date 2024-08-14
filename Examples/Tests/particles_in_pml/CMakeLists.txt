# Add tests (alphabetical order) ##############################################
#

# particles_in_pml_2d #########################################################
#
add_warpx_test(
    particles_in_pml_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/particles_in_pml/inputs_2d_test_particles_in_pml_2d # inputs
    Examples/Tests/particles_in_pml/analysis_particles_in_pml.py # analysis
    diags/diag1000180 # output
)

# particles_in_pml_2d_mr ######################################################
#
add_warpx_test(
    particles_in_pml_2d_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/particles_in_pml/inputs_2d_test_particles_in_pml_2d_mr # inputs
    Examples/Tests/particles_in_pml/analysis_particles_in_pml.py # analysis
    diags/diag1000300 # output
)

# particles_in_pml_3d #########################################################
#
add_warpx_test(
    particles_in_pml_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/particles_in_pml/inputs_3d_test_particles_in_pml_3d # inputs
    Examples/Tests/particles_in_pml/analysis_particles_in_pml.py # analysis
    diags/diag1000120 # output
)

# particles_in_pml_3d_mr ######################################################
#
add_warpx_test(
    particles_in_pml_3d_mr # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/particles_in_pml/inputs_3d_test_particles_in_pml_3d_mr # inputs
    Examples/Tests/particles_in_pml/analysis_particles_in_pml.py # analysis
    diags/diag1000200 # output
)
