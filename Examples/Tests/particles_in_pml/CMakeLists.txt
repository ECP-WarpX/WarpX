# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_particles_in_pml # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/particles_in_pml/inputs_test_2d_particles_in_pml # inputs
    Examples/Tests/particles_in_pml/analysis_particles_in_pml.py # analysis
    diags/diag1000180 # output
)

add_warpx_test(
    test_2d_particles_in_pml_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/particles_in_pml/inputs_test_2d_particles_in_pml_mr # inputs
    Examples/Tests/particles_in_pml/analysis_particles_in_pml.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_3d_particles_in_pml # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/particles_in_pml/inputs_test_3d_particles_in_pml # inputs
    Examples/Tests/particles_in_pml/analysis_particles_in_pml.py # analysis
    diags/diag1000120 # output
)

add_warpx_test(
    test_3d_particles_in_pml_mr # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/particles_in_pml/inputs_test_3d_particles_in_pml_mr # inputs
    Examples/Tests/particles_in_pml/analysis_particles_in_pml.py # analysis
    diags/diag1000200 # output
)
