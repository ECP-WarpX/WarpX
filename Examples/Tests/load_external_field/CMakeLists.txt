# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_load_external_field_grid_picmi # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_test_3d_load_external_field_grid_picmi.py # inputs
    Examples/Tests/load_external_field/analysis_3d.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_3d_load_external_field_particle_picmi # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_test_3d_load_external_field_particle_picmi.py # inputs
    Examples/Tests/load_external_field/analysis_3d.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_rz_load_external_field_grid # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_test_rz_load_external_field_grid # inputs
    Examples/Tests/load_external_field/analysis_rz.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_rz_load_external_field_grid_restart # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_test_rz_load_external_field_grid_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_rz_load_external_field_particles # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_test_rz_load_external_field_particles # inputs
    Examples/Tests/load_external_field/analysis_rz.py # analysis
    diags/diag1000300 # output
)

add_warpx_test(
    test_rz_load_external_field_particles_restart # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_test_rz_load_external_field_particles_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000300 # output
)
