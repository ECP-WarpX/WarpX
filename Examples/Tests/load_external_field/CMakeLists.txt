# Add tests (alphabetical order) ##############################################
#

# load_external_field_3d_grid_picmi ###########################################
#
add_warpx_test(
    load_external_field_3d_grid_picmi # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_3d_test_load_external_field_3d_grid_picmi.py # inputs
    Examples/Tests/load_external_field/analysis_3d.py # analysis
    diags/diag1000300 # output
)

# load_external_field_3d_particle_picmi #######################################
#
add_warpx_test(
    load_external_field_3d_particle_picmi # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_3d_test_load_external_field_3d_particle_picmi.py # inputs
    Examples/Tests/load_external_field/analysis_3d.py # analysis
    diags/diag1000300 # output
)

# load_external_field_rz_grid #################################################
#
add_warpx_test(
    load_external_field_rz_grid # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_rz_test_load_external_field_rz_grid # inputs
    Examples/Tests/load_external_field/analysis_rz.py # analysis
    diags/diag1000300 # output
)

# load_external_field_rz_grid_restart #########################################
#
add_warpx_test(
    load_external_field_rz_grid_restart # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_rz_test_load_external_field_rz_grid_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000300 # output
)

# load_external_field_rz_particles ############################################
#
add_warpx_test(
    load_external_field_rz_particles # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_rz_test_load_external_field_rz_particles # inputs
    Examples/Tests/load_external_field/analysis_rz.py # analysis
    diags/diag1000300 # output
)

# load_external_field_rz_particles_restart ####################################
#
add_warpx_test(
    load_external_field_rz_particles_restart # name
    RZ # dims
    1  # nprocs
    OFF # eb
    Examples/Tests/load_external_field/inputs_rz_test_load_external_field_rz_particles_restart # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000300 # output
)
