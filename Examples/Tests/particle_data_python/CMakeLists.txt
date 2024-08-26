# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_prev_positions_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/particle_data_python/inputs_test_2d_prev_positions_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000010 # output
    OFF # dependency
)

add_warpx_test(
    test_2d_particle_attr_access_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/particle_data_python/inputs_test_2d_particle_attr_access_picmi.py # inputs
    Examples/Tests/particle_data_python/analysis.py # analysis
    diags/diag1000010 # output
    OFF # dependency
)

# FIXME
#add_warpx_test(
#    test_2d_particle_attr_access_unique_picmi # name
#    2 # dims
#    2 # nprocs
#    OFF # eb
#    Examples/Tests/particle_data_python/inputs_test_2d_particle_attr_access_picmi.py # inputs
#    Examples/Tests/particle_data_python/analysis.py # analysis
#    diags/diag1000010 # output
#    OFF # dependency
#)
