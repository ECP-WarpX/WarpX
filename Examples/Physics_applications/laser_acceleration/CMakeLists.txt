# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    Python_LaserAcceleration_1d # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_acceleration/PICMI_inputs_1d.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000100 # output
)
