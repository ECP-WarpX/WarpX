# Add tests (alphabetical order) ##############################################
#

# laser_ion_acc_2d ############################################################
#
add_warpx_test(
    laser_ion_acc_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_ion/inputs_2d_test_laser_ion_acc_2d # inputs
    Examples/analysis_default_openpmd_regression.py # analysis
    diags/diag1/ # output
)

# laser_ion_acc_2d_picmi ######################################################
#
add_warpx_test(
    laser_ion_acc_2d_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Physics_applications/laser_ion/inputs_2d_test_laser_ion_acc_2d_picmi.py # inputs
    Examples/analysis_default_openpmd_regression.py # analysis
    diags/diag1/ # output
)
