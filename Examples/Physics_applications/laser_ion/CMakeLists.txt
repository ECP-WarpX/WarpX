# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_laser_ion_acc  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_laser_ion_acc  # inputs
    analysis_default_openpmd_regression.py  # analysis
    diags/diagInst/  # output
    OFF  # dependency
)

add_warpx_test(
    test_2d_laser_ion_acc_time_avg  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_laser_ion_acc  # inputs; will not be used b/c dependent test
    analysis_time_avg_diags.py  # analysis
    diags/diagTimeAvg/  # output
    test_2d_laser_ion_acc  # dependency
)

add_warpx_test(
    test_2d_laser_ion_acc_picmi  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_laser_ion_acc_picmi.py  # inputs
    analysis_default_openpmd_regression.py  # analysis
    diags/diagInst/  # output
    OFF  # dependency
)

# We do not have an intervals parser for PICMI yet which would accept more than a single integer for the output period
#add_warpx_test(
#    test_2d_laser_ion_acc_time_avg_picmi  # name
#    2  # dims
#    2  # nprocs
#    inputs_test_2d_laser_ion_acc_picmi.py  # inputs
#    analysis_time_avg_diags.py  # analysis
#    diags/diagTimeAvg/  # output
#    test_2d_laser_ion_acc_picmi  # dependency
#)
