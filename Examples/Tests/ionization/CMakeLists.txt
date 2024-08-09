# Add tests (alphabetical order) ##############################################
#

# ionization_boost ############################################################
#
add_warpx_test(
    ionization_boost # name
    2  # dims
    ON # mpi
    Examples/Tests/ionization/inputs_2d_test_ionization_boost # inputs
    Examples/Tests/ionization/analysis_ionization.py # analysis
    diags/diag1000420 # output
)

# ionization_lab ##############################################################
#
add_warpx_test(
    ionization_lab # name
    2  # dims
    ON # mpi
    Examples/Tests/ionization/inputs_2d_test_ionization_lab # inputs
    Examples/Tests/ionization/analysis_ionization.py # analysis
    diags/diag1001600 # output
)

# ionization_picmi ##############################################################
#
add_warpx_test(
    ionization_picmi # name
    2  # dims
    ON # mpi
    Examples/Tests/ionization/inputs_2d_test_ionization_picmi.py # inputs
    Examples/Tests/ionization/analysis_ionization.py # analysis
    diags/diag1001600 # output
)
