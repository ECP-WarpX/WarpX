# Add tests (alphabetical order) ##############################################
#

# dive_cleaning_2d ############################################################
#
add_warpx_test(
    dive_cleaning_2d # name
    2  # dims
    ON # mpi
    Examples/Tests/dive_cleaning/inputs_2d_test_dive_cleaning_2d # inputs
    Examples/Tests/dive_cleaning/analysis.py # analysis
    diags/diag1000128 # output
)

# dive_cleaning_3d ############################################################
#
add_warpx_test(
    dive_cleaning_3d # name
    3  # dims
    ON # mpi
    Examples/Tests/dive_cleaning/inputs_3d_test_dive_cleaning_3d # inputs
    Examples/Tests/dive_cleaning/analysis.py # analysis
    diags/diag1000128 # output
)
