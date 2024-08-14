# Add tests (alphabetical order) ##############################################
#

# divb_cleaning_3d ############################################################
#
add_warpx_test(
    divb_cleaning_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/divb_cleaning/inputs_3d_test_divb_cleaning_3d # inputs
    Examples/Tests/divb_cleaning/analysis.py # analysis
    diags/diag1000400 # output
)
