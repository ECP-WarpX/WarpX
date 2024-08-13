# Add tests (alphabetical order) ##############################################
#

# plasma_lens #################################################################
#
add_warpx_test(
    plasma_lens # name
    3 # dims
    2 # nprocs
    Examples/Tests/plasma_lens/inputs_3d_test_plasma_lens # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)

# plasma_lens_boosted #########################################################
#
add_warpx_test(
    plasma_lens_boosted # name
    3 # dims
    2 # nprocs
    Examples/Tests/plasma_lens/inputs_3d_test_plasma_lens_boosted # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)

# plasma_lens_hard_edged ######################################################
#
add_warpx_test(
    plasma_lens_hard_edged # name
    3 # dims
    2 # nprocs
    Examples/Tests/plasma_lens/inputs_3d_test_plasma_lens_hard_edged # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)

# plasma_lens_picmi ###########################################################
#
add_warpx_test(
    plasma_lens_picmi # name
    3 # dims
    2 # nprocs
    Examples/Tests/plasma_lens/inputs_3d_test_plasma_lens_picmi.py # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)

# plasma_lens_short ###########################################################
#
add_warpx_test(
    plasma_lens_short # name
    3 # dims
    2 # nprocs
    Examples/Tests/plasma_lens/inputs_3d_test_plasma_lens_short # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)
