# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_plasma_lens # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/plasma_lens/inputs_test_3d_plasma_lens # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)

add_warpx_test(
    test_3d_plasma_lens_boosted # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/plasma_lens/inputs_test_3d_plasma_lens_boosted # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)

add_warpx_test(
    test_3d_plasma_lens_hard_edged # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/plasma_lens/inputs_test_3d_plasma_lens_hard_edged # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)

add_warpx_test(
    test_3d_plasma_lens_picmi # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/plasma_lens/inputs_test_3d_plasma_lens_picmi.py # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)

add_warpx_test(
    test_3d_plasma_lens_short # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/plasma_lens/inputs_test_3d_plasma_lens_short # inputs
    Examples/Tests/plasma_lens/analysis.py # analysis
    diags/diag1000084 # output
)
