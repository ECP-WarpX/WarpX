# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_ionization_boost # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/ionization/inputs_test_2d_ionization_boost # inputs
    Examples/Tests/ionization/analysis.py # analysis
    diags/diag1000420 # output
    OFF # dependency
)

add_warpx_test(
    test_2d_ionization_lab # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/ionization/inputs_test_2d_ionization_lab # inputs
    Examples/Tests/ionization/analysis.py # analysis
    diags/diag1001600 # output
    OFF # dependency
)

add_warpx_test(
    test_2d_ionization_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/ionization/inputs_test_2d_ionization_picmi.py # inputs
    Examples/Tests/ionization/analysis.py # analysis
    diags/diag1001600 # output
    OFF # dependency
)
