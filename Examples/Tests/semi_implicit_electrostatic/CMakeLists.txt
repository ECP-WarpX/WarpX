# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_semi_implicit_electrostatic_picmi  # name
    2  # dims
    2  # nprocs
    inputs_test_2d_semi_implicit_electrostatic_picmi.py  # inputs
    analysis.py  # analysis
    diags/field_diag000190  # output
    OFF  # dependency
)
