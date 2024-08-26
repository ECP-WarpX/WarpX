# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_rz_btd # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/btd_rz/inputs_test_rz_btd # inputs
    Examples/Tests/btd_rz/analysis.py # analysis
    diags/diag1000289 # output
    OFF # dependency
)
