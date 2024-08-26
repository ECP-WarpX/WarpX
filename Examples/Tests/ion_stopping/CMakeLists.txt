# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_ion_stopping # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/ion_stopping/inputs_test_3d_ion_stopping # inputs
    Examples/Tests/ion_stopping/analysis.py # analysis
    diags/diag1000010 # output
    OFF # dependency
)
