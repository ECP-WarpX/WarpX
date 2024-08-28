# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_diff_lumi_diag # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/diff_lumi_diag/inputs_test_3d_diff_lumi_diag # inputs
    Examples/Tests/diff_lumi_diag/analysis.py # analysis
    diags/diag1000080 # output
    OFF # dependency
)
