# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_diff_lumi_diag  # name
    3  # dims
    2  # nprocs
    inputs_test_3d_diff_lumi_diag  # inputs
    analysis.py  # analysis
    diags/diag1000080  # output
    OFF  # dependency
)
