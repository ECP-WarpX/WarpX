# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_collider_diagnostics # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/collider_relevant_diags/inputs_test_3d_collider_diagnostics # inputs
    Examples/Tests/collider_relevant_diags/analysis.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)
