# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_radiation_reaction # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/radiation_reaction/inputs_test_3d_radiation_reaction # inputs
    Examples/Tests/radiation_reaction/analysis.py # analysis
    diags/diag1000064 # output
)
