# Add tests (alphabetical order) ##############################################
#

# radiation_reaction ##########################################################
#
add_warpx_test(
    radiation_reaction # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/radiation_reaction/inputs_3d_test_radiation_reaction # inputs
    Examples/Tests/radiation_reaction/analysis.py # analysis
    diags/diag1000064 # output
)
