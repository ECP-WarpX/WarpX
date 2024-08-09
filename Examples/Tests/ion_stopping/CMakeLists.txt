# Add tests (alphabetical order) ##############################################
#

# ion_stopping ################################################################
#
add_warpx_test(
    ion_stopping # name
    3   # dims
    OFF # mpi
    Examples/Tests/ion_stopping/inputs_2d_test_ion_stopping # inputs
    Examples/Tests/ion_stopping/analysis_ion_stopping.py # analysis
    diags/diag1000010 # output
)
