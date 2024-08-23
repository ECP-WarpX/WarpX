# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_nci_corrector # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/nci_fdtd_stability/inputs_test_2d_nci_corrector # inputs
    Examples/Tests/nci_fdtd_stability/analysis_ncicorr.py # analysis
    diags/diag1000600 # output
)

add_warpx_test(
    test_2d_nci_corrector_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/nci_fdtd_stability/inputs_test_2d_nci_corrector_mr # inputs
    Examples/Tests/nci_fdtd_stability/analysis_ncicorr.py # analysis
    diags/diag1000600 # output
)
