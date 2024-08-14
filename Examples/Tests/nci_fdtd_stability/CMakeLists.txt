# Add tests (alphabetical order) ##############################################
#

# nci_corrector ###############################################################
#
add_warpx_test(
    nci_corrector # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nci_fdtd_stability/inputs_2d_test_nci_corrector # inputs
    Examples/Tests/nci_fdtd_stability/analysis_ncicorr.py # analysis
    diags/diag1000600 # output
)

# nci_corrector_mr ############################################################
#
add_warpx_test(
    nci_corrector_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nci_fdtd_stability/inputs_2d_test_nci_corrector_mr # inputs
    Examples/Tests/nci_fdtd_stability/analysis_ncicorr.py # analysis
    diags/diag1000600 # output
)
