# Add tests (alphabetical order) ##############################################
#

# flux_injection_3d ###########################################################
#
add_warpx_test(
    flux_injection_3d # name
    3 # dims
    2 # nprocs
    Examples/Tests/flux_injection/inputs_3d_test_flux_injection_3d # inputs
    Examples/Tests/flux_injection/analysis_flux_injection_3d.py # analysis
    diags/diag1000002 # output
)

# flux_injection_rz ###########################################################
#
add_warpx_test(
    flux_injection_rz # name
    RZ # dims
    2  # nprocs
    Examples/Tests/flux_injection/inputs_rz_test_flux_injection_rz # inputs
    Examples/Tests/flux_injection/analysis_flux_injection_rz.py # analysis
    diags/diag1000120 # output
)
