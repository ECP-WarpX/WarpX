# Add tests (alphabetical order) ##############################################
#

# embedded_boundary_python_api ################################################
#
add_warpx_test(
    embedded_boundary_picmi # name
    3 # dims
    1 # nprocs
    Examples/Tests/embedded_boundary_python_api/inputs_3d_test_embedded_boundary_picmi.py # inputs
    Examples/Tests/embedded_boundary_python_api/analysis.py # analysis
    diags/diag1000002
)
