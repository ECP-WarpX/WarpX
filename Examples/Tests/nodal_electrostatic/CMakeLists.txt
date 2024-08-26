# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_nodal_electrostatic_solver # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/nodal_electrostatic/inputs_test_3d_nodal_electrostatic_solver # inputs
    Examples/Tests/nodal_electrostatic/analysis.py # analysis
    diags/diag1000010 # output
    OFF # depend on this test
)
