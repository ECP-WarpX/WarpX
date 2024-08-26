# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_dirichlet_bc # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/electrostatic_dirichlet_bc/inputs_test_2d_dirichlet_bc # inputs
    Examples/Tests/electrostatic_dirichlet_bc/analysis.py # analysis
    diags/diag1000100 # output
    OFF # depend on this test
)

add_warpx_test(
    test_2d_dirichlet_bc_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/electrostatic_dirichlet_bc/inputs_test_2d_dirichlet_bc_picmi.py # inputs
    Examples/Tests/electrostatic_dirichlet_bc/analysis.py # analysis
    diags/diag1000100 # output
    OFF # depend on this test
)
