# Add tests (alphabetical order) ##############################################
#

# dirichlet_bc ################################################################
#
add_warpx_test(
    dirichlet_bc # name
    2 # dims
    2 # nprocs
    Examples/Tests/electrostatic_dirichlet_bc/inputs_2d_test_dirichlet_bc # inputs
    Examples/Tests/electrostatic_dirichlet_bc/analysis.py # analysis
    diags/diag1000100 # output
)

# dirichlet_bc_picmi ##########################################################
#
add_warpx_test(
    dirichlet_bc_picmi # name
    2 # dims
    2 # nprocs
    Examples/Tests/electrostatic_dirichlet_bc/inputs_2d_test_dirichlet_bc_picmi.py # inputs
    Examples/Tests/electrostatic_dirichlet_bc/analysis.py # analysis
    diags/diag1000100 # output
)
