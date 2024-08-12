# Add tests (alphabetical order) ##############################################
#

# nodal_electrostatic_solver ##################################################
#
add_warpx_test(
    nodal_electrostatic_solver # name
    3   # dims
    OFF # mpi
    Examples/Tests/nodal_electrostatic/inputs_3d_test_nodal_electrostatic_solver # inputs
    Examples/Tests/nodal_electrostatic/analysis_3d.py # analysis
    diags/diag1000010 # output
)
