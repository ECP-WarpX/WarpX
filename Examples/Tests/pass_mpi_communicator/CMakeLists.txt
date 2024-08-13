# Add tests (alphabetical order) ##############################################
#

# pass_mpi_comm_picmi #########################################################
#
add_warpx_test(
    pass_mpi_comm_picmi # name
    2 # dims
    2 # nprocs
    Examples/Tests/pass_mpi_communicator/inputs_2d_test_pass_mpi_comm_picmi.py # inputs
    Examples/Tests/pass_mpi_communicator/analysis.py # analysis
    diags/diag1000010 # output
)
