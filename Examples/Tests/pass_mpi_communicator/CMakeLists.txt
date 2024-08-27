# Add tests (alphabetical order) ##############################################
#

# TODO
# - Enable in pyAMReX (https://github.com/AMReX-Codes/pyamrex/issues/163)
# - Enable related lines in input script
# - Enable analysis script
add_warpx_test(
    test_2d_pass_mpi_comm_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pass_mpi_communicator/inputs_test_2d_pass_mpi_comm_picmi.py # inputs
    OFF #Examples/Tests/pass_mpi_communicator/analysis.py # analysis
    OFF # output
    OFF # dependency
)
