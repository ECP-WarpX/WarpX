# Add tests (alphabetical order) ##############################################
#

# open_bc_poisson_solver ######################################################
#
if(WarpX_FFT)
    add_warpx_test(
        open_bc_poisson_solver # name
        3  # dims
        ON # mpi
        Examples/Tests/open_bc_poisson_solver/inputs_3d_test_open_bc_poisson_solver # inputs
        Examples/Tests/open_bc_poisson_solver/analysis.py # analysis
        diags/diag1000001 # output
    )
endif()
