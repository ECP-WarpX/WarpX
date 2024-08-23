# Add tests (alphabetical order) ##############################################
#

if(WarpX_FFT)
    add_warpx_test(
        test_3d_open_bc_poisson_solver # name
        3 # dims
        2 # nprocs
        OFF # eb
    OFF # preparation script
        Examples/Tests/open_bc_poisson_solver/inputs_test_3d_open_bc_poisson_solver # inputs
        Examples/Tests/open_bc_poisson_solver/analysis.py # analysis
        diags/diag1000001 # output
    )
endif()
