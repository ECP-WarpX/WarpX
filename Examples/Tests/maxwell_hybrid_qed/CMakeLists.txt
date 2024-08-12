# Add tests (alphabetical order) ##############################################
#

# maxwell_hybrid_qed_solver ###################################################
#
if(WarpX_FFT)
    add_warpx_test(
        maxwell_hybrid_qed_solver # name
        2  # dims
        ON # mpi
        Examples/Tests/maxwell_hybrid_qed/inputs_2d_test_maxwell_hybrid_qed_solver # inputs
        Examples/Tests/maxwell_hybrid_qed/analysis_2d.py # analysis
        diags/diag1000300 # output
    )
endif()
