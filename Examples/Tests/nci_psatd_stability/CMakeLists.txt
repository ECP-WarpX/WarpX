# Add tests (alphabetical order) ##############################################
#

# averaged_galilean_2d_psatd ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_2d_psatd # name
        2  # dims
        ON # mpi
        1  # mpi_np
        Examples/Tests/nci_psatd_stability/inputs_avg_2d_test_averaged_galilean_2d_psatd # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
    )
endif()

# averaged_galilean_2d_psatd_hybrid ###########################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_2d_psatd_hybrid # name
        2  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/nci_psatd_stability/inputs_avg_2d_test_averaged_galilean_2d_psatd_hybrid # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
    )
endif()

# averaged_galilean_3d_psatd ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_3d_psatd # name
        3  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/nci_psatd_stability/inputs_avg_3d_test_averaged_galilean_3d_psatd # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000160 # output
    )
endif()

# averaged_galilean_3d_psatd_hybrid ###########################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_3d_psatd_hybrid # name
        3  # dims
        ON # mpi
        2  # mpi_np
        Examples/Tests/nci_psatd_stability/inputs_avg_3d_test_averaged_galilean_3d_psatd_hybrid # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000160 # output
    )
endif()
