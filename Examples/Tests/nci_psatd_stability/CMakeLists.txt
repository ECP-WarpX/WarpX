# Add tests (alphabetical order) ##############################################
#

# averaged_galilean_2d_psatd ##################################################
#
if(WarpX_FFT)
    add_warpx_test(
        averaged_galilean_2d_psatd # name
        2   # dims
        OFF # mpi
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
        Examples/Tests/nci_psatd_stability/inputs_avg_3d_test_averaged_galilean_3d_psatd_hybrid # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000160 # output
    )
endif()

# comoving_2d_psatd_hybrid ####################################################
#
if(WarpX_FFT)
    add_warpx_test(
        comoving_2d_psatd_hybrid # name
        2  # dims
        ON # mpi
        Examples/Tests/nci_psatd_stability/inputs_2d_test_comoving_2d_psatd_hybrid # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000400 # output
    )
endif()

# galilean_2d_psatd ###########################################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_2d_psatd # name
        2   # dims
        OFF # mpi
        Examples/Tests/nci_psatd_stability/inputs_2d_test_galilean_2d_psatd # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
    )
endif()

# galilean_2d_psatd_current_correction ########################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_2d_psatd_current_correction # name
        2  # dims
        ON # mpi
        Examples/Tests/nci_psatd_stability/inputs_2d_test_galilean_2d_psatd_current_correction # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
    )
endif()

# galilean_2d_psatd_current_correction_psb ####################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_2d_psatd_current_correction_psb # name
        2   # dims
        OFF # mpi
        Examples/Tests/nci_psatd_stability/inputs_2d_test_galilean_2d_psatd_current_correction_psb # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
    )
endif()

# galilean_2d_psatd_hybrid ####################################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_2d_psatd_hybrid # name
        2  # dims
        ON # mpi
        Examples/Tests/nci_psatd_stability/inputs_2d_test_galilean_2d_psatd_hybrid # inputs
        Examples/analysis_default_regression.py # analysis
        diags/diag1000400 # output
    )
endif()

# galilean_3d_psatd ###########################################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_3d_psatd # name
        3  # dims
        ON # mpi
        Examples/Tests/nci_psatd_stability/inputs_3d_test_galilean_3d_psatd # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000300 # output
    )
endif()

# galilean_3d_psatd_current_correction ########################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_3d_psatd_current_correction # name
        3  # dims
        ON # mpi
        Examples/Tests/nci_psatd_stability/inputs_3d_test_galilean_3d_psatd_current_correction # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000300 # output
    )
endif()

# galilean_3d_psatd_current_correction_psb ####################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_3d_psatd_current_correction_psb # name
        3   # dims
        OFF # mpi
        Examples/Tests/nci_psatd_stability/inputs_3d_test_galilean_3d_psatd_current_correction_psb # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000300 # output
    )
endif()

# galilean_rz_psatd ###########################################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_rz_psatd # name
        RZ  # dims
        OFF # mpi
        Examples/Tests/nci_psatd_stability/inputs_rz_test_galilean_rz_psatd # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
    )
endif()

# galilean_rz_psatd_current_correction ########################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_rz_psatd_current_correction # name
        RZ # dims
        ON # mpi
        Examples/Tests/nci_psatd_stability/inputs_rz_test_galilean_rz_psatd_current_correction # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
    )
endif()

# galilean_rz_psatd_current_correction_psb ####################################
#
if(WarpX_FFT)
    add_warpx_test(
        galilean_rz_psatd_current_correction_psb # name
        RZ  # dims
        OFF # mpi
        Examples/Tests/nci_psatd_stability/inputs_rz_test_galilean_rz_psatd_current_correction_psb # inputs
        Examples/Tests/nci_psatd_stability/analysis_galilean.py # analysis
        diags/diag1000400 # output
    )
endif()
