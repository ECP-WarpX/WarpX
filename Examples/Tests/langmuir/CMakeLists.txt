# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_1d_langmuir_multi # name
    1 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_1d_langmuir_multi # inputs
    Examples/Tests/langmuir/analysis_1d.py # analysis
    diags/diag1000080 # output
)

add_warpx_test(
    test_2d_langmuir_multi_mr # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_mr # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

add_warpx_test(
    test_2d_langmuir_multi_mr_anisotropic # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_mr_anisotropic # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

add_warpx_test(
    test_2d_langmuir_multi_mr_momentum_conserving # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_mr_momentum_conserving # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_mr_psatd # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_mr_psatd # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

add_warpx_test(
    test_2d_langmuir_multi_nodal # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_nodal # inputs
    Examples/Tests/langmuir/analysis_2d.py # analysis
    diags/diag1000080 # output
)

add_warpx_test(
    test_2d_langmuir_multi_picmi # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd_current_correction # name
        2 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd_current_correction_nodal # name
        2 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_current_correction_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd_momentum_conserving # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_momentum_conserving # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd_multiJ # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd_multiJ_nodal # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_multiJ_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd_nodal # name
        2 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd_vay_deposition # name
        2 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_vay_deposition # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

# FIXME
#if(WarpX_FFT)
#    add_warpx_test(
#        test_2d_langmuir_multi_psatd_vay_deposition_particle_shape_4 # name
#        2 # dims
#        1 # nprocs
#        OFF # eb
#        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_vay_deposition_particle_shape_4 # inputs
#        Examples/Tests/langmuir/analysis_2d.py # analysis
#        diags/diag1000080 # output
#    )
#endif()

if(WarpX_FFT)
    add_warpx_test(
        test_2d_langmuir_multi_psatd_vay_deposition_nodal # name
        2 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_2d_langmuir_multi_psatd_vay_deposition_nodal # inputs
        Examples/Tests/langmuir/analysis_2d.py # analysis
        diags/diag1000080 # output
    )
endif()

add_warpx_test(
    test_3d_langmuir_multi # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_3d_langmuir_multi # inputs
    Examples/Tests/langmuir/analysis_3d.py # analysis
    diags/diag1000040 # output
)

add_warpx_test(
    test_3d_langmuir_multi_nodal # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_nodal # inputs
    Examples/Tests/langmuir/analysis_3d.py # analysis
    diags/diag1000040 # output
)

add_warpx_test(
    test_3d_langmuir_multi_picmi # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_current_correction # name
        3 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_current_correction_nodal # name
        3 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_current_correction_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_div_cleaning # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_div_cleaning # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_momentum_conserving # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_momentum_conserving # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_multiJ # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_multiJ_nodal # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_multiJ_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_nodal # name
        3 # dims
        2 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_vay_deposition # name
        3 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_vay_deposition # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_3d_langmuir_multi_psatd_vay_deposition_nodal # name
        3 # dims
        1 # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_3d_langmuir_multi_psatd_vay_deposition_nodal # inputs
        Examples/Tests/langmuir/analysis_3d.py # analysis
        diags/diag1000040 # output
    )
endif()

add_warpx_test(
    test_rz_langmuir_multi # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_rz_langmuir_multi # inputs
    Examples/Tests/langmuir/analysis_rz.py # analysis
    diags/diag1000080 # output
)

add_warpx_test(
    test_rz_langmuir_multi_picmi # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/langmuir/inputs_test_rz_langmuir_multi_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000040 # output
)

if(WarpX_FFT)
    add_warpx_test(
        test_rz_langmuir_multi_psatd # name
        RZ # dims
        2  # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_rz_langmuir_multi_psatd # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_rz_langmuir_multi_psatd_current_correction # name
        RZ # dims
        1  # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_rz_langmuir_multi_psatd_current_correction # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()

if(WarpX_FFT)
    add_warpx_test(
        test_rz_langmuir_multi_psatd_multiJ # name
        RZ # dims
        2  # nprocs
        OFF # eb
        Examples/Tests/langmuir/inputs_test_rz_langmuir_multi_psatd_multiJ # inputs
        Examples/Tests/langmuir/analysis_rz.py # analysis
        diags/diag1000080 # output
    )
endif()
