# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
    add_warpx_test(
        test_rz_embedded_boundary_diffraction # name
        RZ # dims
        2  # nprocs
        ON # eb
        Examples/Tests/embedded_boundary_diffraction/inputs_test_rz_embedded_boundary_diffraction # inputs
        Examples/Tests/embedded_boundary_diffraction/analysis_fields.py # analysis
        diags/diag1/ # output
    )
endif()
