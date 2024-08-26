# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
add_warpx_test(
    test_rz_spacecraft_charging # name
    RZ # dims
    2  # nprocs
    ON # eb
    Examples/Physics_applications/spacecraft_charging/inputs_test_rz_spacecraft_charging_picmi.py # inputs
    Examples/Physics_applications/spacecraft_charging/analysis.py # analysis
    diags/diag1/ # output
    OFF # dependency
)
endif()
