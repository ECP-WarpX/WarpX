# Add tests (alphabetical order) ##############################################
#

# spacecraft_charging #########################################################
#
if(WarpX_EB)
add_warpx_test(
    spacecraft_charging # name
    RZ # dims
    2  # nprocs
    ON # eb
    Examples/Physics_applications/spacecraft_charging/inputs_rz_test_spacecraft_charging_picmi.py # inputs
    Examples/Physics_applications/spacecraft_charging/analysis.py # analysis
    diags/diag1/ # output
)
endif()
