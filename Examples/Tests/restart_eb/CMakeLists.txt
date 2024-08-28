# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_eb_picmi # name
    3 # dims
    1 # nprocs
    ON # eb
    Examples/Tests/restart_eb/inputs_test_3d_eb_picmi.py # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000060 # output
    OFF # dependency
)

# FIXME
#add_warpx_test(
#    test_3d_eb_picmi_restart # name
#    3 # dims
#    1 # nprocs
#    ON # eb
#    "Examples/Tests/restart_eb/inputs_test_3d_eb_picmi.py amr.restart=../test_3d_eb_picmi/diags/chk000030" # inputs
#    Examples/analysis_default_regression.py # analysis
#    diags/diag1000060 # output
#    test_3d_eb_picmi # dependency
#)
