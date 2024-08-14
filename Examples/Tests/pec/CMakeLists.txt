# Add tests (alphabetical order) ##############################################
#

# pec_field ###################################################################
#
add_warpx_test(
    pec_field # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pec/inputs_3d_test_pec_field # inputs
    Examples/Tests/pec/analysis_pec.py # analysis
    diags/diag1000125 # output
)

# pec_field_mr ################################################################
#
add_warpx_test(
    pec_field_mr # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pec/inputs_3d_test_pec_field_mr # inputs
    Examples/Tests/pec/analysis_pec_mr.py # analysis
    diags/diag1000125 # output
)

# pec_particle ################################################################
#
add_warpx_test(
    pec_particle # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/pec/inputs_3d_test_pec_particle # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000020 # output
)
