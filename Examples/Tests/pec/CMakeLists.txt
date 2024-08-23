# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_pec_field # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pec/inputs_test_3d_pec_field # inputs
    Examples/Tests/pec/analysis_pec.py # analysis
    diags/diag1000125 # output
)

add_warpx_test(
    test_3d_pec_field_mr # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pec/inputs_test_3d_pec_field_mr # inputs
    Examples/Tests/pec/analysis_pec_mr.py # analysis
    diags/diag1000125 # output
)

add_warpx_test(
    test_3d_pec_particle # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/pec/inputs_test_3d_pec_particle # inputs
    Examples/analysis_default_regression.py # analysis
    diags/diag1000020 # output
)
