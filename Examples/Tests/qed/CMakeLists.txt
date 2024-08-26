# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_qed_breit_wheeler # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_2d_qed_breit_wheeler # inputs
    Examples/Tests/qed/analysis_breit_wheeler_yt.py # analysis
    diags/diag1000002 # output
    OFF # dependency
)

add_warpx_test(
    test_2d_qed_breit_wheeler_opmd # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_2d_qed_breit_wheeler_opmd # inputs
    Examples/Tests/qed/analysis_breit_wheeler_opmd.py # analysis
    diags/diag1/ # output
    OFF # dependency
)

add_warpx_test(
    test_3d_qed_breit_wheeler # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_3d_qed_breit_wheeler # inputs
    Examples/Tests/qed/analysis_breit_wheeler_yt.py # analysis
    diags/diag1000002 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_qed_breit_wheeler_opmd # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_3d_qed_breit_wheeler_opmd # inputs
    Examples/Tests/qed/analysis_breit_wheeler_opmd.py # analysis
    diags/diag1/ # output
    OFF # dependency
)

add_warpx_test(
    test_2d_qed_quantum_sync # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_2d_qed_quantum_sync # inputs
    Examples/Tests/qed/analysis_quantum_sync.py # analysis
    diags/diag1000002 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_qed_quantum_sync # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_3d_qed_quantum_sync # inputs
    Examples/Tests/qed/analysis_quantum_sync.py # analysis
    diags/diag1000002 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_qed_schwinger_1 # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_3d_qed_schwinger_1 # inputs
    Examples/Tests/qed/analysis_schwinger.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_qed_schwinger_2 # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_3d_qed_schwinger_2 # inputs
    Examples/Tests/qed/analysis_schwinger.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_qed_schwinger_3 # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_3d_qed_schwinger_3 # inputs
    Examples/Tests/qed/analysis_schwinger.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)

add_warpx_test(
    test_3d_qed_schwinger_4 # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_test_3d_qed_schwinger_4 # inputs
    Examples/Tests/qed/analysis_schwinger.py # analysis
    diags/diag1000001 # output
    OFF # dependency
)
