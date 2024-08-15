# Add tests (alphabetical order) ##############################################
#

# qed_breit_wheeler_2d ########################################################
#
add_warpx_test(
    qed_breit_wheeler_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_2d_test_qed_breit_wheeler_2d # inputs
    Examples/Tests/qed/analysis_breit_wheeler_yt.py # analysis
    diags/diag1000002 # output
)

# qed_breit_wheeler_2d_opmd ###################################################
#
add_warpx_test(
    qed_breit_wheeler_2d_opmd # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_2d_test_qed_breit_wheeler_2d_opmd # inputs
    Examples/Tests/qed/analysis_breit_wheeler_opmd.py # analysis
    diags/diag1/ # output
)

# qed_breit_wheeler_3d ########################################################
#
add_warpx_test(
    qed_breit_wheeler_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_3d_test_qed_breit_wheeler_3d # inputs
    Examples/Tests/qed/analysis_breit_wheeler_yt.py # analysis
    diags/diag1000002 # output
)

# qed_breit_wheeler_3d_opmd ###################################################
#
add_warpx_test(
    qed_breit_wheeler_3d_opmd # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_3d_test_qed_breit_wheeler_3d_opmd # inputs
    Examples/Tests/qed/analysis_breit_wheeler_opmd.py # analysis
    diags/diag1/ # output
)

# qed_quantum_sync_2d #########################################################
#
add_warpx_test(
    qed_quantum_sync_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_2d_test_qed_quantum_sync_2d # inputs
    Examples/Tests/qed/analysis_quantum_sync.py # analysis
    diags/diag1000002 # output
)

# qed_quantum_sync_3d #########################################################
#
add_warpx_test(
    qed_quantum_sync_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_3d_test_qed_quantum_sync_3d # inputs
    Examples/Tests/qed/analysis_quantum_sync.py # analysis
    diags/diag1000002 # output
)

# qed_schwinger_1 #############################################################
#
add_warpx_test(
    qed_schwinger_1 # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_3d_test_qed_schwinger_1 # inputs
    Examples/Tests/qed/analysis_schwinger.py # analysis
    diags/diag1000001 # output
)

# qed_schwinger_2 #############################################################
#
add_warpx_test(
    qed_schwinger_2 # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_3d_test_qed_schwinger_2 # inputs
    Examples/Tests/qed/analysis_schwinger.py # analysis
    diags/diag1000001 # output
)

# qed_schwinger_3 #############################################################
#
add_warpx_test(
    qed_schwinger_3 # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_3d_test_qed_schwinger_3 # inputs
    Examples/Tests/qed/analysis_schwinger.py # analysis
    diags/diag1000001 # output
)

# qed_schwinger_4 #############################################################
#
add_warpx_test(
    qed_schwinger_4 # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/qed/inputs_3d_test_qed_schwinger_4 # inputs
    Examples/Tests/qed/analysis_schwinger.py # analysis
    diags/diag1000001 # output
)
