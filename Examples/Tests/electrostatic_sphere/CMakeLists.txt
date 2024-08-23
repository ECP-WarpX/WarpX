# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_electrostatic_sphere # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/electrostatic_sphere/inputs_test_3d_electrostatic_sphere # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000030 # output
)

add_warpx_test(
    test_3d_electrostatic_sphere_lab_frame # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/electrostatic_sphere/inputs_test_3d_electrostatic_sphere_lab_frame # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000030 # output
)

add_warpx_test(
    test_3d_electrostatic_sphere_lab_frame_mr_emass_10 # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/electrostatic_sphere/inputs_test_3d_electrostatic_sphere_lab_frame_mr_emass_10 # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000002 # output
)

add_warpx_test(
    test_3d_electrostatic_sphere_rel_nodal # name
    3 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/electrostatic_sphere/inputs_test_3d_electrostatic_sphere_rel_nodal # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000030 # output
)

add_warpx_test(
    test_rz_electrostatic_sphere # name
    RZ # dims
    2  # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/electrostatic_sphere/inputs_test_rz_electrostatic_sphere # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000030 # output
)
