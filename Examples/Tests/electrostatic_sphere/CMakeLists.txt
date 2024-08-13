# Add tests (alphabetical order) ##############################################
#

# electrostatic_sphere ########################################################
#
add_warpx_test(
    electrostatic_sphere # name
    3 # dims
    2 # nprocs
    Examples/Tests/electrostatic_sphere/inputs_3d_test_electrostatic_sphere # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000030 # output
)

# electrostatic_sphere_lab_frame ##############################################
#
add_warpx_test(
    electrostatic_sphere_lab_frame # name
    3 # dims
    2 # nprocs
    Examples/Tests/electrostatic_sphere/inputs_3d_test_electrostatic_sphere_lab_frame # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000030 # output
)

# electrostatic_sphere_lab_frame_MR_emass_10 ##################################
#
add_warpx_test(
    electrostatic_sphere_lab_frame_MR_emass_10 # name
    3 # dims
    2 # nprocs
    Examples/Tests/electrostatic_sphere/inputs_3d_test_electrostatic_sphere_lab_frame_MR_emass_10 # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000002 # output
)

# electrostatic_sphere_rel_nodal ##############################################
#
add_warpx_test(
    electrostatic_sphere_rel_nodal # name
    3 # dims
    2 # nprocs
    Examples/Tests/electrostatic_sphere/inputs_3d_test_electrostatic_sphere_rel_nodal # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000030 # output
)

# electrostatic_sphere_rz #####################################################
#
add_warpx_test(
    electrostatic_sphere_rz # name
    RZ # dims
    2  # nprocs
    Examples/Tests/electrostatic_sphere/inputs_rz_test_electrostatic_sphere_rz # inputs
    Examples/Tests/electrostatic_sphere/analysis_electrostatic_sphere.py # analysis
    diags/diag1000030 # output
)
