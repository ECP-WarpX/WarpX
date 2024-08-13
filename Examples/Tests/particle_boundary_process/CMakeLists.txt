# Add tests (alphabetical order) ##############################################
#

# particle_absorption #########################################################
#
add_warpx_test(
    particle_absorption # name
    3 # dims
    2 # nprocs
    Examples/Tests/particle_boundary_process/inputs_3d_test_particle_absorption # inputs
    Examples/Tests/particle_boundary_process/analysis_absorption.py # analysis
    diags/diag1000060 # output
)

# particle_reflection_picmi ###################################################
#
add_warpx_test(
    particle_reflection_picmi # name
    2 # dims
    1 # nprocs
    Examples/Tests/particle_boundary_process/inputs_2d_test_particle_reflection_picmi.py # inputs
    Examples/Tests/particle_boundary_process/analysis_reflection.py # analysis
    diags/diag1000010 # output
)
