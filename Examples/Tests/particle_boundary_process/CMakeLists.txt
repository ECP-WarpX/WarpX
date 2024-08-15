# Add tests (alphabetical order) ##############################################
#

# particle_absorption #########################################################
#
if(WarpX_EB)
    add_warpx_test(
        particle_absorption # name
        3 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/particle_boundary_process/inputs_3d_test_particle_absorption # inputs
        Examples/Tests/particle_boundary_process/analysis_absorption.py # analysis
        diags/diag1000060 # output
    )
endif()

# particle_reflection_picmi ###################################################
#
add_warpx_test(
    particle_reflection_picmi # name
    2 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/particle_boundary_process/inputs_2d_test_particle_reflection_picmi.py # inputs
    Examples/Tests/particle_boundary_process/analysis_reflection.py # analysis
    diags/diag1000010 # output
)
