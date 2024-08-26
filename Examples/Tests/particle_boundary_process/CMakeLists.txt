# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
    add_warpx_test(
        test_3d_particle_absorption # name
        3 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/particle_boundary_process/inputs_test_3d_particle_absorption # inputs
        Examples/Tests/particle_boundary_process/analysis_absorption.py # analysis
        diags/diag1000060 # output
    OFF # depend on this test
    )
endif()

add_warpx_test(
    test_2d_particle_reflection_picmi # name
    2 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/particle_boundary_process/inputs_test_2d_particle_reflection_picmi.py # inputs
    Examples/Tests/particle_boundary_process/analysis_reflection.py # analysis
    diags/diag1000010 # output
    OFF # depend on this test
)
