# Add tests (alphabetical order) ##############################################
#

# particle_boundary_interaction_picmi #########################################
#
if(WarpX_EB)
    add_warpx_test(
        particle_boundary_interaction_picmi # name
        RZ # dims
        2  # nprocs
        ON # eb
        Examples/Tests/particle_boundary_interaction/inputs_rz_test_particle_boundary_interaction_picmi.py # inputs
        Examples/Tests/particle_boundary_interaction/analysis.py # analysis
        diags/diag1/ # output
    )
endif()
