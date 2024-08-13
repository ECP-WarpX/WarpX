# Add tests (alphabetical order) ##############################################
#

# particle_fields_diags #######################################################
#
add_warpx_test(
    particle_fields_diags # name
    3 # dims
    2 # nprocs
    Examples/Tests/particle_fields_diags/inputs_3d_test_particle_fields_diags # inputs
    Examples/Tests/particle_fields_diags/analysis_particle_diags.py # analysis
    diags/diag1000200 # output
)

# FIXME
# particle_fields_diags_single_precision ######################################
#
#add_warpx_test(
#    particle_fields_diags_single_precision # name
#    3 # dims
#    2 # nprocs
#    Examples/Tests/particle_fields_diags/inputs_3d_test_particle_fields_diags_single_precision # inputs
#    Examples/Tests/particle_fields_diags/analysis_particle_diags_single.py # analysis
#    diags/diag1000200 # output
#)
