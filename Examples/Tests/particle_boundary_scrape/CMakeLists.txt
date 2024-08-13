# Add tests (alphabetical order) ##############################################
#

# particle_scrape #############################################################
#
add_warpx_test(
    particle_scrape # name
    3 # dims
    2 # nprocs
    Examples/Tests/particle_boundary_scrape/inputs_3d_test_particle_scrape # inputs
    Examples/Tests/particle_boundary_scrape/analysis_scrape.py # analysis
    diags/diag1000060 # output
)

# particle_scrape_picmi #######################################################
#
add_warpx_test(
    particle_scrape_picmi # name
    3 # dims
    2 # nprocs
    Examples/Tests/particle_boundary_scrape/inputs_3d_test_particle_scrape_picmi.py # inputs
    Examples/Tests/particle_boundary_scrape/analysis_scrape.py # analysis
    diags/diag1000060 # output
)
