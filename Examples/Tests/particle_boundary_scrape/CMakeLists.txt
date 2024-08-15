# Add tests (alphabetical order) ##############################################
#

# particle_scrape #############################################################
#
if(WarpX_EB)
    add_warpx_test(
        particle_scrape # name
        3 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/particle_boundary_scrape/inputs_3d_test_particle_scrape # inputs
        Examples/Tests/particle_boundary_scrape/analysis_scrape.py # analysis
        diags/diag1000060 # output
    )
endif()

# particle_scrape_picmi #######################################################
#
if(WarpX_EB)
    add_warpx_test(
        particle_scrape_picmi # name
        3 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/particle_boundary_scrape/inputs_3d_test_particle_scrape_picmi.py # inputs
        Examples/Tests/particle_boundary_scrape/analysis_scrape.py # analysis
        diags/diag1000060 # output
    )
endif()
