# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
    add_warpx_test(
        test_3d_particle_scrape # name
        3 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/particle_boundary_scrape/inputs_test_3d_particle_scrape # inputs
        Examples/Tests/particle_boundary_scrape/analysis_scrape.py # analysis
        diags/diag1000060 # output
    OFF # depend on this test
    )
endif()

if(WarpX_EB)
    add_warpx_test(
        test_3d_particle_scrape_picmi # name
        3 # dims
        2 # nprocs
        ON # eb
        Examples/Tests/particle_boundary_scrape/inputs_test_3d_particle_scrape_picmi.py # inputs
        Examples/Tests/particle_boundary_scrape/analysis_scrape.py # analysis
        diags/diag1000060 # output
    OFF # depend on this test
    )
endif()
