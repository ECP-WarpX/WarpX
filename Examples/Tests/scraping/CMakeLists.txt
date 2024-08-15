# Add tests (alphabetical order) ##############################################
#

# scraping ####################################################################
#
if(WarpX_EB)
    add_warpx_test(
        scraping # name
        RZ # dims
        2  # nprocs
        ON # eb
        Examples/Tests/scraping/inputs_rz_test_scraping # inputs
        Examples/Tests/scraping/analysis_rz.py # analysis
        diags/diag1000037 # output
    )
endif()

# scraping_filter #############################################################
#
if(WarpX_EB)
    add_warpx_test(
        scraping_filter # name
        RZ # dims
        2  # nprocs
        ON # eb
        Examples/Tests/scraping/inputs_rz_test_scraping_filter # inputs
        Examples/Tests/scraping/analysis_rz_filter.py # analysis
        diags/diag1000037 # output
    )
endif()
