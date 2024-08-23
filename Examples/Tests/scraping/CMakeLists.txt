# Add tests (alphabetical order) ##############################################
#

if(WarpX_EB)
    add_warpx_test(
        test_rz_scraping # name
        RZ # dims
        2  # nprocs
        ON # eb
    OFF # preparation script
        Examples/Tests/scraping/inputs_test_rz_scraping # inputs
        Examples/Tests/scraping/analysis_rz.py # analysis
        diags/diag1000037 # output
    )
endif()

if(WarpX_EB)
    add_warpx_test(
        test_rz_scraping_filter # name
        RZ # dims
        2  # nprocs
        ON # eb
    OFF # preparation script
        Examples/Tests/scraping/inputs_test_rz_scraping_filter # inputs
        Examples/Tests/scraping/analysis_rz_filter.py # analysis
        diags/diag1000037 # output
    )
endif()
