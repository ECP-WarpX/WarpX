# Add tests (alphabetical order) ##############################################
#

# deuterium_deuterium_fusion_3d ###############################################
#
add_warpx_test(
    deuterium_deuterium_fusion_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_3d_test_deuterium_deuterium_fusion_3d # inputs
    Examples/Tests/nuclear_fusion/analysis_two_product_fusion.py # analysis
    diags/diag1000001 # output
)

# deuterium_deuterium_fusion_3d_intraspecies ##################################
#
add_warpx_test(
    deuterium_deuterium_fusion_3d_intraspecies # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_3d_test_deuterium_deuterium_fusion_3d_intraspecies # inputs
    Examples/Tests/nuclear_fusion/analysis_deuterium_deuterium_3d_intraspecies.py # analysis
    diags/diag1000010 # output
)

# deuterium_tritium_fusion_3d #################################################
#
add_warpx_test(
    deuterium_tritium_fusion_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_3d_test_deuterium_tritium_fusion_3d # inputs
    Examples/Tests/nuclear_fusion/analysis_two_product_fusion.py # analysis
    diags/diag1000001 # output
)

# deuterium_tritium_fusion_rz #################################################
#
add_warpx_test(
    deuterium_tritium_fusion_rz # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_rz_test_deuterium_tritium_fusion_rz # inputs
    Examples/Tests/nuclear_fusion/analysis_two_product_fusion.py # analysis
    diags/diag1000001 # output
)

# proton_boron_fusion_2d ######################################################
#
add_warpx_test(
    proton_boron_fusion_2d # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_2d_test_proton_boron_fusion_2d # inputs
    Examples/Tests/nuclear_fusion/analysis_proton_boron_fusion.py # analysis
    diags/diag1000001 # output
)

# proton_boron_fusion_3d ######################################################
#
add_warpx_test(
    proton_boron_fusion_3d # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_3d_test_proton_boron_fusion_3d # inputs
    Examples/Tests/nuclear_fusion/analysis_proton_boron_fusion.py # analysis
    diags/diag1000001 # output
)
