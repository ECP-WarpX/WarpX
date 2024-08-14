# Add tests (alphabetical order) ##############################################
#

# Deuterium_Deuterium_Fusion_3D ###############################################
#
add_warpx_test(
    Deuterium_Deuterium_Fusion_3D # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_3d_test_Deuterium_Deuterium_Fusion_3D # inputs
    Examples/Tests/nuclear_fusion/analysis_two_product_fusion.py # analysis
    diags/diag1000001 # output
)

# Deuterium_Deuterium_Fusion_3D_intraspecies ##################################
#
add_warpx_test(
    Deuterium_Deuterium_Fusion_3D_intraspecies # name
    3 # dims
    1 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_3d_test_Deuterium_Deuterium_Fusion_3D_intraspecies # inputs
    Examples/Tests/nuclear_fusion/analysis_deuterium_deuterium_3d_intraspecies.py # analysis
    diags/diag1000010 # output
)

# Deuterium_Tritium_Fusion_3D #################################################
#
add_warpx_test(
    Deuterium_Tritium_Fusion_3D # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_3d_test_Deuterium_Tritium_Fusion_3D # inputs
    Examples/Tests/nuclear_fusion/analysis_two_product_fusion.py # analysis
    diags/diag1000001 # output
)

# Deuterium_Tritium_Fusion_RZ #################################################
#
add_warpx_test(
    Deuterium_Tritium_Fusion_RZ # name
    RZ # dims
    2  # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_rz_test_Deuterium_Tritium_Fusion_RZ # inputs
    Examples/Tests/nuclear_fusion/analysis_two_product_fusion.py # analysis
    diags/diag1000001 # output
)

# Proton_Boron_Fusion_2D ######################################################
#
add_warpx_test(
    Proton_Boron_Fusion_2D # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_2d_test_Proton_Boron_Fusion_2D # inputs
    Examples/Tests/nuclear_fusion/analysis_proton_boron_fusion.py # analysis
    diags/diag1000001 # output
)

# Proton_Boron_Fusion_3D ######################################################
#
add_warpx_test(
    Proton_Boron_Fusion_3D # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/nuclear_fusion/inputs_3d_test_Proton_Boron_Fusion_3D # inputs
    Examples/Tests/nuclear_fusion/analysis_proton_boron_fusion.py # analysis
    diags/diag1000001 # output
)
