# Add tests (alphabetical order) ##############################################
#

# Langmuir_fluid_1d ###########################################################
#
add_warpx_test(
    Langmuir_fluid_1d # name
    1  # dims
    ON # mpi
    Examples/Tests/langmuir_fluid/inputs_1d_test_Langmuir_fluid_1d # inputs
    Examples/Tests/langmuir_fluid/analysis_1d.py # analysis
    diags/diag1000080 # output
)

# Langmuir_fluid_2d ###########################################################
#
add_warpx_test(
    Langmuir_fluid_2d # name
    2  # dims
    ON # mpi
    Examples/Tests/langmuir_fluid/inputs_2d_test_Langmuir_fluid_2d # inputs
    Examples/Tests/langmuir_fluid/analysis_2d.py # analysis
    diags/diag1000080 # output
)

# Langmuir_fluid_3d ###########################################################
#
add_warpx_test(
    Langmuir_fluid_3d # name
    3  # dims
    ON # mpi
    Examples/Tests/langmuir_fluid/inputs_3d_test_Langmuir_fluid_3d # inputs
    Examples/Tests/langmuir_fluid/analysis_3d.py # analysis
    diags/diag1000040 # output
)

# Langmuir_fluid_rz ###########################################################
#
add_warpx_test(
    Langmuir_fluid_rz # name
    RZ # dims
    ON # mpi
    Examples/Tests/langmuir_fluid/inputs_rz_test_Langmuir_fluid_rz # inputs
    Examples/Tests/langmuir_fluid/analysis_rz.py # analysis
    diags/diag1000080 # output
)
