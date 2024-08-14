# Add tests (alphabetical order) ##############################################
#

# energy_conserving_thermal_plasma ############################################
#
add_warpx_test(
    energy_conserving_thermal_plasma # name
    2 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/energy_conserving_thermal_plasma/inputs_2d_test_energy_conserving_thermal_plasma # inputs
    Examples/Tests/energy_conserving_thermal_plasma/analysis.py # analysis
    diags/diag1000500 # output
)
