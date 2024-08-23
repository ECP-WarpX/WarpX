# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_2d_energy_conserving_thermal_plasma # name
    2 # dims
    2 # nprocs
    OFF # eb
    OFF # preparation script
    Examples/Tests/energy_conserving_thermal_plasma/inputs_test_2d_energy_conserving_thermal_plasma # inputs
    Examples/Tests/energy_conserving_thermal_plasma/analysis.py # analysis
    diags/diag1000500 # output
)
