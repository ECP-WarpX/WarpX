# Add tests (alphabetical order) ##############################################
#

# hard_edged_quadrupoles ######################################################
#
add_warpx_test(
    hard_edged_quadrupoles # name
    3 # dims
    2 # nprocs
    Examples/Tests/accelerator_lattice/inputs_3d_test_hard_edged_quadrupoles # inputs
    Examples/Tests/accelerator_lattice/analysis.py # analysis
    diags/diag1000050 # output
)

# hard_edged_quadrupoles_boosted ##############################################
#
add_warpx_test(
    hard_edged_quadrupoles_boosted # name
    3 # dims
    2 # nprocs
    Examples/Tests/accelerator_lattice/inputs_3d_test_hard_edged_quadrupoles_boosted # inputs
    Examples/Tests/accelerator_lattice/analysis.py # analysis
    diags/diag1000050 # output
)

# hard_edged_quadrupoles_moving ###############################################
#
add_warpx_test(
    hard_edged_quadrupoles_moving # name
    3 # dims
    2 # nprocs
    Examples/Tests/accelerator_lattice/inputs_3d_test_hard_edged_quadrupoles_moving # inputs
    Examples/Tests/accelerator_lattice/analysis.py # analysis
    diags/diag1000050 # output
)
