# Add tests (alphabetical order) ##############################################
#

add_warpx_test(
    test_3d_hard_edged_quadrupoles # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/accelerator_lattice/inputs_test_3d_hard_edged_quadrupoles # inputs
    Examples/Tests/accelerator_lattice/analysis.py # analysis
    diags/diag1000050 # output
    OFF # depend on this test
)

add_warpx_test(
    test_3d_hard_edged_quadrupoles_boosted # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/accelerator_lattice/inputs_test_3d_hard_edged_quadrupoles_boosted # inputs
    Examples/Tests/accelerator_lattice/analysis.py # analysis
    diags/diag1000050 # output
    OFF # depend on this test
)

add_warpx_test(
    test_3d_hard_edged_quadrupoles_moving # name
    3 # dims
    2 # nprocs
    OFF # eb
    Examples/Tests/accelerator_lattice/inputs_test_3d_hard_edged_quadrupoles_moving # inputs
    Examples/Tests/accelerator_lattice/analysis.py # analysis
    diags/diag1000050 # output
    OFF # depend on this test
)
