
#!/usr/bin/env bash
#
# Copyright 2022 Luca Fedeli
#
# License: BSD-3-Clause-LBNL

set -eu -o pipefail

export OMP_NUM_THREADS=2

#
# Generate QED lookup tables using external tool
#
./build_qed/bin/qed_tools/qed_table_generator \
    --table BW --mode DP --dndt_chi_min 0.01 \
    --dndt_chi_max 100 --dndt_how_many 64 \
    --pair_chi_min 0.01 --pair_chi_max 100 \
    --pair_chi_how_many 64 --pair_frac_how_many 64 \
    -o bw_table_tool
./build_qed/bin/qed_tools/qed_table_generator \
    --table QS --mode DP --dndt_chi_min 0.001 \
    --dndt_chi_max 100 --dndt_how_many 64 \
    --em_chi_min 0.001 --em_chi_max 100 \
    --em_chi_how_many 64 --em_frac_how_many 64 \
    --em_frac_min 1e-12 -o qs_table_tool

#
# Generate QED lookup tables using WarpX
#
./build_qed/bin/warpx \
    ./Examples/Tests/qed/quantum_synchrotron/inputs_3d \
    qed_bw.lookup_table_mode = "generate" \
    qed_bw.tab_dndt_chi_min = 0.01 \
    qed_bw.tab_dndt_chi_max = 100.0 \
    qed_bw.tab_dndt_how_many = 64 \
    qed_bw.tab_pair_chi_min = 0.01 \
    qed_bw.tab_pair_chi_max = 100.0 \
    qed_bw.tab_pair_chi_how_many = 64 \
    qed_bw.tab_pair_frac_how_many = 64 \
    qed_bw.save_table_in = "bw_table" \
    qed_qs.lookup_table_mode = "generate" \
    qed_qs.tab_dndt_chi_min = 0.001 \
    qed_qs.tab_dndt_chi_max = 100.0 \
    qed_qs.tab_dndt_how_many = 64 \
    qed_qs.tab_em_chi_min = 0.001 \
    qed_qs.tab_em_frac_min = 1.0e-12 \
    qed_qs.tab_em_chi_max = 100.0 \
    qed_qs.tab_em_chi_how_many = 64 \
    qed_qs.tab_em_frac_how_many = 64 \
    qed_qs.save_table_in = "qs_table"

#
# Convert lookup tables (generated with WarpX and with the external tool) in human-readable format
#
./build_qed/bin/qed_tools/qed_table_reader -i qs_table --table QS --mode DP -o qs_table
./build_qed/bin/qed_tools/qed_table_reader -i qs_table_tool --table QS --mode DP -o qs_table_tool
./build_qed/bin/qed_tools/qed_table_reader -i bw_table --table BW --mode DP -o bw_table
./build_qed/bin/qed_tools/qed_table_reader -i bw_table_tool --table BW --mode DP -o bw_table_tool

#
# Compare the generated lookup tables
#
diff bw_table_dndt bw_table_tool_dndt
diff bw_table_pair bw_table_tool_pair
diff qs_table_phot_em qs_table_tool_phot_em
diff qs_table_dndt qs_table_tool_dndt

#
# Run a WarpX simulation using the lookup tables generated by the external tool
#
./build_qed/bin/warpx.3d \
    ./Examples/Tests/qed/quantum_synchrotron/inputs_3d \
    qed_bw.lookup_table_mode = "load" \
    qed_bw.load_table_from = "bw_table_tool" \
    qed_qs.lookup_table_mode = "load" \
    qed_qs.load_table_from = "qs_table_tool"
