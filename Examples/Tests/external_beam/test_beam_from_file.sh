#! /usr/bin/env bash

rm -rf diags
rm -rf ORIGINAL_diags

mpirun -np 2 ~/warpx/build/bin/warpx inputs_3d

mv diags/diag1/openpmd_000000.h5 ./
mv diags ORIGINAL_diags

mpirun -np 2 ~/warpx/build/bin/warpx inputs_3d \
       driver.injection_style = external_file \
       driver.injection_file = "openpmd_%T.h5"

python compare.py
