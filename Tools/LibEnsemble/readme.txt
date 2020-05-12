#####################################################
# This readme presents results of a convergence scan
# with the input file in sim/inputs on Summit.
#####################################################

A convergence study, starting from an input file similar to ./sim/inputs (with a different set of physical parameters generated the APOSMM generator in a LibEnsemble run) is presented below.
The only parameter changed to go from one resolution to the next is amr.n_cell.

resolution   amr.n_cell   duration    f
res0          16  1024      10s       2.6340181298685122e-05
res1          32  2948      20s       5.074892524831451e-07 <- the one we've been running so far
res2          64  4096      54s       4.557730935993e-06
res3         128  8192     214s       5.683672880520962e-06
res4         256 16384    1092s       5.667607648747023e-06

I ran the last run on 2 GPUs instead of 1, and it took 780s instead of 1092s, meaning we can run it on either 1 GPU or 2.
That way, we can explore it either having all resolutions on the same number of GPUs, or on different number of GPUs.
In the end, we would probably be interested in the former, which can also be further investigating my going to higher resolution.

The simulations can be run with the same input file, with the following changes:
- The timeout in the batch script should be larger for higher resolutions (see numbers above)
- The jsrun command should be modified as, e.g., "jsrun ... warpx.exe inputs amr.n_cell=128 8192" for resolution res3
