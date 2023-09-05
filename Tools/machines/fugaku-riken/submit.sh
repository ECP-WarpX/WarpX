#!/bin/bash
#PJM -L "node=48"
#PJM -L "rscgrp=small"
#PJM -L "elapse=0:30:00"
#PJM -s
#PJM -L "freq=2200,eco_state=2"
#PJM --mpi "max-proc-per-node=12"
#PJM -x PJM_LLIO_GFSCACHE=/vol0004:/vol0003
#PJM --llio localtmp-size=10Gi
#PJM --llio sharedtmp-size=10Gi

export NODES=48
export MPI_RANKS=$((NODES * 12))
export OMP_NUM_THREADS=4

export EXE="./warpx"
export INPUT="i.3d"

export XOS_MMM_L_PAGING_POLICY=demand:demand:demand

llio_transfer ${EXE}

mpiexec -stdout-proc ./output.%j/%/1000r/stdout -stderr-proc ./output.%j/%/1000r/stderr -n ${MPI_RANKS} ${EXE} ${INPUT}

llio_transfer --purge ${EXE}
