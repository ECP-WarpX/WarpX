#!/bin/bash -l

# Copyright 2023 The WarpX Community
#
# Author: Axel Huebl
# License: BSD-3-Clause-LBNL

#SBATCH -t 00:10:00
#SBATCH -N 2
#SBATCH --job-name=WarpX
#SBATCH --account=<proj>
#SBATCH --qos=es_normal
# 2xV100 nodes
#SBATCH --partition=es1
#SBATCH --constraint=es1_v100
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH -o WarpX.o%j
#SBATCH -e WarpX.e%j
#S BATCH --mail-type=all
#S BATCH --mail-user=yourmail@lbl.gov

# executable & inputs file or python interpreter & PICMI script here
EXE=./warpx
INPUTS=inputs_3d

srun ${EXE} ${INPUTS} \
  > output_${SLURM_JOB_ID}.txt
