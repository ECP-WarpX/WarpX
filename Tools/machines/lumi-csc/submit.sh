#!/bin/bash -l

#SBATCH -A <project id>
#SBATCH --job-name=warpx
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=standard-g
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --gpus-per-node=8
#SBATCH --time=00:10:00

export MPICH_GPU_SUPPORT_ENABLED=1

# note (12-12-22)
# this environment setting is currently needed on LUMI to work-around a
# known issue with Libfabric
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

# note (9-2-22, OLCFDEV-1079)
# this environment setting is needed to avoid that rocFFT writes a cache in
# the home directory, which does not scale.
export ROCFFT_RTC_CACHE_PATH=/dev/null

export OMP_NUM_THREADS=1

cat << EOF > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

CPU_BIND="map_cpu:48,56,16,24,1,8,32,40"

srun --cpu-bind=${CPU_BIND} ./select_gpu ./warpx inputs
rm -rf ./select_gpu
