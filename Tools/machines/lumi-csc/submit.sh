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

# Seen since August 2023 seen on OLCF (not yet seen on LUMI?)
# OLCFDEV-1597: OFI Poll Failed UNDELIVERABLE Errors
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1597-ofi-poll-failed-undeliverable-errors
#export MPICH_SMP_SINGLE_COPY_MODE=NONE
#export FI_CXI_RX_MATCH_MODE=software

# note (9-2-22, OLCFDEV-1079)
# this environment setting is needed to avoid that rocFFT writes a cache in
# the home directory, which does not scale.
export ROCFFT_RTC_CACHE_PATH=/dev/null

export OMP_NUM_THREADS=1

# LUMI documentation suggests using the following wrapper script
# to set the ROCR_VISIBLE_DEVICES to the value of SLURM_LOCALID
# see https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
cat << EOF > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=\$SLURM_LOCALID
exec \$*
EOF

chmod +x ./select_gpu

sleep 1

# LUMI documentation suggests using the following CPU bind
# so that the node local rank and GPU ID match
# see https://docs.lumi-supercomputer.eu/runjobs/scheduled-jobs/lumig-job/
CPU_BIND="map_cpu:48,56,16,24,1,8,32,40"

srun --cpu-bind=${CPU_BIND} ./select_gpu ./warpx inputs | tee outputs.txt
rm -rf ./select_gpu
