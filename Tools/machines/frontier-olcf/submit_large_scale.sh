#!/usr/bin/env bash

#SBATCH -A <project id>
#SBATCH -J warpx
#SBATCH -o %x-%j.out
#SBATCH -t 00:10:00
#SBATCH -p batch
#SBATCH --ntasks-per-node=8
# Due to Frontier's Low-Noise Mode Layout only 7 instead of 8 cores are available per process
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#low-noise-mode-layout
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH -N 20
# In order to broadcast the WarpX executable and the linked libraries to each node we need to
# use the Non-Volatile Memory (NVMe) on node storage devices
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#sbcasting-a-binary-with-libraries-stored-on-shared-file-systems
#SBATCH -C nvme

# Replace this with the name of the WarpX executable
exe="warpx"

# From the documentation:
# SBCAST executable from Orion to NVMe -- NOTE: ``-C nvme`` is needed in SBATCH headers to use the NVMe drive
# NOTE: dlopen'd files will NOT be picked up by sbcast
sbcast --send-libs --exclude=NONE -pf ${exe} /mnt/bb/$USER/${exe}
if [ ! "$?" == "0" ]; then
    # CHECK EXIT CODE. When SBCAST fails, it may leave partial files on the compute nodes, and if you continue to launch srun,
    # your application may pick up partially complete shared library files, which would give you confusing errors.
    echo "SBCAST failed!"
    exit 1
fi
# SBCAST sends all libraries detected by `ld` (minus any excluded), and stores them in the same directory in each node's node-local storage
# Any libraries opened by `dlopen` are NOT sent, since they are not known by the linker at run-time.
# All required libraries now reside in /mnt/bb/$USER/${exe}_libs
export LD_LIBRARY_PATH="/mnt/bb/$USER/${exe}_libs"
# libfabric dlopen's several libraries:
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$(pkg-config --variable=libdir libfabric)"
# cray-mpich dlopen's libhsa-runtime64.so and libamdhip64.so (non-versioned), so symlink on each node:
srun -N ${SLURM_NNODES} -n ${SLURM_NNODES} --ntasks-per-node=1 --label -D /mnt/bb/$USER/${exe}_libs \
    bash -c "if [ -f libhsa-runtime64.so.1 ]; then ln -s libhsa-runtime64.so.1 libhsa-runtime64.so; fi;
    if [ -f libamdhip64.so.5 ]; then ln -s libamdhip64.so.5 libamdhip64.so; fi"
# You may notice that some libraries are still linked from /sw/crusher, even after SBCASTing.
# This is because the Spack-build modules use RPATH to find their dependencies. This behavior cannot be changed.
echo "*****ldd /mnt/bb/$USER/${exe}*****"
ldd /mnt/bb/$USER/${exe}
echo "*************************************"

# load cray libs and ROCm libs
#export LD_LIBRARY_PATH=${CRAY_LD_LIBRARY_PATH}:${LD_LIBRARY_PATH}

# From the documentation:
# Each Frontier compute node consists of [1x] 64-core AMD EPYC 7A53
# "Optimized 3rd Gen EPYC" CPU (with 2 hardware threads per physical core) with
# access to 512 GB of DDR4 memory.
# Each node also contains [4x] AMD MI250X, each with 2 Graphics Compute Dies
# (GCDs) for a total of 8 GCDs per node. The programmer can think of the 8 GCDs
# as 8 separate GPUs, each having 64 GB of high-bandwidth memory (HBM2E).

# note (5-16-22 and 7-12-22)
# this environment setting is currently needed on Frontier to work-around a
# known issue with Libfabric (both in the May and June PE)
#export FI_MR_CACHE_MAX_COUNT=0  # libfabric disable caching
# or, less invasive:
export FI_MR_CACHE_MONITOR=memhooks  # alternative cache monitor

# Seen since August 2023
# OLCFDEV-1597: OFI Poll Failed UNDELIVERABLE Errors
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#olcfdev-1597-ofi-poll-failed-undeliverable-errors
export MPICH_SMP_SINGLE_COPY_MODE=NONE
export FI_CXI_RX_MATCH_MODE=software

# note (9-2-22, OLCFDEV-1079)
# this environment setting is needed to avoid that rocFFT writes a cache in
# the home directory, which does not scale.
export ROCFFT_RTC_CACHE_PATH=/dev/null

export OMP_NUM_THREADS=1
export WARPX_NMPI_PER_NODE=8
export TOTAL_NMPI=$(( ${SLURM_JOB_NUM_NODES} * ${WARPX_NMPI_PER_NODE} ))
srun -N${SLURM_JOB_NUM_NODES} -n${TOTAL_NMPI} --ntasks-per-node=${WARPX_NMPI_PER_NODE} \
    /mnt/bb/$USER/${exe} inputs > output.txt
