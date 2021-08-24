#!/bin/bash -l
#SBATCH -C gpu
#SBATCH -t 01:00:00
#SBATCH -J AMReX_CNS
#SBATCH -o AMReX_CNS.o%j
#SBATCH -A <proj>
#SBATCH -N 4
#SBATCH -c 32
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=single:1

# ============
# -N =                 nodes
# -n =                 tasks (MPI ranks, usually = G)
# -G =                 GPUs (full Perlmutter node, 4)
# -c =                 CPU per task (128 total threads on CPU, 32 per GPU)
#
# --ntasks-per-node=   number of tasks (MPI ranks) per node (full node, 4)
# --gpus-per-task=     number of GPUs per task (MPI rank) (full node, 4)
# --gpus-per-node=     number of GPUs per node (full node, 4)
#
# --gpu-bind=single:1  sets only one GPU to be visible to each MPI rank
#                         (quiets AMReX init warnings)
#
# Recommend using --ntasks-per-node=4, --gpus-per-task=1 and --gpu-bind=single:1,
# as they are fixed values and allow for easy scaling with less adjustments.
#
# ============
# e.g.
# For one node:  -N 1, -n 4, -c 32, --gpus-per-node=4, --ntasks-per-node=4
# For two nodes: -N 2, -n 8, -c 32, --gpus-per-node=4, --ntasks-per-node=4

# salloc commands:
# ================
# Single node:
# salloc -N 1 --ntasks-per-node=4 -t 2:00:00 -C gpu -c 32 -G 4 -A <proj>
# Multi node:
# salloc -N 2 --ntasks-per-node=4 -t 2:00:00 -C gpu -c 32 -G 8 -A <proj>

EXE=./warpx
#EXE=../WarpX/build/bin/warpx.3d.MPI.CUDA.DP.OPMD.QED
#EXE=./main3d.gnu.TPROF.MPI.CUDA.ex
INPUTS=inputs_small

# Basic job submissions:
# =============================
# Run inside the current salloc session using available resources.
# Change parameters to match available resources & run with "./run.corigpu"
# srun -n 16 --ntasks-per-node=4 --gpus-per-task=1 --gpu-bind=single:1 ${EXE} ${INPUTS}

# Submit with the SBATCH configuration above to the gpu queue: "sbatch run.corigpu"
# Can also be ran with "./run.corigpu" to run with 1 CPU and 1 GPU.
srun ${EXE} ${INPUTS}





# ==============================
#  NEEDS TESTING AND ADJUSTMENT
# ==============================

# NSight Systems
# ==============

# @@ Simple Example:
#srun nsys profile -o nsys_out.%q{SLURM_PROCID}.%q{SLURM_JOBID} ${EXE} ${INPUTS}

# @@ Recommended Example:
#srun nsys profile -c nvtx -p "<TINY_PROFILER_NAME>@*" -e NSYS_NVTX_PROFILER_REGISTER_ONLY=0 -o nsys_out.%q{SLURM_PROCID}.%q{SLURM_JOBID} ${EXE} ${INPUTS}

# @@ Discussion:
#   This will run nsys profile and store performance data in a qdrep file named after '-o'
#   Open using nsight-sys $(pwd)/nsys_out.#.######.qdrep

#   To capture the NVTX ranges, included in TINY_PROFILE objects, use:
#       "-e NSYS_NVTX_PROFILER_REGISTER_ONLY=0"
#   (TINY_PROFILE's NVTX regions do not use registered strings at this time.)

#   Nsight systems creates a timeline over a single, contiguous block of time.
#   The start of the timeline can be selected using TINY_PROFILER's NVTX markers with:
#     -c nvtx -p "region_name@*"
#   This will turn on the profiling analysis at the first instance of the TINY_PROFILER region
#   and run to the end of the program. To stop the analysis at the end of the same region, add:
#     -x true
#   Note: This will only analyze the first instance of the region, so "-x true" should be used
#   for specific analyses, or on more inclusive timers, e.g. a timer around a full timestep.

# @@ Documentation:
#   For NSight System profiling flags:
#      https://docs.nvidia.com/nsight-systems/profiling/index.html#cli-profile-command-switch-options
#   For NSight examples to launch profiling, including region limiting:
#      https://docs.nvidia.com/nsight-systems/profiling/index.html#example-interactive-cli-command-sequences

# Running NSight Systems on multiple ranks
# ========================================

# Run Nsight Systems only profiling on $PROFILE_RANK rank on a multi-rank job
#    **** Preferred for most basic use cases
#srun ./profile_1rank.sh ${EXE} ${INPUTS}

# Uncomment and copy the following lines into profile_1rank.sh
# Adjust the nsys command line as needed for your test case.
# #!/bin/bash
# PROFILE_RANK=0
# if [ $SLURM_PROCID == $PROFILE_RANK ]; then
#   nsys profile -o nsys_out.%q{SLURM_PROCID}.%q{SLURM_JOBID} "$@"
# else
#   "$@"
# fi

# NSight Compute
# ==============

# Run Nsight Compute:
#    **** This will do a A LOT of analysis. Unless you want the entire job ran 7 times
#    **** with full profiling, limit the kernels profiled with additional flags:
#    For filtering examples, see:
#    https://docs.nvidia.com/nsight-compute/NsightComputeCli/index.html#nvtx-filtering
#    For full list of profile options, see:
#    https://docs.nvidia.com/nsight-compute/NsightComputeCli/index.html#command-line-options-profile
#    Recommended: limit kernels tested within a given BL_PROFILER timer with "--nvtx-include <configuration>"
#      Note: Must use TINY_PROFILE=TRUE and nvtx region names are equal to BL_PROFILER timer names.
#srun nv-nsight-cu-cli -o cucli_out.%q{SLURM_PROCID}.%q{SLURM_JOBID} ${EXE} ${INPUTS}
