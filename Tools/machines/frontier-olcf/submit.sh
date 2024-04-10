#!/usr/bin/env bash

#SBATCH -A <project id>
#SBATCH -J warpx
#SBATCH -o %x-%j.out
#SBATCH -t 00:10:00
#SBATCH -p batch
#SBATCH --ntasks-per-node=8
# Due to Frontier's Low-Noise Mode Layout only
# 7 instead of 8 cores are available per process
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#low-noise-mode-layout
#SBATCH --cpus-per-task=7
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest
#SBATCH -N 20
# In order to broadcast the WarpX executable and
# the linked libraries to each node we need to
# use the Non-Volatile Memory (NVMe) on node storage devices
# https://docs.olcf.ornl.gov/systems/frontier_user_guide.html#sbcasting-a-binary-with-libraries-stored-on-shared-file-systems
#SBATCH -C nvme

# In case you are using pywarpx, please uncomment the following line.
# This is needed to select the correct strategy to broadcast the
# warpx executable and the associated libraries to all the nodes.
#is_pywarpx=true

if [ ${is_pywarpx} ]; then

    # Replace this with the name of the WarpX input script
    input_script="inputs_3d.py"

    # Select WarpX dimensionality ("1", "2", "3", or "RZ")
    warpx_dim="3"

else
    # Replace this with the name of the WarpX executable
    exe="warpx"

    # Replace this with the name of the WarpX inputfile
    inputfile="inputs"
fi

# File to redirect stdandard output to
output="output.txt"

# Time and date when the job script started (system time)
date

########################################################################
########## INSTRUCTIONS BELOW SHOULD NORMALLY NOT BE MODIFIED ##########
########################################################################

# ___ Broadcast WarpX and associated libraries to compute nodes_________

if [ ${is_pywarpx} ]; then

    # Find out the path of the pywarpx module
    PYWARPX_PATH=$(python3 -c "import pywarpx; print(pywarpx.__path__[0])")

    if [ ${warpx_dim} == "1" ]; then
        export WARPX_SO="warpx_pybind_1d.cpython-39-x86_64-linux-gnu.so"
    elif [ ${warpx_dim} == "2" ]; then
        export WARPX_SO="warpx_pybind_2d.cpython-39-x86_64-linux-gnu.so"
    elif [ ${warpx_dim} == "3" ]; then
        export WARPX_SO="warpx_pybind_3d.cpython-39-x86_64-linux-gnu.so"
    elif [ ${warpx_dim} == "RZ" ]; then
        export WARPX_SO="warpx_pybind_rz.cpython-39-x86_64-linux-gnu.so"
    else
        echo "${warpx_dim} is not a valid WarpX dimensionality!"
        exit 1
    fi

    srun -N ${SLURM_NNODES} -n ${SLURM_NNODES} --ntasks-per-node=1 \
        bash -c "mkdir /mnt/bb/$USER/pywarpx/"

    # Broadcast all the files of the pywarpx module to the nodes
    for FILE in $PYWARPX_PATH/*
    do
        if [ ! -d "$FILE" ]; then
                # Broadcasting the libraries linked by WarpX.
            if [ $(basename $FILE) == $WARPX_SO ]; then
                sbcast --send-libs --exclude=NONE -pf $FILE /mnt/bb/$USER/pywarpx/$(basename $FILE)
                if [ ! "$?" == "0" ]; then
                    echo "SBCAST failed!"
                    exit 1
                fi
            else
                sbcast -pf $FILE /mnt/bb/$USER/pywarpx/$(basename $FILE)
                if [ ! "$?" == "0" ]; then
                    echo "SBCAST failed!"
                    exit 1
                fi
            fi
        fi
    done

    # All required libraries should be here now
    export LD_LIBRARY_PATH="/mnt/bb/$USER/pywarpx/${WARPX_SO}_libs"

    # libfabric dlopen's several libraries:
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$(pkg-config --variable=libdir libfabric)"
    # cray-mpich dlopen's libhsa-runtime64.so and
    # libamdhip64.so (non-versioned), so symlink on each node:
    srun -N ${SLURM_NNODES} -n ${SLURM_NNODES} --ntasks-per-node=1 --label -D /mnt/bb/$USER/pywarpx/${WARPX_SO}_libs \
        bash -c "if [ -f libhsa-runtime64.so.1 ]; then ln -s libhsa-runtime64.so.1 libhsa-runtime64.so; fi;
        if [ -f libamdhip64.so.5 ]; then ln -s libamdhip64.so.5 libamdhip64.so; fi"

    # Add broadcasted python libraries to PYTHONPATH
    export PYTHONPATH=/mnt/bb/$USER/:$PYTHONPATH

    echo "*************************************"
    echo "ldd /mnt/bb/$USER/pywarpx/${WARPX_SO} :"
    ldd /mnt/bb/$USER/pywarpx/${WARPX_SO}
    echo ""
    python3 -c "import pywarpx; print('pywarpx.__path__: ', pywarpx.__path__)"
    echo "*************************************"

else

    # From the documentation:
    # SBCAST executable from Orion to NVMe -- NOTE: ``-C nvme``
    # is needed in SBATCH headers to use the NVMe drive
    # NOTE: dlopen'd files will NOT be picked up by sbcast
    sbcast --send-libs --exclude=NONE -pf ${exe} /mnt/bb/$USER/warpx
    if [ "$?" != "0" ]; then
        # CHECK EXIT CODE. When SBCAST fails, it may leave partial
        # files on the compute nodes, and if you continue to launch srun,
        # your application may pick up partially complete shared library files,
        # which would give you confusing errors.
        echo "SBCAST failed!"
        exit 1
    fi
    # SBCAST sends all libraries detected by `ld` (minus any excluded),
    # and stores them in the same directory in each node's node-local storage
    # Any libraries opened by `dlopen` are NOT sent, since they are not
    # known by the linker at run-time.
    # All required libraries now reside in /mnt/bb/$USER/warpx_libs
    export LD_LIBRARY_PATH="/mnt/bb/$USER/warpx_libs"
    # libfabric dlopen's several libraries:
    export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:$(pkg-config --variable=libdir libfabric)"
    # cray-mpich dlopen's libhsa-runtime64.so and
    # libamdhip64.so (non-versioned), so symlink on each node:
    srun -N ${SLURM_NNODES} -n ${SLURM_NNODES} --ntasks-per-node=1 --label -D /mnt/bb/$USER/warpx_libs \
        bash -c "if [ -f libhsa-runtime64.so.1 ]; then ln -s libhsa-runtime64.so.1 libhsa-runtime64.so; fi;
        if [ -f libamdhip64.so.5 ]; then ln -s libamdhip64.so.5 libamdhip64.so; fi"
    # You may notice that some libraries are still linked
    # from /sw/frontier, even after SBCASTing.
    # This is because the Spack-build modules use RPATH to find their dependencies.
    # This behavior cannot be changed.
    echo "*****ldd /mnt/bb/$USER/warpx*****"
    ldd /mnt/bb/$USER/warpx
    echo "*************************************"

fi

# _____________________________________________________________________________


# ___ Set MPI and OMP environment variables____________________________________

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
# OLCFDEV-1597: OFI Poll Failed UNDELIVERABLE Errors
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

# _____________________________________________________________________________


# ___ Run WarpX simulation ____________________________________________________

if [ ${is_pywarpx} ]; then
    srun -N${SLURM_JOB_NUM_NODES} -n${TOTAL_NMPI} --ntasks-per-node=${WARPX_NMPI_PER_NODE} \
        ./${input_script}> ${output}
else
    srun -N${SLURM_JOB_NUM_NODES} -n${TOTAL_NMPI} --ntasks-per-node=${WARPX_NMPI_PER_NODE} \
        /mnt/bb/$USER/warpx ./${inputfile} > ${output}
fi
# _____________________________________________________________________________
