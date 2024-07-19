# please set your project account
export proj=""  # change me!

# remembers the location of this script
export MY_PROFILE=$(cd $(dirname $BASH_SOURCE) && pwd)"/"$(basename $BASH_SOURCE)
if [ -z ${proj-} ]; then echo "WARNING: The 'proj' variable is not yet set in your $MY_PROFILE file! Please edit its line 2 to continue!"; return; fi

module purge
# load these modules for optimization (todo: make this compile)
# module load cmake/3.25.2
# module load intel/19.0.5
# module load cuda/11.6.1
# module load mvapich2-gdr/2.3.7
# export CC=$(which icc)
# export CXX=$(which icpc)
# export FC=$(which ifort)

# load these modules for stability
module load cmake/3.25.2
module load gnu/12.3.0
module load cuda/12.3.0
module load openmpi/4.1.5-hpcx
export CC=$(which gcc)
export CXX=$(which g++)
export FC=$(which gfortran)

export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=${CXX}
export SW_DIR="${HOME}/sw/osc/pitzer/v100" # add subdirectories if you want to load a different set of modules

# optional: for building python binding
module load miniconda3/24.1.2-py310

if [ -d "${SW_DIR}/venvs/warpx-pitzer-v100" ]
then
  source ${SW_DIR}/venvs/warpx-pitzer-v100/bin/activate
fi

# Uncomment this section to set the CPU version of WarpX Python bindings as the default
# if [ -d "${SW_DIR}/venvs/warpx-pitzer" ]
# then
#   source ${SW_DIR}/venvs/warpx-pitzer/bin/activate
# fi

# an alias to request an interactive batch node for one hour
#   for parallel execution, start on the batch node: srun <command>
alias getNode="salloc -N 1 --ntasks-per-node=2 --cpus-per-task=20 --gpus-per-task=v100:1 -t 1:00:00 -A $proj"
# an alias to run a command on a batch node for up to 30min
#   usage: runNode <command>
alias runNode="srun -N 1 --ntasks-per-node=2 --cpus-per-task=20 --gpus-per-task=v100:1 -t 1:00:00 -A $proj"

# optional: for PSATD in RZ geometry support
export CMAKE_PREFIX_PATH=${SW_DIR}/blaspp-2024.05.31:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/lapackpp-2024.05.31:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/blaspp-2024.05.31/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/lapackpp-2024.05.31/lib64:$LD_LIBRARY_PATH

# optional: for QED lookup table generation support (note: it will have compatibility issue with gnu compiler)
# module load boost/1.75.0

# optional: for openPMD support (hdf5 and adios2)
module load hdf5/1.12.2
# module load hdf5/1.12.0)
export CMAKE_PREFIX_PATH=${SW_DIR}/c-blosc-1.21.6:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=${SW_DIR}/adios2-2.10.1:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=${SW_DIR}/c-blosc-1.21.6/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=${SW_DIR}/adios2-2.10.1/lib64:$LD_LIBRARY_PATH
export PATH=${SW_DIR}/adios2-2.10.1/bin:${PATH}
