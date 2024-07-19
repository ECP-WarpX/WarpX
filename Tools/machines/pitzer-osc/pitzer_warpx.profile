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
