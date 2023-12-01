# Compiler settings
export CC=$(which clang)
export CXX=$(which clang++)
export CUDACXX=$(which nvcc)
export CUDAHOSTCXX=$(which clang++)

# cmake -S . -B build -G Ninja \
#     -DCMAKE_WARN_DEPRECATED=OFF \
#     -DWarpX_MPI=OFF \
#     -DWarpX_LIB=ON \
#     -DWarpX_APP=ON \
#     -DWarpX_COMPUTE=CUDA \
#     -DWarpX_DIMS="1;2;RZ;3" \
#     -DWarpX_PYTHON=ON \
#     -DPYINSTALLOPTIONS="--user" \
#     -DWarpX_amrex_src="$HOME/src/amrex" \
#     -DWarpX_amrex_internal=OFF  \
#     -DWarpX_pyamrex_internal=OFF \
#     -DWarpX_pyamrex_src="$HOME/src/pyamrex" \
#     -DCMAKE_INSTALL_PREFIX="$HOME/.local"   \
#     -DCMAKE_PREFIX_PATH="$HOME/src/pyamrex"

cmake --build build --target pip_install -j 8
