#!/usr/bin/env bash
#
# Copyright 2024 Luca Fedeli
#
# This file is part of WarpX.
#

# This script is a developer's tool to perform the
# checks done by the clang-tidy CI test locally.
#
# Note: this script is only tested on Linux

echo "============================================="
echo
echo "This script is a developer's tool to perform the"
echo "checks done by the clang-tidy CI test locally"
echo "_____________________________________________"

# Check source dir
REPO_DIR=$(cd $(dirname ${BASH_SOURCE})/../../ && pwd)
echo
echo "Your current source directory is: ${REPO_DIR}"
echo "_____________________________________________"

# Set number of jobs to use for compilation
PARALLEL="${WARPX_TOOLS_LINTER_PARALLEL:-4}"
echo
echo "${PARALLEL} jobs will be used for compilation."
echo "This can be overridden by setting the environment"
echo "variable WARPX_TOOLS_LINTER_PARALLEL, e.g.: "
echo
echo "$ export WARPX_TOOLS_LINTER_PARALLEL=8"
echo "$ ./Tools/Linter/runClangTidy.sh"
echo "_____________________________________________"

# Check clang version
export CC="${CLANG:-"clang"}"
export CXX="${CLANGXX:-"clang++"}"
export CTIDY="${CLANGTIDY:-"clang-tidy"}"
echo
echo "The following versions of the clang compiler and"
echo "of the clang-tidy linter will be used:"
echo
echo "clang version:"
which ${CC}
${CC} --version
echo
echo "clang++ version:"
which ${CXX}
${CXX} --version
echo
echo "clang-tidy version:"
which ${CTIDY}
${CTIDY} --version
echo
echo "This can be overridden by setting the environment"
echo "variables CLANG, CLANGXX, and CLANGTIDY e.g.: "
echo "$ export CLANG=clang-15"
echo "$ export CLANGXX=clang++-15"
echo "$ export CTIDCLANGTIDYY=clang-tidy-15"
echo "$ ./Tools/Linter/runClangTidy.sh"
echo
echo "*****************************************************"
echo "* Warning: clang v15 is currently used in CI tests. *"
echo "* It is therefore recommended to use this version.  *"
echo "* Otherwise, a newer version may find issues not    *"
echo "* currently covered by CI tests (checks are opt-in) *"
echo "* while older versions may not find all the issues. *"
echo "*****************************************************"
echo "_____________________________________________"

# Install BLAS++
echo
echo "Installing BLAS++"
echo
if [ -d ${REPO_DIR}/src/blaspp ]
then
  cd ${REPO_DIR}/src/blaspp
  git fetch --prune
  git checkout master
  git pull
  cd -
else
  git clone https://github.com/icl-utk-edu/blaspp.git ${REPO_DIR}/src/blaspp
fi
rm -rf ${REPO_DIR}/src/blaspp-build
cmake -S ${REPO_DIR}/src/blaspp -B ${REPO_DIR}/src/blaspp-build -Duse_openmp=ON -DCMAKE_CXX_STANDARD=17 -DCMAKE_INSTALL_PREFIX=${REPO_DIR}/blaspp-install
cmake --build ${REPO_DIR}/src/blaspp-build --target install --parallel ${PARALLEL}
rm -rf ${REPO_DIR}/src/blaspp-build
export CMAKE_PREFIX_PATH=${REPO_DIR}/blaspp-install:${CMAKE_PREFIX_PATH}
echo "_____________________________________________"

# Install LAPACK++
echo
echo "Installing LAPACK++"
echo
if [ -d ${REPO_DIR}/src/lapackpp ]
then
  cd ${REPO_DIR}/src/lapackpp
  git fetch --prune
  git checkout master
  git pull
  cd -
else
  git clone https://github.com/icl-utk-edu/lapackpp.git ${REPO_DIR}/src/lapackpp
fi
rm -rf ${REPO_DIR}/src/lapackpp-build
cmake -S ${REPO_DIR}/src/lapackpp -B ${REPO_DIR}/src/lapackpp-build -DCMAKE_CXX_STANDARD=17 -Dbuild_tests=OFF -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON -DCMAKE_INSTALL_PREFIX=${REPO_DIR}/lapackpp-install
cmake --build ${REPO_DIR}/src/lapackpp-build --target install --parallel ${PARALLEL}
rm -rf ${REPO_DIR}/src/lapackpp-build
export CMAKE_PREFIX_PATH=${REPO_DIR}/lapackpp-install:${CMAKE_PREFIX_PATH}
echo "_____________________________________________"

# Prepare clang-tidy wrapper
echo
echo "Prepare clang-tidy wrapper"
echo "The following wrapper ensures that only source files"
echo "in WarpX/Source/* are actually processed by clang-tidy"
echo
cat > ${REPO_DIR}/clang_tidy_wrapper << EOF
#!/bin/bash
REGEX="[a-z_A-Z0-9\/]*WarpX\/Source[a-z_A-Z0-9\/]+.cpp"
if [[ \$4 =~ \$REGEX ]];then
  ${CTIDY} \$@
fi
EOF
chmod +x ${REPO_DIR}/clang_tidy_wrapper
echo "clang_tidy_wrapper: "
cat ${REPO_DIR}/clang_tidy_wrapper
echo "_____________________________________________"

# Compile Warpx using clang-tidy
echo
echo "Compile Warpx using clang-tidy"
echo "The compilation requires libraries such as"
echo "fftw (for PSATD) and boost (for the QED module)"
echo

rm -rf ${REPO_DIR}/build_clang_tidy

cmake -S ${REPO_DIR} -B ${REPO_DIR}/build_clang_tidy \
  -DCMAKE_CXX_CLANG_TIDY="${REPO_DIR}/clang_tidy_wrapper;--system-headers=0;--config-file=${REPO_DIR}/.clang-tidy" \
  -DCMAKE_VERBOSE_MAKEFILE=ON  \
  -DWarpX_DIMS="1;2;3;RZ"      \
  -DWarpX_MPI=ON               \
  -DWarpX_COMPUTE=OMP          \
  -DWarpX_PSATD=ON             \
  -DWarpX_QED=ON               \
  -DWarpX_QED_TABLE_GEN=ON     \
  -DWarpX_OPENPMD=ON           \
  -DWarpX_PRECISION=SINGLE

cmake --build ${REPO_DIR}/build_clang_tidy -j ${PARALLEL} 2> ${REPO_DIR}/build_clang_tidy/clang-tidy.log

cat ${REPO_DIR}/build_clang_tidy/clang-tidy.log
echo
echo "============================================="
