.. _developers-local-compile:

Fast, Local Compilation
=======================

For simplicity, WarpX :ref:`compilation with CMake <building-cmake>` by default downloads, configures and compiles compatible versions of :ref:`central dependencies <install-dependencies>` such as:

* `AMReX <https://amrex-codes.github.io>`__
* `PICSAR <https://github.com/ECP-WarpX/picsar>`__
* `openPMD-api <https://github.com/openPMD/openPMD-api>`__
* `pyAMReX <https://github.com/AMReX-Codes/pyamrex>`__
* `pybind11 <https://github.com/pybind/pybind11>`__

on-the-fly, which is called a *superbuild*.

In some scenarios, e.g., when compiling without internet, with slow internet access, or when working on WarpX and its dependencies, other strategies might be preferable.
In the below workflows, you as the developer need to make sure to use compatible versions of the dependencies you provide.


.. _developers-local-compile-src:

Compiling From Local Sources
----------------------------

This workflow is best for developers that make changes to WarpX, AMReX, PICSAR, openPMD-api and/or pyAMReX at the same time.
For instance, use this if you add a feature in AMReX and want to try it in WarpX before it is proposed as a pull request for inclusion in AMReX.

Instead of downloading the source code of the above dependencies, one can also use an already cloned source copy.
For instance, clone these dependencies to ``$HOME/src``:

.. code-block:: bash

   cd $HOME/src

   git clone https://github.com/ECP-WarpX/WarpX.git warpx
   git clone https://github.com/AMReX-Codes/amrex.git
   git clone https://github.com/openPMD/openPMD-api.git
   git clone https://github.com/ECP-WarpX/picsar.git
   git clone https://github.com/AMReX-Codes/pyamrex.git
   git clone https://github.com/pybind/pybind11.git

Now modify the dependencies as needed in their source locations, update sources if you cloned them earlier, etc.
When building WarpX, :ref:`the following CMake flags <building-cmake-options>` will use the respective local sources:

.. code-block:: bash

   cd src/warpx

   rm -rf build

   cmake -S . -B build  \
     -DWarpX_PYTHON=ON  \
     -DWarpX_amrex_src=$HOME/src/amrex          \
     -DWarpX_openpmd_src=$HOME/src/openPMD-api  \
     -DWarpX_picsar_src=$HOME/src/picsar        \
     -DWarpX_pyamrex_src=$HOME/src/pyamrex      \
     -DWarpX_pybind11_src=$HOME/src/pybind11

   cmake --build build -j 8


.. _developers-local-compile-findpackage:

Compiling With Pre-Compiled Dependencies
----------------------------------------

This workflow is the best and fastest to compile WarpX, when you just want to change code in WarpX and have the above central dependencies already made available *in the right configurations* (e.g., w/ or w/o MPI or GPU support) from a :ref:`module system <install-hpc>` or :ref:`package manager <install-dependencies>`.

Instead of downloading the source code of the above central dependencies, or using a local copy of their source, we can compile and install those once and instruct CMake to `find their install locations and configurations <https://hsf-training.github.io/hsf-training-cmake-webpage/09-findingpackages/index.html>`__.

WarpX supports this with :ref:`the following CMake flags <building-cmake-options>`:

.. code-block:: bash

   cd src/warpx

   rm -rf build

   cmake -S . -B build  \
     -DWarpX_PYTHON=ON  \
     -DWarpX_amrex_internal=OFF    \
     -DWarpX_openpmd_internal=OFF  \
     -DWarpX_picsar_internal=OFF   \
     -DWarpX_pyamrex_internal=OFF  \
     -DWarpX_pybind11_internal=OFF

   cmake --build build -j 8

As a background, this is also the workflow how WarpX is built in :ref:`package managers such as Spack and Conda-Forge <install-dependencies>`.
