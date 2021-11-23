.. _development-page:

==============================
Python package development use
==============================

This describes how to compile and install WarpX, pywarpx, and mewarpx in order
to make development changes require minimal recompilation.

System Setup - CCache version
-----------------------------

If building for a CUDA GPU, versions of CCache before 4.2 have very
weak support, and so provide minimal benefit. Thus, systems being used
for development should install a later release of CCache. This is most
easily done from its source code.

.. code-block:: bash

    mkdir ccache-tmp
    cd ccache-tmp
    git clone https://github.com/ccache/ccache/
    mkdir ccache-build
    cd ccache-build
    cmake ../ccache
    make -j
    sudo make install

The latter will install in ``/usr/local/bin``

Then, do clean CMake configurations to pick up the path the to newer
``ccache`` following the steps in the next section.

Setup
-----

.. code-block:: bash

    # Create a dedicated build directory, outside the source tree
    mkdir -p ~/build/build_tree_name
    cd  ~/build/build_tree_name

    # Setting this environment variable indicates to pywarpx that it should refer
    # to an already built WarpX library rather than compiling it for itself
    export PYWARPX_LIB_DIR=~/build/build_tree_name/lib/

    # Set WarpX_COMPUTE as desired
    # Set WarpX_LIBS=ON to generate the shared library that pywarpx will reference
    # Pass -G Ninja to use the Ninja build system since it's somewhat faster than GNU Make
    cmake ~/repos/WarpX/ -DWarpX_DIMS=2 -DWarpX_EB=ON -DWarpX_OPENPMD=ON \
        -DWarpX_COMPUTE=CUDA \
        -DWarpX_LIB=ON \
        -G Ninja

    # Run ninja to build
    ninja


    # Create a dedicated Python virtual environment
    mkdir -p ~/venvs/venv_tree_name
    cd  ~/venvs/venv_tree_name

    python3 -m venv .
    source bin/activate

    # Install all of the pre-requisite Python packages
    bash ~/repos/WarpX/mewarpx/medocker/docker_shared/python_packages.sh

    # Install pywarpx
    pip install -e ~/repos/WarpX/Python/

    # Install mewarpx with full set of support libraries
    pip install -e ~/repos/WarpX/mewarpx/[complete]


Modification of C++ code
------------------------

.. code-block:: bash

    # ... edit source
    ninja -C ~/build/build_tree_name


Modification of Python code
---------------------------

No manual steps needed after editing - just run the code!
