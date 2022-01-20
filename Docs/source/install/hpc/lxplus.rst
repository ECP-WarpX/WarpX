.. _building-lxplus:

LXPLUS (CERN)
=============

The LXPLUS cluster is located at CERN.

* `Lxplus documentation <https://lxplusdoc.web.cern.ch>`__
* Batch system: `HTCondor <https://batchdocs.web.cern.ch/index.html>`__
* Filesystem locations:
    * User folder: ``/afs/cern.ch/user/<a>/<account>`` (10GByte)
    * Work folder: ``/afs/cern.ch/work/<a>/<account>`` (100GByte)
    * Eos storage: ``/eos/home-<a>/<account>`` (1T)

Through LXPLUS we have access to CPU and GPU nodes (the latter equipped with NVIDIA A100 or T4 GPUs).

Installation
------------
Only very little software is pre-installed on Lxplus so we show how to install from scratch all the dependencies using `Spack <https://spack.io>`__.

For size reasons it is not advisable to install WarpX in the ``$HOME`` directory, while it should be installed in the "work directory". For this purpose we set an environment variable with the path to the "work directory"

.. code-block:: bash

    export WORK=/afs/cern.ch/work/<a>/<account>/

We clone WarpX in ``$WORK``:

.. code-block:: bash

    cd $WORK
    git clone https://github.com/ECP-WarpX/WarpX.git

GCC
^^^
The pre-installed GNU compiler is outdated so we need a more recent compiler. Here we use the gcc 9.2.0 from the LCG project, but other options are possible.

We activate it by doing

.. code-block:: bash

    source /cvmfs/sft.cern.ch/lcg/releases/gcc/9.2.0/x86_64-centos7/setup.sh

In order to avoid using different compilers this line could be added directly into the ``$HOME/.bashrc`` file.

Spack
^^^^^
We download and activate Spack in ``$WORK``:

.. code-block:: bash

    cd $WORK
    git clone -c feature.manyFiles=true https://github.com/spack/spack.git
    source spack/share/spack/setup-env.sh

When installing packages Spack will try to set permissions in a way which is forbidden on the LXPLUS file system (AFS), but we can avoid this by setting the Spack option:

.. code-block:: bash

    spack config add "config:allow_sgid:false"

Now we add our gcc 9.2.0 compiler to spack:

.. code-block:: bash

    spack compiler find /cvmfs/sft.cern.ch/lcg/releases/gcc/9.2.0-afc57/x86_64-centos7/bin/

Installing the Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^

To install the depndencies we create an anonymous environment.
First of all we copy the envornmnet file in a local folder with the desired name of the environment (e.g. ``warpx-lxplus``):

.. code-block:: bash

    cd $WORK
    mkdir warpx-lxplus
    cp $WORK/WarpX/Tools/machines/lxplus-cern/spack.yaml $WORK/warpx-lxplus/.

Then we activate the environment:

.. code-block:: bash

    spack env activate -d warpx-lxplus

If we are planning on running with GPU support then we must set the environment variable ``SPACK_STACK_USE_CUDA``

.. code-block:: bash

    export SPACK_STACK_USE_CUDA=1

and if we want to use the python interface we must set the environment variable ``SPACK_STACK_USE_PYTHON``

.. code-block:: bash

    export SPACK_STACK_USE_PYTHON=1

Then we can install the required packages:

.. code-block:: bash

    spack install

The environment ``warpx-lxplus`` must be reactivated everytime that we log in so it could be a good idea to add the following lines to the ``.bashrc``:

.. code-block:: bash

    cd $WORK
    spack env activate -d warpx-lxplus
    cd $HOME    

Building WarpX
^^^^^^^^^^^^^^

We prepare and load the Spack software environment as above.
Then we build WarpX:

.. code-block:: bash

    cmake -S . -B build
    cmake --build build

Or if we need to compile with CUDA:

.. code-block:: bash

    cmake -S . -B build -DWarpX_COMPUTE=CUDA -DAMReX_CUDA_ARCH='7.0;7.5'
    cmake --build build

Python Bindings
^^^^^^^^^^^^^^^

Here we assume that a Python interpreter has been set up as explained previously.

Then we compile WarpX as in the previous section (with or without CUDA) adding ``-DWarpX_LIB=ON`` and then we install it into our Python:

.. code-block:: bash

    PYWARPX_LIB_DIR=$PWD/build/lib python3 -m pip wheel .
    python3 -m pip install pywarpx-*whl
