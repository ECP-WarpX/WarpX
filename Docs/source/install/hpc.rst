.. _install-hpc:

HPC
===

On selected high-performance computing (HPC) systems, WarpX has documented or even pre-build installation routines.
Follow the guide here instead of the generic installation routines for optimal stability and best performance.

warpx.profile
-------------

Use a ``warpx.profile`` file to set up your software environment without colliding with other software.
Ideally, store that file directly in your ``$HOME/`` and source it after connecting to the machine:

.. code-block:: bash

   source $HOME/warpx.profile

We list example ``warpx.profile`` files below, which can be used to set up WarpX on various HPC systems.

HPC Systems
-----------

.. toctree::
   :maxdepth: 1

   hpc/adastra
   hpc/cori
   hpc/crusher
   hpc/frontier
   hpc/fugaku
   hpc/hpc3
   hpc/juwels
   hpc/lassen
   hpc/lawrencium
   hpc/lumi
   hpc/lxplus
   hpc/ookami
   hpc/perlmutter
   hpc/quartz
   hpc/summit
   hpc/spock
   hpc/taurus

.. tip::

   Your HPC system is not in the list?
   `Open an issue <https://github.com/ECP-WarpX/WarpX/issues>`__ and together we can document it!
