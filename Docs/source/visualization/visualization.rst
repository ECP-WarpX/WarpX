Visualizing the simulation results
==================================

Output formats
--------------

WarpX can write data either in `plotfile` format (AMReX's native format), or
in `openPMD format <https://www.openpmd.org/>`_ (a common data schema on top of
HDF5 or ADIOS files for particle-in-cell codes).

.. note::

    This is controlled by the parameters ``warpx.plot_int`` (AMReX) or the
    ``warpx.openpmd_int`` & ``warpx.openpmd_backend`` options in the section
    :doc:`../running_cpp/parameters`.

Asynchronous IO
---------------

When using the AMReX `plotfile` format, users can set the ``amrex.async_out=1``
option to perform the IO in a non-blocking fashion, meaning that the simulation
will continue to run while an IO thread controls writing the data to disk.
This can significantly reduce the overall time spent in IO. This is primarily intended for
large runs on supercomputers such as Summit and Cori; depending on the MPI
implementation you are using, you may not see a benefit on your workstation.

When writing plotfiles, each rank will write to a separate file, up to some maximum number
(by default, 64). This maximum can be adjusted using the ``amrex.async_out_nfiles`` inputs
parameter. To use asynchronous IO with than ``amrex.async_out_nfiles`` MPI ranks, WarpX
WarpX must be compiled with the ``MPI_THREAD_MULTIPLE=TRUE`` flag.
Please see :doc:`../building/building` for building instructions.

This section describes some of the tools available to visualize the data:

.. toctree::
   :maxdepth: 1

   yt
   backtransformed_diags
   reduced_diags
   plot_distribution_mapping
   visit
   picviewer
   openpmdviewer
   advanced
   plot_parallel


In addition, WarpX also has In-Situ Visualization capabilities (i.e.
visualizing the data directly from the simulation, without dumping data
files to disk).

.. toctree::
    :maxdepth: 1

    sensei
    ascent

If you like the 3D rendering of laser wakefield acceleration
on the WarpX documentation front page (which is
also the avatar of the ECP-WarpX organization), you can find the serial
analysis script :download:`video_yt.py<../../../Tools/PostProcessing/video_yt.py>` as well
as a parallel analysis script
:download:`video_yt.py<../../../Tools/PostProcessing/yt3d_mpi.py>` used to make a similar
rendering for a beam-driven wakefield simulation, running parallel.

Warning: currently, quantities in the output file for iteration ``n`` are not all defined at the same physical time due to the staggering in time in WarpX.
The table below provides the physical time at which each quantity in the output file is written, in units of time step, for time step n.
"E part" means the electric field at each particle position.

======== ===== =====
quantity staggering
-------- -----------
         FDTD  PSATD
======== ===== =====
E        n     n
B        n     n
j        n-1/2 n-1/2
rho      n     n
position n     n
momentum n-1/2 n-1/2
E part   n     n
B part   n     n
======== ===== =====
