.. _parallelization_warpx:

Parallelization in WarpX
========================

When running a simulation, the domain is split into independent
rectangular sub-domains (called **grids**). This is the way AMReX, a core
component of WarpX, handles parallelization and/or mesh refinement. Furthermore,
this decomposition makes load balancing possible: each MPI rank typically computes
a few grids, and a rank with a lot of work can transfer one or several **grids**
to their neighbors.

A user
does not specify this decomposition explicitly. Instead, the user gives hints to
the code, and the actual decomposition is determined at runtime, depending on
the parallelization. The main user-defined parameters are
``amr.max_grid_size`` and ``amr.blocking_factor``.

AMReX ``max_grid_size`` and ``blocking_factor``
-----------------------------------------------

* ``amr.max_grid_size`` is the maximum number of cells per **grid** along each
  direction (default ``amr.max_grid_size=32`` in 3D).

* ``amr.blocking_factor``: is the minimum number of cells per **grid** along each
  direction (default ``amr.blocking_factor=8``).
  Note that both the domain (at each level) and ``max_grid_size`` must be divisible by ``blocking_factor``.

  .. note::

     You can use the parameters above if you want the same number of cells in all directions.
     Or you can set ``amr.max_grid_size_x``, ``amr.max_grid_size_y`` and ``amr.max_grid_size_z``;
     ``amr.blocking_factor_x``, ``amr.blocking_factor_x`` and ``amr.blocking_factor_x`` to different numbers of cells.

The total number of **grids** is determined using those two restrictions and the number of
ranks used to run the simulation. You can visit `AMReX <https://amrex-codes.github.io/amrex/docs_html/GridCreation.html?highlight=blocking_factor>`_
documentation for more information on the two parameters.

These parameters can have a dramatic impact on the code performance. Each
**grid** in the decomposition is surrounded by guard cells, thus increasing the
amount of data, computation and communication. Hence having a too small
``max_grid_size``, may ruin the code performance.

On the other hand, a too-large ``max_grid_size`` is likely to result in a single
grid per MPI rank, thus preventing load balancing. By setting these two
parameters, the user wants to give some flexibility to the code while avoiding
pathological behaviors.

For more information on this decomposition, see the
`Gridding and Load Balancing <https://amrex-codes.github.io/amrex/docs_html/ManagingGridHierarchy_Chapter.html>`__
page on AMReX documentation.

For specific information on the dynamic load balancer used in WarpX, visit the
`Load Balancing <https://amrex-codes.github.io/amrex/docs_html/LoadBalancing.html>`__
page on the AMReX documentation.

The best values for these parameters strongly depends on a number of parameters,
among which numerical parameters:

* Algorithms used (Maxwell/spectral field solver, filters, order of the
  particle shape factor)

* Number of guard cells (that depends on the particle shape factor and
  the type and order of the Maxwell solver, the filters used, `etc.`)

* Number of particles per cell, and the number of species

and MPI decomposition and computer architecture used for the run:

* GPU or CPU

* Number of OpenMP threads

* Amount of high-bandwidth memory.

Because these parameters put additional constraints on the domain size for a
simulation, it can be cumbersome to calculate the number of cells and the
physical size of the computational domain for a given resolution. This
:download:`Python script <../../../../Tools/DevUtils/compute_domain.py>` does it
automatically.

When using the RZ spectral solver, the values of ``amr.max_grid_size`` and ``amr.blocking_factor`` are constrained since the solver
requires that the full radial extent be within a each block.
For the radial values, any input is ignored and the max grid size and blocking factor are both set equal to the number of radial cells.
For the longitudinal values, the blocking factor has a minimum size of 8, allowing the computational domain of each block to be large enough relative to the guard cells for reasonable performance, but the max grid size and blocking factor must also be small enough so that there will be at least one block per processor.
If max grid size and/or blocking factor are too large, they will be silently reduced as needed.
If there are too many processors so that there is not enough blocks for the number processors, WarpX will abort.
