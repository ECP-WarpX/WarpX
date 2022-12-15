.. _usage_domain_decomposition:

Domain Decomposition
====================

WarpX relies on a spatial domain decomposition for MPI parallelization. It provides two different ways for users to specify this decomposition, a `simple` way recommended for most users, and a `flexible` way recommended if more control is desired. The `flexible` method is required for dynamic load balancing to be useful.

1. Simple Method
----------------

The first and simplest method is to provide the ``warpx.numprocs = nx ny nz`` parameter, either at the command line or somewhere in your inputs deck. In this case, WarpX will split up the overall problem domain into exactly the specified number of subdomains, or `Boxes <https://amrex-codes.github.io/amrex/docs_html/Basics.html#box-intvect-and-indextype>`__ in the AMReX terminology, with the data defined on each `Box` having its own guard cells. The product of ``nx, ny, and nz`` must be exactly the desired number of MPI ranks. Note that, because there is exactly one `Box` per MPI rank when run this way, dynamic load balancing will not be possible, as there is no way of shifting `Boxes` around to achieve a more even load. This is the approach recommended for new users as it is the easiest to use.

.. note::

   If ``warpx.numprocs`` is *not* specified, WarpX will fall back on using the ``amr.max_grid_size`` and ``amr.blocking_factor`` parameters, described below.

2. More General Method
----------------------

The second way of specifying the domain decomposition provides greater flexibility and enables dynamic load balancing, but is not as easy to use. In this method, the user specifies inputs parameters ``amr.max_grid_size`` and ``amr.blocking_factor``, which can be thought of as the maximum and minimum allowed `Box` sizes. Now, the overall problem domain (specified by the ``amr.ncell`` input parameter) will be broken up into some number of `Boxes` with the specified characteristics. By default, WarpX will make the `Boxes` as big as possible given the constraints.

For example, if ``amr.ncell = 768 768 768``, ``amr.max_grid_size =  128``, and ``amr.blocking_factor = 32``, then AMReX will make 6 `Boxes` in each direction, for a total of 216 (the ``amr.blocking_factor`` does not factor in yet; however, see the section on mesh refinement below). If this problem is then run on 54 MPI ranks, there will be 4 boxes per rank initially. This problem could be run on as many as 216 ranks without performing any splitting.

When WarpX is run using this approach to domain decomposition, the number of MPI ranks does not need to be exactly equal to the number of ``Boxes``. Note also that if you run WarpX with more MPI ranks than there are boxes on the base level, WarpX will attempt to split the available ``Boxes`` until there is at least one for each rank to work on; this may cause it violate the constraints of ``amr.max_grid_size`` and ``amr.blocking_factor``.

.. note::

   The AMReX documentation on `Grid Creation <https://amrex-codes.github.io/amrex/docs_html/GridCreation.html#sec-grid-creation>`__ may also be helpful.

3. Performance Considerations
-----------------------------

In terms of performance, in general there is a trade off. Having many small boxes provides flexibility in terms of load balancing; however, the cost is increased time spent in communication due to surface-to-volume effects and increased kernel launch overhead when running on the GPUs. The ideal number of boxes per rank depends on how important dynamic load balancing is on your problem. If your problem is intrinsically well-balanced, like in a uniform plasma, then having a few, large boxes is best. But, if the problem is non-uniform and achieving a good load balance is critical for performance, having more, smaller `Boxes` can be worth it. In general, we find that running with something in the range of 4-8 `Boxes` per process is a good compromise for most problems.

You can also specify a separate `max_grid_size` and `blocking_factor` for each direction, using the parameters ``amr.max_grid_size_x``, ``amr.max_grid_size_y``, etc... . This allows you to request, for example, a "pencil" type domain decomposition that is long in one direction.


4. Mesh Refinement
------------------

With mesh refinement, the above picture is more complicated, as in general the number of boxes can not be predicted at the start of the simulation. The decomposition of the base level will proceed as outlined above. The refined region, however, will be covered by some number of Boxes whose sizes are consistent with ``amr.max_grid_size`` and ``amr.blocking_factor``. With mesh refinement, the blocking factor is important, as WarpX may decide to use `Boxes` smaller than ``amr.max_grid_size`` so as not to over-refine outside of the requested area. Note that you can specify a vector of values to make these parameters vary by level. For example, ``amr.max_grid_size = 128 64`` will make the max grid size be 128 on level 0 and 64 on level 1.

In general, the above performance considerations apply - varying these values such that there are 4-8 Boxes per rank on each level is a good guideline.
