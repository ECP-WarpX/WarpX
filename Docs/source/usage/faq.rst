.. _usage-faq:

FAQ
===

This section lists frequently asked usage questions.


What is MPI thread support level?
---------------------------------

We report this in output on startup together with other information.

That is the `MPI support for threaded execution <https://www.mpich.org/static/docs/v3.1/www3/MPI_Init_thread.html>`__, e.g., with OpenMP or system threads.

We currently only use this for optional, :ref:`async IO with AMReX plotfiles <running-cpp-parameters-diagnostics>`.
In the past, requesting MPI threading support had performance penalties, but we have not seen such anymore on recent systems.
Thus, we request it by default but you can overwrite it with a compile time option if it ever becomes needed.


How do I suppress tiny profiler output if I do not care to see it?
------------------------------------------------------------------

Via ``AMReX_TINY_PROFILE=OFF`` (see: :ref:`build options <building-cmake-options>` and then `AMReX build options <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html#customization-options>`__).
We change the default in ``cmake/dependencies/AMReX.cmake``.

Note that the tiny profiler adds literally no overhead to the simulation runtime, thus we enable it by default.


What do I need to know about using the boosted frame?
-----------------------------------------------------
The boost transformations are calculated internally to WarpX, so the input deck can be designed in the lab frame.  Little modification is needed beyond setting ``warpx.gamma_boost``.
There are new instabilities that arise (see :ref:`boosted frame theory <theory-boostedframe>`) that require different algorithms for stability; please see the :ref:`boosted examples <usage-examples>` for more help.
That said, there are some practical items worth considering that we list here to assist in designing and analyzing boosted frame simulations:

- The boosted frame simulation begins at boosted time :math:`t'=0`.
- Best practice is to set the upper :math:`z` limit of the problem domain to 0 with ``geometry.prob_hi``.
Otherwise it is possible that the boosted simulation begins with downstream objects unexpectedly in the problem domain before there is a chance to correctly calculate all interactions.


What kinds of RZ output do you support?
---------------------------------------

In RZ, supported detail of RZ output depends on the :ref:`output format <dataanalysis-formats>` that is configured in the :ref:`inputs file <running-cpp-parameters-diagnostics>`.

openPMD supports output of the detailed RZ modes and reconstructs representations on-the-fly in post-processing, e.g, in ``openPMD-viewer`` or other tools.
For some tools, this is in-development.

AMReX plotfiles and other in situ methods output a 2D reconstructed Cartesian slice at :math:`\theta=0` by default (and can opt-in to dump raw modes).
