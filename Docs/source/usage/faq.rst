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


What kinds of RZ output do you support?
---------------------------------------

In RZ, supported detail of RZ output depends on the :ref:`output format <dataanalysis-formats>` that is configured in the :ref:`inputs file <running-cpp-parameters-diagnostics>`.

openPMD supports output of the detailed RZ modes and reconstructs representations on-the-fly in post-processing, e.g, in ``openPMD-viewer`` or other tools.
For some tools, this is in-development.

AMReX plotfiles and other in situ methods output a 2D reconstructed Cartesian slice at :math:`\theta=0` by default (and can opt-in to dump raw modes).
