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


What design principles should I keep in mind when creating an input file?
-------------------------------------------------------------------------

Leave a cushion between lasers, particles, and the edge of computational domain.
The laser antenna and plasma species ``zmin`` can be less than or greater than  the ``geometry.prob_hi``,
but not exactly equal.


What do I need to know about using the boosted frame?
-----------------------------------------------------

The input deck can be designed in the lab frame and little modification to the physical set-up is needed --
most of the work is done internally.
Here are a few practical items to assist in designing boosted frame simulations:

- Ions must be explicitly included
- Best practice is to separate counter-propagating objects; things moving to the right should start with :math:`z <= 0` and things stationary or moving to the left (moving to the left in the boosted frame) should start with :math:`z > 0`
- Don't forget the general design principles listed above
- The boosted frame simulation begins at boosted time :math:`t'=0`
- Numerics and algorithms need to be adjusted, as there are numerical instabilities that arise in the boosted frame. For example, setting ``particles.use_fdtd_nci_corr=1`` for an FDTD simulation or setting ``psatd.use_default_v_galilean=1`` for a PSATD simulation. Be careful as this is overly simplistic and these options will not work in all cases.  Please see the :ref:`input parameters documentation <running-cpp-parameters>` and the :ref:`examples <usage-examples>` for more information

An in-depth discussion of the boosted frame is provided in the :ref:`moving window and optimal Lorentz boosted frame <theory-boostedframe>` section.


What kinds of RZ output do you support?
---------------------------------------

In RZ, supported detail of RZ output depends on the :ref:`output format <dataanalysis-formats>` that is configured in the :ref:`inputs file <running-cpp-parameters-diagnostics>`.

openPMD supports output of the detailed RZ modes and reconstructs representations on-the-fly in post-processing, e.g, in ``openPMD-viewer`` or other tools.
For some tools, this is in-development.

AMReX plotfiles and other in situ methods output a 2D reconstructed Cartesian slice at :math:`\theta=0` by default (and can opt-in to dump raw modes).
