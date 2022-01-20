.. _debugging_warpx:

Debugging the code
==================

Sometimes, the code does not give you the result that you are expecting.
This can be due to a variety of reasons, from misunderstandings or changes in the :ref:`input parameters <running-cpp-parameters>`, system specific quirks, or bugs.
You might also want to debug your code as you implement new features in WarpX during development.

This section gives a step-by-step guidance on how to systematically check what might be going wrong.

Debugging Workflow
------------------

Try the following steps to debug a simulation:

#. Check the output text file, usually called ``output.txt``: are there warnings or errors present?
#. On an HPC system, look for the job output and error files, usually called ``WarpX.e...`` and ``WarpX.o...``.
   Read long messages from the top and follow potential guidance.
#. If your simulation already created output data files:
   Check if they look reasonable before the problem occurred; are the initial conditions of the simulation as you expected?
   Do you spot numerical artifacts or instabilities that could point to missing resolution or unexpected/incompatible numerical parameters?
#. Did the job output files indicate a crash? Check the ``Backtrace.<mpirank>`` files for the location of the code that triggered the crash.
   Backtraces are read from bottom (high-level) to top (most specific line that crashed).
#. In case of a crash, Backtraces can be more detailed if you :ref:`re-compile <install-developers>` with debug flags: for example, try compiling with ``-DCMAKE_BUILD_TYPE=RelWithDebInfo`` (some slowdown) or even ``-DCMAKE_BUILD_TYPE=Debug`` (this will make the simulation way slower) and rerun.
#. If debug builds are too costly, try instead compiling with ``-DAMReX_ASSERTIONS=ON`` to activate more checks and rerun.
#. If the problem looks like a memory violation, this could be from an invalid field or particle index access.
   Try compiling with ``-DAMReX_BOUND_CHECK=ON`` (this will make the simulation very slow), and rerun.
#. If the problem looks like a random memory might be used, try initializing memory with signaling Not-a-Number (NaN) values through the runtime option ``fab.init_snan = 1``.
   Further useful runtime options are ``amrex.fpe_trap_invalid``, ``amrex.fpe_trap_zero`` and ``amrex.fpe_trap_overflow`` (see details in the AMReX link below).
#. On Nvidia GPUs, if you suspect the problem might be a race condition due to a missing host / device synchronization, set the environment variable ``export CUDA_LAUNCH_BLOCKING=1`` and rerun.
#. Consider simplifying your input options and re-adding more options after having found a working baseline.

Fore more information, see also the `AMReX Debugging Manual <https://amrex-codes.github.io/amrex/docs_html/Basics.html#debugging>`__.

Last but not least: the community of WarpX developers and users can help if you get stuck.
Collect your above findings, describe where and what you are running and how you installed the code, describe the issue you are seeing with details and input files used and what you already tried.
Can you reproduce the problem with a smaller setup (less parallelism and/or less resolution)?
Report these details in a :ref:`WarpX GitHub issue <contact>`.

Debuggers
---------

See the `AMReX debugger section <https://amrex-codes.github.io/amrex/docs_html/Basics.html#breaking-into-debuggers>`__ on additional runtime parameters to

* disable backtraces
* rethrow exceptions
* avoid AMReX-level signal handling

You will need to set those runtime options to work directly with debuggers.
