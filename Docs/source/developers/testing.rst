.. _developers-testing:

Testing the code
================

When adding a new feature, you want to make sure that (i) you did not break the existing code and (ii) your contribution gives correct results. While existing capabilities are tested regularly remotely (when commits are pushed to an open PR on CI, and every night on local clusters), it can also be useful to run tests on your custom input file. This section details how to use both automated and custom tests.

Continuous Integration in WarpX
-------------------------------

Configuration
^^^^^^^^^^^^^

Our regression tests are using the suite published and documented at `AMReX-Codes/regression_testing <https://github.com/AMReX-Codes/regression_testing>`__.

Most of the configuration of our regression tests happens in ``Regression/Warpx-tests.ini``.
We slightly modify this file in ``Regression/prepare_file_ci.py``.

For example, if you like to change the compiler to compilation to build on Nvidia GPUs, modify this block to add ``-DWarpX_COMPUTE=CUDA``:

.. code-block:: toml

   [source]
   dir = /home/regtester/AMReX_RegTesting/warpx
   branch = development
   cmakeSetupOpts = -DAMReX_ASSERTIONS=ON -DAMReX_TESTING=ON -DWarpX_COMPUTE=CUDA

We also support changing compilation options via the usual :ref:`build enviroment variables <building-cmake-envvars>`.
For instance, compiling with ``clang++ -Werror`` would be:

.. code-block:: sh

   export CXX=$(which clang++)
   export CXXFLAGS="-Werror"


Run the test suite locally
--------------------------

Once your new feature is ready, there are ways to check that you did not break anything.
WarpX has automated tests running every time a commit is added to an open pull request.
The list of automated tests is defined in `./Regression/WarpX-tests.ini <https://github.com/ECP-WarpX/WarpX/blob/development/Regression/WarpX-tests.ini>`__.

For easier debugging, it can be convenient to run the tests on your local machine by executing the script
`./run_test.sh <https://github.com/ECP-WarpX/WarpX/blob/development/run_test.sh>`__ from WarpX's root folder, as illustrated in the examples below:

.. code-block:: sh

    # Example:
    # run all tests defined in ./Regression/WarpX-tests.ini
    ./run_test.sh

    # Example:
    # run only the test named 'pml_x_yee'
    ./run_test.sh pml_x_yee

    # Example:
    # run only the tests named 'pml_x_yee', 'pml_x_ckc' and 'pml_x_psatd'
    ./run_test.sh pml_x_yee pml_x_ckc pml_x_psatd

Note that the script `./run_test.sh <https://github.com/ECP-WarpX/WarpX/blob/development/run_test.sh>`__ runs the tests with the exact same compile-time options and runtime options used to run the tests remotely.

Moreover, the script `./run_test.sh <https://github.com/ECP-WarpX/WarpX/blob/development/run_test.sh>`__ compiles all the executables that are necessary in order to run the chosen tests.
The default number of threads allotted for compiling is set with ``numMakeJobs = 8`` in `./Regression/WarpX-tests.ini <https://github.com/ECP-WarpX/WarpX/blob/ad74bcbdd131a8797339ba38370b1195d0aecffb/Regression/WarpX-tests.ini#L20>`__.
However, when running the tests on a local machine, it is usually possible and convenient to allot more threads for compiling, in order to speed up the builds.
This can be accomplished by setting the environment variable ``WARPX_CI_NUM_MAKE_JOBS``, with the preferred number of threads that fits your local machine, e.g. ``export WARPX_CI_NUM_MAKE_JOBS=16`` (or less if your machine is smaller).
On public CI, we overwrite the value to ``WARPX_CI_NUM_MAKE_JOBS=2``, in order to avoid overloading the available remote resources.
Note that this will not change the number of threads used to run each test, but only the number of threads used to compile each executable necessary to run the tests.

Once the execution of `./run_test.sh <https://github.com/ECP-WarpX/WarpX/blob/development/run_test.sh>`__ is completed, you can find all the relevant files associated with each test in one single directory.
For example, if you run the single test ``pml_x_yee``, as shown above, on 04/30/2021,  you can find all relevant files in the directory ``./test_dir/rt-WarpX/WarpX-tests/2021-04-30/pml_x_yee/``.
The content of this directory will look like the following (possibly including backtraces if the test crashed at runtime):

.. code-block:: sh

    $ ls ./test_dir/rt-WarpX/WarpX-tests/2021-04-30/pml_x_yee/
    analysis_pml_yee.py     # Python analysis script
    inputs_2d               # input file
    main2d.gnu.TEST.TPROF.MTMPI.OMP.QED.GPUCLOCK.ex  # executable
    pml_x_yee.analysis.out  # Python analysis output
    pml_x_yee.err.out       # error output
    pml_x_yee.make.out      # build output
    pml_x_yee_plt00000/     # data output (initialization)
    pml_x_yee_plt00300/     # data output (last time step)
    pml_x_yee.run.out       # test output


Add a test to the suite
-----------------------

There are three steps to follow to add a new automated test (illustrated here for PML boundary conditions):

* An input file for your test, in folder `Example/Tests/...`. For the PML test, the input file is at ``Examples/Tests/PML/inputs_2d``. You can also re-use an existing input file (even better!) and pass specific parameters at runtime (see below).
* A Python script that reads simulation output and tests correctness versus theory or calibrated results. For the PML test, see ``Examples/Tests/PML/analysis_pml_yee.py``. It typically ends with Python statement ``assert( error<0.01 )``.
* If you need a new Python package dependency for testing, add it in ``Regression/requirements.txt``
* Add an entry to ``Regression/WarpX-tests.ini``, so that a WarpX simulation runs your test in the continuous integration process, and the Python script is executed to assess the correctness. For the PML test, the entry is

.. code-block::

   [pml_x_yee]
   buildDir = .
   inputFile = Examples/Tests/PML/inputs2d
   runtime_params = warpx.do_dynamic_scheduling=0 algo.maxwell_solver=yee
   dim = 2
   addToCompileString =
   cmakeSetupOpts = -DWarpX_DIMS=2
   restartTest = 0
   useMPI = 1
   numprocs = 2
   useOMP = 1
   numthreads = 1
   compileTest = 0
   doVis = 0
   analysisRoutine = Examples/Tests/PML/analysis_pml_yee.py

If you re-use an existing input file, you can add arguments to ``runtime_params``, like ``runtime_params = amr.max_level=1 amr.n_cell=32 512 max_step=100 plasma_e.zmin=-200.e-6``.

.. note::

   If you added ``analysisRoutine = Examples/analysis_default_regression.py``, then run the new test case locally and add the :ref:`checksum <developers-checksum>` file for the expected output.

.. note::

   We run those tests on our continuous integration services, which at the moment only have 2 virtual CPU cores.
   Thus, make sure that the product of ``numprocs`` and ``numthreads`` for a test is ``<=2``.


Useful tool for plotfile comparison: ``fcompare``
-------------------------------------------------

AMReX provides ``fcompare``, an executable that takes two ``plotfiles`` as input and returns the absolute and relative difference for each field between these two plotfiles. For some changes in the code, it is very convenient to run the same input file with an old and your current version, and ``fcompare`` the plotfiles at the same iteration. To use it:

.. code-block:: sh

   # Compile the executable
   cd <path to AMReX>/Tools/Plotfile/ # This may change
   make -j 8
   # Run the executable to compare old and new versions
   <path to AMReX>/Tools/Plotfile/fcompare.gnu.ex old/plt00200 new/plt00200

which should return something like

.. code-block:: sh

             variable name             absolute error            relative error
                                          (||A - B||)         (||A - B||/||A||)
   ----------------------------------------------------------------------------
   level = 0
   jx                                 1.044455105e+11               1.021651316
   jy                                  4.08631977e+16               7.734299273
   jz                                 1.877301764e+14               1.073458933
   Ex                                 4.196315448e+10               1.253551615
   Ey                                 3.330698083e+12               6.436470137
   Ez                                 2.598167798e+10              0.6804387128
   Bx                                     273.8687473               2.340209782
   By                                     152.3911863                1.10952567
   Bz                                     37.43212767                 2.1977289
   part_per_cell                                   15                    0.9375
   Ex_fp                              4.196315448e+10               1.253551615
   Ey_fp                              3.330698083e+12               6.436470137
   Ez_fp                              2.598167798e+10              0.6804387128
   Bx_fp                                  273.8687473               2.340209782
   By_fp                                  152.3911863                1.10952567
   Bz_fp                                  37.43212767                 2.1977289
