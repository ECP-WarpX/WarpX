.. _developers-testing:

Testing the code
================

When adding a new feature, you want to make sure that (i) you did not break the existing code and (ii) your contribution gives correct results. While the code is tested regularly remotely (on the cloud when commits are pushed to an open PR, and every night on local clusters), it can also be useful to run tests on your custom input file. This section details how to use both automated and custom tests.

Continuous Integration in WarpX
-------------------------------

Configuration
^^^^^^^^^^^^^

Our regression tests are run with `CTest <https://cmake.org/cmake/help/book/mastering-cmake/chapter/Testing%20With%20CMake%20and%20CTest.html#>`__, an executable that comes with CMake.

The test suite is ready to run once you have configured and built WarpX with CMake, following the instructions that you find in our :ref:`Users <install-cmake>` or :ref:`Developers <building-cmake>` sections.

A test that requires a build option that was not configured and built will be skipped automatically. For example, if you configure and build WarpX in 1D only, any test of dimensionality other than 1D, which would require WarpX to be configured and built in the corresponding dimensionality, will be skipped automatically.

How to run pre-commit tests locally
-----------------------------------

When proposing code changes to Warpx, we perform a couple of automated stylistic and correctness checks on the code change.
You can run those locally before you push to save some time, install them once like this:

.. code-block:: sh

   python -m pip install -U pre-commit
   pre-commit install

See `pre-commit.com <https://pre-commit.com>`__ and our ``.pre-commit-config.yaml`` file in the repository for more details.

How to run automated tests locally
----------------------------------

Once your new feature is ready, there are ways to check that you did not break anything.
WarpX has automated tests running every time a commit is pushed to an open pull request.
The input files and scripts used by the automated tests can be found in the `Examples <https://github.com/ECP-WarpX/WarpX/tree/development/Examples>`__ directory, either under `Physics_applications <https://github.com/ECP-WarpX/WarpX/tree/development/Examples/Physics_applications>`__ or `Tests <https://github.com/ECP-WarpX/WarpX/tree/development/Examples/Tests>`__.

For easier debugging, it can be convenient to run the tests on your local machine by executing CTest as illustrated in the examples below (where we assume that WarpX was configured and built in the directory ``build``):

* List tests available for the current build options:

  .. code-block:: sh

        ctest --test-dir build -N

* Run tests available for the current build options:

  .. code-block:: sh

        ctest --test-dir build

* Run tests available for the current build options in parallel (while preserving existing dependencies between tests):

  .. code-block:: sh

        ctest --test-dir build -j 2

* Run tests available for the current build options and output anything outputted by the test program if the test should fail:

  .. code-block:: sh

        ctest --test-dir build --output-on-failure

* Run tests available for the current build options with verbose output:

  .. code-block:: sh

        ctest --test-dir build --verbose

* Run tests matching the regular expression ``laser_acceleration``:

  .. code-block:: sh

        ctest --test-dir build -R laser_acceleration

* Run tests except those matching the regular expression ``laser_acceleration``:

  .. code-block:: sh

        ctest --test-dir build -E laser_acceleration

Once the execution of CTest is completed, you can find all files associated with each test in its corresponding directory under ``build/bin/``.
For example, if you run the single test ``test_3d_laser_acceleration``, you can find all files associated with this test in the directory ``build/bin/test_3d_laser_acceleration/``.

If you modify the code base locally and want to assess the effects of your code changes on the automated tests, you need to first rebuild WarpX including your code changes and then rerun CTest.

How to add automated tests
--------------------------

As mentioned above, the input files and scripts used by the automated tests can be found in the `Examples <https://github.com/ECP-WarpX/WarpX/tree/development/Examples>`__ directory, either under `Physics_applications <https://github.com/ECP-WarpX/WarpX/tree/development/Examples/Physics_applications>`__ or `Tests <https://github.com/ECP-WarpX/WarpX/tree/development/Examples/Tests>`__.

Each test directory must contain a file named ``CMakeLists.txt`` where all tests associated with the input files and scripts in that directory must be listed.

A new test can be added by adding a corresponding entry in ``CMakeLists.txt`` as illustrated in the examples below:

* Add the **regular test** ``test_1d_laser_acceleration``:

  .. code-block:: cmake

       add_warpx_test(
           test_1d_laser_acceleration  # name
           1  # dims
           2  # nprocs
           OFF  # eb
           inputs_test_1d_laser_acceleration  # inputs
           analysis.py  # analysis
           diags/diag1000100  # output (plotfile)
           OFF  # dependency
       )

* Add the **PICMI test** ``test_2d_laser_acceleration_picmi``:

  .. code-block:: cmake

       add_warpx_test(
           test_2d_laser_acceleration_picmi  # name
           2  # dims
           2  # nprocs
           OFF  # eb
           inputs_test_2d_laser_acceleration_picmi.py  # inputs
           analysis.py  # analysis
           diags/diag1000100  # output (plotfile)
           OFF  # dependency
       )

* Add the **restart test** ``test_3d_laser_acceleration_restart``:

  .. code-block:: cmake

       add_warpx_test(
           test_3d_laser_acceleration_restart  # name
           3  # dims
           2  # nprocs
           OFF  # eb
           inputs_test_3d_laser_acceleration_restart  # inputs
           analysis_default_restart.py  # analysis
           diags/diag1000100  # output (plotfile)
           test_3d_laser_acceleration  # dependency
       )

  Note that the restart has an explicit dependency, namely it can run only provided that the original test, from which the restart checkpoint files will be read, runs first.

* A more complex example. Add the **PICMI test** ``test_rz_laser_acceleration_picmi``, with custom command-line arguments ``--test`` and ``dir``, and openPMD time series output:

  .. code-block:: cmake

       add_warpx_test(
           test_rz_laser_acceleration_picmi  # name
           RZ  # dims
           2   # nprocs
           OFF  # eb
           "inputs_test_rz_laser_acceleration_picmi.py --test --dir 1"  # inputs
           analysis.py  # analysis
           diags/diag1/  # output (openPMD time series)
           OFF  # dependency
       )

If you need a new Python package dependency for testing, please add it in `Regression/requirements.txt <https://github.com/ECP-WarpX/WarpX/blob/development/Regression/requirements.txt>`__.

Sometimes, two tests share a large number of input parameters. The shared input parameters can be collected in a "base" input file that can be passed as a runtime parameter in the actual test input files through the parameter ``FILE``.

Naming conventions for automated tests
--------------------------------------

Note that we currently obey the following snake\_case naming conventions for test names and test input files (which make automation tasks easier, e.g., parsing visually, parsing through code, sorting alphabetically, filtering tests in CTest via ``-R``, etc.):

#. **Regular test names** start with the string ``test_1d_``, ``test_2d_``, ``test_3d_`` or ``test_rz_``, followed by a string that is descriptive of the test. For example, ``test_3d_laser_acceleration``.

#. **PICMI test names** start with the string ``test_1d_``, ``test_2d_``, ``test_3d_`` or ``test_rz_``, followed by a string that is descriptive of the test, and end with the string ``_picmi``. For example, ``test_3d_laser_acceleration_picmi``.

#. **Restart test names** end with the string ``_restart``. For example, ``test_3d_laser_acceleration_restart``.

#. **Test input files** start with the string ``inputs_`` followed by the test name. For example, ``inputs_test_3d_laser_acceleration`` or ``inputs_test_3d_laser_acceleration_picmi.py`` or ``inputs_test_3d_laser_acceleration_restart``.

#. **Base input files** (that is, files collecting input parameters shared between two or more tests) are typically named ``inputs_base_1d``, ``inputs_base_2d``, ``inputs_base_3d`` or ``inputs_base_rz``, possibly followed by additional strings if need be.

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
