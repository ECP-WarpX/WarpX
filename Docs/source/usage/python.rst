.. _usage-picmi:

Python (PICMI)
==============

WarpX uses the PICMI standard for its Python input files.
See `PICMI information and source code <https://github.com/picmi-standard/picmi>`__.

WarpX can be run in one of two modes. It can run as a preprocessor, using the
Python input file to generate an input file to be used by the C++ version, or
it can be run directly from Python.

In either mode, if using a `virtual environment <https://docs.python.org/3/tutorial/venv.html>`__, be sure to activate it before compiling and running WarpX.

Example input files can be found in :ref:`the examples section <usage-examples>`.
The examples support running in both modes by commenting and uncommenting the appropriate lines.

WarpX requires Python version 3.6 or newer.


Using Python input as a preprocessor
------------------------------------

In this case, only the pure Python version needs to be installed, as described :ref:`here <developers-gnumake-python>`.

In order to run a new simulation:

* Create a **new directory**, where the simulation will be run.

* Add a **Python script** in the directory.

The input file should have the line like ``sim.write_input_file(file_name = 'inputs_from_PICMI')``
which runs the preprocessor, generating the AMReX inputs file.

* **Run** the script with Python:

.. code-block:: bash

   python <python_script>

where ``<python_script>`` is the name of the script.
This creates the WarpX input file that you can run as normal with the WarpX executable.


Running WarpX directly from Python
----------------------------------

For this, a full Python installation of WarpX is required, as described in :ref:`the install documentation <install-users>` (:ref:`developers <install-developers>`).

In order to run a new simulation:

* Create a **new directory**, where the simulation will be run.

* Add a **Python script** in the directory.

The input file should have the line ``sim.step()`` which runs the simulation.

* **Run** the script with Python:

.. code-block:: bash

   mpirun -np <n_ranks> python <python_script>

where ``<n_ranks>`` is the number of MPI ranks used, and ``<python_script>``
is the name of the script.
