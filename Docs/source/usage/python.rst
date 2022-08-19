.. _usage-picmi:

Python (PICMI)
==============

WarpX uses the `PICMI standard <https://github.com/picmi-standard/picmi>`__ for its Python input files.
Python version 3.6 or newer is required.

Example input files can be found in :ref:`the examples section <usage-examples>`.
The examples support running in both modes by commenting and uncommenting the appropriate lines.

.. _usage-picmi-parameters:

Parameters
----------

Simulation and grid setup
^^^^^^^^^^^^^^^^^^^^^^^^^

The `Simulation` object is the central object in a PICMI script.
It defines the simulation time, field solver, registered species, etc.

.. autoclass:: pywarpx.picmi.Simulation
    :members: step, add_species, add_laser, write_input_file

Field solvers define the updates of electric and magnetic fields.

.. autoclass:: pywarpx.picmi.ElectromagneticSolver

.. autoclass:: pywarpx.picmi.ElectrostaticSolver

Grid define the geometry and discretization.

.. autoclass:: pywarpx.picmi.Cartesian3DGrid

.. autoclass:: pywarpx.picmi.Cartesian2DGrid

.. autoclass:: pywarpx.picmi.Cartesian1DGrid

.. autoclass:: pywarpx.picmi.CylindricalGrid

.. autoclass:: pywarpx.picmi.EmbeddedBoundary

For convenience, the PICMI interface defines the following constants,
which can be used directly inside any PICMI script. The values are in SI units.

- ``picmi.constants.c``: The speed of light in vacuum.
- ``picmi.constants.ep0``: The vacuum permittivity :math:`\epsilon_0`
- ``picmi.constants.mu0``: The vacuum permeability :math:`\mu_0`
- ``picmi.constants.q_e``: The elementary charge (absolute value of the charge of an electron).
- ``picmi.constants.m_e``: The electron mass
- ``picmi.constants.m_p``: The proton mass

Additionally to self-consistent fields from the field solver, external fields can be applied.

.. autoclass:: pywarpx.picmi.ConstantAppliedField

.. autoclass:: pywarpx.picmi.AnalyticAppliedField

.. autoclass:: pywarpx.picmi.PlasmaLens

.. autoclass:: pywarpx.picmi.Mirror

Diagnostics can be used to output data.

.. autoclass:: pywarpx.picmi.ParticleDiagnostic

.. autoclass:: pywarpx.picmi.FieldDiagnostic

.. autoclass:: pywarpx.picmi.ElectrostaticFieldDiagnostic

Lab-frame diagnostics diagnostics are used when running boosted-frame simulations.

.. autoclass:: pywarpx.picmi.LabFrameParticleDiagnostic

.. autoclass:: pywarpx.picmi.LabFrameFieldDiagnostic

Checkpoints can be controlled

.. autoclass:: pywarpx.picmi.Checkpoint

Particles
^^^^^^^^^

Species objects are a collection of particles with similar properties.
For instance, background plasma electrons, background plasma ions and an externally injected beam could each be their own particle species.

.. autoclass:: pywarpx.picmi.Species

.. autoclass:: pywarpx.picmi.MultiSpecies

Particle distributions can be used for to initialize particles in a particle species.

.. autoclass:: pywarpx.picmi.GaussianBunchDistribution

.. autoclass:: pywarpx.picmi.UniformDistribution

.. autoclass:: pywarpx.picmi.AnalyticDistribution

.. autoclass:: pywarpx.picmi.ParticleListDistribution

Particle layouts determine how to microscopically place macro particles in a grid cell.

.. autoclass:: pywarpx.picmi.GriddedLayout

.. autoclass:: pywarpx.picmi.PseudoRandomLayout

Other operations related to particles

.. autoclass:: pywarpx.picmi.CoulombCollisions

.. autoclass:: pywarpx.picmi.MCCCollisions

Lasers
^^^^^^

Laser profiles can be used to initialize laser pulses in the simulation.

.. autoclass:: pywarpx.picmi.GaussianLaser

.. autoclass:: pywarpx.picmi.AnalyticLaser

Laser injectors control where to initialize laser pulses on the simulation grid.

.. autoclass:: pywarpx.picmi.LaserAntenna


.. _usage-picmi-run:

Running
-------

WarpX can be run in one of two modes. It can run as a preprocessor, using the
Python input file to generate an input file to be used by the C++ version, or
it can be run directly from Python.

In either mode, if using a `virtual environment <https://docs.python.org/3/tutorial/venv.html>`__, be sure to activate it before compiling and running WarpX.


Running WarpX directly from Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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

Using Python input as a preprocessor
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
