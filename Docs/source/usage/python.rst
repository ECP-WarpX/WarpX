.. _usage-picmi:

Python (PICMI)
==============

WarpX uses the `PICMI standard <https://github.com/picmi-standard/picmi>`__ for its Python input files.
Python version 3.7 or newer is required.

Example input files can be found in :ref:`the examples section <usage-examples>`.
In the input file, instances of classes are created defining the various aspects of the simulation.
The `Simulation` object is the central object, where the instances are passed,
defining the simulation time, field solver, registered species, etc.

.. _usage-picmi-parameters:

Classes
----------

Simulation and grid setup
^^^^^^^^^^^^^^^^^^^^^^^^^

Simulation
""""""""""
.. autoclass:: pywarpx.picmi.Simulation
    :members: step, add_species, add_laser, write_input_file

Constants
"""""""""
For convenience, the PICMI interface defines the following constants,
which can be used directly inside any PICMI script. The values are in SI units.

- ``picmi.constants.c``: The speed of light in vacuum.
- ``picmi.constants.ep0``: The vacuum permittivity :math:`\epsilon_0`
- ``picmi.constants.mu0``: The vacuum permeability :math:`\mu_0`
- ``picmi.constants.q_e``: The elementary charge (absolute value of the charge of an electron).
- ``picmi.constants.m_e``: The electron mass
- ``picmi.constants.m_p``: The proton mass

Field solvers define the updates of electric and magnetic fields.

ElectromagneticSolver
"""""""""""""""""""""
.. autoclass:: pywarpx.picmi.ElectromagneticSolver

ElectrostaticSolver
"""""""""""""""""""
.. autoclass:: pywarpx.picmi.ElectrostaticSolver

Cartesian3DGrid
"""""""""""""""
.. autoclass:: pywarpx.picmi.Cartesian3DGrid

Cartesian2DGrid
"""""""""""""""
.. autoclass:: pywarpx.picmi.Cartesian2DGrid

Cartesian1DGrid
"""""""""""""""
.. autoclass:: pywarpx.picmi.Cartesian1DGrid

CylindricalGrid
"""""""""""""""
.. autoclass:: pywarpx.picmi.CylindricalGrid

EmbeddedBoundary
""""""""""""""""
.. autoclass:: pywarpx.picmi.EmbeddedBoundary

Applied fields
^^^^^^^^^^^^^^

ConstantAppliedField
""""""""""""""""""""
.. autoclass:: pywarpx.picmi.ConstantAppliedField

AnalyticAppliedField
""""""""""""""""""""
.. autoclass:: pywarpx.picmi.AnalyticAppliedField

PlasmaLens
""""""""""
.. autoclass:: pywarpx.picmi.PlasmaLens

Mirror
""""""
.. autoclass:: pywarpx.picmi.Mirror

Diagnostics
^^^^^^^^^^^

ParticleDiagnostic
""""""""""""""""""
.. autoclass:: pywarpx.picmi.ParticleDiagnostic

FieldDiagnostic
"""""""""""""""
.. autoclass:: pywarpx.picmi.FieldDiagnostic

ElectrostaticFieldDiagnostic
""""""""""""""""""""""""""""
.. autoclass:: pywarpx.picmi.ElectrostaticFieldDiagnostic

Lab-frame diagnostics diagnostics are used when running boosted-frame simulations.

LabFrameFieldDiagnostic
"""""""""""""""""""""""
.. autoclass:: pywarpx.picmi.LabFrameFieldDiagnostic

Checkpoint
""""""""""
.. autoclass:: pywarpx.picmi.Checkpoint

Particles
^^^^^^^^^

Species objects are a collection of particles with similar properties.
For instance, background plasma electrons, background plasma ions and an externally injected beam could each be their own particle species.

Species
"""""""
.. autoclass:: pywarpx.picmi.Species

MultiSpecies
""""""""""""
.. autoclass:: pywarpx.picmi.MultiSpecies

Particle distributions can be used for to initialize particles in a particle species.

GaussianBunchDistribution
"""""""""""""""""""""""""
.. autoclass:: pywarpx.picmi.GaussianBunchDistribution

UniformDistribution
"""""""""""""""""""
.. autoclass:: pywarpx.picmi.UniformDistribution

AnalyticDistribution
""""""""""""""""""""
.. autoclass:: pywarpx.picmi.AnalyticDistribution

ParticleListDistribution
""""""""""""""""""""""""
.. autoclass:: pywarpx.picmi.ParticleListDistribution

Particle layouts determine how to microscopically place macro particles in a grid cell.

GriddedLayout
"""""""""""""
.. autoclass:: pywarpx.picmi.GriddedLayout

PseudoRandomLayout
""""""""""""""""""
.. autoclass:: pywarpx.picmi.PseudoRandomLayout

Other operations related to particles

CoulombCollisions
"""""""""""""""""
.. autoclass:: pywarpx.picmi.CoulombCollisions

MCCCollisions
"""""""""""""
.. autoclass:: pywarpx.picmi.MCCCollisions

Lasers
^^^^^^

Laser profiles can be used to initialize laser pulses in the simulation.

GaussianLaser
"""""""""""""
.. autoclass:: pywarpx.picmi.GaussianLaser

AnalyticLaser
"""""""""""""
.. autoclass:: pywarpx.picmi.AnalyticLaser

Laser injectors control where to initialize laser pulses on the simulation grid.

LaserAntenna
""""""""""""
.. autoclass:: pywarpx.picmi.LaserAntenna


.. _usage-picmi-run:

Running
-------

WarpX can be run in one of two modes. It can run as a preprocessor, using the
Python input file to generate an input file to be used by the C++ version, or
it can be run directly from Python.
The examples support running in both modes by commenting and uncommenting the appropriate lines.

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


Extending a Simulation from Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When running WarpX directly from Python it is possible to interact with the simulation
by installing ``CallbackFunctions``, which will execute a given Python function at a
specific location in the WarpX simulation loop.

.. autoclass:: pywarpx.callbacks.CallbackFunctions

Places in the WarpX loop where callbacks are available include:
``afterinit``, ``beforecollisions``, ``aftercollisions``, ``beforeEsolve``, ``afterEsolve``,
``beforedeposition``, ``afterdeposition``, ``beforestep``, ``afterstep``, ``afterdiagnostics``,
``afterrestart`` and ``oncheckpointsignal``.
See the examples in *Examples/Tests/ParticleDataPython* for references on how to use
``callbacks``.

There are several "hooks" available via the ``libwarpx`` shared library to access and manipulate
simulation objects (particles, fields and memory buffers) as well as general properties
(such as processor number). These "hooks" are accessible through the `Simulation.extension` object.

.. autofunction:: pywarpx.picmi.Simulation.extension.getNProcs

.. autofunction:: pywarpx.picmi.Simulation.extension.getMyProc

.. autofunction:: pywarpx.picmi.Simulation.extension.get_nattr

.. autofunction:: pywarpx.picmi.Simulation.extension.get_nattr_species

.. autofunction:: pywarpx.picmi.Simulation.extension.getistep

.. autofunction:: pywarpx.picmi.Simulation.extension.gett_new

.. autofunction:: pywarpx.picmi.Simulation.extension.evolve

.. autofunction:: pywarpx.picmi.Simulation.extension.finalize

.. autofunction:: pywarpx.picmi.Simulation.extension.getProbLo

.. autofunction:: pywarpx.picmi.Simulation.extension.getProbHi

.. autofunction:: pywarpx.picmi.Simulation.extension.getCellSize

Particles can be added to the simulation at specific positions and with specific
attribute values:

.. autofunction:: pywarpx.picmi.Simulation.extension.add_particles

Properties of the particles already in the simulation can be obtained with various
functions.

.. autofunction:: pywarpx.picmi.Simulation.extension.get_particle_count

.. autofunction:: pywarpx.picmi.Simulation.extension.get_particle_structs

.. autofunction:: pywarpx.picmi.Simulation.extension.get_particle_arrays

The ``get_particle_structs()`` and ``get_particle_arrays()`` functions are called
by several utility functions of the form ``get_particle_{comp_name}`` where
``comp_name`` is one of ``x``, ``y``, ``z``, ``r``, ``theta``, ``id``, ``cpu``,
``weight``, ``ux``, ``uy`` or ``uz``.

The index of some specific component of the particle data can be obtained.

.. autofunction:: pywarpx.picmi.Simulation.extension.get_particle_comp_index

New components can be added via Python.

.. autofunction:: pywarpx.picmi.Simulation.extension.add_real_comp

Various diagnostics are also accessible from Python.
This includes getting the deposited or total charge density from a given species
as well as accessing the scraped particle buffer. See the example in
*Examples/Tests/ParticleBoudaryScrape* for a reference on how to interact
with scraped particle data.

.. autofunction:: pywarpx.picmi.Simulation.extension.get_species_charge_sum

.. autofunction:: pywarpx.picmi.Simulation.extension.depositChargeDensity

.. autofunction:: pywarpx.picmi.Simulation.extension.get_particle_boundary_buffer_size

.. autofunction:: pywarpx.picmi.Simulation.extension.get_particle_boundary_buffer_structs

.. autofunction:: pywarpx.picmi.Simulation.extension.get_particle_boundary_buffer

.. autofunction:: pywarpx.picmi.Simulation.extension.clearParticleBoundaryBuffer

The embedded boundary conditions can be modified when using the electrostatic solver.

.. autofunction:: pywarpx.picmi.Simulation.extension.set_potential_EB

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
