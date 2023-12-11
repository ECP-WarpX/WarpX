.. _usage-picmi:
.. _usage-picmi-run:

Parameters: Python (PICMI)
==========================

This documents on how to use WarpX as a Python script (e.g., ``python3 PICMI_script.py``).

WarpX uses the `PICMI standard <https://github.com/picmi-standard/picmi>`__ for its Python input files.
Complete example input files can be found in :ref:`the examples section <usage-examples>`.

In the input file, instances of classes are created defining the various aspects of the simulation.
A variable of type :py:class:`pywarpx.picmi.Simulation` is the central object to which all other options are passed, defining the simulation time, field solver, registered species, etc.

Once the simulation is fully configured, it can be used in one of two modes.
**Interactive** use is the most common and can be :ref:`extended with custom runtime functionality <usage-picmi-extend>`:

.. tab-set::

   .. tab-item:: Interactive

      :py:meth:`~pywarpx.picmi.Simulation.step`: run directly from Python

   .. tab-item:: Preprocessor

      :py:meth:`~pywarpx.picmi.Simulation.write_input_file`: create an :ref:`inputs file for a WarpX executable <running-cpp-parameters>`


.. _usage-picmi-parameters:

Simulation and Grid Setup
-------------------------

.. autoclass:: pywarpx.picmi.Simulation
    :members: step, add_species, add_laser, write_input_file

.. autoclass:: pywarpx.picmi.Cartesian3DGrid

.. autoclass:: pywarpx.picmi.Cartesian2DGrid

.. autoclass:: pywarpx.picmi.Cartesian1DGrid

.. autoclass:: pywarpx.picmi.CylindricalGrid

.. autoclass:: pywarpx.picmi.EmbeddedBoundary

Field solvers define the updates of electric and magnetic fields.

.. autoclass:: pywarpx.picmi.ElectromagneticSolver

.. autoclass:: pywarpx.picmi.ElectrostaticSolver

Constants
---------

For convenience, the PICMI interface defines the following constants,
which can be used directly inside any PICMI script. The values are in SI units.

- ``picmi.constants.c``: The speed of light in vacuum.
- ``picmi.constants.ep0``: The vacuum permittivity :math:`\epsilon_0`
- ``picmi.constants.mu0``: The vacuum permeability :math:`\mu_0`
- ``picmi.constants.q_e``: The elementary charge (absolute value of the charge of an electron).
- ``picmi.constants.m_e``: The electron mass
- ``picmi.constants.m_p``: The proton mass

Applied fields
--------------

.. autoclass:: pywarpx.picmi.AnalyticInitialField

.. autoclass:: pywarpx.picmi.ConstantAppliedField

.. autoclass:: pywarpx.picmi.AnalyticAppliedField

.. autoclass:: pywarpx.picmi.PlasmaLens

.. autoclass:: pywarpx.picmi.Mirror

Diagnostics
-----------

.. autoclass:: pywarpx.picmi.ParticleDiagnostic

.. autoclass:: pywarpx.picmi.FieldDiagnostic

.. autoclass:: pywarpx.picmi.ElectrostaticFieldDiagnostic

.. autoclass:: pywarpx.picmi.Checkpoint

Lab-frame diagnostics diagnostics are used when running boosted-frame simulations.

.. autoclass:: pywarpx.picmi.LabFrameFieldDiagnostic

Particles
---------

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

Other operations related to particles:

.. autoclass:: pywarpx.picmi.CoulombCollisions

.. autoclass:: pywarpx.picmi.MCCCollisions

.. autoclass:: pywarpx.picmi.FieldIonization

Laser Pulses
------------

Laser profiles can be used to initialize laser pulses in the simulation.

.. autoclass:: pywarpx.picmi.GaussianLaser

.. autoclass:: pywarpx.picmi.AnalyticLaser

Laser injectors control where to initialize laser pulses on the simulation grid.

.. autoclass:: pywarpx.picmi.LaserAntenna


.. _usage-picmi-extend:

Extending a Simulation from Python
----------------------------------

When running WarpX directly from Python it is possible to interact with the simulation.

For instance, with the :py:meth:`~pywarpx.picmi.Simulation.step` method of the simulation class, one could run ``sim.step(nsteps=1)`` in a loop:

.. code-block:: python3

   # Preparation: set up the simulation
   #   sim = picmi.Simulation(...)
   #   ...

   steps = 1000
   for _ in range(steps):
       sim.step(nsteps=1)

       # do something custom with the sim object

As a more flexible alternative, one can install `callback functions <https://en.wikipedia.org/wiki/Callback_(computer_programming)>`__, which will execute a given Python function at a
specific location in the WarpX simulation loop.

.. automodule:: pywarpx.callbacks
   :members: installcallback, uninstallcallback, isinstalled

Data Access
^^^^^^^^^^^

While the simulation is running, callbacks can access the WarpX simulation data *in situ*.

An important object for data access is ``Simulation.extension.warpx``, which is available only during the simulation run.
This object is the Python equivalent to the C++ ``WarpX`` simulation class and provides access to field ``MultiFab`` and ``ParticleContainer`` data.

.. py:function:: pywarpx.picmi.Simulation.extension.warpx.getistep()

.. py:function:: pywarpx.picmi.Simulation.extension.warpx.gett_new()

.. py:function:: pywarpx.picmi.Simulation.extension.warpx.evolve()

.. autofunction:: pywarpx.picmi.Simulation.extension.finalize()

These and other classes are provided through `pyAMReX <https://github.com/AMReX-Codes/pyamrex>`__.
After the simulation is initialized, pyAMReX can be accessed via

.. code-block:: python

   from pywarpx import picmi, libwarpx

   # ... simulation definition ...

   # equivalent to
   #   import amrex.space3d as amr
   # for a 3D simulation
   amr = libwarpx.amr  # picks the right 1d, 2d or 3d variant

.. py:function:: amr.ParallelDescriptor.NProcs()

.. py:function:: amr.ParallelDescriptor.MyProc()

.. py:function:: amr.ParallelDescriptor.IOProcessor()

.. py:function:: amr.ParallelDescriptor.IOProcessorNumber()

Particles can be added to the simulation at specific positions and with specific attribute values:

.. code-block:: python

   from pywarpx import particle_containers, picmi

   # ...

   electron_wrapper = particle_containers.ParticleContainerWrapper("electrons")

.. autofunction:: pywarpx.particle_containers.ParticleContainerWrapper.add_particles

Properties of the particles already in the simulation can be obtained with various functions.

.. autofunction:: pywarpx.particle_containers.ParticleContainerWrapper.get_particle_count

.. autofunction:: pywarpx.particle_containers.ParticleContainerWrapper.get_particle_structs

.. autofunction:: pywarpx.particle_containers.ParticleContainerWrapper.get_particle_arrays

The ``get_particle_structs()`` and ``get_particle_arrays()`` functions are called
by several utility functions of the form ``get_particle_{comp_name}`` where
``comp_name`` is one of ``x``, ``y``, ``z``, ``r``, ``theta``, ``id``, ``cpu``,
``weight``, ``ux``, ``uy`` or ``uz``.

New components can be added via Python.

.. autofunction:: pywarpx.particle_containers.ParticleContainerWrapper.add_real_comp

Various diagnostics are also accessible from Python.
This includes getting the deposited or total charge density from a given species as well as accessing the scraped particle buffer.
See the example in ``Examples/Tests/ParticleBoundaryScrape`` for a reference on how to interact with scraped particle data.

.. autofunction:: pywarpx.particle_containers.ParticleContainerWrapper.get_species_charge_sum

.. autofunction:: pywarpx.particle_containers.ParticleContainerWrapper.deposit_charge_density

.. autofunction:: pywarpx.particle_containers.ParticleBoundaryBufferWrapper.get_particle_boundary_buffer_size

.. autofunction:: pywarpx.particle_containers.ParticleBoundaryBufferWrapper.get_particle_boundary_buffer_structs

.. autofunction:: pywarpx.particle_containers.ParticleBoundaryBufferWrapper.get_particle_boundary_buffer

.. autofunction:: pywarpx.particle_containers.ParticleBoundaryBufferWrapper.clear_buffer

The embedded boundary conditions can be modified when using the electrostatic solver.

.. py:function:: pywarpx.picmi.Simulation.extension.warpx.set_potential_on_eb()
