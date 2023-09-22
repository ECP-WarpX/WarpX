.. _usage-examples:

Examples
========

This section allows you to **download input files** that correspond to different physical situations.

We provide two kinds of inputs:

* AMReX ``inputs`` files, :ref:`with parameters described here <running-cpp-parameters>`,
* PICMI python input files, `with parameters described here <https://picmi-standard.github.io>`__.

For a complete list of all example input files, have a look at our ``Examples/`` directory.
It contains folders and subfolders with self-describing names that you can try. All these input files are automatically tested, so they should always be up-to-date.

Beam-driven electron acceleration
---------------------------------

AMReX ``inputs``:

* :download:`2D case <../../../Examples/Physics_applications/plasma_acceleration/inputs_2d>`
* :download:`2D case in boosted frame <../../../Examples/Physics_applications/plasma_acceleration/inputs_2d_boost>`
* :download:`3D case in boosted frame <../../../Examples/Physics_applications/plasma_acceleration/inputs_3d_boost>`

PICMI:

* :download:`Without mesh refinement <../../../Examples/Physics_applications/plasma_acceleration/PICMI_inputs_plasma_acceleration.py>`
* :download:`With mesh refinement <../../../Examples/Physics_applications/plasma_acceleration/PICMI_inputs_plasma_acceleration_mr.py>`

Laser-driven electron acceleration
----------------------------------

AMReX ``inputs``:

* :download:`1D case <../../../Examples/Physics_applications/laser_acceleration/inputs_1d>`
* :download:`2D case <../../../Examples/Physics_applications/laser_acceleration/inputs_2d>`
* :download:`2D case in boosted frame <../../../Examples/Physics_applications/laser_acceleration/inputs_2d_boost>`
* :download:`3D case <../../../Examples/Physics_applications/laser_acceleration/inputs_3d>`
* :download:`RZ case <../../../Examples/Physics_applications/laser_acceleration/inputs_rz>`

PICMI (Python) scripts:

* :download:`1D case <../../../Examples/Physics_applications/laser_acceleration/PICMI_inputs_1d.py>`
* :download:`2D case with mesh refinement <../../../Examples/Physics_applications/laser_acceleration/PICMI_inputs_2d.py>`
* :download:`3D case <../../../Examples/Physics_applications/laser_acceleration/PICMI_inputs_3d.py>`
* :download:`RZ case <../../../Examples/Physics_applications/laser_acceleration/PICMI_inputs_rz.py>`

Plasma mirror
-------------

:download:`2D case <../../../Examples/Physics_applications/plasma_mirror/inputs_2d>`

Laser-ion acceleration
----------------------

:download:`2D case <../../../Examples/Physics_applications/laser_ion/inputs>`

.. note::

   The resolution of this 2D case is extremely low by default.
   You will need a computing cluster for adequate resolution of the target density, see comments in the input file.

.. warning::

   It is strongly advised to set the parameters ``<species>.zmin / zmax / xmin / ...`` when working with highly dense targets that are limited in one or multiple dimensions.
   The particle creation routine will first create particles everywhere between these limits (`defaulting to box size if unset`), setting particles to invalid only afterwards based on the density profile.
   Not setting these parameters can quickly lead to memory overflows.

Uniform plasma
--------------

:download:`2D case <../../../Examples/Physics_applications/uniform_plasma/inputs_2d>`
:download:`3D case <../../../Examples/Physics_applications/uniform_plasma/inputs_3d>`

Capacitive discharge
--------------------

The Monte-Carlo collision (MCC) model can be used to simulate electron and ion collisions with a neutral background gas. In particular this can be used to study capacitive discharges between parallel plates. The implementation has been tested against the benchmark results from :cite:t:`ex-Turner2013`. The figure below shows a comparison of the ion density as calculated in WarpX (in June 2022 with `PR #3118 <https://github.com/ECP-WarpX/WarpX/pull/3118>`_) compared to the literature results (which can be found `here <https://aip.scitation.org/doi/suppl/10.1063/1.4775084>`__).

.. figure:: https://user-images.githubusercontent.com/40245517/171573007-f7d733c7-c0de-490c-9ed6-ff4c02154358.png
   :alt: MCC benchmark against Turner et. al. (2013).
   :width: 80%

An input file to reproduce the benchmark calculations is linked below.
To run a given case ``-n``, from 1 to 4, execute:

   .. code-block:: bash

      python3 PICMI_inputs_1d.py -n 1

Once the simulation completes an output file ``avg_ion_density.npy`` will be created which can be compared to the literature results as in the plot above. Running case 1 on 4 processors takes roughly 20 minutes to complete.

* :download:`input file <../../../Examples/Physics_applications/capacitive_discharge/PICMI_inputs_1d.py>`

.. note::

   This example needs `additional calibration data for cross sections <https://github.com/ECP-WarpX/warpx-data>`__.
   Download this data alongside your inputs file and update the paths in the inputs file:

   .. code-block:: bash

      git clone https://github.com/ECP-WarpX/warpx-data.git

Test cases
----------

PICMI (Python) test cases included that can be used as a reference:

* :download:`Gaussian beam <../../../Examples/Tests/gaussian_beam/PICMI_inputs_gaussian_beam.py>`
* :download:`Langmuir plasma wave test in 3d <../../../Examples/Tests/langmuir/PICMI_inputs_langmuir_rt.py>`
* :download:`Langmuir plasma wave test in RZ <../../../Examples/Tests/langmuir/PICMI_inputs_langmuir_rz_multimode_analyze.py>`
* :download:`Langmuir plasma wave test in 2D <../../../Examples/Tests/langmuir/PICMI_inputs_langmuir2d.py>`

Manipulating fields via Python
------------------------------

An example of using Python to access the simulation charge density, solve the Poisson equation (using ``superLU``) and write the resulting electrostatic potential back to the simulation is given in the input file below. This example uses the ``fields.py`` module included in the ``pywarpx`` library.

* :download:`Direct Poisson solver example <../../../Examples/Physics_applications/capacitive_discharge/PICMI_inputs_2d.py>`

An example of initializing the fields by accessing their data through Python, advancing the simulation for a chosen number of time steps, and plotting the fields again through Python. The simulation runs with 128 regular cells, 8 guard cells, and 10 PML cells, in each direction. Moreover, it uses div(E) and div(B) cleaning both in the regular grid and in the PML and initializes all available electromagnetic fields (E,B,F,G) identically.

* :download:`Unit pulse with PML <../../../Examples/Tests/python_wrappers/PICMI_inputs_2d.py>`

.. _examples-hybrid-model:

Kinetic-fluid hybrid model
--------------------------

Several examples and benchmarks of the kinetic-fluid hybrid model are shown below. The first few examples are replications
of the verification tests described in :cite:t:`ex-MUNOZ2018`.
The hybrid-PIC model was added to WarpX in `PR #3665 <https://github.com/ECP-WarpX/WarpX/pull/3665>`_ - the figures below
were generated at that time.

Electromagnetic modes
^^^^^^^^^^^^^^^^^^^^^

In this example a simulation is seeded with a thermal plasma while an initial magnetic field is applied in either the
:math:`z` or :math:`x` direction. The simulation is progressed for a large number of steps and the resulting fields are
analyzed for mode excitations.

Right and left circularly polarized electromagnetic waves are supported through the cyclotron motion of the ions, except
in a region of thermal resonances as indicated on the plot below.

.. figure:: https://user-images.githubusercontent.com/40245517/216207688-9c39374a-9e69-45b8-a588-35b087b83d27.png
   :alt: Parallel EM modes in thermal ion plasma
   :width: 70%

Perpendicularly propagating modes are also supported, commonly referred to as ion-Bernstein modes.

.. figure:: https://user-images.githubusercontent.com/40245517/231217944-7d12b8d4-af4b-44f8-a1b9-a2b59ce3a1c2.png
   :alt: Perpendicular EM modes in thermal ion plasma
   :width: 50%

The input file for these examples and the corresponding analysis can be found at:

* :download:`EM modes input <../../../Examples/Tests/ohm_solver_EM_modes/PICMI_inputs.py>`
* :download:`Analysis script <../../../Examples/Tests/ohm_solver_EM_modes/analysis.py>`

The same input script can be used for 1d, 2d or 3d simulations as well as replicating either the parallel propagating or
ion-Bernstein modes as indicated below.

   .. code-block:: bash

      python3 PICMI_inputs.py -dim {1/2/3} --bdir {x/y/z}

Ion beam R instability
^^^^^^^^^^^^^^^^^^^^^^

In this example a low density ion beam interacts with a "core" plasma population which induces an instability.
Based on the relative density between the beam and the core plasma a resonant or non-resonant condition can
be accessed. The figures below show the evolution of the y-component of the magnetic field as the beam and
core plasma interact.

.. figure:: https://user-images.githubusercontent.com/40245517/217923933-6bdb65cb-7d26-40d8-8687-7dd75274bd48.png
   :alt: Resonant ion beam R instability
   :width: 70%

.. figure:: https://user-images.githubusercontent.com/40245517/217925983-b91d6482-69bc-43c1-8c7d-23ebe7c69d49.png
   :alt: Non-resonant ion beam R instability
   :width: 70%

The growth rates of the strongest growing modes for the resonant case are compared
to theory (dashed lines) in the figure below.

.. figure:: https://github.com/ECP-WarpX/WarpX/assets/40245517/a94bb6e5-30e9-4d8f-9e6b-844dc8f51d17
   :alt: Resonant ion beam R instability growth rates
   :width: 50%

The input file for these examples and the corresponding analysis can be found at:

* :download:`Ion beam R instability input <../../../Examples/Tests/ohm_solver_ion_beam_instability/PICMI_inputs.py>`
* :download:`Analysis script <../../../Examples/Tests/ohm_solver_ion_beam_instability/analysis.py>`

The same input script can be used for 1d, 2d or 3d simulations as well as replicating either the resonant or non-resonant
condition as indicated below.

   .. code-block:: bash

      python3 PICMI_inputs.py -dim {1/2/3} --resonant

Ion Landau damping
^^^^^^^^^^^^^^^^^^

Landau damping is a well known process in which electrostatic (acoustic) waves
are damped by transferring energy to particles satisfying a resonance condition.
The process can be simulated by seeding a plasma with a specific acoustic mode
(density perturbation) and tracking the strength of the mode as a function of
time. The figure below shows a set of such simulations with parameters matching
those described in section 4.5 of :cite:t:`ex-MUNOZ2018`. The straight lines show
the theoretical damping rate for the given temperature ratios.

.. figure:: https://user-images.githubusercontent.com/40245517/230523935-3c8d63bd-ee69-4639-b111-f06dad5587f6.png
   :alt: Ion Landau damping
   :width: 70%

The input file for these examples and the corresponding analysis can be found at:

* :download:`Ion Landau damping input <../../../Examples/Tests/ohm_solver_ion_Landau_damping/PICMI_inputs.py>`
* :download:`Analysis script <../../../Examples/Tests/ohm_solver_ion_Landau_damping/analysis.py>`

The same input script can be used for 1d, 2d or 3d simulations and to sweep different
temperature ratios.

   .. code-block:: bash

      python3 PICMI_inputs.py -dim {1/2/3} --temp_ratio {value}

Magnetic reconnection
^^^^^^^^^^^^^^^^^^^^^

Hybrid-PIC codes are often used to simulate magnetic reconnection in space
plasmas. An example of magnetic reconnection from a force-free sheet is
provided, based on the simulation described in :cite:t:`ex-Le2016`.

.. figure:: https://user-images.githubusercontent.com/40245517/229639784-b5d3b596-3550-4570-8761-8d9a67aa4b3b.gif
   :alt: Magnetic reconnection
   :width: 70%

The input file for this example and corresponding analysis can be found at:

* :download:`Magnetic reconnection input <../../../Examples/Tests/ohm_solver_magnetic_reconnection/PICMI_inputs.py>`
* :download:`Analysis script <../../../Examples/Tests/ohm_solver_magnetic_reconnection/analysis.py>`

Many Further Examples, Demos and Tests
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

WarpX runs over 200 integration tests on a variety of modeling cases, which validate and demonstrate its functionality.
Please see the `Examples/Tests/ <https://github.com/ECP-WarpX/WarpX/tree/development/Examples/Tests>`__ directory for many more examples.


.. bibliography::
    :keyprefix: ex-
