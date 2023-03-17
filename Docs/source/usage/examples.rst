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

The Monte-Carlo collision (MCC) model can be used to simulate electron and ion collisions with a neutral background gas. In particular this can be used to study capacitive discharges between parallel plates. The implementation has been tested against the benchmark results from Turner et al. in `Phys. Plasmas 20, 013507, 2013 <https://aip.scitation.org/doi/abs/10.1063/1.4775084>`_. The figure below shows a comparison of the ion density as calculated in WarpX (in June 2022 with `PR #3118 <https://github.com/ECP-WarpX/WarpX/pull/3118>`_) compared to the literature results (which can be found `here <https://aip.scitation.org/doi/suppl/10.1063/1.4775084>`__).

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
