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

* :download:`2D case <../../../Examples/Physics_applications/laser_acceleration/inputs_2d>`
* :download:`2D case in boosted frame <../../../Examples/Physics_applications/laser_acceleration/inputs_2d_boost>`
* :download:`3D case <../../../Examples/Physics_applications/laser_acceleration/inputs_3d>`

PICMI files:

* :download:`Without mesh refinement<../../../Examples/Physics_applications/laser_acceleration/PICMI_inputs_laser_acceleration.py>`

Plasma mirror
-------------

:download:`2D case <../../../Examples/Physics_applications/plasma_mirror/inputs_2d>`

Laser-ion acceleration
----------------------

:download:`2D case <../../../Examples/Physics_applications/laser_ion/inputs>`

.. note::

   The resolution of this 2D case is extremely low by default.
   You will need a computing cluster for adequate resolution of the target density, see comments in the input file.

Uniform plasma
--------------

:download:`2D case <../../../Examples/Physics_applications/uniform_plasma/inputs_2d>`
:download:`3D case <../../../Examples/Physics_applications/uniform_plasma/inputs_3d>`

Capacitive discharge
--------------------

The Monte-Carlo collision (MCC) model can be used to simulate capacitive discharges between parallel plates. The implementation has been tested against the benchmark results from Turner et. al. in `Phys. Plasmas 20, 013507, 2013 <https://aip.scitation.org/doi/abs/10.1063/1.4775084>`_. The figure below shows a comparison of the ion density as calculated in WarpX (in June 2021) compared to the literature results (which can be found `here <https://aip.scitation.org/doi/suppl/10.1063/1.4775084>`__). An input file with parameters for the 1st case is given below except the total simulation steps have been reduced.

.. figure:: https://user-images.githubusercontent.com/40245517/123883650-4f614680-d8fe-11eb-9302-92a282191e8f.png
   :alt: MCC benchmark against Turner et. al. (2013).
   :width: 80%

* :download:`test case 1 <../../../Examples/Physics_applications/capacitive_discharge/inputs_2d>`

.. note::

   This example needs `additional calibration data for cross sections <https://github.com/ECP-WarpX/warpx-data>`__.
   Download this data alongside your inputs file and update the paths in the inputs file:

   .. code-block:: bash

      git clone https://github.com/ECP-WarpX/warpx-data.git

Test cases
----------

PICMI test cases included that can be used as a reference:

* :download:`Gaussian beam <../../../Examples//Modules/gaussian_beam/PICMI_inputs_gaussian_beam.py>`
* :download:`Langmuir plasma wave test in 3d <../../../Examples//Tests/Langmuir/PICMI_inputs_langmuir_rt.py>`
* :download:`Langmuir plasma wave test in RZ <../../../Examples//Tests/Langmuir/PICMI_inputs_langmuir_rz_multimode_analyze.py>`
* :download:`Langmuir plasma wave test in 2D <../../../Examples//Tests/Langmuir/PICMI_inputs_langmuir2d.py>`
