.. _examples-ohm-solver-em-modes:

Ohm Solver: Electromagnetic modes
=================================

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

* :download:`EM modes input <PICMI_inputs.py>`
* :download:`Analysis script <analysis.py>`

The same input script can be used for 1d, 2d or 3d simulations as well as replicating either the parallel propagating or
ion-Bernstein modes as indicated below.

   .. code-block:: bash

      python3 PICMI_inputs.py -dim {1/2/3} --bdir {x/y/z}

A RZ-geometry example case for normal modes propagating along an applied magnetic field in a cylinder is also available.
The analytical solution for these modes are described in :cite:t:`ex-Stix1992` Chapter 6, Sec. 2.

.. figure:: https://user-images.githubusercontent.com/40245517/259251824-33e78375-81d8-410d-a147-3fa0498c66be.png
   :alt: Normal EM modes in a metallic cylinder
   :width: 90%

The input file for this example and corresponding analysis can be found at:

* :download:`Cylindrical modes input <PICMI_inputs_rz.py>`
* :download:`Analysis script <analysis_rz.py>`
