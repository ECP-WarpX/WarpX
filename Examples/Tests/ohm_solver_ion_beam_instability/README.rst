.. _examples-ohm-solver-ion-beam-instability:

Ohm Solver: Ion Beam R Instability
==================================

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

* :download:`Ion beam R instability input <PICMI_inputs.py>`
* :download:`Analysis script <analysis.py>`

The same input script can be used for 1d, 2d or 3d simulations as well as replicating either the resonant or non-resonant
condition as indicated below.

   .. code-block:: bash

      python3 PICMI_inputs.py -dim {1/2/3} --resonant
