.. _examples-ohm-solver-ion-landau-damping:

Ohm Solver: Ion Landau Damping
==============================

Landau damping is a well known process in which electrostatic (acoustic) waves are damped by transferring energy to particles satisfying a resonance condition.
The process can be simulated by seeding a plasma with a specific acoustic mode (density perturbation) and tracking the strength of the mode as a function of time.
The figure below shows a set of such simulations with parameters matching those described in section 4.5 of :cite:t:`ex-MUNOZ2018`.
The straight lines show the theoretical damping rate for the given temperature ratios.

.. figure:: https://user-images.githubusercontent.com/40245517/230523935-3c8d63bd-ee69-4639-b111-f06dad5587f6.png
   :alt: Ion Landau damping
   :width: 70%

The input file for these examples and the corresponding analysis can be found at:

* :download:`Ion Landau damping input <PICMI_inputs.py>`
* :download:`Analysis script <analysis.py>`

The same input script can be used for 1d, 2d or 3d simulations and to sweep different
temperature ratios.

   .. code-block:: bash

      python3 PICMI_inputs.py -dim {1/2/3} --temp_ratio {value}
