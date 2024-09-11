.. _multiphysics-ionization:

Ionization
==========

Field Ionization
----------------

Under the influence of a sufficiently strong external electric field atoms become ionized.
Particularly the dynamics of interactions between ultra-high intensity laser pulses and matter, e.g., Laser-Plasma Acceleration (LPA) with ionization injection, or Laser-Plasma Interactions with solid density targets (LPI) can depend on field ionization dynamics as well.

WarpX models field ionization based on a description of the Ammosov-Delone-Krainov model:cite:p:`mpion-Ammosov1986` following :cite:t:`mpion-ChenPRSTAB13`.

Implementation Details and Assumptions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

    The current implementation makes the following assumptions

    * Energy for ionization processes is not removed from the electromagnetic fields
    * Only one single-level ionization process can occur per macroparticle and time step
    * Ionization happens at the beginning of the PIC loop before the field solve
    * Angular momentum quantum number :math:`l = 0` and magnetic quantum number :math:`m = 0`

.. bibliography::
    :keyprefix: mpion-