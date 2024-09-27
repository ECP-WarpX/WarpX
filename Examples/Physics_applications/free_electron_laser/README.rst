.. _examples-free-electron-laser:

Free-electron laser
===================

This example shows how to simulate the physics of a free-electron laser (FEL) using WarpX.
In this example, a relativistic electron beam is sent through an undulator (represented by an external,
oscillating magnetic field). The radiation emitted by the beam grows exponentially
as the beam travels through the undulator, due to the Free-Electron-Laser instability.

The parameters of the simulation are taken from section 5.1 of :cite:t:`ex-Fallahi2020`.

The simulation is performed in 1D, and uses the boosted-frame technique as described in
:cite:t:`ex-VayFELA2009` and :cite:t:`ex-VayFELB2009`, to reduce the computational cost of the simulation.
(The Lorentz frame of the simulation is moving at the average speed of the beam in the undulator.)
Even though the simulation is run in this boosted frame, the results are reconstructed in the
laboratory frame, using WarpX's ``BackTransformedDiagnostic``.

The effect of space-charge is intentionally turned off in this example, as it may not be properly modeled in 1D.
This is achieved by initializing two species of opposite charge (electrons and positrons) to
represent the physical electron beam, as discussed in :cite:t:`ex-VayFELB2009`.

Run
---

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. literalinclude:: inputs_test_1d_fel
   :language: ini
   :caption: You can copy this file from ``Examples/Physics_applications/beam_beam_collision/inputs_test_1d_fel``.
