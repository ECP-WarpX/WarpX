.. _examples-pwfa:

Beam-Driven Wakefield Acceleration of Electrons
===============================================

This example shows how to model a beam-driven plasma-wakefield accelerator (PWFA) :cite:p:`ex-TajimaDawson1982,ex-Esarey1996`.

PWFA is best performed in 3D or quasi-cylindrical (RZ) geometry, in order to correctly capture some of the key physics (structure of the space-charge fields, beamloading, shape of the accelerating bubble in the blowout regime, etc.).
For physical situations that have close-to-cylindrical symmetry, simulations in RZ geometry capture the relevant physics at a fraction of the computational cost of a 3D simulation.
On the other hand, for physical situation with strong asymmetries (e.g., non-round driver, strong hosing of the accelerated beam, etc.), only 3D simulations are suitable.

Additionally, to speed up computation, this example uses the :ref:`boosted frame method <theory-boostedframe>` to effectively model long acceleration lengths.

Alternatively, an other common approximation for PWFAs is quasi-static modeling, e.g., if effects such as self-injection can be ignored.
In the Beam, Plasma & Accelerator Simulation Toolkit (BLAST), `HiPACE++ <https://hipace.readthedocs.io>`__ provides such methods.

.. note::

    There are two FDTD Maxwell field solvers that compute the field push implemented in WarpX: the Yee and Cole-Karkkainen solver with Cowan coefficients (CKC) solvers.
    The later includes a modification that allows the numerical dispersion of light in vacuum to be exact, and that is why we choose CKC for the example.

Lorentz boosted frame
---------------------

    WarpX simulations can be done in the laboratory or `Lorentz-boosted <https://warpx.readthedocs.io/en/latest/theory/boosted_frame.html>`_ frames.
    In the laboratory frame, there is typically no need to model the plasma ions species, since they are mainly stationary during the short time scales associated with the motion of plasma electrons.
    In the boosted frame, that argument is no longer valid, as ions have relativistic velocities.

    To avoid spurious effects in the boosted frame, we want the longitudinal cell size to be larger than the transverse one.
    Translating this condition to the cell transverse (:math:`d_{x}`) and longitudinal dimensions (:math:`d_{z}`) in the laboratory frame leads to: :math:`d_{x} > (d_{z} (1+\beta_{b}) \gamma_{b})`, 
    where :math:`\beta_{b}` is the boosted frame velocity in units of :math:`c`.

    To determine the total number of time steps of the simulation, we could either set the `<zmax_plasma_to_compute_max_step>` parameter to the end of the plasma (:math:`z_{\textrm{end}}`), or compute it using:

    * boosted frame edge of the simulation box, :math:`\textrm{corner} = l_{e}/ ((1-\beta_{b}) \gamma_{b})`
    * time of interaction in the boosted frame, :math:`T = \frac{z_{\textrm{end}}/\gamma_{b}-\textrm{corner}}{c (1+\beta_{b})}`
    * total number of iterations, :math:`i_{\textrm{max}} = T/dt`

    where :math:`l_{e}` is the position of the left edge of the simulation box (in respect to propagation direction).

.. note::

   TODO: The Python (PICMI) input file should use the boosted frame method, like the ``inputs_3d_boost`` file.


Run
---

This example can be run **either** as:

* **Python** script: ``python3 PICMI_inputs_plasma_acceleration_3d.py`` or
* WarpX **executable** using an input file: ``warpx.3d inputs_3d_boost``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

      .. note::

         TODO: This input file should use the boosted frame method, like the ``inputs_3d_boost`` file.

      .. literalinclude:: PICMI_inputs_plasma_acceleration.py
         :language: python3
         :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/PICMI_inputs_plasma_acceleration.py``.

   .. tab-item:: Executable: Input File

      .. literalinclude:: inputs_3d_boost
         :language: ini
         :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/inputs_3d_boost``.

Visualize
---------

The plasma has a density profile that consists of a long up-ramp followed by a uniform plateau.

With this configuration the driver excites a nonlinear plasma wake and drives the bubble depleted of plasma electrons where the beam accelerates, as can be seen in Fig
.. figure:: https://user-images.githubusercontent.com/1353258/287800852-f994a020-4ecc-4987-bffc-2cb7df6144a9.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDE3MTYyNjksIm5iZiI6MTcwMTcxNTk2OSwicGF0aCI6Ii8xMzUzMjU4LzI4NzgwMDg1Mi1mOTk0YTAyMC00ZWNjLTQ5ODctYmZmYy0yY2I3ZGY2MTQ0YTkucG5nP1gtQW16LUFsZ29yaXRobT1BV1M0LUhNQUMtU0hBMjU2JlgtQW16LUNyZWRlbnRpYWw9QUtJQUlXTkpZQVg0Q1NWRUg1M0ElMkYyMDIzMTIwNCUyRnVzLWVhc3QtMSUyRnMzJTJGYXdzNF9yZXF1ZXN0JlgtQW16LURhdGU9MjAyMzEyMDRUMTg1MjQ5WiZYLUFtei1FeHBpcmVzPTMwMCZYLUFtei1TaWduYXR1cmU9NDkyNWJkMTg2NWM3ZjcwZjVkMjlmNDE1NmRjNWEyZWM5MzgxMWJhZTVjMGMxNjdkZDg1Zjk0NmQ1NGEwMjNiMiZYLUFtei1TaWduZWRIZWFkZXJzPWhvc3QmYWN0b3JfaWQ9MCZrZXlfaWQ9MCZyZXBvX2lkPTAifQ.C_NQceQcqiCDzBoSzIjm3c8QdTLNDdtjJmkQjkhW4c8
   
   :alt: Electron beam, in black, accelerated by wakefield, in red-blue colormap, driven by drive electron beam, in purple.  Set max_steps to 800 in WarpX simulation to generate data for this figure.

   Electron beam, in black, accelerated by wakefield, in red-blue colormap, driven by drive electron beam, in purple.  Set ``max_steps`` to 800 in a WarpX simulation to generate data for this figure.

The beam increases in energy as it propagates in the wake, as shown in Fig. 

.. figure:: https://user-images.githubusercontent.com/1353258/287800852-f994a020-4ecc-4987-bffc-2cb7df6144a9.png?jwt=eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJpc3MiOiJnaXRodWIuY29tIiwiYXVkIjoicmF3LmdpdGh1YnVzZXJjb250ZW50LmNvbSIsImtleSI6ImtleTEiLCJleHAiOjE3MDE3MTYyNjksIm5iZiI6MTcwMTcxNTk2OSwicGF0aCI6Ii8xMzUzMjU4LzI4NzgwMDg1Mi1mOTk0YTAyMC00ZWNjLTQ5ODctYmZmYy0yY2I3ZGY2MTQ0YTkucG5nP1gtQW16LUFsZ29yaXRobT1BV1M0LUhNQUMtU0hBMjU2JlgtQW16LUNyZWRlbnRpYWw9QUtJQUlXTkpZQVg0Q1NWRUg1M0ElMkYyMDIzMTIwNCUyRnVzLWVhc3QtMSUyRnMzJTJGYXdzNF9yZXF1ZXN0JlgtQW16LURhdGU9MjAyMzEyMDRUMTg1MjQ5WiZYLUFtei1FeHBpcmVzPTMwMCZYLUFtei1TaWduYXR1cmU9NDkyNWJkMTg2NWM3ZjcwZjVkMjlmNDE1NmRjNWEyZWM5MzgxMWJhZTVjMGMxNjdkZDg1Zjk0NmQ1NGEwMjNiMiZYLUFtei1TaWduZWRIZWFkZXJzPWhvc3QmYWN0b3JfaWQ9MCZrZXlfaWQ9MCZyZXBvX2lkPTAifQ.C_NQceQcqiCDzBoSzIjm3c8QdTLNDdtjJmkQjkhW4c8
   :alt: Energy gain of electron beam in wake through this simulation (this is the mean beam energy in the boosted frame, NOT the lab frame). Set ``max_steps`` to 800 in a WarpX simulation to generate data for this figure.

   Energy gain of electron beam in wake through this simulation (this is the mean beam energy in the boosted frame, NOT the lab frame). Set ``max_steps`` to 800 in a WarpX simulation to generate data for this figure.

The plots can be generated with the following script:

.. dropdown:: Script ``plot_3d.py``

   .. literalinclude:: plot_3d.py
      :language: python3
      :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/plot_3d.py``.


