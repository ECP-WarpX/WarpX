.. _examples-pwfa:

Beam-Driven Wakefield Acceleration of Electrons
===============================================

This example shows how to model a beam-driven plasma-wakefield accelerator (PWFA) :cite:p:`ex-TajimaDawson1982,ex-Esarey1996`.

PWFA is best performed in 3D or quasi-cylindrical (RZ) geometry, in order to correctly capture some of the key physics (structure of the space-charge fields, beamloading, shape of the accelerating bubble in the blowout regime, etc.).
For physical situations that have close-to-cylindrical symmetry, simulations in RZ geometry capture the relevant physics at a fraction of the computational cost of a 3D simulation.
On the other hand, for physical situation with strong asymmetries (e.g., non-round driver, strong hosing of the accelerated beam, etc.), only 3D simulations are suitable.

Additionally, to speed up computation, this example uses the :ref:`boosted frame method <theory-boostedframe>` to effectively model long acceleration lengths.

Alternatively, another common approximation for PWFAs is quasi-static modeling, useful if effects such as self-injection can be ignored.
In the Beam, Plasma & Accelerator Simulation Toolkit (BLAST), `HiPACE++ <https://hipace.readthedocs.io>`__ provides such methods.

.. note::

    There are two FDTD Maxwell field solvers that compute the field push implemented in WarpX: the Yee and Cole-Karkkainen solver with Cowan coefficients (CKC) solvers.
    The later includes a modification that allows the numerical dispersion of light in vacuum to be exact, and that is why we choose CKC for the example.

Lorentz boosted frame
---------------------

    WarpX simulations can be done in the laboratory or a :ref:`Lorentz-boosted <theory-boostedframe>` frame.
    In the laboratory frame, there is typically no need to model the plasma ions species, since they are mainly stationary during the short time scales associated with the motion of plasma electrons.
    In the boosted frame, that argument is no longer valid, as ions have relativistic velocities.

    To avoid spurious effects in the boosted frame, we want the longitudinal cell size to be larger than the transverse one.
    Translating this condition to the cell transverse (:math:`d_{x}`) and longitudinal dimensions (:math:`d_{z}`) in the laboratory frame leads to: :math:`d_{x} > (d_{z} (1+\beta_{b}) \gamma_{b})`,
    where :math:`\beta_{b}` is the boosted frame velocity in units of :math:`c`.

    To determine the total number of time steps of the simulation, we could either set the ``<zmax_plasma_to_compute_max_step>`` parameter to the end of the plasma (:math:`z_{\textrm{end}}`), or compute it using:

    * boosted frame edge of the simulation box, :math:`\textrm{corner} = l_{e}/ ((1-\beta_{b}) \gamma_{b})`
    * time of interaction in the boosted frame, :math:`T = \frac{z_{\textrm{end}}/\gamma_{b}-\textrm{corner}}{c (1+\beta_{b})}`
    * total number of iterations, :math:`i_{\textrm{max}} = T/dt`

    where :math:`l_{e}` is the position of the left edge of the simulation box (in respect to propagation direction).


Run
---

This example can be run **either** as:

* **Python** script: ``python3 PICMI_inputs_plasma_acceleration_3d.py`` or
* WarpX **executable** using an input file: ``warpx.3d inputs_3d_boost``

For `MPI-parallel <https://www.mpi-forum.org>`__ runs, prefix these lines with ``mpiexec -n 4 ...`` or ``srun -n 4 ...``, depending on the system.

.. tab-set::

   .. tab-item:: Python: Script

      .. literalinclude:: PICMI_inputs_plasma_acceleration_3d.py
         :language: python3
         :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/PICMI_inputs_plasma_acceleration.py``.

   .. tab-item:: Executable: Input File

      .. literalinclude:: inputs_3d_boost
         :language: ini
         :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/inputs_3d_boost``.

Visualize
---------

The plasma has a density profile that consists of a long up-ramp followed by a uniform plateau.

With this configuration the driver expels the plasma electrons, creating a bubble depleted of plasma electrons.
This bubble is the characteristic shape of the nonlinear 3D plasma wake.
This creates an accelerating field :math:`E_z`, as can be seen in :numref:`fig_pwfa_wake`, and transverse focusing forces (not shown).

.. _fig_pwfa_wake:

.. figure:: https://user-images.githubusercontent.com/10621396/292573839-8d6d7296-48f8-4942-94bd-548683fd35fe.png
   :alt: Electron beam, in black, accelerated by wakefield, in red-blue colormap, driven by drive electron beam, in purple. Plasma particles are plotted in green. Set max_steps to 800 in the input to generate data for this figure.
   :width: 100%

   Electron beam, in black, accelerated by wakefield, in red-blue colormap, driven by drive electron beam, in purple. Plasma particles are plotted in green. Set ``max_steps`` to 800 in the input to generate data for this figure. Note that the :math:`E_z` field is plotted with an interpolation option and appears smoother than the numerical resolution would suggest.

The beam increases in energy as it propagates in the wake, as shown in :numref:`fig_pwfa_energy`.
There is a slight inflection point in the evolution of the beam energy, showing how the beam accelerates as the wake forms around the electron beam and then gains energy linearly in time once the wake stabilizes.
Note that this is the mean energy in the boosted frame.
:ref:`Lab Frame <python_lab_frame_diag>`  or :ref:`back-transformed <running-cpp-parameters-diagnostics-btd>` diagnostics can be used to analyze beam energy in the lab frame. See also the :ref:`FAQ <faq-btd>`.

.. _fig_pwfa_energy:

.. figure:: https://user-images.githubusercontent.com/10621396/290962801-97d994f9-d48d-4f76-a37e-d14f6781d680.png
   :alt: Energy gain of electron beam in wake through this simulation (this is the mean beam energy in the boosted frame, NOT the lab frame). Set ``max_steps`` to 800 in the input to generate data for this figure.
   :width: 80%

   Energy gain of electron beam in wake through this simulation (this is the mean beam energy in the boosted frame, NOT the lab frame). Set ``max_steps`` to 800 in the input to generate data for this figure.

The plots can be generated with the following script:

.. dropdown:: Script ``plot_3d.py``

   .. literalinclude:: plot_3d.py
      :language: python3
      :caption: You can copy this file from ``Examples/Physics_applications/plasma_acceleration/plot_3d.py``.
