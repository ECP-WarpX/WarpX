In-depth explanation of a PWFA simulation
=========================================

As described in the :doc:`intro`, one of the key applications of the WarpX exascale computing platform is in modelling future, compact and economic plasma-based accelerators.
In this section we describe the simulation setup of a realistic electron beam driven plasma wakefield accelerator (PWFA) configuration.
For illustration porpuses the setup can be explored with **WarpX** using the example input file :download:`PWFA <../../../Examples/Physics_applications/plasma_acceleration/inputs_2d_boost>`.

The simulation setup consists of 4 particle species: drive beam (driver), witness beam (beam), plasma electrons (plasma_e), and plasma ions (plasma_p).
The species physical parameters are summarized in the following table.

======== ============================================================
Species  Parameters
======== ============================================================
driver   :math:`\gamma` = 48923; N = 2x10^8; σz = 4.0 μm; σx = 2.0 μm
beam     :math:`\gamma` = 48923; N = 6x10^5; σz = 1.0 mm; σx = 0.5 μm
plasma_e n = 1x10^23 m^-3; w = 70 μm; lr = 8 mm; L = 200 mm
plasma_p n = 1x10^23 m^-3; w = 70 μm; lr = 8 mm; L = 200 mm
======== ============================================================

Where :math:`\gamma` is the beam relativisitc Lorentz factor, N is the number of particles, and σx, σy, σz are the beam widths (root-mean-squares of particle positions) in the transverse (x,y) and longitudinal directions.
The plasma, of total lenght L, has a density profile that consists of a lr long linear up-ramp, ranging from 0 to peak value n, is uniform within a transverse width of w and after the up-ramp.

Listed below are the key arguments and best-practices relevant for chosing the pwfa simulation parameters used in the example.


2D Geometry
-----------

    2D cartesian (with longitudinal direction z and transverse x) geometry simulaions can give valuable physical and numerical insight into the simulation requirements and evolution.
    At the same time it is much less time consuming than the full 3D cartesian or cylindrical geometries.


Finite Difference Time Domain
-----------------------------

    For standard plasma wakefield configurations, it is possible to model the physics correctly using the Particle-In-Cell (PIC) Finite Difference Time Domain (FDTD) algorithms (:doc:`picsar_theory`).
    If the simulation contains localised extremely high intensity fields, however, numerical instabilities might arise, such as the numerical Cherenkov instability (:doc:`boosted_frame`).
    In that case, it is recommended to use the Pseudo Spectral Analytical Time Domain (PSATD) or the Pseudo-Spectral Time-Domain (PSTD) algorithms.
    In the example we are describing, it is sufficient to use FDTD.
    

Cole-Karkkainen solver with Cowan coefficients
----------------------------------------------

    There are two FDTD Maxwell field solvers that compute the field push implemented in WarpX: the Yee and Cole-Karkkainen solver with Cowan coefficients (CKC) solvers.
    The later includes a modifitication that allows the numerical dispersion of light in vaccum to be exact, and that is why we choose CKC for the example.
    

Lorentz boosted frame
---------------------

    WarpX simulations can be done in the laboratory or `Lorentz-boosted <https://warpx.readthedocs.io/en/latest/theory/boosted_frame.html>`_ frames.
    In the laboratory frame, there is typically no need to model the palsma ions species, since they are mainly stationary during the short time scales associated with the motion of plasma electrons.
    In the boosted frame, that argument is no longer valid, as ions have relativistic velocities.
    The boosted frame still results in a substantial reduction to the simulation computational cost.

.. note::
   Regardless of the frame that is chosen for the simulation, the input file parameters are defined in respect to the laboratory frame.


Resolution
----------

    Longitudinal and transverse resolutions (i.e. number and dimensions of the PIC grid cells) should be chosen to accurately describe the physical processes taking place in the simulation.
    Convergence scans, where resolution in both directions is gradually increased, should be used to determine the optimal configuration.
    Multiple cells per beam length and width are recommended (our illustrative example resolution is coarse).

.. note::
    To avoid spurious effects, in the boosted frame, we consider the contrain that the transverse cell size should be larger than the transverse one.
    Traslating this condition to the cell transverse (:math:`d_{x}`) and longitudinal dimensions (:math:`d_{z}`) in the laboratory frame leads to:
    .. math:: d_{x} > (d_{z}*(1+\beta_{b})*\gamma_{b})
    where :math:`\beta_{b}` is the boosted frame velocity in units of speed-of-ligth (:math:`c`).


Time step
---------

    The time step (:math:`dt`) is used to iterated over the main PIC loop and is computed by WarpX differently depending on the Maxwell field FDTD solvers used:
    
    * **For Yee** is equal to the CFL parameter chosen in the input file (:doc:`parameters`) times the Courant–Friedrichs–Lewy condition (CFL) that follows the analytical expression in :doc:`picsar_theory`
    * **For CKC:** is equal to CFL times the minimum between the boosted frame cell dimensions

    where CFL is choosen to be below unity and set an optimal trade-off between making the simulation faster and avoiding NCI and other spurious effects.


Duration of the simulation
--------------------------

    To determine the total number of time steps of the simulation, we could either set the `<zmax_plasma_to_compute_max_step>` parameter to the end of the plasma (:math:`z_{end}`), or compute it using:
    * boosted frame edge of the simulation box - :math:`\textrm{corner} = l_{e}/ ((1-\beta_{b}) \gamma_{b})`
    * time of interaction in the boosted frame - :math:`T = \frac{z_end/\gamma_{b}-\textrm{corner}}{c (1+\beta_{b})}`
    * total number of iterations - :math:`i_{\textrm{max}} = T/dt`
    
       





In the example input, all the simulation parameters are defined in the lab frame
regardless of if the boosted frame is used. This is true also for the
longitudinal and transverse dimensions of the simulation box, the diagnostic
time snapshots for back-transformed data to the lab frame from a boosted-frame
simulation. Thus, when one defines the grid size in the lab frame, the
longitudinal resolution remains the same, but the transverse grid sizes need to
be adjusted approximately in the boosted frame with the following relation

The time step in the boosted frame is increased as

Here γ is the Lorentz factor of the boosted frame. In the boosted frame with β close to 1 in the forward direction of the beam propagation, the beam length and plasma length change, respectively, according to

Define the total run time of a simulation by the full transit time of the beam through the plasma, and they are given by, respectively in the lab and boosted frame



assuming the plasma moving at c opposite to the beam direction. Thus the number of time steps in the lab and boosted frame are

It should be pointed out that this example is performed in 2D x-y geometry, which is not equivalent to the realistic simulation. However, the fast turnaround time in 2D simulation helps determine the numerical requirements and the optimized boosted frame, which can then be used in 3D simulations.

