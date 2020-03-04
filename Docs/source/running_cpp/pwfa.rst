In-depth explanation of a PWFA simulation
=========================================

As described in the Introduction
`Introduction <https://warpx.readthedocs.io/en/latest/theory/intro.html>`_, one
of the key applications of the WarpX exascale computing platform is in modelling
future, compact and economic plasma-based accelerators. In this section we
describe the simulation setup of a realistic electron beam driven plasma
wakefield accelerator (PWFA) configuration. For illustration porpuses the setup
can be explored with **WarpX** using the example input file :download:`PWFA
<../../../Examples/Physics_applications/plasma_acceleration/inputs_2d_boost>`.

The simulation setup consists of 4 particle species: drive
beam (driver), witness beam (beam), plasma electrons (plasma_e), and plasma
ions (plasma_p). The species parameters are summarized in the following table.

======== ===============================================
Species  Parameters
======== ===============================================
driver   γ = 48923; N = 2x10^8; σz = 4.0 μm; σx = 2.0 μm
beam     γ = 48923; N = 6x10^5; σz = 1.0 mm; σx = 0.5 μm
plasma_e n = 1x10^23 m^-3
plasma_p n = 1x10^23 m^-3
======== ===============================================

Where γ is the beam relativisitc Lorentz factor, N is the number of particles,
and σx, σy, σz are the beam widths (root-mean-squares of particle positions) in
the transverse (x,y) and longitudinal directions.

Listed below are the key arguments and best-practices relevant for chosing the
pwfa simulation parameters used in the example.


2D Geometry
-----------

    2D cartesian (with longitudinal direction z and transverse x) geometry
    simulaions can give valuable physical and numerical insight into the
    simulation requirements and evolution while being much less time consuming
    than the full 3D cartesian or cylindrical geometries.


Lorentz boosted frame
---------------------

    WarpX simulations can be done in the laboratory or `Lorentz-boosted
    <https://warpx.readthedocs.io/en/latest/theory/boosted_frame.html>`_ frames.
    In the laboratory frame, there is typically no need to model the palsma ions
    species, since they are mainly stationary during the short time scales
    associated with the motion of plasma electrons. In the boosted frame, that
    argument is no longer valid, as ions have relativistic velocities. The
    boosted frame still results in a substantial reduction to the simulation
    computational cost.

.. note::
Regardless of the frame that is chosen for the simulation, the
input file parameters are defined in respect to the laboratory frame.


Resolution
----------

    Longitudinal and transverse resolutions should be chosen to accurately
    describe the physical processes taking place in the simulation. Convergence
    scans, where resolution in both directions is gradually increased, should be
    used to determine the optimal configuration. Multiple cells per beam length
    and width are recommended (our illustrative example resolution is coarse).

.. note::
    To avoid spurious effects, in the boosted frame, we consider the contrain
    that the transverse cell size (dx) should be larger than the transverse one
    (dz), i.e. dx > (dz*2*γb).


Time step
---------

    The time step is computed differently depending on the field solver, being
    computed  (if we
    would use Yee instead of CKC)


(γb=10)


The example we use employs the
boosted frame, with relativistic factor γb=10, to substantially reduce the
computational costs of the simulation and that is why it includes the definition
of plasma_p.

The separation between the driver and witness beams is set to 125 μm so that the
witness beam falls in the accelerating and focusing region of the plasma bubble.




Longitudinal and transverse resolutions should be chosen to accurately describe
the physical processes taking place in the simulation. To avoid spurious effects
 we consider the contrain that, in the boosted frame, the transverse cell size
 (dx=dy) should be larger than the transverse one (dz), i.e. dx > (dz*2*γb).

* the simulation box resolution ratio in the boosted frame should be lower than
1. In other words








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

