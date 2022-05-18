.. _glossary:

Glossary
========

In daily communication, we tend to abbreviate a lot of terms.
It is important to us to make it easy to interact with the WarpX community and thus, this list shall help to clarify often used terms.

Abbreviations
-------------

* **2FA:** `Two-factor-authentication <https://en.wikipedia.org/wiki/Multi-factor_authentication>`__
* **ABLASTR:** Accelerated BLAST Recipes, the library inside WarpX to share functionality with other BLAST codes
* **ALCF:** `Argonne Leadership Computing Facility <https://www.alcf.anl.gov/>`__, a supercomputing center located near Chicago, IL (USA)
* **BLAST:** `Beam, Plasma & Accelerator Simulation Toolkit <https://blast.lbl.gov>`__
* **AMR:** adaptive mesh-refinement
* **BC:** boundary condition (of a simulation)
* **BCK:** `Benkler-Chavannes-Kuster <https://ieeexplore.ieee.org/document/1638381>`__ method, a stabilization technique for small cells in the electromagnetic solver
* **BTD:** backtransformed diagnosics, a method to collect data for analysis from a *boosted frame* simulation
* **CFL:** the Courant-Friedrichs-Lewy condition, a numerical parameter for the numerical convergence of PDE solvers
* **CI:** continuous integration, automated tests that we perform before a proposed code-change is accepted; see PR
* **CPU:** `central processing unit <https://en.wikipedia.org/wiki/Central_processing_unit>`__, we usual mean a socket or generally the host-side of a computer (compared to the accelerator, e.g. GPU)
* **DOE:** `The United States Department of Energy <https://en.wikipedia.org/wiki/United_States_Department_of_Energy>`__, the largest sponsor of national laboratory research in the United States of America
* **ECP:** `Exascale Computing Project <https://www.exascaleproject.org>`__, a :ref:`U.S. DOE funding source <funding>` that supports WarpX development
* **ECT:** `Enlarged Cell Technique <https://ieeexplore.ieee.org/document/4463918>`__, an electromagnetic solver with accurate resolution of perfectly conducting embedded boundaries
* **EB:** embedded boundary, boundary conditions inside the simulation box, e.g. following material surfaces
* **EM:** electromagnetic, e.g. EM PIC
* **ES:** electrostatic, e.g. ES PIC
* **FDTD:** `Finite-difference time-domain or Yee's method <https://en.wikipedia.org/wiki/Finite-difference_time-domain_method>`__, a class of grid-based finite-difference field solvers
* **GPU:** originally graphics processing unit, now used for fast `general purpose computing (GPGPU) <https://en.wikipedia.org/wiki/Graphics_processing_unit#Stream_processing_and_general_purpose_GPUs_(GPGPU)>`__; also called (hardware) accelerator
* **LDRD:** Laboratory Directed Research and Development, a :ref:`funding program in U.S. DOE laboratories <funding>` that kick-started ABLASTR development
* **LPA:** laser-plasma acceleration, historically used for laser-electron acceleration
* **LPI:** laser-plasma interaction
* **LWFA:** laser-wakefield acceleration (of electrons/leptons)
* **MR:** mesh-refinement
* **MVA:** magnetic-vortex acceleration (of protons/ions)
* **NERSC:** `National Energy Research Scientific Computing Center <https://www.nersc.gov/>`__, a supercomputing center located in Berkeley, CA (USA)
* **NSF:** the `National Science Foundation <https://en.wikipedia.org/wiki/National_Science_Foundation>`__, a large public agency in the United States of America, supporting research and education
* **OLCF:** `Oak Ridge Leadership Computing Facility <https://www.olcf.ornl.gov/>`__, a supercomputing center located in Oak Ridge, TN (USA)
* **OTP:** `One-Time-Password <https://en.wikipedia.org/wiki/One-time_password>`__; see 2FA
* **PDE:** `partial differential equation <https://en.wikipedia.org/wiki/Partial_differential_equation>`__, an equation which imposes relations between the various partial derivatives of a multivariable function
* **PIC:** :ref:`particle-in-cell <theory-pic>`, the method implemented in WarpX
* **PR:** github pull request, a proposed change to the WarpX code base
* **PSATD:** pseudo-spectral analytical time-domain method, a spectral field solver with better numerical properties than FDTD solvers
* **PWFA:** plasma-wakefield acceleration
* **QED:** `quantum electrodynamics <https://en.wikipedia.org/wiki/Quantum_electrodynamics>`__
* **RPA:** radiation-pressure acceleration (of protons/ions), e.g. hole-boring (HB) or light-sail (LS) acceleration
* **RZ:** for the coordinate system ``r-z`` in cylindrical geometry; we use "RZ" when we refer to quasi-cylindrical geometry, decomposed in azimuthal modes (see details `here <https://fbpic.github.io/overview/pic_algorithm.html#cylindrical-grid-with-azimuthal-decomposition>`__)
* **TNSA:** target-normal sheet acceleration (of protons/ions)

Terms
-----

* **accelerator:** depending on context, either a *particle accelerator* in physics or a *hardware accelerator* (e.g. GPU) in computing
* **AMReX:** `C++ library for block-structured adaptive mesh-refinement <https://amrex-codes.github.io/>`__, a primary dependency of WarpX
* **boosted frame:** a :ref:`Lorentz-boosted frame of reference <theory-boostedframe>` for a simulation
* **evolve:** this is a generic term to advance a quantity (same nomenclature in AMReX).
              For instance, ``WarpX::EvolveE(dt)`` advances the electric field for duration ``dt``, ``PhysicalParticleContainer::Evolve(...)`` does field gather + particle push + current deposition for all particles in ``PhysicalParticleContainer``, and ``WarpX::EvolveEM`` is the central ``WarpX`` function that performs 1 PIC iteration.
* **Frontier:** an `Exascale supercomputer at OLCF <https://www.olcf.ornl.gov/frontier/>`__
* **laser:** most of the time, we mean a `laser pulse <https://en.wikipedia.org/wiki/Ultrashort_pulse>`__
* **openPMD**: `Open Standard for Particle-Mesh Data Files <https://www.openPMD.org>`__, a community meta-data project for scientific data
* **Perlmutter:** a Berkeley Lab nobel laureate and a `Pre-Exascale supercomputer at NERSC <https://www.nersc.gov/systems/perlmutter/>`__
* **plotfiles**: the internal binary format for data files in *AMReX*
* **Python:** a popular scripted `programming language <https://www.python.org>`__
