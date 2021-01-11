.. _running-cpp-parameters:

Input parameters
================

.. warning::

   This section is currently in development.

.. note::
   The WarpXParser (see :ref:`running-cpp-parameters-parser`) is used for the right-hand-side of all input parameters that consist in a single real number, so expressions like ``<species_name>.density_max = "2.+1."`` and/or using user-defined constants are accepted. See below for more detail.

.. _running-cpp-parameters-overall:

Overall simulation parameters
-----------------------------

* ``authors`` (`string`: e.g. ``"Jane Doe <jane@example.com>, Jimmy Joe <jimmy@example.com>"``)
    Authors of an input file / simulation setup.
    When provided, this information is added as metadata to (openPMD) output files.

* ``max_step`` (`integer`)
    The number of PIC cycles to perform.

* ``warpx.gamma_boost`` (`float`)
    The Lorentz factor of the boosted frame in which the simulation is run.
    (The corresponding Lorentz transformation is assumed to be along ``warpx.boost_direction``.)

    When using this parameter, some of the input parameters are automatically
    converted to the boosted frame. (See the corresponding documentation of each
    input parameters.)

    .. note::

        For now, only the laser parameters will be converted.

* ``warpx.boost_direction`` (string: ``x``, ``y`` or ``z``)
    The direction of the Lorentz-transform for boosted-frame simulations
    (The direction ``y`` cannot be used in 2D simulations.)

* ``warpx.zmax_plasma_to_compute_max_step`` (`float`) optional
    Can be useful when running in a boosted frame. If specified, automatically
    calculates the number of iterations required in the boosted frame for the
    lower `z` end of the simulation domain to reach
    ``warpx.zmax_plasma_to_compute_max_step`` (typically the plasma end,
    given in the lab frame). The value of ``max_step`` is overwritten, and
    printed to standard output. Currently only works if the Lorentz boost and
    the moving window are along the z direction.

* ``warpx.verbose`` (``0`` or ``1``; default is ``1`` for true)
    Controls how much information is printed to the terminal, when running WarpX.

* ``warpx.random_seed`` (`string` or `int` > 0) optional
    If provided ``warpx.random_seed = random``, the random seed will be determined
    using `std::random_device` and `std::clock()`,
    thus every simulation run produces different random numbers.
    If provided ``warpx.random_seed = n``, and it is required that `n > 0`,
    the random seed for each MPI rank is `(mpi_rank+1) * n`,
    where `mpi_rank` starts from 0.
    `n = 1` and ``warpx.random_seed = default``
    produce the default random seed.
    Note that when GPU threading is used,
    one should not expect to obtain the same random numbers,
    even if a fixed ``warpx.random_seed`` is provided.

* ``warpx.do_electrostatic`` (`string`) optional (default `none`)
    Specifies the electrostatic mode. When turned on, instead of updating
    the fields at each iteration with the full Maxwell equations, the fields
    are recomputed at each iteration from the Poisson equation.
    There is no limitation on the timestep in this case, but
    electromagnetic effects (e.g. propagation of radiation, lasers, etc.)
    are not captured. There are two options:

    * ``labframe``: Poisson's equation is solved in the lab frame with
      the charge density of all species combined. There will only be E
      fields.

    * ``relativistic``: Poisson's equation is solved for each species
      seperately taking into account their averaged velocities. The field
      is mapped to the simulation frame and will produce both E and B
      fields.

* ``self_fields_required_precision`` (`float`, default: 1.e-11)
    The relative precision with which the electrostatic space-charge fields should
    be calculated. More specifically, the space-charge fields are
    computed with an iterative Multi-Level Multi-Grid (MLMG) solver.
    This solver can fail to reach the default precision within a reasonable
    This only applies when warpx.do_electrostatic = labframe.

* ``self_fields_max_iters`` (`integer`, default: 200)
    Maximum number of iterations used for MLMG solver for space-charge
    fields calculation. In case if MLMG converges but fails to reach the desired
    ``self_fields_required_precision``, this parameter may be increased.
    This only applies when warpx.do_electrostatic = labframe.

* ``amrex.abort_on_out_of_gpu_memory``  (``0`` or ``1``; default is ``1`` for true)
    When running on GPUs, memory that does not fit on the device will be automatically swapped to host memory when this option is set to ``0``.
    This will cause severe performance drops.
    Note that even with this set to ``1`` WarpX will not catch all out-of-memory events yet when operating close to maximum device memory.
    `Please also see the documentation in AMReX <https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters>`_.

.. _running-cpp-parameters-box:

Setting up the field mesh
-------------------------

* ``amr.n_cell`` (`2 integers in 2D`, `3 integers in 3D`)
    The number of grid points along each direction (on the **coarsest level**)

* ``amr.max_level`` (`integer`, default: ``0``)
    When using mesh refinement, the number of refinement levels that will be used.

    Use 0 in order to disable mesh refinement.
    Note: currently, ``0`` and ``1`` are supported.

* ``amr.ref_ratio`` (`integer` per refined level, default: ``2``)
    When using mesh refinement, this is the refinement ratio per level.
    With this option, all directions are fined by the same ratio.

    Note: in development; currently, ``2`` is supported.

* ``amr.ref_ratio_vect`` (3 `integer`s for x,y,z per refined level)
    When using mesh refinement, this can be used to set the refinement ratio per direction and level, relative to the previous level.

    Example: for three levels, a value of ``2 2 4 8 8 16`` refines the first level by 2-fold in x and y and 4-fold in z compared to the coarsest level (level 0/mother grid); compared to the first level, the second level is refined 8-fold in x and y and 16-fold in z.

    Note: in development; currently allowed value: ``2 2 2``.

* ``geometry.is_periodic`` (`2 integers in 2D`, `3 integers in 3D`)
    Whether the boundary conditions are periodic, in each direction.

    For each direction, use 1 for periodic conditions, 0 otherwise.

* ``geometry.coord_sys`` (`integer`) optional (default `0`)
    Coordinate system used by the simulation. 0 for Cartesian, 1 for cylindrical.

* ``geometry.prob_lo`` and ``geometry.prob_hi`` (`2 floats in 2D`, `3 integers in 3D`; in meters)
    The extent of the full simulation box. This box is rectangular, and thus its
    extent is given here by the coordinates of the lower corner (``geometry.prob_lo``) and
    upper corner (``geometry.prob_hi``). The first axis of the coordinates is x (or r with cylindrical)
    and the last is z.

* ``warpx.fine_tag_lo`` and ``warpx.fine_tag_hi`` (`2 floats in 2D`, `3 integers in 3D`; in meters) optional
    **When using static mesh refinement with 1 level**, the extent of the refined patch.
    This patch is rectangular, and thus its extent is given here by the coordinates
    of the lower corner (``warpx.fine_tag_lo``) and upper corner (``warpx.fine_tag_hi``).

* ``warpx.n_current_deposition_buffer`` (`integer`)
    When using mesh refinement: the particles that are located inside
    a refinement patch, but within ``n_current_deposition_buffer`` cells of
    the edge of this patch, will deposit their charge and current to the
    lower refinement level, instead of depositing to the refinement patch
    itself. See the section :doc:`../../theory/amr` for more details.
    If this variable is not explicitly set in the input script,
    ``n_current_deposition_buffer`` is automatically set so as to be large
    enough to hold the particle shape, on the fine grid

* ``warpx.n_field_gather_buffer`` (`integer`; 0 by default)
    When using mesh refinement: the particles that are located inside
    a refinement patch, but within ``n_field_gather_buffer`` cells of
    the edge of this patch, will gather the fields from the lower refinement
    level, instead of gathering the fields from the refinement patch itself.
    This avoids some of the spurious effects that can occur inside the
    refinement patch, close to its edge. See the section
    :doc:`../../theory/amr` for more details. If this variable is not
    explicitly set in the input script, ``n_field_gather_buffer`` is
    automatically set so that it is one cell larger than
    ``n_current_deposition_buffer``, on the fine grid.

* ``particles.deposit_on_main_grid`` (`list of strings`)
    When using mesh refinement: the particle species whose name are included
    in the list will deposit their charge/current directly on the main grid
    (i.e. the coarsest level), even if they are inside a refinement patch.

* ``particles.gather_from_main_grid`` (`list of strings`)
    When using mesh refinement: the particle species whose name are included
    in the list will gather their fields from the main grid
    (i.e. the coarsest level), even if they are inside a refinement patch.

* ``warpx.n_rz_azimuthal_modes`` (`integer`; 1 by default)
    When using the RZ version, this is the number of azimuthal modes.

.. _running-cpp-parameters-parallelization:

Distribution across MPI ranks and parallelization
-------------------------------------------------

* ``warpx.numprocs`` (`2 ints` for 2D, `3 ints` for 3D) optional (default `none`)
    This optional parameter can be used to control the domain decomposition on the
    coarsest level. The domain will be chopped into the exact number of pieces in each
    dimension as specified by this parameter. If it's not specified, the domain
    decomposition will be determined by the parameters that will be discussed below.  If
    specified, the product of the numbers must be equal to the number of MPI processes.

* ``amr.max_grid_size`` (`integer`) optional (default `128`)
    Maximum allowable size of each **subdomain**
    (expressed in number of grid points, in each direction).
    Each subdomain has its own ghost cells, and can be handled by a
    different MPI rank ; several OpenMP threads can work simultaneously on the
    same subdomain.

    If ``max_grid_size`` is such that the total number of subdomains is
    **larger** that the number of MPI ranks used, than some MPI ranks
    will handle several subdomains, thereby providing additional flexibility
    for **load balancing**.

    When using mesh refinement, this number applies to the subdomains
    of the coarsest level, but also to any of the finer level.

* ``warpx.load_balance_int`` (`string`) optional (default `0`)
    Using the `Intervals parser`_ syntax, this string defines the timesteps at which
    WarpX should try to redistribute the work across MPI ranks, in order to have
    better load balancing.
    Use 0 to disable load_balancing.

    When performing load balancing, WarpX measures the wall time for
    computational parts of the PIC cycle. It then uses this data to decide
    how to redistribute the subdomains across MPI ranks. (Each subdomain
    is unchanged, but its owner is changed in order to have better performance.)
    This relies on each MPI rank handling several (in fact many) subdomains
    (see ``max_grid_size``).

* ``warpx.load_balance_with_sfc`` (`0` or `1`) optional (default `0`)
    If this is `1`: use a Space-Filling Curve (SFC) algorithm in order to
    perform load-balancing of the simulation.
    If this is `0`: the Knapsack algorithm is used instead.

* ``warpx.load_balance_efficiency_ratio_threshold`` (`float`) optional (default `1.1`)
    Controls whether to adopt a proposed distribution mapping computed during a load balance.
    If the the ratio of the proposed to current distribution mapping *efficiency* (i.e.,
    average cost per MPI process; efficiency is a number in the range [0, 1]) is greater
    than the threshold value, the proposed distribution mapping is adopted.  The suggested
    range of values is ``warpx.load_balance_efficiency_ratio_threshold >= 1``, which ensures
    that the new distribution mapping is adopted only if doing so would improve the load
    balance efficiency. The higher the threshold value, the more conservative is the criterion
    for adoption of a proposed distribution; for example, with
    ``warpx.load_balance_efficiency_ratio_threshold = 1``, the proposed distribution is
    adopted *any* time the proposed distribution improves load balancing; if instead
    ``warpx.load_balance_efficiency_ratio_threshold = 2``, the proposed distribution is
    adopted only if doing so would yield a 100% to the load balance efficiency (with this
    threshold value, if the  current efficiency is ``0.45``, the new distribution would only be
    adopted if the proposed efficiency were greater than ``0.9``).

* ``algo.load_balance_costs_update`` (`Heuristic` or `Timers`) optional (default `Timers`)
    If this is `Heuristic`: load balance costs are updated according to a measure of
    particles and cells assigned to each box of the domain.  The cost :math:`c` is
    computed as

    .. math::

            c = n_{\text{particle}} \cdot w_{\text{particle}} + n_{\text{cell}} \cdot w_{\text{cell}},

    where
    :math:`n_{\text{particle}}` is the number of particles on the box,
    :math:`w_{\text{particle}}` is the particle cost weight factor (controlled by ``algo.costs_heuristic_particles_wt``),
    :math:`n_{\text{cell}}` is the number of cells on the box, and
    :math:`w_{\text{cell}}` is the cell cost weight factor (controlled by ``algo.costs_heuristic_cells_wt``).

    If this is `Timers`: costs are updated according to in-code timers.

* ``algo.costs_heuristic_particles_wt`` (`float`) optional
    Particle weight factor used in `Heuristic` strategy for costs update; if running on GPU,
    the particle weight is set to a value determined from single-GPU tests on Summit,
    depending on the choice of solver (FDTD or PSATD) and order of the particle shape.
    If running on CPU, the default value is `0.9`.

* ``algo.costs_heuristic_cells_wt`` (`float`) optional
    Cell weight factor used in `Heuristic` strategy for costs update; if running on GPU,
    the cell weight is set to a value determined from single-GPU tests on Summit,
    depending on the choice of solver (FDTD or PSATD) and order of the particle shape.
    If running on CPU, the default value is `0.1`.

* ``warpx.do_dynamic_scheduling`` (`0` or `1`) optional (default `1`)
    Whether to activate OpenMP dynamic scheduling.

* ``warpx.safe_guard_cells`` (`0` or `1`) optional (default `0`)
    For developers: run in safe mode, exchanging more guard cells, and more often in the PIC loop (for debugging).

.. _running-cpp-parameters-parser:

Math parser and user-defined constants
--------------------------------------

WarpX provides a math parser that reads expressions in the input file.
It can be used in all input parameters that consist in one real number.

WarpX constants
###############

WarpX provides a few pre-defined constants, that can be used for any parameter that consists in one real number.

======== ===================
q_e      elementary charge
m_e      electron mass
m_p      proton mass
epsilon0 vacuum permittivity
clight   speed of light
pi       math constant pi
======== ===================

See ``Source/Utils/WarpXConst.H`` for the values.

User-defined constants
######################

Users can define their own constants in the input file.
These constants can be used for any parameter that consists in one real number.
User-defined constants can contain only letters, numbers and the character ``_``.
The name of each constant has to begin with a letter. The following names are used
by WarpX, and cannot be used as user-defined constants: ``x``, ``y``, ``z``, ``X``, ``Y``, ``t``.
For example, parameters ``a0`` and ``z_plateau`` can be specified with:

* ``my_constants.a0 = 3.0``
* ``my_constants.z_plateau = 150.e-6``

Coordinates
###########

Besides, for profiles that depend on spatial coordinates (the plasma momentum distribution or the laser field, see below `Particle initialization` and `Laser initialization`), the parser will interpret some variables as spatial coordinates. These are specified in the input parameter, i.e., ``density_function(x,y,z)`` and ``field_function(X,Y,t)``.

The parser reads python-style expressions between double quotes, for instance
``"a0*x**2 * (1-y*1.e2) * (x>0)"`` is a valid expression where ``a0`` is a
user-defined constant (see below) and ``x`` and ``y`` are spatial coordinates. The names are case sensitive. The factor
``(x>0)`` is ``1`` where ``x>0`` and ``0`` where ``x<=0``. It allows the user to
define functions by intervals.
The parser reads mathematical functions into an `abstract syntax tree (AST) <https://en.wikipedia.org/wiki/Abstract_syntax_tree>`_, which supports a maximum depth (see :ref:`build options <building-cmake>`_).
Additional terms in a function can create a level of depth in the AST, e.g. ``a+b+c+d`` is parsed in groups of ``[+ a [+ b [+ c [+ d]]]]`` (depth: 4).
A trick to reduce this depth for the parser, e.g. when reaching the limit, is to group expliclity, e.g. via ``(a+b)+(c+d)``, which is parsed in groups of ``[+ [+ a b] [+ c d]]`` (depth: 2).

.. _running-cpp-parameters-particle:

Particle initialization
-----------------------

* ``particles.species_names`` (`strings`, separated by spaces)
    The name of each species. This is then used in the rest of the input deck ;
    in this documentation we use `<species_name>` as a placeholder.

* ``particles.use_fdtd_nci_corr`` (`0` or `1`) optional (default `0`)
    Whether to activate the FDTD Numerical Cherenkov Instability corrector.
    Not currently available in the RZ configuration.

* ``particles.boundary_conditions`` (`string`) optional (default `none`)
    Boundary conditions applied to particles. Options are:
    * ``none``: the boundary conditions applied to particles is determined by ``geometry.is_periodic``.
    * ``absorbing``: particles exiting the simulation domain are discarded.

* ``particles.rigid_injected_species`` (`strings`, separated by spaces)
    List of species injected using the rigid injection method. The rigid injection
    method is useful when injecting a relativistic particle beam, in boosted-frame
    simulation ; see the section :doc:`../../theory/input_output` for more details.
    For species injected using this method, particles are translated along the `+z`
    axis with constant velocity as long as their ``z`` coordinate verifies
    ``z<zinject_plane``. When ``z>zinject_plane``,
    particles are pushed in a standard way, using the specified pusher.
    (see the parameter ``<species_name>.zinject_plane`` below)

* ``<species_name>.species_type`` (`string`) optional (default `unspecified`)
    Type of physical species, ``"electron"``, ``"positron"``, ``"photon"``, ``"hydrogen"``.
    Either this or both ``mass`` and ``charge`` have to be specified.

* ``<species_name>.charge`` (`float`) optional (default `NaN`)
    The charge of one `physical` particle of this species.
    If ``species_type`` is specified, the charge will be set to the physical value and ``charge`` is optional.
    When ``<species>.do_field_ionization = 1``, the physical particle charge is equal to ``ionization_initial_level * charge``, so latter parameter should be equal to q_e (which is defined in WarpX as the elementary charge in coulombs).

* ``<species_name>.mass`` (`float`) optional (default `NaN`)
    The mass of one `physical` particle of this species.
    If ``species_type`` is specified, the mass will be set to the physical value and ``mass`` is optional.

* ``<species_name>.xmin,ymin,zmin`` and ``<species_name>.xmax,ymax,zmax`` (`float`) optional (default unlimited)
    When ``<species_name>.xmin`` and ``<species_name>.xmax`` are set, they delimit the region within which particles are injected.
    If periodic boundary conditions are used in direction ``i``, then the default (i.e. if the range is not specified) range will be the simulation box, ``[geometry.prob_hi[i], geometry.prob_lo[i]]``.

* ``<species_name>.injection_style`` (`string`)
    Determines how the particles will be injected in the simulation.
    The options are:

    * ``NUniformPerCell``: injection with a fixed number of evenly-spaced particles per cell.
      This requires the additional parameter ``<species_name>.num_particles_per_cell_each_dim``.

    * ``NRandomPerCell``: injection with a fixed number of randomly-distributed particles per cell.
      This requires the additional parameter ``<species_name>.num_particles_per_cell``.

    * ``SingleParticle``: Inject a single macroparticle.
      This requires the additional parameters:
      ``<species_name>.single_particle_pos`` (`3 doubles`, particle 3D position [meter])
      ``<species_name>.single_particle_vel`` (`3 doubles`, particle 3D normalized momentum, i.e. :math:`\gamma \beta`)
      ``<species_name>.single_particle_weight`` ( `double`, macroparticle weight, i.e. number of physical particles it represents)

    * ``gaussian_beam``: Inject particle beam with gaussian distribution in
      space in all directions. This requires additional parameters:
      ``<species_name>.q_tot`` (beam charge) optional (default is ``q_tot=0``),
      ``<species_name>.npart`` (number of particles in the beam),
      ``<species_name>.x/y/z_m`` (average position in `x/y/z`),
      ``<species_name>.x/y/z_rms`` (standard deviation in `x/y/z`),
      ``<species_name>.x/y/z_rms`` (standard deviation in `x/y/z`),
      ``<species_name>.x/y/z_cut`` (optional, particles with ``abs(x-x_m) > x_cut*x_rms`` are not injected, same for y and z. ``<species_name>.q_tot`` is the charge of the un-cut beam, so that cutting the distribution is likely to result in a lower total charge),
      and optional argument ``<species_name>.do_symmetrize`` (whether to
      symmetrize the beam in the x and y directions).

    * ``external_file``: Inject macroparticles with properties (mass, charge, position, and momentum - :math:`\gamma \beta m c`) read from an external openPMD file.
      With it users can specify the additional arguments:
      ``<species_name>.injection_file`` (`string`) openPMD file name and
      ``<species_name>.q_tot`` (`double`) optional (default is ``q_tot=0`` and no re-scaling is done, ``weight=q_p``) when specified it is used to re-scale the weight of externally loaded ``N`` physical particles, each of charge ``q_p``, to inject macroparticles of ``weight=<species_name>.q_tot/q_p/N``.
      ``<species_name>.charge`` (`double`) optional (default is read from openPMD file) when set this will be the charge of the physical particle represented by the injected macroparticles.
      ``<species_name>.mass`` (`double`) optional (default is read from openPMD file) when set this will be the charge of the physical particle represented by the injected macroparticles.
      ``<species_name>.z_shift`` (`double`) optional (default is no shift) when set this value will be added to the longitudinal, ``z``, position of the particles.
      The external file must include the species ``openPMD::Record``s labeled ``position`` and ``momentum`` (`double` arrays), with dimensionality and units set via ``openPMD::setUnitDimension`` and ``setUnitSI``.
      If the external file also contains ``openPMD::Records``s for ``mass`` and ``charge`` (constant `double` scalars) then the species will use these, unless overwritten in the input file (see ``<species_name>.mass``, ```<species_name>.charge`` or ```<species_name>.species_type``).
      The ``external_file`` option is currently implemented for 2D, 3D and RZ geometries, with record components in the cartesian coordinates ``(x,y,z)`` for 3D and RZ, and ``(x,z)`` for 2D.
      For more information on the `openPMD format <https://github.com/openPMD>`__ and how to build WarpX with it, please visit :doc:`../building/openpmd`.

* ``<species_name>.num_particles_per_cell_each_dim`` (`3 integers in 3D and RZ, 2 integers in 2D`)
    With the NUniformPerCell injection style, this specifies the number of particles along each axis
    within a cell. Note that for RZ, the three axis are radius, theta, and z and that the recommended
    number of particles per theta is at least two times the number of azimuthal modes requested.
    (It is recommended to do a convergence scan of the number of particles per theta)

* ``<species_name>.do_continuous_injection`` (`0` or `1`)
    Whether to inject particles during the simulation, and not only at
    initialization. This can be required with a moving window and/or when
    running in a boosted frame.

* ``<species_name>.initialize_self_fields`` (`0` or `1`)
    Whether to calculate the space-charge fields associated with this species
    at the beginning of the simulation.
    The fields are calculated for the mean gamma of the species.

* ``<species_name>.self_fields_required_precision`` (`float`, default: 1.e-11)
    The relative precision with which the initial space-charge fields should
    be calculated. More specifically, the initial space-charge fields are
    computed with an iterative Multi-Level Multi-Grid (MLMG) solver.
    For highly-relativistic beams, this solver can fail to reach the default
    precision within a reasonable time ; in that case, users can set a
    relaxed precision requirement through ``self_fields_required_precision``.

* ``<species_name>.self_fields_max_iters`` (`integer`, default: 200)
    Maximum number of iterations used for MLMG solver for initial space-charge
    fields calculation. In case if MLMG converges but fails to reach the desired
    ``self_fields_required_precision``, this parameter may be increased.

* ``<species_name>.profile`` (`string`)
    Density profile for this species. The options are:

    * ``constant``: Constant density profile within the box, or between ``<species_name>.xmin``
      and ``<species_name>.xmax`` (and same in all directions). This requires additional
      parameter ``<species_name>.density``. i.e., the plasma density in :math:`m^{-3}`.

    * ``parse_density_function``: the density is given by a function in the input file.
      It requires additional argument ``<species_name>.density_function(x,y,z)``, which is a
      mathematical expression for the density of the species, e.g.
      ``electrons.density_function(x,y,z) = "n0+n0*x**2*1.e12"`` where ``n0`` is a
      user-defined constant, see above. WARNING: where ``density_function(x,y,z)`` is close to zero, particles will still be injected between ``xmin`` and ``xmax`` etc., with a null weight. This is undesirable because it results in useless computing. To avoid this, see option ``density_min`` below.

* ``<species_name>.density_min`` (`float`) optional (default `0.`)
    Minimum plasma density. No particle is injected where the density is below this value.

* ``<species_name>.density_max`` (`float`) optional (default `infinity`)
    Maximum plasma density. The density at each point is the minimum between the value given in the profile, and `density_max`.

* ``<species_name>.radially_weighted`` (`bool`) optional (default `true`)
    Whether particle's weight is varied with their radius. This only applies to cylindrical geometry.
    The only valid value is true.

    * ``predefined``: use one of WarpX predefined plasma profiles. It requires additional
      arguments ``<species_name>.predefined_profile_name`` and
      ``<species_name>.predefined_profile_params`` (see below).

* ``<species_name>.momentum_distribution_type`` (`string`)
    Distribution of the normalized momentum (`u=p/mc`) for this species. The options are:

    * ``constant``: constant momentum profile. This requires additional parameters
      ``<species_name>.ux``, ``<species_name>.uy`` and ``<species_name>.uz``, the normalized
      momenta in the x, y and z direction respectively.

    * ``gaussian``: gaussian momentum distribution in all 3 directions. This requires
      additional arguments for the average momenta along each direction
      ``<species_name>.ux_m``, ``<species_name>.uy_m`` and ``<species_name>.uz_m`` as
      well as standard deviations along each direction ``<species_name>.ux_th``,
      ``<species_name>.uy_th`` and ``<species_name>.uz_th``.

    * ``maxwell_boltzmann``: Maxwell-Boltzmann distribution that takes a dimensionless
      temperature parameter ``<species_name>.theta`` as an input, where theta is kb*T/(m*c^2),
      kb is the Boltzmann constant, c is the speed of light, and m is the mass of the species.
      It also includes the optional parameter ``<species_name>.beta`` where beta is equal to v/c.
      The plasma will be initialized to move at bulk velocity beta*c in the
      ``<species_name>.bulk_vel_dir = (+/-) 'x', 'y', 'z'`` direction. Please leave no whitespace
      between the sign and the character on input. A direction without a sign will be treated as
      positive. The MB distribution is initialized in the drifting frame by sampling three Gaussian
      distributions in each dimension using, the Box Mueller method, and then the distribution is
      transformed to the simulation frame using the flipping method. The flipping method can be
      found in Zenitani 2015 section III. B. (Phys. Plasmas 22, 042116).

      Note that though the particles may move at relativistic speeds in the simulation frame,
      they are not relativistic in the drift frame. This is as opposed to the Maxwell Juttner
      setting, which initializes particles with relativistic momentums in their drifting frame.

    * ``maxwell_juttner``: Maxwell-Juttner distribution for high temperature plasma. This mode
      requires a dimensionless temperature parameter ``<species_name>.theta``, where theta is equal
      to kb*T/(m*c^2), where kb is the Boltzmann constant, and m is the mass of the species. It also
      includes the optional parameter ``<species_name>.beta`` where beta is equal to v/c. The plasma
      will be initialized to move at velocity beta*c in the
      ``<species_name>.bulk_vel_dir = (+/-) 'x', 'y', 'z'`` direction. Please leave no whitespace
      between the sign and the character on input. A direction without a sign will be treated as
      positive. The MJ distribution will be initialized in the moving frame using the Sobol method,
      and then the distribution will be transformed to the simulation frame using the flipping method.
      Both the Sobol and the flipping method can be found in Zenitani 2015 (Phys. Plasmas 22, 042116).

      Please take notice that particles initialized with this setting can be relativistic in two ways.
      In the simulation frame, they can drift with a relativistic speed beta. Then, in the drifting
      frame they are still moving with relativistic speeds due to high temperature. This is as opposed
      to the Maxwell Boltzmann setting, which initializes non-relativistic plasma in their relativistic
      drifting frame.

    * ``radial_expansion``: momentum depends on the radial coordinate linearly. This
      requires additional parameter ``u_over_r`` which is the slope.

    * ``parse_momentum_function``: the momentum is given by a function in the input
      file. It requires additional arguments ``<species_name>.momentum_function_ux(x,y,z)``,
      ``<species_name>.momentum_function_uy(x,y,z)`` and ``<species_name>.momentum_function_uz(x,y,z)``,
      which gives the distribution of each component of the momentum as a function of space.

* ``<species_name>.zinject_plane`` (`float`)
    Only read if  ``<species_name>`` is in ``particles.rigid_injected_species``.
    Injection plane when using the rigid injection method.
    See ``particles.rigid_injected_species`` above.

* ``<species_name>.rigid_advance`` (`bool`)
    Only read if ``<species_name>`` is in ``particles.rigid_injected_species``.

    * If ``false``, each particle is advanced with its
      own velocity ``vz`` until it reaches ``zinject_plane``.

    * If ``true``, each particle is advanced with the average speed of the species
      ``vzbar`` until it reaches ``zinject_plane``.

* ``species_name.predefined_profile_name`` (`string`)
    Only read of ``<species_name>.electrons.profile`` is `predefined`.

    * If ``parabolic_channel``, the plasma profile is a parabolic profile with
      cosine-like ramps at the beginning and the end of the profile.
      The density is given by

      .. math::

          n = n_0 n(x,y) n(z)

      with

      .. math::

          n(x,y) = 1 + 4\frac{x^2+y^2}{k_p^2 R_c^4}

      where :math:`k_p` is the plasma wavenumber associated with density :math:`n_0`.
      Here, :math:`n(z)` is a cosine-like up-ramp from :math:`0` to :math:`L_{ramp,up}`,
      constant to :math:`1` from :math:`L_{ramp,up}` to :math:`L_{ramp,up} + L_{plateau}`
      and a cosine-like down-ramp from :math:`L_{ramp,up} + L_{plateau}` to
      :math:`L_{ramp,up} + L_{plateau}+L_{ramp,down}`. All parameters are given
      in ``predefined_profile_params``.

* ``<species_name>.predefined_profile_params`` (list of `float`)
    Parameters for the predefined profiles.

    * If ``species_name.predefined_profile_name`` is ``parabolic_channel``,
      ``predefined_profile_params`` contains a space-separated list of the
      following parameters, in this order: :math:`L_{ramp,up}` :math:`L_{plateau}`
      :math:`L_{ramp,down}` :math:`R_c` :math:`n_0`

* ``<species_name>.do_backward_propagation`` (`bool`)
    Inject a backward-propagating beam to reduce the effect of charge-separation
    fields when running in the boosted frame. See examples.

* ``<species_name>.do_splitting`` (`bool`) optional (default `0`)
    Split particles of the species when crossing the boundary from a lower
    resolution domain to a higher resolution domain.

* ``<species_name>.split_type`` (`int`) optional (default `0`)
    Splitting technique. When `0`, particles are split along the simulation
    axes (4 particles in 2D, 6 particles in 3D). When `1`, particles are split
    along the diagonals (4 particles in 2D, 8 particles in 3D).

* ``<species_name>.do_not_deposit`` (`0` or `1` optional; default `0`)
    If `1` is given, both charge deposition and current deposition will
    not be done, thus that species does not contribute to the fields.

* ``<species_name>.do_not_gather`` (`0` or `1` optional; default `0`)
    If `1` is given, field gather from grids will not be done,
    thus that species will not be affected by the field on grids.

* ``<species_name>.do_not_push`` (`0` or `1` optional; default `0`)
    If `1` is given, this species will not be pushed
    by any pusher during the simulation.

* ``<species>.do_back_transformed_diagnostics`` (`0` or `1` optional, default `1`)
    Only used when ``warpx.do_back_transformed_diagnostics=1``. When running in a
    boosted frame, whether or not to plot back-transformed diagnostics for
    this species.

* ``warpx.serialize_ics`` (`0 or 1`)
    Whether or not to use OpenMP threading for particle initialization.

* ``<species>.do_field_ionization`` (`0` or `1`) optional (default `0`)
    Do field ionization for this species (using the ADK theory).

* ``<species>.physical_element`` (`string`)
    Only read if `do_field_ionization = 1`. Symbol of chemical element for
    this species. Example: for Helium, use ``physical_element = He``.
    Elements up to atomic number Z=86 (Radon) are supported, let us know if you need higher Z.

* ``<species>.ionization_product_species`` (`string`)
    Only read if `do_field_ionization = 1`. Name of species in which ionized
    electrons are stored. This species must be created as a regular species
    in the input file (in particular, it must be in `particles.species_names`).

* ``<species>.ionization_initial_level`` (`int`) optional (default `0`)
    Only read if `do_field_ionization = 1`. Initial ionization level of the
    species (must be smaller than the atomic number of chemical element given
    in `physical_element`).

* ``<species>.do_classical_radiation_reaction`` (`int`) optional (default `0`)
    Enables Radiation Reaction (or Radiation Friction) for the species. Species
    must be either electrons or positrons. Boris pusher must be used for the
    simulation

* ``<species>.do_qed`` (`int`) optional (default `0`)
    If `<species>.do_qed = 0` all the QED effects are disabled for this species.
    If `<species>.do_qed = 1` QED effects can be enabled for this species (see below).
    **This feature requires to compile with QED=TRUE**

* ``<species>.do_qed_quantum_sync`` (`int`) optional (default `0`)
    It only works if `<species>.do_qed = 1`. Enables Quantum synchrotron emission for this species.
    Quantum synchrotron lookup table should be either generated or loaded from disk to enable
    this process (see "Lookup tables for QED modules" section below).
    `<species>` must be either an electron or a positron species.
    **This feature requires to compile with QED=TRUE**

* ``<species>.do_qed_breit_wheeler`` (`int`) optional (default `0`)
    It only works if `<species>.do_qed = 1`. Enables non-linear Breit-Wheeler process for this species.
    Breit-Wheeler lookup table should be either generated or loaded from disk to enable
    this process (see "Lookup tables for QED modules" section below).
    `<species>` must be a photon species.
    **This feature requires to compile with QED=TRUE**

* ``<species>.qed_quantum_sync_phot_product_species`` (`string`)
    If an electron or a positron species has the Quantum synchrotron process, a photon product species must be specified
    (the name of an existing photon species must be provided)
    **This feature requires to compile with QED=TRUE**

* ``<species>.qed_breit_wheeler_ele_product_species`` (`string`)
    If a photon species has the Breit-Wheeler process, an electron product species must be specified
    (the name of an existing electron species must be provided)
    **This feature requires to compile with QED=TRUE**

* ``<species>.qed_breit_wheeler_pos_product_species`` (`string`)
    If a photon species has the Breit-Wheeler process, a positron product species must be specified
    (the name of an existing positron species must be provided).
    **This feature requires to compile with QED=TRUE**

* ``<species>.do_resampling`` (`0` or `1`) optional (default `0`)
    If `1` resampling is performed for this species. This means that the number of macroparticles
    will be reduced at specific timesteps while preserving the distribution function as much as
    possible (in particular the weight of the remaining particles will be increased on average).
    This can be useful in situations with continuous creation of particles (e.g. with ionization
    or with QED effects). At least one resampling trigger (see below) must be specified to actually
    perform resampling.

* ``<species>.resampling_algorithm`` (`string`) optional (default `leveling_thinning`)
    The algorithm used for resampling. Currently there is only one option, which is already set by
    default:

    * ``leveling_thinning`` This algorithm is defined in `Muraviev et al., arXiv:2006.08593 (2020) <https://arxiv.org/abs/2006.08593>`_.
      It has two parameters:

        * ``<species>.resampling_algorithm_target_ratio`` (`float`) optional (default `1.5`)
            This **roughly** corresponds to the ratio between the number of particles before and
            after resampling.

        * ``<species>.resampling_algorithm_min_ppc`` (`int`) optional (default `1`)
            Resampling is not performed in cells with a number of macroparticles strictly smaller
            than this parameter.

* ``<species>.resampling_trigger_intervals`` (`string`) optional (default `0`)
    Using the `Intervals parser`_ syntax, this string defines timesteps at which resampling is
    performed.

* ``<species>.resampling_trigger_max_avg_ppc`` (`float`) optional (default `infinity`)
    Resampling is performed everytime the number of macroparticles per cell of the species
    averaged over the whole simulation domain exceeds this parameter.

.. _running-cpp-parameters-laser:

Laser initialization
--------------------

* ``lasers.names`` (list of `string`)
    Name of each laser. This is then used in the rest of the input deck ;
    in this documentation we use `<laser_name>` as a placeholder. The parameters below
    must be provided for each laser pulse.

* ```<laser_name>`.position`` (`3 floats in 3D and 2D` ; in meters)
    The coordinates of one of the point of the antenna that will emit the laser.
    The plane of the antenna is entirely defined by ``<laser_name>.position``
    and ``<laser_name>.direction``.

    ```<laser_name>`.position`` also corresponds to the origin of the coordinates system
    for the laser tranverse profile. For instance, for a Gaussian laser profile,
    the peak of intensity will be at the position given by ``<laser_name>.position``.
    This variable can thus be used to shift the position of the laser pulse
    transversally.

    .. note::
        In 2D, ```<laser_name>`.position`` is still given by 3 numbers,
        but the second number is ignored.

    When running a **boosted-frame simulation**, provide the value of
    ``<laser_name>.position`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame. Note that,
    in this case, the laser antenna will be moving, in the boosted frame.

* ``<laser_name>.polarization`` (`3 floats in 3D and 2D`)
    The coordinates of a vector that points in the direction of polarization of
    the laser. The norm of this vector is unimportant, only its direction matters.

    .. note::
        Even in 2D, all the 3 components of this vectors are important (i.e.
        the polarization can be orthogonal to the plane of the simulation).

*  ``<laser_name>.direction`` (`3 floats in 3D`)
    The coordinates of a vector that points in the propagation direction of
    the laser. The norm of this vector is unimportant, only its direction matters.

    The plane of the antenna that will emit the laser is orthogonal to this vector.

    .. warning::

        When running **boosted-frame simulations**, ``<laser_name>.direction`` should
        be parallel to ``warpx.boost_direction``, for now.

* ``<laser_name>.e_max`` (`float` ; in V/m)
    Peak amplitude of the laser field.

    For a laser with a wavelength :math:`\lambda = 0.8\,\mu m`, the peak amplitude
    is related to :math:`a_0` by:

    .. math::

        E_{max} = a_0 \frac{2 \pi m_e c}{e\lambda} = a_0 \times (4.0 \cdot 10^{12} \;V.m^{-1})

    When running a **boosted-frame simulation**, provide the value of ``<laser_name>.e_max``
    in the laboratory frame, and use ``warpx.gamma_boost`` to automatically
    perform the conversion to the boosted frame.

* ``<laser_name>.wavelength`` (`float`; in meters)
    The wavelength of the laser in vacuum.

    When running a **boosted-frame simulation**, provide the value of
    ``<laser_name>.wavelength`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

* ``<laser_name>.profile`` (`string`)
    The spatio-temporal shape of the laser. The options that are currently
    implemented are:

    - ``"Gaussian"``: The transverse and longitudinal profiles are Gaussian.
    - ``"Harris"``: The transverse profile is Gaussian, but the longitudinal profile
      is given by the Harris function (see ``<laser_name>.profile_duration`` for more details)
    - ``"parse_field_function"``: the laser electric field is given by a function in the
      input file. It requires additional argument ``<laser_name>.field_function(X,Y,t)``, which
      is a mathematical expression , e.g.
      ``<laser_name>.field_function(X,Y,t) = "a0*X**2 * (X>0) * cos(omega0*t)"`` where
      ``a0`` and ``omega0`` are a user-defined constant, see above. The profile passed
      here is the full profile, not only the laser envelope. ``t`` is time and ``X``
      and ``Y`` are coordinates orthogonal to ``<laser_name>.direction`` (not necessarily the
      x and y coordinates of the simulation). All parameters above are required, but
      none of the parameters below are used when ``<laser_name>.parse_field_function=1``. Even
      though ``<laser_name>.wavelength`` and ``<laser_name>.e_max`` should be included in the laser
      function, they still have to be specified as they are used for numerical purposes.
    - ``"from_txye_file"``: the electric field of the laser is read from an external binary file
      whose format is explained below. It requires to provide the name of the binary file
      setting the additional parameter ``<laser_name>.txye_file_name`` (string). It accepts an
      optional parameter ``<laser_name>.time_chunk_size`` (int). This allows to read only
      time_chunk_size timesteps from the binary file. New timesteps are read as soon as they are needed.
      The default value is automatically set to the number of timesteps contained in the binary file
      (i.e. only one read is performed at the beginning of the simulation).
      The external binary file should provide E(x,y,t) on a rectangular (but non necessarily uniform)
      grid. The code performs a bi-linear (in 2D) or tri-linear (in 3D) interpolation to set the field
      values. x,y,t are meant to be in S.I. units, while the field value is meant to be multiplied by
      ``<laser_name>.e_max`` (i.e. in most cases the maximum of abs(E(x,y,t)) should be 1,
      so that the maximum field intensity can be set straightforwardly with ``<laser_name>.e_max``).
      The binary file has to respect the following format:

        * flag to indicate if the grid is uniform or not (1 byte, 0 means non-uniform, !=0 means uniform)

        * np, number of timesteps (uint32_t, must be >=2)

        * nx, number of points along x (uint32_t, must be >=2)

        * ny, number of points along y (uint32_t, must be 1 for 2D simulations and >=2 for 3D simulations)

        * timesteps (double[2] if grid is uniform, double[np] otherwise)

        * x_coords (double[2] if grid is uniform, double[nx] otherwise)

        * y_coords (double[1] if 2D, double[2] if 3D & uniform grid, double[ny] if 3D & non uniform grid)

        * field_data (double[nt * nx * ny], with nt being the slowest coordinate).

      A file at this format can be generated from Python, see an example at ``Examples/Modules/laser_injection_from_file``


*  ``<laser_name>.profile_t_peak`` (`float`; in seconds)
    The time at which the laser reaches its peak intensity, at the position
    given by ``<laser_name>.position`` (only used for the ``"gaussian"`` profile)

    When running a **boosted-frame simulation**, provide the value of
    ``<laser_name>.profile_t_peak`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

*  ``<laser_name>.profile_duration`` (`float` ; in seconds)

    The duration of the laser pulse, defined as :math:`\tau` below:

    - For the ``"gaussian"`` profile:

    .. math::

        E(\boldsymbol{x},t) \propto \exp\left( -\frac{(t-t_{peak})^2}{\tau^2} \right)

    Note that :math:`\tau` relates to the full width at half maximum (FWHM) of *intensity*, which is closer to pulse length measurements in experiments, as :math:`\tau = \mathrm{FWHM}_I / \sqrt{2\ln(2)}` :math:`\approx \mathrm{FWHM}_I / 1.174`.

    - For the ``"harris"`` profile:

    .. math::

        E(\boldsymbol{x},t) \propto \frac{1}{32}\left[10 - 15 \cos\left(\frac{2\pi t}{\tau}\right) + 6 \cos\left(\frac{4\pi t}{\tau}\right) - \cos\left(\frac{6\pi t}{\tau}\right) \right]\Theta(\tau - t)

    When running a **boosted-frame simulation**, provide the value of
    ``<laser_name>.profile_duration`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

* ``<laser_name>.profile_waist`` (`float` ; in meters)
    The waist of the transverse Gaussian laser profile, defined as :math:`w_0` :

    .. math::

        E(\boldsymbol{x},t) \propto \exp\left( -\frac{\boldsymbol{x}_\perp^2}{w_0^2} \right)

* ``<laser_name>.profile_focal_distance`` (`float`; in meters)
    The distance from ``laser_position`` to the focal plane.
    (where the distance is defined along the direction given by ``<laser_name>.direction``.)

    Use a negative number for a defocussing laser instead of a focussing laser.

    When running a **boosted-frame simulation**, provide the value of
    ``<laser_name>.profile_focal_distance`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

*  ``<laser_name>.phi0`` (`float`; in radians)
    The Carrier Envelope Phase, i.e. the phase of the laser oscillation, at the
    position where the laser enveloppe is maximum (only used for the ``"gaussian"`` profile)

* ``<laser_name>.stc_direction`` (`3 floats`) optional (default `1. 0. 0.`)
    Direction of laser spatio-temporal couplings.
    See definition in Akturk et al., Opt Express, vol 12, no 19 (2004).

* ``<laser_name>.zeta`` (`float`; in meters.seconds) optional (default `0.`)
    Spatial chirp at focus in direction ``<laser_name>.stc_direction``. See definition in
    Akturk et al., Opt Express, vol 12, no 19 (2004).

* ``<laser_name>.beta`` (`float`; in seconds) optional (default `0.`)
    Angular dispersion (or angular chirp) at focus in direction ``<laser_name>.stc_direction``.
    See definition in Akturk et al., Opt Express, vol 12, no 19 (2004).

* ``<laser_name>.phi2`` (`float`; in seconds**2) optional (default `0.`)
    Temporal chirp at focus.
    See definition in Akturk et al., Opt Express, vol 12, no 19 (2004).

* ``<laser_name>.do_continuous_injection`` (`0` or `1`) optional (default `0`).
    Whether or not to use continuous injection.
    If the antenna starts outside of the simulation domain but enters it
    at some point (due to moving window or moving antenna in the boosted
    frame), use this so that the laser antenna is injected when it reaches
    the box boundary. If running in a boosted frame, this requires the
    boost direction, moving window direction and laser propagation direction
    to be along `z`. If not running in a boosted frame, this requires the
    moving window and laser propagation directions to be the same (`x`, `y`
    or `z`)

* ``<laser_name>.min_particles_per_mode`` (`int`) optional (default `4`)
    When using the RZ version, this specifies the minimum number of particles
    per angular mode. The laser particles are loaded into radial spokes, with
    the number of spokes given by min_particles_per_mode*(warpx.n_rz_azimuthal_modes-1).

* ``warpx.num_mirrors`` (`int`) optional (default `0`)
    Users can input perfect mirror condition inside the simulation domain.
    The number of mirrors is given by ``warpx.num_mirrors``. The mirrors are
    orthogonal to the `z` direction. The following parameters are required
    when ``warpx.num_mirrors`` is >0.

* ``warpx.mirror_z`` (list of `float`) required if ``warpx.num_mirrors>0``
    ``z`` location of the front of the mirrors.

* ``warpx.mirror_z_width`` (list of `float`) required if ``warpx.num_mirrors>0``
    ``z`` width of the mirrors.

* ``warpx.mirror_z_npoints`` (list of `int`) required if ``warpx.num_mirrors>0``
    In the boosted frame, depending on `gamma_boost`, ``warpx.mirror_z_width``
    can be smaller than the cell size, so that the mirror would not work. This
    parameter is the minimum number of points for the mirror. If
    ``mirror_z_width < dz/cell_size``, the upper bound of the mirror is increased
    so that it contains at least ``mirror_z_npoints``.

* ``warpx.B_ext_grid_init_style`` (string) optional (default is "default")
    This parameter determines the type of initialization for the external
    magnetic field. The "default" style initializes the
    external magnetic field (Bx,By,Bz) to (0.0, 0.0, 0.0).
    The string can be set to "constant" if a constant magnetic field is
    required to be set at initialization. If set to "constant", then an
    additional parameter, namely, ``warpx.B_external_grid`` must be specified.
    If set to ``parse_B_ext_grid_function``, then a mathematical expression can
    be used to initialize the external magnetic field on the grid. It
    requires additional parameters in the input file, namely,
    ``warpx.Bx_external_grid_function(x,y,z)``,
    ``warpx.By_external_grid_function(x,y,z)``,
    ``warpx.Bz_external_grid_function(x,y,z)`` to initialize the external
    magnetic field for each of the three components on the grid.
    Constants required in the expression can be set using ``my_constants``.
    For example, if ``warpx.Bx_external_grid_function(x,y,z)=Bo*x + delta*(y + z)``
    then the constants `Bo` and `delta` required in the above equation
    can be set using ``my_constants.Bo=`` and ``my_constants.delta=`` in the
    input file. For a two-dimensional simulation, it is assumed that the first dimension     is `x` and the second dimension in `z`, and the value of `y` is set to zero.
    Note that the current implementation of the parser for external B-field
    does not work with RZ and the code will abort with an error message.

* ``warpx.E_ext_grid_init_style`` (string) optional (default is "default")
    This parameter determines the type of initialization for the external
    electric field. The "default" style initializes the
    external electric field (Ex,Ey,Ez) to (0.0, 0.0, 0.0).
    The string can be set to "constant" if a constant electric field is
    required to be set at initialization. If set to "constant", then an
    additional parameter, namely, ``warpx.E_external_grid`` must be specified
    in the input file.
    If set to ``parse_E_ext_grid_function``, then a mathematical expression can
    be used to initialize the external magnetic field on the grid. It
    required additional parameters in the input file, namely,
    ``warpx.Ex_external_grid_function(x,y,z)``,
    ``warpx.Ey_external_grid_function(x,y,z)``,
    ``warpx.Ez_external_grid_function(x,y,z)`` to initialize the external
    electric field for each of the three components on the grid.
    Constants required in the expression can be set using ``my_constants``.
    For example, if ``warpx.Ex_external_grid_function(x,y,z)=Eo*x + delta*(y + z)``
    then the constants `Bo` and `delta` required in the above equation
    can be set using ``my_constants.Eo=`` and ``my_constants.delta=`` in the
    input file. For a two-dimensional simulation, it is assumed that the first
    dimension is `x` and the second dimension in `z`,
    and the value of `y` is set to zero.
    Note that the current implementation of the parser for external E-field
    does not work with RZ and the code will abort with an error message.

* ``warpx.E_external_grid`` & ``warpx.B_external_grid`` (list of `3 floats`)
    required when ``warpx.E_ext_grid_init_style="constant"``
    and when ``warpx.B_ext_grid_init_style="constant"``, respectively.
    External uniform and constant electrostatic and magnetostatic field added
    to the grid at initialization. Use with caution as these fields are used for
    the field solver. In particular, do not use any other boundary condition
    than periodic.

*  ``particles.B_ext_particle_init_style`` (string) optional (default is "default")
     This parameter determines the type of initialization for the external
     magnetic field that is applied directly to the particles at every timestep.
     The "default" style sets the external B-field (Bx,By,Bz) to (0.0,0.0,0.0).
     The string can be set to "constant" if a constant external B-field is applied
     every timestep. If this parameter is set to "constant", then an additional
     parameter, namely, ``particles.B_external_particle`` must be specified in
     the input file.
     To parse a mathematical function for the external B-field, use the option
     ``parse_B_ext_particle_function``. This option requires additional parameters
     in the input file, namely,
     ``particles.Bx_external_particle_function(x,y,z,t)``,
     ``particles.By_external_particle_function(x,y,z,t)``,
     ``particles.Bz_external_particle_function(x,y,z,t)`` to apply the external B-field
     on the particles. Constants required in the mathematical expression can be set
     using ``my_constants``. For a two-dimensional simulation, it is assumed that
     the first and second dimensions are `x` and `z`, respectively, and the
     value of the `By` component is set to zero.
     Note that the current implementation of the parser for B-field on particles
     is applied in cartesian co-ordinates as a function of (x,y,z) even for RZ.

*    ``particles.E_ext_particle_init_style`` (string) optional (default is "default")
     This parameter determines the type of initialization for the external
     electric field that is applied directly to the particles at every timestep.
     The "default" style set the external E-field (Ex,Ey,Ez) to (0.0,0.0,0.0).
     The string can be set to "constant" if a constant external E-field is to be
     used in the simulation at every timestep. If this parameter is set to "constant",
     then an additional parameter, namely, ``particles.E_external_particle`` must be
     specified in the input file.
     To parse a mathematical function for the external E-field, use the option
     ``parse_E_ext_particle_function``. This option requires additional
     parameters in the input file, namely,
     ``particles.Ex_external_particle_function(x,y,z,t)``,
     ``particles.Ey_external_particle_function(x,y,z,t)``,
     ``particles.Ez_external_particle_function(x,y,z,t)`` to apply the external E-field
     on the particles. Constants required in the mathematical expression can be set
     using ``my_constants``. For a two-dimensional simulation, similar to the B-field,
     it is assumed that the first and second dimensions are `x` and `z`, respectively,
     and the value of the `Ey` component is set to zero.
     The current implementation of the parser for B-field on particles
     is applied in cartesian co-ordinates as a function of (x,y,z) even for RZ.

* ``particles.E_external_particle`` & ``particles.B_external_particle`` (list of `float`) optional (default `0. 0. 0.`)
    Two separate parameters which add an externally applied uniform E-field or
    B-field to each particle which is then added to the field values gathered
    from the grid in the PIC cycle.

.. _running-cpp-parameters-collision:

Collision initialization
------------------------

WarpX provides a relativistic elastic Monte Carlo binary collision model,
following the algorithm given by `Perez et al. (Phys. Plasmas 19, 083104, 2012) <https://doi.org/10.1063/1.4742167>`_.

* ``collisions.collision_names`` (`strings`, separated by spaces)
    The name of each collision type.
    This is then used in the rest of the input deck;
    in this documentation we use ``<collision_name>`` as a placeholder.

* ``<collision_name>.species`` (`strings`, two species names separated by spaces)
    The names of two species, between which the collision will be considered.
    The number of provided ``<collision_name>.species`` should match
    the number of collision names, i.e. ``collisions.collision_names``.

* ``<collision_name>.CoulombLog`` (`float`) optional
    A provided fixed Coulomb logarithm of the collision type
    ``<collision_name>``.
    For example, a typical Coulomb logarithm has a form of
    :math:`\ln(\lambda_D/R)`,
    where :math:`\lambda_D` is the Debye length,
    :math:`R\approx1.4A^{1/3}` is the effective Coulombic radius of the nucleus,
    :math:`A` is the mass number.
    If this is not provided, or if a non-positive value is provided,
    a Coulomb logarithm will be computed automatically according to the algorithm.
    a Coulomb logarithm will be computed automatically according to the algorithm in
    `Perez et al. (Phys. Plasmas 19, 083104, 2012) <https://doi.org/10.1063/1.4742167>`_.

* ``<collision_name>.ndt`` (`int`) optional
    Execute collision every # time steps.
    The default value is 1.

.. _running-cpp-parameters-numerics:

Numerics and algorithms
-----------------------

* ``warpx.cfl`` (`float`)
    The ratio between the actual timestep that is used in the simulation
    and the Courant-Friedrichs-Lewy (CFL) limit. (e.g. for `warpx.cfl=1`,
    the timestep will be exactly equal to the CFL limit.)

* ``warpx.use_filter`` (`0 or 1`)
    Whether to smooth the charge and currents on the mesh, after depositing
    them from the macroparticles. This uses a bilinear filter
    (see the sub-section **Filtering** in :doc:`../theory/theory`).
    When using the RZ spectral solver, the filtering is done in k-space.

* ``warpx.filter_npass_each_dir`` (`3 int`) optional (default `1 1 1`)
    Number of passes along each direction for the bilinear filter.
    In 2D simulations, only the first two values are read.

* ``warpx.use_filter_compensation`` (`0` or `1`; default: `0`)
    Whether to add compensation when applying filtering.
    This is only supported with the RZ spectral solver.

* ``warpx.use_damp_fields_in_z_guard`` (`0` or `1`)
    When using the RZ spectrol solver, specifies whether to apply a
    damping factor to the E and B fields in the guard cells
    along z that extend beyond the edge of the domain.
    When the boundary conditions along z are not periodic, this defaults to
    true, otherwise false. The damping profile is
    a sine squared and is applied to the fields on the outer half of the guards.
    This damping is useful for damping high frequency numerical artifacts that
    occur when there is parallel decomposition along z with non-periodic boundary
    conditions.

* ``algo.current_deposition`` (`string`, optional)
    This parameter selects the algorithm for the deposition of the current density.
    Available options are: ``direct``, ``esirkepov``, and ``vay``. The default choice
    is ``esirkepov`` for FDTD maxwell solvers and ``direct`` for standard or
    Galilean PSATD solver (that is, with ``algo.maxwell_solver = psatd``).

    1. ``direct``

       The current density is deposited as described in the section :ref:`current_deposition`.
       This deposition scheme does not conserve charge.

    2. ``esirkepov``

       The current density is deposited as described in
       `(Esirkepov, CPC, 2001) <https://www.sciencedirect.com/science/article/pii/S0010465500002289>`_.
       This deposition scheme guarantees charge conservation for shape factors of arbitrary order.

    3. ``vay``

       The current density is deposited as described in `(Vay et al, 2013) <https://doi.org/10.1016/j.jcp.2013.03.010>`_ (see section :ref:`current_deposition` for more details).
       This option guarantees charge conservation only when used in combination
       with ``psatd.periodic_single_box_fft=1``, that is, only for periodic single-box
       simulations with global FFTs without guard cells. The implementation for domain
       decomposition with local FFTs over guard cells is planned but not yet completed.

* ``algo.charge_deposition`` (`string`, optional)
    The algorithm for the charge density deposition. Available options are:

     - ``standard``: standard charge deposition algorithm, described in
       the section :doc:`../theory/picsar_theory`.

* ``algo.field_gathering`` (`string`, optional)
    The algorithm for field gathering. Available options are:

     - ``energy-conserving``: gathers directly from the grid points (either staggered
       or nodal gridpoints depending on ``warpx.do_nodal``).
     - ``momentum-conserving``: first average the fields from the grid points to
       the nodes, and then gather from the nodes.

     If ``algo.field_gathering`` is not specified, the default is ``energy-conserving``.
     If ``warpx.do_nodal`` is ``true``, then ``energy-conserving`` and ``momentum-conserving``
     are equivalent.


* ``algo.particle_pusher`` (`string`, optional)
    The algorithm for the particle pusher. Available options are:

     - ``boris``: Boris pusher.
     - ``vay``: Vay pusher (see `Vay, Phys. Plasmas (2008) <https://aip.scitation.org/doi/10.1063/1.2837054>`__)
     - ``higuera``: Higuera-Cary pusher (see `Higuera and Cary, Phys. Plasmas (2017) <https://aip.scitation.org/doi/10.1063/1.4979989>`__)

     If ``algo.particle_pusher`` is not specified, ``boris`` is the default.

* ``algo.maxwell_solver`` (`string`, optional)
    The algorithm for the Maxwell field solver.
    Available options are:

     - ``yee``: Yee FDTD solver.
     - ``ckc``: (not available in ``RZ`` geometry) Cole-Karkkainen solver with Cowan
       coefficients (see `Cowan, PRSTAB 16 (2013) <https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.16.041303>`__)
     - ``psatd``: Pseudo-spectral solver (see :ref:`theory <theory-pic-mwsolve-psatd>`)

     If ``algo.maxwell_solver`` is not specified, ``yee`` is the default.

* ``algo.em_solver_medium`` (`string`, optional)
    The medium for evaluating the Maxwell solver. Available options are :

    - ``vacuum``: vacuum properties are used in the Maxwell solver.
    - ``macroscopic``: macroscopic Maxwell equation is evaluated. If this option is selected, then the corresponding properties of the medium must be provided using ``macroscopic.sigma``, ``macroscopic.epsilon``, and ``macroscopic.mu`` for each case where the initialization style is ``constant``.  Otherwise if the initialization style uses the parser, ``macroscopic.sigma_function(x,y,z)``, ``macroscopic.epsilon_function(x,y,z)`` and/or ``macroscopic.mu_function(x,y,z)`` must be provided using the parser initialization style for spatially varying macroscopic properties.

    If ``algo.em_solver_medium`` is not specified, ``vacuum`` is the default.

* ``algo.macroscopic_sigma_method`` (`string`, optional)
    The algorithm for updating electric field when ``algo.em_solver_medium`` is macroscopic. Available options are:

    - ``backwardeuler`` is a fully-implicit, first-order in time scheme for E-update (default).
    - ``laxwendroff`` is the semi-implicit, second order in time scheme for E-update.
    Comparing the two methods, Lax-Wendroff is more prone to developing oscillations and requires a smaller timestep for stability. On the other hand, Backward Euler is more robust but it is first-order accurate in time compared to the second-order Lax-Wendroff method.

* ``macroscopic.sigma_function(x,y,z)``, ``macroscopic.epsilon_function(x,y,z)``, ``macroscopic.mu_function(x,y,z)`` (`string`)
     To initialize spatially varying conducitivy, permittivity, and permeability, respectively,
     using a mathematical function in the input. Constants required in the
     mathematical expression can be set using ``my_constants``. These parameters are parsed
     if ``algo.em_solver_medium=macroscopic``.

* ``macroscopic.sigma``, ``macroscopic.epsilon``, ``macroscopic.mu`` (`double`)
    To initialize a constant conductivity, permittivity, and permeability of the
    computational medium, respectively. The default values are the corresponding values
    in vacuum.

* ``interpolation.nox``, ``interpolation.noy``, ``interpolation.noz`` (`1`, `2`, or `3` ; default: 1)
    The order of the shape factors for the macroparticles, for the 3 dimensions of space.
    Lower-order shape factors result in faster simulations, but more noisy results,

    Note that in the current implementation in WarpX these 3 numbers must be equal.

* ``interpolation.galerkin_scheme`` (`0` or `1`)
    Whether to use a Galerkin scheme when gathering fields to particles.
    When set to `1`, the interpolation orders used for field-gathering are reduced for certain field components along certain directions.
    For example, `E_z` is gathered using ``interpolation.nox``, ``interpolation.noy``, and ``interpolation.noz - 1``.
    See equations 21-23 of (`Godfrey and Vay, 2013 <https://doi.org/10.1016/j.jcp.2013.04.006>`_) and associated references for details.
    Defaults to `1` unless ``warpx.do_nodal = 1`` and/or ``algo.field_gathering = momentum-conserving``.

* ``interpolation.field_gathering_nox``, ``interpolation.field_gathering_noy``, ``interpolation.field_gathering_noz`` (default: ``2`` in all directions)
    The order of the interpolation used with staggered grids (``warpx.do_nodal = 0``) and momentum-conserving field gathering (``algo.field_gathering = momentum-conserving``) to interpolate the electric and magnetic fields from the cell centers to the cell nodes, before gathering the fields from the cell nodes to the particle positions. High-order interpolation (order 8 in each direction, at least) is necessary to ensure stability in typical LWFA boosted-frame simulations using the Galilean PSATD or comoving PSATD schemes. This arbitrary-order interpolation is used only when the PSATD solver is used for Maxwell's equations. With the FDTD solver, basic linear interpolation is used instead.

* ``warpx.do_dive_cleaning`` (`0` or `1` ; default: 0)
    Whether to use modified Maxwell equations that progressively eliminate
    the error in :math:`div(E)-\rho`. This can be useful when using a current
    deposition algorithm which is not strictly charge-conserving, or when
    using mesh refinement. These modified Maxwell equation will cause the error
    to propagate (at the speed of light) to the boundaries of the simulation
    domain, where it can be absorbed.

* ``warpx.do_nodal`` (`0` or `1` ; default: 0)
    Whether to use a nodal grid (i.e. all fields are defined at the
    same points in space) or a staggered grid (i.e. Yee grid ; different
    fields are defined at different points in space)

* ``warpx.do_subcycling`` (`0` or `1`; default: 0)
    Whether or not to use sub-cycling. Different refinement levels have a
    different cell size, which results in different Courant–Friedrichs–Lewy
    (CFL) limits for the time step. By default, when using mesh refinement,
    the same time step is used for all levels. This time step is
    taken as the CFL limit of the finest level. Hence, for coarser
    levels, the timestep is only a fraction of the CFL limit for this
    level, which may lead to numerical artifacts. With sub-cycling, each level
    evolves with its own time step, set to its own CFL limit. In practice, it
    means that when level 0 performs one iteration, level 1 performs two
    iterations. Currently, this option is only supported when
    ``amr.max_level = 1``. More information can be found at
    https://ieeexplore.ieee.org/document/8659392.

* ``psatd.nox``, ``psatd.noy``, ``pstad.noz`` (`integer`) optional (default `16` for all)
    The order of accuracy of the spatial derivatives, when using the code compiled with a PSATD solver.
    If ``psatd.periodic_single_box_fft`` is used, these can be set to ``inf`` for infinite-order PSATD.

* ``psatd.nx_guard`, ``psatd.ny_guard``, ``psatd.nz_guard`` (`integer`) optional
    The number of guard cells to use with PSATD solver.
    If not set by users, these values are calculated automatically and determined *empirically* and
    would be equal the order of the solver for nodal grid, and half the order of the solver for staggered.

* ``psatd.periodic_single_box_fft`` (`0` or `1`; default: 0)
    If true, this will *not* incorporate the guard cells into the box over which FFTs are performed.
    This is only valid when WarpX is run with periodic boundaries and a single box.
    In this case, using `psatd.periodic_single_box_fft` is equivalent to using a global FFT over the whole domain.
    Therefore, all the approximations that are usually made when using local FFTs with guard cells
    (for problems with multiple boxes) become exact in the case of the periodic, single-box FFT without guard cells.

* ``psatd.fftw_plan_measure`` (`0` or `1`)
    Defines whether the parameters of FFTW plans will be initialized by
    measuring and optimizing performance (``FFTW_MEASURE`` mode; activated by default here).
    If ``psatd.fftw_plan_measure`` is set to ``0``, then the best parameters of FFTW
    plans will simply be estimated (``FFTW_ESTIMATE`` mode).
    See `this section of the FFTW documentation <http://www.fftw.org/fftw3_doc/Planner-Flags.html>`__
    for more information.

* ``psatd.current_correction`` (`0` or `1`; default: `0`)
    If true, a current correction scheme in Fourier space is applied in order to guarantee charge conservation.

    If ``psatd.v_galilean`` is zero, the spectral solver used is the standard PSATD scheme described in (`Vay et al, JCP 243, 2013 <https://doi.org/10.1016/j.jcp.2013.03.010>`_) and the current correction reads

    .. math::
       \widehat{\boldsymbol{J}}^{\,n+1/2}_{\mathrm{correct}} = \widehat{\boldsymbol{J}}^{\,n+1/2}
       - \bigg(\boldsymbol{k}\cdot\widehat{\boldsymbol{J}}^{\,n+1/2}
       - i \frac{\widehat{\rho}^{n+1} - \widehat{\rho}^{n}}{\Delta{t}}\bigg) \frac{\boldsymbol{k}}{k^2}

    If ``psatd.v_galilean`` is non-zero, the spectral solver used is the Galilean PSATD scheme described in (`Lehe et al, PRE 94, 2016 <https://doi.org/10.1103/PhysRevE.94.053305>`_) and the current correction reads

    .. math::
       \widehat{\boldsymbol{J}}^{\,n+1/2}_{\mathrm{correct}} = \widehat{\boldsymbol{J}}^{\,n+1/2}
       - \bigg(\boldsymbol{k}\cdot\widehat{\boldsymbol{J}}^{\,n+1/2} - (\boldsymbol{k}\cdot\boldsymbol{v}_G)
       \,\frac{\widehat\rho^{n+1} - \widehat\rho^{n}\theta^2}{1 - \theta^2}\bigg) \frac{\boldsymbol{k}}{k^2}

    where :math:`\theta=\exp(i\,\boldsymbol{k}\cdot\boldsymbol{v}_G\,\Delta{t}/2)`.

    This option is currently implemented only for the standard PSATD and Galilean PSATD schemes, while it is not yet available for the averaged Galilean PSATD scheme (activated by the input parameter ``psatd.do_time_averaging``).

    This option guarantees charge conservation only when used in combination with ``psatd.periodic_single_box_fft=1``, namely for periodic single-box simulations with global FFTs without guard cells.
    The implementation for domain decomposition with local FFTs over guard cells is planned but not yet completed.

* ``psatd.update_with_rho`` (`0` or `1`)
    If true, the update equation for the electric field is expressed in terms of both the current density and the charge density, namely :math:`\widehat{\boldsymbol{J}}^{\,n+1/2}`, :math:`\widehat\rho^{n}`, and :math:`\widehat\rho^{n+1}`.
    If false, instead, the update equation for the electric field is expressed in terms of the current density :math:`\widehat{\boldsymbol{J}}^{\,n+1/2}` only.
    If charge is expected to be conserved (by setting, for example, ``psatd.current_correction=1``), then the two formulations are expected to be equivalent.

    This option is currently implemented only for the standard PSATD and Galilean PSATD schemes, while it is not yet available for the averaged Galilean PSATD scheme (activated by the input parameter ``psatd.do_time_averaging``).

    If ``psatd.v_galilean`` is zero, the spectral solver used is the standard PSATD scheme described in (`Vay et al, JCP 243, 2013 <https://doi.org/10.1016/j.jcp.2013.03.010>`_):

    1. if ``psatd.update_with_rho=0``, the update equation for the electric field reads

    .. math::
       \begin{split}
       \widehat{\boldsymbol{E}}^{\,n+1}= & \:
       C \widehat{\boldsymbol{E}}^{\,n} + i \, \frac{S c}{k} \boldsymbol{k}\times\widehat{\boldsymbol{B}}^{\,n}
       - \frac{S}{\epsilon_0 c \, k} \widehat{\boldsymbol{J}}^{\,n+1/2} \\[0.2cm]
       & +\frac{1-C}{k^2} (\boldsymbol{k}\cdot\widehat{\boldsymbol{E}}^{\,n}) \boldsymbol{k}
       + \frac{1}{\epsilon_0 k^2} \left(\frac{S}{c \, k}-\Delta{t}\right)
       (\boldsymbol{k}\cdot\widehat{\boldsymbol{J}}^{\,n+1/2}) \boldsymbol{k}
       \end{split}

    2. if ``psatd.update_with_rho=1``, the update equation for the electric field reads

    .. math::
       \begin{split}
       \widehat{\boldsymbol{E}}^{\,n+1}= & \:
       C\widehat{\boldsymbol{E}}^{\,n} + i \, \frac{S c}{k} \boldsymbol{k}\times\widehat{\boldsymbol{B}}^{\,n}
       - \frac{S}{\epsilon_0 c \, k} \widehat{\boldsymbol{J}}^{\,n+1/2} \\[0.2cm]
       & + \frac{i}{\epsilon_0 k^2} \left(C-\frac{S}{c\,k}\frac{1}{\Delta{t}}\right)
       \widehat{\rho}^{n} \boldsymbol{k} - \frac{i}{\epsilon_0 k^2} \left(1-\frac{S}{c \, k}
       \frac{1}{\Delta{t}}\right)\widehat{\rho}^{n+1} \boldsymbol{k}
       \end{split}

    The coefficients :math:`C` and :math:`S` are defined in (`Vay et al, JCP 243, 2013 <https://doi.org/10.1016/j.jcp.2013.03.010>`_).

    If ``psatd.v_galilean`` is non-zero, the spectral solver used is the Galilean PSATD scheme described in (`Lehe et al, PRE 94, 2016 <https://doi.org/10.1103/PhysRevE.94.053305>`_):

    1. if ``psatd.update_with_rho=0``, the update equation for the electric field reads

    .. math::
       \begin{split}
       \widehat{\boldsymbol{E}}^{\,n+1} = & \:
       \theta^{2} C \widehat{\boldsymbol{E}}^{\,n} + i \, \theta^{2} \frac{S c}{k}
       \boldsymbol{k}\times\widehat{\boldsymbol{B}}^{\,n}
       + \frac{i \, \nu \, \theta \, \chi_1 - \theta^{2} S}{\epsilon_0 c \, k}
       \widehat{\boldsymbol{J}}^{\,n+1/2} \\[0.2cm]
       & + \theta^{2} \frac{\chi_2-\chi_3}{k^{2}}
       (\boldsymbol{k}\cdot\widehat{\boldsymbol{E}}^{\,n}) \boldsymbol{k}
       + i \, \frac{\chi_2\left(\theta^{2}-1\right)}{\epsilon_0 c \, k^{3} \nu}
       (\boldsymbol{k}\cdot\widehat{\boldsymbol{J}}^{\,n+1/2}) \boldsymbol{k}
       \end{split}

    2. if ``psatd.update_with_rho=1``, the update equation for the electric field reads

    .. math::
       \begin{split}
       \widehat{\boldsymbol{E}}^{\,n+1} = & \:
       \theta^{2} C \widehat{\boldsymbol{E}}^{\,n} + i \, \theta^{2} \frac{S c}{k}
       \boldsymbol{k}\times\widehat{\boldsymbol{B}}^{\,n}
       + \frac{i \, \nu \, \theta \, \chi_1 - \theta^{2} S}{\epsilon_0 c \, k}
       \widehat{\boldsymbol{J}}^{\,n+1/2} \\[0.2cm]
       & + i \, \frac{\theta^{2} \chi_3}{\epsilon_0 k^{2}} \widehat{\rho}^{\,n} \boldsymbol{k}
       - i \, \frac{\chi_2}{\epsilon_0 k^{2}} \widehat{\rho}^{\,n+1} \boldsymbol{k}
       \end{split}

    The coefficients :math:`C`, :math:`S`, :math:`\theta`, :math:`\nu`, :math:`\chi_1`, :math:`\chi_2`, and :math:`\chi_3` are defined in (`Lehe et al, PRE 94, 2016 <https://doi.org/10.1103/PhysRevE.94.053305>`_).

    The default value for ``psatd.update_with_rho`` is ``1`` if ``psatd.v_galilean`` is non-zero and ``0`` otherwise.

    Note that the update with and without rho is also supported in RZ geometry.

* ``pstad.v_galilean`` (`3 floats`, in units of the speed of light; default `0. 0. 0.`)
    Defines the galilean velocity.
    Non-zero `v_galilean` activates Galilean algorithm, which suppresses the Numerical Cherenkov instability
    in boosted-frame simulation. This requires the code to be compiled with `USE_PSATD=TRUE`.
    (see the sub-section Numerical Stability and alternate formulation
    in a Galilean frame in :doc:`../theory/boosted-frame`).
    It also requires the use of the `direct` current deposition option
    `algo.current_deposition = direct` (does not work with Esirkepov algorithm).

* ``psatd.v_comoving`` (3 floating-point values, in units of the speed of light; default ``0. 0. 0.``)
    Defines the comoving velocity in the comoving PSATD scheme.
    A non-zero comoving velocity selects the comoving PSATD algorithm, which suppresses the numerical Cherenkov instability (NCI) in boosted-frame simulations, under certain assumptions. This option requires that WarpX is compiled with ``USE_PSATD = TRUE``. It also requires the use of direct current deposition (``algo.current_deposition = direct``) and has not been neither implemented nor tested with other current deposition schemes.

* ``psatd.do_time_averaging`` (`0` or `1`; default: 0)
    Whether to use an averaged Galilean PSATD algorithm or standard Galilean PSATD.

* ``warpx.override_sync_int`` (`string`) optional (default `1`)
    Using the `Intervals parser`_ syntax, this string defines the timesteps at which
    synchronization of sources (`rho` and `J`) on grid nodes at box boundaries is performed.
    Since the grid nodes at the interface between two neighbor boxes are duplicated in both
    boxes, an instability can occur if they have too different values.
    This option makes sure that they are synchronized periodically.

* ``warpx.use_hybrid_QED`` ('bool'; default: 0)
    Will use the Hybird QED Maxwell solver when pushing fields: a QED correction is added to the
    field solver to solve non-linear Maxwell's equations, according to [Quantum Electrodynamics
    vacuum polarization solver, P. Carneiro et al., `ArXiv 2016 <https://arxiv.org/abs/1607.04224>`__].
    Note that this option can only be used with the PSATD build. Furthermore,
    warpx.do_nodal must be set to `1` which is not its default value.

 * ``warpx.quantum_xi`` ('float'; default: 1.3050122.e-52)
     Overwrites the actual quantum parameter used in Maxwell's QED equations. Assigning a
     value here will make the simulation unphysical, but will allow QED effects to become more apparent.
     Note that this option will only have an effect if the ``warpx.use_Hybrid_QED`` flag is also triggered.

 * ``warpx.do_device_synchronize_before_profile`` (`bool`) optional (default `1`)
    When running in an accelerated platform, whether to call a deviceSynchronize around profiling regions.
    This allows the profiler to give meaningful timers, but (hardly) slows down the simulation.

 * ``warpx.sort_int`` (`string`) optional (defaults: ``-1`` on CPU; ``4`` on GPU)
     Using the `Intervals parser`_ syntax, this string defines the timesteps at which particles are
     sorted by bin.
     If ``<=0``, do not sort particles.
     It is turned on on GPUs for performance reasons (to improve memory locality).

 * ``warpx.sort_bin_size`` (list of `int`) optional (default ``4 4 4``)
     If ``sort_int`` is activated particles are sorted in bins of ``sort_bin_size`` cells.
     In 2D, only the first two elements are read.

.. _running-cpp-parameters-boundary:

Boundary conditions
-------------------

* ``warpx.do_pml`` (`0` or `1`; default: 1)
    Whether to add Perfectly Matched Layers (PML) around the simulation box,
    and around the refinement patches. See the section :doc:`../../theory/PML`
    for more details.

* ``warpx.pml_ncell`` (`int`; default: 10)
    The depth of the PML, in number of cells.

* ``warpx.pml_delta`` (`int`; default: 10)
    The characteristic depth, in number of cells, over which
    the absorption coefficients of the PML increases.

* ``warpx.do_pml_in_domain`` (`int`; default: 0)
    Whether to create the PML inside the simulation area or outside. If inside,
    it allows the user to propagate particles in PML and to use extended PML

* ``warpx.do_pml_has_particles`` (`int`; default: 0)
    Whether to propagate particles in PML or not. Can only be done if PML are in simulation domain,
    i.e. if `warpx.do_pml_in_domain = 1`.

* ``warpx.do_pml_j_damping`` (`int`; default: 0)
    Whether to damp current in PML. Can only be used if particles are propagated in PML,
    i.e. if `warpx.do_pml_has_particles = 1`.

* ``warpx.do_pml_Lo`` (`2 ints in 2D`, `3 ints in 3D`; default: `1 1 1`)
    The directions along which one wants a pml boundary condition for lower boundaries on mother grid.

* ``warpx.do_pml_Hi`` (`2 floats in 2D`, `3 floats in 3D`; default: `1 1 1`)
    The directions along which one wants a pml boundary condition for upper boundaries on mother grid.

.. _running-cpp-parameters-diagnostics:

Diagnostics and output
----------------------

In-situ visualization
^^^^^^^^^^^^^^^^^^^^^

WarpX has three types of diagnostics:
``FullDiagnostics`` consist in dumps of fields and particles at given iterations,
``BackTransformedDiagnostics`` are used when running a simulation in a boosted frame, to reconstruct output data to the lab frame, and
``ReducedDiags`` allow the user to compute some reduced quantity (particle temperature, max of a field) and write a small amount of data to text files.
Similar to what is done for physical species, WarpX has a class Diagnostics that allows users to initialize different diagnostics, each of them with different fields, resolution and period.
This currently applies to standard diagnostics, but should be extended to back-transformed diagnostics and reduced diagnostics (and others) in a near future.

Full Diagnostics
^^^^^^^^^^^^^^^^

``FullDiagnostics`` consist in dumps of fields and particles at given iterations.
Similar to what is done for physical species, WarpX has a class Diagnostics that allows users to initialize different diagnostics, each of them with different fields, resolution and period.
The user specifies the number of diagnostics and the name of each of them, and then specifies options for each of them separately.
Note that some parameter (those that do not start with a ``<diag_name>.`` prefix) apply to all diagnostics.
This should be changed in the future.
In-situ capabilities can be used by turning on Sensei or Ascent (provided they are installed) through the output format, see below.

* ``diagnostics.enable`` (`0` or `1`, optional, default `1`)
    Whether to enable or disable diagnostics. This flag overwrites all other diagnostics input parameters.

* ``diagnostics.diags_names`` (list of `string` optional, default `empty`)
    Name of each diagnostics.
    example: ``diagnostics.diags_names = diag1 my_second_diag``.

* ``<diag_name>.period`` (`string` optional, default `0`)
    Using the `Intervals parser`_ syntax, this string defines the timesteps at which data is dumped.
    Use a negative number or 0 to disable data dumping.
    This is ``0`` (disabled) by default.
    example: ``diag1.period = 10,20:25:1``.

* ``<diag_name>.diag_type`` (`string`)
    Type of diagnostics. So far, only ``Full`` is supported.
    example: ``diag1.diag_type = Full``.

* ``<diag_name>.format`` (`string` optional, default ``plotfile``)
    Flush format. Possible values are:

    * ``plotfile`` for native AMReX format.

    * ``checkpoint`` for a checkpoint file, only works with ``<diag_name>.diag_type = Full``.

    * ``openpmd`` for OpenPMD format `openPMD <https://www.openPMD.org>`_.
      Requires to build WarpX with ``USE_OPENPMD=TRUE`` (see :ref:`instructions <building-openpmd>`).

    * ``ascent`` for in-situ visualization using Ascent.

    * ``sensei`` for in-situ visualization using Sensei.

    example: ``diag1.format = openpmd``.

* ``<diag_name>.sensei_config`` (`string`)
  Only read if ``<diag_name>.format = sensei``.
  Points to the SENSEI XML file which selects and configures the desired back end.

* ``<diag_name>.sensei_pin_mesh`` (`integer`; 0 by default)
  Only read if ``<diag_name>.format = sensei``.
  When 1 lower left corner of the mesh is pinned to 0.,0.,0.

* ``<diag_name>.openpmd_backend`` (``bp``, ``h5`` or ``json``) optional, only used if ``<diag_name>.format = openpmd``
    `I/O backend <https://openpmd-api.readthedocs.io/en/latest/backends/overview.html>`_ for `openPMD <https://www.openPMD.org>`_ data dumps.
    ``bp`` is the `ADIOS I/O library <https://csmd.ornl.gov/adios>`_, ``h5`` is the `HDF5 format <https://www.hdfgroup.org/solutions/hdf5/>`_, and ``json`` is a `simple text format <https://en.wikipedia.org/wiki/JSON>`_.
    ``json`` only works with serial/single-rank jobs.
    When WarpX is compiled with openPMD support, the first available backend in the order given above is taken.

* ``<diag_name>.openpmd_tspf`` (`bool`, optional, default ``true``) only read if ``<diag_name>.format = openpmd``.
    Whether to write one file per timestep.

* ``<diag_name>.fields_to_plot`` (list of `strings`, optional)
    Fields written to output.
    Possible values: ``Ex`` ``Ey`` ``Ez`` ``Bx`` ``By`` ``Bz`` ``jx`` ``jy`` ``jz`` ``part_per_cell`` ``rho`` ``phi`` ``F`` ``part_per_grid`` ``divE`` ``divB`` and ``rho_<species_name>``, where ``<species_name>`` must match the name of one of the available particle species. Note that ``phi`` will only be written out when do_electrostatic==labframe.
    Default is ``<diag_name>.fields_to_plot = Ex Ey Ez Bx By Bz jx jy jz``.
    Note that the fields are averaged on the cell centers before they are written to file.

* ``<diag_name>.plot_raw_fields`` (`0` or `1`) optional (default `0`)
    By default, the fields written in the plot files are averaged on the cell centers.
    When ```warpx.plot_raw_fields`` is `1`, then the raw (i.e. unaveraged)
    fields are also saved in the output files.
    Only works with ``<diag_name>.format = plotfile``.
    See `this section <https://yt-project.org/doc/examining/loading_data.html#viewing-raw-fields-in-warpx>`_
    in the yt documentation for more details on how to view raw fields.

* ``<diag_name>.plot_raw_fields_guards`` (`0` or `1`) optional (default `0`)
    Only used when ``warpx.plot_raw_fields`` is ``1``.
    Whether to include the guard cells in the output of the raw fields.
    Only works with ``<diag_name>.format = plotfile``.

* ``<diag_name>.plot_finepatch`` (`0` or `1`) optional (default `0`)
    Only used when mesh refinement is activated and ``warpx.plot_raw_fields`` is ``1``.
    Whether to output the data of the fine patch, in the plot files.
    Only works with ``<diag_name>.format = plotfile``.

* ``<diag_name>.plot_crsepatch`` (`0` or `1`) optional (default `0`)
    Only used when mesh refinement is activated and ``warpx.plot_raw_fields`` is ``1``.
    Whether to output the data of the coarse patch, in the plot files.
    Only works with ``<diag_name>.format = plotfile``.

* ``<diag_name>.coarsening_ratio`` (list of `int`) optional (default `1 1 1`)
    Reduce size of the field output by this ratio in each dimension.
    (This is done by averaging the field over 1 or 2 points along each direction, depending on the staggering).
    If ``blocking_factor`` and ``max_grid_size`` are used for the domain decomposition, as detailed in
    the :ref:`parallelization <parallelization_warpx>` section, ``coarsening_ratio`` should be an integer
    divisor of ``blocking_factor``. If ``warpx.numprocs`` is used instead, the total number of cells in a given
    dimension must be a multiple of the ``coarsening_ratio`` multiplied by ``numprocs`` in that dimension.

* ``<diag_name>.file_prefix`` (`string`) optional (default `diags/plotfiles/plt`)
    Root for output file names. Supports sub-directories.

* ``<diag_name>.diag_lo`` (list `float`, 1 per dimension) optional (default `-infinity -infinity -infinity`)
    Lower corner of the output fields (if smaller than ``warpx.dom_lo``, then set to ``warpx.dom_lo``). Currently, when the ``diag_lo`` is different from ``warpx.dom_lo``, particle output is disabled.

* ``<diag_name>.diag_hi`` (list `float`, 1 per dimension) optional (default `+infinity +infinity +infinity`)
    Higher corner of the output fields (if larger than ``warpx.dom_hi``, then set to ``warpx.dom_hi``). Currently, when the ``diag_hi`` is different from ``warpx.dom_hi``, particle output i
s disabled.

* ``<diag_name>.write_species`` (`0` or `1`) optional (default `1`)
   Whether to write species output or not. For checkpoint format, always set this parameter to 1.

* ``<diag_name>.species`` (list of `string`, default all physical species in the simulation)
    Which species dumped in this diagnostics.

* ``<diag_name>.<species_name>.variables`` (list of `strings` separated by spaces, optional)
    List of particle quantities to write to output.
    Choices are ``w`` for the particle weight and ``ux`` ``uy`` ``uz`` for the particle momenta.
    By default, all particle quantities are written.
    If ``<diag_name>.<species_name>.variables = none``, no particle data are written, except for particle positions, which are always included.

* ``<diag_name>.<species_name>.random_fraction`` (`float`) optional
    If provided ``<diag_name>.<species_name>.random_fraction = a``, only `a` fraction of the particle data of this species will be dumped randomly in diag ``<diag_name>``, i.e. if `rand() < a`, this particle will be dumped, where `rand()` denotes a random number generator.
    The value `a` provided should be between 0 and 1.

* ``<diag_name>.<species_name>.uniform_stride`` (`int`) optional
    If provided ``<diag_name>.<species_name>.uniform_stride = n``,
    every `n` particle of this species will be dumped, selected uniformly.
    The value provided should be an integer greater than or equal to 0.

* ``<diag_name>.<species_name>.plot_filter_function(t,x,y,z,ux,uy,uz)`` (`string`) optional
    Users can provide an expression returning a boolean for whether a particle is dumped (the exact test is whether the return value is `> 0.5`).
    `t` represents the physical time in seconds during the simulation.
    `x, y, z` represent particle positions in the unit of meter.
    `ux, uy, uz` represent particle velocities in the unit of
    :math:`\gamma v/c`, where
    :math:`\gamma` is the Lorentz factor,
    :math:`v/c` is the particle velocity normalized by the speed of light.
    E.g. If provided `(x>0.0)*(uz<10.0)` only those particles located at
    positions `x` greater than `0`, and those having velocity `uz` less than 10,
    will be dumped.

* ``amrex.async_out`` (`0` or `1`) optional (default `0`)
    Whether to use asynchronous IO when writing plotfiles. This only has an effect
    when using the AMReX plotfile format. Please see :doc:`../visualization/visualization`
    for more information.

* ``amrex.async_out_nfiles`` (`int`) optional (default `64`)
    The maximum number of files to write to when using asynchronous IO.
    To use asynchronous IO with more than ``amrex.async_out_nfiles`` MPI ranks,
    WarpX must be compiled with the ``MPI_THREAD_MULTIPLE=TRUE`` flag.
    Please see :doc:`../visualization/visualization` for more information.

Back-Transformed Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``BackTransformedDiagnostics`` are used when running a simulation in a boosted frame, to reconstruct output data to the lab frame, and

* ``warpx.do_back_transformed_diagnostics`` (`0` or `1`)
    Whether to use the **back-transformed diagnostics** (i.e. diagnostics that
    perform on-the-fly conversion to the laboratory frame, when running
    boosted-frame simulations)

* ``warpx.lab_data_directory`` (`string`)
    The directory in which to save the lab frame data when using the
    **back-transformed diagnostics**. If not specified, the default is
    is `lab_frame_data`.

* ``warpx.num_snapshots_lab`` (`integer`)
    Only used when ``warpx.do_back_transformed_diagnostics`` is ``1``.
    The number of lab-frame snapshots that will be written.

* ``warpx.dt_snapshots_lab`` (`float`, in seconds)
    Only used when ``warpx.do_back_transformed_diagnostics`` is ``1``.
    The time interval inbetween the lab-frame snapshots (where this
    time interval is expressed in the laboratory frame).

* ``warpx.dz_snapshots_lab`` (`float`, in meters)
    Only used when ``warpx.do_back_transformed_diagnostics`` is ``1``.
    Distance between the lab-frame snapshots (expressed in the laboratory
    frame). ``dt_snapshots_lab`` is then computed by
    ``dt_snapshots_lab = dz_snapshots_lab/c``. Either `dt_snapshots_lab`
    or `dz_snapshot_lab` is required.

* ``warpx.do_back_transformed_fields`` (`0 or 1`)
    Whether to use the **back-transformed diagnostics** for the fields.

* ``warpx.back_transformed_diag_fields`` (space-separated list of `string`)
    Which fields to dumped in back-transformed diagnostics. Choices are
    'Ex', 'Ey', Ez', 'Bx', 'By', Bz', 'jx', 'jy', jz' and 'rho'. Example:
    ``warpx.back_transformed_diag_fields = Ex Ez By``. By default, all fields
    are dumped.

* ``slice.num_slice_snapshots_lab`` (`integer`)
    Only used when ``warpx.do_back_transformed_diagnostics`` is ``1``.
    The number of back-transformed field and particle data that
    will be written for the reduced domain defined by ``slice.dom_lo``
    and ``slice.dom_hi``. Note that the 'slice' is a reduced
    diagnostic which could be 1D, 2D, or 3D, aligned with the co-ordinate axes.
    These slices can be visualized using read_raw_data.py and the HDF5 format can
    be visualized using the h5py library. Please see the documentation on visualization
    for further details.

* ``slice.dt_slice_snapshots_lab`` (`float`, in seconds)
    Only used when ``warpx.do_back_transformed_diagnostics`` is ``1``.
    The time interval between the back-transformed reduced diagnostics (where this
    time interval is expressed in the laboratory frame).

* ``slice.particle_slice_width_lab`` (`float`, in meters)
    Only used when ``warpx.do_back_transformed_diagnostics`` is ``1`` and
    ``slice.num_slice_snapshots_lab`` is non-zero. Particles are
    copied from the full back-transformed diagnostic to the reduced
    slice diagnostic if there are within the user-defined width from
    the slice region defined by ``slice.dom_lo`` and ``slice.dom_hi``.


Reduced Diagnostics
^^^^^^^^^^^^^^^^^^^

``ReducedDiags`` allow the user to compute some reduced quantity (particle temperature, max of a field) and write a small amount of data to text files.

* ``warpx.reduced_diags_names`` (`strings`, separated by spaces)
    The names given by the user of simple reduced diagnostics.
    Also the names of the output `.txt` files.
    This reduced diagnostics aims to produce simple outputs
    of the time history of some physical quantities.
    If ``warpx.reduced_diags_names`` is not provided in the input file,
    no reduced diagnostics will be done.
    This is then used in the rest of the input deck;
    in this documentation we use `<reduced_diags_name>` as a placeholder.

* ``<reduced_diags_name>.type`` (`string`)
    The type of reduced diagnostics associated with this `<reduced_diags_name>`.
    For example, ``ParticleEnergy`` and ``FieldEnergy``.
    All available types will be described below in detail.
    For all reduced diagnostics,
    the first and the second columns in the output file are
    the time step and the corresponding physical time in seconds, respectively.

    * ``ParticleEnergy``
        This type computes both the total and the mean
        relativistic particle kinetic energy among all species.

        .. math::

            E_p = \sum_{i=1}^N ( \sqrt{ p_i^2 c^2 + m_0^2 c^4 } - m_0 c^2 ) w_i

        where :math:`p` is the relativistic momentum,
        :math:`c` is the speed of light,
        :math:`m_0` is the rest mass,
        :math:`N` is the number of particles,
        :math:`w` is the individual particle weight.

        The output columns are
        total :math:`E_p` of all species,
        :math:`E_p` of each species,
        total mean energy :math:`E_p / \sum w_i`,
        mean energy of each species.

    * ``FieldEnergy``
        This type computes the electric and magnetic field energy.

        .. math::

            E_f = \sum [ \varepsilon_0 E^2 / 2 + B^2 / ( 2 \mu_0 ) ] \Delta V

        where
        :math:`E` is the electric field,
        :math:`B` is the magnetic field,
        :math:`\varepsilon_0` is the vacuum permittivity,
        :math:`\mu_0` is the vacuum permeability,
        :math:`\Delta V` is the cell volume (or area for 2D),
        the sum is over all cells.

        The output columns are
        total field energy :math:`E_f`,
        :math:`E` field energy,
        :math:`B` field energy, at mesh refinement levels from 0 to :math:`n`.

    * ``FieldMaximum``
        This type computes the maximum value of each component of the electric and magnetic fields
        and of the norm of the electric and magnetic field vectors.
        Measuring maximum fields in a plasma might be very noisy in PIC, use this instead
        for analysis of scenarios such as an electromagnetic wave propagating in vacuum.

        The output columns are
        the maximum value of the :math:`E_x` field,
        the maximum value of the :math:`E_y` field,
        the maximum value of the :math:`E_z` field,
        the maximum value of the norm :math:`|E|` of the electric field,
        the maximum value of the :math:`B_x` field,
        the maximum value of the :math:`B_y` field,
        the maximum value of the :math:`B_z` field and
        the maximum value of the norm :math:`|B|` of the magnetic field,
        at mesh refinement levels from  0 to :math:`n`.

        Note that the fields are averaged on the cell centers before their maximum values are
        computed.

    * ``RhoMaximum``
        This type computes the maximum and minimum values of the total charge density as well as
        the maximum absolute value of the charge density of each charged species.
        Please be aware that measuring maximum charge densities might be very noisy in PIC simulations.

        The output columns are
        the maximum value of the :math:`rho` field,
        the minimum value of the :math:`rho` field,
        the maximum value of the absolute :math:`|rho|` field of each charged species.

        Note that the charge densities are averaged on the cell centers before their maximum values
        are computed.

    * ``ParticleNumber``
        This type computes the total number of macroparticles and of physical particles (i.e. the
        sum of their weights) in the whole simulation domain (for each species and summed over all
        species). It can be useful in particular for simulations with creation (ionization, QED
        processes) or removal (resampling) of particles.

        The output columns are
        total number of macroparticles summed over all species,
        total number of macroparticles of each species,
        sum of the particles' weight summed over all species,
        sum of the particles' weight of each species.

    * ``BeamRelevant``
        This type computes properties of a particle beam relevant for particle accelerators,
        like position, momentum, emittance, etc.

        ``<reduced_diags_name>.species`` must be provided,
        such that the diagnostics are done for this (beam-like) species only.

        The output columns (for 3D-XYZ) are the following, where the average is done over
        the whole species (typical usage: the particle beam is in a separate species):

        [1], [2], [3]: The mean values of beam positions (m)
        :math:`\langle x \rangle`, :math:`\langle y \rangle`,
        :math:`\langle z \rangle`.

        [4], [5], [6]: The mean values of beam relativistic momenta (kg m/s)
        :math:`\langle p_x \rangle`, :math:`\langle p_y \rangle`,
        :math:`\langle p_z \rangle`.

        [7]: The mean Lorentz factor :math:`\langle \gamma \rangle`.

        [8], [9], [10]: The RMS values of beam positions (m)
        :math:`\delta_x = \sqrt{ \langle (x - \langle x \rangle)^2 \rangle }`,
        :math:`\delta_y = \sqrt{ \langle (y - \langle y \rangle)^2 \rangle }`,
        :math:`\delta_z = \sqrt{ \langle (z - \langle z \rangle)^2 \rangle }`.

        [11], [12], [13]: The RMS values of beam relativistic momenta (kg m/s)
        :math:`\delta_{px} = \sqrt{ \langle (p_x - \langle p_x \rangle)^2 \rangle }`,
        :math:`\delta_{py} = \sqrt{ \langle (p_y - \langle p_y \rangle)^2 \rangle }`,
        :math:`\delta_{pz} = \sqrt{ \langle (p_z - \langle p_z \rangle)^2 \rangle }`.

        [14]: The RMS value of the Lorentz factor
        :math:`\sqrt{ \langle (\gamma - \langle \gamma \rangle)^2 \rangle }`.

        [15], [16], [17]: beam projected transverse RMS normalized emittance (m)
        :math:`\epsilon_x = \dfrac{1}{mc} \sqrt{\delta_x^2 \delta_{px}^2 -
        \Big\langle (x-\langle x \rangle) (p_x-\langle p_x \rangle) \Big\rangle^2}`,
        :math:`\epsilon_y = \dfrac{1}{mc} \sqrt{\delta_y^2 \delta_{py}^2 -
        \Big\langle (y-\langle y \rangle) (p_y-\langle p_y \rangle) \Big\rangle^2}`,
        :math:`\epsilon_z = \dfrac{1}{mc} \sqrt{\delta_z^2 \delta_{pz}^2 -
        \Big\langle (z-\langle z \rangle) (p_z-\langle p_z \rangle) \Big\rangle^2}`.

        [18]: The charge of the beam (C).

        For 2D-XZ,
        :math:`\langle y \rangle`,
        :math:`\delta_y`, and
        :math:`\epsilon_y` will not be outputed.

    * ``LoadBalanceCosts``
        This type computes the cost, used in load balancing, for each box on the domain.
        The cost :math:`c` is computed as

        .. math::

            c = n_{\text{particle}} \cdot w_{\text{particle}} + n_{\text{cell}} \cdot w_{\text{cell}},

        where
        :math:`n_{\text{particle}}` is the number of particles on the box,
        :math:`w_{\text{particle}}` is the particle cost weight factor (controlled by ``algo.costs_heuristic_particles_wt``),
        :math:`n_{\text{cell}}` is the number of cells on the box, and
        :math:`w_{\text{cell}}` is the cell cost weight factor (controlled by ``algo.costs_heuristic_cells_wt``).

    * ``ParticleHistogram``
        This type computes a user defined particle histogram.

        * ``<reduced_diags_name>.species`` (`string`)
            A species name must be provided,
            such that the diagnostics are done for this species.

        * ``<reduced_diags_name>.histogram_function(t,x,y,z,ux,uy,uz)`` (`string`)
            A histogram function must be provided.
            `t` represents the physical time in seconds during the simulation.
            `x, y, z` represent particle positions in the unit of meter.
            `ux, uy, uz` represent the particle velocities in the unit of
            :math:`\gamma v/c`, where
            :math:`\gamma` is the Lorentz factor,
            :math:`v/c` is the particle velocity normalized by the speed of light.
            E.g.
            ``x`` produces the position (density) distribution in `x`.
            ``ux`` produces the velocity distribution in `x`,
            ``sqrt(ux*ux+uy*uy+uz*uz)`` produces the speed distribution.
            The default value of the histogram without normalization is
            :math:`f = \sum\limits_{i=1}^N w_i`, where
            :math:`\sum\limits_{i=1}^N` is the sum over :math:`N` particles
            in that bin,
            :math:`w_i` denotes the weight of the ith particle.

        * ``<reduced_diags_name>.bin_number`` (`int` > 0)
            This is the number of bins used for the histogram.

        * ``<reduced_diags_name>.bin_max`` (`float`)
            This is the maximum value of the bins.

        * ``<reduced_diags_name>.bin_min`` (`float`)
            This is the minimum value of the bins.

        * ``<reduced_diags_name>.normalization`` (optional)
            This provides options to normalize the histogram:

            ``unity_particle_weight``
            uses unity particle weight to compute the histogram,
            such that the values of the histogram are
            the number of counted macroparticles in that bin,
            i.e.  :math:`f = \sum\limits_{i=1}^N 1`,
            :math:`N` is the number of particles in that bin.

            ``max_to_unity`` will normalize the histogram such that
            its maximum value is one.

            ``area_to_unity`` will normalize the histogram such that
            the area under the histogram is one,
            so the histogram is also the probability density function.

            If nothing is provided,
            the macroparticle weight will be used to compute
            the histogram, and no normalization will be done.

        The output columns are
        values of the 1st bin, the 2nd bin, ..., the nth bin.
        An example input file and a loading pything script of
        using the histogram reduced diagnostics
        are given in ``Examples/Tests/initial_distribution/``.

    * ``ParticleExtrema``
        This type computes the minimum and maxmium values of
        particle position, momentum, gamma, weight,
        and the :math:`\chi` parameter for QED species.

        ``<reduced_diags_name>.species`` must be provided,
        such that the diagnostics are done for this species only.

        The output columns are
        minimum and maximum position :math:`x`, :math:`y`, :math:`z`;
        minimum and maximum momentum :math:`p_x`, :math:`p_y`, :math:`p_z`;
        minimum and maximum gamma :math:`\gamma`;
        minimum and maximum weight :math:`w`;
        minimum and maximum :math:`\chi`.

        Note that when the QED parameter :math:`\chi` is computed,
        field gather is carried out at every output,
        so the time of the diagnostic may be long
        depending on the simulation size.

* ``<reduced_diags_name>.frequency`` (`string`) optional (default ``1``)
    Using the `Intervals Parser`_ syntax, this string defines the timesteps at which reduced
    diagnostics are written to file.

* ``<reduced_diags_name>.path`` (`string`) optional (default `./diags/reducedfiles/`)
    The path that the output file will be stored.

* ``<reduced_diags_name>.extension`` (`string`) optional (default `txt`)
    The extension of the output file.

* ``<reduced_diags_name>.separator`` (`string`) optional (default a `whitespace`)
    The separator between row values in the output file.
    The default separator is a whitespace.

Lookup tables and other settings for QED modules
-----------------------------------------------------------------------------

Lookup tables store pre-computed values for functions used by the QED modules.
**This feature requires to compile with QED=TRUE (and also with QED_TABLE_GEN=TRUE for table generation) **

* ``qed_bw.lookup_table_mode`` (`string`)
    There are three options to prepare the lookup table required by the Breit-Wheeler module:

    * ``builtin``:  a built-in table is used (Warning: the table gives reasonable results but its resolution
    is quite low).

    * ``generate``: a new table is generated. This option requires Boost math library
      (version >= 1.66) and to compile with ``QED_TABLE_GEN=TRUE``. All
      the following parameters must be specified (table 1 is used to evolve the optical depth
      of the photons, while table 2 is used for pair generation):

        * ``qed_bw.tab_dndt_chi_min`` (`float`): minimum chi parameter for lookup table 1 (
          used for the evolution of the optical depth of the photons)

        * ``qed_bw.tab_dndt_chi_max`` (`float`): maximum chi parameter for lookup table 1

        * ``qed_bw.tab_dndt_how_many`` (`int`): number of points to be used for lookup table 1

        * ``qed_bw.tab_pair_chi_min`` (`float`): minimum chi parameter for lookup table 2 (
          used for pair generation)

        * ``qed_bw.tab_pair_chi_max`` (`float`): maximum chi parameter for lookup table 2

        * ``qed_bw.tab_pair_chi_how_many`` (`int`): number of points to be used for chi axis in lookup table 2

        * ``qed_bw.tab_pair_frac_how_many`` (`int`): number of points to be used for the second axis in lookup table 2
          (the second axis is the ratio between the quantum parameter of the less energetic particle of the pair and the
          quantum parameter of the photon).

        * ``qed_bw.save_table_in`` (`string`): where to save the lookup table

    * ``load``: a lookup table is loaded from a pre-generated binary file. The following parameter
      must be specified:

        * ``qed_bw.load_table_from`` (`string`): name of the lookup table file to read from.

* ``qed_qs.lookup_table_mode`` (`string`)
    There are three options to prepare the lookup table required by the Quantum Synchrotron module:

    * ``builtin``: a built-in table is used (Warning: the table gives reasonable results but its resolution
    is quite low).

    * ``generate``: a new table is generated. This option requires Boost math library
      (version >= 1.66) and to compile with ``QED_TABLE_GEN=TRUE``. All
      the following parameters must be specified (table 1 is used to evolve the optical depth
      of the particles, while table 2 is used for photon emission):

        * ``qed_qs.tab_dndt_chi_min`` (`float`): minimum chi parameter for lookup table 1 (
          used for the evolution of the optical depth of electrons and positrons)

        * ``qed_qs.tab_dndt_chi_max`` (`float`): maximum chi parameter for lookup table 1

        * ``qed_qs.tab_dndt_how_many`` (`int`): number of points to be used for lookup table 1

        * ``qed_qs.tab_em_chi_min`` (`float`): minimum chi parameter for lookup table 2 (
          used for photon emission)

        * ``qed_qs.tab_em_chi_max`` (`float`): maximum chi parameter for lookup table 2

        * ``qed_qs.tab_em_chi_how_many`` (`int`): number of points to be used for chi axis in lookup table 2

        * ``qed_qs.tab_em_frac_how_many`` (`int`): number of points to be used for the second axis in lookup table 2
          (the second axis is the ratio between the quantum parameter of the photon and the
          quantum parameter of the charged particle).

        * ``qed_qs.tab_em_frac_min`` (`float`): minimum value to be considered for the second axis of lookup table 2

        * ``qed_bw.save_table_in`` (`string`): where to save the lookup table

    * ``load``: a lookup table is loaded from a pre-generated binary file. The following parameter
      must be specified:

        * ``qed_qs.load_table_from`` (`string`): name of the lookup table file to read from.

* ``qed_bw.chi_min`` (`float`): minimum chi parameter to be considered by the Breit-Wheeler engine
    (suggested value : 0.01)

* ``qed_qs.chi_min`` (`float`): minimum chi parameter to be considered by the Quantum Synchrotron engine
    (suggested value : 0.001)

* ``qed_qs.photon_creation_energy_threshold`` (`float`) optional (default `2`)
    Energy threshold for photon particle creation in `*me*c^2` units.

* ``warpx.do_qed_schwinger`` (`bool`) optional (default `0`)
    If this is 1, Schwinger electron-positron pairs can be generated in vacuum in the cells where the EM field is high enough.
    Activating the Schwinger process requires the code to be compiled with ``QED=TRUE`` and ``PICSAR``.
    If ``warpx.do_qed_schwinger = 1``, Schwinger product species must be specified with
    ``qed_schwinger.ele_product_species`` and ``qed_schwinger.pos_product_species``.
    **Note: implementation of this feature is in progress.**
    So far it requires ``warpx.do_nodal=1`` and does not support mesh refinement, cylindrical coordinates or single precision.

* ``qed_schwinger.ele_product_species`` (`string`)
    If Schwinger process is activated, an electron product species must be specified
    (the name of an existing electron species must be provided).

* ``qed_schwinger.pos_product_species`` (`string`)
    If Schwinger process is activated, a positron product species must be specified
    (the name of an existing positron species must be provided).

* ``qed_schwinger.y_size`` (`float`; in meters)
    If Schwinger process is activated with ``DIM=2D``, a transverse size must be specified.
    It is used to convert the pair production rate per unit volume into an actual number of created particles.
    This value should correspond to the typical transverse extent for which the EM field has a very high value
    (e.g. the beam waist for a focused laser beam).

* ``qed_schwinger.xmin,ymin,zmin`` and ``qed_schwinger.xmax,ymax,zmax`` (`float`) optional (default unlimited)
    When ``qed_schwinger.xmin`` and ``qed_schwinger.xmax`` are set, they delimit the region within
    which Schwinger pairs can be created.
    The same is applicable in the other directions.

* ``qed_schwinger.threshold_poisson_gaussian`` (`integer`) optional (default `25`)
    If the expected number of physical pairs created in a cell at a given timestep is smaller than this threshold,
    a Poisson distribution is used to draw the actual number of physical pairs created.
    Otherwise a Gaussian distribution is used.
    Note that, regardless of this parameter, the number of macroparticles created is at most one per cell
    per timestep per species (with a weight corresponding to the number of physical pairs created).

Checkpoints and restart
-----------------------
WarpX supports checkpoints/restart via AMReX.
The checkpoint capability can be turned with regular diagnostics: ``<diag_name>.format = checkpoint``.

* ``amr.restart`` (`string`)
    Name of the checkpoint file to restart from. Returns an error if the folder does not exist
    or if it is not properly formatted.

Intervals parser
----------------

WarpX can parse time step interval expressions of the form ``start:stop:period``, e.g.
``1:2:3, 4::, 5:6, :, ::10``.
A comma is used as a separator between groups of intervals, which we call slices.
The resulting time steps are the `union set <https://en.wikipedia.org/wiki/Union_(set_theory)>`_ of all given slices.
White spaces are ignored.
A single slice can have 0, 1 or 2 colons ``:``, just as `numpy slices <https://numpy.org/doc/stable/reference/generated/numpy.s_.html>`_, but with inclusive upper bound for ``stop``.

* For 0 colon the given value is the period

* For 1 colon the given string is of the type ``start:stop``

* For 2 colons the given string is of the type ``start:stop:period``

Any value that is not given is set to default.
Default is ``0`` for the start, ``std::numeric_limits<int>::max()`` for the stop and ``1`` for the
period.
For the 1 and 2 colon syntax, actually having the integers in the string is optional
(this means that ``::5``, ``100 ::10`` and ``100 :`` are all valid syntaxes).

**Examples**

* ``something_int = 50`` -> do something at timesteps 0, 50, 100, 150, etc.
  (equivalent to ``something_int = ::50``)

* ``something_int = 300:600:100`` -> do something at timesteps 300, 400, 500 and 600.

* ``something_int = 300::50`` -> do something at timesteps 300, 350, 400, 450, etc.

* ``something_int = 105:108,205:208`` -> do something at timesteps 105, 106, 107, 108,
  205, 206, 207 and 208. (equivalent to ``something_int = 105 : 108 : , 205 : 208 :``)

* ``something_int = :`` or  ``something_int = ::`` -> do something at every timestep.

* ``something_int = 167:167,253:253,275:425:50`` do something at timesteps 167, 253, 275,
  325, 375 and 425.

This is essentially the python slicing syntax except that the stop is inclusive
(``0:100`` contains 100) and that no colon means that the given value is the period.

Note that if a given period is zero or negative, the correspoding slice is disregarded.
For example, ``something_int = -1`` deactivates ``something`` and
``something_int = ::-1,100:1000:25`` is equivalent to ``something_int = 100:1000:25``.
