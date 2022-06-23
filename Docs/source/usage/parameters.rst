.. _running-cpp-parameters:

Input Parameters
================

.. note::

   The AMReX parser (see :ref:`running-cpp-parameters-parser`) is used for the right-hand-side of all input parameters that consist of one or more integers or floats, so expressions like ``<species_name>.density_max = "2.+1."`` and/or using user-defined constants are accepted.

.. _running-cpp-parameters-overall:

Overall simulation parameters
-----------------------------

* ``authors`` (`string`: e.g. ``"Jane Doe <jane@example.com>, Jimmy Joe <jimmy@example.com>"``)
    Authors of an input file / simulation setup.
    When provided, this information is added as metadata to (openPMD) output files.

* ``max_step`` (`integer`)
    The number of PIC cycles to perform.

* ``stop_time`` (`float`; in seconds)
    The maximum physical time of the simulation. Can be provided instead of ``max_step``. If both
    ``max_step`` and ``stop_time`` are provided, both criteria are used and the simulation stops
    when the first criterion is hit.

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

* ``warpx.always_warn_immediately`` (``0`` or ``1``; default is ``0`` for false)
    If set to ``1``, WarpX immediately prints every warning message as soon as
    it is generated. It is mainly intended for debug purposes, in case a simulation
    crashes before a global warning report can be printed.

* ``warpx.abort_on_warning_threshold`` (string: ``low``, ``medium`` or ``high``) optional
    Optional threshold to abort as soon as a warning is raised.
    If the threshold is set, warning messages with priority greater than or
    equal to the threshold trigger an immediate abort.
    It is mainly intended for debug purposes, and is best used with
    ``warpx.always_warn_immediately=1``.

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
      separately taking into account their averaged velocities. The field
      is mapped to the simulation frame and will produce both E and B
      fields.

    See the `AMReX documentation <https://amrex-codes.github.io/amrex/docs_html/LinearSolvers.html#>`_
    for details of the MLMG solver (the default solver used with electrostatic
    simulations). The default behavior of the code is to check whether there is
    non-zero charge density in the system and if so force the MLMG solver to
    use the solution max norm when checking convergence. If there is no charge
    density, the MLMG solver will switch to using the initial guess max norm
    error when evaluating convergence and an absolute error tolerance of
    :math:`10^{-6}` :math:`\mathrm{V/m}^2` will be used (unless a different
    non-zero value is specified by the user via
    ``warpx.self_fields_absolute_tolerance``).

* ``warpx.self_fields_required_precision`` (`float`, default: 1.e-11)
    The relative precision with which the electrostatic space-charge fields should
    be calculated. More specifically, the space-charge fields are
    computed with an iterative Multi-Level Multi-Grid (MLMG) solver.
    This solver can fail to reach the default precision within a reasonable time.
    This only applies when warpx.do_electrostatic = labframe.

* ``warpx.self_fields_absolute_tolerance`` (`float`, default: 0.0)
    The absolute tolerance with which the space-charge fields should be
    calculated in units of :math:`\mathrm{V/m}^2`. More specifically, the acceptable
    residual with which the solution can be considered converged. In general
    this should be left as the default, but in cases where the simulation state
    changes very little between steps it can occur that the initial guess for
    the MLMG solver is so close to the converged value that it fails to improve
    that solution sufficiently to reach the ``self_fields_required_precision``
    value.

* ``warpx.self_fields_max_iters`` (`integer`, default: 200)
    Maximum number of iterations used for MLMG solver for space-charge
    fields calculation. In case if MLMG converges but fails to reach the desired
    ``self_fields_required_precision``, this parameter may be increased.
    This only applies when warpx.do_electrostatic = labframe.

* ``warpx.self_fields_verbosity`` (`integer`, default: 2)
    The vebosity used for MLMG solver for space-charge fields calculation. Currently
    MLMG solver looks for verbosity levels from 0-5. A higher number results in more
    verbose output.

* ``amrex.abort_on_out_of_gpu_memory``  (``0`` or ``1``; default is ``1`` for true)
    When running on GPUs, memory that does not fit on the device will be automatically swapped to host memory when this option is set to ``0``.
    This will cause severe performance drops.
    Note that even with this set to ``1`` WarpX will not catch all out-of-memory events yet when operating close to maximum device memory.
    `Please also see the documentation in AMReX <https://amrex-codes.github.io/amrex/docs_html/GPU.html#inputs-parameters>`_.

Signal Handling
^^^^^^^^^^^^^^^

WarpX can handle Unix (Linux/macOS) `process signals <https://en.wikipedia.org/wiki/Signal_(IPC)>`__.
This can be useful to configure jobs on HPC and cloud systems to shut down cleanly when they are close to reaching their allocated walltime or to steer the simulation behavior interactively.

Allowed signal names are documented in the `C++ standard <https://en.cppreference.com/w/cpp/utility/program/SIG_types>`__ and `POSIX <https://pubs.opengroup.org/onlinepubs/9699919799/basedefs/signal.h.html>`__.
We follow the same naming, but remove the ``SIG`` prefix, e.g., the WarpX signal configuration name for ``SIGINT`` is ``INT``.

* ``warpx.break_signals`` (array of `string`, separated by spaces) optional
    A list of signal names or numbers that the simulation should
    handle by cleanly terminating at the next timestep

* ``warpx.checkpoint_signals`` (array of `string`, separated by spaces) optional
    A list of signal names or numbers that the simulation should
    handle by outputting a checkpoint at the next timestep. A
    diagnostic of type `checkpoint` must be configured.

.. note::

   Certain signals are only available on specific platforms, please see the links above for details.
   Typically supported on Linux and macOS are ``HUP``, ``INT``, ``QUIT``, ``ABRT``, ``USR1``, ``USR2``, ``TERM``, ``TSTP``, ``URG``, and ``IO`` among others.

   Signals to think about twice before overwriting in *interactive simulations*:
   Note that ``INT`` (interupt) is the signal that ``Ctrl+C`` sends on the terminal, which most people use to abort a process; once overwritten you need to abort interactive jobs with, e.g., ``Ctrl+\`` (``QUIT``) or sending the ``KILL`` signal.
   The ``TSTP`` (terminal stop) command is sent interactively from ``Ctrl+Z`` to temporarily send a process to sleep (until send in the background with commands such as ``bg`` or continued with ``fg``), overwriting it would thus disable that functionality.
   The signals ``KILL`` and ``STOP`` cannot be used.

   The ``FPE`` signal should not be overwritten in WarpX, as it is `controlled by AMReX <https://amrex-codes.github.io/amrex/docs_html/Debugging.html#breaking-into-debuggers>`__ for :ref:`debug workflows that catch invalid floating-point operations <debugging_warpx>`.


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

* ``amr.ref_ratio_vect`` (`3 integers for x,y,z per refined level`)
    When using mesh refinement, this can be used to set the refinement ratio per direction and level, relative to the previous level.

    Example: for three levels, a value of ``2 2 4 8 8 16`` refines the first level by 2-fold in x and y and 4-fold in z compared to the coarsest level (level 0/mother grid); compared to the first level, the second level is refined 8-fold in x and y and 16-fold in z.

    Note: in development; currently allowed value: ``2 2 2``.

* ``geometry.dims`` (`string`)
    The dimensions of the simulation geometry.
    Supported values are ``1``, ``2``, ``3``, ``RZ``.
    For ``3``, a cartesian geometry of ``x``, ``y``, ``z`` is modeled.
    For ``2``, the axes are ``x`` and ``z`` and all physics in ``y`` is assumed to be translation symmetric.
    For ``1``, the only axis is ``z`` and the dimensions ``x`` and ``y`` are translation symmetric.
    For ``RZ``, we apply an azimuthal mode decomposition, with ``warpx.n_rz_azimuthal_modes`` providing further control.

    Note that this value has to match the :ref:`WarpX_DIMS <building-cmake-options>` compile-time option.
    If you installed WarpX from a :ref:`package manager <install-users>`, then pick the right executable by name.

* ``geometry.n_rz_azimuthal_modes`` (`integer`; 1 by default)
    When using the RZ version, this is the number of azimuthal modes.
    The default is ``1``, which corresponds to a perfectly axisymmetric simulation.

* ``geometry.prob_lo`` and ``geometry.prob_hi`` (`2 floats in 2D`, `3 floats in 3D`; in meters)
    The extent of the full simulation box. This box is rectangular, and thus its
    extent is given here by the coordinates of the lower corner (``geometry.prob_lo``) and
    upper corner (``geometry.prob_hi``). The first axis of the coordinates is x
    (or r with cylindrical) and the last is z.

* ``warpx.do_moving_window`` (`integer`; 0 by default)
    Whether to use a moving window for the simulation

* ``warpx.moving_window_dir`` (either ``x``, ``y`` or ``z``)
    The direction of the moving window.

* ``warpx.moving_window_v`` (`float`)
    The speed of moving window, in units of the speed of light
    (i.e. use ``1.0`` for a moving window that moves exactly at the speed of light)

* ``warpx.start_moving_window_step`` (`integer`; 0 by default)
    The timestep at which the moving window starts.

* ``warpx.end_moving_window_step`` (`integer`; default is ``-1`` for false)
    The timestep at which the moving window ends.

* ``warpx.fine_tag_lo`` and ``warpx.fine_tag_hi`` (`2 floats in 2D`, `3 floats in 3D`; in meters) optional
    **When using static mesh refinement with 1 level**, the extent of the refined patch.
    This patch is rectangular, and thus its extent is given here by the coordinates
    of the lower corner (``warpx.fine_tag_lo``) and upper corner (``warpx.fine_tag_hi``).

* ``warpx.refine_plasma`` (`integer`) optional (default `0`)
    Increase the number of macro-particles that are injected "ahead" of a mesh
    refinement patch in a moving window simulation.

    Note: in development; only works with static mesh-refinement, specific
    to moving window plasma injection, and requires a single refined level.

* ``warpx.n_current_deposition_buffer`` (`integer`)
    When using mesh refinement: the particles that are located inside
    a refinement patch, but within ``n_current_deposition_buffer`` cells of
    the edge of this patch, will deposit their charge and current to the
    lower refinement level, instead of depositing to the refinement patch
    itself. See the :ref:`mesh-refinement section <theory-amr>` for more details.
    If this variable is not explicitly set in the input script,
    ``n_current_deposition_buffer`` is automatically set so as to be large
    enough to hold the particle shape, on the fine grid

* ``warpx.n_field_gather_buffer`` (`integer`; 0 by default)
    When using mesh refinement: the particles that are located inside
    a refinement patch, but within ``n_field_gather_buffer`` cells of
    the edge of this patch, will gather the fields from the lower refinement
    level, instead of gathering the fields from the refinement patch itself.
    This avoids some of the spurious effects that can occur inside the
    refinement patch, close to its edge. See the
    :ref:`mesh-refinement section <theory-amr>` for more details.
    If this variable is not
    explicitly set in the input script, ``n_field_gather_buffer`` is
    automatically set so that it is one cell larger than
    ``n_current_deposition_buffer``, on the fine grid.

* ``warpx.do_single_precision_comms`` (`integer`; 0 by default)
    Perform MPI communications for field guard regions in single precision.
    Only meaningful for ``WarpX_PRECISION=DOUBLE``.

* ``particles.deposit_on_main_grid`` (`list of strings`)
    When using mesh refinement: the particle species whose name are included
    in the list will deposit their charge/current directly on the main grid
    (i.e. the coarsest level), even if they are inside a refinement patch.

* ``particles.gather_from_main_grid`` (`list of strings`)
    When using mesh refinement: the particle species whose name are included
    in the list will gather their fields from the main grid
    (i.e. the coarsest level), even if they are inside a refinement patch.

.. _running-cpp-parameters-bc:

Domain Boundary Conditions
--------------------------

* ``boundary.field_lo`` and ``boundary.field_hi`` (`2 strings` for 2D, `3 strings` for 3D, `pml` by default)
    Boundary conditions applied to fields at the lower and upper domain boundaries.
    Options are:

    * ``Periodic``: This option can be used to set periodic domain boundaries. Note that if the fields for lo in a certain dimension are set to periodic, then the corresponding upper boundary must also be set to periodic. If particle boundaries are not specified in the input file, then particles boundaries by default will be set to periodic. If particles boundaries are specified, then they must be set to periodic corresponding to the periodic field boundaries.

    * ``pml`` (default): This option can be used to add Perfectly Matched Layers (PML) around the simulation domain. See the :ref:`PML theory section <theory-bc>` for more details.
      Additional pml algorithms can be explored using the parameters ``warpx.do_pml_in_domain``, ``warpx.pml_has_particles``, and ``warpx.do_pml_j_damping``.

    * ``absorbing_silver_mueller``: This option can be used to set the Silver-Mueller absorbing boundary conditions. These boundary conditions are simpler and less computationally expensive than the pml, but are also less effective at absorbing the field. They only work with the Yee Maxwell solver.

    * ``damped``: This is the recommended option in the moving direction when using the spectral solver with moving window (currently only supported along z). This boundary condition applies a damping factor to the electric and magnetic fields in the outer half of the guard cells, using a sine squared profile. As the spectral solver is by nature periodic, the damping prevents fields from wrapping around to the other end of the domain when the periodicity is not desired. This boundary condition is only valid when using the spectral solver.

    * ``pec``: This option can be used to set a Perfect Electric Conductor at the simulation boundary. For the electromagnetic solve, at PEC, the tangential electric field and the normal magnetic field are set to 0. This boundary can be used to model a dielectric or metallic surface. In the guard-cell region, the tangential electric field is set equal and opposite to the respective field component in the mirror location across the PEC boundary, and the normal electric field is set equal to the field component in the mirror location in the domain across the PEC boundary. Similarly, the tangential (and normal) magnetic field components are set equal (and opposite) to the respective magnetic field components in the mirror locations across the PEC boundary. Note that PEC boundary is invalid at `r=0` for the RZ solver. Please use ``none`` option. This boundary condition does not work with the spectral solver.
      If an electrostatic field solve is used the boundary potentials can also be set through ``boundary.potential_lo_x/y/z`` and ``boundary.potential_hi_x/y/z`` (default `0`).

    * ``none``: No boundary condition is applied to the fields with the electromagnetic solver. This option must be used for the RZ-solver at `r=0`. If the electrostatic solver is used, a Neumann boundary condition (with gradient equal to 0) will be applied on the specified boundary.

* ``boundary.particle_lo`` and ``boundary.particle_hi`` (`2 strings` for 2D, `3 strings` for 3D, `absorbing` by default)
    Options are:

    * ``Absorbing``: Particles leaving the boundary will be deleted.

    * ``Periodic``: Particles leaving the boundary will re-enter from the opposite boundary. The field boundary condition must be consistenly set to periodic and both lower and upper boundaries must be periodic.

    * ``Reflecting``: Particles leaving the boundary are reflected from the boundary back into the domain.
      When ``boundary.reflect_all_velocities`` is false, the sign of only the normal velocity is changed, otherwise the sign of all velocities are changed.

* ``boundary.reflect_all_velocities`` (`bool`) optional (default `false`)
    For a reflecting boundary condition, this flags whether the sign of only the normal velocity is changed or all velocities.


Additional PML parameters
-------------------------

* ``warpx.pml_ncell`` (`int`; default: 10)
    The depth of the PML, in number of cells.

* ``do_similar_dm_pml`` (`int`; default: 1)
    Whether or not to use an amrex::DistributionMapping for the PML grids that is `similar` to the mother grids, meaning that the
    mapping will be computed to minimize the communication costs between the PML and the mother grids.

* ``warpx.pml_delta`` (`int`; default: 10)
    The characteristic depth, in number of cells, over which
    the absorption coefficients of the PML increases.

* ``warpx.do_pml_in_domain`` (`int`; default: 0)
    Whether to create the PML inside the simulation area or outside. If inside,
    it allows the user to propagate particles in PML and to use extended PML

* ``warpx.pml_has_particles`` (`int`; default: 0)
    Whether to propagate particles in PML or not. Can only be done if PML are in simulation domain,
    i.e. if `warpx.do_pml_in_domain = 1`.

* ``warpx.do_pml_j_damping`` (`int`; default: 0)
    Whether to damp current in PML. Can only be used if particles are propagated in PML,
    i.e. if `warpx.pml_has_particles = 1`.

* ``warpx.v_particle_pml`` (`float`; default: 1)
    When ``warpx.do_pml_j_damping = 1``, the assumed velocity of the particles to be absorbed in the PML, in units of the speed of light `c`.

* ``warpx.do_pml_dive_cleaning`` (`bool`; default: 1)
    Whether to use divergence cleaning for E in the PML region.
    The value must match ``warpx.do_pml_divb_cleaning`` (either both false or both true).
    This option seems to be necessary in order to avoid strong Nyquist instabilities in 3D simulations with the PSATD solver, open boundary conditions and PML in all directions. 2D simulations and 3D simulations with open boundary conditions and PML only in one direction might run well even without divergence cleaning.
    This option is implemented only for the PSATD solver.

* ``warpx.do_pml_divb_cleaning`` (`bool`; default: 1)
    Whether to use divergence cleaning for B in the PML region.
    The value must match ``warpx.do_pml_dive_cleaning`` (either both false or both true).
    This option seems to be necessary in order to avoid strong Nyquist instabilities in 3D simulations with the PSATD solver, open boundary conditions and PML in all directions. 2D simulations and 3D simulations with open boundary conditions and PML only in one direction might run well even without divergence cleaning.

.. _running-cpp-parameters-eb:

Embedded Boundary Conditions
----------------------------

* ``warpx.eb_implicit_function`` (`string`)
    A function of `x`, `y`, `z` that defines the surface of the embedded
    boundary. That surface lies where the function value is 0 ;
    the physics simulation area is where the function value is negative ;
    the interior of the embeddded boundary is where the function value is positive.

* ``warpx.eb_potential(x,y,z,t)`` (`string`)
    Only used when ``warpx.do_electrostatic=labframe``. Gives the value of
    the electric potential at the surface of the embedded boundary,
    as a function of  `x`, `y`, `z` and time. This function is also evaluated
    inside the embedded boundary. For this reason, it is important to define
    this function in such a way that it is constant inside the embedded boundary.

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

* ``algo.load_balance_intervals`` (`string`) optional (default `0`)
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

* ``algo.load_balance_efficiency_ratio_threshold`` (`float`) optional (default `1.1`)
    Controls whether to adopt a proposed distribution mapping computed during a load balance.
    If the the ratio of the proposed to current distribution mapping *efficiency* (i.e.,
    average cost per MPI process; efficiency is a number in the range [0, 1]) is greater
    than the threshold value, the proposed distribution mapping is adopted.  The suggested
    range of values is ``algo.load_balance_efficiency_ratio_threshold >= 1``, which ensures
    that the new distribution mapping is adopted only if doing so would improve the load
    balance efficiency. The higher the threshold value, the more conservative is the criterion
    for adoption of a proposed distribution; for example, with
    ``algo.load_balance_efficiency_ratio_threshold = 1``, the proposed distribution is
    adopted *any* time the proposed distribution improves load balancing; if instead
    ``algo.load_balance_efficiency_ratio_threshold = 2``, the proposed distribution is
    adopted only if doing so would yield a 100% to the load balance efficiency (with this
    threshold value, if the  current efficiency is ``0.45``, the new distribution would only be
    adopted if the proposed efficiency were greater than ``0.9``).

* ``algo.load_balance_with_sfc`` (`0` or `1`) optional (default `0`)
    If this is `1`: use a Space-Filling Curve (SFC) algorithm in order to
    perform load-balancing of the simulation.
    If this is `0`: the Knapsack algorithm is used instead.

* ``algo.load_balance_knapsack_factor`` (`float`) optional (default `1.24`)
    Controls the maximum number of boxes that can be assigned to a rank during
    load balance when using the 'knapsack' policy for update of the distribution
    mapping; the maximum is
    `load_balance_knapsack_factor*(average number of boxes per rank)`.
    For example, if there are 4 boxes per rank and `load_balance_knapsack_factor=2`,
    no more than 8 boxes can be assigned to any rank.

* ``algo.load_balance_costs_update`` (`heuristic` or `timers` or `gpuclock`) optional (default `timers`)
    If this is `heuristic`: load balance costs are updated according to a measure of
    particles and cells assigned to each box of the domain.  The cost :math:`c` is
    computed as

    .. math::

            c = n_{\text{particle}} \cdot w_{\text{particle}} + n_{\text{cell}} \cdot w_{\text{cell}},

    where
    :math:`n_{\text{particle}}` is the number of particles on the box,
    :math:`w_{\text{particle}}` is the particle cost weight factor (controlled by ``algo.costs_heuristic_particles_wt``),
    :math:`n_{\text{cell}}` is the number of cells on the box, and
    :math:`w_{\text{cell}}` is the cell cost weight factor (controlled by ``algo.costs_heuristic_cells_wt``).

    If this is `timers`: costs are updated according to in-code timers.

    If this is `gpuclock`: [**requires to compile with option** ``-DWarpX_GPUCLOCK=ON``]
    costs are measured as (max-over-threads) time spent in current deposition
    routine (only applies when running on GPUs).

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

WarpX uses AMReX's math parser that reads expressions in the input file.
It can be used in all input parameters that consist of one or more integers or floats.
Integer input expecting boolean, 0 or 1, are not parsed.
Note that when multiple values are expected, the expressions are space delimited.
For integer input values, the expressions are evaluated as real numbers and the final result rounded to the nearest integer.
See `this section <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parser>`__ of the AMReX documentation for a complete list of functions supported by the math parser.

WarpX constants
^^^^^^^^^^^^^^^

WarpX provides a few pre-defined constants, that can be used for any parameter that consists of one or more floats.

======== ===================
q_e      elementary charge
m_e      electron mass
m_p      proton mass
m_u      unified atomic mass unit (Dalton)
epsilon0 vacuum permittivity
mu0      vacuum permeability
clight   speed of light
kb       Boltzmann's constant (J/K)
pi       math constant pi
======== ===================

See ``Source/Utils/WarpXConst.H`` for the values.

User-defined constants
^^^^^^^^^^^^^^^^^^^^^^

Users can define their own constants in the input file.
These constants can be used for any parameter that consists of one or more integers or floats.
User-defined constant names can contain only letters, numbers and the character ``_``.
The name of each constant has to begin with a letter. The following names are used
by WarpX, and cannot be used as user-defined constants: ``x``, ``y``, ``z``, ``X``, ``Y``, ``t``.
The values of the constants can include the predefined WarpX constants listed above as well as other user-defined constants.
For example:

* ``my_constants.a0 = 3.0``
* ``my_constants.z_plateau = 150.e-6``
* ``my_constants.n0 = 1.e22``
* ``my_constants.wp = sqrt(n0*q_e**2/(epsilon0*m_e))``

Coordinates
^^^^^^^^^^^

Besides, for profiles that depend on spatial coordinates (the plasma momentum distribution or the laser field, see below `Particle initialization` and `Laser initialization`), the parser will interpret some variables as spatial coordinates. These are specified in the input parameter, i.e., ``density_function(x,y,z)`` and ``field_function(X,Y,t)``.

The parser reads python-style expressions between double quotes, for instance
``"a0*x**2 * (1-y*1.e2) * (x>0)"`` is a valid expression where ``a0`` is a
user-defined constant (see above) and ``x`` and ``y`` are spatial coordinates. The names are case sensitive. The factor
``(x>0)`` is ``1`` where ``x>0`` and ``0`` where ``x<=0``. It allows the user to
define functions by intervals.
Alternatively the expression above can be written as ``if(x>0, a0*x**2 * (1-y*1.e2), 0)``.

.. _running-cpp-parameters-particle:

Particle initialization
-----------------------

* ``particles.species_names`` (`strings`, separated by spaces)
    The name of each species. This is then used in the rest of the input deck ;
    in this documentation we use `<species_name>` as a placeholder.

* ``particles.use_fdtd_nci_corr`` (`0` or `1`) optional (default `0`)
    Whether to activate the FDTD Numerical Cherenkov Instability corrector.
    Not currently available in the RZ configuration.

* ``particles.rigid_injected_species`` (`strings`, separated by spaces)
    List of species injected using the rigid injection method. The rigid injection
    method is useful when injecting a relativistic particle beam, in boosted-frame
    simulation ; see the :ref:`input-output section <theory-io>` for more details.
    For species injected using this method, particles are translated along the `+z`
    axis with constant velocity as long as their ``z`` coordinate verifies
    ``z<zinject_plane``. When ``z>zinject_plane``,
    particles are pushed in a standard way, using the specified pusher.
    (see the parameter ``<species_name>.zinject_plane`` below)

* ``particles.do_tiling`` (`bool`) optional (default `false` if WarpX is compiled for GPUs, `true` otherwise)
    Controls whether tiling ('cache blocking') transformation is used for particles.
    Tiling should be on when using OpenMP and off when using GPUs.

* ``<species_name>.species_type`` (`string`) optional (default `unspecified`)
    Type of physical species.
    Currently, the accepted species are ``"electron"``, ``"positron"``, ``"photon"``, ``"hydrogen"`` (or equivalently ``"proton"``), ``"helium"`` (or equivalently ``"alpha"``), ``"boron"``, ``"carbon"``, ``"oxygen"``, ``"nitrogen"``, ``"argon"``, ``"copper"`` and ``"xenon"``.
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

* ``<species_name>.injection_style`` (`string`; default: ``none``)
    Determines how the (macro-)particles will be injected in the simulation.
    The number of particles per cell is always given with respect to the coarsest level (level 0/mother grid), even if particles are immediately assigned to a refined patch.

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

    * ``MultipleParticles``: Inject multiple macroparticles.
      This requires the additional parameters:
      ``<species_name>.multiple_particles_pos_x`` (list of `doubles`, X positions of the particles [meter])
      ``<species_name>.multiple_particles_pos_y`` (list of `doubles`, Y positions of the particles [meter])
      ``<species_name>.multiple_particles_pos_z`` (list of `doubles`, Z positions of the particles [meter])
      ``<species_name>.multiple_particles_vel_x`` (list of `doubles`, X normalized momenta of the particles, i.e. :math:`\gamma \beta_x`)
      ``<species_name>.multiple_particles_vel_y`` (list of `doubles`, Y normalized momenta of the particles, i.e. :math:`\gamma \beta_y`)
      ``<species_name>.multiple_particles_vel_z`` (list of `doubles`, Z normalized momenta of the particles, i.e. :math:`\gamma \beta_z`)
      ``<species_name>.multiple_particles_weight`` (list of `doubles`, macroparticle weights, i.e. number of physical particles each represents)

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
      The external file must include the species ``openPMD::Record`` labeled ``position`` and ``momentum`` (`double` arrays), with dimensionality and units set via ``openPMD::setUnitDimension`` and ``setUnitSI``.
      If the external file also contains ``openPMD::Records`` for ``mass`` and ``charge`` (constant `double` scalars) then the species will use these, unless overwritten in the input file (see ``<species_name>.mass``, ``<species_name>.charge`` or ``<species_name>.species_type``).
      The ``external_file`` option is currently implemented for 2D, 3D and RZ geometries, with record components in the cartesian coordinates ``(x,y,z)`` for 3D and RZ, and ``(x,z)`` for 2D.
      For more information on the `openPMD format <https://github.com/openPMD>`__ and how to build WarpX with it, please visit :ref:`the install section <install-developers>`.

    * ``NFluxPerCell``: Continuously inject a flux of macroparticles from a planar surface.
      The density specified by the density profile is interpreted to have the units of #/m^2/s.
      This requires the additional parameters:
      ``<species_name>.surface_flux_pos`` (`double`, location of the injection plane [meter])
      ``<species_name>.flux_normal_axis`` (`x`, `y`, or `z` for 3D, `x` or `z` for 2D, or `r`, `t`, or `z` for RZ. When `flux_normal_axis` is `r` or `t`, the `x` and `y` components of the user-specified momentum distribution are interpreted as the `r` and `t` components respectively)
      ``<species_name>.flux_direction`` (`-1` or `+1`, direction of flux relative to the plane)
      ``<species_name>.num_particles_per_cell`` (`double`)
      ``<species_name>.flux_tmin`` (`double`, Optional time at which the flux will be turned on. Ignored when negative.)
      ``<species_name>.flux_tmax`` (`double`, Optional time at which the flux will be turned off. Ignored when negative.)

    * ``none``: Do not inject macro-particles (for example, in a simulation that starts with neutral, ionizable atoms, one may want to create the electrons species -- where ionized electrons can be stored later on -- without injecting electron macro-particles).

* ``<species_name>.num_particles_per_cell_each_dim`` (`3 integers in 3D and RZ, 2 integers in 2D`)
    With the NUniformPerCell injection style, this specifies the number of particles along each axis
    within a cell. Note that for RZ, the three axis are radius, theta, and z and that the recommended
    number of particles per theta is at least two times the number of azimuthal modes requested.
    (It is recommended to do a convergence scan of the number of particles per theta)

* ``<species_name>.random_theta`` (`bool`) optional (default `1`)
    When using RZ geometry, whether to randomize the azimuthal position of particles.
    This is used when ``<species_name>.injection_style = NUniformPerCell``.

* ``<species_name>.do_splitting`` (`bool`) optional (default `0`)
    Split particles of the species when crossing the boundary from a lower
    resolution domain to a higher resolution domain.

    Currently implemented on CPU only.

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

* ``<species_name>.self_fields_absolute_tolerance`` (`float`, default: 0.0)
    The absolute tolerance with which the space-charge fields should be
    calculated in units of :math:`\mathrm{V/m}^2`. More specifically, the acceptable
    residual with which the solution can be considered converged. In general
    this should be left as the default, but in cases where the simulation state
    changes very little between steps it can occur that the initial guess for
    the MLMG solver is so close to the converged value that it fails to improve
    that solution sufficiently to reach the ``self_fields_required_precision``
    value.

* ``<species_name>.self_fields_max_iters`` (`integer`, default: 200)
    Maximum number of iterations used for MLMG solver for initial space-charge
    fields calculation. In case if MLMG converges but fails to reach the desired
    ``self_fields_required_precision``, this parameter may be increased.

* ``<species_name>.profile`` (`string`)
    Density profile for this species. The options are:

    * ``constant``: Constant density profile within the box, or between ``<species_name>.xmin``
      and ``<species_name>.xmax`` (and same in all directions). This requires additional
      parameter ``<species_name>.density``. i.e., the plasma density in :math:`m^{-3}`.

    * ``predefined``: Predefined density profile.
      This requires additional parameters ``<species_name>.predefined_profile_name`` and ``<species_name>.predefined_profile_params``.
      Currently, only a parabolic channel density profile is implemented.

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

    * ``at_rest``: Particles are initialized with zero momentum.

    * ``constant``: constant momentum profile. This can be controlled with the additional parameters
      ``<species_name>.ux``, ``<species_name>.uy`` and ``<species_name>.uz``, the normalized
      momenta in the x, y and z direction respectively, which are all ``0.`` by default.

    * ``gaussian``: gaussian momentum distribution in all 3 directions. This can be controlled with the
      additional arguments for the average momenta along each direction
      ``<species_name>.ux_m``, ``<species_name>.uy_m`` and ``<species_name>.uz_m`` as
      well as standard deviations along each direction ``<species_name>.ux_th``,
      ``<species_name>.uy_th`` and ``<species_name>.uz_th``.
      These 6 parameters are all ``0.`` by default.

    * ``gaussianflux``: Gaussian momentum flux distribution, which is Gaussian in the plane and v*Gaussian normal to the plane.
      It can only be used when ``injection_style = NFluxPerCell``.
      This can be controlled with the additional arguments to specify the plane's orientation, ``<species_name>.flux_normal_axis`` and
      ``<species_name>.flux_direction``, for the average momenta along each direction
      ``<species_name>.ux_m``, ``<species_name>.uy_m`` and ``<species_name>.uz_m``, as
      well as standard deviations along each direction ``<species_name>.ux_th``,
      ``<species_name>.uy_th`` and ``<species_name>.uz_th``. Note that the average momenta normal to the plane is not used.
      ``ux_m``, ``uy_m``, ``uz_m``, ``ux_th``, ``uy_th`` and ``uz_th`` are all ``0.`` by default.

    * ``maxwell_boltzmann``: Maxwell-Boltzmann distribution that takes a dimensionless
      temperature parameter :math:`\theta` as an input, where :math:`\theta = \frac{k_\mathrm{B} \cdot T}{m \cdot c^2}`,
      "math:`T` is the temperature in Kelvin, :math:`k_\mathrm{B}` is the Boltzmann constant, :math:`c` is the speed of light, and :math:`m` is the mass of the species.
      Theta is specified by a combination of ``<species_name>.theta_distribution_type``, ``<species_name>.theta``, and ``<species_name>.theta_function(x,y,z)`` (see below).
      For values of :math:`\theta > 0.01`, errors due to ignored relativistic terms exceed 1%.
      Temperatures less than zero are not allowed.
      The plasma can be initialized to move at a bulk velocity :math:`\beta = v/c`.
      The speed is specified by the parameters ``<species_name>.beta_distribution_type``, ``<species_name>.beta``, and ``<species_name>.beta_function(x,y,z)`` (see below).
      :math:`\beta` can be positive or negative and is limited to the range :math:`-1 < \beta < 1`.
      The direction of the velocity field is given by ``<species_name>.bulk_vel_dir = (+/-) 'x', 'y', 'z'``, and must be the same across the domain.
      Please leave no whitespace
      between the sign and the character on input. A direction without a sign will be treated as
      positive. The MB distribution is initialized in the drifting frame by sampling three Gaussian
      distributions in each dimension using, the Box Mueller method, and then the distribution is
      transformed to the simulation frame using the flipping method. The flipping method can be
      found in Zenitani 2015 section III. B. (Phys. Plasmas 22, 042116).
      By default, ``beta`` is equal to ``0.`` and ``bulk_vel_dir`` is ``+x``.

      Note that though the particles may move at relativistic speeds in the simulation frame,
      they are not relativistic in the drift frame. This is as opposed to the Maxwell Juttner
      setting, which initializes particles with relativistic momentums in their drifting frame.

    * ``maxwell_juttner``: Maxwell-Juttner distribution for high temperature plasma that takes a dimensionless temperature parameter :math:`\theta` as an input, where :math:`\theta = \frac{k_\mathrm{B} \cdot T}{m \cdot c^2}`,
      :math:`T` is the temperature in Kelvin, :math:`k_\mathrm{B}` is the Boltzmann constant, and :math:`m` is the mass of the species.
      Theta is specified by a combination of ``<species_name>.theta_distribution_type``, ``<species_name>.theta``, and ``<species_name>.theta_function(x,y,z)`` (see below).
      The Sobol method used to generate the distribution will not terminate for :math:`\theta \lesssim 0.1`, and the code will abort if it encounters a temperature below that threshold.
      The Maxwell-Boltzmann distribution is recommended for temperatures in the range :math:`0.01 < \theta < 0.1`.
      Errors due to relativistic effects can be expected to approximately between 1% and 10%.
      The plasma can be initialized to move at a bulk velocity :math:`\beta = v/c`.
      The speed is specified by the parameters ``<species_name>.beta_distribution_type``, ``<species_name>.beta``, and ``<species_name>.beta_function(x,y,z)`` (see below).
      :math:`\beta` can be positive or negative and is limited to the range :math:`-1 < \beta < 1`.
      The direction of the velocity field is given by ``<species_name>.bulk_vel_dir = (+/-) 'x', 'y', 'z'``, and must be the same across the domain.
      Please leave no whitespace
      between the sign and the character on input. A direction without a sign will be treated as
      positive. The MJ distribution will be initialized in the moving frame using the Sobol method,
      and then the distribution will be transformed to the simulation frame using the flipping method.
      Both the Sobol and the flipping method can be found in Zenitani 2015 (Phys. Plasmas 22, 042116).
      By default, ``beta`` is equal to ``0.`` and ``bulk_vel_dir`` is ``+x``.

      Please take notice that particles initialized with this setting can be relativistic in two ways.
      In the simulation frame, they can drift with a relativistic speed beta. Then, in the drifting
      frame they are still moving with relativistic speeds due to high temperature. This is as opposed
      to the Maxwell Boltzmann setting, which initializes non-relativistic plasma in their relativistic
      drifting frame.

    * ``radial_expansion``: momentum depends on the radial coordinate linearly. This
      can be controlled with additional parameter ``u_over_r`` which is the slope (``0.`` by default).

    * ``parse_momentum_function``: the momentum is given by a function in the input
      file. It requires additional arguments ``<species_name>.momentum_function_ux(x,y,z)``,
      ``<species_name>.momentum_function_uy(x,y,z)`` and ``<species_name>.momentum_function_uz(x,y,z)``,
      which gives the distribution of each component of the momentum as a function of space.

* ``<species_name>.theta_distribution_type`` (`string`) optional (default ``constant``)
    Only read if ``<species_name>.momentum_distribution_type`` is ``maxwell_boltzmann`` or ``maxwell_juttner``.
    See documentation for these distributions (above) for constraints on values of theta. Temperatures less than zero are not allowed.

    * If ``constant``, use a constant temperature, given by the required float parameter ``<species_name>.theta``.

    * If ``parser``, use a spatially-dependent analytic parser function, given by the required parameter ``<species_name>.theta_function(x,y,z)``.

* ``<species_name>.beta_distribution_type`` (`string`) optional (default ``constant``)
    Only read if ``<species_name>.momentum_distribution_type`` is ``maxwell_boltzmann`` or ``maxwell_juttner``.
    See documentation for these distributions (above) for constraints on values of beta.

    * If ``constant``, use a constant speed, given by the required float parameter ``<species_name>.beta``.

    * If ``parser``, use a spatially-dependent analytic parser function, given by the required parameter ``<species_name>.beta_function(x,y,z)``.

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
    Only read if ``<species_name>.profile`` is ``predefined``.

    * If ``parabolic_channel``, the plasma profile is a parabolic profile with
      cosine-like ramps at the beginning and the end of the profile.
      The density is given by

      .. math::

          n = n_0 n(x,y) n(z-z_0)

      with

      .. math::

          n(x,y) = 1 + 4\frac{x^2+y^2}{k_p^2 R_c^4}

      where :math:`k_p` is the plasma wavenumber associated with density :math:`n_0`.
      Here, with :math:`z_0` as the start of the plasma, :math:`n(z-z_0)` is a cosine-like up-ramp from :math:`0` to :math:`L_{ramp,up}`,
      constant to :math:`1` from :math:`L_{ramp,up}` to :math:`L_{ramp,up} + L_{plateau}`
      and a cosine-like down-ramp from :math:`L_{ramp,up} + L_{plateau}` to
      :math:`L_{ramp,up} + L_{plateau}+L_{ramp,down}`. All parameters are given
      in ``predefined_profile_params``.

* ``<species_name>.predefined_profile_params`` (list of `float`)
    Parameters for the predefined profiles.

    * If ``species_name.predefined_profile_name`` is ``parabolic_channel``,
      ``predefined_profile_params`` contains a space-separated list of the
      following parameters, in this order: :math:`z_0` :math:`L_{ramp,up}` :math:`L_{plateau}`
      :math:`L_{ramp,down}` :math:`R_c` :math:`n_0`

* ``<species_name>.do_backward_propagation`` (`bool`)
    Inject a backward-propagating beam to reduce the effect of charge-separation
    fields when running in the boosted frame. See examples.

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

* ``<species_name>.addIntegerAttributes`` (list of `string`)
    User-defined integer particle attribute for species, ``species_name``.
    These integer attributes will be initialized with user-defined functions
    when the particles are generated.
    If the user-defined integer attribute is ``<int_attrib_name>`` then the
    following required parameter must be specified to initialize the attribute.
    * ``<species_name>.attribute.<int_attrib_name>(x,y,z,ux,uy,uz,t)`` (`string`)
    ``t`` represents the physical time in seconds during the simulation.
    ``x``, ``y``, ``z`` represent particle positions in the unit of meter.
    ``ux``, ``uy``, ``uz`` represent the particle velocities in the unit of
    :math:`\gamma v/c`, where
    :math:`\gamma` is the Lorentz factor,
    :math:`v/c` is the particle velocity normalized by the speed of light.
    E.g. If ``electrons.addIntegerAttributes = upstream``
    and ``electrons.upstream(x,y,z,ux,uy,uz,t) = (x>0.0)*1`` is provided
    then, an integer attribute ``upstream`` is added to all electron particles
    and when these particles are generated, the particles with position less than ``0``
    are assigned a value of ``1``.

* ``<species_name>.addRealAttributes`` (list of `string`)
    User-defined real particle attribute for species, ``species_name``.
    These real attributes will be initialized with user-defined functions
    when the particles are generated.
    If the user-defined real attribute is ``<real_attrib_name>`` then the
    following required parameter must be specified to initialize the attribute.

   * ``<species_name>.attribute.<real_attrib_name>(x,y,z,ux,uy,uz,t)`` (`string`)
     ``t`` represents the physical time in seconds during the simulation.
     ``x``, ``y``, ``z`` represent particle positions in the unit of meter.
     ``ux``, ``uy``, ``uz`` represent the particle velocities in the unit of
     :math:`\gamma v/c`, where
     :math:`\gamma` is the Lorentz factor,
     :math:`v/c` is the particle velocity normalized by the speed of light.

* ``<species>.save_particles_at_xlo/ylo/zlo``,  ``<species>.save_particles_at_xhi/yhi/zhi`` and ``<species>.save_particles_at_eb`` (`0` or `1` optional, default `0`)
    If `1` particles of this species will be copied to the scraped particle
    buffer for the specified boundary if they leave the simulation domain in
    the specified direction. **If USE_EB=TRUE** the ``save_particles_at_eb``
    flag can be set to `1` to also save particle data for the particles of this
    species that impact the embedded boundary.
    The scraped particle buffer can be used to track particle fluxes out of the
    simulation but is currently only accessible via the Python interface. The
    function ``get_particle_boundary_buffer``, found in the
    ``picmi.Simulation`` class as
    ``sim.extension.get_particle_boundary_buffer()``, can be
    used to access the scraped particle buffer. An entry is included for every
    particle in the buffer of the timestep at which the particle was scraped.
    This can be accessed by passing the argument ``comp_name="step_scraped"`` to
    the above mentioned function.

    .. note::
        Currently the scraped particle buffer relies on the user to access the
        data in the buffer for processing and periodically clear the buffer. The
        buffer will grow unbounded as particles are scraped and therefore could
        lead to memory issues if not periodically cleared. To clear the buffer
        call ``warpx_clearParticleBoundaryBuffer()``.

* ``<species>.do_back_transformed_diagnostics`` (`0` or `1` optional, default `1`)
    Only used when ``warpx.do_back_transformed_diagnostics=1``. When running in a
    boosted frame, whether or not to plot back-transformed diagnostics for
    this species.

* ``warpx.serialize_initial_conditions`` (`0` or `1`) optional (default `0`)
    Serialize the initial conditions for reproducible testing.
    Mainly whether or not to use OpenMP threading for particle initialization.

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

* ``<species>.do_qed_quantum_sync`` (`int`) optional (default `0`)
    Enables Quantum synchrotron emission for this species.
    Quantum synchrotron lookup table should be either generated or loaded from disk to enable
    this process (see "Lookup tables for QED modules" section below).
    `<species>` must be either an electron or a positron species.
    **This feature requires to compile with QED=TRUE**

* ``<species>.do_qed_breit_wheeler`` (`int`) optional (default `0`)
    Enables non-linear Breit-Wheeler process for this species.
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

* ``<laser_name>.position`` (`3 floats in 3D and 2D` ; in meters)
    The coordinates of one of the point of the antenna that will emit the laser.
    The plane of the antenna is entirely defined by ``<laser_name>.position``
    and ``<laser_name>.direction``.

    ``<laser_name>.position`` also corresponds to the origin of the coordinates system
    for the laser tranverse profile. For instance, for a Gaussian laser profile,
    the peak of intensity will be at the position given by ``<laser_name>.position``.
    This variable can thus be used to shift the position of the laser pulse
    transversally.

    .. note::
        In 2D, ``<laser_name>.position`` is still given by 3 numbers,
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

* ``<laser_name>.direction`` (`3 floats in 3D`)
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

        E_{max} = a_0 \frac{2 \pi m_e c^2}{e\lambda} = a_0 \times (4.0 \cdot 10^{12} \;V.m^{-1})

    When running a **boosted-frame simulation**, provide the value of ``<laser_name>.e_max``
    in the laboratory frame, and use ``warpx.gamma_boost`` to automatically
    perform the conversion to the boosted frame.

* ``<laser_name>.a0`` (`float` ; dimensionless)
    Peak normalized amplitude of the laser field (given in the lab frame, just as ``e_max`` above).
    See the description of ``<laser_name>.e_max`` for the conversion between ``a0`` and ``e_max``.
    Exactly one of ``a0`` and ``e_max`` must be specified.

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
      setting the additional parameter ``<laser_name>.txye_file_name`` (`string`). It accepts an
      optional parameter ``<laser_name>.time_chunk_size`` (`int`). This allows to read only
      time_chunk_size timesteps from the binary file. New timesteps are read as soon as they are needed.
      The default value is automatically set to the number of timesteps contained in the binary file
      (i.e. only one read is performed at the beginning of the simulation).
      It also accepts the optional parameter ``<laser_name>.delay`` (`float`; in seconds), which allows
      delaying (``delay > 0``) or anticipating (``delay < 0``) the laser by the specified amount of time.
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


* ``<laser_name>.profile_t_peak`` (`float`; in seconds)
    The time at which the laser reaches its peak intensity, at the position
    given by ``<laser_name>.position`` (only used for the ``"gaussian"`` profile)

    When running a **boosted-frame simulation**, provide the value of
    ``<laser_name>.profile_t_peak`` in the laboratory frame, and use ``warpx.gamma_boost``
    to automatically perform the conversion to the boosted frame.

* ``<laser_name>.profile_duration`` (`float` ; in seconds)
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

* ``<laser_name>.phi0`` (`float`; in radians) optional (default `0.`)
    The Carrier Envelope Phase, i.e. the phase of the laser oscillation, at the
    position where the laser envelope is maximum (only used for the ``"gaussian"`` profile)

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

External fields
---------------

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
    input file. For a two-dimensional simulation, it is assumed that the first dimension
    is `x` and the second dimension in `z`, and the value of `y` is set to zero.
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
    be used to initialize the external electric field on the grid. It
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

* ``particles.E_ext_particle_init_style`` & ``particles.B_ext_particle_init_style`` (string) optional (default "none")
    These parameters determine the type of the external electric and
    magnetic fields respectively that are applied directly to the particles at every timestep.
    The field values are specified in the lab frame.
    With the default ``none`` style, no field is applied.
    Possible values are ``constant``, ``parse_E_ext_particle_function`` or ``parse_B_ext_particle_function``, or
    ``repeated_plasma_lens``.

    * ``constant``: a constant field is applied, given by the input parameters
      ``particles.E_external_particle`` or ``particles.B_external_particle``, which are lists of the field components.

    * ``parse_E_ext_particle_function`` or ``parse_B_ext_particle_function``: the field is specified as an analytic
      expression that is a function of space (x,y,z) and time (t), relative to the lab frame.
      The E-field is specified by the input parameters:

        * ``particles.Ex_external_particle_function(x,y,z,t)``

        * ``particles.Ey_external_particle_function(x,y,z,t)``

        * ``particles.Ez_external_particle_function(x,y,z,t)``

      The B-field is specified by the input parameters:

        * ``particles.Bx_external_particle_function(x,y,z,t)``

        * ``particles.By_external_particle_function(x,y,z,t)``

        * ``particles.Bz_external_particle_function(x,y,z,t)``

      Note that the position is defined in Cartesian coordinates, as a function of (x,y,z), even for RZ.

    * ``repeated_plasma_lens``: apply a series of plasma lenses. The properties of the lenses are defined in the
      lab frame by the input parameters:

        * ``repeated_plasma_lens_period``, the period length of the repeat, a single float number,

        * ``repeated_plasma_lens_starts``, the start of each lens relative to the period, an array of floats,

        * ``repeated_plasma_lens_lengths``, the length of each lens, an array of floats,

        * ``repeated_plasma_lens_strengths_E``, the electric focusing strength of each lens, an array of floats, when
          ``particles.E_ext_particle_init_style`` is set to ``repeated_plasma_lens``.

        * ``repeated_plasma_lens_strengths_B``, the magnetic focusing strength of each lens, an array of floats, when
          ``particles.B_ext_particle_init_style`` is set to ``repeated_plasma_lens``.

      The applied field is uniform longitudinally (along z) with a hard edge,
      where residence corrections are used for more accurate field calculation. On the time step when a particle enters
      or leaves each lens, the field applied is scaled by the fraction of the time step spent within the lens.
      The fields are of the form :math:`E_x = \mathrm{strength} \cdot x`, :math:`E_y = \mathrm{strength} \cdot y`,
      and :math:`E_z = 0`, and
      :math:`B_x = \mathrm{strength} \cdot y`, :math:`B_y = -\mathrm{strength} \cdot x`, and :math:`B_z = 0`.

.. _running-cpp-parameters-collision:

Collision initialization
------------------------

WarpX provides several particle collions models, using varying degrees of approximation.

* ``collisions.collision_names`` (`strings`, separated by spaces)
    The name of each collision type.
    This is then used in the rest of the input deck;
    in this documentation we use ``<collision_name>`` as a placeholder.

* ``<collision_name>.type`` (`string`) optional
    The type of collsion. The types implemented are:

    - ``pairwisecoulomb`` for pairwise Coulomb collisions, the default if unspecified.
      This provides a pair-wise relativistic elastic Monte Carlo binary Coulomb collision model,
      following the algorithm given by `Perez et al. (Phys. Plasmas 19, 083104, 2012) <https://doi.org/10.1063/1.4742167>`_.
      When the RZ mode is used, `warpx.n_rz_azimuthal_modes` must be set to 1 at the moment,
      since the current implementation of the collision module assumes axisymmetry.
    - ``background_mcc`` for collisions between particles and a neutral background.
      This is a relativistic Monte Carlo treatment for particles colliding
      with a neutral background gas. The implementation follows the so-called
      null collision strategy discussed for example in `Birdsall (IEEE Transactions on
      Plasma Science, vol. 19, no. 2, pp. 65-85, 1991) <https://ieeexplore.ieee.org/document/106800>`_.
      See also :ref:`collisions section <theory-collisions>`.
    - ``background_stopping`` for slowing of ions due to collisions with electrons or ions.
      This implements the approximate formulae as derived in Introduction to Plasma Physics,
      from Goldston and Rutherford, section 14.2.

* ``<collision_name>.species`` (`strings`)
    If using ``pairwisecoulomb`` type this should be the names of two species,
    between which the collision will be considered.
    The number of provided ``<collision_name>.species`` should match
    the number of collision names, i.e. ``collisions.collision_names``.
    If using ``background_mcc`` or ``background_stopping`` type this should be the name of the species for which collisions will be included.
    Only one species name should be given.

* ``<collision_name>.CoulombLog`` (`float`) optional
    Only for ``pairwisecoulomb``. A provided fixed Coulomb logarithm of the
    collision type ``<collision_name>``.
    For example, a typical Coulomb logarithm has a form of
    :math:`\ln(\lambda_D/R)`,
    where :math:`\lambda_D` is the Debye length,
    :math:`R\approx1.4A^{1/3}` is the effective Coulombic radius of the nucleus,
    :math:`A` is the mass number.
    If this is not provided, or if a non-positive value is provided,
    a Coulomb logarithm will be computed automatically according to the algorithm in
    `Perez et al. (Phys. Plasmas 19, 083104, 2012) <https://doi.org/10.1063/1.4742167>`_.

* ``<collision_name>.ndt`` (`int`) optional
    Execute collision every # time steps. The default value is 1.

* ``<collision_name>.background_density`` (`float`)
    Only for ``background_mcc`` and ``background_stopping``. The density of the background in :math:`m^{-3}`.
    Can also provide ``<collision_name>.background_density(x,y,z,t)`` using the parser
    initialization style for spatially and temporally varying density. With ``background_mcc``, if a function
    is used for the background density, the input parameter ``<collision_name>.max_background_density``
    must also be provided to calculate the maximum collision probability.

* ``<collision_name>.background_temperature`` (`float`)
    Only for ``background_mcc`` and ``background_stopping``. The temperature of the background in Kelvin.
    Can also provide ``<collision_name>.background_temperature(x,y,z,t)`` using the parser
    initialization style for spatially and temporally varying temperature.

* ``<collision_name>.background_mass`` (`float`) optional
    Only for ``background_mcc`` and ``background_stopping``. The mass of the background gas in kg.
    With ``background_mcc``, if not given the mass of the colliding species will be used unless ionization is
    included in which case the mass of the product species will be used.
    With ``background_stopping``, and ``background_type`` set to ``electrons``, if not given defaults to the electron mass. With
    ``background_type`` set to ``ions``, the mass must be given.

* ``<collision_name>.background_charge_state`` (`float`)
    Only for ``background_stopping``, where it is required when ``background_type`` is set to ``ions``.
    This specifies the charge state of the background ions.

* ``<collision_name>.background_type`` (`string`)
    Only for ``background_stopping``, where it is required, the type of the background.
    The possible values are ``electrons`` and ``ions``. When ``electrons``, equation 14.12 from Goldston and Rutherford is used.
    This formula is based on Coulomb collisions with the approximations that :math:`M_b >> m_e` and :math:`V << v_{thermal\_e}`,
    and the assumption that the electrons have a Maxwellian distribution with temperature :math:`T_e`.

    .. math::
        \frac{dV}{dt} = - \frac{2^{1/2}n_eZ_b^2e^4m_e^{1/2}\log\Lambda}{12\pi^{3/2}\epsilon_0M_bT_e^{3/2}}V

    where :math:`V` is each velocity component, :math:`n_e` is the background density, :math:`Z_b` is the ion charge state,
    :math:`e` is the electron charge, :math:`m_e` is the background mass, :math:`\log\Lambda=\log((12\pi/Z_b)(n_e\lambda_{de}^3))`,
    :math:`\lambda_{de}` is the DeBye length, and :math:`M_b` is the ion mass.
    The equation is integrated over a time step, giving :math:`V(t+dt) = V(t)*\exp(-\alpha*{dt})`
    where :math:`\alpha` is the factor multiplying :math:`V`.

    When ``ions``, equation 14.20 is used.
    This formula is based on Coulomb collisions with the approximations that :math:`M_b >> M` and :math:`V >> v_{thermal\_i}`.
    The background ion temperature only appears in the :math:`\log\Lambda` term.

    .. math::
        \frac{dW_b}{dt} = - \frac{2^{1/2}n_iZ^2Z_b^2e^4M_b^{1/2}\log\Lambda}{8\pi\epsilon_0MW_b^{1/2}}

    where :math:`W_b` is the ion energy, :math:`n_i` is the background density,
    :math:`Z` is the charge state of the background ions, :math:`Z_b` is the ion charge state,
    :math:`e` is the electron charge, :math:`M_b` is the ion mass, :math:`\log\Lambda=\log((12\pi/Z_b)(n_i\lambda_{di}^3))`,
    :math:`\lambda_{di}` is the DeBye length, and :math:`M` is the background ion mass.
    The equation is integrated over a time step, giving :math:`W_b(t+dt) = ((W_b(t)^{3/2}) - 3/2\beta{dt})^{2/3}`
    where :math:`\beta` is the term on the r.h.s except :math:`W_b`.


* ``<collision_name>.scattering_processes`` (`strings` separated by spaces)
    Only for ``background_mcc``. The MCC scattering processes that should be
    included. Available options are ``elastic``, ``back`` & ``charge_exchange``
    for ions and ``elastic``, ``excitationX`` & ``ionization`` for electrons.
    The ``elastic`` option uses hard-sphere scattering, with a differential
    cross section that is independent of angle.
    With ``charge_exchange``, the ion velocity is replaced with the neutral
    velocity, chosen from a Maxwellian based on the value of
    ``<collision_name>.background_temperature``.
    Multiple excitation events can be included for electrons corresponding to
    excitation to different levels, the ``X`` above can be changed to a unique
    identifier for each excitation process. For each scattering process specified
    a path to a cross-section data file must  also be given. We use
    ``<scattering_process>`` as a placeholder going forward.

* ``<collision_name>.<scattering_process>_cross_section`` (`string`)
    Only for ``background_mcc``. Path to the file containing cross-section data
    for the given scattering processes. The cross-section file must have exactly
    2 columns of data, the first containing equally spaced energies in eV and the
    second the corresponding cross-section in :math:`m^2`.

* ``<collision_name>.<scattering_process>_energy`` (`float`)
    Only for ``background_mcc``. If the scattering process is either
    ``excitationX`` or ``ionization`` the energy cost of that process must be given in eV.

* ``<collision_name>.ionization_species`` (`float`)
    Only for ``background_mcc``. If the scattering process is ``ionization`` the
    produced species must also be given. For example if argon properties is used
    for the background gas, a species of argon ions should be specified here.

.. _running-cpp-parameters-numerics:

Numerics and algorithms
-----------------------

* ``warpx.cfl`` (`float`) optional (default `0.999`)
    The ratio between the actual timestep that is used in the simulation
    and the Courant-Friedrichs-Lewy (CFL) limit. (e.g. for `warpx.cfl=1`,
    the timestep will be exactly equal to the CFL limit.)

* ``warpx.use_filter`` (`0` or `1`; default: `1`, except for RZ FDTD)
    Whether to smooth the charge and currents on the mesh, after depositing them from the macro-particles.
    This uses a bilinear filter (see the :ref:`filtering section <theory-pic-filter>`).
    The default is `1` in all cases, except for simulations in RZ geometry using the FDTD solver.
    With the RZ PSATD solver, the filtering is done in :math:`k`-space.

    .. warning::

       Known bug: filter currently not working with FDTD solver in RZ geometry (see https://github.com/ECP-WarpX/WarpX/issues/1943).

* ``warpx.filter_npass_each_dir`` (`3 int`) optional (default `1 1 1`)
    Number of passes along each direction for the bilinear filter.
    In 2D simulations, only the first two values are read.

* ``warpx.use_filter_compensation`` (`0` or `1`; default: `0`)
    Whether to add compensation when applying filtering.
    This is only supported with the RZ spectral solver.

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
       the :ref:`particle-in-cell theory section <theory-pic>`.

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

* ``algo.particle_shape`` (`integer`; `1`, `2`, or `3`)
    The order of the shape factors (splines) for the macro-particles along all spatial directions: `1` for linear, `2` for quadratic, `3` for cubic.
    Low-order shape factors result in faster simulations, but may lead to more noisy results.
    High-order shape factors are computationally more expensive, but may increase the overall accuracy of the results. For production runs it is generally safer to use high-order shape factors, such as cubic order.

    Note that this input parameter is not optional and must always be set in all input files provided that there is at least one particle species (set in input as ``particles.species_names``) or one laser species (set in input as ``lasers.names``) in the simulation. No default value is provided automatically.

* ``algo.maxwell_solver`` (`string`, optional)
    The algorithm for the Maxwell field solver.
    Available options are:

     - ``yee``: Yee FDTD solver.
     - ``ckc``: (not available in ``RZ`` geometry) Cole-Karkkainen solver with Cowan
       coefficients (see `Cowan, PRSTAB 16 (2013) <https://journals.aps.org/prab/abstract/10.1103/PhysRevSTAB.16.041303>`__)
     - ``psatd``: Pseudo-spectral solver (see :ref:`theory <theory-pic-mwsolve-psatd>`)
     - ``ect``: Enlarged cell technique (conformal finite difference solver. See Xiao and Liu,
                IEEE Antennas and Propagation Society International Symposium (2005),
                <https://ieeexplore.ieee.org/document/1551259>)

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
     To initialize spatially varying conductivity, permittivity, and permeability, respectively,
     using a mathematical function in the input. Constants required in the
     mathematical expression can be set using ``my_constants``. These parameters are parsed
     if ``algo.em_solver_medium=macroscopic``.

* ``macroscopic.sigma``, ``macroscopic.epsilon``, ``macroscopic.mu`` (`double`)
    To initialize a constant conductivity, permittivity, and permeability of the
    computational medium, respectively. The default values are the corresponding values
    in vacuum.

* ``interpolation.galerkin_scheme`` (`0` or `1`)
    Whether to use a Galerkin scheme when gathering fields to particles.
    When set to `1`, the interpolation orders used for field-gathering are reduced for certain field components along certain directions.
    For example, :math:`E_z` is gathered using ``algo.particle_shape`` along :math:`(x,y)` and ``algo.particle_shape - 1`` along :math:`z`.
    See equations (21)-(23) of (`Godfrey and Vay, 2013 <https://doi.org/10.1016/j.jcp.2013.04.006>`_) and associated references for details.
    Defaults to `1` unless ``warpx.do_nodal = 1`` and/or ``algo.field_gathering = momentum-conserving``.

    .. warning::

        The default behavior should not normally be changed.
        At present, this parameter is intended mainly for testing and development purposes.

* ``interpolation.field_centering_nox``, ``interpolation.field_centering_noy``, ``interpolation.field_centering_noz`` (default: ``2`` in all directions)
    The order of interpolation used with staggered grids (``warpx.do_nodal = 0``) and momentum-conserving field gathering (``algo.field_gathering = momentum-conserving``) to interpolate the electric and magnetic fields from the cell centers to the cell nodes, before gathering the fields from the cell nodes to the particle positions. High-order interpolation (order 8 in each direction, at least) is necessary to ensure stability in typical LWFA boosted-frame simulations using the Galilean PSATD or comoving PSATD schemes. This finite-order interpolation is used only when the PSATD solver is used for Maxwell's equations. With the FDTD solver, basic linear interpolation is used instead.

* ``interpolation.current_centering_nox``, ``interpolation.current_centering_noy``, ``interpolation.current_centering_noz`` (default: ``2`` in all directions)
    The order of interpolation used to center the currents from nodal to staggered grids (if ``warpx.do_current_centering = 1``), before pushing the Maxwell fields on staggered grids. This finite-order interpolation is used only when the PSATD solver is used for Maxwell's equations. With the FDTD solver, basic linear interpolation is used instead.

* ``warpx.do_current_centering`` (`0` or `1` ; default: 0)
    If true, the current is deposited on a nodal grid and then centered to a staggered grid (Yee grid), using finite-order interpolation. If ``warpx.do_nodal = 1``, the Maxwell fields are pushed on a nodal grid, it is not necessary to center the currents to a staggered grid, and we set therefore ``warpx.do_current_centering = 0`` automatically, overwriting the user-defined input.

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

* ``warpx.do_multi_J`` (`0` or `1`; default: `0`)
    Whether to use the multi-J algorithm, where current deposition and field update are performed multiple times within each time step. The number of sub-steps is determined by the input parameter ``warpx.do_multi_J_n_depositions``. Unlike sub-cycling, field gathering is performed only once per time step, as in regular PIC cycles. When ``warpx.do_multi_J = 1``, we perform linear interpolation of two distinct currents deposited at the beginning and the end of the time step, instead of using one single current deposited at half time. For simulations with strong numerical Cherenkov instability (NCI), it is recommended to use the multi-J algorithm in combination with ``psatd.do_time_averaging = 1``.

* ``warpx.do_multi_J_n_depositions`` (integer)
    Number of sub-steps to use with the multi-J algorithm, when ``warpx.do_multi_J = 1``.
    Note that this input parameter is not optional and must always be set in all input files where ``warpx.do_multi_J = 1``. No default value is provided automatically.


* ``psatd.nox``, ``psatd.noy``, ``pstad.noz`` (`integer`) optional (default `16` for all)
    The order of accuracy of the spatial derivatives, when using the code compiled with a PSATD solver.
    If ``psatd.periodic_single_box_fft`` is used, these can be set to ``inf`` for infinite-order PSATD.

* ``psatd.nx_guard``, ``psatd.ny_guard``, ``psatd.nz_guard`` (`integer`) optional
    The number of guard cells to use with PSATD solver.
    If not set by users, these values are calculated automatically and determined *empirically* and
    would be equal the order of the solver for nodal grid, and half the order of the solver for staggered.

* ``psatd.periodic_single_box_fft`` (`0` or `1`; default: 0)
    If true, this will *not* incorporate the guard cells into the box over which FFTs are performed.
    This is only valid when WarpX is run with periodic boundaries and a single box.
    In this case, using `psatd.periodic_single_box_fft` is equivalent to using a global FFT over the whole domain.
    Therefore, all the approximations that are usually made when using local FFTs with guard cells
    (for problems with multiple boxes) become exact in the case of the periodic, single-box FFT without guard cells.

* ``psatd.current_correction`` (`0` or `1`; default: `1`, with the exceptions mentioned below)
    If true, a current correction scheme in Fourier space is applied in order to guarantee charge conservation.
    The default value is ``psatd.current_correction=1``, unless a charge-conserving current deposition scheme is used (by setting ``algo.current_deposition=esirkepov`` or ``algo.current_deposition=vay``) or unless the ``div(E)`` cleaning scheme is used (by setting ``warpx.do_dive_cleaning=1``).

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

    This option is currently implemented only for the standard PSATD, Galilean PSATD, and averaged Galilean PSATD schemes, while it is not yet available for the multi-J algorithm.

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
    in a Galilean frame in the :ref:`theory section <theory-boostedframe-galilean>`).
    It also requires the use of the `direct` current deposition option
    `algo.current_deposition = direct` (does not work with Esirkepov algorithm).

* ``psatd.v_comoving`` (3 floating-point values, in units of the speed of light; default ``0. 0. 0.``)
    Defines the comoving velocity in the comoving PSATD scheme.
    A non-zero comoving velocity selects the comoving PSATD algorithm, which suppresses the numerical Cherenkov instability (NCI) in boosted-frame simulations, under certain assumptions. This option requires that WarpX is compiled with ``USE_PSATD = TRUE``. It also requires the use of direct current deposition (``algo.current_deposition = direct``) and has not been neither implemented nor tested with other current deposition schemes.

* ``psatd.do_time_averaging`` (`0` or `1`; default: 0)
    Whether to use an averaged Galilean PSATD algorithm or standard Galilean PSATD.

* ``warpx.override_sync_intervals`` (`string`) optional (default `1`)
    Using the `Intervals parser`_ syntax, this string defines the timesteps at which
    synchronization of sources (`rho` and `J`) and fields (`E` and `B`) on grid nodes at box
    boundaries is performed. Since the grid nodes at the interface between two neighbor boxes are
    duplicated in both boxes, an instability can occur if they have too different values.
    This option makes sure that they are synchronized periodically.
    Note that if Perfectly Matched Layers (PML) are used, synchronization of the `E` and `B` fields
    is performed at every timestep regardless of this parameter.

* ``warpx.use_hybrid_QED`` (`bool`; default: 0)
    Will use the Hybird QED Maxwell solver when pushing fields: a QED correction is added to the
    field solver to solve non-linear Maxwell's equations, according to [Quantum Electrodynamics
    vacuum polarization solver, P. Carneiro et al., `ArXiv 2016 <https://arxiv.org/abs/1607.04224>`__].
    Note that this option can only be used with the PSATD build. Furthermore,
    warpx.do_nodal must be set to `1` which is not its default value.

* ``warpx.quantum_xi`` (`float`; default: 1.3050122.e-52)
     Overwrites the actual quantum parameter used in Maxwell's QED equations. Assigning a
     value here will make the simulation unphysical, but will allow QED effects to become more apparent.
     Note that this option will only have an effect if the ``warpx.use_Hybrid_QED`` flag is also triggered.

* ``warpx.do_device_synchronize`` (`bool`) optional (default `1`)
    When running in an accelerated platform, whether to call a ``amrex::Gpu::synchronize()`` around profiling regions.
    This allows the profiler to give meaningful timers, but (hardly) slows down the simulation.

* ``warpx.sort_intervals`` (`string`) optional (defaults: ``-1`` on CPU; ``4`` on GPU)
     Using the `Intervals parser`_ syntax, this string defines the timesteps at which particles are
     sorted by bin.
     If ``<=0``, do not sort particles.
     It is turned on on GPUs for performance reasons (to improve memory locality).

* ``warpx.sort_bin_size`` (list of `int`) optional (default ``1 1 1``)
     If ``sort_intervals`` is activated particles are sorted in bins of ``sort_bin_size`` cells.
     In 2D, only the first two elements are read.

.. _running-cpp-parameters-diagnostics:

Diagnostics and output
----------------------

.. _running-cpp-parameters-diagnostics-insitu:

In-situ visualization
^^^^^^^^^^^^^^^^^^^^^

WarpX has four types of diagnostics:
``FullDiagnostics`` consist in dumps of fields and particles at given iterations,
``BackTransformedDiagnostics`` are used when running a simulation in a boosted frame, to reconstruct output data to the lab frame,
``BoundaryScrapingDiagnostics`` are used to collect the particles that are absorbed at the boundary, throughout the simulation, and
``ReducedDiags`` allow the user to compute some reduced quantity (particle temperature, max of a field) and write a small amount of data to text files.
Similar to what is done for physical species, WarpX has a class Diagnostics that allows users to initialize different diagnostics, each of them with different fields, resolution and period.
This currently applies to standard diagnostics, but should be extended to back-transformed diagnostics and reduced diagnostics (and others) in a near future.

.. _running-cpp-parameters-diagnostics-full:

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

* ``<diag_name>.intervals`` (`string`)
    Using the `Intervals parser`_ syntax, this string defines the timesteps at which data is dumped.
    Use a negative number or 0 to disable data dumping.
    example: ``diag1.intervals = 10,20:25:1``.
    Note that by default the last timestep is dumped regardless of this parameter. This can be
    changed using the parameter ``<diag_name>.dump_last_timestep`` described below.

* ``<diag_name>.dump_last_timestep`` (`bool` optional, default `1`)
    If this is `1`, the last timestep is dumped regardless of ``<diag_name>.period``.

* ``<diag_name>.diag_type`` (`string`)
    Type of diagnostics. ``Full`` and ``BackTransformed``
    example: ``diag1.diag_type = Full`` or ``diag1.diag_type = BackTransformed``

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

* ``<diag_name>.openpmd_encoding`` (optional, ``v`` (variable based), ``f`` (file based) or ``g`` (group based) ) only read if ``<diag_name>.format = openpmd``.
     openPMD `file output encoding <https://openpmd-api.readthedocs.io/en/0.14.0/usage/concepts.html#iteration-and-series>`__.
     File based: one file per timestep (slower), group/variable based: one file for all steps (faster)).
     ``variable based`` is an `experimental feature with ADIOS2 <https://openpmd-api.readthedocs.io/en/0.14.0/backends/adios2.html#experimental-new-adios2-schema>`__ and not supported for back-transformed diagnostics.
     Default: ``f`` (full diagnostics)

* ``<diag_name>.adios2_operator.type`` (``zfp``, ``blosc``) optional,
    `ADIOS2 I/O operator type <https://openpmd-api.readthedocs.io/en/0.14.0/details/backendconfig.html#adios2>`__ for `openPMD <https://www.openPMD.org>`_ data dumps.

* ``<diag_name>.adios2_operator.parameters.*`` optional,
    `ADIOS2 I/O operator parameters <https://openpmd-api.readthedocs.io/en/0.14.0/details/backendconfig.html#adios2>`__ for `openPMD <https://www.openPMD.org>`_ data dumps.

    A typical example for `ADIOS2 output using lossless compression <https://openpmd-api.readthedocs.io/en/0.14.0/details/backendconfig.html#adios2>`__ with ``blosc`` using the ``zstd`` compressor and 6 CPU treads per MPI Rank (e.g. for a `GPU run with spare CPU resources <https://arxiv.org/abs/1706.00522>`__):

    .. code-block:: text

        <diag_name>.adios2_operator.type = blosc
        <diag_name>.adios2_operator.parameters.compressor = zstd
        <diag_name>.adios2_operator.parameters.clevel = 1
        <diag_name>.adios2_operator.parameters.doshuffle = BLOSC_BITSHUFFLE
        <diag_name>.adios2_operator.parameters.threshold = 2048
        <diag_name>.adios2_operator.parameters.nthreads = 6  # per MPI rank (and thus per GPU)

    or for the lossy ZFP compressor using very strong compression per scalar:

    .. code-block:: text

        <diag_name>.adios2_operator.type = zfp
        <diag_name>.adios2_operator.parameters.precision = 3

* ``<diag_name>.adios2_engine.type`` (``bp4``, ``sst``, ``ssc``, ``dataman``) optional,
    `ADIOS2 Engine type <https://openpmd-api.readthedocs.io/en/0.14.0/details/backendconfig.html#adios2>`__ for `openPMD <https://www.openPMD.org>`_ data dumps.
    See full list of engines at `ADIOS2 readthedocs <https://adios2.readthedocs.io/en/latest/engines/engines.html>`__

* ``<diag_name>.adios2_engine.parameters.*`` optional,
    `ADIOS2 Engine parameters <https://openpmd-api.readthedocs.io/en/0.14.0/details/backendconfig.html#adios2>`__ for `openPMD <https://www.openPMD.org>`_ data dumps.

    An example for parameters for the BP engine are setting the number of writers (``NumAggregators``), transparently redirecting data to burst buffers etc.
    A detailed list of engine-specific parameters are available at the official `ADIOS2 documentation <https://adios2.readthedocs.io/en/latest/engines/engines.html>`__

    .. code-block:: text

        <diag_name>.adios2_engine.parameters.NumAggregators = 2048
        <diag_name>.adios2_engine.parameters.BurstBufferPath="/mnt/bb/username"

* ``<diag_name>.fields_to_plot`` (list of `strings`, optional)
    Fields written to output.
    Possible scalar fields: ``part_per_cell`` ``rho`` ``phi`` ``F`` ``part_per_grid`` ``divE`` ``divB`` and ``rho_<species_name>``, where ``<species_name>`` must match the name of one of the available particle species. Note that ``phi`` will only be written out when do_electrostatic==labframe. Also, note that for ``<diag_name>.diag_type = BackTransformed``, the only scalar field currently supported is ``rho``.
    Possible vector field components in Cartesian geometry: ``Ex`` ``Ey`` ``Ez`` ``Bx`` ``By`` ``Bz`` ``jx`` ``jy`` ``jz``.
    Possible vector field components in RZ geometry: ``Er`` ``Et`` ``Ez`` ``Br`` ``Bt`` ``Bz`` ``jr`` ``jt`` ``jz``.
    Default is ``<diag_name>.fields_to_plot = Ex Ey Ez Bx By Bz jx jy jz``,
    unless in RZ geometry with ``<diag_name>.format == openpmd``,
    then default is ``<diag_name>.fields_to_plot = Er Et Ez Br Bt Bz jr jt jz``.
    Note that the fields are averaged on the cell centers before they are written to file.
    Also, when ``<diag_name>.format = openpmd``, the RZ modes for all fields are written.
    Otherwise, we reconstruct a 2D Cartesian slice of the fields for output at :math:`\theta=0`.

* ``<diag_name>.dump_rz_modes`` (`0` or `1`) optional (default `0`)
    Whether to save all modes when in RZ.  When ``openpmd_backend = openpmd``, this parameter is ignored and all modes are saved.

* ``<diag_name>.particle_fields_to_plot`` (list of `strings`, optional)
   Names of per-cell diagnostics of particle properties to calculate and output as additional fields.
   Note that the deposition onto the grid does not respect the particle shape factor, but instead uses nearest-grid point interpolation.
   Default is none.
   Parser functions for these field names are specified by ``<diag_name>.particle_fields.<field_name>(x,y,z,ux,uy,uz)``.
   Also, note that this option is only available for ``<diag_name>.diag_type = Full``

* ``<diag_name>.particle_fields_species`` (list of `strings`, optional)
         Species for which to calculate ``particle_fields_to_plot``.
         Fields will be calculated separately for each specified species.
         The default is a list of all of the available particle species.

* ``<diag_name>.particle_fields.<field_name>.do_average`` (`0` or `1`) optional (default `1`)
   Whether the diagnostic is an average or a sum. With an average, the sum over the specified function is divided
   by the sum of the particle weights in each cell.

* ``<diag_name>.particle_fields.<field_name>(x,y,z,ux,uy,uz)`` (parser `string`)
   Parser function to be calculated for each particle per cell. The averaged field written is

        .. math::

            \texttt{<field_name>_<species>} = \frac{\sum_{i=1}^N w_i \, f(x_i,y_i,z_i,u_{x,i},u_{y,i},u_{z,i})}{\sum_{i=1}^N w_i}

   where :math:`w_i` is the particle weight, :math:`f()` is the parser function, and :math:`(x_i,y_i,z_i)` are particle positions in units of a meter. The sums are over all particles of type ``<species>`` in a cell (ignoring the particle shape factor) that satisfy ``<diag_name>.particle_fields.<field_name>.filter(x,y,z,ux,uy,uz)``.
   When ``<diag_name>.particle_fields.<field_name>.do_average`` is `0`, the division by the sum over particle weights is not done.
   In 1D or 2D, the particle coordinates will follow the WarpX convention. :math:`(u_{x,i},u_{y,i},u_{z,i})` are components of the particle four-velocity. :math:`u = \gamma v/c`, :math:`\gamma` is the Lorentz factor, :math:`v` is the particle velocity, and :math:`c` is the speed of light.

* ``<diag_name>.particle_fields.<field_name>.filter(x,y,z,ux,uy,uz)`` (parser `string`, optional)
    Parser function returning a boolean for whether to include a particle in the diagnostic.
    If not specified, all particles will be included (see above).
    The function arguments are the same as above.

* ``<diag_name>.plot_raw_fields`` (`0` or `1`) optional (default `0`)
    By default, the fields written in the plot files are averaged on the cell centers.
    When ``<diag_name>.plot_raw_fields = 1``, then the raw (i.e. non-averaged)
    fields are also saved in the output files.
    Only works with ``<diag_name>.format = plotfile``.
    See `this section <https://yt-project.org/doc/examining/loading_data.html#viewing-raw-fields-in-warpx>`_
    in the yt documentation for more details on how to view raw fields.

* ``<diag_name>.plot_raw_fields_guards`` (`0` or `1`) optional (default `0`)
    Only used when ``<diag_name>.plot_raw_fields = 1``.
    Whether to include the guard cells in the output of the raw fields.
    Only works with ``<diag_name>.format = plotfile``.

* ``<diag_name>.coarsening_ratio`` (list of `int`) optional (default `1 1 1`)
    Reduce size of the field output by this ratio in each dimension.
    (This is done by averaging the field over 1 or 2 points along each direction, depending on the staggering).
    If ``blocking_factor`` and ``max_grid_size`` are used for the domain decomposition, as detailed in
    the :ref:`parallelization <parallelization_warpx>` section, ``coarsening_ratio`` should be an integer
    divisor of ``blocking_factor``. If ``warpx.numprocs`` is used instead, the total number of cells in a given
    dimension must be a multiple of the ``coarsening_ratio`` multiplied by ``numprocs`` in that dimension.

* ``<diag_name>.file_prefix`` (`string`) optional (default `diags/<diag_name>`)
    Root for output file names. Supports sub-directories.

* ``<diag_name>.file_min_digits`` (`int`) optional (default `6`)
    The minimum number of digits used for the iteration number appended to the diagnostic file names.

* ``<diag_name>.diag_lo`` (list `float`, 1 per dimension) optional (default `-infinity -infinity -infinity`)
    Lower corner of the output fields (if smaller than ``warpx.dom_lo``, then set to ``warpx.dom_lo``). Currently, when the ``diag_lo`` is different from ``warpx.dom_lo``, particle output is disabled.

* ``<diag_name>.diag_hi`` (list `float`, 1 per dimension) optional (default `+infinity +infinity +infinity`)
    Higher corner of the output fields (if larger than ``warpx.dom_hi``, then set to ``warpx.dom_hi``). Currently, when the ``diag_hi`` is different from ``warpx.dom_hi``, particle output is disabled.

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
    Users can provide an expression returning a boolean for whether a particle is dumped.
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
    when using the AMReX plotfile format.
    Please see the :ref:`data analysis section <dataanalysis-formats>` for more information.

* ``amrex.async_out_nfiles`` (`int`) optional (default `64`)
    The maximum number of files to write to when using asynchronous IO.
    To use asynchronous IO with more than ``amrex.async_out_nfiles`` MPI ranks,
    WarpX must be configured with ``-DWarpX_MPI_THREAD_MULTIPLE=ON``.
    Please see the :ref:`data analysis section <dataanalysis-formats>` for more information.

* ``warpx.field_io_nfiles`` and ``warpx.particle_io_nfiles`` (`int`) optional (default `1024`)
    The maximum number of files to use when writing field and particle data to plotfile directories.

* ``warpx.mffile_nstreams`` (`int`) optional (default `4`)
    Limit the number of concurrent readers per file.

.. _running-cpp-parameters-diagnostics-btd:

BackTransformed Diagnostics (with support for Plotfile/openPMD output)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``BackTransformed`` diag type are used when running a simulation in a boosted frame, to reconstruct output data to the lab frame. This option can be set using ``<diag_name>.diag_type = BackTransformed``. Note that this diagnostic is not currently supported for RZ.  Additional options for this diagnostic include:

* ``<diag_name>.num_snapshots_lab`` (`integer`)
    Only used when ``<diag_name>.diag_type`` is ``BackTransformed``.
    The number of lab-frame snapshots that will be written.

* ``<diag_name>.dt_snapshots_lab`` (`float`, in seconds)
    Only used when ``<diag_name>.diag_type`` is ``BackTransformed``.
    The time interval inbetween the lab-frame snapshots (where this
    time interval is expressed in the laboratory frame).

* ``<diag_name>.dz_snapshots_lab`` (`float`, in meters)
    Only used when ``<diag_name>.diag_type`` is ``BackTransformed``.
    Distance between the lab-frame snapshots (expressed in the laboratory
    frame). ``dt_snapshots_lab`` is then computed by
    ``dt_snapshots_lab = dz_snapshots_lab/c``. Either `dt_snapshots_lab`
    or `dz_snapshot_lab` is required.

* ``<diag_name>.buffer_size`` (`integer`)
    Only used when ``<diag_name>.diag_type`` is ``BackTransformed``.
    The default size of the back transformed diagnostic buffers used to generate lab-frame
    data is 256. That is, when the multifab with lab-frame data has 256 z-slices,
    the data will be flushed out. However, if many lab-frame snapshots are required for
    diagnostics and visualization, the GPU may run out of memory with many large boxes with
    a size of 256 in the z-direction. This input parameter can then be used to set a
    smaller buffer-size, preferably multiples of 8, such that, a large number of
    lab-frame snapshot data can be generated without running out of gpu memory.
    The downside to using a small buffer size, is that the I/O time may increase due
    to frequent flushes of the lab-frame data. The other option is to keep the default
    value for buffer size and use slices to reduce the memory footprint and maintain
    optimum I/O performance.

Back-Transformed Diagnostics (legacy output)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``BackTransformedDiagnostics`` are used when running a simulation in a boosted frame, to reconstruct output data to the lab frame

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

* ``warpx.buffer_size`` (`integer`)
    The default size of the back transformed diagnostic buffers used to generate lab-frame
    data is 256. That is, when the multifab with lab-frame data has 256 z-slices,
    the data will be flushed out. However, if many lab-frame snapshots are required for
    diagnostics and visualization, the GPU may run out of memory with many large boxes with
    a size of 256 in the z-direction. This input parameter can then be used to set a
    smaller buffer-size, preferably multiples of 8, such that, a large number of
    lab-frame snapshot data can be generated without running out of gpu memory.
    The downside to using a small buffer size, is that the I/O time may increase due
    to frequent flushes of the lab-frame data. The other option is to keep the default
    value for buffer size and use slices to reduce the memory footprint and maintain
    optimum I/O performance.

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

Boundary Scraping Diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``BoundaryScrapingDiagnostics`` are used to collect the particles that are absorbed at the embedded boundary, throughout the simulation.
(Note that this diagnostics does not save any field ; it only saves particles.)
Currently, the only supported output format is openPMD, so the user needs to set ``<diag>.format=openpmd``. In addition, the user needs
to set ``<species>.save_particles_at_eb=1`` for each of the species that are to be saved in this diagnostic.

In addition to their usual attributes, the saved particles have an additional integer attribute ``timestamp``, which
indicates the PIC iteration at which each particle was absorbed at the boundary.

.. _running-cpp-parameters-diagnostics-reduced:

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
    in this documentation we use ``<reduced_diags_name>`` as a placeholder.

* ``<reduced_diags_name>.type`` (`string`)
    The type of reduced diagnostics associated with this ``<reduced_diags_name>``.
    For example, ``ParticleEnergy``, ``FieldEnergy``, etc.
    All available types are described below in detail.
    For all reduced diagnostics,
    the first and the second columns in the output file are
    the time step and the corresponding physical time in seconds, respectively.

    * ``ParticleEnergy``
        This type computes the total and mean relativistic particle kinetic energy among all species:

        .. math::

            E_p = \sum_{i=1}^N w_i \, \left( \sqrt{|\boldsymbol{p}_i|^2 c^2 + m_0^2 c^4} - m_0 c^2 \right)

        where :math:`\boldsymbol{p}_i` is the relativistic momentum of the :math:`i`-th particle, :math:`c` is the speed of light, :math:`m_0` is the rest mass, :math:`N` is the number of particles, and :math:`w_i` is the weight of the :math:`i`-th particle.

        The output columns are the total energy of all species, the total energy per species, the total mean energy :math:`E_p / \sum_i w_i` of all species, and the total mean energy per species.

    * ``ParticleMomentum``
        This type computes the total and mean relativistic particle momentum among all species:

        .. math::

            \boldsymbol{P}_p = \sum_{i=1}^N w_i \, \boldsymbol{p}_i

        where :math:`\boldsymbol{p}_i` is the relativistic momentum of the :math:`i`-th particle, :math:`N` is the number of particles, and :math:`w_i` is the weight of the :math:`i`-th particle.

        The output columns are the components of the total momentum of all species, the total momentum per species, the total mean momentum :math:`\boldsymbol{P}_p / \sum_i w_i` of all species, and the total mean momentum per species.

    * ``FieldEnergy``
        This type computes the electromagnetic field energy

        .. math::

            E_f = \frac{1}{2} \sum_{\text{cells}} \left( \varepsilon_0 |\boldsymbol{E}|^2 + \frac{|\boldsymbol{B}|^2}{\mu_0} \right) \Delta V

        where :math:`\boldsymbol{E}` is the electric field, :math:`\boldsymbol{B}` is the magnetic field, :math:`\varepsilon_0` is the vacuum permittivity, :math:`\mu_0` is the vacuum permeability, :math:`\Delta V` is the cell volume (or cell area in 2D), and the sum is over all cells.

        The output columns are the total field energy :math:`E_f`, the :math:`\boldsymbol{E}` field energy, and the :math:`\boldsymbol{B}` field energy, at each mesh refinement level.

    * ``FieldMomentum``
        This type computes the electromagnetic field momentum

        .. math::

            \boldsymbol{P}_f = \varepsilon_0 \sum_{\text{cells}} \left( \boldsymbol{E} \times \boldsymbol{B} \right) \Delta V

        where :math:`\boldsymbol{E}` is the electric field, :math:`\boldsymbol{B}` is the magnetic field, :math:`\varepsilon_0` is the vacuum permittivity, :math:`\Delta V` is the cell volume (or cell area in 2D), and the sum is over all cells.

        The output columns are the components of the total field momentum :math:`\boldsymbol{P}_f` at each mesh refinement level.

        Note that the fields are *not* averaged on the cell centers before their energy is
        computed.

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

    * ``FieldProbe``
        This type computes the value of each component of the electric and magnetic fields
        and of the Poynting vector (a measure of electromagnetic flux) at points in the domain.

        Multiple geometries for point probes can be specified via ``<reduced_diags_name>.probe_geometry = ...``:

        * ``Point`` (default): a single point
        * ``Line``: a line of points with equal spacing
        * ``Plane``: a plane of points with equal spacing

        **Point**: The point where the fields are measured is specified through the input parameters ``<reduced_diags_name>.x_probe``, ``<reduced_diags_name>.y_probe`` and ``<reduced_diags_name>.z_probe``.

        **Line**: probe a 1 dimensional line of points to create a line detector.
        Initial input parameters ``x_probe``, ``y_probe``, and ``z_probe`` designate one end of the line detector, while the far end is specified via ``<reduced_diags_name>.x1_probe``, ``<reduced_diags_name>.y1_probe``, ``<reduced_diags_name>.z1_probe``.
        Additionally, ``<reduced_diags_name>.resolution`` must be defined to give the number of detector points along the line (equally spaced) to probe.

        **Plane**: probe a 2 dimensional plane of points to create a square plane detector.
        Initial input parameters ``x_probe``, ``y_probe``, and ``z_probe`` designate the center of the detector.
        The detector plane is normal to a vector specified by ``<reduced_diags_name>.target_normal_x``, ``<reduced_diags_name>.target_normal_y``, and ``<reduced_diags_name>.target_normal_z``.
        Note that it is not necessary to specify the ``target_normal`` vector in a 2D simulation (the only supported normal is in ``y``).
        The top of the plane is perpendicular to an "up" vector denoted by ``<reduced_diags_name>.target_up_x``, ``<reduced_diags_name>.target_up_y``, and ``<reduced_diags_name>.target_up_z``.
        The detector has a square radius to be determined by ``<reduced_diags_name>.detector_radius``.
        Similarly to the line detector, the plane detector requires a resolution ``<reduced_diags_name>.resolution``, which denotes the number of detector particles along each side of the square detector.

        The output columns are
        the value of the :math:`E_x` field,
        the value of the :math:`E_y` field,
        the value of the :math:`E_z` field,
        the value of the :math:`B_x` field,
        the value of the :math:`B_y` field,
        the value of the :math:`B_z` field and
        the value of the Poynting Vector :math:`|S|` of the electromagnetic fields,
        at mesh refinement levels from  0 to :math:`n`, at point (:math:`x`, :math:`y`, :math:`z`).

        Note: the norms are always interpolated to the measurement point before they are written
        to file. The electromagnetic field components are interpolated to the measurement point
        by default, but can they be saved as non-averaged by setting
        ``<reduced_diags_name>.raw_fields = true``, in which case the raw fields for the cell
        containing the measurement point are saved.
        The interpolation order can be set by specifying ``<reduced_diags_name>.interp_order``,
        otherwise it is set to ``1``.
        Integrated electric and magnetic field components can instead be obtained by specifying
        ``<reduced_diags_name>.integrate == true``.
        In a *moving window* simulation, the FieldProbe can be set to follow the moving frame by specifying ``<reduced_diags_name>.do_moving_window_FP = 1`` (default 0).

        .. warning::

           The FieldProbe reduced diagnostic does not yet add a Lorentz back transformation for boosted frame simulations.
           Thus, it records field data in the boosted frame, not (yet) in the lab frame.

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

    * ``FieldReduction``
        This type computes an arbitrary reduction of the positions and the electromagnetic fields.

        * ``<reduced_diags_name>.reduced_function(x,y,z,Ex,Ey,Ez,Bx,By,Bz)`` (`string`)
            An analytic function to be reduced must be provided, using the math parser.

        * ``<reduced_diags_name>.reduction_type`` (`string`)
            The type of reduction to be performed. It must be either ``Maximum``, ``Minimum`` or
            ``Integral``.
            ``Integral`` computes the spatial integral of the function defined in the parser by
            summing its value on all grid points and multiplying the result by the volume of a
            cell.
            Please be also aware that measuring maximum quantities might be very noisy in PIC
            simulations.

        The only output column is the reduced value.

        Note that the fields are averaged on the cell centers before the reduction is performed.

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

    * ``LoadBalanceEfficiency``
        This type computes the load balance efficiency, given the present costs
        and distribution mapping. Load balance efficiency is computed as the
        mean cost over all ranks, divided by the maximum cost over all ranks.
        Until costs are recorded, load balance efficiency is output as `-1`;
        at earliest, the load balance efficiency can be output starting at step
        `2`, since costs are not recorded until step `1`.

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

        * ``<reduced_diags_name>.filter_function(t,x,y,z,ux,uy,uz)`` (`string`) optional
            Users can provide an expression returning a boolean for whether a particle is taken
            into account when calculating the histogram.
            `t` represents the physical time in seconds during the simulation.
            `x, y, z` represent particle positions in the unit of meter.
            `ux, uy, uz` represent particle velocities in the unit of
            :math:`\gamma v/c`, where
            :math:`\gamma` is the Lorentz factor,
            :math:`v/c` is the particle velocity normalized by the speed of light.
            E.g. If provided `(x>0.0)*(uz<10.0)` only those particles located at
            positions `x` greater than `0`, and those having velocity `uz` less than 10,
            will be taken into account when calculating the histogram.

        The output columns are
        values of the 1st bin, the 2nd bin, ..., the nth bin.
        An example input file and a loading python script of
        using the histogram reduced diagnostics
        are given in ``Examples/Tests/initial_distribution/``.

    * ``ParticleExtrema``
        This type computes the minimum and maximum values of
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

* ``<reduced_diags_name>.intervals`` (`string`)
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
------------------------------------------------

Lookup tables store pre-computed values for functions used by the QED modules.
**This feature requires to compile with QED=TRUE (and also with QED_TABLE_GEN=TRUE for table generation)**

* ``qed_bw.lookup_table_mode`` (`string`)
    There are three options to prepare the lookup table required by the Breit-Wheeler module:

    * ``builtin``:  a built-in table is used (Warning: the table gives reasonable results but its resolution is quite low).

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

    * ``builtin``: a built-in table is used (Warning: the table gives reasonable results but its resolution is quite low).

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
    Schwinger process requires either ``warpx.do_nodal=1`` or
    ``algo.field_gathering=momentum-conserving`` (so that different field components are computed
    at the same location in the grid) and does not currently support mesh refinement, cylindrical
    coordinates or single precision.

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
For the 1 and 2 colon syntax, actually having values in the string is optional
(this means that ``::5``, ``100 ::10`` and ``100 :`` are all valid syntaxes).

All values can be expressions that will be parsed in the same way as other integer input parameters.

**Examples**

* ``something_intervals = 50`` -> do something at timesteps 0, 50, 100, 150, etc.
  (equivalent to ``something_intervals = ::50``)

* ``something_intervals = 300:600:100`` -> do something at timesteps 300, 400, 500 and 600.

* ``something_intervals = 300::50`` -> do something at timesteps 300, 350, 400, 450, etc.

* ``something_intervals = 105:108,205:208`` -> do something at timesteps 105, 106, 107, 108,
  205, 206, 207 and 208. (equivalent to ``something_intervals = 105 : 108 : , 205 : 208 :``)

* ``something_intervals = :`` or  ``something_intervals = ::`` -> do something at every timestep.

* ``something_intervals = 167:167,253:253,275:425:50`` do something at timesteps 167, 253, 275,
  325, 375 and 425.

This is essentially the python slicing syntax except that the stop is inclusive
(``0:100`` contains 100) and that no colon means that the given value is the period.

Note that if a given period is zero or negative, the corresponding slice is disregarded.
For example, ``something_intervals = -1`` deactivates ``something`` and
``something_intervals = ::-1,100:1000:25`` is equivalent to ``something_intervals = 100:1000:25``.
