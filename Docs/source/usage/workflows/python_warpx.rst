Accessing global WarpX functionalities (e.g., extract timestep)
---------------------------------------------------------------

An important object is ``sim.extension.warpx``, which is the Python equivalent to the
C++ ``WarpX`` simulation class and gives access to global functionalities:

.. py:class:: WarpX

   .. py:method:: getistep(lev: int)

      Get the current step on mesh-refinement level ``lev``.

   .. py:method:: gett_new(lev: int)

      Get the current physical time on mesh-refinement level ``lev``.

   .. py:method:: getdt(lev: int)

      Get the current physical time step size on mesh-refinement level ``lev``.

   .. py:method:: multi_particle_container

   .. py:method:: get_particle_boundary_buffer

   .. py:method:: implicit_solver

      Return the :py:class:`ImplicitSolver`, or ``None`` when the evolve scheme is explicit.

   .. py:method:: deposit_mass_matrices

      Zero and deposit the mass matrices from all species.
      The mass matrices are the linear response of the deposited current density to the electric field.
      They are only allocated by the implicit evolve schemes that use them
      (e.g. ``implicit_evolve.use_mass_matrices_jacobian = 1``).

   .. py:method:: sync_mass_matrices

      Sum the guard cells of the mass matrices into the valid cells.

   .. py:method:: set_potential_on_domain_boundary(potential_[lo/hi]_[x/y/z]: str)

      The potential on the domain boundaries can be modified when using the electrostatic solver.
      This function updates the strings and function parsers which set the domain
      boundary potentials during the Poisson solve.

   .. py:method:: set_potential_on_eb(potential: str)

      The embedded boundary (EB) conditions can be modified when using the electrostatic solver.
      This set the EB potential string and updates the function parser.

   .. py:method:: evolve(numsteps=-1)

      Evolve the simulation the specified number of steps.

   .. py:method:: step(numsteps=-1)

      An alias to the evolve method.

   .. autofunction:: pywarpx.picmi.Simulation.extension.finalize

The ``MultiParticleContainer`` returned by ``multi_particle_container`` (also available as ``sim.particles``)
holds all species and exposes, among others:

.. py:class:: MultiParticleContainer

   .. py:method:: get(name: str)

      Return the ``WarpXParticleContainer`` of the species ``name``.

   .. py:method:: push_p(lev: int, dt: float, Ex, Ey, Ez, Bx, By, Bz, momentum_push_type: str = "Full")

      Push the momentum of the particles of all species over ``dt``, leaving their positions unchanged.
      The fields are gathered from the given ``MultiFab`` objects, which must have their guard cells filled
      and the staggering of ``Efield_fp`` and ``Bfield_fp``.

The implicit evolve schemes expose their solver through ``implicit_solver``:

.. py:class:: ImplicitSolver

   .. py:method:: pre_linear_solve

      What the implicit solvers call before every linear solve: when mass matrices are in use, deposit them
      from the current particle state, mirror the symmetric halves of their diagonal blocks and save the
      electric field that they were linearized around.

   .. py:method:: finish_mass_matrices

      Mirror the symmetric half of the diagonal blocks of the mass matrices, which the deposition leaves out.

   .. py:method:: apply_mass_matrices(out_name: str, in_name: str, in_ref_name: str = None, baseline_name: str = None, scale: float = 1.0, zero_out_first: bool = False)

      Compute ``out += scale * S * (in - in_ref) [+ baseline]``, where ``S`` are the mass matrices and every
      field is looked up by name in the ``MultiFab`` register. ``in`` must have its guard cells filled.
      See `tests/unit/test_mass_matrices.py <https://github.com/BLAST-WarpX/warpx/blob/development/tests/unit/test_mass_matrices.py>`__
      for an example that checks the mass matrices against a particle push.
