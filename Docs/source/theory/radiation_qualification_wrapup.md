# Radiation qualification checkpoint: WIP

This update publishes the current qualification work, not a declaration of
production readiness. It does not merge the PR, include circuit-coupling packs,
or replace kinetic ions with a FLASH fluid closure. AMR is not a completion
requirement for this work.

## Implemented scope

- Opt-in stationary implicit diffusion and native electron caloric exchange,
  with private candidate state, residual checks and bounded substep retries.
- Opt-in gray moving moment transport in periodic Cartesian geometry, coupled
  to native kinetic ions and hybrid electron energy. Actual finite ion work and
  deferred particle-owned impulse/work are separately accounted and restarted.
- Native shape-consistent force assignment and adjoint work, with independent
  moving-beam ODE, trapped-pulse transport and oblique-interface checks.
- Particle diagnostics use output copies for unit conversion, preserving live
  particle state. Hybrid deposited moment history is preserved on restart.
- Guarded boundary building blocks: positive M1 mirror flux, independent face
  accounting, reflecting-wall force assignment and transactional native
  particle-carry reflection with a checkpointed pending-only wall ledger.

Existing legacy defaults remain. Experimental options and unsupported-path
guards are documented in the input reference. In particular, these boundary
building blocks do **not** enable the full moving radiation/material wall path.

## Evidence and unresolved gates

The public-branch CPU regression subset passes 328/328 checks after restoring
the existing parent carry-ownership validation ahead of its new wall ledger.
This includes the legacy radiation cases and the new implicit/moving exchange,
boundary, rejection, diagnostic and restart checks. The full-run log is
`build-production-cpu/pr-owner-full-tests.log` in the PR qualification worktree.

The latest carry-boundary tests passed on CPU with one rank, two ranks and a
one-to-two-rank checkpoint continuation. They include exact live-state rejection,
normal/reflect-all reflection, two carry owners, migration, overflow rejection,
and compensated cancellation across restart. CUDA qualification of this latest
increment was outstanding when it was packaged.

Earlier CPU/P40 qualification covers the periodic moving beam, trapped pulse,
shape assignment and several independent boundary kernels. Those results do
not qualify all dimensions or all boundary combinations. The strict full cold
material-wall GPU trajectory and cross-backend bounds still fail, even though
total energy conservation and restored deposited-history checks pass. No
physics tolerance has been loosened to turn those failures into passes.

Remaining work includes the integrated physical-boundary material feedback and
complete represented-plus-deferred wall inventory, resolution of those strict
GPU trajectory gates, moving multigroup frequency redistribution, packet/moment
conversion with moving material, and practical full 3D/RZ qualification.
Green CI is a regression/build milestone, not evidence that these open physics
requirements are complete.

### CI checksum follow-up

The 2D CI matrix also exposed a real initialization defect: native bulk plasma
injection does not call the generic runtime-attribute initializer, so newly
registered radiation carries could contain allocator residue. Both bulk
injection paths now explicitly zero newborn radiation accounts. The flux path
does so after redistribution of its temporary container, whose runtime
components are not communicated. Existing particle accounts are not reset.
The new pre-filled-capacity birth test passes with the fix and fails when the
initializer is removed. The affected 2D moment/particle-carry subset then
passes 77/77 locally, including MPI, rollback and restart cases.
The expanded 1D radiation/coupling subset also passes 329/329 checks
(`build-ci-tools/birth-final-1d-tests.log`).

Single-precision CI compilation additionally requires native-precision flux
literals and explicit ownership/query declarations. Those changes retain the
moving runtime's double-precision guard. Invalid-tolerance rejection now has
explicit zero, negative, out-of-range, infinity and NaN coverage; no accepted
tolerance range or physics assertion was relaxed.

Azure build 6472 passed the RZ pinch simulation and its original physics
analysis, but its existing `test_rz_theta_implicit_dynamic_pinch` checksum
differed by up to 1.58e-6. A local PETSc/MPI A/B comparison changes only the
diagnostic unit-conversion behavior: restoring the old live SI round trip
changes checksums by up to 4.80e-6 relative to the corrected copy-only writer.
Both variants pass the original 1e-12 energy and charge-conservation bounds
and the solver-iteration limits. Corrected/legacy maximum relative energy
errors are 1.97e-15 / 3.10e-14; charge RMS errors are 8.45e-14 / 8.00e-14.

Neither local variant reproduces the exact CI checksum. This establishes
sensitivity to diagnostic mutation, not unique attribution of the CI delta;
the local PETSc/MPI environment differs from CI. No generated benchmark or
physics tolerance was changed. Artifacts are under `build-rz-checksum/` and
`build-ci-tools/rz-diagnostic-comparison.json` in the PR worktree. The temporary
legacy writer was removed immediately after building the control executable.

CodeQL also identified two reference-test products evaluated in `double`
before conversion to `long double`. Promoting the operands before arithmetic
fixes those intermediate-precision issues; the existing moving-flux test
passes without changing its bound.

## Detailed records

```{toctree}
:maxdepth: 1

radiation_moving_runtime
radiation_moving_material_plan
radiation_moving_moment_discretization
radiation_moving_beam_verification
radiation_particle_shape_exchange
radiation_moment_boundary_plan
radiation_material_wall_qualification
radiation_oblique_interface_qualification
radiation_reflecting_cavity_qualification
```

Build/test artifacts referenced in these records are local qualification logs,
not bundled simulation data. No FLASH source or remote-machine credentials are
included in this update.
