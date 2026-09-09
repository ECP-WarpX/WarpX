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
