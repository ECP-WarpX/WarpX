/* Copyright 2022-2024 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: S. Eric Clark (LLNL), Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */
#include "FieldSolver/MagnetostaticSolver/MagnetostaticSolver.H"
#include "Particles/MultiParticleContainer.H"
#include "Python/callbacks.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <ablastr/coarsen/sample.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/fields/VectorPoissonSolver.H>


using namespace amrex;

void
WarpX::ComputeInitialMagnetostaticField()
{
    // First check setup appropriateness
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        this->max_level == 0,
        "Magnetostatic solver not implemented with mesh refinement."
    );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        do_magnetostatic_solve,
        "Magnetostatic solver called but warpx.do_magnetostatic is false."
    );
#ifdef WARPX_DIM_RZ
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        n_rz_azimuthal_modes == 1,
        "Error: RZ magnetostatic only implemented for a single mode"
    );
#endif

    // This function is specially setup for the initial magnetostatic solve
    // where particles and their velocities are synchronized. The goal is to
    // calculate the plasma contribution to the B-field. It is assumed that the
    // B-field is steady in time, i.e. dB/dt = 0 and therefore there is no
    // inductive E-field component.

    // Perform current deposition at t_0 (source of Poisson solver).
    mypc->DepositCurrent(current_fp, dt[0], 0.0_rt);

    // Synchronize J:
    // filter (if used), exchange guard cells, interpolate across MR levels
    // and apply boundary conditions
    SyncCurrent(current_fp, current_cp, current_buf);
    // SyncCurrent does not include a call to FillBoundary
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int idim = 0; idim < 3; ++idim) {
            current_fp[lev][idim]->FillBoundary(Geom(lev).periodicity());
        }
        ApplyJfieldBoundary(
            lev, current_fp[lev][0].get(), current_fp[lev][1].get(),
            current_fp[lev][2].get(), PatchType::fine
        );
    }

    // Reference magnetostatic multifabs
    auto& Afield_fp_nodal = m_magnetostatic_solver->Afield_fp_nodal;
    auto& Afield_fp_old = m_magnetostatic_solver->Afield_fp_old;
    auto& current_fp_temp = m_magnetostatic_solver->current_fp_temp;

    // Copy current_fp values to current_fp_temp which stores the currents at
    // the previous time step. This in line with the assumption of initial
    // steady state.
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab::Copy(
                *current_fp_temp[lev][idim], *current_fp[lev][idim],
                0, 0, 1, current_fp_temp[lev][idim]->nGrowVect()
            );
        }
    }

    // set the boundary and current density potentials
    setVectorPotentialBC(m_magnetostatic_solver->Afield_fp_nodal);

    // Compute the vector potential A, by solving the Poisson equation
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( !IsPythonCallbackInstalled("poissonsolver"),
        "Python Level Poisson Solve not supported for Magnetostatic implementation.");

    computeVectorPotential(
        current_fp_temp, Afield_fp_nodal,
        m_magnetostatic_solver->required_precision,
        m_magnetostatic_solver->absolute_tolerance,
        m_magnetostatic_solver->max_iters,
        m_magnetostatic_solver->verbosity
    );

    // interpolate A from a nodal multifab to an appropriately staggered one
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        m_magnetostatic_solver->InterpolateAfield(Afield_fp[lev], lev);
        Afield_fp[lev][0]->FillBoundary(Geom(lev).periodicity());
        Afield_fp[lev][1]->FillBoundary(Geom(lev).periodicity());
        Afield_fp[lev][2]->FillBoundary(Geom(lev).periodicity());
    }

    // Copy Afield_fp values to Afield_fp_old, such that the initial dA/dt is 0.
    for (int lev = 0; lev <= maxLevel(); lev++)
    {
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab::Copy(
                *Afield_fp_old[lev][idim], *Afield_fp[lev][idim],
                0, 0, 1, Afield_fp_old[lev][idim]->nGrowVect()
            );
        }
    }

    // At this point Afield_fp contains the vector potential at t = 0 and
    // we are ready to obtain B^0.
    ComputeBfromVectorPotential();
}

void
WarpX::ComputeMagnetostaticField()
{
    WARPX_PROFILE("WarpX::ComputeMagnetostaticField");

    // Fields have been reset in Electrostatic solver for this time step, these
    // fields are added into the B and E fields after electrostatic solve
    ComputeMagnetostaticFieldLabFrame();
}

void
WarpX::ComputeMagnetostaticFieldLabFrame()
{
    WARPX_PROFILE("WarpX::ComputeMagnetostaticFieldLabFrame");

    // Perform current deposition at t_{n+1/2} (source of Poisson solver).
    mypc->DepositCurrent(current_fp, dt[0], -0.5_rt * dt[0]);

    // Synchronize J:
    // filter (if used), exchange guard cells, interpolate across MR levels
    // and apply boundary conditions
    SyncCurrent(current_fp, current_cp, current_buf);
    // SyncCurrent does not include a call to FillBoundary
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (int idim = 0; idim < 3; ++idim) {
            current_fp[lev][idim]->FillBoundary(Geom(lev).periodicity());
        }
        ApplyJfieldBoundary(
            lev, current_fp[lev][0].get(), current_fp[lev][1].get(),
            current_fp[lev][2].get(), PatchType::fine
        );
    }

    // Reference magnetostatic multifabs
    auto& Afield_fp_nodal = m_magnetostatic_solver->Afield_fp_nodal;
    auto& Afield_fp_old = m_magnetostatic_solver->Afield_fp_old;
    auto& current_fp_temp = m_magnetostatic_solver->current_fp_temp;

    // At this point J^{n-1/2} is stored in `current_fp_temp` and J^{n+1/2}
    // in `current_fp`. In order to calculate B^{n+1} we need to obtain
    // J^{n+1}. Use extrapolation for this:
    // J^{n+1} = 1/2 * J_i^{n-1/2} + 3/2 * J_i^{n+1/2}.
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int idim = 0; idim < 3; ++idim) {
            MultiFab::LinComb(
                *current_fp_temp[lev][idim],
                0.5_rt, *current_fp_temp[lev][idim], 0,
                1.5_rt, *current_fp[lev][idim], 0,
                0, 1, current_fp_temp[lev][idim]->nGrowVect()
            );
        }
    }

    // set the boundary and current density potentials
    setVectorPotentialBC(m_magnetostatic_solver->Afield_fp_nodal);

    // Compute the vector potential A, by solving the Poisson equation
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( !IsPythonCallbackInstalled("poissonsolver"),
        "Python Level Poisson Solve not supported for Magnetostatic implementation.");

    computeVectorPotential(
        current_fp_temp, Afield_fp_nodal,
        m_magnetostatic_solver->required_precision,
        m_magnetostatic_solver->absolute_tolerance,
        m_magnetostatic_solver->max_iters,
        m_magnetostatic_solver->verbosity
    );

    // interpolate A from a nodal multifab to an appropriately staggered one
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        m_magnetostatic_solver->InterpolateAfield(Afield_fp[lev], lev);
        Afield_fp[lev][0]->FillBoundary(Geom(lev).periodicity());
        Afield_fp[lev][1]->FillBoundary(Geom(lev).periodicity());
        Afield_fp[lev][2]->FillBoundary(Geom(lev).periodicity());

        // Time filter A to avoid high frequency noise in the B-field becoming
        // large E-fields. A low-pass filter is used in which the A-field
        // solution is constructed as a linear combination of the new and old
        // solutions.
        if (m_magnetostatic_solver->t_filter_param < 1.0_rt)
        {
            for (int idim = 0; idim < 3; ++idim) {
                MultiFab::LinComb(
                    *Afield_fp[lev][idim],
                    (1.0_rt - m_magnetostatic_solver->t_filter_param),
                    *Afield_fp_old[lev][idim], 0,
                    m_magnetostatic_solver->t_filter_param,
                    *Afield_fp[lev][idim], 0,
                    0, 1, Afield_fp[lev][idim]->nGrowVect()
                );
            }
        }
    }

    // At this point Afield_fp contains the vector potential at t = n+1 and
    // we are ready to obtain B^{n+1}.
    ComputeBfromVectorPotential();

    // update E-field from the B-field inductance, but only if time-filtering
    // of A is done otherwise the algorithm is unstable
    if (m_magnetostatic_solver->t_filter_param < 1.0_rt)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            m_magnetostatic_solver->UpdateEfromInductance(
                Efield_fp[lev], Afield_fp[lev], getdt(lev), lev
            );
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        for (int idim = 0; idim < 3; ++idim) {
            // Copy J^{n+1/2} to current_fp_temp so that it is available at the
            // next step as J^{n-1/2}
            MultiFab::Copy(*current_fp_temp[lev][idim], *current_fp[lev][idim],
                           0, 0, 1, current_fp_temp[lev][idim]->nGrowVect());
            // Copy A^{n+1/2} to Afield_fp_old so that is is available at the
            // next step as A^{n-1/2}
            MultiFab::Copy(*Afield_fp_old[lev][idim], *Afield_fp[lev][idim],
                           0, 0, 1, Afield_fp_old[lev][idim]->nGrowVect());
        }
    }
}

/* Compute the vector potential `A` by solving the Poisson equation with `J` as
   a source.
   This uses the amrex solver.

    More specifically, this solves the equation
    \f[
        \vec{\nabla}^2 r \vec{A} = - r \mu_0 \vec{J}
 \f]

   \param[in] curr The current density
   \param[out] A The vector potential to be computed by this function
   \param[in] required_precision The relative convergence threshold for the MLMG solver
   \param[in] absolute_tolerance The absolute convergence threshold for the MLMG solver
   \param[in] max_iters The maximum number of iterations allowed for the MLMG solver
   \param[in] verbosity The verbosity setting for the MLMG solver
*/
void
WarpX::computeVectorPotential (const amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>,3> >& curr,
                                amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>,3> >& A,
                                Real const required_precision,
                                Real absolute_tolerance,
                                int const max_iters,
                                int const verbosity) const
{
    // create a vector to our fields, sorted by level
    amrex::Vector<amrex::Array<amrex::MultiFab*,3>> sorted_curr;
    amrex::Vector<amrex::Array<amrex::MultiFab*,3>> sorted_A;
    for (int lev = 0; lev <= finest_level; ++lev) {
        sorted_curr.emplace_back(amrex::Array<amrex::MultiFab*,3> ({curr[lev][0].get(),
                                                                    curr[lev][1].get(),
                                                                    curr[lev][2].get()}));
        sorted_A.emplace_back(amrex::Array<amrex::MultiFab*,3> ({A[lev][0].get(),
                                                                 A[lev][1].get(),
                                                                 A[lev][2].get()}));
    }

#if defined(AMREX_USE_EB)
    amrex::Vector<amrex::EBFArrayBoxFactory const *> factories;
    for (int lev = 0; lev <= finest_level; ++lev) {
        factories.push_back(&WarpX::fieldEBFactory(lev));
    }
    const std::optional<amrex::Vector<amrex::EBFArrayBoxFactory const *> > eb_farray_box_factory({factories});
#else
    const std::optional<amrex::Vector<amrex::FArrayBoxFactory const *> > eb_farray_box_factory;
#endif

    ablastr::fields::computeVectorPotential(
        sorted_curr,
        sorted_A,
        required_precision,
        absolute_tolerance,
        max_iters,
        verbosity,
        this->geom,
        this->dmap,
        this->grids,
        m_magnetostatic_solver->m_boundary_handler,
        WarpX::do_single_precision_comms,
        this->ref_ratio,
        gett_new(0),
        eb_farray_box_factory
    );
}


/* \brief Set Dirichlet/Neumann boundary conditions for the magnetostatic solver.

    The given potential's values are fixed on the boundaries of the given
    dimension according to the desired values from the simulation input file,
    boundary.potential_lo and boundary.potential_hi.

   \param[inout] A The vector potential
*/
void
WarpX::setVectorPotentialBC ( amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>,3>>& A ) const
{
    // check if any dimension has non-periodic boundary conditions
    if (!m_magnetostatic_solver->m_boundary_handler.has_non_periodic) { return; }

    auto dirichlet_flag = m_magnetostatic_solver->m_boundary_handler.dirichlet_flag;

    // loop over all mesh refinement levels and set the boundary values
    for (int lev=0; lev <= max_level; lev++) {

        amrex::Box domain = Geom(lev).Domain();
        domain.surroundingNodes();
        for (int adim=0; adim < 3; adim++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(*A[lev][adim], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                // Extract the vector potential
                auto A_arr = A[lev][adim]->array(mfi);
                // Extract tileboxes for which to loop
                const Box& tb  = mfi.tilebox( A[lev][adim]->ixType().toIntVect());

                // loop over dimensions
                for (int idim=0; idim<AMREX_SPACEDIM; idim++){
                    // check if neither boundaries in this dimension should be set
                    if (!(dirichlet_flag[adim][2*idim] || dirichlet_flag[adim][2*idim+1])) { continue; }

                    // a check can be added below to test if the boundary values
                    // are already correct, in which case the ParallelFor over the
                    // cells can be skipped

                    if (!domain.strictly_contains(tb)) {
                        amrex::ParallelFor( tb,
                            [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                                IntVect iv(AMREX_D_DECL(i,j,k));

                                if (dirichlet_flag[adim][2*idim] && iv[idim] == domain.smallEnd(idim)) {
                                    A_arr(i,j,k) = 0.;
                                }

                                if (dirichlet_flag[adim][2*idim+1] && iv[idim] == domain.bigEnd(idim)) {
                                    A_arr(i,j,k) = 0.;
                                }
                            } // loop ijk
                        );
                    }
                } // idim
            } // MFIter
        }
    } // lev
}


MagnetostaticSolver::MagnetostaticSolver ( int nlevs_max )
{
    ReadParameters();
    AllocateMFs(nlevs_max);
}

void MagnetostaticSolver::ReadParameters ()
{
    const ParmParse pp_warpx("warpx");

    // Grab solver tolerances from input parameters
    utils::parser::queryWithParser(
        pp_warpx, "self_fields_required_precision", required_precision);
    utils::parser::queryWithParser(
        pp_warpx, "self_fields_absolute_tolerance", absolute_tolerance);
    absolute_tolerance *= PhysConst::c;
    utils::parser::queryWithParser(pp_warpx, "self_fields_max_iters", max_iters);

    pp_warpx.query("self_fields_verbosity", verbosity);

    utils::parser::queryWithParser(
        pp_warpx, "magnetostatic_t_filtering_parameter", t_filter_param);

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        t_filter_param > 0.0_rt,
        "Using `t_filter_parameter = 0` with the magnetostatic solver is "
        "equivalent to not using the magnetostatic solver at all."
    );
}

void MagnetostaticSolver::AllocateMFs (int nlevs_max)
{
    Afield_fp_nodal.resize(nlevs_max);
    Afield_fp_old.resize(nlevs_max);
    current_fp_temp.resize(nlevs_max);
}

void MagnetostaticSolver::AllocateLevelMFs (int lev, const BoxArray& ba, const DistributionMapping& dm,
                                            const int ncomps, const IntVect& ngJ, const IntVect& ngEB,
                                            const IntVect& Ex_nodal_flag,
                                            const IntVect& Ey_nodal_flag,
                                            const IntVect& Ez_nodal_flag)
{
    // The "Afield_fp_nodal" multifab is used to store the vector potential
    // as obtained from the linear solver.
    WarpX::AllocInitMultiFab(Afield_fp_nodal[lev][0], amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, ngJ, lev, "Afield_fp_nodal[x]", 0.0_rt);
    WarpX::AllocInitMultiFab(Afield_fp_nodal[lev][1], amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, ngJ, lev, "Afield_fp_nodal[y]", 0.0_rt);
    WarpX::AllocInitMultiFab(Afield_fp_nodal[lev][2], amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, ngJ, lev, "Afield_fp_nodal[z]", 0.0_rt);

    // The "Afield_fp_old" multifab is used to store the vector potential
    // from the previous solve.
    WarpX::AllocInitMultiFab(Afield_fp_old[lev][0], amrex::convert(ba, Ex_nodal_flag),
        dm, ncomps, ngEB, lev, "Afield_fp_old[x]", 0.0_rt);
    WarpX::AllocInitMultiFab(Afield_fp_old[lev][1], amrex::convert(ba, Ey_nodal_flag),
        dm, ncomps, ngEB, lev, "Afield_fp_old[y]", 0.0_rt);
    WarpX::AllocInitMultiFab(Afield_fp_old[lev][2], amrex::convert(ba, Ez_nodal_flag),
        dm, ncomps, ngEB, lev, "Afield_fp_old[z]", 0.0_rt);

    // The "current_fp_temp" multifab is used to store the current density
    // interpolated or extrapolated to appropriate timesteps.
    WarpX::AllocInitMultiFab(current_fp_temp[lev][0], amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, ngJ, lev, "current_fp_temp[x]", 0.0_rt);
    WarpX::AllocInitMultiFab(current_fp_temp[lev][1], amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, ngJ, lev, "current_fp_temp[y]", 0.0_rt);
    WarpX::AllocInitMultiFab(current_fp_temp[lev][2], amrex::convert(ba, IntVect(AMREX_D_DECL(1,1,1))),
        dm, ncomps, ngJ, lev, "current_fp_temp[z]", 0.0_rt);
}

void MagnetostaticSolver::ClearLevel (int lev)
{
    for (int i = 0; i < 3; ++i) {
        Afield_fp_nodal[lev][i].reset();
        Afield_fp_old[lev][i].reset();
        current_fp_temp[lev][i].reset();
    }
}

void MagnetostaticSolver::InitData () {
    // Store the boundary conditions for the field solver
    auto& warpx = WarpX::GetInstance();
    m_boundary_handler.defineVectorPotentialBCs(
        warpx.field_boundary_lo, warpx.field_boundary_hi, warpx.Geom(0).ProbLo(0)
    );

    // get appropriate staggering for interpolation of nodal A-field to the
    // E-field staggering (which already matches the set staggering for
    // Afield_fp_old).
    auto Ex_stag = Afield_fp_old[0][0]->ixType().toIntVect();
    auto Ey_stag = Afield_fp_old[0][1]->ixType().toIntVect();
    auto Ez_stag = Afield_fp_old[0][2]->ixType().toIntVect();

    // copy data to device
    for ( int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        E_IndexType[0][idim] = Ex_stag[idim];
        E_IndexType[1][idim] = Ey_stag[idim];
        E_IndexType[2][idim] = Ez_stag[idim];
    }

    // Below we set all the unused dimensions to have nodal values
    // since these values will be interpolated from a nodal grid - if this is
    // not done the Interp function returns nonsense values.
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ) || defined(WARPX_DIM_1D_Z)
    E_IndexType[0][2] = 1;
    E_IndexType[1][2] = 1;
    E_IndexType[2][2] = 1;
#endif
#if defined(WARPX_DIM_1D_Z)
    E_IndexType[0][1] = 1;
    E_IndexType[1][1] = 1;
    E_IndexType[2][1] = 1;
#endif
}

void MagnetostaticSolver::InterpolateAfield (
    amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>& Afield,
    const int lev
) {
    using namespace ablastr::coarsen::sample;

    // Parameters for `interp` that maps from Yee to nodal mesh and back
    amrex::GpuArray<int, 3> const& nodal = {1, 1, 1};
    // The "coarsening is just 1 i.e. no coarsening"
    amrex::GpuArray<int, 3> const& coarsen = {1, 1, 1};

    // interpolate Afield_fp_nodal[lev][i] to Afield[lev][i]
    for (int idim = 0; idim < 3; idim++) {
        auto dst_stag = E_IndexType[idim];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*Afield[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Array4<amrex::Real const> const& src_arr = Afield_fp_nodal[lev][idim]->const_array(mfi);
            Array4<amrex::Real> const& dst_arr = Afield[idim]->array(mfi);

            const Box bx = mfi.tilebox();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dst_arr(i, j, k) = Interp(src_arr, nodal, dst_stag, coarsen, i, j, k, 0);
            });
        }
    }
}

void MagnetostaticSolver::UpdateEfromInductance (
    amrex::Array<std::unique_ptr<amrex::MultiFab>, 3>& Efield,
    amrex::Array<std::unique_ptr<amrex::MultiFab>, 3> const& Afield,
    Real dt, const int lev
) {
    auto dt_inv = 1.0_rt / dt;
    // The E-field due to inductance (-dA/dt) is added to the E-field, which
    // already contains the capacitive part from the electrostatic potential.
    // Here we assume that the time derivative of A does not change over one
    // step since we only have access to A^{n} and A^{n+1}, but should
    // really find A^{n+3/2} in order to accurately calculate (dA/dt)^{n+1}.
    // Afield contains A^{n+1} and Afield_fp_old contains A^{n}.
    // E_ind = -dA/dt = -(Afield - Afield_fp_old) / dt
    for (int idim = 0; idim < 3; ++idim) {
        // Use Afield_fp_old as an intermediate storage location to calculate E_ind
        MultiFab::LinComb(
            *Afield_fp_old[lev][idim],
            dt_inv, *Afield_fp_old[lev][idim], 0,
            -dt_inv, *Afield[idim], 0,
            0, 1, Afield_fp_old[lev][idim]->nGrowVect()
        );
        // Add E_ind to existing E-field
        MultiFab::Add(
            *Efield[idim], *Afield_fp_old[lev][idim], 0, 0, 1,
            Afield_fp_old[lev][idim]->nGrowVect()
        );
    }
}
