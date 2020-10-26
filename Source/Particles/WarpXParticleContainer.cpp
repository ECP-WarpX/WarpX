/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote
 * Jean-Luc Vay, Luca Fedeli, Maxence Thevenet
 * Michael Rowan, Remi Lehe, Revathi Jambunathan
 * Weiqun Zhang, Yinjian Zhao, levinem
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "MultiParticleContainer.H"
#include "WarpXParticleContainer.H"
#include "WarpX.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/CoarsenMR.H"
// Import low-level single-particle kernels
#include "Pusher/GetAndSetPosition.H"
#include "Pusher/UpdatePosition.H"
#include "Deposition/CurrentDeposition.H"
#include "Deposition/ChargeDeposition.H"

#include <AMReX_AmrParGDB.H>
#include <AMReX.H>

#include <limits>


using namespace amrex;

WarpXParIter::WarpXParIter (ContainerType& pc, int level)
    : amrex::ParIter<0,0,PIdx::nattribs>(pc, level,
             MFItInfo().SetDynamic(WarpX::do_dynamic_scheduling))
{
}

WarpXParIter::WarpXParIter (ContainerType& pc, int level, MFItInfo& info)
    : amrex::ParIter<0,0,PIdx::nattribs>(pc, level,
                   info.SetDynamic(WarpX::do_dynamic_scheduling))
{
}

WarpXParticleContainer::WarpXParticleContainer (AmrCore* amr_core, int ispecies)
    : ParticleContainer<0,0,PIdx::nattribs>(amr_core->GetParGDB())
    , species_id(ispecies)
{
    SetParticleSize();
    ReadParameters();

    // build up the map of string names to particle component numbers
    particle_comps["w"]  = PIdx::w;
    particle_comps["ux"] = PIdx::ux;
    particle_comps["uy"] = PIdx::uy;
    particle_comps["uz"] = PIdx::uz;
#ifdef WARPX_DIM_RZ
    particle_comps["theta"] = PIdx::theta;
#endif

    // Initialize temporary local arrays for charge/current deposition
    int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
    num_threads = omp_get_num_threads();
#endif
    local_rho.resize(num_threads);
    local_jx.resize(num_threads);
    local_jy.resize(num_threads);
    local_jz.resize(num_threads);
}

void
WarpXParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
        ParmParse pp("particles");

#ifdef AMREX_USE_GPU
        do_tiling = false; // By default, tiling is off on GPU
#else
        do_tiling = true;
#endif
        pp.query("do_tiling",  do_tiling);

        initialized = true;
    }
}

void
WarpXParticleContainer::AllocData ()
{
    // have to resize here, not in the constructor because grids have not
    // been built when constructor was called.
    reserveData();
    resizeData();
}

void
WarpXParticleContainer::AddNParticles (int /*lev*/,
                                       int n, const ParticleReal* x, const ParticleReal* y, const ParticleReal* z,
                                       const ParticleReal* vx, const ParticleReal* vy, const ParticleReal* vz,
                                       int nattr, const ParticleReal* attr, int uniqueparticles, amrex::Long id)
{
    // nattr is unused below but needed in the BL_ASSERT
    amrex::ignore_unused(nattr);

    BL_ASSERT(nattr == 1);

    const ParticleReal* weight = attr;

    int ibegin, iend;
    if (uniqueparticles) {
        ibegin = 0;
        iend = n;
    } else {
        int myproc = ParallelDescriptor::MyProc();
        int nprocs = ParallelDescriptor::NProcs();
        int navg = n/nprocs;
        int nleft = n - navg * nprocs;
        if (myproc < nleft) {
            ibegin = myproc*(navg+1);
            iend = ibegin + navg+1;
        } else {
            ibegin = myproc*navg + nleft;
            iend = ibegin + navg;
        }
    }

    //  Add to grid 0 and tile 0
    // Redistribute() will move them to proper places.
    auto& particle_tile = DefineAndReturnParticleTile(0, 0, 0);

    std::size_t np = iend-ibegin;

#ifdef WARPX_DIM_RZ
    Vector<ParticleReal> theta(np);
#endif

    for (int i = ibegin; i < iend; ++i)
    {
        ParticleType p;
        if (id==-1)
        {
            p.id() = ParticleType::NextID();
        } else {
            p.id() = id;
        }
        p.cpu() = ParallelDescriptor::MyProc();
#if (AMREX_SPACEDIM == 3)
        p.pos(0) = x[i];
        p.pos(1) = y[i];
        p.pos(2) = z[i];
#elif (AMREX_SPACEDIM == 2)
        amrex::ignore_unused(y);
#ifdef WARPX_DIM_RZ
        theta[i-ibegin] = std::atan2(y[i], x[i]);
        p.pos(0) = std::sqrt(x[i]*x[i] + y[i]*y[i]);
#else
        p.pos(0) = x[i];
#endif
        p.pos(1) = z[i];
#endif

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ){
            DefineAndReturnParticleTile(0, 0, 0);
        }

        particle_tile.push_back(p);
    }

    if (np > 0)
    {
        particle_tile.push_back_real(PIdx::w , weight + ibegin, weight + iend);
        particle_tile.push_back_real(PIdx::ux,     vx + ibegin,     vx + iend);
        particle_tile.push_back_real(PIdx::uy,     vy + ibegin,     vy + iend);
        particle_tile.push_back_real(PIdx::uz,     vz + ibegin,     vz + iend);

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ){
            DefineAndReturnParticleTile(0, 0, 0);
        }

        for (int comp = PIdx::uz+1; comp < PIdx::nattribs; ++comp)
        {
#ifdef WARPX_DIM_RZ
            if (comp == PIdx::theta) {
                particle_tile.push_back_real(comp, theta.data(), theta.data() + np);
            }
            else {
                particle_tile.push_back_real(comp, np, 0.0);
            }
#else
            particle_tile.push_back_real(comp, np, 0.0);
#endif
        }

        for (int i = PIdx::nattribs; i < NumRealComps(); ++i)
        {
            particle_tile.push_back_real(i, 0.0);
        }
    }

    Redistribute();
}

/* \brief Current Deposition for thread thread_num
 * \param pti         : Particle iterator
 * \param wp          : Array of particle weights
 * \param uxp uyp uzp : Array of particle
 * \param ion_lev      : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param jx jy jz    : Full array of current density
 * \param offset      : Index of first particle for which current is deposited
 * \param np_to_depose: Number of particles for which current is deposited.
                        Particles [offset,offset+np_tp_depose] deposit current
 * \param thread_num  : Thread number (if tiling)
 * \param lev         : Level of box that contains particles
 * \param depos_lev   : Level on which particles deposit (if buffers are used)
 * \param dt          : Time step for particle level
 */
void
WarpXParticleContainer::DepositCurrent(WarpXParIter& pti,
                                       RealVector& wp, RealVector& uxp,
                                       RealVector& uyp, RealVector& uzp,
                                       const int * const ion_lev,
                                       MultiFab* jx, MultiFab* jy, MultiFab* jz,
                                       const long offset, const long np_to_depose,
                                       int thread_num, int lev, int depos_lev,
                                       Real dt)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((depos_lev==(lev-1)) ||
                                     (depos_lev==(lev  )),
                                     "Deposition buffers only work for lev-1");

    // If no particles, do not do anything
    if (np_to_depose == 0) return;

    // If user decides not to deposit
    if (do_not_deposit) return;

    // Number of guard cells for local deposition of J
    WarpX& warpx = WarpX::GetInstance();
    const int ng_J = warpx.get_ng_depos_J().max();

    // Extract deposition order (same order along all directions) and check that
    // particles shape fits within the guard cells.
    // NOTE: In specific situations where the staggering of J and the current
    // deposition algorithm are not trivial, this check might be too relaxed
    // and we might include a particle that should deposit part of its current
    // in a neighboring box. However, this should catch particles traveling many
    // cells away, for example with algorithms that allow for large time steps.
    const int shape_extent = static_cast<int>(WarpX::nox / 2);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::numParticlesOutOfRange(pti, ng_J - shape_extent) == 0,
        "Particles shape does not fit within guard cells used for local current deposition");

    const std::array<Real,3>& dx = WarpX::CellSize(std::max(depos_lev,0));
    Real q = this->charge;

    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCurrent::CurrentDeposition", blp_deposit);
    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCurrent::Accumulate", blp_accumulate);

    // Get tile box where current is deposited.
    // The tile box is different when depositing in the buffers (depos_lev<lev)
    // or when depositing inside the level (depos_lev=lev)
    Box tilebox;
    if (lev == depos_lev) {
        tilebox = pti.tilebox();
    } else {
        const IntVect& ref_ratio = WarpX::RefRatio(depos_lev);
        tilebox = amrex::coarsen(pti.tilebox(),ref_ratio);
    }

#ifndef AMREX_USE_GPU
    // Staggered tile boxes (different in each direction)
    Box tbx = convert( tilebox, jx->ixType().toIntVect() );
    Box tby = convert( tilebox, jy->ixType().toIntVect() );
    Box tbz = convert( tilebox, jz->ixType().toIntVect() );
#endif

    tilebox.grow(ng_J);

#ifdef AMREX_USE_GPU
    amrex::ignore_unused(thread_num);
    // GPU, no tiling: j<xyz>_arr point to the full j<xyz> arrays
    auto & jx_fab = jx->get(pti);
    auto & jy_fab = jy->get(pti);
    auto & jz_fab = jz->get(pti);
    Array4<Real> const& jx_arr = jx->array(pti);
    Array4<Real> const& jy_arr = jy->array(pti);
    Array4<Real> const& jz_arr = jz->array(pti);
#else
    tbx.grow(ng_J);
    tby.grow(ng_J);
    tbz.grow(ng_J);

    // CPU, tiling: j<xyz>_arr point to the local_j<xyz>[thread_num] arrays
    local_jx[thread_num].resize(tbx, jx->nComp());
    local_jy[thread_num].resize(tby, jy->nComp());
    local_jz[thread_num].resize(tbz, jz->nComp());

    // local_jx[thread_num] is set to zero
    local_jx[thread_num].setVal(0.0);
    local_jy[thread_num].setVal(0.0);
    local_jz[thread_num].setVal(0.0);

    auto & jx_fab = local_jx[thread_num];
    auto & jy_fab = local_jy[thread_num];
    auto & jz_fab = local_jz[thread_num];
    Array4<Real> const& jx_arr = local_jx[thread_num].array();
    Array4<Real> const& jy_arr = local_jy[thread_num].array();
    Array4<Real> const& jz_arr = local_jz[thread_num].array();
#endif

    const auto GetPosition = GetParticlePosition(pti, offset);

    // Lower corner of tile box physical domain
    // Note that this includes guard cells since it is after tilebox.ngrow
    const Dim3 lo = lbound(tilebox);
    // Take into account Galilean shift
    Real cur_time = warpx.gett_new(lev);
    const auto& time_of_last_gal_shift = warpx.time_of_last_gal_shift;
    Real time_shift = (cur_time + 0.5*dt - time_of_last_gal_shift);
    amrex::Array<amrex::Real,3> galilean_shift = {
        m_v_galilean[0]* time_shift,
        m_v_galilean[1]*time_shift,
        m_v_galilean[2]*time_shift };
    const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, galilean_shift, depos_lev);

    if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Esirkepov) {
        if (WarpX::do_nodal==1) {
          amrex::Abort("The Esirkepov algorithm cannot be used with a nodal grid.");
        }
        if ( (m_v_galilean[0]!=0) or (m_v_galilean[1]!=0) or (m_v_galilean[2]!=0)){
            amrex::Abort("The Esirkepov algorithm cannot be used with the Galilean algorithm.");
        }
    }

    WARPX_PROFILE_VAR_START(blp_deposit);
    if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Esirkepov) {
        if        (WarpX::nox == 1){
            doEsirkepovDepositionShapeN<1>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_arr, jy_arr, jz_arr, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 2){
            doEsirkepovDepositionShapeN<2>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_arr, jy_arr, jz_arr, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 3){
            doEsirkepovDepositionShapeN<3>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_arr, jy_arr, jz_arr, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes);
        }
    } else if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Vay) {
        if        (WarpX::nox == 1){
            doVayDepositionShapeN<1>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes );
        } else if (WarpX::nox == 2){
            doVayDepositionShapeN<2>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes );
        } else if (WarpX::nox == 3){
            doVayDepositionShapeN<3>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx, xyzmin, lo, q,
                WarpX::n_rz_azimuthal_modes );
        }
    } else {
        if        (WarpX::nox == 1){
            doDepositionShapeN<1>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx,
                xyzmin, lo, q, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 2){
            doDepositionShapeN<2>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx,
                xyzmin, lo, q, WarpX::n_rz_azimuthal_modes);
        } else if (WarpX::nox == 3){
            doDepositionShapeN<3>(
                GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                jx_fab, jy_fab, jz_fab, np_to_depose, dt, dx,
                xyzmin, lo, q, WarpX::n_rz_azimuthal_modes);
        }
    }
    WARPX_PROFILE_VAR_STOP(blp_deposit);

#ifndef AMREX_USE_GPU
    // CPU, tiling: atomicAdd local_j<xyz> into j<xyz>
    WARPX_PROFILE_VAR_START(blp_accumulate);
    (*jx)[pti].atomicAdd(local_jx[thread_num], tbx, tbx, 0, 0, jx->nComp());
    (*jy)[pti].atomicAdd(local_jy[thread_num], tby, tby, 0, 0, jy->nComp());
    (*jz)[pti].atomicAdd(local_jz[thread_num], tbz, tbz, 0, 0, jz->nComp());
    WARPX_PROFILE_VAR_STOP(blp_accumulate);
#endif
}

/* \brief Charge Deposition for thread thread_num
 * \param pti         : Particle iterator
 * \param wp          : Array of particle weights
 * \param ion_lev     : Pointer to array of particle ionization level. This is
                         required to have the charge of each macroparticle
                         since q is a scalar. For non-ionizable species,
                         ion_lev is a null pointer.
 * \param rho         : Full array of charge density
 * \param icomp       : Component of rho into which charge is deposited.
                        0: old value (before particle push).
                        1: new value (after particle push).
 * \param offset      : Index of first particle for which charge is deposited
 * \param np_to_depose: Number of particles for which charge is deposited.
                        Particles [offset,offset+np_tp_depose] deposit charge
 * \param thread_num  : Thread number (if tiling)
 * \param lev         : Level of box that contains particles
 * \param depos_lev   : Level on which particles deposit (if buffers are used)
 */
void
WarpXParticleContainer::DepositCharge (WarpXParIter& pti, RealVector& wp,
                                       const int * const ion_lev,
                                       amrex::MultiFab* rho, int icomp,
                                       const long offset, const long np_to_depose,
                                       int thread_num, int lev, int depos_lev)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE((depos_lev==(lev-1)) ||
                                     (depos_lev==(lev  )),
                                     "Deposition buffers only work for lev-1");

    // If no particles, do not do anything
    if (np_to_depose == 0) return;

    // If user decides not to deposit
    if (do_not_deposit) return;

    // Number of guard cells for local deposition of rho
    WarpX& warpx = WarpX::GetInstance();
    const int ng_rho = warpx.get_ng_depos_rho().max();

    // Extract deposition order (same order along all directions) and check that
    // particles shape fits within the guard cells.
    // NOTE: In specific situations where the staggering of rho and the charge
    // deposition algorithm are not trivial, this check might be too strict and
    // we might need to relax it, as currently done for the current deposition.
    const int shape_extent = static_cast<int>(WarpX::nox / 2 + 1);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        amrex::numParticlesOutOfRange(pti, ng_rho - shape_extent) == 0,
        "Particles shape does not fit within guard cells used for local charge deposition");

    const std::array<Real,3>& dx = WarpX::CellSize(std::max(depos_lev,0));
    const Real q = this->charge;

    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCharge::ChargeDeposition", blp_ppc_chd);
    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCharge::Accumulate", blp_accumulate);

    // Get tile box where charge is deposited.
    // The tile box is different when depositing in the buffers (depos_lev<lev)
    // or when depositing inside the level (depos_lev=lev)
    Box tilebox;
    if (lev == depos_lev) {
        tilebox = pti.tilebox();
    } else {
        const IntVect& ref_ratio = WarpX::RefRatio(depos_lev);
        tilebox = amrex::coarsen(pti.tilebox(),ref_ratio);
    }

#ifndef AMREX_USE_GPU
    // Staggered tile box
    Box tb = amrex::convert( tilebox, rho->ixType().toIntVect() );
#endif

    tilebox.grow(ng_rho);

    const int nc = WarpX::ncomps;

#ifdef AMREX_USE_GPU
    amrex::ignore_unused(thread_num);
    // GPU, no tiling: rho_fab points to the full rho array
    MultiFab rhoi(*rho, amrex::make_alias, icomp*nc, nc);
    auto & rho_fab = rhoi.get(pti);
#else
    tb.grow(ng_rho);

    // CPU, tiling: rho_fab points to local_rho[thread_num]
    local_rho[thread_num].resize(tb, nc);

    // local_rho[thread_num] is set to zero
    local_rho[thread_num].setVal(0.0);

    auto & rho_fab = local_rho[thread_num];
#endif

    const auto GetPosition = GetParticlePosition(pti, offset);

    // Lower corner of tile box physical domain
    // Note that this includes guard cells since it is after tilebox.ngrow
    Real cur_time = warpx.gett_new(lev);
    Real dt = warpx.getdt(lev);
    const auto& time_of_last_gal_shift = warpx.time_of_last_gal_shift;
    // Take into account Galilean shift
    Real time_shift_rho_old = (cur_time - time_of_last_gal_shift);
    Real time_shift_rho_new = (cur_time + dt - time_of_last_gal_shift);
    amrex::Array<amrex::Real,3> galilean_shift;
    if (icomp==0){
        galilean_shift = {
            m_v_galilean[0]*time_shift_rho_old,
            m_v_galilean[1]*time_shift_rho_old,
            m_v_galilean[2]*time_shift_rho_old };
    } else{
        galilean_shift = {
            m_v_galilean[0]*time_shift_rho_new,
            m_v_galilean[1]*time_shift_rho_new,
            m_v_galilean[2]*time_shift_rho_new };
    }
    const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, galilean_shift, depos_lev);

    // Indices of the lower bound
    const Dim3 lo = lbound(tilebox);

    WARPX_PROFILE_VAR_START(blp_ppc_chd);
    if        (WarpX::nox == 1){
        doChargeDepositionShapeN<1>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                    rho_fab, np_to_depose, dx, xyzmin, lo, q,
                                    WarpX::n_rz_azimuthal_modes);
    } else if (WarpX::nox == 2){
        doChargeDepositionShapeN<2>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                    rho_fab, np_to_depose, dx, xyzmin, lo, q,
                                    WarpX::n_rz_azimuthal_modes);
    } else if (WarpX::nox == 3){
        doChargeDepositionShapeN<3>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                    rho_fab, np_to_depose, dx, xyzmin, lo, q,
                                    WarpX::n_rz_azimuthal_modes);
    }
    WARPX_PROFILE_VAR_STOP(blp_ppc_chd);

#ifndef AMREX_USE_GPU
    // CPU, tiling: atomicAdd local_rho into rho
    WARPX_PROFILE_VAR_START(blp_accumulate);
    (*rho)[pti].atomicAdd(local_rho[thread_num], tb, tb, 0, icomp*nc, nc);
    WARPX_PROFILE_VAR_STOP(blp_accumulate);
#endif
}

void
WarpXParticleContainer::DepositCharge (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                                        bool local, bool reset,
                                        bool do_rz_volume_scaling)
{
#ifdef WARPX_DIM_RZ
    (void)do_rz_volume_scaling;
#endif
    // Loop over the refinement levels
    int const finest_level = rho.size() - 1;
    for (int lev = 0; lev <= finest_level; ++lev) {

        // Reset the `rho` array if `reset` is True
        if (reset) rho[lev]->setVal(0.0, rho[lev]->nGrow());

        // Loop over particle tiles and deposit charge on each level
#ifdef _OPENMP
#pragma omp parallel
        {
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const long np = pti.numParticles();
            auto& wp = pti.GetAttribs(PIdx::w);

            int* AMREX_RESTRICT ion_lev;
            if (do_field_ionization){
                ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
            } else {
                ion_lev = nullptr;
            }

            DepositCharge(pti, wp, ion_lev, rho[lev].get(), 0, 0, np, thread_num, lev, lev);
        }
#ifdef _OPENMP
        }
#endif

#ifdef WARPX_DIM_RZ
        if (do_rz_volume_scaling) {
            WarpX::GetInstance().ApplyInverseVolumeScalingToChargeDensity(rho[lev].get(), lev);
        }
#else
        ignore_unused(do_rz_volume_scaling);
#endif

        // Exchange guard cells
        if (!local) rho[lev]->SumBoundary( m_gdb->Geom(lev).periodicity() );
    }

    // Now that the charge has been deposited at each level,
    // we average down from fine to crse
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        const DistributionMapping& fine_dm = rho[lev+1]->DistributionMap();
        BoxArray coarsened_fine_BA = rho[lev+1]->boxArray();
        coarsened_fine_BA.coarsen(m_gdb->refRatio(lev));
        MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, rho[lev+1]->nComp(), 0);
        coarsened_fine_data.setVal(0.0);

        int const refinement_ratio = 2;

        CoarsenMR::Coarsen( coarsened_fine_data, *rho[lev+1], IntVect(refinement_ratio) );
        rho[lev]->ParallelAdd( coarsened_fine_data, m_gdb->Geom(lev).periodicity() );
    }
}

std::unique_ptr<MultiFab>
WarpXParticleContainer::GetChargeDensity (int lev, bool local)
{
    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);
    const auto& dm = m_gdb->DistributionMap(lev);
    BoxArray nba = ba;
#if (!defined WARPX_DIM_RZ) || (!defined WARPX_USE_PSATD)
    nba.surroundingNodes();
#endif

    // Number of guard cells for local deposition of rho
    WarpX& warpx = WarpX::GetInstance();
    const int ng_rho = warpx.get_ng_depos_rho().max();

    auto rho = std::make_unique<MultiFab>(nba,dm,WarpX::ncomps,ng_rho);
    rho->setVal(0.0);

#ifdef _OPENMP
#pragma omp parallel
    {
#endif
#ifdef _OPENMP
        int thread_num = omp_get_thread_num();
#else
        int thread_num = 0;
#endif

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const long np = pti.numParticles();
            auto& wp = pti.GetAttribs(PIdx::w);

            int* AMREX_RESTRICT ion_lev;
            if (do_field_ionization){
                ion_lev = pti.GetiAttribs(particle_icomps["ionization_level"]).dataPtr();
            } else {
                ion_lev = nullptr;
            }

            DepositCharge(pti, wp, ion_lev, rho.get(), 0, 0, np,
                          thread_num, lev, lev);
        }
#ifdef _OPENMP
    }
#endif

#ifdef WARPX_DIM_RZ
    WarpX::GetInstance().ApplyInverseVolumeScalingToChargeDensity(rho.get(), lev);
#endif

    if (!local) rho->SumBoundary(gm.periodicity());

    return rho;
}

Real WarpXParticleContainer::sumParticleCharge(bool local) {

    amrex::Real total_charge = 0.0;

    const int nLevels = finestLevel();
    for (int lev = 0; lev < nLevels; ++lev)
    {

#ifdef _OPENMP
#pragma omp parallel reduction(+:total_charge)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& wp = pti.GetAttribs(PIdx::w);
            for (unsigned long i = 0; i < wp.size(); i++) {
                total_charge += wp[i];
            }
        }
    }

    if (!local) ParallelDescriptor::ReduceRealSum(total_charge);
    total_charge *= this->charge;
    return total_charge;
}

std::array<Real, 3> WarpXParticleContainer::meanParticleVelocity(bool local) {

    amrex::Real vx_total = 0.0;
    amrex::Real vy_total = 0.0;
    amrex::Real vz_total = 0.0;

    long np_total = 0;

    amrex::Real inv_clight_sq = 1.0/PhysConst::c/PhysConst::c;

    const int nLevels = finestLevel();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<Real, Real, Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (int lev = 0; lev <= nLevels; ++lev) {
            for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
                const auto uxp = pti.GetAttribs(PIdx::ux).data();
                const auto uyp = pti.GetAttribs(PIdx::uy).data();
                const auto uzp = pti.GetAttribs(PIdx::uz).data();

                const long np = pti.numParticles();
                np_total += np;

                reduce_op.eval(np, reduce_data,
                               [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                               {
                                   Real usq = (uxp[i]*uxp[i] +
                                               uyp[i]*uyp[i] +
                                               uzp[i]*uzp[i])*inv_clight_sq;
                                   Real gaminv = 1.0_rt/std::sqrt(1.0_rt + usq);
                                   return {uxp[i]*gaminv,  uyp[i]*gaminv, uzp[i]*gaminv};
                               });
            }
        }

        ReduceTuple hv = reduce_data.value();
        vx_total = amrex::get<0>(hv);
        vy_total = amrex::get<1>(hv);
        vz_total = amrex::get<2>(hv);
    }
    else
#endif
    {
        for (int lev = 0; lev <= nLevels; ++lev) {
#ifdef _OPENMP
#pragma omp parallel reduction(+:vx_total, vy_total, vz_total, np_total)
#endif
            for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
                auto& ux = pti.GetAttribs(PIdx::ux);
                auto& uy = pti.GetAttribs(PIdx::uy);
                auto& uz = pti.GetAttribs(PIdx::uz);

                np_total += pti.numParticles();

                for (unsigned long i = 0; i < ux.size(); i++) {
                    Real usq = (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i])*inv_clight_sq;
                    Real gaminv = 1.0_rt/std::sqrt(1.0_rt + usq);
                    vx_total += ux[i]*gaminv;
                    vy_total += uy[i]*gaminv;
                    vz_total += uz[i]*gaminv;
                }
            }
        }
    }

    if (!local) {
        ParallelDescriptor::ReduceRealSum(vx_total);
        ParallelDescriptor::ReduceRealSum(vy_total);
        ParallelDescriptor::ReduceRealSum(vz_total);
        ParallelDescriptor::ReduceLongSum(np_total);
    }

    std::array<Real, 3> mean_v;
    if (np_total > 0) {
        mean_v[0] = vx_total / np_total;
        mean_v[1] = vy_total / np_total;
        mean_v[2] = vz_total / np_total;
    }

    return mean_v;
}

Real WarpXParticleContainer::maxParticleVelocity(bool local) {

    amrex::ParticleReal max_v = 0.0;

    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev)
    {

#ifdef _OPENMP
#pragma omp parallel reduction(max:max_v)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& ux = pti.GetAttribs(PIdx::ux);
            auto& uy = pti.GetAttribs(PIdx::uy);
            auto& uz = pti.GetAttribs(PIdx::uz);
            for (unsigned long i = 0; i < ux.size(); i++) {
                max_v = std::max(max_v, std::sqrt(ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i]));
            }
        }
    }

    if (!local) ParallelAllReduce::Max(max_v, ParallelDescriptor::Communicator());
    return max_v;
}

void
WarpXParticleContainer::PushX (amrex::Real dt)
{
    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev) {
        PushX(lev, dt);
    }
}

void
WarpXParticleContainer::PushX (int lev, amrex::Real dt)
{
    WARPX_PROFILE("WarpXParticleContainer::PushX()");

    if (do_not_push) return;

    amrex::LayoutData<amrex::Real>* cost = WarpX::getCosts(lev);

#ifdef _OPENMP
#pragma omp parallel
#endif
    {

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
            }
            Real wt = amrex::second();

            //
            // Particle Push
            //

            const auto GetPosition = GetParticlePosition(pti);
                  auto SetPosition = SetParticlePosition(pti);

            // - momenta are stored as a struct of array, in `attribs`
            auto& attribs = pti.GetAttribs();
            ParticleReal* AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr();
            ParticleReal* AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr();
            ParticleReal* AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr();

            // Loop over the particles and update their position
            amrex::ParallelFor( pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                                    ParticleReal x, y, z;
                                    GetPosition(i, x, y, z);
                                    UpdatePosition(x, y, z, ux[i], uy[i], uz[i], dt);
                                    SetPosition(i, x, y, z);
                }
            );

            if (cost && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add( &(*cost)[pti.index()], wt);
            }
        }
    }
}

// When using runtime components, AMReX requires to touch all tiles
// in serial and create particles tiles with runtime components if
// they do not exist (or if they were defined by default, i.e.,
// without runtime component).
void WarpXParticleContainer::defineAllParticleTiles () noexcept
{
    tmp_particle_data.resize(finestLevel()+1);
    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
        for (auto mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            tmp_particle_data[lev][std::make_pair(grid_id,tile_id)];
            DefineAndReturnParticleTile(lev, grid_id, tile_id);
        }
    }
}

// This function is called in Redistribute, just after locate
void
WarpXParticleContainer::particlePostLocate(ParticleType& p,
                                           const ParticleLocData& pld,
                                           const int lev)
{
    if (not do_splitting) return;

    // Tag particle if goes to higher level.
    // It will be split later in the loop
    if (pld.m_lev == lev+1
        and p.id() != NoSplitParticleID
        and p.id() >= 0)
    {
        p.id() = DoSplitParticleID;
    }

    if (pld.m_lev == lev-1){
        // For the moment, do not do anything if particles goes
        // to lower level.
    }
}

void
WarpXParticleContainer::ApplyBoundaryConditions (ParticleBC boundary_conditions){
    WARPX_PROFILE("WarpXParticleContainer::ApplyBoundaryConditions()");
    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto GetPosition = GetParticlePosition(pti);
            const Real xmin = Geom(lev).ProbLo(0);
            const Real xmax = Geom(lev).ProbHi(0);
#ifdef WARPX_DIM_3D
            const Real ymin = Geom(lev).ProbLo(1);
            const Real ymax = Geom(lev).ProbHi(1);
#endif
            const Real zmin = Geom(lev).ProbLo(AMREX_SPACEDIM-1);
            const Real zmax = Geom(lev).ProbHi(AMREX_SPACEDIM-1);

            ParticleTileType& ptile = ParticlesAt(lev, pti);
            ParticleType * const pp = ptile.GetArrayOfStructs()().data();

            // Loop over particles and apply BC to each particle
            amrex::ParallelFor(
                pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i) {
                    ParticleType& p = pp[i];
                    ParticleReal x, y, z;
                    GetPosition(i, x, y, z);
#ifdef WARPX_DIM_3D
                    if (x < xmin || x > xmax || y < ymin || y > ymax || z < zmin || z > zmax){
                        if (boundary_conditions == ParticleBC::absorbing) p.id() = -1;
                    }
#else
                    if (x < xmin || x > xmax || z < zmin || z > zmax){
                        if (boundary_conditions == ParticleBC::absorbing) p.id() = -1;
                    }
#endif
                }
            );
        }
    }
}
