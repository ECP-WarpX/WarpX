/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote
 * Jean-Luc Vay, Luca Fedeli, Maxence Thevenet
 * Michael Rowan, Remi Lehe, Revathi Jambunathan
 * Weiqun Zhang, Yinjian Zhao, levinem
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXParticleContainer.H"

#include "ablastr/particles/DepositCharge.H"
#include "Deposition/ChargeDeposition.H"
#include "Deposition/CurrentDeposition.H"
#include "Deposition/SharedDepositionUtils.H"
#include "Pusher/GetAndSetPosition.H"
#include "Pusher/UpdatePosition.H"
#include "ParticleBoundaries_K.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <ablastr/coarsen/average.H>
#include <ablastr/utils/Communication.H>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_AmrParGDB.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_Dim3.H>
#include <AMReX_Extension.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuAllocators.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParGDB.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleContainerBase.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParticleTransformation.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_Random.H>
#include <AMReX_TinyProfiler.H>
#include <AMReX_Utility.H>


#ifdef AMREX_USE_OMP
#   include <omp.h>
#endif

#include <algorithm>
#include <cmath>

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
    : NamedComponentParticleContainer<DefaultAllocator>(amr_core->GetParGDB())
    , species_id(ispecies)
{
    SetParticleSize();
    ReadParameters();

    // Initialize temporary local arrays for charge/current deposition
#ifdef AMREX_USE_OMP
    int num_threads = 1;
#pragma omp parallel
#pragma omp single
    num_threads = omp_get_num_threads();
#else
    const int num_threads = 1;
#endif

    local_rho.resize(num_threads);
    local_jx.resize(num_threads);
    local_jy.resize(num_threads);
    local_jz.resize(num_threads);

    // The boundary conditions are read in in ReadBCParams but a child class
    // can allow these value to be overwritten if different boundary
    // conditions are desired for a specific species
#ifndef WARPX_DIM_1D_Z
    m_boundary_conditions.SetBoundsX(WarpX::particle_boundary_lo[0], WarpX::particle_boundary_hi[0]);
#endif
#ifdef WARPX_DIM_3D
    m_boundary_conditions.SetBoundsY(WarpX::particle_boundary_lo[1], WarpX::particle_boundary_hi[1]);
    m_boundary_conditions.SetBoundsZ(WarpX::particle_boundary_lo[2], WarpX::particle_boundary_hi[2]);
#elif WARPX_DIM_XZ || WARPX_DIM_RZ
    m_boundary_conditions.SetBoundsZ(WarpX::particle_boundary_lo[1], WarpX::particle_boundary_hi[1]);
#else
    m_boundary_conditions.SetBoundsZ(WarpX::particle_boundary_lo[0], WarpX::particle_boundary_hi[0]);
#endif
    m_boundary_conditions.BuildReflectionModelParsers();
}

void
WarpXParticleContainer::ReadParameters ()
{
    static bool initialized = false;
    if (!initialized)
    {
        const ParmParse pp_particles("particles");
        pp_particles.query("do_tiling", do_tiling);
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
                                       int n, const amrex::ParticleReal* x,
                                       const amrex::ParticleReal* y,
                                       const amrex::ParticleReal* z,
                                       const amrex::ParticleReal* ux,
                                       const amrex::ParticleReal* uy,
                                       const amrex::ParticleReal* uz,
                                       const int nattr_real, const amrex::ParticleReal* attr_real,
                                       const int nattr_int, const int* attr_int,
                                       int uniqueparticles, amrex::Long id)
{
    using namespace amrex::literals;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE((PIdx::nattribs + nattr_real - 1) <= NumRealComps(),
                                     "Too many real attributes specified");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(nattr_int <= NumIntComps(),
                                     "Too many integer attributes specified");

    int ibegin, iend;
    if (uniqueparticles) {
        ibegin = 0;
        iend = n;
    } else {
        const int myproc = amrex::ParallelDescriptor::MyProc();
        const int nprocs = amrex::ParallelDescriptor::NProcs();
        const int navg = n/nprocs;
        const int nleft = n - navg * nprocs;
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

    using PinnedTile = amrex::ParticleTile<Particle<NStructReal, NStructInt>,
                                           NArrayReal, NArrayInt,
                                           amrex::PinnedArenaAllocator>;
    PinnedTile pinned_tile;
    pinned_tile.define(NumRuntimeRealComps(), NumRuntimeIntComps());

    const std::size_t np = iend-ibegin;

    // treat weight as a special attr since it will always be specified
    amrex::Vector<amrex::ParticleReal> weight(np);

#ifdef WARPX_DIM_RZ
    amrex::Vector<amrex::ParticleReal> theta(np);
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
        p.cpu() = amrex::ParallelDescriptor::MyProc();
#if defined(WARPX_DIM_3D)
        p.pos(0) = x[i];
        p.pos(1) = y[i];
        p.pos(2) = z[i];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        amrex::ignore_unused(y);
#ifdef WARPX_DIM_RZ
        theta[i-ibegin] = std::atan2(y[i], x[i]);
        p.pos(0) = std::sqrt(x[i]*x[i] + y[i]*y[i]);
#else
        p.pos(0) = x[i];
#endif
        p.pos(1) = z[i];
#else //AMREX_SPACEDIM == 1
        amrex::ignore_unused(x,y);
        p.pos(0) = z[i];
#endif

        pinned_tile.push_back(p);

        // grab weight from the attr array
        weight[i-ibegin] = attr_real[i*nattr_real];
    }

    if (np > 0)
    {
        pinned_tile.push_back_real(PIdx::w , weight.data(), weight.data() + np);
        pinned_tile.push_back_real(PIdx::ux,     ux + ibegin,     ux + iend);
        pinned_tile.push_back_real(PIdx::uy,     uy + ibegin,     uy + iend);
        pinned_tile.push_back_real(PIdx::uz,     uz + ibegin,     uz + iend);

        if ( (NumRuntimeRealComps()>0) || (NumRuntimeIntComps()>0) ){
            DefineAndReturnParticleTile(0, 0, 0);
        }

        for (int comp = PIdx::uz+1; comp < PIdx::nattribs; ++comp)
        {
#ifdef WARPX_DIM_RZ
            if (comp == PIdx::theta) {
                pinned_tile.push_back_real(comp, theta.data(), theta.data() + np);
            }
            else {
                pinned_tile.push_back_real(comp, np, 0.0_prt);
            }
#else
            pinned_tile.push_back_real(comp, np, 0.0_prt);
#endif
        }

        // Initialize nattr_real - 1 runtime real attributes from data in the attr_real array
        for (int j = PIdx::nattribs; j < PIdx::nattribs + nattr_real - 1; ++j)
        {
            // get the next attribute from attr_real array
            amrex::Vector<amrex::ParticleReal> attr_vals(np);
            for (int i = ibegin; i < iend; ++i)
            {
                attr_vals[i-ibegin] = attr_real[j - PIdx::nattribs + 1 + i*nattr_real];
            }
            pinned_tile.push_back_real(j, attr_vals.data(), attr_vals.data() + np);
        }

        // Initialize nattr_int runtime integer attributes from data in the attr_int array
        for (int j = 0; j < nattr_int; ++j)
        {
            // get the next attribute from attr_int array
            amrex::Vector<int> attr_vals(np);
            for (int i = ibegin; i < iend; ++i)
            {
                attr_vals[i-ibegin] = attr_int[j + i*nattr_int];
            }
            pinned_tile.push_back_int(j, attr_vals.data(), attr_vals.data() + np);
        }

        // Default initialize the other real and integer runtime attributes
        DefaultInitializeRuntimeAttributes(pinned_tile, nattr_real - 1, nattr_int,
                                           amrex::RandomEngine{});

        auto old_np = particle_tile.numParticles();
        auto new_np = old_np + pinned_tile.numParticles();
        particle_tile.resize(new_np);
        amrex::copyParticles(
            particle_tile, pinned_tile, 0, old_np, pinned_tile.numParticles()
        );
    }

    Redistribute();
}

/* \brief Current Deposition for thread thread_num
 * \param pti         : Particle iterator
 * \param wp          : Array of particle weights
 * \param uxp uyp uzp : Array of particle momenta
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
 * \param relative_time: Time at which to deposit J, relative to the time of the
 *                       current positions of the particles. When different than 0,
 *                       the particle position will be temporarily modified to match
 *                       the time of the deposition.
 */
void
WarpXParticleContainer::DepositCurrent (WarpXParIter& pti,
                                        RealVector const & wp, RealVector const & uxp,
                                        RealVector const & uyp, RealVector const & uzp,
                                        int const * const ion_lev,
                                        amrex::MultiFab * const jx, amrex::MultiFab * const jy, amrex::MultiFab * const jz,
                                        long const offset, long const np_to_depose,
                                        int const thread_num, const int lev, int const depos_lev,
                                        amrex::Real const dt, amrex::Real const relative_time)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE((depos_lev==(lev-1)) ||
                                     (depos_lev==(lev  )),
                                     "Deposition buffers only work for lev-1");

    // If no particles, do not do anything
    if (np_to_depose == 0) return;

    // If user decides not to deposit
    if (do_not_deposit) return;

    // Number of guard cells for local deposition of J
    const WarpX& warpx = WarpX::GetInstance();

    const amrex::IntVect& ng_J = warpx.get_ng_depos_J();

    // Extract deposition order and check that particles shape fits within the guard cells.
    // NOTE: In specific situations where the staggering of J and the current deposition algorithm
    // are not trivial, this check might be too relaxed and we might include a particle that should
    // deposit part of its current in a neighboring box. However, this should catch particles
    // traveling many cells away, for example with algorithms that allow for large time steps.

#if   defined(WARPX_DIM_1D_Z)
    const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(WarpX::noz/2));
#elif   defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(WarpX::nox/2),
                                                       static_cast<int>(WarpX::noz/2));
#elif defined(WARPX_DIM_3D)
    const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(WarpX::nox/2),
                                                       static_cast<int>(WarpX::noy/2),
                                                       static_cast<int>(WarpX::noz/2));
#endif

    // On CPU: particles deposit on tile arrays, which have a small number of guard cells ng_J
    // On GPU: particles deposit directly on the J arrays, which usually have a larger number of guard cells
#ifndef AMREX_USE_GPU
    const amrex::IntVect range = ng_J - shape_extent;
#else
    // Jx, Jy and Jz have the same number of guard cells, hence it is sufficient to check for Jx
    const amrex::IntVect range = jx->nGrowVect() - shape_extent;
#endif
    amrex::ignore_unused(range); // for release builds
    AMREX_ASSERT_WITH_MESSAGE(
        amrex::numParticlesOutOfRange(pti, range) == 0,
        "Particles shape does not fit within tile (CPU) or guard cells (GPU) used for current deposition");

    const std::array<Real,3>& dx = WarpX::CellSize(std::max(depos_lev,0));
    const amrex::ParticleReal q = this->charge;

    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCurrent::Sorting", blp_sort);
    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCurrent::FindMaxTilesize",
            blp_get_max_tilesize);
    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCurrent::DirectCurrentDepKernel",
            direct_current_dep_kernel);
    WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCurrent::EsirkepovCurrentDepKernel",
            esirkepov_current_dep_kernel);
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
    const std::array<amrex::Real, 3>& xyzmin = WarpX::LowerCorner(tilebox, depos_lev, 0.5_rt*dt);

    if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Esirkepov) {
        if (WarpX::grid_type == GridType::Collocated) {
          WARPX_ABORT_WITH_MESSAGE("The Esirkepov algorithm cannot be used with a collocated grid.");
        }
    }

    WARPX_PROFILE_VAR_START(blp_deposit);
    amrex::LayoutData<amrex::Real> * const costs = WarpX::getCosts(lev);
    amrex::Real * const cost = costs ? &((*costs)[pti.index()]) : nullptr;

    // If doing shared mem current deposition, get tile info
    if (WarpX::do_shared_mem_current_deposition) {
        const Geometry& geom = Geom(lev);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();
        const auto domain = geom.Domain();

        Box box = pti.validbox();
        box.grow(ng_J);
        const amrex::IntVect bin_size = WarpX::shared_tilesize;

        //sort particles by bin
        WARPX_PROFILE_VAR_START(blp_sort);
        amrex::DenseBins<ParticleType> bins;
        {
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos = ptile.GetArrayOfStructs();
            auto pstruct_ptr = aos().dataPtr();

            const int ntiles = numTilesInBox(box, true, bin_size);

            bins.build(ptile.numParticles(), pstruct_ptr, ntiles,
                    [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) -> unsigned int
                    {
                        Box tbox;
                        auto iv = getParticleCell(p, plo, dxi, domain);
                        AMREX_ASSERT(box.contains(iv));
                        auto tid = getTileIndex(iv, box, true, bin_size, tbox);
                        return static_cast<unsigned int>(tid);
                    });
        }
        WARPX_PROFILE_VAR_STOP(blp_sort);
        WARPX_PROFILE_VAR_START(blp_get_max_tilesize);
            //get the maximum size necessary for shared mem
            // get tile boxes
        //get the maximum size necessary for shared mem
#if AMREX_SPACEDIM > 0
        const int sizeX = getMaxTboxAlongDim(box.size()[0], WarpX::shared_tilesize[0]);
#endif
#if AMREX_SPACEDIM > 1
        const int sizeZ = getMaxTboxAlongDim(box.size()[1], WarpX::shared_tilesize[1]);
#endif
#if AMREX_SPACEDIM > 2
        const int sizeY = getMaxTboxAlongDim(box.size()[2], WarpX::shared_tilesize[2]);
#endif
        const amrex::IntVect max_tbox_size( AMREX_D_DECL(sizeX,sizeZ,sizeY) );
        WARPX_PROFILE_VAR_STOP(blp_get_max_tilesize);

        // Now pick current deposition algorithm
        if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Esirkepov) {
            WARPX_ABORT_WITH_MESSAGE("Cannot do shared memory deposition with Esirkepov algorithm");
        }
        else if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Vay) {
            WARPX_ABORT_WITH_MESSAGE("Cannot do shared memory deposition with Vay algorithm");
        }
        else {
            WARPX_PROFILE_VAR_START(direct_current_dep_kernel);
            if        (WarpX::nox == 1){
                doDepositionSharedShapeN<1>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, relative_time, dx,
                    xyzmin, lo, q, WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo, bins, box, geom, max_tbox_size);
            } else if (WarpX::nox == 2){
                doDepositionSharedShapeN<2>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, relative_time, dx,
                    xyzmin, lo, q, WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo, bins, box, geom, max_tbox_size);
            } else if (WarpX::nox == 3){
                doDepositionSharedShapeN<3>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, relative_time, dx,
                    xyzmin, lo, q, WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo, bins, box, geom, max_tbox_size);
            }
            WARPX_PROFILE_VAR_STOP(direct_current_dep_kernel);
        }
    }
    // If not doing shared memory deposition, call normal kernels
    else {
        if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Esirkepov) {
            if        (WarpX::nox == 1){
                doEsirkepovDepositionShapeN<1>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_arr, jy_arr, jz_arr, np_to_depose, dt, relative_time, dx, xyzmin, lo, q,
                    WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            } else if (WarpX::nox == 2){
                doEsirkepovDepositionShapeN<2>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_arr, jy_arr, jz_arr, np_to_depose, dt, relative_time, dx, xyzmin, lo, q,
                    WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            } else if (WarpX::nox == 3){
                doEsirkepovDepositionShapeN<3>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_arr, jy_arr, jz_arr, np_to_depose, dt, relative_time, dx, xyzmin, lo, q,
                    WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            }
        } else if (WarpX::current_deposition_algo == CurrentDepositionAlgo::Vay) {
            if        (WarpX::nox == 1){
                doVayDepositionShapeN<1>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, dt, relative_time, dx, xyzmin, lo, q,
                    WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            } else if (WarpX::nox == 2){
                doVayDepositionShapeN<2>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, dt, relative_time, dx, xyzmin, lo, q,
                    WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            } else if (WarpX::nox == 3){
                doVayDepositionShapeN<3>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, dt, relative_time, dx, xyzmin, lo, q,
                    WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            }
        } else { // Direct deposition
            if        (WarpX::nox == 1){
                doDepositionShapeN<1>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, relative_time, dx,
                    xyzmin, lo, q, WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            } else if (WarpX::nox == 2){
                doDepositionShapeN<2>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, relative_time, dx,
                    xyzmin, lo, q, WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            } else if (WarpX::nox == 3){
                doDepositionShapeN<3>(
                    GetPosition, wp.dataPtr() + offset, uxp.dataPtr() + offset,
                    uyp.dataPtr() + offset, uzp.dataPtr() + offset, ion_lev,
                    jx_fab, jy_fab, jz_fab, np_to_depose, relative_time, dx,
                    xyzmin, lo, q, WarpX::n_rz_azimuthal_modes, cost,
                    WarpX::load_balance_costs_update_algo);
            }
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

void
WarpXParticleContainer::DepositCurrent (
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& J,
    const amrex::Real dt, const amrex::Real relative_time)
{
    // Loop over the refinement levels
    int const finest_level = J.size() - 1;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Loop over particle tiles and deposit current on each level
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
        {
        const int thread_num = omp_get_thread_num();
#else
        const int thread_num = 0;
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const long np = pti.numParticles();
            const auto & wp = pti.GetAttribs(PIdx::w);
            const auto & uxp = pti.GetAttribs(PIdx::ux);
            const auto & uyp = pti.GetAttribs(PIdx::uy);
            const auto & uzp = pti.GetAttribs(PIdx::uz);

            int* AMREX_RESTRICT ion_lev = nullptr;
            if (do_field_ionization)
            {
                ion_lev = pti.GetiAttribs(particle_icomps["ionizationLevel"]).dataPtr();
            }

            DepositCurrent(pti, wp, uxp, uyp, uzp, ion_lev,
                           J[lev][0].get(), J[lev][1].get(), J[lev][2].get(),
                           0, np, thread_num, lev, lev, dt, relative_time);
        }
#ifdef AMREX_USE_OMP
        }
#endif
    }
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
                        Particles [offset,offset+np_tp_depose) deposit charge
 * \param thread_num  : Thread number (if tiling)
 * \param lev         : Level of box that contains particles
 * \param depos_lev   : Level on which particles deposit (if buffers are used)
 */
void
WarpXParticleContainer::DepositCharge (WarpXParIter& pti, RealVector const& wp,
                                       const int * const ion_lev,
                                       amrex::MultiFab* rho, const int icomp,
                                       const long offset, const long np_to_depose,
                                       const int thread_num, const int lev, const int depos_lev)
{
    if (WarpX::do_shared_mem_charge_deposition)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE((depos_lev==(lev-1)) ||
                                         (depos_lev==(lev  )),
                                         "Deposition buffers only work for lev-1");

        // If no particles, do not do anything
        if (np_to_depose == 0) return;

        // Number of guard cells for local deposition of rho
        const WarpX& warpx = WarpX::GetInstance();
        const amrex::IntVect& ng_rho = warpx.get_ng_depos_rho();

        // Extract deposition order and check that particles shape fits within the guard cells.
        // NOTE: In specific situations where the staggering of rho and the charge deposition algorithm
        // are not trivial, this check might be too strict and we might need to relax it, as currently
        // done for the current deposition.

#if   defined(WARPX_DIM_1D_Z)
        const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(WarpX::noz/2));
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
        const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(WarpX::nox/2+1),
                                                           static_cast<int>(WarpX::noz/2+1));
#elif defined(WARPX_DIM_3D)
        const amrex::IntVect shape_extent = amrex::IntVect(static_cast<int>(WarpX::nox/2+1),
                                                           static_cast<int>(WarpX::noy/2+1),
                                                           static_cast<int>(WarpX::noz/2+1));
#endif

        // On CPU: particles deposit on tile arrays, which have a small number of guard cells ng_rho
        // On GPU: particles deposit directly on the rho array, which usually have a larger number of guard cells
#ifndef AMREX_USE_GPU
        const amrex::IntVect range = ng_rho - shape_extent;
#else
        const amrex::IntVect range = rho->nGrowVect() - shape_extent;
#endif

        AMREX_ASSERT_WITH_MESSAGE(
                                  amrex::numParticlesOutOfRange(pti, range) == 0,
                                  "Particles shape does not fit within tile (CPU) or guard cells (GPU) used for charge deposition");
        amrex::ignore_unused(range); // In case the assertion isn't compiled

        const std::array<Real,3>& dx = WarpX::CellSize(std::max(depos_lev,0));
        const Real q = this->charge;

        WARPX_PROFILE_VAR_NS("WarpXParticleContainer::DepositCharge::Sorting", blp_sort);
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

        const auto ix_type = rho->ixType().toIntVect();
#ifndef AMREX_USE_GPU
        // Staggered tile box
        Box tb = amrex::convert( tilebox, ix_type );
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

        // Lower corner of tile box physical domain
        // Note that this includes guard cells since it is after tilebox.ngrow
        // Take into account Galilean shift
        const Real dt = warpx.getdt(lev);
        const amrex::Real time_shift_delta = (icomp == 0 ? 0.0_rt : dt);
        const std::array<amrex::Real,3>& xyzmin = WarpX::LowerCorner(
                tilebox, depos_lev, time_shift_delta);

        // Indices of the lower bound
        const Dim3 lo = lbound(tilebox);

        // HACK - sort particles by bin here.
        WARPX_PROFILE_VAR_START(blp_sort);
        amrex::DenseBins<ParticleType> bins;
        {
            const Geometry& geom = Geom(lev);
            const auto dxi = geom.InvCellSizeArray();
            const auto plo = geom.ProbLoArray();
            const auto domain = geom.Domain();

            auto& ptile = ParticlesAt(lev, pti);
            auto& aos   = ptile.GetArrayOfStructs();
            auto pstruct_ptr = aos().dataPtr();

            Box box = pti.validbox();
            box.grow(ng_rho);
            const amrex::IntVect bin_size = WarpX::shared_tilesize;
            const int ntiles = numTilesInBox(box, true, bin_size);

            bins.build(ptile.numParticles(), pstruct_ptr, ntiles,
                       [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) -> unsigned int
                       {
                           Box tbx;
                           auto iv = getParticleCell(p, plo, dxi, domain);
                           AMREX_ASSERT(box.contains(iv));
                           auto tid = getTileIndex(iv, box, true, bin_size, tbx);
                           return static_cast<unsigned int>(tid);
                       });
        }
        WARPX_PROFILE_VAR_STOP(blp_sort);

        // get tile boxes
        amrex::Gpu::DeviceVector<Box> tboxes(bins.numBins(), amrex::Box());
        {
            const Geometry& geom = Geom(lev);
            const auto dxi = geom.InvCellSizeArray();
            const auto plo = geom.ProbLoArray();
            const auto domain = geom.Domain();

            auto& ptile = ParticlesAt(lev, pti);
            auto& aos   = ptile.GetArrayOfStructs();
            auto pstruct_ptr = aos().dataPtr();

            Box box = pti.validbox();
            box.grow(ng_rho);
            const amrex::IntVect bin_size = WarpX::shared_tilesize;

            const auto offsets_ptr = bins.offsetsPtr();
            auto tbox_ptr = tboxes.dataPtr();
            auto permutation = bins.permutationPtr();
            amrex::ParallelFor(bins.numBins(),
                               [=] AMREX_GPU_DEVICE (int ibin) {
                                   const int bin_start = offsets_ptr[ibin];
                                   const int bin_stop = offsets_ptr[ibin+1];
                                   if (bin_start < bin_stop) {
                                       auto p = pstruct_ptr[permutation[bin_start]];
                                       Box tbx;
                                       auto iv = getParticleCell(p, plo, dxi, domain);
                                       AMREX_ASSERT(box.contains(iv));
                                       [[maybe_unused]] auto tid = getTileIndex(iv, box, true, bin_size, tbx);
                                       AMREX_ASSERT(tid == ibin);
                                       AMREX_ASSERT(tbx.contains(iv));
                                       tbox_ptr[ibin] = tbx;
                                   }
                               });
        }

        // compute max tile box size in each direction
        amrex::IntVect max_tbox_size;
        {
            ReduceOps<AMREX_D_DECL(ReduceOpMax, ReduceOpMax, ReduceOpMax)> reduce_op;
            ReduceData<AMREX_D_DECL(int, int, int)> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            const auto boxes_ptr = tboxes.dataPtr();
            reduce_op.eval(tboxes.size(), reduce_data,
                           [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                           {
                               const Box& box = boxes_ptr[i];
                               if (box.ok()) {
                                   IntVect si = box.length();
                                   return {AMREX_D_DECL(si[0], si[1], si[2])};
                               } else {
                                   return {AMREX_D_DECL(0, 0, 0)};
                               }
                           });

            ReduceTuple hv = reduce_data.value();

            max_tbox_size = IntVect(AMREX_D_DECL(amrex::get<0>(hv),
                                                 amrex::get<1>(hv),
                                                 amrex::get<2>(hv)));
        }

        WARPX_PROFILE_VAR_START(blp_ppc_chd);
        amrex::LayoutData<amrex::Real>* costs = WarpX::getCosts(lev);
        amrex::Real* cost = costs ? &((*costs)[pti.index()]) : nullptr;

        const auto GetPosition = GetParticlePosition(pti, offset);
        const Geometry& geom = Geom(lev);
        Box box = pti.validbox();
        box.grow(ng_rho);

        if (WarpX::nox == 1){
            doChargeDepositionSharedShapeN<1>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                              rho_fab, ix_type, np_to_depose, dx, xyzmin, lo, q,
                                              WarpX::n_rz_azimuthal_modes, cost,
                                              WarpX::load_balance_costs_update_algo, bins, box, geom, max_tbox_size);
        } else if (WarpX::nox == 2){
            doChargeDepositionSharedShapeN<2>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                              rho_fab, ix_type, np_to_depose, dx, xyzmin, lo, q,
                                              WarpX::n_rz_azimuthal_modes, cost,
                                              WarpX::load_balance_costs_update_algo, bins, box, geom, max_tbox_size);
        } else if (WarpX::nox == 3){
            doChargeDepositionSharedShapeN<3>(GetPosition, wp.dataPtr()+offset, ion_lev,
                                              rho_fab, ix_type, np_to_depose, dx, xyzmin, lo, q,
                                              WarpX::n_rz_azimuthal_modes, cost,
                                              WarpX::load_balance_costs_update_algo, bins, box, geom, max_tbox_size);
        }
#ifndef AMREX_USE_GPU
        // CPU, tiling: atomicAdd local_rho into rho
        WARPX_PROFILE_VAR_START(blp_accumulate);
        (*rho)[pti].atomicAdd(local_rho[thread_num], tb, tb, 0, icomp*nc, nc);
        WARPX_PROFILE_VAR_STOP(blp_accumulate);
#endif
    } else {

        const WarpX& warpx = WarpX::GetInstance();

        // deposition guards
        //   note: this is smaller than rho->nGrowVect() for PSATD
        const amrex::IntVect& ng_rho = warpx.get_ng_depos_rho();

        const std::array<amrex::Real,3>& dx = WarpX::CellSize(std::max(depos_lev,0));
        amrex::IntVect ref_ratio;
        if (lev == depos_lev) {
            ref_ratio = IntVect(AMREX_D_DECL(1, 1, 1 ));
        } else {
            ref_ratio = WarpX::RefRatio(depos_lev);
        }
        const int nc = WarpX::ncomps;

        // Get tile box where charge is deposited.
        // The tile box is different when depositing in the buffers (depos_lev<lev)
        // or when depositing inside the level (depos_lev=lev)
        amrex::Box tilebox;
        if (lev == depos_lev) {
            tilebox = pti.tilebox();
        } else {
            tilebox = amrex::coarsen(pti.tilebox(), ref_ratio);
        }
        tilebox.grow(ng_rho);

        // Lower corner of tile box physical domain
        // Note that this includes guard cells since it is after tilebox.ngrow
        // Take into account Galilean shift
        const amrex::Real dt = warpx.getdt(lev);
        const amrex::Real time_shift_delta = (icomp == 0 ? 0.0_rt : dt);
        const std::array<amrex::Real,3>& xyzmin = WarpX::LowerCorner(tilebox, depos_lev, time_shift_delta);

        // pointer to costs data
        amrex::LayoutData<amrex::Real>* costs = WarpX::getCosts(lev);
        amrex::Real* cost = costs ? &((*costs)[pti.index()]) : nullptr;

        AMREX_ALWAYS_ASSERT(WarpX::nox == WarpX::noy);
        AMREX_ALWAYS_ASSERT(WarpX::nox == WarpX::noz);

        ablastr::particles::deposit_charge<WarpXParticleContainer>(
            pti, wp, this->charge, ion_lev,
            rho, local_rho[thread_num],
            WarpX::noz, dx, xyzmin, WarpX::n_rz_azimuthal_modes,
            ng_rho, depos_lev, ref_ratio,
            offset, np_to_depose,
            icomp, nc,
            cost, WarpX::load_balance_costs_update_algo, WarpX::do_device_synchronize
        );
    }
}

void
WarpXParticleContainer::DepositCharge (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                                       const bool local, const bool reset,
                                       const bool apply_boundary_and_scale_volume,
                                       const bool interpolate_across_levels,
                                       const int icomp)
{
    WARPX_PROFILE("WarpXParticleContainer::DepositCharge");

    // Loop over the refinement levels
    int const finest_level = rho.size() - 1;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        // Reset the rho array if reset is True
        int const nc = WarpX::ncomps;
        if (reset) rho[lev]->setVal(0., icomp*nc, nc, rho[lev]->nGrowVect());

        // Loop over particle tiles and deposit charge on each level
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
        {
        const int thread_num = omp_get_thread_num();
#else
        const int thread_num = 0;
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const long np = pti.numParticles();
            auto const & wp = pti.GetAttribs(PIdx::w);

            int* AMREX_RESTRICT ion_lev = nullptr;
            if (do_field_ionization)
            {
                ion_lev = pti.GetiAttribs(particle_icomps["ionizationLevel"]).dataPtr();
            }

            DepositCharge(pti, wp, ion_lev, rho[lev].get(), icomp, 0, np, thread_num, lev, lev);
        }
#ifdef AMREX_USE_OMP
        }
#endif

#ifdef WARPX_DIM_RZ
        if (apply_boundary_and_scale_volume)
        {
            WarpX::GetInstance().ApplyInverseVolumeScalingToChargeDensity(rho[lev].get(), lev);
        }
#endif

        // Exchange guard cells
        if (local == false) {
            ablastr::utils::communication::SumBoundary(
                *rho[lev],
                WarpX::do_single_precision_comms,
                m_gdb->Geom(lev).periodicity()
            );
        }

#ifndef WARPX_DIM_RZ
        if (apply_boundary_and_scale_volume)
        {
            // Reflect density over PEC boundaries, if needed.
            WarpX::GetInstance().ApplyRhofieldBoundary(lev, rho[lev].get(), PatchType::fine);
        }
#endif
    }

    // Now that the charge has been deposited at each level,
    // we average down from fine to crse
    if (interpolate_across_levels)
    {
        for (int lev = finest_level - 1; lev >= 0; --lev) {
            const DistributionMapping& fine_dm = rho[lev+1]->DistributionMap();
            BoxArray coarsened_fine_BA = rho[lev+1]->boxArray();
            coarsened_fine_BA.coarsen(m_gdb->refRatio(lev));
            const IntVect ngrow = (rho[lev+1]->nGrowVect()+1)/m_gdb->refRatio(lev);
            MultiFab coarsened_fine_data(coarsened_fine_BA, fine_dm, rho[lev+1]->nComp(), ngrow );
            coarsened_fine_data.setVal(0.0);

            ablastr::coarsen::average::Coarsen(coarsened_fine_data, *rho[lev + 1], m_gdb->refRatio(lev) );
            ablastr::utils::communication::ParallelAdd(*rho[lev], coarsened_fine_data, 0, 0,
                                                       rho[lev]->nComp(),
                                                       amrex::IntVect::TheZeroVector(),
                                                       amrex::IntVect::TheZeroVector(),
                                                       WarpX::do_single_precision_comms,
                                                       m_gdb->Geom(lev).periodicity());
        }
    }
}

std::unique_ptr<MultiFab>
WarpXParticleContainer::GetChargeDensity (int lev, bool local)
{
    const auto& gm = m_gdb->Geom(lev);
    const auto& ba = m_gdb->ParticleBoxArray(lev);
    const auto& dm = m_gdb->DistributionMap(lev);
    BoxArray nba = ba;

#ifdef WARPX_DIM_RZ
    const bool is_PSATD_RZ =
        (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD);
#else
    const bool is_PSATD_RZ = false;
#endif
    if( !is_PSATD_RZ )
        nba.surroundingNodes();

    // Number of guard cells for local deposition of rho
    const WarpX& warpx = WarpX::GetInstance();
    const int ng_rho = warpx.get_ng_depos_rho().max();

    auto rho = std::make_unique<MultiFab>(nba,dm,WarpX::ncomps,ng_rho);
    rho->setVal(0.0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
    {
#endif
#ifdef AMREX_USE_OMP
        const int thread_num = omp_get_thread_num();
#else
        const int thread_num = 0;
#endif

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const long np = pti.numParticles();
            auto& wp = pti.GetAttribs(PIdx::w);

            int* AMREX_RESTRICT ion_lev;
            if (do_field_ionization){
                ion_lev = pti.GetiAttribs(particle_icomps["ionizationLevel"]).dataPtr();
            } else {
                ion_lev = nullptr;
            }

            DepositCharge(pti, wp, ion_lev, rho.get(), 0, 0, np,
                          thread_num, lev, lev);
        }
#ifdef AMREX_USE_OMP
    }
#endif

#ifdef WARPX_DIM_RZ
    WarpX::GetInstance().ApplyInverseVolumeScalingToChargeDensity(rho.get(), lev);
#endif

    if (local == false) { ablastr::utils::communication::SumBoundary(*rho, WarpX::do_single_precision_comms, gm.periodicity()); }

#ifndef WARPX_DIM_RZ
    // Reflect density over PEC boundaries, if needed.
    WarpX::GetInstance().ApplyRhofieldBoundary(lev, rho.get(), PatchType::fine);
#endif

    return rho;
}

amrex::ParticleReal WarpXParticleContainer::sumParticleCharge(bool local) {

    amrex::ParticleReal total_charge = 0.0;

    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev)
    {

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(+:total_charge)
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto& wp = pti.GetAttribs(PIdx::w);
            for (const auto& ww : wp) {
                total_charge += ww;
            }
        }
    }

    if (local == false) ParallelDescriptor::ReduceRealSum(total_charge);
    total_charge *= this->charge;
    return total_charge;
}

std::array<ParticleReal, 3> WarpXParticleContainer::meanParticleVelocity(bool local) {

    amrex::ParticleReal vx_total = 0.0_prt;
    amrex::ParticleReal vy_total = 0.0_prt;
    amrex::ParticleReal vz_total = 0.0_prt;

    amrex::Long np_total = 0;

    const amrex::ParticleReal inv_clight_sq = 1.0_prt/PhysConst::c/PhysConst::c;

    const int nLevels = finestLevel();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<ParticleReal, ParticleReal, ParticleReal> reduce_data(reduce_op);
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
                                   amrex::ParticleReal usq = (uxp[i]*uxp[i] +
                                                              uyp[i]*uyp[i] +
                                                              uzp[i]*uzp[i])*inv_clight_sq;
                                   const amrex::ParticleReal gaminv = 1.0_prt/std::sqrt(1.0_prt + usq);
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
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(+:vx_total, vy_total, vz_total, np_total)
#endif
            for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
            {
                auto& ux = pti.GetAttribs(PIdx::ux);
                auto& uy = pti.GetAttribs(PIdx::uy);
                auto& uz = pti.GetAttribs(PIdx::uz);

                np_total += pti.numParticles();

                for (unsigned long i = 0; i < ux.size(); i++) {
                    const amrex::ParticleReal usq = (ux[i]*ux[i] + uy[i]*uy[i] + uz[i]*uz[i])*inv_clight_sq;
                    const amrex::ParticleReal gaminv = 1.0_prt/std::sqrt(1.0_prt + usq);
                    vx_total += ux[i]*gaminv;
                    vy_total += uy[i]*gaminv;
                    vz_total += uz[i]*gaminv;
                }
            }
        }
    }

    if (local == false) {
        ParallelDescriptor::ReduceRealSum(vx_total);
        ParallelDescriptor::ReduceRealSum(vy_total);
        ParallelDescriptor::ReduceRealSum(vz_total);
        ParallelDescriptor::ReduceLongSum(np_total);
    }

    std::array<amrex::ParticleReal, 3> mean_v;
    if (np_total > 0) {
        mean_v[0] = vx_total / np_total;
        mean_v[1] = vy_total / np_total;
        mean_v[2] = vz_total / np_total;
    }

    return mean_v;
}

amrex::ParticleReal WarpXParticleContainer::maxParticleVelocity(bool local) {

    amrex::ParticleReal max_v = 0.0;

    const int nLevels = finestLevel();
    for (int lev = 0; lev <= nLevels; ++lev)
    {

#ifdef AMREX_USE_OMP
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

    if (local == false) ParallelAllReduce::Max(max_v, ParallelDescriptor::Communicator());
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

    amrex::LayoutData<amrex::Real>* costs = WarpX::getCosts(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {

        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            if (costs && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
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

            if (costs && WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
            {
                amrex::Gpu::synchronize();
                wt = amrex::second() - wt;
                amrex::HostDevice::Atomic::Add( &(*costs)[pti.index()], wt);
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
WarpXParticleContainer::ApplyBoundaryConditions (){
    WARPX_PROFILE("WarpXParticleContainer::ApplyBoundaryConditions()");

    // Periodic boundaries are handled in AMReX code
    if (m_boundary_conditions.CheckAll(ParticleBoundaryType::Periodic)) return;

    auto boundary_conditions = m_boundary_conditions.data;

    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (WarpXParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            auto GetPosition = GetParticlePosition(pti);
            auto SetPosition = SetParticlePosition(pti);
#ifndef WARPX_DIM_1D_Z
            const Real xmin = Geom(lev).ProbLo(0);
            const Real xmax = Geom(lev).ProbHi(0);
#endif
#ifdef WARPX_DIM_3D
            const Real ymin = Geom(lev).ProbLo(1);
            const Real ymax = Geom(lev).ProbHi(1);
#endif
            const Real zmin = Geom(lev).ProbLo(WARPX_ZINDEX);
            const Real zmax = Geom(lev).ProbHi(WARPX_ZINDEX);

            ParticleTileType& ptile = ParticlesAt(lev, pti);
            ParticleType * const pp = ptile.GetArrayOfStructs()().data();

            auto& soa = ptile.GetStructOfArrays();
            amrex::ParticleReal * const AMREX_RESTRICT ux = soa.GetRealData(PIdx::ux).data();
            amrex::ParticleReal * const AMREX_RESTRICT uy = soa.GetRealData(PIdx::uy).data();
            amrex::ParticleReal * const AMREX_RESTRICT uz = soa.GetRealData(PIdx::uz).data();

            // Loop over particles and apply BC to each particle
            amrex::ParallelForRNG(
                pti.numParticles(),
                [=] AMREX_GPU_DEVICE (long i, amrex::RandomEngine const& engine) {
                    ParticleType& p = pp[i];

                    // skip particles that are already flagged for removal
                    if (p.id() < 0) return;

                    ParticleReal x, y, z;
                    GetPosition.AsStored(i, x, y, z);
                    // Note that for RZ, (x, y, z) is actually (r, theta, z).

                    bool particle_lost = false;
                    ApplyParticleBoundaries::apply_boundaries(
#ifndef WARPX_DIM_1D_Z
                                                              x, xmin, xmax,
#endif
#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
                                                              y,
#endif
#if (defined WARPX_DIM_3D)
                                                              ymin, ymax,
#endif
                                                              z, zmin, zmax,
                                                              ux[i], uy[i], uz[i], particle_lost,
                                                              boundary_conditions, engine);

                    if (particle_lost) {
                        p.id() = -p.id();
                    } else {
                        SetPosition.AsStored(i, x, y, z);
                    }
                }
            );
        }
    }
}
