/* Copyright 2019 David Grote, Luca Fedeli, Maxence Thevenet
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PhotonParticleContainer.H"

#ifdef WARPX_QED
#   include "Particles/ElementaryProcess/QEDInternals/BreitWheelerEngineWrapper.H"
#endif
#include "Particles/Gather/FieldGather.H"
#include "Particles/Gather/GetExternalFields.H"
#include "Particles/PhysicalParticleContainer.H"
#include "Particles/Pusher/CopyParticleAttribs.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/Pusher/UpdatePositionPhoton.H"
#include "Particles/WarpXParticleContainer.H"
#include "WarpX.H"

#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Dim3.H>
#include <AMReX_Extension.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_StructOfArrays.H>

#include <algorithm>
#include <array>
#include <map>
#include <memory>

using namespace amrex;

PhotonParticleContainer::PhotonParticleContainer (AmrCore* amr_core, int ispecies,
                                                  const std::string& name)
    : PhysicalParticleContainer(amr_core, ispecies, name)
{
    ParmParse pp_species_name(species_name);

#ifdef WARPX_QED
        //Find out if Breit Wheeler process is enabled
        pp_species_name.query("do_qed_breit_wheeler", m_do_qed_breit_wheeler);

        //If Breit Wheeler process is enabled, look for the target electron and positron
        //species
        if(m_do_qed_breit_wheeler){
            pp_species_name.get("qed_breit_wheeler_ele_product_species", m_qed_breit_wheeler_ele_product_name);
            pp_species_name.get("qed_breit_wheeler_pos_product_species", m_qed_breit_wheeler_pos_product_name);
        }

        //Check for processes which do not make sense for photons
        bool test_quantum_sync = false;
        pp_species_name.query("do_qed_quantum_sync", test_quantum_sync);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        test_quantum_sync == 0,
        "ERROR: do_qed_quantum_sync can be 1 for species NOT listed in particles.photon_species only!");
        //_________________________________________________________
#endif

}

void PhotonParticleContainer::InitData()
{
    AddParticles(0); // Note - add on level 0

    Redistribute();  // We then redistribute

}

void
PhotonParticleContainer::PushPX (WarpXParIter& pti,
                                 amrex::FArrayBox const * exfab,
                                 amrex::FArrayBox const * eyfab,
                                 amrex::FArrayBox const * ezfab,
                                 amrex::FArrayBox const * bxfab,
                                 amrex::FArrayBox const * byfab,
                                 amrex::FArrayBox const * bzfab,
                                 const amrex::IntVect ngEB, const int /*e_is_nodal*/,
                                 const long offset,
                                 const long np_to_push,
                                 int lev, int gather_lev,
                                 amrex::Real dt, ScaleFields /*scaleFields*/, DtType a_dt_type)
{
    // Get cell size on gather_lev
    const std::array<Real,3>& dx = WarpX::CellSize(std::max(gather_lev,0));

    // Get box from which field is gathered.
    // If not gathering from the finest level, the box is coarsened.
    amrex::Box box;
    if (lev == gather_lev) {
        box = pti.tilebox();
    } else {
        const IntVect& ref_ratio = WarpX::RefRatio(gather_lev);
        box = amrex::coarsen(pti.tilebox(),ref_ratio);
    }

    // Add guard cells to the box.
    box.grow(ngEB);

    auto& attribs = pti.GetAttribs();

    // Extract pointers to the different particle quantities
    ParticleReal* const AMREX_RESTRICT ux = attribs[PIdx::ux].dataPtr() + offset;
    ParticleReal* const AMREX_RESTRICT uy = attribs[PIdx::uy].dataPtr() + offset;
    ParticleReal* const AMREX_RESTRICT uz = attribs[PIdx::uz].dataPtr() + offset;

#ifdef WARPX_QED
    BreitWheelerEvolveOpticalDepth evolve_opt;
    amrex::ParticleReal* AMREX_RESTRICT p_optical_depth_BW = nullptr;
    const bool local_has_breit_wheeler = has_breit_wheeler();
    if (local_has_breit_wheeler) {
        evolve_opt = m_shr_p_bw_engine->build_evolve_functor();
        p_optical_depth_BW = pti.GetAttribs(particle_comps["opticalDepthBW"]).dataPtr() + offset;
    }
#endif

    auto copyAttribs = CopyParticleAttribs(pti, tmp_particle_data, offset);
    int do_copy = (WarpX::do_back_transformed_diagnostics &&
                   do_back_transformed_diagnostics && a_dt_type!=DtType::SecondHalf);

    const auto GetPosition = GetParticlePosition(pti, offset);
    auto SetPosition = SetParticlePosition(pti, offset);

    const auto getExternalEB = GetExternalEBField(pti, offset);

    // Lower corner of tile box physical domain (take into account Galilean shift)
    amrex::Real cur_time = WarpX::GetInstance().gett_new(lev);
    const auto& time_of_last_gal_shift = WarpX::GetInstance().time_of_last_gal_shift;
    amrex::Real time_shift = (cur_time - time_of_last_gal_shift);
    amrex::Array<amrex::Real,3> galilean_shift = {
        m_v_galilean[0]*time_shift,
        m_v_galilean[1]*time_shift,
        m_v_galilean[2]*time_shift };
    const std::array<Real, 3>& xyzmin = WarpX::LowerCorner(box, galilean_shift, gather_lev);

    const Dim3 lo = lbound(box);

    bool galerkin_interpolation = WarpX::galerkin_interpolation;
    int nox = WarpX::nox;
    int n_rz_azimuthal_modes = WarpX::n_rz_azimuthal_modes;

    amrex::GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
    amrex::GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};

    amrex::Array4<const amrex::Real> const& ex_arr = exfab->array();
    amrex::Array4<const amrex::Real> const& ey_arr = eyfab->array();
    amrex::Array4<const amrex::Real> const& ez_arr = ezfab->array();
    amrex::Array4<const amrex::Real> const& bx_arr = bxfab->array();
    amrex::Array4<const amrex::Real> const& by_arr = byfab->array();
    amrex::Array4<const amrex::Real> const& bz_arr = bzfab->array();

    amrex::IndexType const ex_type = exfab->box().ixType();
    amrex::IndexType const ey_type = eyfab->box().ixType();
    amrex::IndexType const ez_type = ezfab->box().ixType();
    amrex::IndexType const bx_type = bxfab->box().ixType();
    amrex::IndexType const by_type = byfab->box().ixType();
    amrex::IndexType const bz_type = bzfab->box().ixType();

    const auto t_do_not_gather = do_not_gather;

    amrex::ParallelFor(
        np_to_push,
        [=] AMREX_GPU_DEVICE (long i) {
            if (do_copy) copyAttribs(i);
            ParticleReal x, y, z;
            GetPosition(i, x, y, z);

            amrex::ParticleReal Exp=0, Eyp=0, Ezp=0;
            amrex::ParticleReal Bxp=0, Byp=0, Bzp=0;

            if(!t_do_not_gather){
                // first gather E and B to the particle positions
                doGatherShapeN(x, y, z, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                               ex_arr, ey_arr, ez_arr, bx_arr, by_arr, bz_arr,
                               ex_type, ey_type, ez_type, bx_type, by_type, bz_type,
                               dx_arr, xyzmin_arr, lo, n_rz_azimuthal_modes,
                               nox, galerkin_interpolation);
            }
            getExternalEB(i, Exp, Eyp, Ezp, Bxp, Byp, Bzp);

#ifdef WARPX_QED
            if (local_has_breit_wheeler) {
                evolve_opt(ux[i], uy[i], uz[i], Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                    dt, p_optical_depth_BW[i]);
            }
#endif

            UpdatePositionPhoton( x, y, z, ux[i], uy[i], uz[i], dt );
            SetPosition(i, x, y, z);
        }
    );
}

void
PhotonParticleContainer::Evolve (int lev,
                                 const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                                 const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                                 MultiFab& jx, MultiFab& jy, MultiFab& jz,
                                 MultiFab* cjx, MultiFab* cjy, MultiFab* cjz,
                                 MultiFab* rho, MultiFab* crho,
                                 const MultiFab* cEx, const MultiFab* cEy, const MultiFab* cEz,
                                 const MultiFab* cBx, const MultiFab* cBy, const MultiFab* cBz,
                                 Real t, Real dt, DtType a_dt_type, bool skip_deposition)
{
    // This does gather, push and depose.
    // Push and depose have been re-written for photons
    PhysicalParticleContainer::Evolve (lev,
                                       Ex, Ey, Ez,
                                       Bx, By, Bz,
                                       jx, jy, jz,
                                       cjx, cjy, cjz,
                                       rho, crho,
                                       cEx, cEy, cEz,
                                       cBx, cBy, cBz,
                                       t, dt, a_dt_type, skip_deposition);

}
