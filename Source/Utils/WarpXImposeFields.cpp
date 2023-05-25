/* Copyright 2023 Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>
#include <AMReX_VisMF.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>

using namespace amrex;

namespace
{
    std::unique_ptr<MultiFab> read_snapshot_field (BoxArray& ba,
                                                   DistributionMapping& dm,
                                                   IndexType itype,
                                                   std::string const& filename)
    {
        if (ba.empty()) {
            auto mf = std::make_unique<MultiFab>();
            VisMF::Read(*mf, filename);
            ba = mf->boxArray();
            dm = mf->DistributionMap();
            return mf;
        } else {
            VisMF vismf(filename);
            auto mf = std::make_unique<MultiFab>(amrex::convert(ba,itype),
                                                 dm, vismf.nComp(),
                                                 vismf.nGrowVect());
            VisMF::Read(*mf, filename);
            return mf;
        }
    }

    template <typename ND, typename CC>
    void impose_snapshot_field_direct (MultiFab& mf, MultiFab const& mf_ext,
                                       BoxArray const& baext,
                                       DistributionMapping const& dmext,
                                       Box const& planebox,
                                       std::map<int,int> const& index_map,
                                       ND const& nd_interp_info, CC const& cc_interp_info)
    {
        const int zdir = AMREX_SPACEDIM-1;
        auto typ = mf.ixType();
        bool is_nodal = typ.nodeCentered(zdir);
        MultiFab mf_src(amrex::convert(baext,typ), dmext, 1, 0);
        mf_src.setVal(0.0_rt);
        mf_src.ParallelCopy(mf_ext, 0, 0, 1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            bx.setBig(zdir, std::min(bx.bigEnd(zdir),planebox.smallEnd(zdir)-1));
            if (bx.ok()) {
                mf[mfi].template setVal<RunOn::Device>(0.0_rt, bx);
            }

            bx = mfi.tilebox() & amrex::convert(planebox,typ);
            if (bx.ok()) {
                auto const& dst = mf.array(mfi);
                auto const& src = mf_src.array(index_map.at(mfi.index()));
                AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM > 1);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
#if (AMREX_SPACEDIM == 2)
                    auto [jj,w] = is_nodal ? nd_interp_info(j)
                                           : cc_interp_info(j);
                    dst(i,j,k) = src(i,jj,k)*w + src(i,jj+1,k)*(1.0_rt-w);
#else
                    auto [kk,w] = is_nodal ? nd_interp_info(k)
                                           : cc_interp_info(k);
                    dst(i,j,k) = src(i,j,kk)*w + src(i,j,kk+1)*(1.0_rt-w);
#endif
                });
            }
        }
    }

    // mf = c1*mfext1 + c2*mfext2
    // mfext1 has the same index type as mf.
    // mfext2 has a different index type in the boosted direction.
    template <typename ND, typename CC, typename N2C, typename C2N>
    void impose_snapshot_field_boost (MultiFab& mf,
                                      MultiFab const& mfext1, Real c1,
                                      MultiFab const& mfext2, Real c2,
                                      BoxArray const& baext,
                                      DistributionMapping const& dmext,
                                      Box const& planebox,
                                      std::map<int,int> const& index_map,
                                      ND const& nd_interp_info, CC const& cc_interp_info,
                                      N2C const& nd2cc_interp_info, C2N const& cc2nd_interp_info)
    {
        const int zdir = AMREX_SPACEDIM-1;
        auto typ1 = mfext1.ixType();
        auto typ2 = mfext2.ixType();
        typ2.flip(zdir);
        AMREX_ALWAYS_ASSERT(typ1 == typ2);
        typ2.flip(zdir);
        bool is_nodal = typ1.nodeCentered(zdir);
        MultiFab mfsrc1(amrex::convert(baext,typ1), dmext, 1, 0);
        mfsrc1.setVal(0.0_rt);
        mfsrc1.ParallelCopy(mfext1, 0, 0, 1);
        MultiFab mfsrc2(amrex::convert(baext,typ2), dmext, 1, 0);
        mfsrc2.setVal(0.0_rt);
        mfsrc2.ParallelCopy(mfext2, 0, 0, 1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box bx = mfi.tilebox();
            bx.setBig(zdir, std::min(bx.bigEnd(zdir),planebox.smallEnd(zdir)-1));
            if (bx.ok()) {
                mf[mfi].template setVal<RunOn::Device>(0.0_rt, bx);
            }

            bx = mfi.tilebox() & amrex::convert(planebox,typ1);
            if (bx.ok()) {
                auto const& dst = mf.array(mfi);
                auto const idx = index_map.at(mfi.index());
                auto const& src1 = mfsrc1.array(idx);
                auto const& src2 = mfsrc2.array(idx);
                AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM > 1);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
#if (AMREX_SPACEDIM == 2)
                    auto [jj1,w1] = is_nodal ? nd_interp_info(j)
                                             : cc_interp_info(j);
                    auto v1 = src1(i,jj1,k)*w1 + src1(i,jj1+1,k)*(1.0_rt-w1);
                    auto [jj2,w2] = is_nodal ? cc2nd_interp_info(j)
                                             : nd2cc_interp_info(j);
                    auto v2 = src2(i,jj2,k)*w2 + src2(i,jj2+1,k)*(1.0_rt-w2);
                    // Lorentz transformation
                    dst(i,j,k) = v1*c1 + v2*c2;
#else
                    auto [kk1,w1] = is_nodal ? nd_interp_info(k)
                                             : cc_interp_info(k);
                    auto v1 = src1(i,j,kk1)*w1 + src1(i,j,kk1+1)*(1.0_rt-w1);
                    auto [kk2,w2] = is_nodal ? cc2nd_interp_info(k)
                                             : nd2cc_interp_info(k);
                    auto v2 = src2(i,j,kk2)*w2 + src2(i,j,kk2+1)*(1.0_rt-w2);
                    // Lorentz transformation
                    dst(i,j,k) = v1*c1 + v2*c2;
#endif
                });
            }
        }
    }
}

void
WarpX::InitImposeFieldsGeom ()
{
    if (m_impose_field_type == ImposeFieldType::station) {
        Vector<char> headerfile;
        ParallelDescriptor::ReadAndBcastFile(m_impose_field_file_path+"/StationHeader", headerfile);
        std::istringstream is(std::string(headerfile.data()));
        is >> m_impose_station_info.location;
        std::string tt;
        std::getline(is,tt); // eat \n
        while (std::getline(is, tt)) {
            std::istringstream istt(tt);
            Real t0, t1;
            istt >> t0;
            istt >> t1;
            m_impose_station_info.time.emplace_back(t0,t1);
        }
        m_impose_station_info.ibuffer1 = -1;
        m_impose_station_info.ibuffer2 = -1;
    } else if (m_impose_field_type == ImposeFieldType::snapshot) {
        PlotFileData pf(m_impose_field_file_path);
        RealBox rb(pf.probLo(), pf.probHi());
        {
            const int zdir = AMREX_SPACEDIM-1;
            const auto zlen = rb.length(zdir);
            rb.setHi(zdir, rb.lo(zdir));
            rb.setLo(zdir, rb.hi(zdir)-zlen);
        }
        {
            ParmParse pp("geometry");
            pp.addarr("prob_lo", Vector<Real>{AMREX_D_DECL(rb.lo(0),rb.lo(1),rb.lo(2))});
            pp.addarr("prob_hi", Vector<Real>{AMREX_D_DECL(rb.hi(0),rb.hi(1),rb.hi(2))});
        }
        for (int lev = 0; lev <= pf.finestLevel(); ++lev) {
            m_impose_snapshot_info.geom.emplace_back(pf.probDomain(lev), rb, 0,
                                                     Array<int,AMREX_SPACEDIM>{AMREX_D_DECL(0,0,0)});
        }
        m_impose_snapshot_info.t_lab = pf.time();
        // This is only correct for gamma_boost=1. It will get updated for
        // gamma_boost>1 in ConvertLabParamsToBoost().
        m_t_boost_offset = m_impose_snapshot_info.t_lab;
    }
}

void
WarpX::ImposeFieldsInPlane ()
{
    if (!impose_E_field_in_plane || !impose_B_field_in_plane) return;

    AMREX_ALWAYS_ASSERT(boost_direction[0] != 1 && boost_direction[1] != 1);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(max_level == 0,
                                     "Imposing field in a plane is not implemented for more than one level");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!add_external_E_field && !add_external_B_field,
                                     "Cannot have both add-external-fields and impose-fields-in-plane");

    if (m_impose_field_type == ImposeFieldType::station) {
        ImposeStationFieldsInPlane();
    } else {
        ImposeSnapshotFieldsInPlane();
    }
}

void
WarpX::ImposeStationFieldsInPlane ()
{
    const int zdir = WARPX_ZINDEX;
    const Real t_boost = t_new[0];
    const Real z_lab = m_impose_station_info.location;
    const Real t_lab_min = m_impose_station_info.time.front().first;
    const Real t_lab_max = m_impose_station_info.time.back().second;
    const Real beta = beta_boost;
    const Real gamma = gamma_boost;

    // Give z_lab and t_boost, compute the impose location in boosted frame
    const Real zmid = z_lab/gamma - beta*PhysConst::c*t_boost;

    // Given z_boost, and t_boost, compute lab

    AMREX_ALWAYS_ASSERT(finestLevel() == 0 &&
                        Geom(0).Domain().smallEnd(zdir) == 0);
    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
        const auto ngz = std::max(Efield_fp[lev][0]->nGrowVect()[zdir],
                                  Bfield_fp[lev][0]->nGrowVect()[zdir]);
        const auto zlo  = Geom(lev).ProbLo(zdir);
        const auto dz   = Geom(lev).CellSize(zdir);
        const Real zmin = zmid - Real(0.5)*dz*ngz;
        const Real zmax = zmid - Real(0.5)*dz*ngz;
        const auto izmin = int(std::floor((zmin-zlo)/dz));
        const auto izmax = int(std::floor((zmax-zlo)/dz))+1;
        auto t_lab_min = gamma*(t_boost+beta/PhysConst::c*(zlo+dz*Real(izmin)));
        auto t_lab_max = gamma*(t_boost+beta/PhysConst::c*(zlo+dz*Real(izmax)));

        Box const nddom = amrex::surroundingNodes(Geom(lev).Domain());
        Box imposebox = nddom;
        imposebox.setSmall(zdir,izmin).setBig(zdir,izmax);

        if (! nddom.intersects(imposebox)) {
            return;
        }

        BoxArray ba = amrex::convert(Efield_fp[lev][0]->boxArray(),nddom.ixType());
        DistributionMapping const& dm = Efield_fp[lev][0]->DistributionMap();

        BoxList bl(ba.ixType());
        Vector<int> procmap;
    }

#if 0
    const Real t_lab = t_boost / gamma_boost + beta_boost * z_lab / PhysConst::c;
    const Real z_boost = gamma_boost * (z_lab - beta_boost*PhysConst::c*t_lab);

    if (t_lab < m_impose_station_info.time.front().first) {
        return;
    }

    if (t_lab >= m_impose_station_info.time.back().second) {
        return;
    }

    int ibuffer;
    const auto nbuffers = static_cast<int>(m_impose_station_info.time.size());
    for (ibuffer = 0; ibuffer < nbuffers; ++ibuffer) {
        if (t_lab >= m_impose_station_info.time[ibuffer].first &&
            t_lab <  m_impose_station_info.time[ibuffer].second) {
            break;
        }
    }
#if 0
    AMREX_ALWAYS_ASSERT(ibuffer != nbuffers);

    if (ibuffer != m_impose_station_info.ibuffer) {
        m_impose_station_info.ibuffer = ibuffer;
        amrex::Print() << "xxxxx ibuffer = " << ibuffer
                       << " " << t_lab << " " << m_impose_station_info.time[ibuffer].first << " " << m_impose_station_info.time[ibuffer].second << "\n";
        amrex::Abort("xxxxx need to read data");
    }
#endif
#endif
}

void
WarpX::ImposeSnapshotFieldsInPlane ()
{
    if (Efield_fp_external[0][0] == nullptr &&
        Bfield_fp_external[0][0] == nullptr) {
        for (int lev = 0; lev <= finestLevel(); ++lev) {
            std::string raw_field_path = m_impose_field_file_path;
            raw_field_path.append("/raw_fields/Level_")
                .append(std::to_string(lev)).append("/");

            BoxArray ba0;
            DistributionMapping dm0;

            if (impose_E_field_in_plane) {
                std::array<std::string,3> Exyz{"Ex_fp","Ey_fp","Ez_fp"};
                for (int idim = 0; idim < 3; ++idim) {
                    Efield_fp_external[lev][idim] = read_snapshot_field
                        (ba0, dm0, Efield_fp[lev][idim]->ixType(),
                         raw_field_path+Exyz[idim]);
                }
            }

            if (impose_B_field_in_plane) {
                std::array<std::string,3> Bxyz{"Bx_fp","By_fp","Bz_fp"};
                for (int idim = 0; idim < 3; ++idim) {
                    Bfield_fp_external[lev][idim] = read_snapshot_field
                        (ba0, dm0, Bfield_fp[lev][idim]->ixType(),
                         raw_field_path+Bxyz[idim]);
                }
            }
        }
    }

    const Real t = t_new[0];
    const Real ct = PhysConst::c * t;
    const Real beta = beta_boost;
    const Real gamma = gamma_boost;
    const Real betact = beta_boost * ct;

    auto z_lab_to_boost = [=] (Real zlab) -> Real
    {
        return zlab / gamma - betact;
    };

    auto z_boost_to_lab = [=] AMREX_GPU_HOST_DEVICE (Real zboost) -> Real
    {
        return gamma * (zboost + betact);
    };

    const Real ct_impose = PhysConst::c * m_impose_snapshot_info.t_lab;
    auto ct_lab = [=] AMREX_GPU_HOST_DEVICE (Real zboost) -> Real
    {
        return gamma * (ct + beta * zboost) - ct_impose;
    };

    AMREX_ALWAYS_ASSERT(finestLevel() == 0);
    for (int lev = 0; lev <= finestLevel(); ++lev)
    {
        const int zdir = AMREX_SPACEDIM-1;
        const Real zmax = z_lab_to_boost(m_impose_snapshot_info.geom[lev].ProbHi(zdir));
        const auto ngz = std::max(Efield_fp[lev][0]->nGrowVect()[zdir],
                                  Bfield_fp[lev][0]->nGrowVect()[zdir]);
        const auto zlo  = Geom(lev).ProbLo(zdir);
        const auto dz   = Geom(lev).CellSize(zdir);
        const Real zmin = zmax - dz*ngz;
        const auto izmin = int(std::floor((zmin-zlo)/dz));
        const auto izmax = int(std::floor((zmax-zlo)/dz))+1;

        const auto zlo_ext  = m_impose_snapshot_info.geom[lev].ProbLo(zdir);
        const auto dz_ext   = m_impose_snapshot_info.geom[lev].CellSize(zdir);

        AMREX_ALWAYS_ASSERT(Geom(lev).Domain().smallEnd(zdir) == 0 &&
                            m_impose_snapshot_info.geom[lev].Domain().smallEnd(zdir) == 0);

        // Given a cell index, this returns information for interpolation,
        // the lower index in the external field and its weight
        auto cc_interp_info = [=] AMREX_GPU_HOST_DEVICE (int k)
            -> std::pair<int,Real>
        {
            Real z = zlo + (k+0.5_rt)*dz;
            Real zext = z_boost_to_lab(z) - ct_lab(z);
            int kext = int(std::floor((zext-zlo_ext)/dz_ext));
            Real zcext = zlo_ext + (kext+0.5_rt)*dz_ext;
            if (zcext <= zext) {
                Real w = std::max(1.0_rt - (zext-zcext)/dz_ext, 0.0_rt);
                return {kext, w};
            } else {
                Real w = std::min((zcext-zext)/dz_ext, 1.0_rt);
                return {kext-1, w};
            }
        };

        // Given a nodal index, this returns information for interpolation,
        // the lower index in the external field and its weight
        auto nd_interp_info = [=] AMREX_GPU_HOST_DEVICE (int k)
            -> std::pair<int,Real>
        {
            Real z = zlo + k*dz;
            Real zext = z_boost_to_lab(z) - ct_lab(z);
            int kext = int(std::floor((zext-zlo_ext)/dz_ext));
            Real w = 1.0_rt - (zext - (zlo_ext+kext*dz_ext)) / dz_ext;
            w = std::max(std::min(w,1.0_rt),0.0_rt);
            return {kext, w};
        };

        // In boosted frame, we need to interpolate from cc to nd and nd to cc.
        auto nd2cc_interp_info = [=] AMREX_GPU_HOST_DEVICE (int k) // k is cc
            -> std::pair<int,Real>
        {
            Real z = zlo + (k+0.5_rt)*dz;
            Real zext = z_boost_to_lab(z) - ct_lab(z);
            int kext = int(std::floor((zext-zlo_ext)/dz_ext));
            Real w = 1.0_rt - (zext - (zlo_ext+kext*dz_ext)) / dz_ext;
            w = std::max(std::min(w,1.0_rt),0.0_rt);
            return {kext, w};
        };

        auto cc2nd_interp_info = [=] AMREX_GPU_HOST_DEVICE (int k) // k is nd
            -> std::pair<int,Real>
        {
            Real z = zlo + k*dz;
            Real zext = z_boost_to_lab(z) - ct_lab(z);
            int kext = int(std::floor((zext-zlo_ext)/dz_ext));
            Real zcext = zlo_ext + (kext+0.5_rt)*dz_ext;
            if (zcext <= zext) {
                Real w = std::max(1.0_rt - (zext-zcext)/dz_ext, 0.0_rt);
                return {kext, w};
            } else {
                Real w = std::min((zcext-zext)/dz_ext, 1.0_rt);
                return {kext-1, w};
            }
        };

        Box planebox = Geom(lev).Domain();
        planebox.surroundingNodes(zdir).setSmall(zdir, izmin).setBig(zdir, izmax);

        BoxArray ba = Efield_fp[lev][0]->boxArray();
        ba.enclosedCells().surroundingNodes(zdir);

        DistributionMapping const& dm = Efield_fp[lev][0]->DistributionMap();

        BoxList bl(ba.ixType());
        Vector<int> procmap;
        std::map<int,int> index_map;

        int nboxes = ba.size();
        int kmin_ext = std::numeric_limits<int>::max();
        int kmax_ext = std::numeric_limits<int>::lowest();
        for (int ibox = 0; ibox < nboxes; ++ibox) {
            Box b = ba[ibox] & planebox;
            if (b.ok()) {
                auto cclo = cc_interp_info(b.smallEnd(zdir));
                auto cchi = cc_interp_info(b.bigEnd(zdir)-1);
                auto ndlo = nd_interp_info(b.smallEnd(zdir));
                auto ndhi = nd_interp_info(b.bigEnd(zdir));
                int lo = std::min(cclo.first  ,ndlo.first  );
                int hi = std::max(cchi.first+2,ndhi.first+1);
                if (beta_boost != 1.0_rt) {
                    cclo = cc2nd_interp_info(b.smallEnd(zdir));
                    cchi = cc2nd_interp_info(b.bigEnd(zdir));
                    ndlo = nd2cc_interp_info(b.smallEnd(zdir));
                    ndhi = nd2cc_interp_info(b.bigEnd(zdir)-1);
                    lo = std::min({lo, cclo.first  , ndlo.first  });
                    hi = std::max({hi, cchi.first+2, ndhi.first+1});
                }
                b.setSmall(zdir,lo);
                b.setBig  (zdir,hi);
                bl.push_back(b);
                procmap.push_back(dm[ibox]);
                index_map[ibox] = bl.size() - 1;
                kmin_ext = std::min(kmin_ext, b.smallEnd(zdir));
                kmax_ext = std::max(kmax_ext, b.bigEnd  (zdir));
            }
        }

        if (kmax_ext < kmin_ext) {
            // The plane is now out of the simulation domain.
            return;
        }

        BoxArray baext(std::move(bl));
        DistributionMapping dmext(std::move(procmap));

        if (beta_boost != 0.0_rt) {
            AMREX_ALWAYS_ASSERT(impose_E_field_in_plane && impose_B_field_in_plane);
        }

        if (impose_E_field_in_plane) {
            if (beta_boost == 0.0_rt) {
                for (int idim = 0; idim < 3; ++idim) {
                    impose_snapshot_field_direct(*Efield_fp         [lev][idim],
                                                 *Efield_fp_external[lev][idim],
                                                 baext, dmext, planebox, index_map,
                                                 nd_interp_info, cc_interp_info);
                }
            } else {
                impose_snapshot_field_boost(*Efield_fp         [lev][0],
                                            *Efield_fp_external[lev][0],
                                            gamma_boost,
                                            *Bfield_fp_external[lev][1],
                                            -gamma_boost*beta_boost*PhysConst::c,
                                            baext, dmext, planebox, index_map,
                                            nd_interp_info, cc_interp_info,
                                            nd2cc_interp_info, cc2nd_interp_info);

                impose_snapshot_field_boost(*Efield_fp         [lev][1],
                                            *Efield_fp_external[lev][1],
                                            gamma_boost,
                                            *Bfield_fp_external[lev][0],
                                            gamma_boost*beta_boost*PhysConst::c,
                                            baext, dmext, planebox, index_map,
                                            nd_interp_info, cc_interp_info,
                                            nd2cc_interp_info, cc2nd_interp_info);

                impose_snapshot_field_direct(*Efield_fp         [lev][2],
                                             *Efield_fp_external[lev][2],
                                             baext, dmext, planebox, index_map,
                                             nd_interp_info, cc_interp_info);
            }
        }

        if (impose_B_field_in_plane) {
            if (beta_boost == 0.0_rt) {
                for (int idim = 0; idim < 3; ++idim) {
                    impose_snapshot_field_direct(*Bfield_fp         [lev][idim],
                                                 *Bfield_fp_external[lev][idim],
                                                 baext, dmext, planebox, index_map,
                                                 nd_interp_info, cc_interp_info);
                }
            } else {
                impose_snapshot_field_boost(*Bfield_fp         [lev][0],
                                            *Bfield_fp_external[lev][0],
                                            gamma_boost,
                                            *Efield_fp_external[lev][1],
                                            gamma_boost*beta_boost/PhysConst::c,
                                            baext, dmext, planebox, index_map,
                                            nd_interp_info, cc_interp_info,
                                            nd2cc_interp_info, cc2nd_interp_info);

                impose_snapshot_field_boost(*Bfield_fp         [lev][1],
                                            *Bfield_fp_external[lev][1],
                                            gamma_boost,
                                            *Efield_fp_external[lev][0],
                                            -gamma_boost*beta_boost/PhysConst::c,
                                            baext, dmext, planebox, index_map,
                                            nd_interp_info, cc_interp_info,
                                            nd2cc_interp_info, cc2nd_interp_info);

                impose_snapshot_field_direct(*Bfield_fp         [lev][2],
                                             *Bfield_fp_external[lev][2],
                                             baext, dmext, planebox, index_map,
                                             nd_interp_info, cc_interp_info);
            }
        }
    }
}
