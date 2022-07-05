/* Copyright 2021 Lorenzo Giacomel
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpX.H"

#ifdef AMREX_USE_EB
#  include "Utils/TextMsg.H"
#  include "Utils/WarpXUtil.H"

#  include <AMReX.H>
#  include <AMReX_Array.H>
#  include <AMReX_Array4.H>
#  include <AMReX_BLProfiler.H>
#  include <AMReX_Box.H>
#  include <AMReX_BoxArray.H>
#  include <AMReX_BoxList.H>
#  include <AMReX_Config.H>
#  include <AMReX_EB2.H>
#  include <AMReX_EB_utils.H>
#  include <AMReX_FabArray.H>
#  include <AMReX_FabFactory.H>
#  include <AMReX_GpuControl.H>
#  include <AMReX_GpuDevice.H>
#  include <AMReX_GpuQualifiers.H>
#  include <AMReX_IntVect.H>
#  include <AMReX_Loop.H>
#  include <AMReX_MFIter.H>
#  include <AMReX_MultiFab.H>
#  include <AMReX_iMultiFab.H>
#  include <AMReX_ParmParse.H>
#  include <AMReX_Parser.H>
#  include <AMReX_REAL.H>
#  include <AMReX_SPACE.H>
#  include <AMReX_Vector.H>

#  include <cstdlib>
#  include <string>

#endif

#ifdef AMREX_USE_EB
namespace {
    class ParserIF
        : public amrex::GPUable
    {
    public:
        ParserIF (const amrex::ParserExecutor<3>& a_parser)
            : m_parser(a_parser)
            {}

        ParserIF (const ParserIF& rhs) noexcept = default;
        ParserIF (ParserIF&& rhs) noexcept = default;
        ParserIF& operator= (const ParserIF& rhs) = delete;
        ParserIF& operator= (ParserIF&& rhs) = delete;

        AMREX_GPU_HOST_DEVICE inline
        amrex::Real operator() (AMREX_D_DECL(amrex::Real x, amrex::Real y,
                                             amrex::Real z)) const noexcept {
#if defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
            return m_parser(x,amrex::Real(0.0),y);
#else
            return m_parser(x,y,z);
#endif
        }

        inline amrex::Real operator() (const amrex::RealArray& p) const noexcept {
            return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
        }

    private:
        amrex::ParserExecutor<3> m_parser; //! function parser with three arguments (x,y,z)
    };
}
#endif

void
WarpX::InitEB ()
{
#ifdef AMREX_USE_EB
    BL_PROFILE("InitEB");

    amrex::ParmParse pp_warpx("warpx");
    std::string impf;
    pp_warpx.query("eb_implicit_function", impf);
    if (! impf.empty()) {
        auto eb_if_parser = makeParser(impf, {"x", "y", "z"});
        ParserIF pif(eb_if_parser.compile<3>());
        auto gshop = amrex::EB2::makeShop(pif, eb_if_parser);
         // The last argument of amrex::EB2::Build is the maximum coarsening level
         // to which amrex should try to coarsen the EB.  It will stop after coarsening
         // as much as it can, if it cannot coarsen to that level.  Here we use a big
         // number (e.g., maxLevel()+20) for multigrid solvers.  Because the coarse
         // level has only 1/8 of the cells on the fine level, the memory usage should
         // not be an issue.
        amrex::EB2::Build(gshop, Geom(maxLevel()), maxLevel(), maxLevel()+20);
    } else {
        amrex::ParmParse pp_eb2("eb2");
        if (!pp_eb2.contains("geom_type")) {
            std::string geom_type = "all_regular";
            pp_eb2.add("geom_type", geom_type); // use all_regular by default
        }
        // See the comment above on amrex::EB2::Build for the hard-wired number 20.
        amrex::EB2::Build(Geom(maxLevel()), maxLevel(), maxLevel()+20);
    }

#endif
}

#ifdef AMREX_USE_EB
void
WarpX::ComputeEdgeLengths (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& edge_lengths,
                           const amrex::EBFArrayBoxFactory& eb_fact) {
#ifndef WARPX_DIM_RZ
    BL_PROFILE("ComputeEdgeLengths");

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &edge_centroid = eb_fact.getEdgeCent();
#ifdef WARPX_DIM_XZ
    edge_lengths[1]->setVal(0.);
#endif
    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi){
#ifdef WARPX_DIM_XZ
        for (int idim = 0; idim < 3; ++idim){
            if(idim == 1) continue;
#elif defined(WARPX_DIM_3D)
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
#else
        amrex::Abort(Utils::TextMsg::Err(
            "ComputeEdgeLengths: Only implemented in 2D3V and 3D3V"));
#endif
            amrex::Box box = mfi.tilebox(edge_lengths[idim]->ixType().toIntVect(),
                                         edge_lengths[idim]->nGrowVect());
            amrex::FabType fab_type = flags[mfi].getType(box);
            auto const &edge_lengths_dim = edge_lengths[idim]->array(mfi);

            if (fab_type == amrex::FabType::regular) {
                // every cell in box is all regular
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    edge_lengths_dim(i, j, k) = 1.;
                });
            } else if (fab_type == amrex::FabType::covered) {
                // every cell in box is all covered
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    edge_lengths_dim(i, j, k) = 0.;
                });
            } else {
#ifdef WARPX_DIM_XZ
                int idim_amrex = idim;
                if(idim == 2) idim_amrex = 1;
                auto const &edge_cent = edge_centroid[idim_amrex]->const_array(mfi);
#elif defined(WARPX_DIM_3D)
                auto const &edge_cent = edge_centroid[idim]->const_array(mfi);
#else
                amrex::Abort(Utils::TextMsg::Err(
                    "ComputeEdgeLengths: Only implemented in 2D3V and 3D3V"));
#endif
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (edge_cent(i, j, k) == amrex::Real(-1.0)) {
                        // This edge is all covered
                        edge_lengths_dim(i, j, k) = 0.;
                    } else if (edge_cent(i, j, k) == amrex::Real(1.0)) {
                        // This edge is all open
                        edge_lengths_dim(i, j, k) = 1.;
                    } else {
                        // This edge is cut.
                        edge_lengths_dim(i, j, k) = 1 - amrex::Math::abs(amrex::Real(2.0)
                                                                        * edge_cent(i, j, k));
                    }

                });
            }
        }
    }
#endif
}


void
WarpX::ComputeFaceAreas (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& face_areas,
                         const amrex::EBFArrayBoxFactory& eb_fact) {
#ifndef WARPX_DIM_RZ
    BL_PROFILE("ComputeFaceAreas");

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
#ifdef WARPX_DIM_XZ
    //In 2D the volume frac is actually the area frac.
    auto const &area_frac = eb_fact.getVolFrac();
#elif defined(WARPX_DIM_3D)
    auto const &area_frac = eb_fact.getAreaFrac();
#else
    amrex::Abort(Utils::TextMsg::Err(
        "ComputeFaceAreas: Only implemented in 2D3V and 3D3V"));
#endif

#ifdef WARPX_DIM_XZ
    face_areas[0]->setVal(0.);
    face_areas[2]->setVal(0.);
#endif
    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
#ifdef WARPX_DIM_XZ
        // In 2D we change the extrema of the for loop so that we only have the case idim=1
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
#elif defined(WARPX_DIM_3D)
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
#else
        amrex::Abort(Utils::TextMsg::Err(
            "ComputeFaceAreas: Only implemented in 2D3V and 3D3V"));
#endif
            amrex::Box box = mfi.tilebox(face_areas[idim]->ixType().toIntVect(),
                                         face_areas[idim]->nGrowVect());
            amrex::FabType fab_type = flags[mfi].getType(box);
            auto const &face_areas_dim = face_areas[idim]->array(mfi);
            if (fab_type == amrex::FabType::regular) {
                // every cell in box is all regular
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    face_areas_dim(i, j, k) = amrex::Real(1.);
                });
            } else if (fab_type == amrex::FabType::covered) {
                // every cell in box is all covered
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    face_areas_dim(i, j, k) = amrex::Real(0.);
                });
            } else {
#ifdef WARPX_DIM_XZ
                auto const &face = area_frac.const_array(mfi);
#elif defined(WARPX_DIM_3D)
                auto const &face = area_frac[idim]->const_array(mfi);
#else
                amrex::Abort(Utils::TextMsg::Err(
                    "ComputeFaceAreas: Only implemented in 2D3V and 3D3V"));
#endif
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    face_areas_dim(i, j, k) = face(i, j, k);
                });
            }
        }
    }
#endif
}


void
WarpX::ScaleEdges (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& edge_lengths,
                   const std::array<amrex::Real,3>& cell_size) {
#ifndef WARPX_DIM_RZ
    BL_PROFILE("ScaleEdges");

    for (amrex::MFIter mfi(*edge_lengths[0]); mfi.isValid(); ++mfi) {
#ifdef WARPX_DIM_XZ
        for (int idim = 0; idim < 3; ++idim){
            if(idim == 1) continue;
#elif defined(WARPX_DIM_3D)
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
#else
        amrex::Abort(Utils::TextMsg::Err(
            "ScaleEdges: Only implemented in 2D3V and 3D3V"));
#endif
            const amrex::Box& box = mfi.tilebox(edge_lengths[idim]->ixType().toIntVect(),
                                                edge_lengths[idim]->nGrowVect() );
            auto const &edge_lengths_dim = edge_lengths[idim]->array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edge_lengths_dim(i, j, k) *= cell_size[idim];
            });
        }
    }
#endif
}

void
WarpX::ScaleAreas(std::array< std::unique_ptr<amrex::MultiFab>, 3 >& face_areas,
                  const std::array<amrex::Real,3>& cell_size) {
#ifndef WARPX_DIM_RZ
    BL_PROFILE("ScaleAreas");

    amrex::Real full_area;

    for (amrex::MFIter mfi(*face_areas[0]); mfi.isValid(); ++mfi) {
#ifdef WARPX_DIM_XZ
        // In 2D we change the extrema of the for loop so that we only have the case idim=1
        for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
#elif defined(WARPX_DIM_3D)
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
#else
        amrex::Abort(Utils::TextMsg::Err(
            "ScaleAreas: Only implemented in 2D3V and 3D3V"));
#endif
            const amrex::Box& box = mfi.tilebox(face_areas[idim]->ixType().toIntVect(),
                                                face_areas[idim]->nGrowVect() );
#ifdef WARPX_DIM_XZ
            full_area = cell_size[0]*cell_size[2];
#elif defined(WARPX_DIM_3D)
            if (idim == 0) {
                full_area = cell_size[1]*cell_size[2];
            } else if (idim == 1) {
                full_area = cell_size[0]*cell_size[2];
            } else {
                full_area = cell_size[0]*cell_size[1];
            }
#else
            amrex::Abort(Utils::TextMsg::Err(
                "ScaleAreas: Only implemented in 2D3V and 3D3V"));
#endif
            auto const &face_areas_dim = face_areas[idim]->array(mfi);

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                face_areas_dim(i, j, k) *= full_area;
            });

        }
    }
#endif
}


void
WarpX::MarkCells(){
#ifndef WARPX_DIM_RZ
    auto const &cell_size = CellSize(maxLevel());

#ifdef WARPX_DIM_3D
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
#elif defined(WARPX_DIM_XZ)
    m_flag_info_face[maxLevel()][0]->setVal(0.);
    m_flag_info_face[maxLevel()][2]->setVal(0.);
    m_flag_ext_face[maxLevel()][0]->setVal(0.);
    m_flag_ext_face[maxLevel()][2]->setVal(0.);
    // In 2D we change the extrema of the for loop so that we only have the case idim=1
    for (int idim = 1; idim < AMREX_SPACEDIM; ++idim) {
#else
    amrex::Abort(Utils::TextMsg::Err(
        "MarkCells: Only implemented in 2D3V and 3D3V"));
#endif
        for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {
            //amrex::Box const &box = mfi.tilebox(m_face_areas[maxLevel()][idim]->ixType().toIntVect());
            const amrex::Box& box = mfi.tilebox(m_face_areas[maxLevel()][idim]->ixType().toIntVect(),
                                                m_face_areas[maxLevel()][idim]->nGrowVect() );

            auto const &S = m_face_areas[maxLevel()][idim]->array(mfi);
            auto const &flag_info_face = m_flag_info_face[maxLevel()][idim]->array(mfi);
            auto const &flag_ext_face = m_flag_ext_face[maxLevel()][idim]->array(mfi);
            const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
            const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
            const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
            auto const &mod_areas_dim = m_area_mod[maxLevel()][idim]->array(mfi);

            const amrex::Real dx = cell_size[0];
            const amrex::Real dy = cell_size[1];
            const amrex::Real dz = cell_size[2];

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Minimal area for this cell to be stable
                mod_areas_dim(i, j, k) = S(i, j, k);
                double S_stab;
                if(idim == 0){
                    S_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                    lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
                }else if(idim == 1){
#ifdef WARPX_DIM_XZ
                    S_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j + 1, k) * dz,
                                             lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
#elif defined(WARPX_DIM_3D)
                    S_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k + 1) * dz,
                                             lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
#else
                    amrex::Abort(Utils::TextMsg::Err(
                        "MarkCells: Only implemented in 2D3V and 3D3V"));
#endif
                }else {
                    S_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                             ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
                }

                // Does this face need to be extended?
                // The difference between flag_info_face and flag_ext_face is that:
                //     - for every face flag_info_face contains a:
                //          * 0 if the face needs to be extended
                //          * 1 if the face is large enough to lend area to other faces
                //          * 2 if the face is actually intruded by other face
                //       Here we only take care of the first two cases. The entries corresponding
                //       to the intruded faces are going to be set in the function ComputeFaceExtensions
                //     - for every face flag_ext_face contains a:
                //          * 1 if the face needs to be extended
                //          * 0 otherwise
                //       In the function ComputeFaceExtensions, after the cells are extended, the
                //       corresponding entries in flag_ext_face are set to zero. This helps to keep
                //       track of which cells could not be extended
                flag_ext_face(i, j, k) = int(S(i, j, k) < S_stab && S(i, j, k) > 0);
                if(flag_ext_face(i, j, k)){
                    flag_info_face(i, j, k) = 0;
                }
                // Is this face available to lend area to other faces?
                // The criterion is that the face has to be interior and not already unstable itself
                if(int(S(i, j, k) > 0 && !flag_ext_face(i, j, k))) {
                    flag_info_face(i, j, k) = 1;
                }
            });
        }
    }
#endif
}
#endif

void
WarpX::ComputeDistanceToEB () {
#ifdef AMREX_USE_EB
    BL_PROFILE("ComputeDistanceToEB");
    const amrex::EB2::IndexSpace& eb_is = amrex::EB2::IndexSpace::top();
    for (int lev=0; lev<=maxLevel(); lev++) {
        const amrex::EB2::Level& eb_level = eb_is.getLevel(Geom(lev));
        auto const eb_fact = fieldEBFactory(lev);
        amrex::FillSignedDistance(*m_distance_to_eb[lev], eb_level, eb_fact, 1);
    }
#endif
}
