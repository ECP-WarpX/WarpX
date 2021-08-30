/* Copyright 2021 Lorenzo Giacomel
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpX.H"

#ifdef AMREX_USE_EB
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
#if (AMREX_SPACEDIM == 2)
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

#ifdef WARPX_DIM_RZ
    amrex::Abort("Embedded Boundaries not implemented in RZ geometry");
#endif

    amrex::ParmParse pp_warpx("warpx");
    std::string impf;
    pp_warpx.query("eb_implicit_function", impf);
    if (! impf.empty()) {
        auto eb_if_parser = makeParser(impf, {"x", "y", "z"});
        ParserIF pif(eb_if_parser.compile<3>());
        auto gshop = amrex::EB2::makeShop(pif, eb_if_parser);
        amrex::EB2::Build(gshop, Geom(maxLevel()), maxLevel(), maxLevel());
    } else {
        amrex::ParmParse pp_eb2("eb2");
        if (!pp_eb2.contains("geom_type")) {
            std::string geom_type = "all_regular";
            pp_eb2.add("geom_type", geom_type); // use all_regular by default
        }
        amrex::EB2::Build(Geom(maxLevel()), maxLevel(), maxLevel());
    }

#endif
}


void
WarpX::ComputeEdgeLengths () {
#ifdef AMREX_USE_EB
    BL_PROFILE("ComputeEdgeLengths");

    auto const eb_fact = fieldEBFactory(maxLevel());

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &edge_centroid = eb_fact.getEdgeCent();
    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi){
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            const amrex::Box& box = mfi.tilebox(m_edge_lengths[maxLevel()][idim]->ixType().toIntVect(),
                                                m_edge_lengths[maxLevel()][idim]->nGrowVect() );
            amrex::FabType fab_type = flags[mfi].getType(box);
            auto const &edge_lengths_dim = m_edge_lengths[maxLevel()][idim]->array(mfi);
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
                auto const &edge_cent = edge_centroid[idim]->const_array(mfi);
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
WarpX::ComputeFaceAreas () {
#ifdef AMREX_USE_EB
    BL_PROFILE("ComputeFaceAreas");

    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &area_frac = eb_fact.getAreaFrac();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const amrex::Box& box = mfi.tilebox(m_face_areas[maxLevel()][idim]->ixType().toIntVect(),
                                                m_face_areas[maxLevel()][idim]->nGrowVect() );
            amrex::FabType fab_type = flags[mfi].getType(box);
            auto const &face_areas_dim = m_face_areas[maxLevel()][idim]->array(mfi);
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
                auto const &face = area_frac[idim]->const_array(mfi);
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    face_areas_dim(i, j, k) = face(i, j, k);
                });
            }
        }
    }
#endif
}


void
WarpX::ScaleEdges () {
#ifdef AMREX_USE_EB
    BL_PROFILE("ScaleEdges");

    auto const &cell_size = CellSize(maxLevel());
    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &flags = eb_fact.getMultiEBCellFlagFab();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const amrex::Box& box = mfi.tilebox(m_edge_lengths[maxLevel()][idim]->ixType().toIntVect(),
                                                m_edge_lengths[maxLevel()][idim]->nGrowVect() );
            auto const &edge_lengths_dim = m_edge_lengths[maxLevel()][idim]->array(mfi);
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                edge_lengths_dim(i, j, k) *= cell_size[idim];
            });
        }
    }
#endif
}


void
WarpX::ScaleAreas() {
#ifdef AMREX_USE_EB
    BL_PROFILE("ScaleAreas");

    auto const& cell_size = CellSize(maxLevel());
    amrex::Real full_area;

    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &flags = eb_fact.getMultiEBCellFlagFab();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const amrex::Box& box = mfi.tilebox(m_face_areas[maxLevel()][idim]->ixType().toIntVect(),
                                                m_face_areas[maxLevel()][idim]->nGrowVect() );
            if (idim == 0) {
                full_area = cell_size[1]*cell_size[2];
            } else if (idim == 1) {
                full_area = cell_size[0]*cell_size[2];
            } else {
                full_area = cell_size[0]*cell_size[1];
            }
            auto const &face_areas_dim = m_face_areas[maxLevel()][idim]->array(mfi);

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                face_areas_dim(i, j, k) *= full_area;
            });

            if(WarpX::maxwell_solver_id==MaxwellSolverAlgo::ECT) {
                auto const &mod_areas_dim = m_area_mod[maxLevel()][idim]->array(mfi);
                amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        mod_areas_dim(i, j, k) = face_areas_dim(i, j, k);
                });
            }
        }
    }
#endif
}


void
WarpX::MarkCells(){
#ifdef AMREX_USE_EB
    auto const &cell_size = CellSize(maxLevel());

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
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
            const amrex::Real dx = cell_size[0];
            const amrex::Real dy = cell_size[1];
            const amrex::Real dz = cell_size[2];

            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                // Minimal area for this cell to be stable
                double S_stab;
                if(idim == 0){
                    S_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                    lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
                }else if(idim == 1){
                    S_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k + 1) * dz,
                                             lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
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


void
WarpX::ComputeDistanceToEB () {
#ifdef AMREX_USE_EB
    BL_PROFILE("ComputeDistanceToEB");

    amrex::ParmParse pp_warpx("warpx");
    std::string impf;
    pp_warpx.query("eb_implicit_function", impf);
    if (! impf.empty()) {
        auto eb_if_parser = makeParser(impf, {"x", "y", "z"});
        ParserIF pif(eb_if_parser.compile<3>());
        auto gshop = amrex::EB2::makeShop(pif);
        amrex::FillImpFunc(*m_distance_to_eb[maxLevel()], gshop, Geom(maxLevel()));
        m_distance_to_eb[maxLevel()]->negate(m_distance_to_eb[maxLevel()]->nGrow()); // signed distance f = - imp. f.
    } else {
        m_distance_to_eb[maxLevel()]->setVal(100.0); // some positive value
    }
#endif
}
