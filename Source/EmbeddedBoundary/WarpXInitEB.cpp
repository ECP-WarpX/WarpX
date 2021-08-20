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
#  include <AMReX_ParmParse.H>
#  include <AMReX_Parser.H>
#  include <AMReX_REAL.H>
#  include <AMReX_SPACE.H>
#  include <AMReX_Vector.H>

#  include <array>
#  include <cstdlib>
#  include <memory>
#  include <string>
#  include <vector>

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

/**
 * \brief Compute the length of the mesh edges. Here the length is a value in [0, 1].
 *        An edge of length 0 is fully covered.
 */
void
WarpX::ComputeEdgeLengths () {
#ifdef AMREX_USE_EB
    BL_PROFILE("ComputeEdgeLengths");

    auto const eb_fact = fieldEBFactory(maxLevel());

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &edge_centroid = eb_fact.getEdgeCent();
    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi){
        amrex::Box const &box = mfi.validbox();
        amrex::FabType fab_type = flags[mfi].getType(box);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
            auto const &edge_lengths_dim = m_edge_lengths[maxLevel()][idim]->array(mfi);
            if (fab_type == amrex::FabType::regular) {
                // every cell in box is all regular
                amrex::LoopOnCpu(amrex::convert(box, amrex::Box(edge_lengths_dim).ixType()),
                                 [=](int i, int j, int k) {
                    edge_lengths_dim(i, j, k) = 1.;
                });
            } else if (fab_type == amrex::FabType::covered) {
                // every cell in box is all covered
                amrex::LoopOnCpu(amrex::convert(box, amrex::Box(edge_lengths_dim).ixType()),
                                 [=](int i, int j, int k) {
                    edge_lengths_dim(i, j, k) = 0.;
                });
            } else {
                auto const &edge_cent = edge_centroid[idim]->const_array(mfi);
                amrex::LoopOnCpu(amrex::convert(box, amrex::Box(edge_cent).ixType()),
                                [=](int i, int j, int k) {
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

/**
 * \brief Compute the area of the mesh faces. Here the area is a value in [0, 1].
 *        An edge of area 0 is fully covered.
 */
void
WarpX::ComputeFaceAreas () {
#ifdef AMREX_USE_EB
    BL_PROFILE("ComputeFaceAreas");

    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &area_frac = eb_fact.getAreaFrac();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        amrex::FabType fab_type = flags[mfi].getType(box);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            auto const &face_areas_dim = m_face_areas[maxLevel()][idim]->array(mfi);
            if (fab_type == amrex::FabType::regular) {
                // every cell in box is all regular
                amrex::LoopOnCpu(amrex::convert(box, amrex::Box(face_areas_dim).ixType()),
                                [=](int i, int j, int k) {
                    face_areas_dim(i, j, k) = amrex::Real(1.);
                });
            } else if (fab_type == amrex::FabType::covered) {
                // every cell in box is all covered
                amrex::LoopOnCpu(amrex::convert(box, amrex::Box(face_areas_dim).ixType()),
                                 [=](int i, int j, int k) {
                    face_areas_dim(i, j, k) = amrex::Real(0.);
                });
            } else {
                auto const &face = area_frac[idim]->const_array(mfi);
                amrex::LoopOnCpu(amrex::convert(box, amrex::Box(face).ixType()),
                                [=](int i, int j, int k) {
                    face_areas_dim(i, j, k) = face(i, j, k);
                });
            }
        }
    }
#endif
}

/**
 * \brief Scale the edges lengths by the mesh width to obtain the real lengths.
 */
void
WarpX::ScaleEdges () {
#ifdef AMREX_USE_EB
    BL_PROFILE("ScaleEdges");

    auto const &cell_size = CellSize(maxLevel());
    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &flags = eb_fact.getMultiEBCellFlagFab();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            auto const &edge_lengths_dim = m_edge_lengths[maxLevel()][idim]->array(mfi);
            amrex::LoopOnCpu(amrex::convert(box, amrex::Box(edge_lengths_dim).ixType()),
                             [=](int i, int j, int k) {
                edge_lengths_dim(i, j, k) *= cell_size[idim];
            });
        }
    }
#endif
}

/**
 * \brief Scale the edges areas by the mesh width to obtain the real areas.
 */
void
WarpX::ScaleAreas() {
#ifdef AMREX_USE_EB
    BL_PROFILE("ScaleAreas");

    auto const& cell_size = CellSize(maxLevel());
    amrex::Real full_area;

    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &flags = eb_fact.getMultiEBCellFlagFab();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (idim == 0) {
                full_area = cell_size[1]*cell_size[2];
            } else if (idim == 1) {
                full_area = cell_size[0]*cell_size[2];
            } else {
                full_area = cell_size[0]*cell_size[1];
            }
            auto const &face_areas_dim = m_face_areas[maxLevel()][idim]->array(mfi);
            amrex::LoopOnCpu(amrex::convert(box, amrex::Box(face_areas_dim).ixType()),
                             [=](int i, int j, int k) {
                                face_areas_dim(i, j, k) *= full_area;
            });
        }
    }
#endif
}

/**
 * \brief Compute the level set function used for particle-boundary interaction.
 */
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
