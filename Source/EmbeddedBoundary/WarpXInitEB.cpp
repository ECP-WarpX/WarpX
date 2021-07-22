#include "WarpX.H"
#include <math.h>

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
WarpX::InitEB (int lev)
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
        amrex::EB2::Build(gshop, Geom(lev), maxLevel(), lev);
    } else {
        amrex::ParmParse pp_eb2("eb2");
        if (!pp_eb2.contains("geom_type")) {
            std::string geom_type = "all_regular";
            pp_eb2.add("geom_type", geom_type); // use all_regular by default
        }
        amrex::EB2::Build(Geom(lev), maxLevel(), lev);
    }
#else
    amrex::ignore_unused(lev);
#endif
}

/**
 * \brief Compute the length of the mesh edges. Here the length is a value in [0, 1].
 *        An edge of length 0 is fully covered.
 */
void
WarpX::ComputeEdgeLengths (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& edge_lengths,
                           int lev, bool flag_cp) {
#ifdef AMREX_USE_EB
    BL_PROFILE("ComputeEdgeLengths");

    //This variable is equal to lev if this is a fine patch or
    // equal to lev -1 if this is a coarse patch
    int lev_loc = lev;
    if(flag_cp and lev > 0){
        lev_loc = lev -1;
    }
    auto const eb_fact = fieldEBFactory(lev_loc);

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &edge_centroid = eb_fact.getEdgeCent();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim){
        for (amrex::MFIter mfi(*edge_lengths[idim]); mfi.isValid(); ++mfi){
            amrex::Box const &box = mfi.tilebox(edge_lengths[idim]->ixType().toIntVect());
            amrex::FabType fab_type = flags[mfi].getType(box);
            auto const &edge_lengths_dim = edge_lengths[idim]->array(mfi);
            if (fab_type == amrex::FabType::regular) {
                amrex::ParallelFor(box, [=](int i, int j, int k) {
                    edge_lengths_dim(i, j, k) = 1.;
                });
            } else if (fab_type == amrex::FabType::covered) {
                amrex::ParallelFor(box, [=](int i, int j, int k) {
                    edge_lengths_dim(i, j, k) = 0.;
                });
            } else {
                auto const &edge_cent = edge_centroid[idim]->const_array(mfi);
                amrex::ParallelFor(box, [=](int i, int j, int k) {
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
#else
    amrex::ignore_unused(edge_lengths, lev, flag_cp);
#endif
}

/**
 * \brief Compute the area of the mesh faces. Here the area is a value in [0, 1].
 *        An edge of area 0 is fully covered.
 */
void
WarpX::ComputeFaceAreas (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& face_areas,
                         int lev, bool flag_cp) {
#ifdef AMREX_USE_EB
    BL_PROFILE("ComputeFaceAreas");

    //This variable is equal to lev if this is a fine patch or
    // equal to lev -1 if this is a coarse patch
    int lev_loc = lev;
    if(flag_cp and lev > 0){
        lev_loc = lev -1;
    }

    auto const eb_fact = fieldEBFactory(lev_loc);
    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &area_frac = eb_fact.getAreaFrac();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (amrex::MFIter mfi(*face_areas[idim]); mfi.isValid(); ++mfi) {
            amrex::Box const &box = mfi.tilebox(face_areas[idim]->ixType().toIntVect());
            amrex::FabType fab_type = flags[mfi].getType(box);
            auto const &face_areas_dim = face_areas[idim]->array(mfi);
            if (fab_type == amrex::FabType::regular) {
                // every cell in box is all regular
                amrex::ParallelFor(box, [=](int i, int j, int k) {
                    face_areas_dim(i, j, k) = amrex::Real(1.);
                });
            } else if (fab_type == amrex::FabType::covered) {
                // every cell in box is all covered
                amrex::ParallelFor(box, [=](int i, int j, int k) {
                    face_areas_dim(i, j, k) = amrex::Real(0.);
                });
            } else {
                auto const &face = area_frac[idim]->const_array(mfi);
                amrex::ParallelFor(box, [=](int i, int j, int k) {
                    face_areas_dim(i, j, k) = face(i, j, k);
                });
            }
        }
    }
#else
    amrex::ignore_unused(face_areas, lev, flag_cp);
#endif
}

/**
 * \brief Scale the edges lengths by the mesh width to obtain the real lengths.
 */
void
WarpX::ScaleEdges (std::array< std::unique_ptr<amrex::MultiFab>, 3 >& edge_lengths,
                   int lev, bool flag_cp) {
#ifdef AMREX_USE_EB
    BL_PROFILE("ScaleEdges");

    //This variable is equal to lev if this is a fine patch or
    // equal to lev -1 if this is a coarse patch
    int lev_loc = lev;
    if(flag_cp and lev > 0){
        lev_loc = lev -1;
    }

    auto const &cell_size = CellSize(lev_loc);
    auto const eb_fact = fieldEBFactory(lev_loc);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (amrex::MFIter mfi(*edge_lengths[idim]); mfi.isValid(); ++mfi) {
            amrex::Box const &box = mfi.tilebox(edge_lengths[idim]->ixType().toIntVect());
            auto const &edge_lengths_dim = edge_lengths[idim]->array(mfi);
            amrex::ParallelFor(box, [=](int i, int j, int k) {
                edge_lengths_dim(i, j, k) *= cell_size[idim];
            });
        }
    }
#else
    amrex::ignore_unused(edge_lengths, lev, flag_cp);
#endif
}

/**
 * \brief Scale the edges areas by the mesh width to obtain the real areas.
 */
void
WarpX::ScaleAreas(std::array< std::unique_ptr<amrex::MultiFab>, 3 >& face_areas,
                  int lev, bool flag_cp) {
#ifdef AMREX_USE_EB
    BL_PROFILE("ScaleAreas");

    //This variable is equal to lev if this is a fine patch or
    // equal to lev -1 if this is a coarse patch
    int lev_loc = lev;
    if(flag_cp and lev > 0){
        lev_loc = lev -1;
    }

    auto const& cell_size = CellSize(lev_loc);
    amrex::Real full_area;

    auto const eb_fact = fieldEBFactory(lev_loc);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (amrex::MFIter mfi(*face_areas[idim]); mfi.isValid(); ++mfi) {
            amrex::Box const &box = mfi.tilebox(face_areas[idim]->ixType().toIntVect());
            if (idim == 0) {
                full_area = cell_size[1]*cell_size[2];
            } else if (idim == 1) {
                full_area = cell_size[0]*cell_size[2];
            } else {
                full_area = cell_size[0]*cell_size[1];
            }
            auto const &face_areas_dim = face_areas[idim]->array(mfi);
            amrex::ParallelFor(box, [=](int i, int j, int k) {
                face_areas_dim(i, j, k) *= full_area;
            });
        }
    }
#else
    amrex::ignore_unused(face_areas, lev, flag_cp);
#endif
}

//0 for unst, 1 for stable and available, 2 for stable, available and intruded

/**
 * \brief Initialize information for cell extensions.
 *        The flags convention for m_flag_info_face is as follows
 *          - 0 for unstable cells
 *          - 1 for stable cells which have not been intruded
 *          - 2 for stable cells which have been intruded
 *        Here we cannot know if a cell is intruded or not so we initialize all stable cells with 1
 */
// TODO: add the inputs here
void
WarpX::MarkCells(std::array< std::unique_ptr<amrex::MultiFab>, 3 >& edge_lengths,
                 std::array< std::unique_ptr<amrex::MultiFab>, 3 >& face_areas,
                 std::array< std::unique_ptr<amrex::iMultiFab>, 3 >& flag_info_face,
                 std::array< std::unique_ptr<amrex::iMultiFab>, 3 >& flag_ext_face,
                 int lev, bool flag_cp){
#ifdef AMREX_USE_EB

    int lev_loc = lev;
    if(flag_cp and lev > 0){
        lev_loc = lev -1;
    }

    auto const &cell_size = CellSize(lev_loc);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (amrex::MFIter mfi(*Bfield_fp[lev][idim]); mfi.isValid(); ++mfi) {
            amrex::Box const &box = mfi.tilebox(face_areas[idim]->ixType().toIntVect());

            auto const &S = face_areas[idim]->array(mfi);
            auto const &flag_info_face_dim = flag_info_face[idim]->array(mfi);
            auto const &flag_ext_face_dim = flag_ext_face[idim]->array(mfi);
            const auto &lx = edge_lengths[0]->array(mfi);
            const auto &ly = edge_lengths[1]->array(mfi);
            const auto &lz = edge_lengths[2]->array(mfi);
            amrex::Real dx = cell_size[0];
            amrex::Real dy = cell_size[1];
            amrex::Real dz = cell_size[2];

            amrex::ParallelFor(box, [=](int i, int j, int k) {
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

                // Does this face need to be extended? This is the same as
                // flag_info_face(i, j, k) = 0 here but it is modified later to keep track o which
                // faces need multi-cell extension
                flag_ext_face_dim(i, j, k) = int(S(i, j, k) < S_stab and !amrex::isnan(S(i, j, k))
                                              and S(i, j, k) > 0);
                if(flag_ext_face_dim(i, j, k)){
                    flag_info_face_dim(i, j, k) = 0;
                }
                // Is this face available to lend area to other faces?
                // The criterion is that the face has to be interior and not already unstable itself
                if(int(S(i, j, k) > 0 and !flag_ext_face_dim(i, j, k))) {
                    flag_info_face_dim(i, j, k) = 1;
                }
            });
        }
    }
#else
    amrex::ignore_unused(edge_lengths, face_areas, flag_info_face, flag_ext_face, lev, flag_cp);
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
