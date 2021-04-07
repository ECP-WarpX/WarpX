#include "WarpX.H"

#include <AMReX_Config.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EB2.H>
#   include <AMReX_ParmParse.H>
#endif


void
WarpX::InitEB ()
{
#ifdef AMREX_USE_EB
    BL_PROFILE("InitEB");

    amrex::ParmParse pp_eb2("eb2");
    if (!pp_eb2.contains("geom_type")) {
        pp_eb2.add("geom_type", "all_regular"); // use all_regular by default
    }
    amrex::EB2::Build(Geom(maxLevel()), maxLevel(), maxLevel());

#endif
}

void
WarpX::ComputeEdgeLengths() {
    BL_PROFILE("ComputeEdgeLengths");

    auto const eb_fact = fieldEBFactory(maxLevel());

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &edge_centroid = eb_fact.getEdgeCent();
//  	for ( amrex::MFIter mfi(*Bfield_fp[maxLevel()][0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        amrex::FabType fab_type = flags[mfi].getType(box);
        if (fab_type == amrex::FabType::regular) {
            // every cell in box is all regular
        } else if (fab_type == amrex::FabType::covered) {
            // every cell in box is all covered
        } else {
            // mix of regular, cut and covered
            // Note that edge_centroid[idim] is a MultiCutFab that does not have
            // data if the fab type regular or covered.
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const &edge = edge_centroid[idim]->const_array(mfi);
                auto const &edge_lengths_dim = edge_lengths[maxLevel()][idim]->array(mfi);
                amrex::LoopOnCpu(amrex::convert(box, amrex::Box(edge).ixType()),
                                 [=](int i, int j, int k) {
                                     if (edge(i, j, k) == amrex::Real(-1.0)) {
                                         // This edge is all covered
                                         edge_lengths_dim(i, j, k) = 0.;
                                     } else if (edge(i, j, k) == amrex::Real(1.0)) {
                                         // This edge is all open
                                         edge_lengths_dim(i, j, k) = 1.;
                                     } else {
                                         // This edge is cut.
                                         amrex::Real edge_cent = edge(i, j, k); // edge centroid: (-0.5,0.5)
                                         if (edge_cent < amrex::Real(0.0)) {
                                             // The right side is covered.
                                             edge_lengths_dim(i, j, k)  =
                                                     amrex::Real(2.0) * edge_cent + amrex::Real(0.5) + 0.5; // (0, 1)
                                         } else {
                                             // The left side is covered
                                             // TODO: rewrite well
                                             edge_lengths_dim(i, j, k)  = 1 -
                                                     (amrex::Real(2.0) * edge_cent - amrex::Real(0.5) + 0.5); // (0, 1)
                                         }
                                     }
                                 });
            }
        }

    }

}

void
WarpX::ComputeFaceAreas() {
    BL_PROFILE("ComputeFaceAreas");

    auto const eb_fact = fieldEBFactory(maxLevel());

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &area_frac = eb_fact.getAreaFrac();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        amrex::FabType fab_type = flags[mfi].getType(box);
        if (fab_type == amrex::FabType::regular) {
            // every cell in box is all regular
        } else if (fab_type == amrex::FabType::covered) {
            // every cell in box is all covered
        } else {
            // mix of regular, cut and covered
            // Note that edge_centroid[idim] is a MultiCutFab that does not have
            // data if the fab type regular or covered.
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const &face = area_frac[idim]->const_array(mfi);
                auto const &face_areas_dim = face_areas[maxLevel()][idim]->array(mfi);
                amrex::LoopOnCpu(amrex::convert(box, amrex::Box(face).ixType()),
                                 [=](int i, int j, int k) {
                                    face_areas_dim(i, j, k) = face(i, j, k);
                                 });
            }
        }

    }

}

void
WarpX::ScaleEdges() {
    auto const& cell_size = CellSize(maxLevel());

    auto const eb_fact = fieldEBFactory(maxLevel());

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &edge_centroid = eb_fact.getEdgeCent();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        amrex::FabType fab_type = flags[mfi].getType(box);
        if (fab_type == amrex::FabType::regular) {
            // every cell in box is all regular
        } else if (fab_type == amrex::FabType::covered) {
            // every cell in box is all covered
        } else {
            // mix of regular, cut and covered
            // Note that edge_centroid[idim] is a MultiCutFab that does not have
            // data if the fab type regular or covered.
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                auto const &edge = edge_centroid[1]->const_array(mfi);
			    auto const &edge_lengths_dim = edge_lengths[maxLevel()][idim]->array(mfi);
                for(int i = lo.x; i <= hi.x; ++i) {
                    for(int j = lo.y; j <= hi.y; ++j) {
                        edge_lengths_dim(i,j,hi.z/2) *= cell_size[idim];
                    }
                }
            }
        }
    }
}

void
WarpX::ScaleAreas() {
    auto const& cell_size = CellSize(maxLevel());
    amrex::Real full_area;

    auto const eb_fact = fieldEBFactory(maxLevel());

    auto const &flags = eb_fact.getMultiEBCellFlagFab();
    auto const &area_frac = eb_fact.getAreaFrac();

    for (amrex::MFIter mfi(flags); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        amrex::FabType fab_type = flags[mfi].getType(box);
        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        if (fab_type == amrex::FabType::regular) {
            // every cell in box is all regular
        } else if (fab_type == amrex::FabType::covered) {
            // every cell in box is all covered
        } else {
            // mix of regular, cut and covered
            // Note that edge_centroid[idim] is a MultiCutFab that does not have
            // data if the fab type regular or covered.
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                if (idim == 0) {
                    full_area = cell_size[1]*cell_size[2];
                } else if (idim == 1) {
                    full_area = cell_size[0]*cell_size[2];
                } else {
                    full_area = cell_size[0]*cell_size[1];
                }
                auto const &face = area_frac[idim]->const_array(mfi);
			    auto const &face_areas_dim = face_areas[maxLevel()][idim]->array(mfi);
                for(int i = lo.x; i <= hi.x; ++i) {
                    for(int j = lo.y; j <= hi.y; ++j) {
                        face_areas_dim(i,j,hi.z/2) *= full_area;
                    }
                }
            }
        }
    }
}
