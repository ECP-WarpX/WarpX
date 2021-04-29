//
// Created by Lorenzo Giacomel on 14/04/2021.
//

#ifndef WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H
#define WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H

struct FaceInfoBox {
  amrex::Gpu::DeviceVector<int> i_face, j_face, k_face;
  amrex::Gpu::DeviceVector<amrex::Real> area, rho_face;
  amrex::BaseFab<std::vector<int>> inds;
};

#endif //WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H
