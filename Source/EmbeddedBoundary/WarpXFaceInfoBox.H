//
// Created by Lorenzo Giacomel on 14/04/2021.
//

#ifndef WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H
#define WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H

struct FaceInfoBox {
  amrex::Gpu::DeviceVector<int> i_face, j_face, k_face;
  amrex::Gpu::DeviceVector<amrex::Real> area;
  // this vector stores the indices to access in the previous vectors
  amrex::Gpu::DeviceVector<int> inds;
  //number of entries of inds which correspond to this element
  amrex::BaseFab<int> size;
  //each entry points an element inside the inds vector all the elements. After that element
  // #size entries indicate the indices to access the vectors
  amrex::BaseFab<int*> inds_pointer;
};

#endif //WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H
