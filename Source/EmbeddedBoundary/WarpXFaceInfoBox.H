/* Copyright 2021 Lorenzo Giacomel
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H
#define WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H

#include "WarpXFaceInfoBox_fwd.H"

#include <AMReX_Gpu.H>
#include <AMReX_BaseFab.H>

#include <bitset>
#include <map>
#include <utility>

struct FaceInfoBox {
    enum class Neighbours {n, s, e, w, nw, ne, sw, se};

    amrex::Gpu::DeviceVector<Neighbours> neigh_faces;
    amrex::Gpu::DeviceVector<amrex::Real> area;
    // this vector stores the indices to access in the previous vectors
    amrex::Gpu::DeviceVector<int> inds;
    //number of entries of inds which correspond to this element
    amrex::BaseFab<int> size;
    // Each entry points to an element inside the inds vector
    amrex::BaseFab<int*> inds_pointer;

    int vecs_size;

    /**
    * \brief add the neighbor i, j to the list of intruded neighbors.
    */
    AMREX_GPU_HOST_DEVICE
    static void addConnectedNeighbor(int i, int j, int ind, Neighbours* neigh_face_ptr)
    {
        static const auto m_i_j_to_neig = std::map<std::pair<int,int>, Neighbours>{
            {{-1,-1}, Neighbours::nw}, {{-1,0}, Neighbours::w}, {{-1,1}, Neighbours::sw},
            {{0 ,-1}, Neighbours::n }, {{0 ,1}, Neighbours::s},
            {{1 ,-1}, Neighbours::ne}, {{1 ,0}, Neighbours::e}, {{1 ,1}, Neighbours::se}};

        neigh_face_ptr[ind] = m_i_j_to_neig.at(std::make_pair(i,j));
    }

    /**
    * \brief writes into i_face and j_face the intruded neighbors indices;
    */
    AMREX_GPU_HOST_DEVICE
    static amrex::Array1D<int, 0, 1> uint8_to_inds(Neighbours mask)
    {
        static const auto m_neig_to_i_j = std::map<Neighbours,std::pair<int,int>>{
            {Neighbours::nw, {-1,-1}}, {Neighbours::w, {-1,0}}, {Neighbours::sw, {-1,1}},
            {Neighbours::n , {0 ,-1}}, {Neighbours::s, {0 ,1}},
            {Neighbours::ne, {1 ,-1}}, {Neighbours::e, {1 ,0}}, {Neighbours::se, {1 ,1}}};

        const auto [i,j] = m_neig_to_i_j.at(mask);
        return amrex::Array1D<int, 0, 1>{i, j};
    }
};

#endif //WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H
