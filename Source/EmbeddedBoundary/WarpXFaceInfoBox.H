//
// Created by Lorenzo Giacomel on 14/04/2021.
//
#ifndef WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H
#define WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H

#include<bitset>

struct FaceInfoBox {
    //amrex::Gpu::DeviceVector<int> i_face, j_face, k_face;
    amrex::Gpu::DeviceVector<uint8_t> neigh_faces;
    amrex::Gpu::DeviceVector<amrex::Real> area;
    // this vector stores the indices to access in the previous vectors
    amrex::Gpu::DeviceVector<int> inds;
    //number of entries of inds which correspond to this element
    amrex::BaseFab<int> size;
    //each entry points an element inside the inds vector all the elements. After that element
    // #size entries indicate the indices to access the vectors
    amrex::BaseFab<int*> inds_pointer;

    static const uint8_t ul = 1;
    static const uint8_t l = 2;
    static const uint8_t bl = 4;
    static const uint8_t u = 8;
    static const uint8_t d = 16;
    static const uint8_t ur = 32;
    static const uint8_t r = 64;
    static const uint8_t br = 128;

    /**
    * \brief add the neighbor i, j to the list of intruded neighbors.
    */
    static void addConnectedNeighbor(int i, int j, int ind, uint8_t* neigh_face_ptr){

        if(i == -1 and j == -1){
            *(neigh_face_ptr + ind) = ul;
        }else if(i == -1 and j == 0){
            *(neigh_face_ptr + ind) = l;
        }else if(i == -1 and j == 1){
            *(neigh_face_ptr + ind) = bl;
        }else if(i == 0 and j == -1){
            *(neigh_face_ptr + ind) = u;
        }else if(i == 0 and j == 1){
            *(neigh_face_ptr + ind) = d;
        }else if(i == 1 and j == -1){
            *(neigh_face_ptr + ind) = ur;
        }else if(i == 1 and j == 0){
            *(neigh_face_ptr + ind) = r;
        }else if(i == 1 and j == 1){
            *(neigh_face_ptr + ind) = br;
        }
    }

/**
* \brief writes into i_face and j_face the intruded neighbors indices;
*/
    static amrex::Array1D<int, 0, 1> uint8_to_inds(uint8_t mask){
        amrex::Array1D<int, 0, 1> res;

        if(mask & ul){
            res(0) = -1;
            res(1) = -1;
        }else if(mask & l){
            res(0) = -1;
            res(1) = 0;
        }else if(mask & bl){
            res(0) = -1;
            res(1) = 1;
        }else if(mask & u){
            res(0) = 0;
            res(1) = -1;
        }else if(mask & d){
            res(0) = -1;
            res(1) = 1;
        }else if(mask & ur){
            res(0) = 1;
            res(1) = -1;
        }else if(mask & r){
            res(0) = 1;
            res(1) = 0;
        }else if(mask & br){
            res(0) = 1;
            res(1) = 1;
        }

        return res;
    }
};

#endif //WARPX_SOURCE_EMBEDDEDBOUNDARY_WARPXFACEINFOBOX_H
