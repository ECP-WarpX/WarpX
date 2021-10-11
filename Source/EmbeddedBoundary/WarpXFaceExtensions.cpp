/* Copyright 2021 Lorenzo Giacomel
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "WarpX.H"
#include <AMReX_Scan.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>


amrex::Array1D<int, 0, 2>
WarpX::CountExtFaces() {
    amrex::Array1D<int, 0, 2> sums{0, 0, 0};
#ifdef AMREX_USE_EB

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        amrex::ReduceOps<amrex::ReduceOpSum> reduce_ops;
        amrex::ReduceData<int> reduce_data(reduce_ops);
        for (amrex::MFIter mfi(*m_flag_ext_face[maxLevel()][idim]); mfi.isValid(); ++mfi) {
            amrex::Box const &box = mfi.validbox();
            auto const &flag_ext_face = m_flag_ext_face[maxLevel()][idim]->array(mfi);
            reduce_ops.eval(box, reduce_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> amrex::GpuTuple<int> {
                    return flag_ext_face(i, j, k);
                });
        }

        auto r = reduce_data.value();
        sums(idim) = amrex::get<0>(r);
    }

    amrex::ParallelDescriptor::ReduceIntSum(&(sums(0)), AMREX_SPACEDIM);
#endif
    return sums;
}


void
WarpX::ComputeFaceExtensions(){
#ifdef AMREX_USE_EB
    if(WarpX::verbose) {
        amrex::Array1D<int, 0, 2> N_ext_faces = CountExtFaces();
        amrex::Print() << "Faces to be extended in x:\t" << N_ext_faces(0) << std::endl;
        amrex::Print() << "Faces to be extended in y:\t" << N_ext_faces(1) << std::endl;
        amrex::Print() << "Faces to be extended in z:\t" << N_ext_faces(2) << std::endl;
    }

    InitBorrowing();
    ComputeOneWayExtensions();

    if(WarpX::verbose) {
        amrex::Array1D<int, 0, 2> N_ext_faces_after_one_way = CountExtFaces();
        amrex::Print() << "Faces to be extended after one way extension in x:\t" <<
                       N_ext_faces_after_one_way(0) << std::endl;
        amrex::Print() << "Faces to be extended after one way extension in y:\t" <<
                       N_ext_faces_after_one_way(1) << std::endl;
        amrex::Print() << "Faces to be extended after one way extension in z:\t" <<
                       N_ext_faces_after_one_way(2) << std::endl;
    }

    ComputeEightWaysExtensions();
    ShrinkBorrowing();

    amrex::Array1D<int, 0, 2> N_ext_faces_after_eight_ways = CountExtFaces();
    if(WarpX::verbose) {
        amrex::Print() << "Faces to be extended after eight ways extension in x:\t" <<
                       N_ext_faces_after_eight_ways(0) << std::endl;
        amrex::Print() << "Faces to be extended after eight ways extension in y:\t" <<
                       N_ext_faces_after_eight_ways(1) << std::endl;
        amrex::Print() << "Faces to be extended after eight ways extension ins z:\t" <<
                       N_ext_faces_after_eight_ways(2) << std::endl;
    }
    if (N_ext_faces_after_eight_ways(0) > 0) {
        amrex::Abort("Some x faces could not be extended");
    }
    if (N_ext_faces_after_eight_ways(1) > 0) {
        amrex::Abort("Some y faces could not be extended");
    }
    if (N_ext_faces_after_eight_ways(2) > 0) {
        amrex::Abort("Some z faces could not be extended");
    }
#endif
}


void
WarpX::InitBorrowing() {
    int idim = 0;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        auto &borrowing_x = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_x.inds_pointer.resize(box);
        borrowing_x.size.resize(box);
        borrowing_x.size.setVal<amrex::RunOn::Device>(0);
        amrex::Long ncells = box.numPts();
        // inds, neigh_faces and area are extended to their largest possible size here, but they are
        // resized to a much smaller size later on, based on the actual number of neighboring
        // intruded faces for each unstable face.
        borrowing_x.inds.resize(8*ncells);
        borrowing_x.neigh_faces.resize(8*ncells);
        borrowing_x.area.resize(8*ncells);
    }

    idim = 1;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        auto &borrowing_y = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_y.inds_pointer.resize(box);
        borrowing_y.size.resize(box);
        borrowing_y.size.setVal<amrex::RunOn::Device>(0);
        amrex::Long ncells = box.numPts();
        borrowing_y.inds.resize(8*ncells);
        borrowing_y.neigh_faces.resize(8*ncells);
        borrowing_y.area.resize(8*ncells);
    }

    idim = 2;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        auto &borrowing_z = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_z.inds_pointer.resize(box);
        borrowing_z.size.resize(box);
        borrowing_z.size.setVal<amrex::RunOn::Device>(0);
        amrex::Long ncells = box.numPts();
        borrowing_z.inds.resize(8*ncells);
        borrowing_z.neigh_faces.resize(8*ncells);
        borrowing_z.area.resize(8*ncells);
    }
}


AMREX_GPU_DEVICE
int
ComputeNBorrowOneFaceExtension(const amrex::Dim3 cell, const amrex::Real S_ext,
                                      const amrex::Array4<amrex::Real>& S_red,
                                      const amrex::Array4<int>& flag_info_face,
                                      const amrex::Array4<int>& flag_ext_face, const int idim) {
    const int i = cell.x;
    const int j = cell.y;
    const int k = cell.z;
    int n_borrow = 0;
    bool stop = false;
    if(idim == 0){
        for (int j_n = -1; j_n < 2; j_n++) {
            for (int k_n = -1; k_n < 2; k_n++) {
                //This if makes sure that we don't visit the "diagonal neighbours"
                if (!(j_n == k_n || j_n == -k_n)) {
                    // Here a face is available if it doesn't need to be extended itself and if its
                    // area exceeds Sz_ext. Here we need to take into account if the intruded face
                    // has given away already some area, so we use Sz_red rather than Sz.
                    // If no face is available we don't do anything and we will need to use the
                    // multi-face extensions.
                    if (S_red(i, j + j_n, k + k_n) > S_ext
                        && (flag_info_face(i, j + j_n, k + k_n) == 1 ||
                             flag_info_face(i, j + j_n, k + k_n) == 2)
                        && flag_ext_face(i, j, k) && ! stop) {
                        n_borrow += 1;
                        stop = true;
                    }
                }
            }
        }
    }else if(idim == 1){
        for (int i_n = -1; i_n < 2; i_n++) {
            for (int k_n = -1; k_n < 2; k_n++) {
                //This if makes sure that we don't visit the "diagonal neighbours"
                if (!(i_n == k_n || i_n == -k_n)) {
                    // Here a face is available if it doesn't need to be extended itself and if its
                    // area exceeds Sz_ext. Here we need to take into account if the intruded face
                    // has given away already some area, so we use Sz_red rather than Sz.
                    // If no face is available we don't do anything and we will need to use the
                    // multi-face extensions.
                    if (S_red(i + i_n, j, k + k_n) > S_ext
                        && (flag_info_face(i + i_n, j, k + k_n) == 1
                             || flag_info_face(i + i_n, j, k + k_n) == 2)
                        && flag_ext_face(i, j, k) && ! stop) {
                        n_borrow += 1;
                        stop = true;
                    }
                }
            }
        }
    }else if(idim == 2) {
        for (int i_n = -1; i_n < 2; i_n++) {
            for (int j_n = -1; j_n < 2; j_n++) {
                //This if makes sure that we don't visit the "diagonal neighbours"
                if (!(i_n == j_n || i_n == -j_n)) {
                    // Here a face is available if it doesn't need to be extended itself and if its
                    // area exceeds Sz_ext. Here we need to take into account if the intruded face
                    // has given away already some area, so we use Sz_red rather than Sz.
                    // If no face is available we don't do anything and we will need to use the
                    // multi-face extensions.
                    if (S_red(i + i_n, j + j_n, k) > S_ext
                        && (flag_info_face(i + i_n, j + j_n, k) == 1
                             || flag_info_face(i + i_n, j + j_n, k) == 2)
                        && flag_ext_face(i, j, k) && ! stop) {
                        n_borrow += 1;
                        stop = true;
                    }
                }
            }
        }
    }

    return n_borrow;
}


AMREX_GPU_DEVICE
int
ComputeNBorrowEightFacesExtension(const amrex::Dim3 cell, const amrex::Real S_ext,
                                         const amrex::Array4<amrex::Real>& S_red,
                                         const amrex::Array4<amrex::Real>& S,
                                         const amrex::Array4<int>& flag_info_face, const int idim) {
    const int i = cell.x;
    const int j = cell.y;
    const int k = cell.z;
    int n_borrow = 0;
    amrex::Array2D<int, 0, 2, 0, 2> local_avail{};

    if(idim == 0) {
        // Use a local 3x3 array to keep track of which of the neighboring faces is available to be
        // intruded.
        for (int j_loc = 0; j_loc <= 2; j_loc++) {
            for (int k_loc = 0; k_loc <= 2; k_loc++) {
                local_avail(j_loc, k_loc) = (flag_info_face(i, j + j_loc - 1, k + k_loc - 1) == 1
                                             || flag_info_face(i, j + j_loc - 1, k + k_loc - 1) == 2);
            }
        }

        // denom is the total area that can be borrowed from the neighboring cells
        amrex::Real denom = local_avail(0, 1) * S(i, j - 1, k) +
            local_avail(2, 1) * S(i, j + 1, k) +
            local_avail(1, 0) * S(i, j, k - 1) +
            local_avail(1, 2) * S(i, j, k + 1) +
            local_avail(0, 0) * S(i, j - 1, k - 1) +
            local_avail(2, 0) * S(i, j + 1, k - 1) +
            local_avail(0, 2) * S(i, j - 1, k + 1) +
            local_avail(2, 2) * S(i, j + 1, k + 1);

        // This flag will be used to check if any of the intruded faces would give away too much
        // area and needs to be excluded from the extension process. It is initialized to true to
        // make sure that we perform at least one iteration of the while loop.
        bool neg_face = true;

        // If the total area of the neighboring available faces is greater than the area needed for
        // the extension, then we can try to do the extension. (denom >= S_ext)
        // If any of the surrounding faces would be giving away to much area, then we have to
        // exclude it from the extension process and repeat. We keep trying until either all the faces
        // are big enough to be intruded or until we find out that none of the faces can be intruded
        // (that is denom = 0).
        while (denom >= S_ext && neg_face && denom > 0) {
            neg_face = false;
            for (int j_n = -1; j_n < 2; j_n++) {
                for (int k_n = -1; k_n < 2; k_n++) {
                    if (local_avail(j_n + 1, k_n + 1)) {
                        amrex::Real patch = S_ext * S(i, j + j_n, k + k_n) / denom;
                        if (S_red(i, j + j_n, k + k_n) - patch <= 0) {
                            neg_face = true;
                            // Whenever a face cannot be intruded we set the corresponding entry
                            // in local_avail to False
                            local_avail(j_n + 1, k_n + 1) = false;
                        }
                    }
                }
            }

            denom = local_avail(0, 1) * S(i, j - 1, k) +
                local_avail(2, 1) * S(i, j + 1, k) +
                local_avail(1, 0) * S(i, j, k - 1) +
                local_avail(1, 2) * S(i, j, k + 1) +
                local_avail(0, 0) * S(i, j - 1, k - 1) +
                local_avail(2, 0) * S(i, j + 1, k - 1) +
                local_avail(0, 2) * S(i, j - 1, k + 1) +
                local_avail(2, 2) * S(i, j + 1, k + 1);
        }
    }else if(idim == 1) {
        for(int i_loc = 0; i_loc <= 2; i_loc++){
            for(int k_loc = 0; k_loc <= 2; k_loc++){
                local_avail(i_loc, k_loc) = (flag_info_face(i + i_loc - 1, j, k + k_loc - 1) == 1
                                             || flag_info_face(i + i_loc - 1, j, k + k_loc - 1) == 2);
            }
        }

        amrex::Real denom = local_avail(0, 1) * S(i - 1, j, k) +
            local_avail(2, 1) * S(i + 1, j, k) +
            local_avail(1, 0) * S(i, j, k - 1) +
            local_avail(1, 2) * S(i, j, k + 1) +
            local_avail(0, 0) * S(i - 1, j, k - 1) +
            local_avail(2, 0) * S(i + 1, j, k - 1) +
            local_avail(0, 2) * S(i - 1, j, k + 1) +
            local_avail(2, 2) * S(i + 1, j, k + 1);

        bool neg_face = true;

        while(denom >= S_ext && neg_face && denom > 0){
            neg_face = false;
            for (int i_n = -1; i_n < 2; i_n++) {
                for (int k_n = -1; k_n < 2; k_n++) {
                    if(local_avail(i_n + 1, k_n + 1)){
                        amrex::Real patch = S_ext * S(i + i_n, j, k + k_n) / denom;
                        if(S_red(i + i_n, j, k + k_n) - patch <= 0) {
                            neg_face = true;
                            local_avail(i_n + 1, k_n + 1) = false;
                        }
                    }
                }
            }

            denom = local_avail(0, 1) * S(i - 1, j, k) +
                local_avail(2, 1) * S(i + 1, j, k) +
                local_avail(1, 0) * S(i, j, k - 1) +
                local_avail(1, 2) * S(i, j, k + 1) +
                local_avail(0, 0) * S(i - 1, j, k - 1) +
                local_avail(2, 0) * S(i + 1, j, k - 1) +
                local_avail(0, 2) * S(i - 1, j, k + 1) +
                local_avail(2, 2) * S(i + 1, j, k + 1);
        }

    } else if(idim == 2){
        for(int i_loc = 0; i_loc <= 2; i_loc++){
            for(int j_loc = 0; j_loc <= 2; j_loc++){
                local_avail(i_loc, j_loc) = (flag_info_face(i + i_loc - 1, j + j_loc - 1, k) == 1
                                             || flag_info_face(i + i_loc - 1, j + j_loc - 1, k) == 2);
            }
        }

        amrex::Real denom = local_avail(0, 1) * S(i - 1, j, k) +
            local_avail(2, 1) * S(i + 1, j, k) +
            local_avail(1, 0) * S(i, j - 1, k) +
            local_avail(1, 2) * S(i, j + 1, k) +
            local_avail(0, 0) * S(i - 1, j - 1, k) +
            local_avail(2, 0) * S(i + 1, j - 1, k) +
            local_avail(0, 2) * S(i - 1, j + 1, k) +
            local_avail(2, 2) * S(i + 1, j + 1, k);

        bool neg_face = true;

        while(denom >= S_ext && neg_face && denom > 0){
            neg_face = false;
            for (int i_n = -1; i_n < 2; i_n++) {
                for (int j_n = -1; j_n < 2; j_n++) {
                    if(local_avail(i_n + 1, j_n + 1)){
                        amrex::Real patch = S_ext * S(i + i_n, j + j_n, k) / denom;
                        if(S_red(i + i_n, j + j_n, k) - patch <= 0) {
                            neg_face = true;
                            local_avail(i_n + 1, j_n + 1) = false;
                        }
                    }
                }
            }

            denom = local_avail(0, 1) * S(i - 1, j, k) +
                local_avail(2, 1) * S(i + 1, j, k) +
                local_avail(1, 0) * S(i, j - 1, k) +
                local_avail(1, 2) * S(i, j + 1, k) +
                local_avail(0, 0) * S(i - 1, j - 1, k) +
                local_avail(2, 0) * S(i + 1, j - 1, k) +
                local_avail(0, 2) * S(i - 1, j + 1, k) +
                local_avail(2, 2) * S(i + 1, j + 1, k);
        }

    }

    // We count the number of entries in local_avail which are still True, this is the number of
    // neighboring faces which are intruded
    for(int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
            n_borrow += local_avail(ii, jj);
        }
    }

    return n_borrow;
}


void
WarpX::ComputeOneWayExtensions() {
#ifdef AMREX_USE_EB
    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &cell_size = CellSize(maxLevel());

    // Do the extensions in the x-plane
    int idim = 0;

    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sx = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_x = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_info_face_x = m_flag_info_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_x = (*m_borrowing[maxLevel()][idim])[mfi];
        auto const &borrowing_x_inds_pointer = borrowing_x.inds_pointer.array();
        auto const &borrowing_x_size = borrowing_x.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_x_inds = borrowing_x.inds.data();
        FaceInfoBox::Neighbours* borrowing_x_neigh_faces = borrowing_x.neigh_faces.data();
        amrex::Real* borrowing_x_area = borrowing_x.area.data();
        int& vecs_size_x = borrowing_x.vecs_size;

        auto const &Sx_mod = m_area_mod[maxLevel()][idim]->array(mfi);
        const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
        const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
        const amrex::Real dy = cell_size[1];
        const amrex::Real dz = cell_size[2];

        vecs_size_x = amrex::Scan::PrefixSum<int>(ncells,
        [=] AMREX_GPU_DEVICE (int icell) {
            const amrex::Dim3 cell = box.atOffset(icell).dim3();
            const int i = cell.x;
            const int j = cell.y;
            const int k = cell.z;
            // If the face doesn't need to be extended break the loop
            if (!flag_ext_face_x(i, j, k)) {
                return 0;
            }

            amrex::Real Sx_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                  lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
            amrex::Real Sx_ext = Sx_stab - Sx(i, j, k);
            int n_borrow = ComputeNBorrowOneFaceExtension(cell, Sx_ext, Sx_mod, flag_info_face_x,
                                                          flag_ext_face_x, idim);

            borrowing_x_size(i, j, k) = n_borrow;
            return n_borrow;
        },
        [=] AMREX_GPU_DEVICE (int icell, int ps){
            amrex::Dim3 cell = box.atOffset(icell).dim3();
            const int i = cell.x;
            const int j = cell.y;
            const int k = cell.z;
            int nborrow = borrowing_x_size(i, j, k);
            if (nborrow == 0) {
                borrowing_x_inds_pointer(i, j, k) = nullptr;
            } else{
                borrowing_x_inds_pointer(i, j, k) = borrowing_x_inds + ps;

                amrex::Real Sx_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                      lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
                amrex::Real Sx_ext = Sx_stab - Sx(i, j, k);
                for (int j_n = -1; j_n < 2; j_n++) {
                    for (int k_n = -1; k_n < 2; k_n++) {
                        //This if makes sure that we don't visit the "diagonal neighbours"
                        if( !(j_n == k_n || j_n == -k_n)){
                            // Here a face is available if it doesn't need to be extended itself and if its
                            // area exceeds Sz_ext. Here we need to take into account if the intruded face
                            // has given away already some area, so we use Sz_red rather than Sz.
                            // If no face is available we don't do anything and we will need to use the
                            // multi-face extensions.
                            if (Sx_mod(i, j + j_n, k + k_n) > Sx_ext
                                && ( flag_info_face_x(i, j + j_n, k + k_n) == 1
                                      || flag_info_face_x(i, j + j_n, k + k_n) == 2)
                                && flag_ext_face_x(i, j, k)) {
                                Sx_mod(i, j + j_n, k + k_n) -= Sx_ext;
                                // Insert the index of the face info
                                borrowing_x_inds[ps] = ps;
                                // Store the information about the intruded face in the dataset of the
                                // faces which are borrowing area
                                FaceInfoBox::addConnectedNeighbor(j_n, k_n, ps,
                                                                 borrowing_x_neigh_faces);

                                borrowing_x_area[ps] = Sx_ext;

                                flag_info_face_x(i, j + j_n, k + k_n) = 2;
                                // Add the area to the intruding face.
                                Sx_mod(i, j, k) += Sx_ext;
                                flag_ext_face_x(i, j, k) = false;
                            }
                        }
                    }
                }
            }
        }, amrex::Scan::Type::exclusive);


    }
    // Do the extensions in the y-plane
    idim = 1;

    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sy = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_y = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_info_face_y = m_flag_info_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_y = (*m_borrowing[maxLevel()][idim])[mfi];
        auto const &borrowing_y_inds_pointer = borrowing_y.inds_pointer.array();
        auto const &borrowing_y_size = borrowing_y.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_y_inds = borrowing_y.inds.data();
        FaceInfoBox::Neighbours* borrowing_y_neigh_faces = borrowing_y.neigh_faces.data();
        amrex::Real* borrowing_y_area = borrowing_y.area.data();
        int& vecs_size_y = borrowing_y.vecs_size;

        auto const &Sy_mod = m_area_mod[maxLevel()][idim]->array(mfi);
        const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
        const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
        const amrex::Real dx = cell_size[0];
        const amrex::Real dz = cell_size[2];

        vecs_size_y = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell) {
                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_y(i, j, k)) {
                    return 0;
                }

                amrex::Real Sy_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k + 1) * dz,
                                                      lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
                amrex::Real Sy_ext = Sy_stab - Sy(i, j, k);
                int n_borrow =
                    ComputeNBorrowOneFaceExtension(cell, Sy_ext, Sy_mod, flag_info_face_y,
                                                   flag_ext_face_y, idim);

                borrowing_y_size(i, j, k) = n_borrow;
                return n_borrow;
                },
            [=] AMREX_GPU_DEVICE (int icell, int ps){
                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;
                const int nborrow = borrowing_y_size(i, j, k);

                if (nborrow == 0) {
                    borrowing_y_inds_pointer(i, j, k) = nullptr;
                } else{
                    borrowing_y_inds_pointer(i, j, k) = borrowing_y_inds + ps;

                    amrex::Real Sy_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k+1) * dz,
                                                          lz(i, j, k) * dx, lz(i+1, j, k) * dz});
                    amrex::Real Sy_ext = Sy_stab - Sy(i, j, k);
                    for (int i_n = -1; i_n < 2; i_n++) {
                        for (int k_n = -1; k_n < 2; k_n++) {
                            //This if makes sure that we don't visit the "diagonal neighbours"
                            if( !(i_n == k_n || i_n == -k_n)){
                                // Here a face is available if it doesn't need to be extended itself and if its
                                // area exceeds Sz_ext. Here we need to take into account if the intruded face
                                // has given away already some area, so we use Sz_red rather than Sz.
                                // If no face is available we don't do anything and we will need to use the
                                // multi-face extensions.
                                if (Sy_mod(i + i_n, j, k + k_n) > Sy_ext
                                    && (flag_info_face_y(i + i_n, j, k + k_n) == 1
                                         || flag_info_face_y(i + i_n, j, k + k_n) == 2)
                                    && flag_ext_face_y(i, j, k)) {
                                    Sy_mod(i + i_n, j, k + k_n) -= Sy_ext;
                                    // Insert the index of the face info
                                    borrowing_y_inds[ps] = ps;
                                    // Store the information about the intruded face in the dataset of the
                                    // faces which are borrowing area
                                    FaceInfoBox::addConnectedNeighbor(i_n, k_n, ps,
                                                                      borrowing_y_neigh_faces);
                                    borrowing_y_area[ps] = Sy_ext;

                                    flag_info_face_y(i + i_n, j, k + k_n) = 2;
                                    // Add the area to the intruding face.
                                    Sy_mod(i, j, k) = Sy(i, j, k) + Sy_ext;
                                    flag_ext_face_y(i, j, k) = false;
                                }
                            }
                        }
                    }
                }
            }, amrex::Scan::Type::exclusive);

    }
    // Do the extensions in the z-plane
    idim = 2;

    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sz = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_z = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_info_face_z = m_flag_info_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_z = (*m_borrowing[maxLevel()][idim])[mfi];
        auto const &borrowing_z_inds_pointer = borrowing_z.inds_pointer.array();
        auto const &borrowing_z_size = borrowing_z.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_z_inds = borrowing_z.inds.data();
        FaceInfoBox::Neighbours* borrowing_z_neigh_faces = borrowing_z.neigh_faces.data();
        amrex::Real* borrowing_z_area = borrowing_z.area.data();
        int& vecs_size_z = borrowing_z.vecs_size;

        auto const &Sz_mod = m_area_mod[maxLevel()][idim]->array(mfi);
        const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
        const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
        const amrex::Real dx = cell_size[0];
        const amrex::Real dy = cell_size[1];

        vecs_size_z = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell) {
                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_z(i, j, k)) {
                    return 0;
                }

                amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                      ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
                amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);
                const int n_borrow =
                    ComputeNBorrowOneFaceExtension(cell, Sz_ext, Sz_mod, flag_info_face_z,
                                                   flag_ext_face_z, idim);


                borrowing_z_size(i, j, k) = n_borrow;
                return n_borrow;
            },
            [=] AMREX_GPU_DEVICE (int icell, int ps){
                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;
                const int nborrow = borrowing_z_size(i, j, k);
                if (nborrow == 0) {
                    borrowing_z_inds_pointer(i, j, k) = nullptr;
                } else{
                    borrowing_z_inds_pointer(i, j, k) = borrowing_z_inds + ps;

                    amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                          ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
                    amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);
                    for (int i_n = -1; i_n < 2; i_n++) {
                        for (int j_n = -1; j_n < 2; j_n++) {
                            //This if makes sure that we don't visit the "diagonal neighbours"
                            if( !(i_n == j_n || i_n == -j_n)){
                                // Here a face is available if it doesn't need to be extended itself and if its
                                // area exceeds Sz_ext. Here we need to take into account if the intruded face
                                // has given away already some area, so we use Sz_red rather than Sz.
                                // If no face is available we don't do anything and we will need to use the
                                // multi-face extensions.
                                if (Sz_mod(i + i_n, j + j_n, k) > Sz_ext
                                    && (flag_info_face_z(i + i_n, j + j_n, k) == 1
                                         || flag_info_face_z(i + i_n, j + j_n, k) == 2)
                                    && flag_ext_face_z(i, j, k)) {
                                    Sz_mod(i + i_n, j + j_n, k) -= Sz_ext;
                                    // Insert the index of the face info
                                    borrowing_z_inds[ps] = ps;
                                    // Store the information about the intruded face in the dataset of the
                                    // faces which are borrowing area
                                    FaceInfoBox::addConnectedNeighbor(i_n, j_n, ps,
                                                                      borrowing_z_neigh_faces);
                                    borrowing_z_area[ps] = Sz_ext;

                                    flag_info_face_z(i + i_n, j + j_n, k) = 2;
                                    // Add the area to the intruding face.
                                    Sz_mod(i, j, k) = Sz(i, j, k) + Sz_ext;
                                    flag_ext_face_z(i, j, k) = false;
                                }
                            }
                        }
                    }
                }
            },
        amrex::Scan::Type::exclusive);

    }

#endif
}


void
WarpX::ComputeEightWaysExtensions() {
#ifdef AMREX_USE_EB

    auto const &cell_size = CellSize(maxLevel());
    int idim = 0;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sx = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_x = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_info_face_x = m_flag_info_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_x = (*m_borrowing[maxLevel()][idim])[mfi];
        auto const &borrowing_x_inds_pointer = borrowing_x.inds_pointer.array();
        auto const &borrowing_x_size = borrowing_x.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_x_inds = borrowing_x.inds.data();
        FaceInfoBox::Neighbours* borrowing_x_neigh_faces = borrowing_x.neigh_faces.data();
        amrex::Real* borrowing_x_area = borrowing_x.area.data();
        int& vecs_size_x = borrowing_x.vecs_size;

        auto const &Sx_mod = m_area_mod[maxLevel()][idim]->array(mfi);
        const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
        const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
        const amrex::Real dy = cell_size[1];
        const amrex::Real dz = cell_size[2];

        vecs_size_x += amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell){
                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_x(i, j, k)) {
                    return 0;
                }
                amrex::Real Sx_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                      lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
                amrex::Real Sx_ext = Sx_stab - Sx(i, j, k);
                const int n_borrow = ComputeNBorrowEightFacesExtension(cell, Sx_ext, Sx_mod, Sx,
                                                                 flag_info_face_x, idim);

                borrowing_x_size(i, j, k) = n_borrow;
                return n_borrow;
                },
            [=] AMREX_GPU_DEVICE (int icell, int ps){
                ps += vecs_size_x;

                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;

                if (!flag_ext_face_x(i, j, k)) {
                    return;
                }

                const int nborrow = borrowing_x_size(i, j, k);
                if (nborrow == 0) {
                    borrowing_x_inds_pointer(i, j, k) = nullptr;
                } else{
                    borrowing_x_inds_pointer(i, j, k) = borrowing_x_inds + ps;
                    Sx_mod(i, j, k) = Sx(i, j, k);
                    amrex::Real Sx_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                  lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
                    amrex::Real Sx_ext = Sx_stab - Sx(i, j, k);
                    amrex::Array2D<amrex::Real, 0, 2, 0, 2> local_avail{};
                    for(int j_loc = 0; j_loc <= 2; j_loc++){
                        for(int k_loc = 0; k_loc <= 2; k_loc++){
                            local_avail(j_loc, k_loc) = (flag_info_face_x(i, j  + j_loc - 1, k + k_loc - 1) == 1
                                                         || flag_info_face_x(i, j  + j_loc - 1, k + k_loc - 1) == 2);
                        }
                    }

                    amrex::Real denom = local_avail(0, 1) * Sx(i, j - 1, k) +
                                        local_avail(2, 1) * Sx(i, j + 1, k) +
                                        local_avail(1, 0) * Sx(i, j, k - 1) +
                                        local_avail(1, 2) * Sx(i, j, k + 1) +
                                        local_avail(0, 0) * Sx(i, j - 1, k - 1) +
                                        local_avail(2, 0) * Sx(i, j + 1, k - 1) +
                                        local_avail(0, 2) * Sx(i, j - 1, k + 1) +
                                        local_avail(2, 2) * Sx(i, j + 1, k + 1);

                    bool neg_face = true;

                    while(denom >= Sx_ext && neg_face && denom > 0){
                        neg_face = false;
                        for (int j_n = -1; j_n < 2; j_n++) {
                            for (int k_n = -1; k_n < 2; k_n++) {
                                if(local_avail(j_n + 1, k_n + 1)){
                                    amrex::Real patch = Sx_ext * Sx(i, j + j_n, k+ k_n) / denom;
                                    if(Sx_mod(i, j + j_n, k + k_n) - patch <= 0) {
                                        neg_face = true;
                                        local_avail(j_n + 1, k_n + 1) = false;
                                    }
                                }
                            }
                        }

                        denom = local_avail(0, 1) * Sx(i, j - 1, k) +
                                local_avail(2, 1) * Sx(i, j + 1, k) +
                                local_avail(1, 0) * Sx(i, j, k - 1) +
                                local_avail(1, 2) * Sx(i, j, k + 1) +
                                local_avail(0, 0) * Sx(i, j - 1, k - 1) +
                                local_avail(2, 0) * Sx(i, j + 1, k - 1) +
                                local_avail(0, 2) * Sx(i, j - 1, k + 1) +
                                local_avail(2, 2) * Sx(i, j + 1, k + 1);
                    }

                    if(denom >= Sx_ext){
                        Sx_mod(i, j, k) = Sx(i, j, k);
                        int count = 0;
                        for (int j_n = -1; j_n < 2; j_n++) {
                            for (int k_n = -1; k_n < 2; k_n++) {
                                if(local_avail(j_n + 1, k_n + 1)){
                                    amrex::Real patch = Sx_ext * Sx(i, j + j_n, k + k_n) / denom;
                                    borrowing_x_inds[ps + count] = ps + count;
                                    // Store the information about the intruded face in the dataset of the
                                    // faces which are borrowing area
                                    FaceInfoBox::addConnectedNeighbor(j_n, k_n, ps + count,
                                                                      borrowing_x_neigh_faces);
                                    borrowing_x_area[ps + count] = patch;

                                    flag_info_face_x(i, j + j_n, k + k_n) = 2;
                                    Sx_mod(i, j, k) += patch;
                                    Sx_mod(i, j + j_n, k + k_n) -= patch;
                                    count+=1;
                                }
                            }
                        }
                        flag_ext_face_x(i, j, k) = false;
                    }
                }
            },
            amrex::Scan::Type::exclusive);

    }

    idim = 1;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {

        // Do the extensions in the y-plane
        amrex::Box const &box = mfi.validbox();

        auto const &Sy = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_y = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_info_face_y = m_flag_info_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_y = (*m_borrowing[maxLevel()][idim])[mfi];
        auto const &borrowing_y_inds_pointer = borrowing_y.inds_pointer.array();
        auto const &borrowing_y_size = borrowing_y.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_y_inds = borrowing_y.inds.data();
        FaceInfoBox::Neighbours* borrowing_y_neigh_faces = borrowing_y.neigh_faces.data();
        amrex::Real* borrowing_y_area = borrowing_y.area.data();
        int& vecs_size_y = borrowing_y.vecs_size;

        auto const &Sy_mod = m_area_mod[maxLevel()][idim]->array(mfi);
        const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
        const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dz = cell_size[2];

        vecs_size_y += amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell){
                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_y(i, j, k)){
                    return 0;
                }
                amrex::Real Sy_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k + 1) * dz,
                                                    lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
                amrex::Real Sy_ext = Sy_stab - Sy(i, j, k);
                const int n_borrow = ComputeNBorrowEightFacesExtension(cell, Sy_ext, Sy_mod, Sy,
                                                                 flag_info_face_y, idim);

                borrowing_y_size(i, j, k) = n_borrow;
                return n_borrow;
            },
            [=] AMREX_GPU_DEVICE (int icell, int ps) {

                ps += vecs_size_y;

                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;

                if (!flag_ext_face_y(i, j, k)) {
                    return;
                }

                int nborrow = borrowing_y_size(i, j, k);
                if (nborrow == 0) {
                    borrowing_y_inds_pointer(i, j, k) = nullptr;
                } else {
                    borrowing_y_inds_pointer(i, j, k) = borrowing_y_inds + ps;

                    Sy_mod(i, j, k) = Sy(i, j, k);
                    amrex::Real Sy_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k + 1) * dz,
                                                          lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
                    amrex::Real Sy_ext = Sy_stab - Sy(i, j, k);
                    amrex::Array2D<amrex::Real, 0, 2, 0, 2> local_avail{};
                    for (int i_loc = 0; i_loc <= 2; i_loc++) {
                        for (int k_loc = 0; k_loc <= 2; k_loc++) {
                            local_avail(i_loc, k_loc) = (flag_info_face_y(i + i_loc - 1, j, k + k_loc - 1) == 1
                                                         || flag_info_face_y(i + i_loc - 1, j, k + k_loc - 1) == 2)   ;
                        }
                    }

                    amrex::Real denom = local_avail(0, 1) * Sy(i - 1, j, k) +
                        local_avail(2, 1) * Sy(i + 1, j, k) +
                        local_avail(1, 0) * Sy(i, j, k - 1) +
                        local_avail(1, 2) * Sy(i, j, k + 1) +
                        local_avail(0, 0) * Sy(i - 1, j, k - 1) +
                        local_avail(2, 0) * Sy(i + 1, j, k - 1) +
                        local_avail(0, 2) * Sy(i - 1, j, k + 1) +
                        local_avail(2, 2) * Sy(i + 1, j, k + 1);

                    bool neg_face = true;

                    while (denom >= Sy_ext && neg_face && denom > 0) {
                        neg_face = false;
                        for (int i_n = -1; i_n < 2; i_n++) {
                            for (int k_n = -1; k_n < 2; k_n++) {
                                if (local_avail(i_n + 1, k_n + 1)) {
                                    amrex::Real patch = Sy_ext * Sy(i + i_n, j, k + k_n) / denom;
                                    if (Sy_mod(i + i_n, j, k + k_n) - patch <= 0) {
                                        neg_face = true;
                                        local_avail(i_n + 1, k_n + 1) = false;
                                    }
                                }
                            }
                        }

                        denom = local_avail(0, 1) * Sy(i - 1, j, k) +
                            local_avail(2, 1) * Sy(i + 1, j, k) +
                            local_avail(1, 0) * Sy(i, j, k - 1) +
                            local_avail(1, 2) * Sy(i, j, k + 1) +
                            local_avail(0, 0) * Sy(i - 1, j, k - 1) +
                            local_avail(2, 0) * Sy(i + 1, j, k - 1) +
                            local_avail(0, 2) * Sy(i - 1, j, k + 1) +
                            local_avail(2, 2) * Sy(i + 1, j, k + 1);
                    }

                    if (denom >= Sy_ext) {
                        Sy_mod(i, j, k) = Sy(i, j, k);
                        int count = 0;
                        for (int i_n = -1; i_n < 2; i_n++) {
                            for (int k_n = -1; k_n < 2; k_n++) {
                                if (local_avail(i_n + 1, k_n + 1)) {
                                    amrex::Real patch = Sy_ext * Sy(i + i_n, j, k + k_n) / denom;
                                    borrowing_y_inds[ps + count] = ps + count;
                                    FaceInfoBox::addConnectedNeighbor(i_n, k_n, ps + count,
                                                                      borrowing_y_neigh_faces);
                                    borrowing_y_area[ps + count] = patch;

                                    flag_info_face_y(i + i_n, j, k + k_n) = 2;
                                    Sy_mod(i, j, k) += patch;
                                    Sy_mod(i + i_n, j, k + k_n) -= patch;
                                    count += 1;
                                }
                            }
                        }
                        flag_ext_face_y(i, j, k) = false;
                    }
                }
            },
            amrex::Scan::Type::exclusive);

    }

    // Do the extensions in the z-plane
    idim = 2;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sz = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_z = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_info_face_z = m_flag_info_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_z = (*m_borrowing[maxLevel()][idim])[mfi];
        auto const &borrowing_z_inds_pointer = borrowing_z.inds_pointer.array();
        auto const &borrowing_z_size = borrowing_z.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_z_inds = borrowing_z.inds.data();
        FaceInfoBox::Neighbours* borrowing_z_neigh_faces = borrowing_z.neigh_faces.data();
        amrex::Real* borrowing_z_area = borrowing_z.area.data();
        int& vecs_size_z = borrowing_z.vecs_size;

        auto const &Sz_mod = m_area_mod[maxLevel()][idim]->array(mfi);
        const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
        const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
        const amrex::Real dx = cell_size[0];
        const amrex::Real dy = cell_size[1];

        vecs_size_z += amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell){
                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_z(i, j, k)) {
                    return 0;
                }
                amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                    ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
                amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);
                const int n_borrow = ComputeNBorrowEightFacesExtension(cell, Sz_ext, Sz_mod, Sz,
                                                                 flag_info_face_z, idim);

                borrowing_z_size(i, j, k) = n_borrow;
                return n_borrow;
            },
            [=] AMREX_GPU_DEVICE (int icell, int ps) {

                ps += vecs_size_z;

                const amrex::Dim3 cell = box.atOffset(icell).dim3();
                const int i = cell.x;
                const int j = cell.y;
                const int k = cell.z;

                if (!flag_ext_face_z(i, j, k)) {
                    return;
                }

                const int nborrow = borrowing_z_size(i, j, k);
                if (nborrow == 0) {
                    borrowing_z_inds_pointer(i, j, k) = nullptr;
                } else {
                    borrowing_z_inds_pointer(i, j, k) = borrowing_z_inds + ps;

                    Sz_mod(i, j, k) = Sz(i, j, k);
                    amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                          ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
                    amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);
                    amrex::Array2D<amrex::Real, 0, 2, 0, 2> local_avail{};
                    for(int i_loc = 0; i_loc <= 2; i_loc++){
                        for(int j_loc = 0; j_loc <= 2; j_loc++){
                            local_avail(i_loc, j_loc) = (flag_info_face_z(i + i_loc - 1, j + j_loc - 1, k) == 1
                                                         || flag_info_face_z(i + i_loc - 1, j + j_loc - 1, k) == 2);
                        }
                    }

                    amrex::Real denom = local_avail(0, 1) * Sz(i - 1, j, k) +
                                        local_avail(2, 1) * Sz(i + 1, j, k) +
                                        local_avail(1, 0) * Sz(i, j - 1, k) +
                                        local_avail(1, 2) * Sz(i, j + 1, k) +
                                        local_avail(0, 0) * Sz(i - 1, j - 1, k) +
                                        local_avail(2, 0) * Sz(i + 1, j - 1, k) +
                                        local_avail(0, 2) * Sz(i - 1, j + 1, k) +
                                        local_avail(2, 2) * Sz(i + 1, j + 1, k);

                    bool neg_face = true;

                    while(denom >= Sz_ext && neg_face && denom > 0){
                        neg_face = false;
                        for (int i_n = -1; i_n < 2; i_n++) {
                            for (int j_n = -1; j_n < 2; j_n++) {
                                if(local_avail(i_n + 1, j_n + 1)){
                                    amrex::Real patch = Sz_ext * Sz(i + i_n, j + j_n, k) / denom;
                                    if(Sz_mod(i + i_n, j + j_n, k) - patch <= 0) {
                                        neg_face = true;
                                        local_avail(i_n + 1, j_n + 1) = false;
                                    }
                                }
                            }
                        }

                        denom = local_avail(0, 1) * Sz(i - 1, j, k) +
                                local_avail(2, 1) * Sz(i + 1, j, k) +
                                local_avail(1, 0) * Sz(i, j - 1, k) +
                                local_avail(1, 2) * Sz(i, j + 1, k) +
                                local_avail(0, 0) * Sz(i - 1, j - 1, k) +
                                local_avail(2, 0) * Sz(i + 1, j - 1, k) +
                                local_avail(0, 2) * Sz(i - 1, j + 1, k) +
                                local_avail(2, 2) * Sz(i + 1, j + 1, k);
                    }

                    if(denom >= Sz_ext){
                        Sz_mod(i, j, k) = Sz(i, j, k);
                        int count = 0;
                        for (int i_n = -1; i_n < 2; i_n++) {
                            for (int j_n = -1; j_n < 2; j_n++) {
                                if(local_avail(i_n + 1, j_n + 1)){
                                    amrex::Real patch = Sz_ext * Sz(i + i_n, j + j_n, k) / denom;
                                    borrowing_z_inds[ps + count] = ps + count;
                                    FaceInfoBox::addConnectedNeighbor(i_n, j_n, ps + count,
                                                                      borrowing_z_neigh_faces);
                                    borrowing_z_area[ps + count] = patch;

                                    flag_info_face_z(i + i_n, j + j_n, k) = 2;

                                    Sz_mod(i, j, k) += patch;
                                    Sz_mod(i + i_n, j + j_n, k) -= patch;
                                    count +=1;
                                }
                            }
                        }
                        flag_ext_face_z(i, j, k) = false;
                    }
                }
            },
            amrex::Scan::Type::exclusive);

    }

#endif
}


void
WarpX::ShrinkBorrowing() {
    int idim = 0;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {
        auto &borrowing_x = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_x.inds.resize(borrowing_x.vecs_size);
        borrowing_x.neigh_faces.resize(borrowing_x.vecs_size);
        borrowing_x.area.resize(borrowing_x.vecs_size);
    }

    idim = 1;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {
        auto &borrowing_y = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_y.inds.resize(borrowing_y.vecs_size);
        borrowing_y.neigh_faces.resize(borrowing_y.vecs_size);
        borrowing_y.area.resize(borrowing_y.vecs_size);
    }

    idim = 2;
    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {
        auto &borrowing_z = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_z.inds.resize(borrowing_z.vecs_size);
        borrowing_z.neigh_faces.resize(borrowing_z.vecs_size);
        borrowing_z.area.resize(borrowing_z.vecs_size);
    }
}
