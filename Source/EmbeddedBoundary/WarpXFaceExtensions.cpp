#include "WarpX.H"
#include <AMReX_Scan.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>

/**
* \brief auxiliary function to count the amount of faces which still need to be extended
*/
amrex::Array1D<int, 0, 2>
WarpX::CountExtFaces(std::array<std::unique_ptr<amrex::iMultiFab>, 3>& flag_ext_face) {
    amrex::Array1D<int, 0, 2> sums{0, 0, 0};
#ifdef AMREX_USE_EB
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        int sum = 0;
        for (amrex::MFIter mfi(*flag_ext_face[idim]); mfi.isValid(); ++mfi) {
            amrex::Box const &box = mfi.validbox();
            auto const &flag_ext_face_dim = flag_ext_face[idim]->array(mfi);
            amrex::ParallelFor(box, [=, &sum](int i, int j, int k) {
                // Do I need to flag the reduction somehow?
                sum = sum + flag_ext_face_dim(i, j, k);
            });
        }
        sums(idim) = sum;
    }
#else
    amrex::ignore_unused(flag_ext_face);
#endif
    return sums;
}


/**
* \brief This is the main function computing the cell extension. Where possible it computes one-way
 *       extensions and, when this is not possible, it does eight-ways extensions.
*/
void
WarpX::ComputeFaceExtensions(std::array<std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
                             std::array<std::unique_ptr<amrex::MultiFab>, 3>& area_mod,
                             std::array<std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
                             std::array<std::unique_ptr<amrex::MultiFab>, 3> const& face_areas,
                             std::array<std::unique_ptr<amrex::iMultiFab>, 3>& flag_ext_face,
                             std::array<std::unique_ptr<amrex::iMultiFab>, 3>& flag_info_face,
                             std::array<std::unique_ptr<amrex::LayoutData<FaceInfoBox>>, 3>& borrowing,
                             int lev, bool flag_cp){
#ifdef AMREX_USE_EB
    amrex::Array1D<int, 0, 2> N_ext_faces = CountExtFaces(flag_ext_face);
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces(0));
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces(1));
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces(2));
    amrex::Print()<< "Faces to be extended in x:\t" << N_ext_faces(0) <<std::endl;
    amrex::Print()<< "Faces to be extended in y:\t" << N_ext_faces(1) <<std::endl;
    amrex::Print()<< "Faces to be extended in z:\t" << N_ext_faces(2) <<std::endl;

    InitBorrowing(Bfield, borrowing);
    amrex::Array1D<int, 0, 2> temp_inds = ComputeOneWayExtensions(Bfield, area_mod, edge_lengths,
                                                                  face_areas, flag_ext_face,
                                                                  flag_info_face, borrowing,
                                                                  lev, flag_cp);
    amrex::Array1D<int, 0, 2> N_ext_faces_after_one_way = CountExtFaces(flag_ext_face);
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces_after_one_way(0),
                                            amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces_after_one_way(1),
                                            amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces_after_one_way(2),
                                            amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print()<< "Faces to be extended after one way extension in x:\t" <<
                     N_ext_faces_after_one_way(0) <<std::endl;
    amrex::Print()<< "Faces to be extended after one way extension in y:\t" <<
                     N_ext_faces_after_one_way(1) <<std::endl;
    amrex::Print()<< "Faces to be extended after one way extension in z:\t" <<
                     N_ext_faces_after_one_way(2) <<std::endl;
    ComputeEightWaysExtensions(Bfield, area_mod, edge_lengths, face_areas, flag_ext_face,
                               flag_info_face, borrowing, temp_inds, lev, flag_cp);
    amrex::Array1D<int, 0, 2> N_ext_faces_after_eight_ways = CountExtFaces(flag_ext_face);
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces_after_eight_ways(0),
                                            amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces_after_eight_ways(1),
                                            amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::ParallelDescriptor::ReduceIntSum(N_ext_faces_after_eight_ways(2),
                                            amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print()<< "Faces to be extended after eight ways extension in x:\t" <<
                     N_ext_faces_after_eight_ways(0) <<std::endl;
    amrex::Print()<< "Faces to be extended after eight ways extension in y:\t" <<
                     N_ext_faces_after_eight_ways(1) <<std::endl;
    amrex::Print()<< "Faces to be extended after eight ways extension ins z:\t" <<
                     N_ext_faces_after_eight_ways(2) <<std::endl;

    if (N_ext_faces_after_eight_ways(0) > 0) {
        amrex::Abort("Some x faces could not be extended");
    }
    if (N_ext_faces_after_eight_ways(1) > 0) {
        amrex::Abort("Some y faces could not be extended");
    }
    if (N_ext_faces_after_eight_ways(2) > 0) {
        amrex::Abort("Some z faces could not be extended");
    }
#else
    amrex::ignore_unused(Bfield, area_mod, edge_lengths, face_areas, flag_ext_face, flag_info_face,
                         borrowing, lev, flag_cp);
#endif
}

/**
* \brief this function initializes the memory for the cell FaceInfoBoxes
*/
void
WarpX::InitBorrowing(std::array<std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
                     std::array<std::unique_ptr<amrex::LayoutData<FaceInfoBox>>, 3>& borrowing) {
#ifdef AMREX_USE_EB
    int idim = 0;
    for (amrex::MFIter mfi(*Bfield[0]); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        auto &borrowing_x = (*borrowing[idim])[mfi];
        borrowing_x.inds_pointer.resize(box);
        borrowing_x.size.resize(box);
        borrowing_x.size.setVal<amrex::RunOn::Host>(0);
        amrex::Long ncells = box.numPts();
        borrowing_x.inds.resize(8*ncells);
        borrowing_x.neigh_faces.resize(8*ncells);
        borrowing_x.area.resize(8*ncells);
    }

    idim = 1;
    for (amrex::MFIter mfi(*Bfield[idim]); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        auto &borrowing_y = (*borrowing[idim])[mfi];
        borrowing_y.inds_pointer.resize(box);
        borrowing_y.size.resize(box);
        borrowing_y.size.setVal<amrex::RunOn::Host>(0);
        amrex::Long ncells = box.numPts();
        borrowing_y.inds.resize(8*ncells);
        borrowing_y.neigh_faces.resize(8*ncells);
        borrowing_y.area.resize(8*ncells);
    }

    idim = 2;
    for (amrex::MFIter mfi(*Bfield[idim]); mfi.isValid(); ++mfi) {
        amrex::Box const &box = mfi.validbox();
        auto &borrowing_z = (*borrowing[idim])[mfi];
        borrowing_z.inds_pointer.resize(box);
        borrowing_z.size.resize(box);
        borrowing_z.size.setVal<amrex::RunOn::Host>(0);
        amrex::Long ncells = box.numPts();
        borrowing_z.inds.resize(8*ncells);
        borrowing_z.neigh_faces.resize(8*ncells);
        borrowing_z.area.resize(8*ncells);
    }
#else
    amrex::ignore_unused(Bfield, borrowing);
#endif
}

/**
* \brief For the face of cell pointing in direction idim, compute the number of faces
* we need to intrude. For the one-way extension this function returns only one or zero: one if the
* face can be extended withe the one-way extension, zeros if it can't.
*/
int
WarpX::ComputeNBorrowOneFaceExtension(amrex::Dim3 cell, amrex::Real S_ext,
                                      const amrex::Array4<amrex::Real>& S_red,
                                      const amrex::Array4<int>& flag_info_face,
                                      const amrex::Array4<int>& flag_ext_face, int idim) {
#ifdef AMREX_USE_EB
    int i = cell.x;
    int j = cell.y;
    int k = cell.z;
    int n_borrow = 0;
    bool stop = false;
    if(idim == 0){
        for (int j_n = -1; j_n < 2; j_n++) {
            for (int k_n = -1; k_n < 2; k_n++) {
                //This if makes sure that we don't visit the "diagonal neighbours"
                if (!(j_n == k_n or j_n == -k_n)) {
                    // Here a face is available if it doesn't need to be extended itself and if its
                    // area exceeds Sz_ext. Here we need to take into account if the intruded face
                    // has given away already some area, so we use Sz_red rather than Sz.
                    // If no face is available we don't do anything and we will need to use the
                    // multi-face extensions.
                    if (S_red(i, j + j_n, k + k_n) > S_ext
                        and (flag_info_face(i, j + j_n, k + k_n) == 1 or
                             flag_info_face(i, j + j_n, k + k_n) == 2)
                        and flag_ext_face(i, j, k) and not stop) {
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
                if (!(i_n == k_n or i_n == -k_n)) {
                    // Here a face is available if it doesn't need to be extended itself and if its
                    // area exceeds Sz_ext. Here we need to take into account if the intruded face
                    // has given away already some area, so we use Sz_red rather than Sz.
                    // If no face is available we don't do anything and we will need to use the
                    // multi-face extensions.
                    if (S_red(i + i_n, j, k + k_n) > S_ext
                        and (flag_info_face(i + i_n, j, k + k_n) == 1
                             or flag_info_face(i + i_n, j, k + k_n) == 2)
                        and flag_ext_face(i, j, k) and not stop) {
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
                if (!(i_n == j_n or i_n == -j_n)) {
                    // Here a face is available if it doesn't need to be extended itself and if its
                    // area exceeds Sz_ext. Here we need to take into account if the intruded face
                    // has given away already some area, so we use Sz_red rather than Sz.
                    // If no face is available we don't do anything and we will need to use the
                    // multi-face extensions.
                    if (S_red(i + i_n, j + j_n, k) > S_ext
                        and (flag_info_face(i + i_n, j + j_n, k) == 1
                             or flag_info_face(i + i_n, j + j_n, k) == 2)
                        and flag_ext_face(i, j, k) and not stop) {
                        n_borrow += 1;
                        stop = true;
                    }
                }
            }
        }
    }

    return n_borrow;
#else
    amrex::ignore_unused(cell, S_ext, S_red, flag_info_face, flag_ext_face, idim);
    return 0;
#endif
}

/**
* \brief For the face of cell pointing in direction idim, compute the number of faces
* we need to intrude. For the one-way extension this function returns only one or zero: one if the
* face can be extended withe the one-way extension, zeros if it can't.
*/
int
WarpX::ComputeNBorrowEightFacesExtension(amrex::Dim3 cell, amrex::Real S_ext,
                                         const amrex::Array4<amrex::Real>& S_red,
                                         const amrex::Array4<amrex::Real>& S,
                                         const amrex::Array4<int>& flag_info_face, int idim) {
#ifdef AMREX_USE_EB
    int i = cell.x;
    int j = cell.y;
    int k = cell.z;
    int n_borrow = 0;
    amrex::Array2D<amrex::Real, 0, 2, 0, 2> local_avail{};

    if(idim == 0) {
        for (int j_loc = 0; j_loc <= 2; j_loc++) {
            for (int k_loc = 0; k_loc <= 2; k_loc++) {
                local_avail(j_loc, k_loc) = (flag_info_face(i, j + j_loc - 1, k + k_loc - 1) == 1
                                             or flag_info_face(i, j + j_loc - 1, k + k_loc - 1) == 2);
            }
        }

        amrex::Real denom = local_avail(0, 1) * S(i, j - 1, k) +
            local_avail(2, 1) * S(i, j + 1, k) +
            local_avail(1, 0) * S(i, j, k - 1) +
            local_avail(1, 2) * S(i, j, k + 1) +
            local_avail(0, 0) * S(i, j - 1, k - 1) +
            local_avail(2, 0) * S(i, j + 1, k - 1) +
            local_avail(0, 2) * S(i, j - 1, k + 1) +
            local_avail(2, 2) * S(i, j + 1, k + 1);

        bool neg_face = true;

        while (denom >= S_ext and neg_face and denom > 0) {
            neg_face = false;
            for (int j_n = -1; j_n < 2; j_n++) {
                for (int k_n = -1; k_n < 2; k_n++) {
                    if (local_avail(j_n + 1, k_n + 1)) {
                        amrex::Real patch = S_ext * S(i, j + j_n, k + k_n) / denom;
                        if (S_red(i, j + j_n, k + k_n) - patch <= 0) {
                            neg_face = true;
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
                                             or flag_info_face(i + i_loc - 1, j, k + k_loc - 1) == 2);
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

        while(denom >= S_ext and neg_face and denom > 0){
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
                                             or flag_info_face(i + i_loc - 1, j + j_loc - 1, k) == 2);
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

        while(denom >= S_ext and neg_face and denom > 0){
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

    for(int ii = 0; ii < 3; ii++) {
        for (int jj = 0; jj < 3; jj++) {
            n_borrow += local_avail(ii, jj);
        }
    }

    return n_borrow;
#else
    amrex::ignore_unused(S_ext, S_red, S, flag_info_face, idim);
    return 0;
#endif
}

/**
* \brief Do the one-way extension
*/
amrex::Array1D<int, 0, 2>
WarpX::ComputeOneWayExtensions(std::array<std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
                               std::array<std::unique_ptr<amrex::MultiFab>, 3>& area_mod,
                               std::array<std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
                               std::array<std::unique_ptr<amrex::MultiFab>, 3> const& face_areas,
                               std::array<std::unique_ptr<amrex::iMultiFab>, 3>& flag_ext_face,
                               std::array<std::unique_ptr<amrex::iMultiFab>, 3>& flag_info_face,
                               std::array<std::unique_ptr<amrex::LayoutData<FaceInfoBox>>, 3>& borrowing,
                               int lev, bool flag_cp) {
#ifdef AMREX_USE_EB
    //This variable is equal to lev if this is a fine patch or
    // equal to lev -1 if this is a coarse patch
    int lev_loc = lev;
    if(flag_cp and lev > 0){
        lev_loc = lev -1;
    }

    auto const &cell_size = CellSize(lev_loc);

    int nelems_x;
    int nelems_y;
    int nelems_z;

    // Do the extensions in the x-plane
    int idim = 0;

    for (amrex::MFIter mfi(*Bfield[idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sx = face_areas[idim]->array(mfi);
        auto const &flag_ext_face_x = flag_ext_face[idim]->array(mfi);
        auto const &flag_info_face_x = flag_info_face[idim]->array(mfi);
        auto &borrowing_x = (*borrowing[idim])[mfi];
        auto const &borrowing_x_inds_pointer = borrowing_x.inds_pointer.array();
        auto const &borrowing_x_size = borrowing_x.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_x_inds = borrowing_x.inds.data();
        uint8_t* borrowing_x_neigh_faces = borrowing_x.neigh_faces.data();
        double* borrowing_x_area = borrowing_x.area.data();

        auto const &Sx_mod = area_mod[idim]->array(mfi);
        const auto &ly = edge_lengths[1]->array(mfi);
        const auto &lz = edge_lengths[2]->array(mfi);
        amrex::Real dy = cell_size[1];
        amrex::Real dz = cell_size[2];

        nelems_x = amrex::Scan::PrefixSum<int>(ncells,
        [=] AMREX_GPU_DEVICE (int icell) {
            amrex::Dim3 cell = box.atOffset(icell).dim3();
            int i = cell.x;
            int j = cell.y;
            int k = cell.z;
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
            int i = cell.x;
            int j = cell.y;
            int k = cell.z;
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
                        if( !(j_n == k_n or j_n == -k_n)){
                            // Here a face is available if it doesn't need to be extended itself and if its
                            // area exceeds Sz_ext. Here we need to take into account if the intruded face
                            // has given away already some area, so we use Sz_red rather than Sz.
                            // If no face is available we don't do anything and we will need to use the
                            // multi-face extensions.
                            if (Sx_mod(i, j + j_n, k + k_n) > Sx_ext
                                and ( flag_info_face_x(i, j + j_n, k + k_n) == 1
                                      or flag_info_face_x(i, j + j_n, k + k_n) == 2)
                                and flag_ext_face_x(i, j, k)) {
                                Sx_mod(i, j + j_n, k + k_n) -= Sx_ext;
                                // Insert the index of the face info
                                *(borrowing_x_inds + ps) = ps;
                                // Store the information about the intruded face in the dataset of the
                                // faces which are borrowing area
                                FaceInfoBox::addConnectedNeighbor(j_n, k_n, ps,
                                                                 borrowing_x_neigh_faces);
                                //*(borrowing_x_i_face + ps) = i;
                                //*(borrowing_x_j_face + ps) = j + j_n;
                                //*(borrowing_x_k_face + ps) = k + k_n;
                                *(borrowing_x_area + ps) = Sx_ext;

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

    for (amrex::MFIter mfi(*Bfield[idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sy = face_areas[idim]->array(mfi);
        auto const &flag_ext_face_y = flag_ext_face[idim]->array(mfi);
        auto const &flag_info_face_y = flag_info_face[idim]->array(mfi);
        auto &borrowing_y = (*borrowing[idim])[mfi];
        auto const &borrowing_y_inds_pointer = borrowing_y.inds_pointer.array();
        auto const &borrowing_y_size = borrowing_y.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_y_inds = borrowing_y.inds.data();
        uint8_t* borrowing_y_neigh_faces = borrowing_y.neigh_faces.data();
        double* borrowing_y_area = borrowing_y.area.data();

        auto const &Sy_mod = area_mod[idim]->array(mfi);
        const auto &lx = edge_lengths[0]->array(mfi);
        const auto &lz = edge_lengths[2]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dz = cell_size[2];

        nelems_y = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell) {
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
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
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                int nborrow = borrowing_y_size(i, j, k);

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
                            if( !(i_n == k_n or i_n == -k_n)){
                                // Here a face is available if it doesn't need to be extended itself and if its
                                // area exceeds Sz_ext. Here we need to take into account if the intruded face
                                // has given away already some area, so we use Sz_red rather than Sz.
                                // If no face is available we don't do anything and we will need to use the
                                // multi-face extensions.
                                if (Sy_mod(i + i_n, j, k + k_n) > Sy_ext
                                    and (flag_info_face_y(i + i_n, j, k + k_n) == 1
                                         or flag_info_face_y(i + i_n, j, k + k_n) == 2)
                                    and flag_ext_face_y(i, j, k)) {
                                    Sy_mod(i + i_n, j, k + k_n) -= Sy_ext;
                                    // Insert the index of the face info
                                    *(borrowing_y_inds + ps) = ps;
                                    // Store the information about the intruded face in the dataset of the
                                    // faces which are borrowing area
                                    FaceInfoBox::addConnectedNeighbor(i_n, k_n, ps,
                                                                      borrowing_y_neigh_faces);
                                    *(borrowing_y_area + ps) = Sy_ext;

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

        borrowing_y.inds.resize(nelems_y);

    }
    // Do the extensions in the z-plane
    idim = 2;

    for (amrex::MFIter mfi(*Bfield[idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sz = face_areas[idim]->array(mfi);
        auto const &flag_ext_face_z = flag_ext_face[idim]->array(mfi);
        auto const &flag_info_face_z = flag_info_face[idim]->array(mfi);
        auto &borrowing_z = (*borrowing[idim])[mfi];
        auto const &borrowing_z_inds_pointer = borrowing_z.inds_pointer.array();
        auto const &borrowing_z_size = borrowing_z.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_z_inds = borrowing_z.inds.data();
        uint8_t* borrowing_z_neigh_faces = borrowing_z.neigh_faces.data();
        double* borrowing_z_area = borrowing_z.area.data();

        auto const &Sz_mod = area_mod[idim]->array(mfi);
        const auto &lx = edge_lengths[0]->array(mfi);
        const auto &ly = edge_lengths[1]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dy = cell_size[1];


        nelems_z = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell) {
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_z(i, j, k)) {
                    //borrowing_z_size(i, j, k) = 0;
                    return 0;
                }

                amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                      ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
                amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);
                int n_borrow =
                    ComputeNBorrowOneFaceExtension(cell, Sz_ext, Sz_mod, flag_info_face_z,
                                                   flag_ext_face_z, idim);


              borrowing_z_size(i, j, k) = n_borrow;
                return n_borrow;
            },
            [=] AMREX_GPU_DEVICE (int icell, int ps){
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                int nborrow = borrowing_z_size(i, j, k);
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
                            if( !(i_n == j_n or i_n == -j_n)){
                                // Here a face is available if it doesn't need to be extended itself and if its
                                // area exceeds Sz_ext. Here we need to take into account if the intruded face
                                // has given away already some area, so we use Sz_red rather than Sz.
                                // If no face is available we don't do anything and we will need to use the
                                // multi-face extensions.
                                if (Sz_mod(i + i_n, j + j_n, k) > Sz_ext
                                    and (flag_info_face_z(i + i_n, j + j_n, k) == 1
                                         or flag_info_face_z(i + i_n, j + j_n, k) == 2)
                                    and flag_ext_face_z(i, j, k)) {
                                    Sz_mod(i + i_n, j + j_n, k) -= Sz_ext;
                                    // Insert the index of the face info
                                    *(borrowing_z_inds + ps) = ps;
                                    // Store the information about the intruded face in the dataset of the
                                    // faces which are borrowing area
                                    FaceInfoBox::addConnectedNeighbor(i_n, j_n, ps,
                                                                      borrowing_z_neigh_faces);
                                    *(borrowing_z_area+ps) = Sz_ext;

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

    return {nelems_x, nelems_y, nelems_z};
#else
    amrex::ignore_unused(Bfield, area_mod, edge_lengths, face_areas, flag_ext_face, flag_info_face,
                         borrowing, lev, flag_cp);
    return {0, 0, 0};
#endif
}

/**
* \brief Do the 8-ways extension
*/

void
WarpX::ComputeEightWaysExtensions(std::array<std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
                                  std::array<std::unique_ptr<amrex::MultiFab>, 3>& area_mod,
                                  std::array<std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
                                  std::array<std::unique_ptr<amrex::MultiFab>, 3> const& face_areas,
                                  std::array<std::unique_ptr<amrex::iMultiFab>, 3>& flag_ext_face,
                                  std::array<std::unique_ptr<amrex::iMultiFab>, 3>& flag_info_face,
                                  std::array<std::unique_ptr<amrex::LayoutData<FaceInfoBox>>, 3>& borrowing,
                                  amrex::Array1D<int, 0, 2> temp_inds, int lev, bool flag_cp) {
#ifdef AMREX_USE_EB

    //This variable is equal to lev if this is a fine patch or
    // equal to lev -1 if this is a coarse patch
    int lev_loc = lev;
    if(flag_cp and lev > 0){
        lev_loc = lev -1;
    }

    auto const &cell_size = CellSize(lev_loc);

    int n_borrow_offset_x = temp_inds(0);
    int n_borrow_offset_y = temp_inds(1);
    int n_borrow_offset_z = temp_inds(2);

    int idim = 0;
    for (amrex::MFIter mfi(*Bfield[idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sx = face_areas[idim]->array(mfi);
        auto const &flag_ext_face_x = flag_ext_face[idim]->array(mfi);
        auto const &flag_info_face_x = flag_info_face[idim]->array(mfi);
        auto &borrowing_x = (*borrowing[idim])[mfi];
        auto const &borrowing_x_inds_pointer = borrowing_x.inds_pointer.array();
        auto const &borrowing_x_size = borrowing_x.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_x_inds = borrowing_x.inds.data();
        uint8_t* borrowing_x_neigh_faces = borrowing_x.neigh_faces.data();
        double* borrowing_x_area = borrowing_x.area.data();

        auto const &Sx_mod = area_mod[idim]->array(mfi);
        const auto &ly = edge_lengths[1]->array(mfi);
        const auto &lz = edge_lengths[2]->array(mfi);
        amrex::Real dy = cell_size[1];
        amrex::Real dz = cell_size[2];
        int nelems_x;

        nelems_x = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell){
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_x(i, j, k)) {
                    return 0;
                }
                amrex::Real Sx_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                      lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
                amrex::Real Sx_ext = Sx_stab - Sx(i, j, k);
                int n_borrow = ComputeNBorrowEightFacesExtension(cell, Sx_ext, Sx_mod, Sx,
                                                                 flag_info_face_x, idim);

                borrowing_x_size(i, j, k) = n_borrow;
                return n_borrow;
                },
            [=] AMREX_GPU_DEVICE (int icell, int ps){
                ps += n_borrow_offset_x;

                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;

                if (!flag_ext_face_x(i, j, k)) {
                    return;
                }

                int nborrow = borrowing_x_size(i, j, k);
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
                                                         or flag_info_face_x(i, j  + j_loc - 1, k + k_loc - 1) == 2);
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

                    while(denom >= Sx_ext and neg_face and denom > 0){
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
                                    *(borrowing_x_inds + ps + count) = ps + count;
                                    // Store the information about the intruded face in the dataset of the
                                    // faces which are borrowing area
                                    FaceInfoBox::addConnectedNeighbor(j_n, k_n, ps + count,
                                                                      borrowing_x_neigh_faces);
                                    *(borrowing_x_area + ps + count) = patch;

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
    for (amrex::MFIter mfi(*Bfield[idim]); mfi.isValid(); ++mfi) {

        // Do the extensions in the y-plane
        amrex::Box const &box = mfi.validbox();

        auto const &Sy = face_areas[idim]->array(mfi);
        auto const &flag_ext_face_y = flag_ext_face[idim]->array(mfi);
        auto const &flag_info_face_y = flag_info_face[idim]->array(mfi);
        auto &borrowing_y = (*borrowing[idim])[mfi];
        auto const &borrowing_y_inds_pointer = borrowing_y.inds_pointer.array();
        auto const &borrowing_y_size = borrowing_y.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_y_inds = borrowing_y.inds.data();
        uint8_t* borrowing_y_neigh_faces = borrowing_y.neigh_faces.data();
        double* borrowing_y_area = borrowing_y.area.data();

        auto const &Sy_mod = area_mod[idim]->array(mfi);
        const auto &lx = edge_lengths[0]->array(mfi);
        const auto &lz = edge_lengths[2]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dz = cell_size[2];
        int nelems_y;

        nelems_y = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell){
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_y(i, j, k)){
                    return 0;
                }
                amrex::Real Sy_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k + 1) * dz,
                                                    lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
                amrex::Real Sy_ext = Sy_stab - Sy(i, j, k);
                int n_borrow = ComputeNBorrowEightFacesExtension(cell, Sy_ext, Sy_mod, Sy,
                                                                 flag_info_face_y, idim);

                borrowing_y_size(i, j, k) = n_borrow;
                return n_borrow;
            },
            [=] AMREX_GPU_DEVICE (int icell, int ps) {

                ps += n_borrow_offset_y;

                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;

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
                                                         or flag_info_face_y(i + i_loc - 1, j, k + k_loc - 1) == 2)   ;
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

                    while (denom >= Sy_ext and neg_face and denom > 0) {
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
                                    *(borrowing_y_inds + ps + count) = ps + count;
                                    FaceInfoBox::addConnectedNeighbor(i_n, k_n, ps + count,
                                                                      borrowing_y_neigh_faces);
                                    *(borrowing_y_area + ps + count) = patch;

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
            amrex::Scan::Type::exclusive
        );
    }

    // Do the extensions in the z-plane
    idim = 2;
    for (amrex::MFIter mfi(*Bfield[idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sz = face_areas[idim]->array(mfi);
        auto const &flag_ext_face_z = flag_ext_face[idim]->array(mfi);
        auto const &flag_info_face_z = flag_info_face[idim]->array(mfi);
        auto &borrowing_z = (*borrowing[idim])[mfi];
        auto const &borrowing_z_inds_pointer = borrowing_z.inds_pointer.array();
        auto const &borrowing_z_size = borrowing_z.size.array();
        amrex::Long ncells = box.numPts();
        int* borrowing_z_inds = borrowing_z.inds.data();
        uint8_t* borrowing_z_neigh_faces = borrowing_z.neigh_faces.data();
        double* borrowing_z_area = borrowing_z.area.data();

        auto const &Sz_mod = area_mod[idim]->array(mfi);
        const auto &lx = edge_lengths[0]->array(mfi);
        const auto &ly = edge_lengths[1]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dy = cell_size[1];
        int nelems_z;

        nelems_z = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell){
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_z(i, j, k)) {
                    return 0;
                }
                amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                    ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
                amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);
                int n_borrow = ComputeNBorrowEightFacesExtension(cell, Sz_ext, Sz_mod, Sz,
                                                                 flag_info_face_z, idim);

                borrowing_z_size(i, j, k) = n_borrow;
                return n_borrow;
            },
            [=] AMREX_GPU_DEVICE (int icell, int ps) {

                ps += n_borrow_offset_z;

                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;

                if (!flag_ext_face_z(i, j, k)) {
                    return;
                }

                int nborrow = borrowing_z_size(i, j, k);
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
                                                         or flag_info_face_z(i + i_loc - 1, j + j_loc - 1, k) == 2);
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

                    while(denom >= Sz_ext and neg_face and denom > 0){
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
                                    *(borrowing_z_inds + ps + count) = ps + count;
                                    FaceInfoBox::addConnectedNeighbor(i_n, j_n, ps + count,
                                                                      borrowing_z_neigh_faces);
                                    *(borrowing_z_area + ps + count) = patch;

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
#else
    amrex::ignore_unused(Bfield, area_mod, edge_lengths, face_areas, flag_ext_face, flag_info_face,
                         borrowing, temp_inds, lev, flag_cp);
#endif
}

