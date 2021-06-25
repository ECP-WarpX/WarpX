#include "WarpX.H"
#include <AMReX_Scan.H>

/**
* \brief auxiliary function to count the amount of faces which still need to be extended
*/
int
WarpX::CountExtFaces() {
    int sum = 0;
#ifdef AMREX_USE_EB
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        for (amrex::MFIter mfi(*m_flag_ext_face[maxLevel()][idim]); mfi.isValid(); ++mfi) {
            amrex::Box const &box = mfi.validbox();
            auto const &flag_ext_face = m_flag_ext_face[maxLevel()][idim]->array(mfi);
            amrex::ParallelFor(box, [=, &sum](int i, int j, int k) {
                // Do I need to flag the reduction somehow?
                sum = sum + flag_ext_face(i, j, k);
            });
        }
    }
#endif
    return sum;
}


/**
* \brief Compute the face extensions. For the moment it only uses the one-way extension but at a
* later point we will also have the 4-ways and 8-ways extensions.
*/
void
WarpX::ComputeFaceExtensions(){
#ifdef AMREX_USE_EB
    unsigned int N_ext_faces = CountExtFaces();
    std::cout<< "Cells to be extended before:\t" << N_ext_faces <<std::endl;
    amrex::Array1D<int, 0, 2> temp_inds = ComputeOneWayExtensions();
    unsigned int N_ext_faces_after_one_way = CountExtFaces();
    std::cout<< "Cells to be extended after one way extension:\t" <<
                N_ext_faces_after_one_way <<std::endl;
    //ComputeEightWaysExtensions(temp_inds);
    //unsigned int N_ext_faces_after_eight_ways = CountExtFaces();
    //std::cout<< "Cells to be extended after eight ways extension:\t" <<
    //            N_ext_faces_after_eight_ways <<std::endl;
    if (N_ext_faces_after_one_way > 0) {
        amrex::Abort("Some faces could not be extended");
    }
#endif
}

/**
* \brief For the face of cell pointing in direction idim, compute the number of faces
* we need to intrude. For the one-way extension this function returns only one or zero: one if the
* face can be extended withe the one-way extension, zeros if it can't.
*/
int WarpX::ComputeNBorrowOneCellExtension(amrex::Dim3 cell, amrex::Real S_ext,
                                          const amrex::Array4<amrex::Real>& S_red,
                                          const amrex::Array4<int>& flag_avail_face,
                                          const amrex::Array4<int>& flag_ext_face, int idim) {
    int i = cell.x;
    int j = cell.y;
    int k = cell.z;
    int nelems = 0;
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
                        and flag_avail_face(i, j + j_n, k + k_n)
                        and flag_ext_face(i, j, k) and not stop) {
                        nelems += 1;
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
                        and flag_avail_face(i + i_n, j, k + k_n)
                        and flag_ext_face(i, j, k) and not stop) {
                        nelems += 1;
                        stop = true;
                    }
                }
            }
        }
    }else if(idim == 2){
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
                        and flag_avail_face(i + i_n, j + j_n, k)
                        and flag_ext_face(i, j, k) and not stop) {
                        nelems += 1;
                        stop = true;
                    }
                }
            }
        }
    }else{
        amrex::Abort("Template argument should be 0, 1 or 2");
    }
    return nelems;
}

/**
* \brief Do the one-way extension
*/
amrex::Array1D<int, 0, 2>
WarpX::ComputeOneWayExtensions() {
#ifdef AMREX_USE_EB
    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &cell_size = CellSize(maxLevel());

    int nelems_x;
    int nelems_y;
    int nelems_z;

    // Do the extensions in the x-plane
    int idim = 0;

    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][0]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sx = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_intr_face_x = m_flag_intr_face[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_x = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_x = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_x = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_x.inds_pointer.resize(box);
        borrowing_x.size.resize(box);
        auto const &borrowing_x_inds_pointer = borrowing_x.inds_pointer.array();
        auto const &borrowing_x_size = borrowing_x.size.array();
        amrex::Long ncells = box.numPts();
        borrowing_x.inds.resize(8*ncells);
        borrowing_x.i_face.resize(8*ncells);
        borrowing_x.j_face.resize(8*ncells);
        borrowing_x.k_face.resize(8*ncells);
        borrowing_x.area.resize(8*ncells);
        borrowing_x.rho_face.resize(8*ncells);
        int* borrowing_x_inds = borrowing_x.inds.data();
        int* borrowing_x_i_face = borrowing_x.i_face.data();
        int* borrowing_x_j_face = borrowing_x.j_face.data();
        int* borrowing_x_k_face = borrowing_x.k_face.data();
        double* borrowing_x_rho_area = borrowing_x.rho_face.data();
        double* borrowing_x_area = borrowing_x.area.data();

        auto const &Sx_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sx_enl = m_area_enl[maxLevel()][idim]->array(mfi);
        const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
        const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
        amrex::Real dy = cell_size[1];
        amrex::Real dz = cell_size[2];

        nelems_x = amrex::Scan::PrefixSum<int>(ncells,
        [=] (int icell) {
            amrex::Dim3 cell = box.atOffset(icell).dim3();
            int i = cell.x;
            int j = cell.y;
            int k = cell.z;
            // If the face doesn't need to be extended break the loop
            if (!flag_ext_face_x(i, j, k))
                return 0;
            //one cell extension, therefore the_size_for_this_cell(cell) = 1 or 0
            amrex::Real Sx_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                  lz(i, j, k) * dy, lz(i, j + 1, k) * dy});
            amrex::Real Sx_ext = Sx_stab - Sx(i, j, k);
            int n_borrow = the_size_for_this_cell(cell, Sx_ext, Sx_red, flag_avail_face_x,
                                                  flag_ext_face_x, idim);

            //int n_borrow = 1;
            //int n_lend = 1;
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
                            if (Sx_red(i, j + j_n, k + k_n) > Sx_ext
                                and flag_avail_face_x(i, j + j_n, k + k_n)
                                and flag_ext_face_x(i, j, k)) {
                                Sx_red(i, j + j_n, k + k_n) -= Sx_ext;
                                // Insert the index of the face info
                                *(borrowing_x_inds + ps) = ps;
                                // advance the borrowing index
                                // Store the information about the intruded face in the dataset of the
                                // faces which are borrowing area
                                *(borrowing_x_i_face + ps) = i;
                                *(borrowing_x_j_face + ps) = j + j_n;
                                *(borrowing_x_k_face + ps) = k + k_n;
                                *(borrowing_x_area + ps) = Sx_ext;

                                // rho_face is computed at every time step. Here we just insert a NaN to
                                // allocate the space in the rho_face vector and we will overwrite it later.
                                *(borrowing_x_rho_area + ps) = NAN;
                                flag_intr_face_x(i, j + j_n, k + k_n) = true;
                                // Add the area to the intruding face.
                                // TODO: this is not needed here, but it's needed in the multi-face
                                //  extensions. Should I remove it?
                                Sx_enl(i, j, k) = Sx(i, j, k) + Sx_ext;
                                flag_ext_face_x(i, j, k) = false;
                            }
                        }
                    }
                }
            }
        }, amrex::Scan::Type::exclusive);

        borrowing_x.inds.resize(nelems_x);

    }
    // Do the extensions in the y-plane
    idim = 1;

    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sy = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_intr_face_y = m_flag_intr_face[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_y = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_y = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_y = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_y.inds_pointer.resize(box);
        borrowing_y.size.resize(box);
        auto const &borrowing_inds_pointer_y = borrowing_y.inds_pointer.array();
        auto const &borrowing_size_y = borrowing_y.size.array();
        amrex::Long ncells = box.numPts();
        borrowing_y.inds.resize(8*ncells);
        borrowing_y.i_face.resize(8*ncells);
        borrowing_y.j_face.resize(8*ncells);
        borrowing_y.k_face.resize(8*ncells);
        borrowing_y.area.resize(8*ncells);
        borrowing_y.rho_face.resize(8*ncells);
        int* inds_ptr_borrowing_y = borrowing_y.inds.data();
        int* i_face_ptr_borrowing_y = borrowing_y.i_face.data();
        int* j_face_ptr_borrowing_y = borrowing_y.j_face.data();
        int* k_face_ptr_borrowing_y = borrowing_y.k_face.data();
        double* rho_face_ptr_borrowing_y = borrowing_y.rho_face.data();
        double* area_ptr_borrowing_y = borrowing_y.area.data();
        auto const &Sy_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sy_enl = m_area_enl[maxLevel()][idim]->array(mfi);
        const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
        const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dz = cell_size[2];

        nelems_y = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell) {
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_y(i, j, k))
                    return 0;
                //one cell extension, therefore the_size_for_this_cell(cell) = 1 or 0
                amrex::Real Sy_stab = 0.5 * std::max({lx(i, j, k) * dz, lx(i, j, k + 1) * dz,
                                                      lz(i, j, k) * dx, lz(i + 1, j, k) * dx});
                amrex::Real Sy_ext = Sy_stab - Sy(i, j, k);
                int n_borrow = the_size_for_this_cell(cell, Sy_ext, Sy_red, flag_avail_face_y,
                                                      flag_ext_face_y, idim);

                borrowing_size_y(i, j, k) = n_borrow;
                return n_borrow;
                },
            [=] AMREX_GPU_DEVICE (int icell, int ps){
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                int nborrow = borrowing_size_y(i, j, k);
                if (nborrow == 0) {
                    borrowing_inds_pointer_y(i, j, k) = nullptr;
                } else{
                    borrowing_inds_pointer_y(i, j, k) = inds_ptr_borrowing_y + ps;

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
                                if (Sy_red(i + i_n, j, k + k_n) > Sy_ext
                                    and flag_avail_face_y(i + i_n, j, k + k_n)
                                    and flag_ext_face_y(i, j, k)) {
                                    Sy_red(i + i_n, j, k + k_n) -= Sy_ext;
                                    // Insert the index of the face info
                                    *(inds_ptr_borrowing_y + ps) = ps;
                                    // Store the information about the intruded face in the dataset of the
                                    // faces which are borrowing area
                                    *(i_face_ptr_borrowing_y + ps) = i + i_n;
                                    *(j_face_ptr_borrowing_y + ps) = j;
                                    *(k_face_ptr_borrowing_y + ps) = k + k_n;
                                    *(area_ptr_borrowing_y + ps) = Sy_ext;

                                    // rho_face is computed at every time step. Here we just insert a NaN to
                                    // allocate the space in the rho_face vector and we will overwrite it later.
                                    *(rho_face_ptr_borrowing_y + ps) = NAN;

                                    flag_intr_face_y(i + i_n, j, k + k_n) = true;
                                    // Add the area to the intruding face.
                                    // TODO: this is not needed here, but it's needed in the multi-face
                                    //  extensions. Should I remove it?
                                    Sy_enl(i, j, k) = Sy(i, j, k) + Sy_ext;
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

    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][idim]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        auto const &Sz = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_intr_face_z = m_flag_intr_face[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_z = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_z = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &borrowing_z = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_z.inds_pointer.resize(box);
        borrowing_z.size.resize(box);
        auto const &borrowing_inds_pointer_z = borrowing_z.inds_pointer.array();
        auto const &borrowing_size_z = borrowing_z.size.array();
        amrex::Long ncells = box.numPts();
        borrowing_z.inds.resize(8*ncells);
        borrowing_z.i_face.resize(8*ncells);
        borrowing_z.j_face.resize(8*ncells);
        borrowing_z.k_face.resize(8*ncells);
        borrowing_z.area.resize(8*ncells);
        borrowing_z.rho_face.resize(8*ncells);
        int* inds_ptr_borrowing_z = borrowing_z.inds.data();
        int* i_face_ptr_borrowing_z = borrowing_z.i_face.data();
        int* j_face_ptr_borrowing_z = borrowing_z.j_face.data();
        int* k_face_ptr_borrowing_z = borrowing_z.k_face.data();
        double* rho_face_ptr_borrowing_z = borrowing_z.rho_face.data();
        double* area_ptr_borrowing_z = borrowing_z.area.data();

        auto const &Sz_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sz_enl = m_area_enl[maxLevel()][idim]->array(mfi);
        const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
        const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dy = cell_size[1];


        nelems_z = amrex::Scan::PrefixSum<int>(ncells,
            [=] AMREX_GPU_DEVICE (int icell) {
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                // If the face doesn't need to be extended break the loop
                if (!flag_ext_face_z(i, j, k))
                    return 0;
                //one cell extension, therefore the_size_for_this_cell(cell) = 1 or 0
                amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                      ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
                amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);
                int n_borrow = the_size_for_this_cell(cell, Sz_ext, Sz_red, flag_avail_face_z,
                                                      flag_ext_face_z, idim);


                borrowing_size_z(i, j, k) = n_borrow;
                return n_borrow;
            },
            [=] AMREX_GPU_DEVICE (int icell, int ps){
                amrex::Dim3 cell = box.atOffset(icell).dim3();
                int i = cell.x;
                int j = cell.y;
                int k = cell.z;
                int nborrow = borrowing_size_z(i, j, k);
                if (nborrow == 0) {
                    borrowing_inds_pointer_z(i, j, k) = nullptr;
                } else{
                    borrowing_inds_pointer_z(i, j, k) = inds_ptr_borrowing_z + ps;

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
                                if (Sz_red(i + i_n, j + j_n, k) > Sz_ext
                                    and flag_avail_face_z(i + i_n, j + j_n, k)
                                    and flag_ext_face_z(i, j, k)) {
                                    Sz_red(i + i_n, j + j_n, k) -= Sz_ext;
                                    // Insert the index of the face info
                                    *(inds_ptr_borrowing_z + ps) = ps;
                                    // advance the borrowing index
                                    // Store the information about the intruded face in the dataset of the
                                    // faces which are borrowing area
                                    *(i_face_ptr_borrowing_z+ps) = i + i_n;
                                    *(j_face_ptr_borrowing_z+ps) = j + j_n;
                                    *(k_face_ptr_borrowing_z+ps) = k;

                                    *(area_ptr_borrowing_z+ps) = Sz_ext;
                                    // rho_face is computed at every time step. Here we just insert a NaN to
                                    // allocate the space in the rho_face vector and we will overwrite it later.
                                    *(rho_face_ptr_borrowing_z + ps) = NAN;

                                    flag_intr_face_z(i + i_n, j + j_n, k) = true;
                                    // Add the area to the intruding face.
                                    // TODO: this is not needed here, but it's needed in the multi-face
                                    //  extensions. Should I remove it?
                                    Sz_enl(i, j, k) = Sz(i, j, k) + Sz_ext;
                                    flag_ext_face_z(i, j, k) = false;
                                }
                            }
                        }
                    }
                }
            },
        amrex::Scan::Type::exclusive);

        borrowing_z.inds.resize(nelems_z);
    }

    return {nelems_x, nelems_y, nelems_z};
#else
    return {0, 0, 0};
#endif
}

/**
* \brief Do the 8-ways extension
*/
/*
void
WarpX::ComputeEightWaysExtensions(amrex::Array1D<int, 0, 5> temp_inds) {
#ifdef AMREX_USE_EB
    int temp_ind_borrowing_x = temp_inds(0);
    int temp_ind_lending_x = temp_inds(1);
    int temp_ind_borrowing_y = temp_inds(2);
    int temp_ind_lending_y = temp_inds(3);
    int temp_ind_borrowing_z = temp_inds(4);
    int temp_ind_lending_z = temp_inds(5);
    auto const &cell_size = CellSize(maxLevel());

    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][2]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
        const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
        const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dy = cell_size[1];
        amrex::Real dz = cell_size[2];

        // Do the extensions in the x-plane
        int idim = 0;
        auto const &Sx = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_x = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_x = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &lending_x = (*m_lending[maxLevel()][idim])[mfi];
        auto &borrowing_x = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_x.inds_pointer.resize(box);
        lending_x.inds_pointer.resize(box);
        borrowing_x.size.resize(box);
        lending_x.size.resize(box);
        auto const &lending_inds_pointer_x = lending_x.inds_pointer.array();
        auto const &borrowing_inds_pointer_x = borrowing_x.inds_pointer.array();
        auto const &lending_size_x = lending_x.size.array();
        auto const &borrowing_size_x = borrowing_x.size.array();
        auto const &Sx_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sx_enl = m_area_enl[maxLevel()][idim]->array(mfi);

        amrex::ParallelFor(box,
                         [=, &borrowing_x, &lending_x,
                             &temp_ind_borrowing_x, &temp_ind_lending_x](int i, int j, int k) {

            if( !flag_ext_face_x(i, j, k) ) return;

            Sx_enl(i, j, k) = Sx(i, j, k);
            amrex::Real Sx_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                  lz(i, j, k) * dy, lz(i, j, k + 1) * dy});
            amrex::Real Sx_ext = Sx_stab - Sx(i, j, k);
            amrex::Array2D<amrex::Real, 0, 2, 0, 2> local_avail{};
            for(int j_loc = 0; j_loc <= 2; j_loc++){
                for(int k_loc = 0; k_loc <= 2; k_loc++){
                    local_avail(j_loc, k_loc) = flag_avail_face_x(i, j  + j_loc - 1, k + k_loc - 1);
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

            bool neg_cell = true;

            while(denom >= Sx_ext and neg_cell and denom > 0){
                neg_cell = false;
                for (int j_n = -1; j_n < 2; j_n++) {
                    for (int k_n = -1; k_n < 2; k_n++) {
                        if(local_avail(j_n + 1, k_n + 1)){
                            amrex::Real patch = Sx_ext * Sx(i, j + j_n, k+ k_n) / denom;
                            if(Sx_red(i, j + j_n, k + k_n) - patch <= 0) {
                                neg_cell = true;
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
                Sx_enl(i, j, k) = Sx(i, j, k);
                for (int j_n = -1; j_n < 2; j_n++) {
                    for (int k_n = -1; k_n < 2; k_n++) {
                        if(local_avail(j_n + 1, k_n + 1)){
                            amrex::Real patch = Sx_ext * Sx(i, j + j_n, k + k_n) / denom;
                            borrowing_x.inds.push_back(temp_ind_borrowing_x);
                            // If it has not been done store the pointer
                            if(borrowing_size_x(i, j, k) == 0){
                                borrowing_inds_pointer_x(i, j, k) = &borrowing_x.inds.back();
                            }
                            // increase the size
                            borrowing_size_x(i, j, k) += 1;
                            // advance the borrowing index
                            temp_ind_borrowing_x = temp_ind_borrowing_x + 1;
                            borrowing_x.i_face.push_back(i);
                            borrowing_x.j_face.push_back(j + j_n);
                            borrowing_x.k_face.push_back(k + k_n);
                            borrowing_x.area.push_back(Sx_ext);
                            borrowing_x.rho_face.push_back(NAN);
                            lending_x.inds.push_back(temp_ind_lending_x);
                            if(lending_size_x(i, j+j_n, k+k_n) == 0){
                                //push back and initialize the pointer
                                lending_x.inds.push_back(temp_ind_lending_x);
                                lending_inds_pointer_x(i, j+j_n, k+k_n) = &lending_x.inds.back();
                            }else{
                                //insert it where it has to be.
                                // To do so, first we get the pointer to the first index relative to
                                // this face
                                int* p = lending_inds_pointer_x(i, j+j_n, k+k_n);
                                // Then we get an iterator starting at the position of that index
                                amrex::Gpu::DeviceVector<int>::iterator it = lending_x.inds.begin()
                                    + std::distance(lending_x.inds.data(), p);
                                // And now we insert the index as first index relative to this face
                                //lending_x.inds.
                                lending_x.inds.insert(it, temp_ind_lending_x);

                            }
                            // increase the size
                            lending_size_x(i, j + j_n, k + k_n) += 1;
                            // advance the lending index
                            temp_ind_lending_x = temp_ind_lending_x + 1;
                            lending_x.i_face.push_back(i);
                            lending_x.j_face.push_back(j);
                            lending_x.k_face.push_back(k);
                            lending_x.area.push_back(Sx_ext);
                            lending_x.rho_face.push_back(NAN);
                            Sx_enl(i, j, k) += patch;
                            Sx_red(i, j + j_n, k + k_n) -= patch;
                        }
                    }
                }
                flag_ext_face_x(i, j, k) = false;
            }
        });

        // Do the extensions in the y-plane
        idim = 1;
        auto const &Sy = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_y = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_y = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &lending_y = (*m_lending[maxLevel()][idim])[mfi];
        auto &borrowing_y = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_y.inds_pointer.resize(box);
        lending_y.inds_pointer.resize(box);
        borrowing_y.size.resize(box);
        lending_y.size.resize(box);
        auto const &lending_inds_pointer_y = lending_y.inds_pointer.array();
        auto const &borrowing_inds_pointer_y = borrowing_y.inds_pointer.array();
        auto const &lending_size_y = lending_y.size.array();
        auto const &borrowing_size_y = borrowing_y.size.array();
        auto const &Sy_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sy_enl = m_area_enl[maxLevel()][idim]->array(mfi);

        amrex::ParallelFor(box,
                         [=, &borrowing_y, &lending_y,
                             &temp_ind_borrowing_y, &temp_ind_lending_y](int i, int j, int k) {

            if( !flag_ext_face_y(i, j, k) ) return;

            Sy_enl(i, j, k) = Sy(i, j, k);
            amrex::Real Sy_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                  lz(i, j, k) * dy, lz(i, j, k + 1) * dy});
            amrex::Real Sy_ext = Sy_stab - Sy(i, j, k);
            amrex::Array2D<amrex::Real, 0, 2, 0, 2> local_avail{};
            for(int i_loc = 0; i_loc <= 2; i_loc++){
                for(int k_loc = 0; k_loc <= 2; k_loc++){
                    local_avail(i_loc, k_loc) = flag_avail_face_y(i + i_loc - 1, j, k + k_loc - 1);
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

            bool neg_cell = true;

            while(denom >= Sy_ext and neg_cell and denom > 0){
                neg_cell = false;
                for (int i_n = -1; i_n < 2; i_n++) {
                    for (int k_n = -1; k_n < 2; k_n++) {
                        if(local_avail(i_n + 1, k_n + 1)){
                            amrex::Real patch = Sy_ext * Sy(i + i_n, j, k + k_n) / denom;
                            if(Sy_red(i + i_n, j, k + k_n) - patch <= 0) {
                                neg_cell = true;
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

            if(denom >= Sy_ext){
                Sy_enl(i, j, k) = Sy(i, j, k);
                for (int i_n = -1; i_n < 2; i_n++) {
                    for (int k_n = -1; k_n < 2; k_n++) {
                        if(local_avail(i_n + 1, k_n + 1)){
                            amrex::Real patch = Sy_ext * Sy(i + i_n, j, k + k_n) / denom;
                            borrowing_y.inds.push_back(temp_ind_borrowing_y);
                            // If it has not been done store the pointer
                            if(borrowing_size_y(i, j, k) == 0){
                                borrowing_inds_pointer_y(i, j, k) = &borrowing_y.inds.back();
                            }
                            // increase the size
                            borrowing_size_y(i, j, k) += 1;
                            // advance the lending index
                            temp_ind_borrowing_y = temp_ind_borrowing_y + 1;
                            borrowing_y.i_face.push_back(i + i_n);
                            borrowing_y.j_face.push_back(j);
                            borrowing_y.k_face.push_back(k + k_n);
                            borrowing_y.area.push_back(Sy_ext);
                            borrowing_y.rho_face.push_back(NAN);
                            lending_y.inds.push_back(temp_ind_lending_y);
                            if(lending_size_y(i + i_n, j, k+k_n) == 0){
                                //push back and initialize the pointer
                                lending_y.inds.push_back(temp_ind_lending_x);
                                lending_inds_pointer_y(i+i_n, j, k+k_n) = &lending_y.inds.back();
                            }else{
                                //insert it where it has to be.
                                // To do so, first we get the pointer to the first index relative to
                                // this face
                                int* p = lending_inds_pointer_y(i + i_n, j, k+k_n);
                                // Then we get an iterator starting at the position of that index
                                amrex::Gpu::DeviceVector<int>::iterator it = lending_y.inds.begin()
                                    + std::distance(lending_y.inds.data(), p);
                                // And now we insert the index as first index relative to this face
                                //lending_x.inds.
                                lending_y.inds.insert(it, temp_ind_lending_y);

                            }
                            // increase the size
                            lending_size_y(i + i_n, j, k + k_n) += 1;
                            // advance the lending index
                            temp_ind_lending_y = temp_ind_lending_y + 1;
                            lending_y.i_face.push_back(i);
                            lending_y.j_face.push_back(j);
                            lending_y.k_face.push_back(k);
                            lending_y.area.push_back(Sy_ext);
                            lending_y.rho_face.push_back(NAN);
                            Sy_enl(i, j, k) += patch;
                            Sy_red(i + i_n, j, k + k_n) -= patch;
                        }
                    }
                }
                flag_ext_face_y(i, j, k) = false;
            }
        });

        // Do the extensions in the z-plane
        idim = 2;
        auto const &Sz = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_z = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_z = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &lending_z = (*m_lending[maxLevel()][idim])[mfi];
        auto &borrowing_z = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_z.inds_pointer.resize(box);
        lending_z.inds_pointer.resize(box);
        borrowing_z.size.resize(box);
        lending_z.size.resize(box);
        auto const &borrowing_inds_pointer_z = borrowing_z.inds_pointer.array();
        auto const &lending_inds_pointer_z = lending_z.inds_pointer.array();
        auto const &lending_size_z = lending_z.size.array();
        auto const &borrowing_size_z = borrowing_z.size.array();
        auto const &Sz_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sz_enl = m_area_enl[maxLevel()][idim]->array(mfi);

        amrex::ParallelFor(box,
                         [=, &borrowing_z, &lending_z,
                             &temp_ind_borrowing_z, &temp_ind_lending_z](int i, int j, int k) {

            if( !flag_ext_face_z(i, j, k) ) return;

            Sz_enl(i, j, k) = Sz(i, j, k);
            amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                  ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
            amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);
            amrex::Array2D<amrex::Real, 0, 2, 0, 2> local_avail{};
            for(int i_loc = 0; i_loc <= 2; i_loc++){
                for(int j_loc = 0; j_loc <= 2; j_loc++){
                    local_avail(i_loc, j_loc) = flag_avail_face_z(i + i_loc - 1, j + j_loc - 1, k);
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

            bool neg_cell = true;

            while(denom >= Sz_ext and neg_cell and denom > 0){
                neg_cell = false;
                for (int i_n = -1; i_n < 2; i_n++) {
                    for (int j_n = -1; j_n < 2; j_n++) {
                        if(local_avail(i_n + 1, j_n + 1)){
                            amrex::Real patch = Sz_ext * Sz(i + i_n, j + j_n, k) / denom;
                            if(Sz_red(i + i_n, j + j_n, k) - patch <= 0) {
                                neg_cell = true;
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
                Sz_enl(i, j, k) = Sz(i, j, k);
                for (int i_n = -1; i_n < 2; i_n++) {
                    for (int j_n = -1; j_n < 2; j_n++) {
                        if(local_avail(i_n + 1, j_n + 1)){
                            amrex::Real patch = Sz_ext * Sz(i + i_n, j + j_n, k) / denom;
                            borrowing_z.inds.push_back(temp_ind_borrowing_z);
                            // If it has not been done store the pointer
                            if(borrowing_size_z(i, j, k) == 0){
                                borrowing_inds_pointer_z(i, j, k) = &borrowing_z.inds.back();
                            }
                            // increase the size
                            borrowing_size_z(i, j, k) += 1;
                            // advance the lending index
                            temp_ind_borrowing_z = temp_ind_borrowing_z + 1;
                            borrowing_z.i_face.push_back(i + i_n);
                            borrowing_z.j_face.push_back(j + j_n);
                            borrowing_z.k_face.push_back(k);
                            borrowing_z.area.push_back(Sz_ext);
                            borrowing_z.rho_face.push_back(NAN);
                            lending_z.inds.push_back(temp_ind_lending_z);
                            if(lending_size_y(i + i_n, j + j_n, k) == 0){
                                //push back and initialize the pointer
                                lending_z.inds.push_back(temp_ind_lending_x);
                                lending_inds_pointer_z(i+i_n, j + j_n, k) = &lending_z.inds.back();
                            }else {
                                //insert it where it has to be.
                                // To do so, first we get the pointer to the first index relative to
                                // this face
                                int *p = lending_inds_pointer_y(i + i_n, j + j_n, k);
                                // Then we get an iterator starting at the position of that index
                                amrex::Gpu::DeviceVector<int>::iterator it = lending_z.inds.begin()
                                    + std::distance(lending_z.inds.data(), p);
                                // And now we insert the index as first index relative to this face
                                //lending_x.inds.
                                lending_z.inds.insert(it, temp_ind_lending_z);

                            }
                            // increase the size
                            lending_size_z(i, j, k) += 1;
                            // advance the lending index
                            temp_ind_lending_z = temp_ind_lending_z + 1;
                            lending_z.i_face.push_back(i);
                            lending_z.j_face.push_back(j);
                            lending_z.k_face.push_back(k);
                            lending_z.area.push_back(Sz_ext);
                            lending_z.rho_face.push_back(NAN);
                            Sz_enl(i, j, k) += patch;
                            Sz_red(i + i_n, j + j_n, k) -= patch;
                        }
                    }
                }
                flag_ext_face_z(i, j, k) = false;
            }
        });

    }
#else
    amrex::ignore_unused(temp_inds);
#endif
}
*/
