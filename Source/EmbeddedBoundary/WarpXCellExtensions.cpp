#include "WarpX.H"

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
    amrex::Array1D<int, 0, 5> temp_inds = ComputeOneWayExtensions();
    unsigned int N_ext_faces_after_one_way = CountExtFaces();
    std::cout<< "Cells to be extended after one way extension:\t" <<
                N_ext_faces_after_one_way <<std::endl;
    ComputeEightWaysExtensions(temp_inds);
    unsigned int N_ext_faces_after_eight_ways = CountExtFaces();
    std::cout<< "Cells to be extended after eight ways extension:\t" <<
                N_ext_faces_after_eight_ways <<std::endl;
    if (N_ext_faces_after_eight_ways > 0) {
        amrex::Abort("Some faces could not be extended");
    }
#endif
}

/**
* \brief Do the one-way extension
*/
amrex::Array1D<int, 0, 5>
WarpX::ComputeOneWayExtensions() {
#ifdef AMREX_USE_EB
    auto const eb_fact = fieldEBFactory(maxLevel());
    auto const &cell_size = CellSize(maxLevel());

    int temp_ind_borrowing_x = 0;
    int temp_ind_lending_x = 0;
    int temp_ind_borrowing_y = 0;
    int temp_ind_lending_y = 0;
    int temp_ind_borrowing_z = 0;
    int temp_ind_lending_z = 0;

    for (amrex::MFIter mfi(*Bfield_fp[maxLevel()][2]); mfi.isValid(); ++mfi) {

        amrex::Box const &box = mfi.validbox();

        // Do the extensions in the x-plane
        int idim = 0;
        auto const &Sx = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_x = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_x = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &lending_x = (*m_lending[maxLevel()][idim])[mfi];
        auto &borrowing_x = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_x.inds.resize(box);
        lending_x.inds.resize(box);
        auto const &lending_inds_x = lending_x.inds.array();
        auto const &borrowing_inds_x = borrowing_x.inds.array();
        auto const &Sx_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sx_enl = m_area_enl[maxLevel()][idim]->array(mfi);
        const auto &lx = m_edge_lengths[maxLevel()][0]->array(mfi);
        const auto &ly = m_edge_lengths[maxLevel()][1]->array(mfi);
        const auto &lz = m_edge_lengths[maxLevel()][2]->array(mfi);
        amrex::Real dx = cell_size[0];
        amrex::Real dy = cell_size[1];
        amrex::Real dz = cell_size[2];

        amrex::ParallelFor(box,
                         [=, &borrowing_x, &lending_x,
                             &temp_ind_borrowing_x, &temp_ind_lending_x] AMREX_GPU_HOST(int i, int j, int k) {

            // If the face doesn't need to be extended break the loop
            if( !flag_ext_face_x(i, j, k) ) return;

            amrex::Real Sx_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                  lz(i, j, k) * dy, lz(i, j, k + 1) * dy});
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
                            borrowing_inds_x(i, j, k).push_back(temp_ind_borrowing_x);
                            temp_ind_borrowing_x = temp_ind_borrowing_x + 1;
                            // Store the information about the intruded face in the dataset of the
                            // faces which are borrowing area
                            borrowing_x.i_face.push_back(i);
                            borrowing_x.j_face.push_back(j + j_n);
                            borrowing_x.k_face.push_back(k + k_n);
                            borrowing_x.area.push_back(Sx_ext);
                            // rho_face is computed at every time step. Here we just insert a NaN to
                            // allocate the space in the rho_face vector and we will overwrite it later.
                            borrowing_x.rho_face.push_back(NAN);
                            lending_inds_x(i, j + j_n, k + k_n).push_back(temp_ind_lending_x);
                            // Store the information about the intruding face in the dataset of the
                            // faces which are lending area
                            temp_ind_lending_x = temp_ind_lending_x + 1;
                            lending_x.i_face.push_back(i);
                            lending_x.j_face.push_back(j);
                            lending_x.k_face.push_back(k);
                            lending_x.area.push_back(Sx_ext);
                            lending_x.rho_face.push_back(NAN);
                            // Add the area to the intruding face.
                            // TODO: this is not needed here, but it's needed in the multi-face
                            //  extensions. Should I remove it?
                            Sx_enl(i, j, k) = Sx(i, j, k) + Sx_ext;
                            flag_ext_face_x(i, j, k) = false;
                        }
                    }
                }
            }
        });

        // Do the extensions in the y-plane
        idim = 1;
        auto const &Sy = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_y = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_y = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &lending_y = (*m_lending[maxLevel()][idim])[mfi];
        auto &borrowing_y = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_y.inds.resize(box);
        lending_y.inds.resize(box);
        auto const &borrowing_inds_y = borrowing_y.inds.array();
        auto const &lending_inds_y = lending_y.inds.array();
        auto const &Sy_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sy_enl = m_area_enl[maxLevel()][idim]->array(mfi);

        amrex::LoopOnCpu(box,
                         [=, &borrowing_y, &lending_y,
                             &temp_ind_borrowing_y, &temp_ind_lending_y](int i, int j, int k) {

            if( !flag_ext_face_y(i, j, k) ) return;

            amrex::Real Sy_stab = 0.5 * std::max({ly(i, j, k) * dz, ly(i, j, k + 1) * dz,
                                                  lz(i, j, k) * dy, lz(i, j, k + 1) * dy});
            amrex::Real Sy_ext = Sy_stab - Sy(i, j, k);
            for (int i_n = -1; i_n < 2; i_n++) {
                for (int k_n = -1; k_n < 2; k_n++) {
                    if(!(i_n == k_n or i_n == -k_n) ){
                        if (Sy_red(i + i_n, j, k + k_n) > Sy_ext
                        and flag_avail_face_y(i + i_n, j, k + k_n) and flag_ext_face_y(i, j, k)) {
                            Sy_red(i + i_n, j, k + k_n) -= Sy_ext;
                            borrowing_inds_y(i, j, k).push_back(temp_ind_borrowing_y);
                            temp_ind_borrowing_y = temp_ind_borrowing_y + 1;
                            borrowing_y.i_face.push_back(i + i_n);
                            borrowing_y.j_face.push_back(j);
                            borrowing_y.k_face.push_back(k + k_n);
                            borrowing_y.area.push_back(Sy_ext);
                            borrowing_y.rho_face.push_back(NAN);
                            lending_inds_y(i + i_n, j, k + k_n).push_back(temp_ind_lending_y);
                            temp_ind_lending_y = temp_ind_lending_y + 1;
                            lending_y.i_face.push_back(i);
                            lending_y.j_face.push_back(j);
                            lending_y.k_face.push_back(k);
                            lending_y.area.push_back(Sy_ext);
                            lending_y.rho_face.push_back(NAN);
                            Sy_enl(i, j, k) = Sy(i, j, k) + Sy_ext;
                            flag_ext_face_y(i, j, k) = false;
                        }
                    }
                }
            }
        });

        // Do the extensions in the z-plane
        idim = 2;
        auto const &Sz = m_face_areas[maxLevel()][idim]->array(mfi);
        auto const &flag_ext_face_z = m_flag_ext_face[maxLevel()][idim]->array(mfi);
        auto const &flag_avail_face_z = m_flag_avail_face[maxLevel()][idim]->array(mfi);
        auto &lending_z = (*m_lending[maxLevel()][idim])[mfi];
        auto &borrowing_z = (*m_borrowing[maxLevel()][idim])[mfi];
        borrowing_z.inds.resize(box);
        lending_z.inds.resize(box);
        auto const &borrowing_inds_z = borrowing_z.inds.array();
        auto const &lending_inds_z = lending_z.inds.array();
        auto const &Sz_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sz_enl = m_area_enl[maxLevel()][idim]->array(mfi);

        amrex::LoopOnCpu(amrex::convert(box, amrex::Box(Sz).ixType()),
                         [=, &borrowing_z, &lending_z,
                          &temp_ind_borrowing_z, &temp_ind_lending_z](int i, int j, int k) {

            if( !flag_ext_face_z(i, j, k) ) return;

            amrex::Real Sz_stab = 0.5 * std::max({lx(i, j, k) * dy, lx(i, j + 1, k) * dy,
                                                  ly(i, j, k) * dx, ly(i + 1, j, k) * dx});
            amrex::Real Sz_ext = Sz_stab - Sz(i, j, k);

            for (int i_n = -1; i_n < 2; i_n++) {
                for (int j_n = -1; j_n < 2; j_n++) {
                    if(!(i_n == j_n or i_n == -j_n)){
                        if (Sz_red(i + i_n, j + j_n, k) > Sz_ext
                            and flag_avail_face_z( i + i_n, j + j_n, k)
                            and flag_ext_face_z(i, j, k)) {
                            Sz_red(i + i_n, j + j_n, k) -= Sz_ext;
                            borrowing_inds_z(i, j, k).push_back(temp_ind_borrowing_z);
                            temp_ind_borrowing_z = temp_ind_borrowing_z + 1;
                            borrowing_z.i_face.push_back(i + i_n);
                            borrowing_z.j_face.push_back(j + j_n);
                            borrowing_z.k_face.push_back(k);
                            borrowing_z.area.push_back(Sz_ext);
                            borrowing_z.rho_face.push_back(NAN);
                            lending_inds_z(i + i_n, j + j_n, k).push_back(temp_ind_lending_z);
                            temp_ind_lending_z = temp_ind_lending_z + 1;
                            lending_z.i_face.push_back(i);
                            lending_z.j_face.push_back(j);
                            lending_z.k_face.push_back(k);
                            lending_z.area.push_back(Sz_ext);
                            lending_z.rho_face.push_back(NAN);
                            Sz_enl(i, j, k) = Sz(i, j, k) + Sz_ext;
                            flag_ext_face_z(i, j, k) = false;
                        }
                    }
                }
            }

        });
    }
    return {temp_ind_borrowing_x, temp_ind_lending_x,
            temp_ind_borrowing_y, temp_ind_lending_y,
            temp_ind_borrowing_z, temp_ind_lending_z};
#else
    return {0, 0, 0, 0, 0, 0};
#endif
}

/**
* \brief Do the 8-ways extension
*/
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
        borrowing_x.inds.resize(box);
        lending_x.inds.resize(box);
        auto const &lending_inds_x = lending_x.inds.array();
        auto const &borrowing_inds_x = borrowing_x.inds.array();
        auto const &Sx_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sx_enl = m_area_enl[maxLevel()][idim]->array(mfi);

        amrex::LoopOnCpu(box,
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
                            borrowing_inds_x(i, j, k).push_back(temp_ind_borrowing_x);
                            temp_ind_borrowing_x = temp_ind_borrowing_x + 1;
                            borrowing_x.i_face.push_back(i);
                            borrowing_x.j_face.push_back(j + j_n);
                            borrowing_x.k_face.push_back(k + k_n);
                            borrowing_x.area.push_back(Sx_ext);
                            borrowing_x.rho_face.push_back(NAN);
                            lending_inds_x(i, j + j_n, k  + k_n).push_back(temp_ind_lending_x);
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
        borrowing_y.inds.resize(box);
        lending_y.inds.resize(box);
        auto const &lending_inds_y = lending_y.inds.array();
        auto const &borrowing_inds_y = borrowing_y.inds.array();
        auto const &Sy_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sy_enl = m_area_enl[maxLevel()][idim]->array(mfi);

        amrex::LoopOnCpu(box,
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
                            borrowing_inds_y(i, j, k).push_back(temp_ind_borrowing_y);
                            temp_ind_borrowing_y = temp_ind_borrowing_y + 1;
                            borrowing_y.i_face.push_back(i + i_n);
                            borrowing_y.j_face.push_back(j);
                            borrowing_y.k_face.push_back(k + k_n);
                            borrowing_y.area.push_back(Sy_ext);
                            borrowing_y.rho_face.push_back(NAN);
                            lending_inds_y(i + i_n, j, k  + k_n).push_back(temp_ind_lending_y);
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
        borrowing_z.inds.resize(box);
        lending_z.inds.resize(box);
        auto const &lending_inds_z = lending_z.inds.array();
        auto const &borrowing_inds_z = borrowing_z.inds.array();
        auto const &Sz_red = m_area_red[maxLevel()][idim]->array(mfi);
        auto const &Sz_enl = m_area_enl[maxLevel()][idim]->array(mfi);

        amrex::LoopOnCpu(box,
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
                            borrowing_inds_z(i, j, k).push_back(temp_ind_borrowing_z);
                            temp_ind_borrowing_z = temp_ind_borrowing_z + 1;
                            borrowing_z.i_face.push_back(i + i_n);
                            borrowing_z.j_face.push_back(j + j_n);
                            borrowing_z.k_face.push_back(k);
                            borrowing_z.area.push_back(Sz_ext);
                            borrowing_z.rho_face.push_back(NAN);
                            lending_inds_z(i + i_n, j + j_n, k).push_back(temp_ind_lending_z);
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