/* Copyright 2019 David Grote, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_FIELD_DATA_H_
#define WARPX_SPECTRAL_FIELD_DATA_H_

#include "SpectralFieldData_fwd.H"

#include "AnyFFT.H"
#include "SpectralKSpace.H"
#include "Utils/WarpX_Complex.H"

#include <AMReX_BaseFab.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_FabArray.H>
#include <AMReX_IndexType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#include <vector>

// Declare type for spectral fields
using SpectralField = amrex::FabArray< amrex::BaseFab <Complex> >;

class SpectralFieldIndex
{
    public:

        /**
         * \brief Constructor of the class SpectralFieldIndex
         *
         * Set integer indices to access data in spectral space
         * and total number of fields to be stored.
         *
         * \param[in] update_with_rho whether rho is used in the field update equations
         * \param[in] time_averaging  whether the time averaging algorithm is used
         * \param[in] do_multi_J      whether the multi-J algorithm is used (hence two currents
         *                            computed at the beginning and the end of the time interval
         *                            instead of one current computed at half time)
         * \param[in] dive_cleaning   whether to use div(E) cleaning to account for errors in
         *                            Gauss law (new field F in the update equations)
         * \param[in] divb_cleaning   whether to use div(B) cleaning to account for errors in
         *                            div(B) = 0 law (new field G in the update equations)
         * \param[in] pml             whether the indices are used to access spectral data
         *                            for the PML spectral solver
         * \param[in] pml_rz          whether the indices are used to access spectral data
         *                            for the RZ PML spectral solver
         */
        SpectralFieldIndex (const bool update_with_rho,
                            const bool time_averaging,
                            const bool do_multi_J,
                            const bool dive_cleaning,
                            const bool divb_cleaning,
                            const bool pml,
                            const bool pml_rz = false);

        /**
         * \brief Default constructor
         */
        SpectralFieldIndex () = default;

        /**
         * \brief Default destructor
         */
        ~SpectralFieldIndex () = default;

    // Total number of fields that are actually allocated
    int n_fields;

    // Indices overwritten in the constructor, for the fields that are actually allocated
    // (index -1 will never be used, unless there is some bug in the code implementation,
    // which would result in a runtime crash due to out-of-bound accesses that can be detected
    // by running the code in DEBUG mode)

    // Always
    int Ex = -1, Ey = -1, Ez = -1;
    int Bx = -1, By = -1, Bz = -1;
    int Jx = -1, Jy = -1, Jz = -1;
    int rho_old = -1, rho_new = -1, divE = -1;

    // Time averaging
    int Ex_avg = -1, Ey_avg = -1, Ez_avg = -1;
    int Bx_avg = -1, By_avg = -1, Bz_avg = -1;

    // Multi-J, div(E) and div(B) cleaning
    int Jx_new = -1, Jy_new = -1, Jz_new = -1;
    int F = -1, G = -1;

    // PML
    int Exy = -1, Exz = -1, Eyx = -1, Eyz = -1, Ezx = -1, Ezy = -1;
    int Bxy = -1, Bxz = -1, Byx = -1, Byz = -1, Bzx = -1, Bzy = -1;

    // PML with div(E) and/or div(B) cleaning
    int Exx = -1, Eyy = -1, Ezz = -1, Bxx = -1, Byy = -1, Bzz = -1;
    int Fx  = -1, Fy  = -1, Fz  = -1, Gx  = -1, Gy  = -1, Gz  = -1;

    // PML RZ
    int Er_pml = -1, Et_pml = -1, Br_pml = -1, Bt_pml = -1;
};

/** \brief Class that stores the fields in spectral space, and performs the
 *  Fourier transforms between real space and spectral space
 */
class SpectralFieldData
{

    public:
        SpectralFieldData( const int lev,
                           const amrex::BoxArray& realspace_ba,
                           const SpectralKSpace& k_space,
                           const amrex::DistributionMapping& dm,
                           const int n_field_required,
                           const bool periodic_single_box);
        SpectralFieldData() = default; // Default constructor
        SpectralFieldData& operator=(SpectralFieldData&& field_data) = default;
        ~SpectralFieldData();

        void ForwardTransform (const int lev,
                               const amrex::MultiFab& mf, const int field_index,
                               const int i_comp);

        void BackwardTransform (const int lev, amrex::MultiFab& mf, const int field_index,
                                const amrex::IntVect& fill_guards, const int i_comp);

        // `fields` stores fields in spectral space, as multicomponent FabArray
        SpectralField fields;

    private:
        // tmpRealField and tmpSpectralField store fields
        // right before/after the Fourier transform
        SpectralField tmpSpectralField; // contains Complexs
        amrex::MultiFab tmpRealField; // contains Reals
        AnyFFT::FFTplans forward_plan, backward_plan;
        // Correcting "shift" factors when performing FFT from/to
        // a cell-centered grid in real space, instead of a nodal grid
        SpectralShiftFactor xshift_FFTfromCell, xshift_FFTtoCell,
                            zshift_FFTfromCell, zshift_FFTtoCell;
#if defined(WARPX_DIM_3D)
        SpectralShiftFactor yshift_FFTfromCell, yshift_FFTtoCell;
#endif

        bool m_periodic_single_box;
};

#endif // WARPX_SPECTRAL_FIELD_DATA_H_
