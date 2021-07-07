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

/** Index for the PML fields, when stored in spectral space,
 *  (n_fields is automatically set to the total number of fields)
 *  TODO How to include the diagonal components only when needed?
 */
struct SpectralPMLIndex {
  enum {Exx=0, Exy, Exz, Eyx, Eyy, Eyz, Ezx, Ezy, Ezz,
        Bxx  , Bxy, Bxz, Byx, Byy, Byz, Bzx, Bzy, Bzz, Fx, Fy, Fz, Gx, Gy, Gz, n_fields};
};

class SpectralFieldIndex
{
    public:

        /**
         * \brief Constructor
         */
        SpectralFieldIndex (const bool update_with_rho,
                            const bool time_averaging,
                            const bool J_linear_in_time,
                            const bool dive_cleaning,
                            const bool divb_cleaning);

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
    // (index -1 will never be used, unless there is some bug in the code implementation)
    int Ex = -1;
    int Ey = -1;
    int Ez = -1;
    int Bx = -1;
    int By = -1;
    int Bz = -1;
    int Jx = -1;
    int Jy = -1;
    int Jz = -1;
    int rho_old = -1;
    int rho_new = -1;
    int divE = -1;
    int Ex_avg = -1;
    int Ey_avg = -1;
    int Ez_avg = -1;
    int Bx_avg = -1;
    int By_avg = -1;
    int Bz_avg = -1;
    int Jx_new = -1;
    int Jy_new = -1;
    int Jz_new = -1;
    int F = -1;
    int G = -1;
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
                               const int i_comp, const amrex::IntVect& stag);
        AMREX_FORCE_INLINE
        void ForwardTransform (const int lev,
                               const amrex::MultiFab& mf, const int field_index, const int i_comp)
        {
            ForwardTransform(lev, mf, field_index, i_comp, mf.ixType().toIntVect());
        }

        void BackwardTransform (const int lev, amrex::MultiFab& mf, const int field_index, const int i_comp);

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
#if (AMREX_SPACEDIM==3)
        SpectralShiftFactor yshift_FFTfromCell, yshift_FFTtoCell;
#endif

        bool m_periodic_single_box;
};

#endif // WARPX_SPECTRAL_FIELD_DATA_H_
