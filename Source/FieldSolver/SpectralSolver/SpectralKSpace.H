/* Copyright 2019 David Grote, Maxence Thevenet, Remi Lehe
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_SPECTRAL_K_SPACE_H_
#define WARPX_SPECTRAL_K_SPACE_H_

#include "Utils/WarpX_Complex.H"

#include <AMReX_BoxArray.H>
#include <AMReX_LayoutData.H>


// `KVectorComponent` and `SpectralShiftFactor` hold one 1D array
// ("DeviceVector") for each box ("LayoutData"). The arrays are
// only allocated if the corresponding box is owned by the local MPI rank.
using RealKVector = amrex::Gpu::DeviceVector<amrex::Real>;
using KVectorComponent = amrex::LayoutData< RealKVector >;
using SpectralShiftFactor = amrex::LayoutData<
                           amrex::Gpu::DeviceVector<Complex> >;

// Indicate the type of correction "shift" factor to apply
// when the FFT is performed from/to a cell-centered grid in real space.
struct ShiftType {
    enum{ TransformFromCellCentered=0, TransformToCellCentered=1 };
};

/**
 * \brief Class that represents the spectral space.
 *
 * (Contains info about the size of the spectral space corresponding
 * to each box in `realspace_ba`, as well as the value of the
 * corresponding k coordinates)
 */
class SpectralKSpace
{
    public:
        amrex::BoxArray spectralspace_ba;
        SpectralKSpace() : dx(amrex::RealVect::Zero) {}
        SpectralKSpace( const amrex::BoxArray& realspace_ba,
                        const amrex::DistributionMapping& dm,
                        const amrex::RealVect realspace_dx );
        KVectorComponent getKComponent(
            const amrex::DistributionMapping& dm,
            const amrex::BoxArray& realspace_ba,
            const int i_dim, const bool only_positive_k ) const;
        KVectorComponent getModifiedKComponent(
            const amrex::DistributionMapping& dm, const int i_dim,
            const int n_order, const bool nodal ) const;
        SpectralShiftFactor getSpectralShiftFactor(
            const amrex::DistributionMapping& dm, const int i_dim,
            const int shift_type ) const;

    protected:
        amrex::Array<KVectorComponent, AMREX_SPACEDIM> k_vec;
        // 3D: k_vec is an Array of 3 components, corresponding to kx, ky, kz
        // 2D: k_vec is an Array of 2 components, corresponding to kx, kz
        amrex::RealVect dx;
};

/**
 * \brief Returns an array of coefficients (Fornberg coefficients), corresponding
 * to the weight of each point in a finite-difference approximation of a derivative
 * (up to order \c n_order).
 *
 * \param[in] n_order order of the finite-difference approximation
 * \param[in] nodal   whether the finite-difference approximation is computed
 *                    on a nodal grid or a staggered grid
 */
amrex::Vector<amrex::Real>
getFornbergStencilCoefficients(const int n_order, const bool nodal);

#endif
