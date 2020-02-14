/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include <SpectralKSpace.H>
#include <SpectralSolver.H>
#include <PsatdAlgorithm.H>
#include <PMLPsatdAlgorithm.H>

using namespace amrex;

/**
 * \brief Initialize the spectral Maxwell solver
 *
 * This function selects the spectral algorithm to be used, allocates the
 * corresponding coefficients for the discretized field update equation,
 * and prepares the structures that store the fields in spectral space.
 *
 * \param norder_x Order of accuracy of the spatial derivatives along x
 * \param norder_y Order of accuracy of the spatial derivatives along y
 * \param norder_z Order of accuracy of the spatial derivatives along z
 * \param nodal    Whether the solver is applied to a nodal or staggered grid
 * \param dx       Cell size along each dimension
 * \param dt       Time step
 * \param pml      Whether the boxes in which the solver is applied are PML boxes
 */
SpectralSolver::SpectralSolver(
                const amrex::BoxArray& realspace_ba,
                const amrex::DistributionMapping& dm,
                const int norder_x, const int norder_y,
                const int norder_z, const bool nodal,
                const amrex::RealVect dx, const amrex::Real dt,
                const bool pml ) {

    // Initialize all structures using the same distribution mapping dm

    // - Initialize k space object (Contains info about the size of
    // the spectral space corresponding to each box in `realspace_ba`,
    // as well as the value of the corresponding k coordinates)
    const SpectralKSpace k_space= SpectralKSpace(realspace_ba, dm, dx);

    // - Select the algorithm depending on the input parameters
    //   Initialize the corresponding coefficients over k space
    if (pml) {
        algorithm = std::unique_ptr<PMLPsatdAlgorithm>( new PMLPsatdAlgorithm(
            k_space, dm, norder_x, norder_y, norder_z, nodal, dt ) );
    } else {
        algorithm = std::unique_ptr<PsatdAlgorithm>( new PsatdAlgorithm(
            k_space, dm, norder_x, norder_y, norder_z, nodal, dt ) );
    }

    m_realspace_ba = realspace_ba;
    m_dm = dm;
    m_norder_x = norder_x;
    m_norder_y = norder_y;
    m_norder_z = norder_z;
    m_nodal = nodal;
    m_dx = dx;

    // Initialize arrays for fields in spectral space + FFT plans
    field_data = SpectralFieldData( realspace_ba, k_space, dm, algorithm->getRequiredNumberOfFields() );
};

/**
 * Compute divergence of E
 */
void
SpectralSolver::ComputeSpectralDivE( std::unique_ptr<amrex::MultiFab>& divE,
                                     const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield ) {

    using Idx = SpectralFieldIndex;

    // Forward Fourier transform of E
    ForwardTransform( *Efield[0], Idx::Ex );
    ForwardTransform( *Efield[1], Idx::Ey );
    ForwardTransform( *Efield[2], Idx::Ez );

    const SpectralKSpace k_space= SpectralKSpace( m_realspace_ba, m_dm, m_dx );

    KVectorComponent modified_kx_vec = k_space.getModifiedKComponent( m_dm, 0, m_norder_x, m_nodal );
#if (AMREX_SPACEDIM==3)
    KVectorComponent modified_ky_vec = k_space.getModifiedKComponent( m_dm, 1, m_norder_y, m_nodal );
    KVectorComponent modified_kz_vec = k_space.getModifiedKComponent( m_dm, 2, m_norder_z, m_nodal );
#else
    KVectorComponent modified_kz_vec = k_space.getModifiedKComponent( m_dm, 1, m_norder_z, m_nodal );
#endif

    // Loop over boxes
    for (MFIter mfi(field_data.fields); mfi.isValid(); ++mfi){

        const Box& bx = field_data.fields[mfi].box();

        // Extract arrays for the fields to be updated
        Array4<Complex> fields = field_data.fields[mfi].array();
        // Extract pointers for the k vectors
        const Real* modified_kx_arr = modified_kx_vec[mfi].dataPtr();
#if (AMREX_SPACEDIM==3)
        const Real* modified_ky_arr = modified_ky_vec[mfi].dataPtr();
#endif
        const Real* modified_kz_arr = modified_kz_vec[mfi].dataPtr();

        // Loop over indices within one box
        ParallelFor(bx,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            using Idx = SpectralFieldIndex;
            // Shortcuts for the components of E
            const Complex Ex = fields(i,j,k,Idx::Ex);
            const Complex Ey = fields(i,j,k,Idx::Ey);
            const Complex Ez = fields(i,j,k,Idx::Ez);
            // k vector values
            const Real kx = modified_kx_arr[i];
#if (AMREX_SPACEDIM==3)
            const Real ky = modified_ky_arr[j];
            const Real kz = modified_kz_arr[k];
#else
            constexpr Real ky = 0;
            const Real kz = modified_kz_arr[j];
#endif
            const Complex I = Complex{0,1};

            // div(E) in Fourier space
            fields(i,j,k,Idx::divE) = I*(kx*Ex+ky*Ey+kz*Ez);
        });
    }

    // Backward Fourier transform
    BackwardTransform( *divE, Idx::divE );
};
