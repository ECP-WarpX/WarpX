/* Copyright 2020 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_FINITE_DIFFERENCE_SOLVER_H_
#define WARPX_FINITE_DIFFERENCE_SOLVER_H_

#include <AMReX_MultiFab.H>
#include "MacroscopicProperties/MacroscopicProperties.H"
#include "BoundaryConditions/PML.H"

/**
 * \brief Top-level class for the electromagnetic finite-difference solver
 *
 * Stores the coefficients of the finite-difference stencils,
 * and has member functions to update fields over one time step.
 */
class FiniteDifferenceSolver
{
    public:

        // Constructor
        /** \brief Initialize the finite-difference Maxwell solver (for a given refinement level)
         *
         * This function initializes the stencil coefficients for the chosen finite-difference algorithm
         *
         * \param fdtd_algo Identifies the chosen algorithm, as defined in WarpXAlgorithmSelection.H
         * \param cell_size Cell size along each dimension, for the chosen refinement level
         * \param do_nodal  Whether the solver is applied to a nodal or staggered grid
         */
        FiniteDifferenceSolver (
            int const fdtd_algo,
            std::array<amrex::Real,3> cell_size,
            bool const do_nodal );

        void EvolveB ( std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
                       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
                       int lev, amrex::Real const dt );

        void EvolveE ( std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
                       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
                       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
                       std::unique_ptr<amrex::MultiFab> const& Ffield,
                       int lev, amrex::Real const dt );

        void EvolveF ( std::unique_ptr<amrex::MultiFab>& Ffield,
                       std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
                       std::unique_ptr<amrex::MultiFab> const& rhofield,
                       int const rhocomp,
                       amrex::Real const dt );

        void ApplySilverMuellerBoundary(
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
            amrex::Box domain_box,
            amrex::Real const dt );

        void ComputeDivE ( const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
                           amrex::MultiFab& divE );

        /**
          * \brief Macroscopic E-update for non-vacuum medium using the user-selected
          * finite-difference algorithm and macroscopic sigma-method defined in
          * WarpXAlgorithmSelection.H
          *
          * \param[out] Efield  vector of electric field MultiFabs updated at a given level
          * \param[in] Bfield   vector of magnetic field MultiFabs at a given level
          * \param[in] Jfield   vector of current density MultiFabs at a given level
          * \param[in] dt       timestep of the simulation
          * \param[in] macroscopic_properties contains user-defined properties of the medium.
          */

        void MacroscopicEvolveE ( std::array< std::unique_ptr<amrex::MultiFab>, 3>& Efield,
                            std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
                            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
                            amrex::Real const dt,
                            std::unique_ptr<MacroscopicProperties> const& macroscopic_properties);

        void EvolveBPML ( std::array< amrex::MultiFab*, 3 > Bfield,
                      std::array< amrex::MultiFab*, 3 > const Efield,
                      amrex::Real const dt,
                      const bool dive_cleaning);

       void EvolveEPML ( std::array< amrex::MultiFab*, 3 > Efield,
                      std::array< amrex::MultiFab*, 3 > const Bfield,
                      std::array< amrex::MultiFab*, 3 > const Jfield,
                      amrex::MultiFab* const Ffield,
                      MultiSigmaBox const& sigba,
                      amrex::Real const dt, bool pml_has_particles );

       void EvolveFPML ( amrex::MultiFab* Ffield,
                     std::array< amrex::MultiFab*, 3 > const Efield,
                     amrex::Real const dt );

    private:

        int m_fdtd_algo;
        bool m_do_nodal;

#ifdef WARPX_DIM_RZ
        amrex::Real m_dr, m_rmin;
        int m_nmodes;
        amrex::Gpu::DeviceVector<amrex::Real> m_stencil_coefs_r;
        amrex::Gpu::DeviceVector<amrex::Real> m_stencil_coefs_z;
#else
        amrex::Gpu::DeviceVector<amrex::Real> m_stencil_coefs_x;
        amrex::Gpu::DeviceVector<amrex::Real> m_stencil_coefs_y;
        amrex::Gpu::DeviceVector<amrex::Real> m_stencil_coefs_z;
#endif

    public:
        // The member functions below contain extended __device__ lambda.
        // In order to compile with nvcc, they need to be public.

#ifdef WARPX_DIM_RZ
        template< typename T_Algo >
        void EvolveBCylindrical (
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
            const int lev,
            amrex::Real const dt );

        template< typename T_Algo >
        void EvolveECylindrical (
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
            std::unique_ptr<amrex::MultiFab> const& Ffield,
            const int lev,
            amrex::Real const dt );

        template< typename T_Algo >
        void EvolveFCylindrical (
            std::unique_ptr<amrex::MultiFab>& Ffield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
            std::unique_ptr<amrex::MultiFab> const& rhofield,
            int const rhocomp,
            amrex::Real const dt );

        template< typename T_Algo >
        void ComputeDivECylindrical (
            const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
            amrex::MultiFab& divE );

#else
        template< typename T_Algo >
        void EvolveBCartesian (
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Bfield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
            int lev, amrex::Real const dt );

        template< typename T_Algo >
        void EvolveECartesian (
            std::array< std::unique_ptr<amrex::MultiFab>, 3 >& Efield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Bfield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
            std::unique_ptr<amrex::MultiFab> const& Ffield,
            int lev, amrex::Real const dt );

        template< typename T_Algo >
        void EvolveFCartesian (
            std::unique_ptr<amrex::MultiFab>& Ffield,
            std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Efield,
            std::unique_ptr<amrex::MultiFab> const& rhofield,
            int const rhocomp,
            amrex::Real const dt );

        template< typename T_Algo >
        void ComputeDivECartesian (
            const std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
            amrex::MultiFab& divE );

        template< typename T_Algo, typename T_MacroAlgo >
        void MacroscopicEvolveECartesian (
            std::array< std::unique_ptr< amrex::MultiFab>, 3>& Efield,
            std::array< std::unique_ptr< amrex::MultiFab>, 3> const &Bfield,
            std::array< std::unique_ptr< amrex::MultiFab>, 3> const& Jfield,
            amrex::Real const dt,
            std::unique_ptr<MacroscopicProperties> const& macroscopic_properties);

        template< typename T_Algo >
        void EvolveBPMLCartesian (
            std::array< amrex::MultiFab*, 3 > Bfield,
            std::array< amrex::MultiFab*, 3 > const Efield,
            amrex::Real const dt,
            const bool dive_cleaning);

        template< typename T_Algo >
        void EvolveEPMLCartesian (
            std::array< amrex::MultiFab*, 3 > Efield,
            std::array< amrex::MultiFab*, 3 > const Bfield,
            std::array< amrex::MultiFab*, 3 > const Jfield,
            amrex::MultiFab* const Ffield,
            MultiSigmaBox const& sigba,
            amrex::Real const dt, bool pml_has_particles );

        template< typename T_Algo >
        void EvolveFPMLCartesian ( amrex::MultiFab* Ffield,
                      std::array< amrex::MultiFab*, 3 > const Efield,
                      amrex::Real const dt );

#endif

};

#endif // WARPX_FINITE_DIFFERENCE_SOLVER_H_
