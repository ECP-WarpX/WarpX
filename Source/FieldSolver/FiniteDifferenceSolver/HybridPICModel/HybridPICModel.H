/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Roelof Groenewald (TAE Technologies)
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_HYBRIDPICMODEL_H_
#define WARPX_HYBRIDPICMODEL_H_

#include "HybridPICModel_fwd.H"

#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"

#include <AMReX_Array.H>
#include <AMReX_REAL.H>


/**
 * \brief This class contains the parameters needed to evaluate hybrid field
 * solutions (kinetic ions with fluid electrons).
 */
class HybridPICModel
{
public:
    HybridPICModel (int nlevs_max); // constructor

    /** Read user-defined model parameters. Called in constructor. */
    void ReadParameters ();

    /** Allocate hybrid-PIC specific multifabs. Called in constructor. */
    void AllocateMFs (int nlevs_max);
    void AllocateLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                           int ncomps, const amrex::IntVect& ngJ, const amrex::IntVect& ngRho,
                           const amrex::IntVect& jx_nodal_flag, const amrex::IntVect& jy_nodal_flag,
                           const amrex::IntVect& jz_nodal_flag, const amrex::IntVect& rho_nodal_flag);

    /** Helper function to clear values from hybrid-PIC specific multifabs. */
    void ClearLevel (int lev);

    void InitData ();

    /**
     * \brief
     * Function to evaluate the external current expressions and populate the
     * external current multifab. Note the external current can be a function
     * of time and therefore this should be re-evaluated at every step.
     */
    void GetCurrentExternal (
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& edge_lengths
    );
    void GetCurrentExternal (
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
        int lev
    );

    /**
     * \brief
     * Function to calculate the total current based on Ampere's law while
     * neglecting displacement current (J = curl x B). Used in the Ohm's law
     * solver (kinetic-fluid hybrid model).
     *
     * \param[in] Bfield       Magnetic field from which the current is calculated.
     * \param[in] edge_lengths Length of cell edges taking embedded boundaries into account
     */
    void CalculateCurrentAmpere (
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Bfield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& edge_lengths
    );
    void CalculateCurrentAmpere (
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
        int lev
    );

    /**
     * \brief
     * Function to update the E-field using Ohm's law (hybrid-PIC model).
     */
    void HybridPICSolveE (
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Jfield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Bfield,
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> const& rhofield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& edge_lengths,
        bool include_resistivity_term);

    void HybridPICSolveE (
        std::array< std::unique_ptr<amrex::MultiFab>, 3>& Efield,
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Jfield,
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
        std::unique_ptr<amrex::MultiFab> const& rhofield,
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
        int lev, bool include_resistivity_term);

    void HybridPICSolveE (
        std::array< std::unique_ptr<amrex::MultiFab>, 3>& Efield,
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Jfield,
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
        std::unique_ptr<amrex::MultiFab> const& rhofield,
        std::array< std::unique_ptr<amrex::MultiFab>, 3> const& edge_lengths,
        int lev, PatchType patch_type, bool include_resistivity_term);

    void BfieldEvolveRK (
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Jfield,
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> const& rhofield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& edge_lengths,
        amrex::Real dt, DtType a_dt_type,
        amrex::IntVect ng, std::optional<bool> nodal_sync);

    void BfieldEvolveRK (
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Jfield,
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> const& rhofield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& edge_lengths,
        amrex::Real dt, int lev, DtType dt_type,
        amrex::IntVect ng, std::optional<bool> nodal_sync);

    void FieldPush (
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& Jfield,
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> const& rhofield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>> const& edge_lengths,
        amrex::Real dt, DtType dt_type,
        amrex::IntVect ng, std::optional<bool> nodal_sync);

    /**
     * \brief
     * Function to calculate the electron pressure at a given timestep type
     * using the simulation charge density. Used in the Ohm's law solver
     * (kinetic-fluid hybrid model).
     */
    void CalculateElectronPressure (          DtType a_dt_type);
    void CalculateElectronPressure (int lev,  DtType a_dt_type);

    /**
     * \brief Fill the electron pressure multifab given the kinetic particle
     * charge density (and assumption of quasi-neutrality) using the user
     * specified electron equation of state.
     *
     * \param[out] Pe_filed scalar electron pressure MultiFab at a given level
     * \param[in] rho_field scalar ion charge density Multifab at a given level
     */
    void FillElectronPressureMF (
        std::unique_ptr<amrex::MultiFab> const& Pe_field,
        amrex::MultiFab* const& rho_field ) const;

    // Declare variables to hold hybrid-PIC model parameters
    /** Number of substeps to take when evolving B */
    int m_substeps = 10;

    /** Electron temperature in eV */
    amrex::Real m_elec_temp;
    /** Reference electron density */
    amrex::Real m_n0_ref = 1.0;
    /** Electron pressure scaling exponent */
    amrex::Real m_gamma = 5.0/3.0;

    /** Plasma density floor - if n < n_floor it will be set to n_floor */
    amrex::Real m_n_floor = 1.0;

    /** Plasma resistivity */
    std::string m_eta_expression = "0.0";
    std::unique_ptr<amrex::Parser> m_resistivity_parser;
    amrex::ParserExecutor<2> m_eta;
    bool m_resistivity_has_J_dependence = false;

    /** External current */
    std::string m_Jx_ext_grid_function = "0.0";
    std::string m_Jy_ext_grid_function = "0.0";
    std::string m_Jz_ext_grid_function = "0.0";
    std::array< std::unique_ptr<amrex::Parser>, 3> m_J_external_parser;
    std::array< amrex::ParserExecutor<4>, 3> m_J_external;
    bool m_external_field_has_time_dependence = false;

    // Declare multifabs specifically needed for the hybrid-PIC model
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_fp_temp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp_temp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp_ampere;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp_external;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > electron_pressure_fp;

    // Helper functions to retrieve hybrid-PIC multifabs
    [[nodiscard]] amrex::MultiFab*
    get_pointer_current_fp_ampere  (int lev, int direction) const
    {
        return current_fp_ampere[lev][direction].get();
    }

    [[nodiscard]] amrex::MultiFab*
    get_pointer_electron_pressure_fp  (int lev) const
    {
        return electron_pressure_fp[lev].get();
    }

    /** Gpu Vector with index type of the Jx multifab */
    amrex::GpuArray<int, 3> Jx_IndexType;
    /** Gpu Vector with index type of the Jy multifab */
    amrex::GpuArray<int, 3> Jy_IndexType;
    /** Gpu Vector with index type of the Jz multifab */
    amrex::GpuArray<int, 3> Jz_IndexType;
    /** Gpu Vector with index type of the Bx multifab */
    amrex::GpuArray<int, 3> Bx_IndexType;
    /** Gpu Vector with index type of the By multifab */
    amrex::GpuArray<int, 3> By_IndexType;
    /** Gpu Vector with index type of the Bz multifab */
    amrex::GpuArray<int, 3> Bz_IndexType;
    /** Gpu Vector with index type of the Ex multifab */
    amrex::GpuArray<int, 3> Ex_IndexType;
    /** Gpu Vector with index type of the Ey multifab */
    amrex::GpuArray<int, 3> Ey_IndexType;
    /** Gpu Vector with index type of the Ez multifab */
    amrex::GpuArray<int, 3> Ez_IndexType;
};

/**
 * \brief
 * This struct contains only static functions to compute the electron pressure
 * using the particle density at a given point and the user provided reference
 * density and temperatures.
 */
struct ElectronPressure {

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static amrex::Real get_pressure (amrex::Real const n0,
                                     amrex::Real const T0,
                                     amrex::Real const gamma,
                                     amrex::Real const rho) {
        return n0 * T0 * std::pow((rho/PhysConst::q_e)/n0, gamma);
    }
};

#endif // WARPX_HYBRIDPICMODEL_H_
