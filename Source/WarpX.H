/* Copyright 2016-2020 Andrew Myers, Ann Almgren, Aurore Blelly
 *                     Axel Huebl, Burlen Loring, David Grote
 *                     Glenn Richardson, Junmin Gu, Luca Fedeli
 *                     Mathieu Lobet, Maxence Thevenet, Michael Rowan
 *                     Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 *                     Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_H_
#define WARPX_H_

#include "Evolve/WarpXDtType.H"
#include "Particles/MultiParticleContainer.H"
#include "BoundaryConditions/PML.H"
#include "Diagnostics/BackTransformedDiagnostic.H"
#include "Diagnostics/MultiDiagnostics.H"
#include "Filter/BilinearFilter.H"
#include "Filter/NCIGodfreyFilter.H"
#include "Diagnostics/ReducedDiags/MultiReducedDiags.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/IntervalsParser.H"
#include "FieldSolver/FiniteDifferenceSolver/MacroscopicProperties/MacroscopicProperties.H"

#include "FieldSolver/FiniteDifferenceSolver/FiniteDifferenceSolver.H"
#ifdef WARPX_USE_PSATD
#   ifdef WARPX_DIM_RZ
#       include "FieldSolver/SpectralSolver/SpectralSolverRZ.H"
#   else
#       include "FieldSolver/SpectralSolver/SpectralSolver.H"
#   endif
#endif

#include "Parallelization/GuardCellManager.H"

#ifdef WARPX_USE_OPENPMD
#   include "Diagnostics/WarpXOpenPMD.H"
#endif

#include "Parser/WarpXParserWrapper.H"

#include <AMReX_AmrCore.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_Print.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealVect.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_LayoutData.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>

#ifdef AMREX_USE_OMP
#   include <omp.h>
#endif

#include <iostream>
#include <memory>
#include <array>

enum struct PatchType : int
{
    fine,
    coarse
};

class WarpX
    : public amrex::AmrCore
{
public:

    friend class PML;

    static WarpX& GetInstance ();
    static void ResetInstance ();

    WarpX ();
    ~WarpX ();

    static std::string Version (); //!< Version of WarpX executable
    static std::string PicsarVersion (); //!< Version of PICSAR dependency

    int Verbose () const { return verbose; }


    void InitData ();

    void Evolve (int numsteps = -1);

    MultiParticleContainer& GetPartContainer () { return *mypc; }

    static void shiftMF (amrex::MultiFab& mf, const amrex::Geometry& geom,
                         int num_shift, int dir, amrex::Real external_field=0.0,
                         bool useparser = false, HostDeviceParser<3> const& field_parser={});

    static void GotoNextLine (std::istream& is);

    //! Author of an input file / simulation setup
    static std::string authors;

    // Initial field on the grid.
    static amrex::Vector<amrex::Real> E_external_grid;
    static amrex::Vector<amrex::Real> B_external_grid;

    // Initialization Type for External E and B on grid
    static std::string B_ext_grid_s;
    static std::string E_ext_grid_s;

    // Parser for B_external on the grid
    static std::string str_Bx_ext_grid_function;
    static std::string str_By_ext_grid_function;
    static std::string str_Bz_ext_grid_function;
    // Parser for E_external on the grid
    static std::string str_Ex_ext_grid_function;
    static std::string str_Ey_ext_grid_function;
    static std::string str_Ez_ext_grid_function;

    // ParserWrapper for B_external on the grid
    std::unique_ptr<ParserWrapper<3> > Bxfield_parser;
    std::unique_ptr<ParserWrapper<3> > Byfield_parser;
    std::unique_ptr<ParserWrapper<3> > Bzfield_parser;
    // ParserWrapper for E_external on the grid
    std::unique_ptr<ParserWrapper<3> > Exfield_parser;
    std::unique_ptr<ParserWrapper<3> > Eyfield_parser;
    std::unique_ptr<ParserWrapper<3> > Ezfield_parser;

    // Algorithms
    static long current_deposition_algo;
    static long charge_deposition_algo;
    static long field_gathering_algo;
    static long particle_pusher_algo;
    static int maxwell_solver_id;
    static long load_balance_costs_update_algo;
    static int em_solver_medium;
    static int macroscopic_solver_algo;
    static amrex::Vector<int> field_boundary_lo;
    static amrex::Vector<int> field_boundary_hi;
    static amrex::Vector<int> particle_boundary_lo;
    static amrex::Vector<int> particle_boundary_hi;


    // If true, the current is deposited on a nodal grid and then centered onto a staggered grid
    static bool do_current_centering;

    // PSATD: If true (overwritten by the user in the input file), the current correction
    // defined in equation (19) of https://doi.org/10.1016/j.jcp.2013.03.010 is applied
    bool current_correction = false;

    // PSATD: If true, the update equation for E contains both J and rho (at times n and n+1):
    // default is false for standard PSATD and true for Galilean PSATD (set in WarpX.cpp)
    bool update_with_rho = false;

    // div(E) cleaning
    static bool do_dive_cleaning;

    // Interpolation order
    static long nox;
    static long noy;
    static long noz;

    // Order of finite-order centering of fields (staggered to nodal)
    static int field_centering_nox;
    static int field_centering_noy;
    static int field_centering_noz;

    // Order of finite-order centering of currents (nodal to staggered)
    static int current_centering_nox;
    static int current_centering_noy;
    static int current_centering_noz;

    // Number of modes for the RZ multimode version
    static int n_rz_azimuthal_modes;
    static int ncomps;

    static bool use_fdtd_nci_corr;
    static bool galerkin_interpolation;

    static bool use_filter;
    static bool use_kspace_filter;
    static bool use_filter_compensation;
    static bool use_damp_fields_in_z_guard;
    static bool serialize_ics;

    // Back transformation diagnostic
    static bool do_back_transformed_diagnostics;
    static std::string lab_data_directory;
    static int  num_snapshots_lab;
    static amrex::Real dt_snapshots_lab;
    static bool do_back_transformed_fields;
    static bool do_back_transformed_particles;

    // Boosted frame parameters
    static amrex::Real gamma_boost;
    static amrex::Real beta_boost;
    static amrex::Vector<int> boost_direction;
    static amrex::Real zmax_plasma_to_compute_max_step;
    static int do_compute_max_step_from_zmax;

    static bool do_dynamic_scheduling;
    static bool refine_plasma;

    static IntervalsParser sort_intervals;
    static amrex::IntVect sort_bin_size;

    static int do_subcycling;

    static bool do_device_synchronize_before_profile;
    static bool safe_guard_cells;

    // buffers
    static int n_field_gather_buffer;       //! in number of cells from the edge (identical for each dimension)
    static int n_current_deposition_buffer; //! in number of cells from the edge (identical for each dimension)

    // do nodal
    static int do_nodal;

    std::array<const amrex::MultiFab* const, 3>
    get_array_Bfield_aux  (const int lev) const {
        return {
            Bfield_aux[lev][0].get(),
            Bfield_aux[lev][1].get(),
            Bfield_aux[lev][2].get()
        };
    }
    std::array<const amrex::MultiFab* const, 3>
    get_array_Efield_aux  (const int lev) const {
        return {
            Efield_aux[lev][0].get(),
            Efield_aux[lev][1].get(),
            Efield_aux[lev][2].get()
        };
    }

    amrex::MultiFab * get_pointer_Efield_aux  (int lev, int direction) const { return Efield_aux[lev][direction].get(); }
    amrex::MultiFab * get_pointer_Bfield_aux  (int lev, int direction) const { return Bfield_aux[lev][direction].get(); }

    amrex::MultiFab * get_pointer_Efield_fp  (int lev, int direction) const { return Efield_fp[lev][direction].get(); }
    amrex::MultiFab * get_pointer_Bfield_fp  (int lev, int direction) const { return Bfield_fp[lev][direction].get(); }
    amrex::MultiFab * get_pointer_current_fp  (int lev, int direction) const { return current_fp[lev][direction].get(); }
    amrex::MultiFab * get_pointer_rho_fp  (int lev) const { return rho_fp[lev].get(); }
    amrex::MultiFab * get_pointer_F_fp  (int lev) const { return F_fp[lev].get(); }
    amrex::MultiFab * get_pointer_phi_fp  (int lev) const { return phi_fp[lev].get(); }

    amrex::MultiFab * get_pointer_Efield_cp  (int lev, int direction) const { return Efield_cp[lev][direction].get(); }
    amrex::MultiFab * get_pointer_Bfield_cp  (int lev, int direction) const { return Bfield_cp[lev][direction].get(); }
    amrex::MultiFab * get_pointer_current_cp  (int lev, int direction) const { return current_cp[lev][direction].get(); }
    amrex::MultiFab * get_pointer_rho_cp  (int lev) const { return rho_cp[lev].get(); }
    amrex::MultiFab * get_pointer_F_cp  (int lev) const { return F_cp[lev].get(); }

    const amrex::MultiFab& getcurrent (int lev, int direction) {return *current_fp[lev][direction];}
    const amrex::MultiFab& getEfield  (int lev, int direction) {return *Efield_aux[lev][direction];}
    const amrex::MultiFab& getBfield  (int lev, int direction) {return *Bfield_aux[lev][direction];}

    const amrex::MultiFab& getcurrent_cp (int lev, int direction) {return *current_cp[lev][direction];}
    const amrex::MultiFab& getEfield_cp  (int lev, int direction) {return  *Efield_cp[lev][direction];}
    const amrex::MultiFab& getBfield_cp  (int lev, int direction) {return  *Bfield_cp[lev][direction];}
    const amrex::MultiFab& getrho_cp (int lev) {return  *rho_cp[lev];}

    const amrex::MultiFab& getcurrent_fp (int lev, int direction) {return *current_fp[lev][direction];}
    const amrex::MultiFab& getEfield_fp  (int lev, int direction) {return *Efield_fp[lev][direction];}
    const amrex::MultiFab& getBfield_fp  (int lev, int direction) {return *Bfield_fp[lev][direction];}
    const amrex::MultiFab& getrho_fp (int lev) {return *rho_fp[lev];}
    const amrex::MultiFab& getphi_fp (int lev) {return *phi_fp[lev];}
    const amrex::MultiFab& getF_fp (int lev) {return *F_fp[lev];}

    const amrex::MultiFab& getEfield_avg_fp (int lev, int direction) {return *Efield_avg_fp[lev][direction];}
    const amrex::MultiFab& getBfield_avg_fp (int lev, int direction) {return *Bfield_avg_fp[lev][direction];}

    bool DoPML () const {return do_pml;}

    /** get low-high-low-high-... vector for each direction indicating if mother grid PMLs are enabled */
    std::vector<bool> getPMLdirections() const;

    static amrex::LayoutData<amrex::Real>* getCosts (int lev);

    void setLoadBalanceEfficiency (const int lev, const amrex::Real efficiency)
    {
        if (m_instance)
        {
            m_instance->load_balance_efficiency[lev] = efficiency;
        } else
        {
            return;
        }
    }

    amrex::Real getLoadBalanceEfficiency (const int lev)
    {
        if (m_instance)
        {
            return m_instance->load_balance_efficiency[lev];
        } else
        {
            return -1;
        }
    }

    static amrex::IntVect filter_npass_each_dir;
    BilinearFilter bilinear_filter;
    amrex::Vector< std::unique_ptr<NCIGodfreyFilter> > nci_godfrey_filter_exeybz;
    amrex::Vector< std::unique_ptr<NCIGodfreyFilter> > nci_godfrey_filter_bxbyez;

    amrex::Real time_of_last_gal_shift  = 0;
    amrex::Array<amrex::Real,3> m_v_galilean = {{0}};
    amrex::Array<amrex::Real,3> m_galilean_shift = {{0}};

    amrex::Array<amrex::Real,3> m_v_comoving = {{0.}};

    static int num_mirrors;
    amrex::Vector<amrex::Real> mirror_z;
    amrex::Vector<amrex::Real> mirror_z_width;
    amrex::Vector<int> mirror_z_npoints;

    /// object with all reduced diagnotics, similar to MultiParticleContainer for species.
    MultiReducedDiags* reduced_diags;

    void applyMirrors(amrex::Real time);

    /** Determine the timestep of the simulation. */
    void ComputeDt ();

    // Compute max_step automatically for simulations in a boosted frame.
    void computeMaxStepBoostAccelerator(const amrex::Geometry& geom);
    int  MoveWindow (bool move_j);

    /**
     * \brief
     * This function shifts the boundary of the grid by 'm_v_galilean*dt'.
     * In doding so, only positions attributes are changed while fields remain unchanged.
     */
    void ShiftGalileanBoundary ();
    void UpdatePlasmaInjectionPosition (amrex::Real dt);
    void ResetProbDomain (const amrex::RealBox& rb);
    void EvolveE (         amrex::Real dt);
    void EvolveE (int lev, amrex::Real dt);
    void EvolveB (         amrex::Real dt);
    void EvolveB (int lev, amrex::Real dt);
    void EvolveF (         amrex::Real dt, DtType dt_type);
    void EvolveF (int lev, amrex::Real dt, DtType dt_type);
    void EvolveB (int lev, PatchType patch_type, amrex::Real dt);
    void EvolveE (int lev, PatchType patch_type, amrex::Real dt);
    void EvolveF (int lev, PatchType patch_type, amrex::Real dt, DtType dt_type);
    void ApplySilverMuellerBoundary (amrex::Real dt);

    void MacroscopicEvolveE (         amrex::Real dt);
    void MacroscopicEvolveE (int lev, amrex::Real dt);
    void MacroscopicEvolveE (int lev, PatchType patch_type, amrex::Real dt);

    /** \brief apply QED correction on electric field
     * \param dt vector of time steps (for all levels)
     */
    void Hybrid_QED_Push (         amrex::Vector<amrex::Real> dt);

    /** \brief apply QED correction on electric field for level lev
     * \param lev mesh refinement level
     * \param dt time step
     */
    void Hybrid_QED_Push (int lev, amrex::Real dt);

    /** \brief apply QED correction on electric field for level lev and patch type patch_type
     * \param lev mesh refinement level
     * \param dt patch_type which MR patch: PatchType::fine or PatchType::coarse
     * \param dt time step
     */
    void Hybrid_QED_Push (int lev, PatchType patch_type, amrex::Real dt);

    static amrex::Real quantum_xi_c2;

    /** \brief perform load balance; compute and communicate new `amrex::DistributionMapping`
     */
    void LoadBalance ();
    /** \brief resets costs to zero
     */
    void ResetCosts ();

    /** \brief returns the load balance interval
     */
    IntervalsParser get_load_balance_intervals () const {return load_balance_intervals;}

    /**
     * \brief Private function for spectral solver
     * Applies a damping factor in the guards cells that extend
     * beyond the extent of the domain, reducing fluctuations that
     * can appear in parallel simulations. This will be called
     * when use_damp_fields_in_z_guard is true.
     */
    void DampFieldsInGuards (std::array<std::unique_ptr<amrex::MultiFab>,3>& Efield,
                             std::array<std::unique_ptr<amrex::MultiFab>,3>& Bfield);

#ifdef WARPX_DIM_RZ
    void ApplyInverseVolumeScalingToCurrentDensity(amrex::MultiFab* Jx,
                                                   amrex::MultiFab* Jy,
                                                   amrex::MultiFab* Jz,
                                                   int lev);

    void ApplyInverseVolumeScalingToChargeDensity(amrex::MultiFab* Rho,
                                                  int lev);
#endif

    void DampPML ();
    void DampPML (int lev);
    void DampPML (int lev, PatchType patch_type);

    void DampJPML ();
    void DampJPML (int lev);
    void DampJPML (int lev, PatchType patch_type);

    void CopyJPML ();

    /**
     * \brief Synchronize the nodal points of the PML MultiFabs
     */
    void NodalSyncPML ();

    /**
     * \brief Synchronize the nodal points of the PML MultiFabs for given MR level
     */
    void NodalSyncPML (int lev);

    /**
     * \brief Synchronize the nodal points of the PML MultiFabs for given MR level and patch
     */
    void NodalSyncPML (int lev, PatchType patch_type);

    PML* GetPML (int lev);

    /** Run the ionization module on all species */
    void doFieldIonization ();
    /** Run the ionization module on all species at level lev
     * \param lev level
     */
    void doFieldIonization (int lev);

#ifdef WARPX_QED
    /** Run the QED module on all species */
    void doQEDEvents ();
    /** Run the QED module on all species at level lev
     * \param lev level
     */
    void doQEDEvents (int lev);
#endif

    void PushParticlesandDepose (int lev, amrex::Real cur_time, DtType a_dt_type=DtType::Full, bool skip_current=false);
    void PushParticlesandDepose (         amrex::Real cur_time, bool skip_current=false);

    // This function does aux(lev) = fp(lev) + I(aux(lev-1)-cp(lev)).
    // Caller must make sure fp and cp have ghost cells filled.
    void UpdateAuxilaryData ();
    void UpdateAuxilaryDataStagToNodal ();
    void UpdateAuxilaryDataSameType ();

    /**
     * \brief This function is called if \c warpx.do_current_centering = 1 and
     * it centers the currents from a nodal grid to a staggered grid (Yee) using
     * finite-order interpolation based on the Fornberg coefficients.
     *
     * \param[in,out] dst destination \c MultiFab where the results of the finite-order centering are stored
     * \param[in] src source \c MultiFab that contains the values of the nodal current to be centered
     */
    void UpdateCurrentNodalToStag (amrex::MultiFab& dst, amrex::MultiFab const& src);

    // Fill boundary cells including coarse/fine boundaries
    void FillBoundaryB   (amrex::IntVect ng);
    void FillBoundaryE   (amrex::IntVect ng);
    void FillBoundaryB_avg   (amrex::IntVect ng);
    void FillBoundaryE_avg   (amrex::IntVect ng);

    void FillBoundaryF   (amrex::IntVect ng);
    void FillBoundaryAux (amrex::IntVect ng);
    void FillBoundaryE   (int lev, amrex::IntVect ng);
    void FillBoundaryB   (int lev, amrex::IntVect ng);
    void FillBoundaryE_avg   (int lev, amrex::IntVect ng);
    void FillBoundaryB_avg   (int lev, amrex::IntVect ng);

    void FillBoundaryF   (int lev, amrex::IntVect ng);
    void FillBoundaryAux (int lev, amrex::IntVect ng);

    void SyncCurrent ();
    void SyncRho ();

    amrex::Vector<int> getnsubsteps () const {return nsubsteps;}
    int getnsubsteps (int lev) const {return nsubsteps[lev];}
    amrex::Vector<int> getistep () const {return istep;}
    int getistep (int lev) const {return istep[lev];}
    void setistep (int lev, int ii) {istep[lev] = ii;}
    amrex::Vector<amrex::Real> gett_old () const {return t_old;}
    amrex::Real gett_old (int lev) const {return t_old[lev];}
    amrex::Vector<amrex::Real> gett_new () const {return t_new;}
    amrex::Real gett_new (int lev) const {return t_new[lev];}
    void sett_new (int lev, amrex::Real time) {t_new[lev] = time;}
    amrex::Vector<amrex::Real> getdt () const {return dt;}
    amrex::Real getdt (int lev) const {return dt[lev];}
    amrex::Real getmoving_window_x() const {return moving_window_x;}
    amrex::Real getcurrent_injection_position () const {return current_injection_position;}
    bool getis_synchronized() const {return is_synchronized;}
    void setplot_rho (bool a_plot_rho) {plot_rho = a_plot_rho;}

    int maxStep () const {return max_step;}
    amrex::Real stopTime () const {return stop_time;}

    void AverageAndPackFields( amrex::Vector<std::string>& varnames,
        amrex::Vector<amrex::MultiFab>& mf_avg, const amrex::IntVect ngrow) const;

    void prepareFields( int const step, amrex::Vector<std::string>& varnames,
        amrex::Vector<amrex::MultiFab>& mf_avg,
        amrex::Vector<const amrex::MultiFab*>& output_mf,
        amrex::Vector<amrex::Geometry>& output_geom ) const;

    static std::array<amrex::Real,3> CellSize (int lev);
    static amrex::RealBox getRealBox(const amrex::Box& bx, int lev);
    static std::array<amrex::Real,3> LowerCorner (const amrex::Box& bx,
                                                  std::array<amrex::Real,3> galilean_shift, int lev);
    static std::array<amrex::Real,3> UpperCorner (const amrex::Box& bx, int lev);

    /*
      /brief This computes the lower of the problem domain, taking into account any shift when using the Galilean algorithm.
     */
    std::array<amrex::Real,3> LowerCornerWithGalilean (const amrex::Box& bx, const amrex::Array<amrex::Real,3>& v_galilean, int lev);

    static amrex::IntVect RefRatio (int lev);

    static const amrex::iMultiFab* CurrentBufferMasks (int lev);
    static const amrex::iMultiFab* GatherBufferMasks (int lev);

    static int do_electrostatic;

    // Parameters for lab frame electrostatic
    static amrex::Real self_fields_required_precision;
    static int self_fields_max_iters;

    static int do_moving_window;
    static int moving_window_dir;
    static amrex::Real moving_window_v;
    static bool fft_do_time_averaging;

    // slice generation //
    static int num_slice_snapshots_lab;
    static amrex::Real dt_slice_snapshots_lab;
    static amrex::Real particle_slice_width_lab;
    amrex::RealBox getSliceRealBox() const {return slice_realbox;}

    // these should be private, but can't due to Cuda limitations
    static void ComputeDivB (amrex::MultiFab& divB, int const dcomp,
                             const std::array<const amrex::MultiFab* const, 3>& B,
                             const std::array<amrex::Real,3>& dx);

    static void ComputeDivB (amrex::MultiFab& divB, int const dcomp,
                             const std::array<const amrex::MultiFab* const, 3>& B,
                             const std::array<amrex::Real,3>& dx, amrex::IntVect const ngrow);

    void ComputeDivE(amrex::MultiFab& divE, const int lev);

    const amrex::IntVect getngE() const { return guard_cells.ng_alloc_EB; }
    const amrex::IntVect getngF() const { return guard_cells.ng_alloc_F; }
    const amrex::IntVect getngUpdateAux() const { return guard_cells.ng_UpdateAux; }
    const amrex::IntVect get_ng_depos_J() const {return guard_cells.ng_depos_J;}
    const amrex::IntVect get_ng_depos_rho() const {return guard_cells.ng_depos_rho;}

    /** Coarsest-level Domain Decomposition
     *
     * If specified, the domain will be chopped into the exact number
     * of pieces in each dimension as specified by this parameter.
     *
     * @return the number of MPI processes per dimension if specified, otherwise a 0-vector
     */
    const amrex::IntVect get_numprocs() const {return numprocs;}

    void ComputeSpaceChargeField (bool const reset_fields);
    void AddSpaceChargeField (WarpXParticleContainer& pc);
    void AddSpaceChargeFieldLabFrame ();
    void computePhi (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                     amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                     std::array<amrex::Real, 3> const beta = {{0,0,0}},
                     amrex::Real const required_precision=amrex::Real(1.e-11),
                     const int max_iters=200) const;
    void computePhiRZ (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                       amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                       std::array<amrex::Real, 3> const beta,
                       amrex::Real const required_precision,
                       int const max_iters) const;
    void computePhiCartesian (const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                       amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                       std::array<amrex::Real, 3> const beta,
                       amrex::Real const required_precision,
                       int const max_iters) const;

    void computeE (amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& E,
                   const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<amrex::Real, 3> const beta = {{0,0,0}} ) const;
    void computeB (amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3> >& B,
                   const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi,
                   std::array<amrex::Real, 3> const beta = {{0,0,0}} ) const;

    /**
     * \brief
     * This function initializes the E and B fields on each level
     * using the parser and the user-defined function for the external fields.
     * The subroutine will parse the x_/y_z_external_grid_function and
     * then, the field multifab is initialized based on the (x,y,z) position
     * on the staggered yee-grid or cell-centered grid, in the interior cells
     * and guard cells.
     *
     * \param[in] mfx, x-component of the field to be initialized
     * \param[in] mfy, y-component of the field to be initialized
     * \param[in] mfz, z-component of the field to be initialized
     * \param[in] xfield_parser, parser function to initialize x-field
     * \param[in] yfield_parser, parser function to initialize y-field
     * \param[in] zfield_parser, parser function to initialize z-field
     * \param[in] lev, level of the Multifabs that is initialized
     */
    void InitializeExternalFieldsOnGridUsingParser (
         amrex::MultiFab *mfx, amrex::MultiFab *mfy, amrex::MultiFab *mfz,
         HostDeviceParser<3> const& xfield_parser, HostDeviceParser<3> const& yfield_parser,
         HostDeviceParser<3> const& zfield_parser, const int lev);

    /** \brief adds particle and cell contributions in cells to compute heuristic
     * cost in each box on each level, and records in `costs`
     * @param[in] costs vector of (`unique_ptr` to) vectors; expected to be initialized
     * to correct number of boxes and boxes per level
     */
    void ComputeCostsHeuristic (amrex::Vector<std::unique_ptr<amrex::LayoutData<amrex::Real> > >& costs);

    void ApplyFilterandSumBoundaryRho (int lev, int glev, amrex::MultiFab& rho, int icomp, int ncomp);

#ifdef WARPX_USE_PSATD
    // Device vectors of stencil coefficients used for finite-order centering of fields
    amrex::Gpu::DeviceVector<amrex::Real> device_field_centering_stencil_coeffs_x;
    amrex::Gpu::DeviceVector<amrex::Real> device_field_centering_stencil_coeffs_y;
    amrex::Gpu::DeviceVector<amrex::Real> device_field_centering_stencil_coeffs_z;
    // Device vectors of stencil coefficients used for finite-order centering of currents
    amrex::Gpu::DeviceVector<amrex::Real> device_current_centering_stencil_coeffs_x;
    amrex::Gpu::DeviceVector<amrex::Real> device_current_centering_stencil_coeffs_y;
    amrex::Gpu::DeviceVector<amrex::Real> device_current_centering_stencil_coeffs_z;
#endif

protected:

    /**
     * \brief
     *  This function initializes E, B, rho, and F, at all the levels
     *  of the multifab. rho and F are initialized with 0.
     *  The E and B fields are initialized using user-defined inputs.
     *  The initialization type is set using "B_ext_grid_init_style"
     *  and "E_ext_grid_init_style". The initialization style is set to "default"
     *  if not explicitly defined by the user, and the E and B fields are
     *  initialized with E_external_grid and B_external_grid, respectively, each with
     *  a default value of 0.
     *  If the initialization type for the E and B field is "constant",
     *  then, the E and B fields at all the levels are initialized with
     *  user-defined values for E_external_grid and B_external_grid.
     *  If the initialization type for B-field is set to
     *  "parse_B_ext_grid_function", then, the parser is used to read
     *  Bx_external_grid_function(x,y,z), By_external_grid_function(x,y,z),
     *  and Bz_external_grid_function(x,y,z).
     *  Similarly, if the E-field initialization type is set to
     *  "parse_E_ext_grid_function", then, the parser is used to read
     *  Ex_external_grid_function(x,y,z), Ey_external_grid_function(x,y,z),
     *  and Ex_external_grid_function(x,y,z). The parser for the E and B
     *  initialization assumes that the function has three independent
     *  variables, at max, namely, x, y, z. However, any number of constants
     *  can be used in the function used to define the E and B fields on the grid.
     */
    void InitLevelData (int lev, amrex::Real time);

    //! Tagging cells for refinement
    virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) final;

    //! Use this function to override the Level 0 grids made by AMReX.
    //! This function is called in amrex::AmrCore::InitFromScratch.
    virtual void PostProcessBaseGrids (amrex::BoxArray& ba0) const final;

    //! Make a new level from scratch using provided BoxArray and
    //! DistributionMapping.  Only used during initialization.  Called
    //! by AmrCoreInitFromScratch.
    virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
                                          const amrex::DistributionMapping& dm) final;

    //! Make a new level using provided BoxArray and
    //! DistributionMapping and fill with interpolated coarse level
    //! data.  Called by AmrCore::regrid.
    virtual void MakeNewLevelFromCoarse (int /*lev*/, amrex::Real /*time*/, const amrex::BoxArray& /*ba*/,
                                         const amrex::DistributionMapping& /*dm*/) final
        { amrex::Abort("MakeNewLevelFromCoarse: To be implemented"); }

    //! Remake an existing level using provided BoxArray and
    //! DistributionMapping and fill with existing fine and coarse
    //! data.  Called by AmrCore::regrid.
    virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
                              const amrex::DistributionMapping& dm) final;

    //! Delete level data.  Called by AmrCore::regrid.
    virtual void ClearLevel (int lev) final;

private:

    // Singleton is used when the code is run from python
    static WarpX* m_instance;

    ///
    /// Advance the simulation by numsteps steps, electromagnetic case.
    ///
    void EvolveEM(int numsteps);

    void FillBoundaryB (int lev, PatchType patch_type, amrex::IntVect ng);
    void FillBoundaryE (int lev, PatchType patch_type, amrex::IntVect ng);
    void FillBoundaryF (int lev, PatchType patch_type, amrex::IntVect ng);

    void FillBoundaryB_avg (int lev, PatchType patch_type, amrex::IntVect ng);
    void FillBoundaryE_avg (int lev, PatchType patch_type, amrex::IntVect ng);

    /**
     * \brief Synchronize the nodal points of the electric field MultiFabs
     */
    void NodalSyncE ();

    /**
     * \brief Synchronize the nodal points of the electric field MultiFabs for given MR level
     */
    void NodalSyncE (int lev);

    /**
     * \brief Synchronize the nodal points of the electric field MultiFabs for given MR level and patch
     */
    void NodalSyncE (int lev, PatchType patch_type);

    /**
     * \brief Synchronize the nodal points of the magnetic field MultiFabs
     */
    void NodalSyncB ();

    /**
     * \brief Synchronize the nodal points of the magnetic field MultiFabs for given MR level
     */
    void NodalSyncB (int lev);

    /**
     * \brief Synchronize the nodal points of the magnetic field MultiFabs for given MR level and patch
     */
    void NodalSyncB (int lev, PatchType patch_type);

    void OneStep_nosub (amrex::Real t);
    void OneStep_sub1 (amrex::Real t);

    void RestrictCurrentFromFineToCoarsePatch (int lev);
    void AddCurrentFromFineLevelandSumBoundary (int lev);
    void StoreCurrent (int lev);
    void RestoreCurrent (int lev);
    void ApplyFilterandSumBoundaryJ (int lev, PatchType patch_type);
    void NodalSyncJ (int lev, PatchType patch_type);

    void RestrictRhoFromFineToCoarsePatch (int lev);
    void ApplyFilterandSumBoundaryRho (int lev, PatchType patch_type, int icomp, int ncomp);
    void AddRhoFromFineLevelandSumBoundary (int lev, int icomp, int ncomp);
    void NodalSyncRho (int lev, PatchType patch_type, int icomp, int ncomp);

    /**
     * \brief Private function for current correction in Fourier space
     * (equation (19) of https://doi.org/10.1016/j.jcp.2013.03.010):
     * loops over the MR levels and applies the correction on the fine and coarse
     * patches (calls the virtual method \c CurrentCorrection of the spectral
     * algorithm in use, via the public interface defined in the class SpectralSolver).
     */
    void CurrentCorrection ();

    /**
     * \brief Private function for Vay deposition in Fourier space
     * (equations (20)-(24) of https://doi.org/10.1016/j.jcp.2013.03.010):
     * loops over the MR levels and applies the correction on the fine and coarse
     * patches (calls the virtual method \c VayDeposition of the spectral
     * algorithm in use, via the public interface defined in the class SpectralSolver).
     */
    void VayDeposition ();

    void ReadParameters ();

    /** This function queries deprecated input parameters and abort
     *  the run if one of them is specified. */
    void BackwardCompatibility ();

    void InitFromScratch ();

    void AllocLevelData (int lev, const amrex::BoxArray& new_grids,
                         const amrex::DistributionMapping& new_dmap);

    void InitFromCheckpoint ();
    void PostRestart ();

    void InitPML ();
    void ComputePMLFactors ();

    void InitFilter ();

    void InitDiagnostics ();

    void InitNCICorrector ();

    /**
     * \brief Check that the number of guard cells is smaller than the number of valid cells,
     * for all available MultiFabs, and abort otherwise.
     */
    void CheckGuardCells();

    /**
     * \brief Check that the number of guard cells is smaller than the number of valid cells,
     * for a given MultiFab, and abort otherwise.
     */
    void CheckGuardCells(amrex::MultiFab const& mf);

    /** Check the requested resources and write performance hints */
    void PerformanceHints ();

    std::unique_ptr<amrex::MultiFab> GetCellCenteredData();

    void BuildBufferMasks ();
    void BuildBufferMasksInBox ( const amrex::Box tbx, amrex::IArrayBox &buffer_mask,
                                 const amrex::IArrayBox &guard_mask, const int ng );
    const amrex::iMultiFab* getCurrentBufferMasks (int lev) const {
        return current_buffer_masks[lev].get();
    }
    const amrex::iMultiFab* getGatherBufferMasks (int lev) const {
        return gather_buffer_masks[lev].get();
    }

#ifdef WARPX_USE_PSATD
    /**
     * \brief Re-orders the Fornberg coefficients so that they can be used more conveniently for
     * finite-order centering operations. For example, for finite-order centering of order 6,
     * the Fornberg coefficients \c (c_0,c_1,c_2) are re-ordered as \c (c_2,c_1,c_0,c_0,c_1,c_2).
     *
     * \param[in,out] ordered_coeffs host vector where the re-ordered Fornberg coefficients will be stored
     * \param[in] unordered_coeffs host vector storing the original sequence of Fornberg coefficients
     * \param[in] order order of the finite-order centering along a given direction
     */
    void ReorderFornbergCoefficients (amrex::Vector<amrex::Real>& ordered_coeffs,
                                      amrex::Vector<amrex::Real>& unordered_coeffs,
                                      const int order);
    /**
     * \brief Allocates and initializes the stencil coefficients used for the finite-order centering
     * of fields and currents, and stores them in the given device vectors.
     *
     * \param[in,out] device_centering_stencil_coeffs_x device vector where the stencil coefficients along x will be stored
     * \param[in,out] device_centering_stencil_coeffs_y device vector where the stencil coefficients along y will be stored
     * \param[in,out] device_centering_stencil_coeffs_z device vector where the stencil coefficients along z will be stored
     * \param[in] centering_nox order of the finite-order centering along x
     * \param[in] centering_noy order of the finite-order centering along y
     * \param[in] centering_noz order of the finite-order centering along z
     */
    void AllocateCenteringCoefficients (amrex::Gpu::DeviceVector<amrex::Real>& device_centering_stencil_coeffs_x,
                                        amrex::Gpu::DeviceVector<amrex::Real>& device_centering_stencil_coeffs_y,
                                        amrex::Gpu::DeviceVector<amrex::Real>& device_centering_stencil_coeffs_z,
                                        const int centering_nox,
                                        const int centering_noy,
                                        const int centering_noz);
#endif

    void AllocLevelMFs (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm,
                        const amrex::IntVect& ngE, const amrex::IntVect& ngJ,
                        const amrex::IntVect& ngRho, const amrex::IntVect& ngF,
                        const bool aux_is_nodal);
#ifdef WARPX_USE_PSATD
#   ifdef WARPX_DIM_RZ
    void AllocLevelSpectralSolverRZ (amrex::Vector<std::unique_ptr<SpectralSolverRZ>>& spectral_solver,
                                     const int lev,
                                     const amrex::BoxArray& realspace_ba,
                                     const amrex::DistributionMapping& dm,
                                     const std::array<amrex::Real,3>& dx);
#   else
    void AllocLevelSpectralSolver (amrex::Vector<std::unique_ptr<SpectralSolver>>& spectral_solver,
                                   const int lev,
                                   const amrex::BoxArray& realspace_ba,
                                   const amrex::DistributionMapping& dm,
                                   const std::array<amrex::Real,3>& dx,
                                   const bool pml_flag=false);
#   endif
#endif

    amrex::Vector<int> istep;      // which step?
    amrex::Vector<int> nsubsteps;  // how many substeps on each level?

    amrex::Vector<amrex::Real> t_new;
    amrex::Vector<amrex::Real> t_old;
    amrex::Vector<amrex::Real> dt;

    // Particle container
    std::unique_ptr<MultiParticleContainer> mypc;
    std::unique_ptr<MultiDiagnostics> multi_diags;

    // Boosted Frame Diagnostics
    std::unique_ptr<BackTransformedDiagnostic> myBFD;

    //
    // Fields: First array for level, second for direction
    //

    // Full solution
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_aux;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_aux;

    // Fine patch
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > F_fp;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_fp;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > phi_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_avg_fp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_avg_fp;
    // store fine patch
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_store;

    // Nodal MultiFab for nodal current deposition if warpx.do_current_centering = 1
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>> current_fp_nodal;

    // Coarse patch
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > F_cp;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_avg_cp;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_avg_cp;

    // Copy of the coarse aux
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3 > > Efield_cax;
    amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_cax;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > current_buffer_masks;
    amrex::Vector<std::unique_ptr<amrex::iMultiFab> > gather_buffer_masks;

    // If charge/current deposition buffers are used
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_buf;
    amrex::Vector<std::unique_ptr<amrex::MultiFab> > charge_buf;

    // PML
    int do_pml = 1;
    int do_silver_mueller = 0;
    int pml_ncell = 10;
    int pml_delta = 10;
    int pml_has_particles = 0;
    int do_pml_j_damping = 0;
    int do_pml_in_domain = 0;
    bool do_pml_dive_cleaning; // default set in WarpX.cpp
    bool do_pml_divb_cleaning; // default set in WarpX.cpp
    amrex::IntVect do_pml_Lo = amrex::IntVect::TheUnitVector();
    amrex::IntVect do_pml_Hi = amrex::IntVect::TheUnitVector();
    amrex::Vector<std::unique_ptr<PML> > pml;

    amrex::Real moving_window_x = std::numeric_limits<amrex::Real>::max();
    amrex::Real current_injection_position = 0;

    // Plasma injection parameters
    int warpx_do_continuous_injection = 0;
    int num_injected_species = -1;
    amrex::Vector<int> injected_plasma_species;

    int n_buffer = 4;
    amrex::Real const_dt = amrex::Real(0.5e-11);

    // Macroscopic properties
    std::unique_ptr<MacroscopicProperties> m_macroscopic_properties;


    // Load balancing
    /** Load balancing intervals that reads the "load_balance_intervals" string int the input file
     * for getting steps at which load balancing is performed */
    IntervalsParser load_balance_intervals;
    /** Collection of LayoutData to keep track of weights used in load balancing
     * routines. Contains timer-based or heuristic-based costs depending on input option */
    amrex::Vector<std::unique_ptr<amrex::LayoutData<amrex::Real> > > costs;
    /** Load balance with 'space filling curve' strategy. */
    int load_balance_with_sfc = 0;
    /** Controls the maximum number of boxes that can be assigned to a rank during
     * load balance via the 'knapsack' strategy; e.g., if there are 4 boxes per rank,
     * `load_balance_knapsack_factor=2` limits the maximum number of boxes that can
     * be assigned to a rank to 8. */
    amrex::Real load_balance_knapsack_factor = amrex::Real(1.24);
    /** Threshold value that controls whether to adopt the proposed distribution
     * mapping during load balancing.  The new distribution mapping is adopted
     * if the ratio of proposed distribution mapping efficiency to current
     * distribution mapping efficiency is larger than the threshold; 'efficiency'
     * here means the average cost per MPI rank.  */
    amrex::Real load_balance_efficiency_ratio_threshold = amrex::Real(1.1);
    /** Current load balance efficiency for each level.  */
    amrex::Vector<amrex::Real> load_balance_efficiency;
    /** Weight factor for cells in `Heuristic` costs update.
     * Default values on GPU are determined from single-GPU tests on Summit.
     * The problem setup for these tests is an empty (i.e. no particles) domain
     * of size 256 by 256 by 256 cells, from which the average time per iteration
     * per cell is computed. */
    amrex::Real costs_heuristic_cells_wt = amrex::Real(-1);
    /** Weight factor for particles in `Heuristic` costs update.
     * Default values on GPU are determined from single-GPU tests on Summit.
     * The problem setup for these tests is a high-ppc (27 particles per cell)
     * uniform plasma on a domain of size 128 by 128 by 128, from which the approximate
     * time per iteration per particle is computed. */
    amrex::Real costs_heuristic_particles_wt = amrex::Real(-1);

    // Determines timesteps for override sync
    IntervalsParser override_sync_intervals;

    // Other runtime parameters
    int verbose = 1;

    bool use_hybrid_QED = 0;

    int max_step   = std::numeric_limits<int>::max();
    amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

    int regrid_int = -1;

    amrex::Real cfl = amrex::Real(0.7);

    std::string restart_chkfile;

    bool plot_rho = false;

    amrex::VisMF::Header::Version plotfile_headerversion  = amrex::VisMF::Header::Version_v1;
    amrex::VisMF::Header::Version slice_plotfile_headerversion  = amrex::VisMF::Header::Version_v1;

    bool use_single_read = true;
    bool use_single_write = true;
    int mffile_nstreams = 4;
    int field_io_nfiles = 1024;
    int particle_io_nfiles = 1024;

    amrex::RealVect fine_tag_lo;
    amrex::RealVect fine_tag_hi;

    bool is_synchronized = true;

    guardCellManager guard_cells;

    //Slice Parameters
    int slice_max_grid_size;
    int slice_plot_int = -1;
    amrex::RealBox slice_realbox;
    amrex::IntVect slice_cr_ratio;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > F_slice;
    amrex::Vector<            std::unique_ptr<amrex::MultiFab>      > rho_slice;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > current_slice;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Efield_slice;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > Bfield_slice;

    bool fft_periodic_single_box = false;
    int nox_fft = 16;
    int noy_fft = 16;
    int noz_fft = 16;

    // Domain decomposition on Level 0
    amrex::IntVect numprocs{0};

    //
    // Embeded Boundary
    //

    // Factory for field data
    amrex::Vector<std::unique_ptr<amrex::FabFactory<amrex::FArrayBox> > > m_field_factory;

    amrex::FabFactory<amrex::FArrayBox> const& fieldFactory (int lev) const noexcept {
        return *m_field_factory[lev];
    }
#ifdef AMREX_USE_EB
    amrex::EBFArrayBoxFactory const& fieldEBFactory (int lev) const noexcept {
        return static_cast<amrex::EBFArrayBoxFactory const&>(*m_field_factory[lev]);
    }
#endif

    void InitEB ();

private:
    // void EvolvePSATD (int numsteps);
    void PushPSATD (amrex::Real dt);
    void PushPSATD (int lev, amrex::Real dt);

    int fftw_plan_measure = 1; // used with PSATD

#ifdef WARPX_USE_PSATD
#   ifdef WARPX_DIM_RZ
        amrex::Vector<std::unique_ptr<SpectralSolverRZ>> spectral_solver_fp;
        amrex::Vector<std::unique_ptr<SpectralSolverRZ>> spectral_solver_cp;
#   else
        amrex::Vector<std::unique_ptr<SpectralSolver>> spectral_solver_fp;
        amrex::Vector<std::unique_ptr<SpectralSolver>> spectral_solver_cp;
#   endif

public:

#   ifdef WARPX_DIM_RZ
        SpectralSolverRZ&
#   else
        SpectralSolver&
#   endif
            get_spectral_solver_fp (int lev) {return *spectral_solver_fp[lev];}
#endif

private:
    amrex::Vector<std::unique_ptr<FiniteDifferenceSolver>> m_fdtd_solver_fp;
    amrex::Vector<std::unique_ptr<FiniteDifferenceSolver>> m_fdtd_solver_cp;
};

#endif
