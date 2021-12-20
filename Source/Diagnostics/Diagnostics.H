#ifndef WARPX_DIAGNOSTICS_H_
#define WARPX_DIAGNOSTICS_H_

#include "ParticleDiag/ParticleDiag.H"

#include "ComputeDiagFunctors/ComputeDiagFunctor_fwd.H"
#include "FlushFormats/FlushFormat_fwd.H"

#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#include <memory>
#include <string>
#include <vector>

/**
 * \brief base class for diagnostics.
 * Contains main routines to filter, compute and flush diagnostics.
 *
 * Each specific diagnostics derives from this class.
 */
class Diagnostics
{
public:
    Diagnostics (int i, std::string name);
    /** Virtual Destructor to handle clean destruction of derived classes */
    virtual ~Diagnostics () ;
    /** Pack (stack) all fields in the cell-centered output MultiFab m_mf_output.
     *
     * Fields are computed (e.g., cell-centered or back-transformed)
       on-the-fly using a functor. */
    void ComputeAndPack ();
    /** \brief Flush particle and field buffers to file using the FlushFormat member variable.
     *
     * This function should belong to class Diagnostics and not be virtual, as it flushes
     * the particle buffers (name not found yet) and the field buffers (m_mf_output), both
     * of which are members of Diagnostics. Yet, the implementation is left to derived classes
     * for now because:
     * - The functions underlying FlushFormat::WriteToFile expect a geometry object, which is
     *   WarpX::geom for FullDiagnostics but should be constructed for BTDiagnostics;
     * - The functions underlying FlushFormat::WriteToFile do not support writing a buffer to file
     *   multiple times yet.
     *  When these are fixed, the implementation of Flush should be in Diagnostics.cpp
     * \param[in] i_buffer index of the buffer data to be flushed.
     */
    virtual void Flush (int i_buffer) = 0;
    /** Initialize pointers to main fields and allocate output multifab m_mf_output. */
    void InitData ();
    /** Initialize functors that store pointers to the fields requested by the user.
     *
     * Derived classes MUST implement this function, and it must allocate m_all_field_functors
     * and fill it with ComputeDiagFunctor objects.
     * The functions is called at intiailization and when the domain is decomposed
     * during the simulation to load-balance.
     * \param[in] lev level on which the vector of unique_ptrs to field functors is initialized.
     */
    virtual void InitializeFieldFunctors (int lev) = 0;
    /** whether to compute and pack data in output buffers at this time step
     * \param[in] step current time step
     * \param[in] force_flush if true, return true for any step
     * \return bool, whether to compute and pack data
     */
    virtual bool DoComputeAndPack(int step, bool force_flush=false) = 0;
    /** whether to flush at this time step
     * \param[in] step current time step
     * \param[in] i_buffer index of the buffer data to be flushed.
     * \param[in] force_flush if true, return true for any step
     * \return bool, whether to flush
     */
    virtual bool DoDump (int step, int i_buffer, bool force_flush=false) = 0;
    /** Start a new iteration, i.e., dump has not been done yet. */
    void NewIteration () {m_already_done = false;}
    /** Perform necessary operations with user-defined diagnostic parameters
     *  to filter (coarsen, slice), compute (cell-center, back-transform),
     *  and flush the output data stored in buffers, m_mf_output.
     * \param[in] step current timestep
     * \param[in] force_flush used to force-fully write data stored in buffers.
     */
    void FilterComputePackFlush (int step, bool force_flush=false);
    /** Whether the last timestep is always dumped */
    bool DoDumpLastTimestep () const {return  m_dump_last_timestep;}

protected:
    /** Read Parameters of the base Diagnostics class */
    bool BaseReadParameters ();
    /** Initialize member variables of the base Diagnostics class. */
    void InitBaseData ();
    /** Initialize m_mf_output vectors and data required to construct the buffers
     * \param[in] i_buffer index of buffer for which the output MultiFab is defined.
     * \param[in] lev level on which the output MultiFab is defined
     */
    virtual void InitializeFieldBufferData (int i_buffer, int lev) = 0;
    /** Initialize member variables and arrays specific to the diagnostics in the
     *  derived classes.(FullDiagnostics, BTDiagnostics)
     */
    virtual void DerivedInitData () {}
    /** This function initialized particle buffers (not implemented in diagnostics, yet) */
    virtual void InitializeParticleBuffer () = 0;
    /** Prepare data (either fill-boundary or cell-centered data for
        back-transform diagnostics) to be processed for diagnostics.
     */
    virtual void PrepareFieldDataForOutput () {}
    /** Update the physical extent of the diagnostic domain for moving window and
     *  galilean shift simulations
     *
     * \param[in] step current time step
     */
    virtual void MovingWindowAndGalileanDomainShift (int step) { amrex::ignore_unused(step); }
    /** Name of diagnostics: runtime parameter given in the input file. */
    std::string m_diag_name;
    /** Prefix for output directories */
    std::string m_file_prefix;
    /** Minimum number of digits to iteration number in file name */
    int m_file_min_digits = 5;
    /** Index of diagnostics in MultiDiagnostics::alldiags */
    int m_diag_index;
    /** Names of  each component requested by the user.
     * in cylindrical geometry, this list is appended with
     * automatically-constructed names for all modes of all fields */
    amrex::Vector< std::string > m_varnames;
    /** If true, a dump is performed at the last timestep regardless of the required dump timesteps */
    bool m_dump_last_timestep = true;
    /** format for output files, "plotfile" or "openpmd" or "sensei" or "ascent"
     * The checkpoint format is applicable for FullDiagnostics only.
     */
    std::string m_format = "plotfile";
    /** Whether this iteration has already been dumped, to avoid writing data twice */
    int m_already_done = false;
    /** This class is responsible for flushing the data to file */
    std::unique_ptr<FlushFormat> m_flush_format;
    /** output multifab, where all fields are computed (cell-centered or back-transformed)
     *  and stacked.
     *  The first vector is for total number of snapshots. (=1 for FullDiagnostics)
     *  The second vector is loops over the total number of levels.
     */
    amrex::Vector< amrex::Vector< amrex::MultiFab > >  m_mf_output;

    /** Geometry that defines the domain attributes corresponding to output multifab.
     *  Specifically, the user-defined physical co-ordinates for the diagnostics
     *  is used to construct the geometry information for each MultiFab at
     *  the respective levels. This geometry will be used to write out
     *  plotfile data using the WriteToFile() function
     */
    amrex::Vector< amrex::Vector <amrex::Geometry> > m_geom_output;
    // a particle buffer here?
    int nlev; /**< number of levels to output */
    int nmax_lev; /**< max_level to allocate output multifab and vector of field functors. */
    /** Number of levels to be output*/
    int nlev_output;
    /** Names of species to write to output */
    std::vector< std::string > m_output_species_names;
    /** Names of all species in the simulation */
    std::vector< std::string > m_all_species_names;
    /** Each element of this vector handles output for 1 species */
    amrex::Vector< ParticleDiag > m_output_species;
    /** Vector of (pointers to) functors to compute output fields, per level,
     * per component. This allows for simple operations (averaging to
     * cell-center for standard EB fields) as well as more involved operations
     * (back-transformed diagnostics, filtering, reconstructing cartesian
     * fields in cylindrical). */
    amrex::Vector< amrex::Vector <std::unique_ptr<ComputeDiagFunctor > > > m_all_field_functors;
    /** Coarsening ratio such that, fields are averaged to the coarsened grid.
      * The ratio should render the grid to be coarsenable (as defined by AMReX). */
    amrex::IntVect m_crse_ratio = amrex::IntVect(1);
    /** Lower corner of the diagnostics output, in physical coordinates */
    amrex::Vector< amrex::Real> m_lo;
    /** Higher corner of the diagnostics output, in physical coordinates */
    amrex::Vector< amrex::Real> m_hi;
    /** Number of output buffers. The value is set to 1 for all FullDiagnostics */
    int m_num_buffers;
    /** Array of species indices that dump rho per species */
    amrex::Vector<int> m_rho_per_species_index;
    /** Temporarily clear species output for BTD until particle buffer is added */
    virtual void TMP_ClearSpeciesDataForBTD() {}
};

#endif // WARPX_DIAGNOSTICS_H_
