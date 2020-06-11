#ifndef WARPX_DIAGNOSTICS_H_
#define WARPX_DIAGNOSTICS_H_

#include "FlushFormats/FlushFormat.H"
#include "ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "ParticleDiag/ParticleDiag.H"
#include "Utils/IntervalsParser.H"
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

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
     * Fields are computed (e.g., cell-centered or back-transformed)
       on-the-fly using a functor.*/
    void ComputeAndPack ();
    /** Flush m_mf_output and particles to file. */
    virtual void Flush (){};
    /** Initialize pointers to main fields and allocate output multifab m_mf_output. */
    virtual void InitData (){};
    /** whether to flush at this time step
     * \param[in] step current time step
     * \param[in] force_flush if true, return true for any step
     * \return bool, whether to flush
     */
    virtual bool DoDump (int step, bool force_flush=false){ return false;};
    /** Initialize functors that store pointers to the fields requested by the user.
     * \param[in] lev level on which the vector of unique_ptrs to field functors is initialized.
     */
    virtual void InitializeFieldFunctors (int lev){};
    /** Start a new iteration, i.e., dump has not been done yet. */
    void NewIteration () {m_already_done = false;};

protected:
    /** Read Parameters of the base Diagnostics class */
    bool ReadBaseParameters ();
    /** Initialize member variables of the base Diagnostics class. */
    void InitBaseData ();
    /** Prepare data (either fill-boundary or cell-centered data for
        back-transform diagnostics) to be processed for diagnostics.
     */
    virtual void PrepareFieldDataForOutput () {};
    /** Name of diagnostics: runtime parameter given in the input file. */
    std::string m_diag_name;
    /** Prefix for output directories */
    std::string m_file_prefix;
    /** Index of diagnostics in MultiDiagnostics::alldiags */
    int m_diag_index;
    /** Names of  each component requested by the user.
     * in cylindrical geometry, this list is appended with
     * automatically-constructed names for all modes of all fields */
    amrex::Vector< std::string > m_varnames;
    /** format for output files, "plotfile" or "openpmd" */
    std::string m_format = "plotfile";
    /** Whether this iteration has already been dumped, to avoid writing data twice */
    int m_already_done = false;
    /** This class is responsible for flushing the data to file */
    FlushFormat* m_flush_format;
    /** output multifab, where all fields are computed (cell-centered or back-transformed)
     *  and stacked.
     *  The first vector is for total number of snapshots. (=1 for FullDiagnostics)
     *  The second vector is loops over the total number of levels.
     */
    amrex::Vector< amrex::Vector< amrex::MultiFab > >  m_mf_output;
    // a particle buffer here?
    int nlev; /**< number of levels to output */
    int nmax_lev; /**< max_level to allocate output multifab and vector of field functors. */
    /** Name of species to write to file */
    std::vector< std::string > m_species_names;
    /** Each element of this vector handles output for 1 species */
    amrex::Vector< ParticleDiag > m_all_species;
    /** Vector of (pointers to) functors to compute output fields, per level,
     * per component. This allows for simple operations (averaging to
     * cell-center for standard EB fields) as well as more involved operations
     * (back-transformed diagnostics, filtering, reconstructing cartesian
     * fields in cylindrical). */
    amrex::Vector< amrex::Vector <std::unique_ptr<ComputeDiagFunctor const> > > m_all_field_functors;
    /** Coarsening ratio such that, fields are averaged to the coarsened grid.
      * The ratio should render the grid to be coarsenable (as defined by AMReX). */
    amrex::IntVect m_crse_ratio = amrex::IntVect(1);
    /** Lower corner of the diagnostics output, in physical coordinates */
    amrex::Vector< amrex::Real> m_lo;
    /** Higher corner of the diagnostics output, in physical coordinates */
    amrex::Vector< amrex::Real> m_hi;
    /** Number of output buffers. The value is set to 1 for all FullDiagnostics */
    int m_num_buffers;
};

#endif // WARPX_DIAGNOSTICS_H_
