#ifndef WARPX_FULLDIAGNOSTICS_H_
#define WARPX_FULLDIAGNOSTICS_H_

#include "Diagnostics.H"

class
FullDiagnostics final : public Diagnostics
{
public:
    FullDiagnostics (int i, std::string name);
//    ~FullDiagnostics ();
//    void BTInitData (){};
private:
    void ReadParameters ();
    /** Determines timesteps at which full diagnostics are written to file */
    IntervalsParser m_intervals;
    /** Whether to plot raw (i.e., NOT cell-centered) fields */
    bool m_plot_raw_fields = false;
    /** Whether to plot guard cells of raw fields */
    bool m_plot_raw_fields_guards = false;
    /** Whether to plot charge density rho in raw fields */
    bool m_plot_raw_rho = false;
    /** Whether to plot F (charge conservation error) in raw fields */
    bool m_plot_raw_F = false;
    /** Flush m_mf_output and particles to file. */
    void Flush () override; 
    /** Flush raw data */
    void FlushRaw () override;
    bool DoDump (int step, bool force_flush=false) override;
    void InitData () override;
    /** Append varnames with names for all modes of a field
     * \param[in] field field name (includes component, e.g., Er)
     * \param[in] ncomp number of components (modes, real and imag)
     */
    void AddRZModesToOutputNames (const std::string& field, int ncomp);
    /** Append m_all_field_functors[lev] with all modes of all components of E B and j
     * and add name to varnames.
     * \param[in] lev level on which source multifabs are defined.
     */
    void AddRZModesToDiags (int lev);
    /** Whether to dump the RZ modes */
    bool m_dump_rz_modes = false;
    /** Define the cell-centered multifab m_mf_output depending on user-defined
      * lo and hi and coarsening ratio. This MultiFab may have a different BoxArray and
      * DistributionMap than field MultiFabs in the simulation.
      *
      * \param[in] lev level on which source multifabs are defined
      */
    void DefineDiagMultiFab ( int lev );
    /** Initialize functors that store pointers to the fields requested by the user.
      * \param[in] lev level on which the vector of unique_ptrs to field functors is initialized.
      */
    void InitializeFieldFunctors (int lev) override;
    /** Pack (stack) all fields in the cell-centered output MultiFab m_mf_output.
     * Fields are computed (e.g., cell-centered) on-the-fly using a functor.*/
    void ComputeAndPack () override;
    
};

#endif // WARPX_FULLDIAGNOSTICS_H_
