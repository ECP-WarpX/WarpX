#ifndef WARPX_RECORDINGPLANEDIAGNOSTICS_H
#define WARPX_RECORDINGPLANEDIAGNOSTICS_H

#include "Diagnostics.H"
#include "Diagnostics/ComputeDiagFunctors/ComputeDiagFunctor.H"
#include "Utils/Parser/IntervalsParser.H"

#include <string>


class RecordingPlaneDiagnostics final
    : public Diagnostics
{
public:
    RecordingPlaneDiagnostics (int i, std::string name);
private:
    /** Read user-requested parameters for station diagnostics */
    utils::parser::SliceParser m_intervals;
    void ReadParameters ();

    /** Station location where field values are recorded in time */
    amrex::Real m_station_loc;

    /** The direction normal to the station slice where fields are recorded in time*/
#if defined(WARPX_DIM_1D_Z)
    enum struct StationNormalDir {z=0};
#elif defined(WARPX_DIM_XZ)
    enum struct StationNormalDir {x=0, z};
#elif defined(WARPX_DIM_RZ)
    enum struct StationNormalDir {r=0, z};
#else
    enum struct StationNormalDir {x=0, y, z};
#endif
    StationNormalDir m_station_normal = StationNormalDir::z;
    /** Whether to plot raw fields*/
    bool m_plot_raw_fields = true;
    /** Whether to plot guard cells of raw fields*/
    bool m_plot_raw_fields_guards = true;

    /** whether to flush at this time step
     * \param[in] step current time step
     * \param[in] i_buffer index of a back-transformed snapshot
     * \param[in] force_flush if true, return true for any step
     * \return bool, whether to flush
     */
    bool DoDump (int step, int i_buffer, bool force_flush=false) override;

    /** whether to compute and pack cell-centered data in m_mf_output
     * \param[in] step current time step
     * \param[in] force_flush if true, return true for any step since output must be
                  computed in this case
     * \return bool, whether to flush
     */
    bool DoComputeAndPack (int step, bool force_flush=false) override;

    void InitializeFieldFunctors (int lev) override;

    void InitializeBufferData (int i_buffer, int lev, bool restart=false) override;

    void InitializeParticleBuffer () override;

    void InitializeParticleFunctors () override;

    void Flush (int i_buffer, bool /*force_flush */) override;

    void WriteRecordingPlaneHeader (const std::string& filename);

    void PrepareFieldDataForOutput () override;

    void PrepareParticleDataForOutput () override;

    void UpdateBufferData () override;

    void FlushParticleBuffer (std::string path, std::string species_name, int isp, int i_buffer=0);
    void WriteHeaderFile (std::string pdir, amrex::Vector<std::string> real_names,
                          amrex::Vector<std::string> int_names, int isp);
    void WriteParticleData (std::string pdir, amrex::Vector<std::string> real_names,
                            amrex::Vector<std::string> int_names,
                            amrex::Vector<int> real_flags, amrex::Vector<int> int_flags, int isp);

    void ClearParticlesBuffer (int i_buffer);

    int m_buffer_size = 256;
    int m_slice_counter = 0;
    amrex::Real m_tmin = std::numeric_limits<amrex::Real>::lowest();
    amrex::Real m_tmax = std::numeric_limits<amrex::Real>::lowest();
    amrex::Box m_buffer_box;
    int m_flush_counter = 0;

    bool GetZSliceInDomain (const int lev);
    bool m_last_timeslice_filled = false;
};


#endif // WARPX_RECORDINGPLANEDIAGNOSTICS_H
