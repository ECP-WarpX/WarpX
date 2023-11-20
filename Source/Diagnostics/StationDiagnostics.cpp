#include "StationDiagnostics.H"

#include "ComputeDiagFunctors/StationFunctor.H"
//temporary
#include "ComputeDiagFunctors/CellCenterFunctor.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "WarpX.H"

//temporary
#include <ablastr/coarsen/sample.H>
#include <ablastr/utils/Communication.H>
#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_CoordSys.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_FabFactory.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_MakeType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_VisMF.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using namespace amrex::literals;

StationDiagnostics::StationDiagnostics (int i, std::string name)
    : Diagnostics(i, std::move(name))
{
    ReadParameters();
}

void
StationDiagnostics::ReadParameters ()
{
#ifdef WARPX_DIM_RZ
    amrex::Abort("StationDiagnostics is not implemented for RZ, yet");
#endif

    const amrex::ParmParse pp_diag_name(m_diag_name);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_format == "plotfile", "<diag>.format must be plotfile");

    std::string station_normal_dir;
    pp_diag_name.get("station_normal_dir", station_normal_dir);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        station_normal_dir == "z", "<diag>.station_normal_dir is supported only for z at the moment");
    if (station_normal_dir == "z") {
        m_station_normal = StationNormalDir::z;
    }

    pp_diag_name.get("station_location",m_station_loc);

    // If station recordings are used to restart simulations, then raw fields are needed
    pp_diag_name.query("plot_raw_fields", m_plot_raw_fields);
    pp_diag_name.query("plot_raw_fields", m_plot_raw_fields_guards);

    pp_diag_name.query("buffer_size",m_buffer_size);
    // for now, number of buffers or in this case, z-locations is assumed to be 1
    // This is used to allocate the number of output multi-level multifabs, m_mf_output
    m_num_buffers = 1;

    std::string intervals_string_vec = "0";
    pp_diag_name.query("intervals", intervals_string_vec);
    m_intervals = utils::parser::SliceParser(intervals_string_vec);

    m_varnames = {"Ex", "Ey", "Ez", "Bx", "By", "Bz"};
    m_file_prefix = "diags/" + m_diag_name;
    pp_diag_name.query("file_prefix", m_file_prefix);
}

bool
StationDiagnostics::DoDump (int step, int /* i_buffer*/, bool force_flush)
{
    // Determine criterion to dump output

    return ( (m_slice_counter == m_buffer_size) || force_flush );
}

bool
StationDiagnostics::DoComputeAndPack (int step, bool force_flush)
{
    // we may have to compute and pack everytimestep, but only store at user-defined intervals
    return ( step>=0 );
}

void
StationDiagnostics::InitializeFieldFunctors (int lev)
{
    // need to define a functor that will take a slice of the fields and store it in a given location
    // in the multifab based on current time
    // In this function, we will define m_all_field_functors based on that functor
    auto & warpx = WarpX::GetInstance();
    const int num_stationdiag_functors = 1;
    m_all_field_functors[lev].resize(num_stationdiag_functors);
    int ncomp = 6;
    for (int i = 0; i < num_stationdiag_functors; ++i)
    {
        m_all_field_functors[lev][i] = std::make_unique<StationFunctor>(
                                           warpx.get_array_EBfield_fp(lev), m_station_loc, lev,
                                           m_crse_ratio, ncomp
                                       );
    }
}

void
StationDiagnostics::InitializeBufferData (int i_buffer, int lev, bool restart)
{
    // Define boxArray, dmap, and initialize output multifab
    // This will have the extension in x-y at a given location z, and the third dimension will be time

    auto & warpx = WarpX::GetInstance();

    amrex::IntVect lo(0);
    amrex::IntVect hi(m_buffer_size-1);
    // in the station normal direction, with time, we need to determine number of points
    // hi : (t_max - t_min )/dt
    const amrex::Box diag_box(lo, hi);
    m_buffer_box = diag_box;
    amrex::BoxArray diag_ba;
    diag_ba.define(diag_box);
    amrex::BoxArray ba = diag_ba.maxSize(256);
    amrex::DistributionMapping dmap(ba);
    int ncomps = 6; //  Ex Ey Ez Bx By Bz
    int nghost = 1; //  Ex Ey Ez Bx By Bz
    m_mf_output[0][lev] = amrex::MultiFab( amrex::convert(ba, amrex::IntVect::TheNodeVector()), dmap, ncomps, nghost);
    m_mf_output[0][lev].setVal(0.);
}

void
StationDiagnostics::PrepareFieldDataForOutput ()
{
    const int num_station_functors = 1;
    const int nlev = 1;
    int k_index = m_slice_counter;
    for (int lev = 0; lev < nlev; ++lev) {
        for (int i = 0; i < num_station_functors; ++i) {
            // number of slices = 1
            m_all_field_functors[lev][i]->PrepareFunctorData(0, true, m_station_loc, m_buffer_box,
                                                             m_slice_counter, m_buffer_size, 0);
        }
    }
}


void
StationDiagnostics::UpdateBufferData ()
{
    m_slice_counter++;
}

void
StationDiagnostics::InitializeParticleBuffer ()
{

}

void
StationDiagnostics::Flush (int i_buffer)
{
    auto & warpx = WarpX::GetInstance();
    m_tmax = warpx.gett_new(0);
    std::string filename = m_file_prefix;
    filename = amrex::Concatenate(m_file_prefix, i_buffer, m_file_min_digits);
    constexpr int permission_flag_rwxrxrx = 0755;
    if (! amrex::UtilCreateDirectory(filename, permission_flag_rwxrxrx) ) {
        amrex::CreateDirectoryFailed(filename);
    }
    WriteStationHeader(filename);
    for (int lev = 0; lev < 1; ++lev) {
        const std::string buffer_path = filename + amrex::Concatenate("/Level_",lev,1) + "/";
        if (m_flush_counter == 0) {
            if (! amrex::UtilCreateDirectory(buffer_path, permission_flag_rwxrxrx) ) {
                amrex::CreateDirectoryFailed(buffer_path);
            }
        }
        std::string buffer_string = amrex::Concatenate("buffer-",m_flush_counter,m_file_min_digits);
        const std::string prefix = amrex::MultiFabFileFullPrefix(lev, filename, "Level_", buffer_string);
        amrex::VisMF::Write(m_mf_output[i_buffer][lev], prefix);
    }

    // reset counter
    m_slice_counter = 0;
    // update Flush counter
    m_flush_counter++;
}

void
StationDiagnostics::WriteStationHeader (const std::string& filename)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        const std::string HeaderFileName(filename + "/StationHeader");
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good())
            amrex::FileOpenFailed(HeaderFileName);

        HeaderFile.precision(17);

        HeaderFile << m_tmin << "\n";
        HeaderFile << m_tmax << "\n";
        HeaderFile << m_buffer_size << "\n";
        HeaderFile << m_slice_counter << "\n";
        HeaderFile << m_flush_counter << "\n";


    }
}
