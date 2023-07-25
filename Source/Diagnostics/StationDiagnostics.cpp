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
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
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
    amrex::Print() << " in station diags \n";
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
        //m_all_field_functors[lev][i] = std::make_unique<StationFunctor>(
        //                                   warpx.get_array_EBfield_fp(lev), m_station_loc, lev,
        //                                   m_crse_ratio, ncomp
        //                               );
        // temporary
        m_all_field_functors[lev][i] = std::make_unique<StationFunctor>(
                                           m_cell_centered_data[lev].get(), m_station_loc, lev,
                                           m_crse_ratio, ncomp
                                       );
    }
    m_cell_center_functors[lev].resize(ncomp);
    m_cell_center_functors[lev][0] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev,0),
                                                                         lev, m_crse_ratio);
    m_cell_center_functors[lev][1] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev,1),
                                                                         lev, m_crse_ratio);
    m_cell_center_functors[lev][2] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev,2),
                                                                         lev, m_crse_ratio);
    m_cell_center_functors[lev][3] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev,0),
                                                                         lev, m_crse_ratio);
    m_cell_center_functors[lev][4] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev,1),
                                                                         lev, m_crse_ratio);
    m_cell_center_functors[lev][5] = std::make_unique<CellCenterFunctor>(warpx.get_pointer_Efield_aux(lev,2),
                                                                         lev, m_crse_ratio);
}

void
StationDiagnostics::InitializeBufferData (int i_buffer, int lev, bool restart)
{
    // Define boxArray, dmap, and initialize output multifab
    // This will have the extension in x-y at a given location z, and the third dimension will be time

    auto & warpx = WarpX::GetInstance();

    amrex::IntVect lo(0);
    amrex::IntVect hi(255);
    // in the station normal direction, with time, we need to determine number of points
    // hi : (t_max - t_min )/dt
    const amrex::Box diag_box(lo, hi);
    m_buffer_box = diag_box;
    amrex::BoxArray diag_ba;
    diag_ba.define(diag_box);
    amrex::BoxArray ba = diag_ba.maxSize(256);
    amrex::DistributionMapping dmap(ba);
    int ncomps = 6; //  Ex Ey Ez Bx By Bz
    //m_mf_output[0][lev] = amrex::MultiFab( amrex::convert(ba, amrex::IntVect::TheNodeVector()), dmap, ncomps, 1);
    // temporary
    m_mf_output[0][lev] = amrex::MultiFab( ba, dmap, ncomps, 1);
    m_mf_output[0][lev].setVal(0.);
}

void
StationDiagnostics::PrepareFieldDataForOutput ()
{
    // temporary
    auto & warpx = WarpX::GetInstance();
    for (int lev = 0; lev < 1; ++lev) {
        int icomp_dst = 0;
        for (int icomp = 0; icomp < 6; ++icomp) {
            m_cell_center_functors[lev][icomp]->operator()(*m_cell_centered_data[lev], icomp_dst);
            icomp_dst += m_cell_center_functors[lev][icomp]->nComp();
        }
        ablastr::utils::communication::FillBoundary(*m_cell_centered_data[lev],
                                                    WarpX::do_single_precision_comms,
                                                    warpx.Geom(lev).periodicity());
    }


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
    // reset counter
    m_tmax = warpx.gett_new(0);
    m_slice_counter = 0;
}

// temporary
void
StationDiagnostics::DerivedInitData ()
{
    m_cell_centered_data.resize(nmax_lev);
    m_cell_center_functors.resize(nmax_lev);
    auto & warpx = WarpX::GetInstance();
    const int ngrow = 1;
    const int ncomps = 6;
    for (int lev = 0; lev < nmax_lev; ++lev) {
        amrex::BoxArray ba = warpx.boxArray(lev);
        const amrex::DistributionMapping dmap = warpx.DistributionMap(lev);
        WarpX::AllocInitMultiFab( m_cell_centered_data[lev], ba, dmap, ncomps, amrex::IntVect(ngrow), lev, "cellcentered_Station", 0._rt);
    }
}
