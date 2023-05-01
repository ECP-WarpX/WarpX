/* Copyright 2019-2020 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Laser/LaserProfiles.H"

#include "Utils/Algorithms/LinearInterpolation.H"
#include "Utils/Algorithms/UpperBound.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpX_Complex.H"
#include "Utils/WarpXConst.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iterator>
#include <limits>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// example: data handling & print
#include <vector>   // std::vector
#include <iostream> // std::cout
#include <memory>   // std::shared_ptr

#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
    namespace io = openPMD;
#endif

using namespace amrex;

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::init (
    const amrex::ParmParse& ppl,
    CommonLaserParameters params)
{
#ifdef WARPX_USE_OPENPMD
    if (!std::numeric_limits< double >::is_iec559)
    {
        ablastr::warn_manager::WMRecordWarning("Laser",
            "(Double does not comply with IEEE 754: bad"
            "things will happen parsing the X, Y and T profiles for the laser!)",
            ablastr::warn_manager::WarnPriority::high);
    }
    // Parse the TXYE file
    ppl.get("txye_file_name", m_params.txye_file_name);
    if(m_params.txye_file_name.empty())
    {
        Abort("txye_file_name must be provided for txye_file laser profile!");
    }
    parse_txye_file(m_params.txye_file_name);
    //Set time_chunk_size
    m_params.time_chunk_size = m_params.nt;
    int temp = 1;
    if(utils::parser::queryWithParser(ppl ,"time_chunk_size", temp)){
        m_params.time_chunk_size = min(
            temp, m_params.time_chunk_size);
    }
    if(m_params.time_chunk_size < 2){
        Abort("Error! time_chunk_size must be >= 2!");
    }
    //Reads the (optional) delay
    utils::parser::queryWithParser(ppl, "delay", m_params.t_delay);

    //Allocate memory for E_data Vector
    const int data_size = m_params.time_chunk_size*
            m_params.nx*m_params.ny;
    m_params.E_data.resize(data_size);

    //Read first time chunck
    read_data_t_chuck(0, m_params.time_chunk_size);
    //Copy common params
    m_common_params = params;
#else
    amrex::Abort(Utils::TextMsg::Err("WarpX has to be compiled with option openPMD=ON to read lasy file"));
    amrex::ignore_unused(ppl, params);
#endif
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::update (amrex::Real t)
{
    t += m_params.t_min - m_params.t_delay;
    if(t >= m_params.t_max)
        return;
    const auto idx_times = find_left_right_time_indices(t);
    const auto idx_t_left = idx_times.first;
    const auto idx_t_right = idx_times.second;
    //Load data chunck if needed
    if(idx_t_right >  m_params.last_time_index){
        read_data_t_chuck(idx_t_left, idx_t_left+m_params.time_chunk_size);
    }
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::fill_amplitude (
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    t += m_params.t_min - m_params.t_delay;
    //Amplitude is 0 if time is out of range
    if(t < m_params.t_min ||  t > m_params.t_max){
        amrex::ParallelFor(np,
            [=] AMREX_GPU_DEVICE (int i) {
                amplitude[i] = 0.0_rt;});
        return;
    }
    //Find left and right time indices
    const auto [idx_t_left, idx_t_right] = find_left_right_time_indices(t);

    if(idx_t_left <  m_params.first_time_index){
        Abort("Something bad has happened with the simulation time");
    }
    internal_fill_amplitude_uniform(idx_t_left, np, Xp, Yp, t, amplitude);
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::parse_txye_file(std::string txye_file_name)
{
#ifdef WARPX_USE_OPENPMD
    if(ParallelDescriptor::IOProcessor()){
        auto series = io::Series(txye_file_name, io::Access::READ_ONLY);
        auto i = series.iterations[0];
        auto E = i.meshes["laserEnvelope"];
        auto E_laser = E[io::RecordComponent::SCALAR];
        auto extent = E_laser.getExtent();
        m_params.nt = extent[0];
        m_params.nx = extent[1];
        m_params.ny = extent[2];
        if(m_params.nt <= 1) Abort("nt in txye file must be >=2");
        if(m_params.nx <= 1) Abort("nx in txye file must be >=2");
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
        if(m_params.ny <= 1) Abort("ny in txye file must be >=2 in 3D");
#elif defined(WARPX_DIM_XZ)
        if(m_params.ny != 1) Abort("ny in txye file must be 1 in 2D");
#endif
        // Extract grid offset and grid spacing
        std::vector<double> offset = E.gridGlobalOffset();
        std::vector<double> position = E_laser.position<double>();
        std::vector<double> spacing = E.gridSpacing<double>();
        // Calculate the min and max of the grid
        m_params.t_min = offset[0] + position[0]*spacing[0];
        m_params.t_max = m_params.t_min + (m_params.nt-1)*spacing[0];
        m_params.x_min = offset[1] + position[1]*spacing[1];
        m_params.x_max = m_params.x_min + (m_params.nx-1)*spacing[1];
        m_params.y_min = offset[2] + position[2]*spacing[2];
        m_params.y_max = m_params.y_min + (m_params.ny-1)*spacing[2];
    }
#else
    amrex::ignore_unused(txye_file_name);
#endif
}

std::pair<int,int>
WarpXLaserProfiles::FromTXYEFileLaserProfile::find_left_right_time_indices(amrex::Real t) const
{
    int idx_t_right;
    const auto t_min = m_params.t_min;
    const auto t_max = m_params.t_max;
    const auto temp_idx_t_right = static_cast<int>(std::ceil( (m_params.nt-1)*(t-t_min)/(t_max-t_min)));
    idx_t_right = max(min(temp_idx_t_right, m_params.nt-1),1);
    return std::make_pair(idx_t_right-1, idx_t_right);
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::read_data_t_chuck(int t_begin, int t_end)
{
#ifdef WARPX_USE_OPENPMD
    amrex::Print() << Utils::TextMsg::Info(
        "Reading [" + std::to_string(t_begin) + ", " + std::to_string(t_end) +
        ") data chunk from " + m_params.txye_file_name);
    //Indices of the first and last timestep to read
    std::uint64_t const i_first = max(0, t_begin);
    std::uint64_t const i_last = min(t_end-1, m_params.nt-1);
    auto const time_chunk_size = static_cast<std::uint64_t>( m_params.time_chunk_size);
    if((i_last-i_first+1)*m_params.nx*m_params.ny > static_cast<std::uint64_t>(m_params.E_data.size()))
        Abort("Data chunk to read from file is too large");
    Vector<Complex> h_E_data(m_params.E_data.size());
    if(ParallelDescriptor::IOProcessor()){
        auto series = io::Series(m_params.txye_file_name, io::Access::READ_ONLY);
        auto i = series.iterations[0];
        auto E = i.meshes["laserEnvelope"];
        auto E_laser = E[io::RecordComponent::SCALAR];
        openPMD:: Extent full_extent = E_laser.getExtent();
        openPMD::Extent read_extent = {full_extent[0], full_extent[1],(i_last - i_first + 1)};
        auto x_data = E_laser.loadChunk< std::complex<double> >(io::Offset{0, 0, i_first}, read_extent);
        const int read_size = (i_last - i_first + 1)*m_params.nx*m_params.ny;
        series.flush();
        for (int j=0; j<read_size; j++) {
            h_E_data[j] = Complex{ x_data.get()[j].real(), x_data.get()[j].imag() };
        }
    }
    //Broadcast E_data
    ParallelDescriptor::Bcast(h_E_data.dataPtr(),
        h_E_data.size(), ParallelDescriptor::IOProcessorNumber());

    Gpu::copyAsync(Gpu::hostToDevice,h_E_data.begin(),h_E_data.end(),m_params.E_data.begin());
    Gpu::synchronize();
    //Update first and last indices
    m_params.first_time_index = i_first;
    m_params.last_time_index = i_last;
#else
    amrex::ignore_unused(t_begin, t_end);
#endif
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::internal_fill_amplitude_uniform(
    const int idx_t_left,
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    // Copy member variables to tmp copies
    // and get pointers to underlying data for GPU.
    const amrex::Real omega_t = 2.*MathConst::pi*PhysConst::c*t/m_common_params.wavelength;
    const Complex exp_omega_t = Complex{ std::cos(-omega_t), std::sin(-omega_t) };
    const auto tmp_x_min = m_params.x_min;
    const auto tmp_x_max = m_params.x_max;
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    const auto tmp_y_min = m_params.y_min;
    const auto tmp_y_max = m_params.y_max;
#endif
    const auto tmp_nx = m_params.nx;
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    const auto tmp_ny = m_params.ny;
#endif
    const auto p_E_data = m_params.E_data.dataPtr();
    const auto tmp_idx_first_time = m_params.first_time_index;
    const int idx_t_right = idx_t_left+1;
    const auto t_left = idx_t_left*
        (m_params.t_max-m_params.t_min)/(m_params.nt-1) +
        m_params.t_min;
    const auto t_right = idx_t_right*
        (m_params.t_max-m_params.t_min)/(m_params.nt-1) +
        m_params.t_min;
#if (defined WARPX_DIM_1D_Z)
    amrex::Abort(Utils::TextMsg::Err(
        "WarpXLaserProfiles::FromTXYEFileLaserProfile Not implemented for 1D"));
#endif
    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
    np,
    [=] AMREX_GPU_DEVICE (int i) {
        //Amplitude is zero if we are out of bounds
        if (Xp[i] <= tmp_x_min || Xp[i] >= tmp_x_max){
            amplitude[i] = 0.0_rt;
            return;
        }
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
        if (Yp[i] <= tmp_y_min || Yp[i] >= tmp_y_max){
            amplitude[i] = 0.0_rt;
            return;
        }
#endif
        //Find indices and coordinates along x
        const int temp_idx_x_right = static_cast<int>(
            std::ceil((tmp_nx-1)*(Xp[i]- tmp_x_min)/(tmp_x_max-tmp_x_min)));
        const int idx_x_right =
            max(min(temp_idx_x_right,tmp_nx-1),static_cast<int>(1));
        const int idx_x_left = idx_x_right - 1;
        const auto x_0 =
            idx_x_left*(tmp_x_max-tmp_x_min)/(tmp_nx-1) + tmp_x_min;
        const auto x_1 =
            idx_x_right*(tmp_x_max-tmp_x_min)/(tmp_nx-1) + tmp_x_min;
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
        //Find indices and coordinates along y
        const int temp_idx_y_right = static_cast<int>(
            std::ceil((tmp_ny-1)*(Yp[i]- tmp_y_min)/(tmp_y_max-tmp_y_min)));
        const int idx_y_right =
            max(min(temp_idx_y_right,tmp_ny-1),static_cast<int>(1));
        const int idx_y_left = idx_y_right - 1;
        const auto y_0 =
            idx_y_left*(tmp_y_max-tmp_y_min)/(tmp_ny-1) + tmp_y_min;
        const auto y_1 =
            idx_y_right*(tmp_y_max-tmp_y_min)/(tmp_ny-1) + tmp_y_min;

        //Interpolate amplitude
        const auto idx = [=](int i_interp, int j_interp, int k_interp){
            return
                (i_interp-tmp_idx_first_time)*tmp_nx*tmp_ny+
                j_interp*tmp_ny + k_interp;
        };
        Complex val = utils::algorithms::trilinear_interp(
            t_left, t_right,
            x_0, x_1,
            y_0, y_1,
            p_E_data[idx(idx_t_left, idx_x_left, idx_y_left)],
            p_E_data[idx(idx_t_left, idx_x_left, idx_y_right)],
            p_E_data[idx(idx_t_left, idx_x_right, idx_y_left)],
            p_E_data[idx(idx_t_left, idx_x_right, idx_y_right)],
            p_E_data[idx(idx_t_right, idx_x_left, idx_y_left)],
            p_E_data[idx(idx_t_right, idx_x_left, idx_y_right)],
            p_E_data[idx(idx_t_right, idx_x_right, idx_y_left)],
            p_E_data[idx(idx_t_right, idx_x_right, idx_y_right)],
            t, Xp[i], Yp[i]);
#elif defined(WARPX_DIM_XZ)
        //Interpolate amplitude
        const auto idx = [=](int i_interp, int j_interp){
            return (i_interp-tmp_idx_first_time) * tmp_nx + j_interp;
        };
        Complex val = utils::algorithms::bilinear_interp(
            t_left, t_right,
            x_0, x_1,
            p_E_data[idx(idx_t_left, idx_x_left)],
            p_E_data[idx(idx_t_left, idx_x_right)],
            p_E_data[idx(idx_t_right, idx_x_left)],
            p_E_data[idx(idx_t_right, idx_x_right)],
            t, Xp[i]);
        amrex::ignore_unused(Yp);
#else
        // TODO: implement WARPX_DIM_1D_Z
        Complex val = 0;
        amrex::ignore_unused(x_0, x_1, p_E_data, tmp_idx_first_time,
                             t_left, t_right, Xp, Yp, t, idx_x_left);
#endif
        // The interpolated amplitude was only the envelope.
        // Here we add the laser oscillations.
        amplitude[i] = (val*exp_omega_t).real();
        }
    );
}
