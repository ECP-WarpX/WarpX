/* Copyright 2019-2020 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "Laser/LaserProfiles.H"

#include "Utils/TextMsg.H"
#include "Utils/WarpXUtil.H"
#include "Utils/WarpX_Complex.H"

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

using namespace amrex;

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::init (
    const amrex::ParmParse& ppl,
    const amrex::ParmParse& /* ppc */,
    CommonLaserParameters params)
{
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
    if(queryWithParser(ppl ,"time_chunk_size", temp)){
        m_params.time_chunk_size = min(
            temp, m_params.time_chunk_size);
    }
    if(m_params.time_chunk_size < 2){
        Abort("Error! time_chunk_size must be >= 2!");
    }

    //Reads the (optional) delay
    queryWithParser(ppl, "delay", m_params.t_delay);

    //Allocate memory for E_data Vector
    const int data_size = m_params.time_chunk_size*
            m_params.nx*m_params.ny;
    m_params.E_data.resize(data_size);

    //Read first time chunck
    read_data_t_chuck(0, m_params.time_chunk_size);

    //Copy common params
    m_common_params = params;
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::update (amrex::Real t)
{
    t -= m_params.t_delay;

    if(t >= m_params.t_coords.back())
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
    t -= m_params.t_delay;

    //Amplitude is 0 if time is out of range
    if(t < m_params.t_coords.front() ||  t > m_params.t_coords.back()){
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

    if(m_params.is_grid_uniform){
        internal_fill_amplitude_uniform(
            idx_t_left, np, Xp, Yp, t, amplitude);
    }
    else{
        internal_fill_amplitude_nonuniform(
            idx_t_left, np, Xp, Yp, t, amplitude);
    }
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::parse_txye_file(std::string txye_file_name)
{
    if(ParallelDescriptor::IOProcessor()){
        std::ifstream inp(txye_file_name, std::ios::binary);
        if(!inp) Abort("Failed to open txye file");
        inp.exceptions(std::ios_base::failbit | std::ios_base::badbit);

        //Uniform grid flag
        char flag;
        inp.read(&flag, 1);
        if(!inp) Abort("Failed to read grid type from txye file");
        m_params.is_grid_uniform=flag;

        //Grid points along t, x and y
        inp.read(reinterpret_cast<char*>(&m_params.nt), sizeof(uint32_t));
        inp.read(reinterpret_cast<char*>(&m_params.nx), sizeof(uint32_t));
        inp.read(reinterpret_cast<char*>(&m_params.ny), sizeof(uint32_t));
        if(!inp) Abort("Failed to read sizes from txye file");

        if(m_params.nt <= 1) Abort("nt in txye file must be >=2");
        if(m_params.nx <= 1) Abort("nx in txye file must be >=2");
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
        if(m_params.ny <= 1) Abort("ny in txye file must be >=2 in 3D");
#elif defined(WARPX_DIM_XZ)
        if(m_params.ny != 1) Abort("ny in txye file must be 1 in 2D");
#endif

        //Coordinates
        Vector<double> dbuf_t, dbuf_x, dbuf_y;
        if(m_params.is_grid_uniform){
            dbuf_t.resize(2);
            dbuf_x.resize(2);
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
            dbuf_y.resize(2);
#elif defined(WARPX_DIM_XZ)
            dbuf_y.resize(1);
#endif
        }
        else{
            dbuf_t.resize(m_params.nt);
            dbuf_x.resize(m_params.nx);
            dbuf_y.resize(m_params.ny);
        }

        inp.read(reinterpret_cast<char*>(dbuf_t.dataPtr()),
            dbuf_t.size()*sizeof(double));
        inp.read(reinterpret_cast<char*>(dbuf_x.dataPtr()),
            dbuf_x.size()*sizeof(double));
        inp.read(reinterpret_cast<char*>(dbuf_y.dataPtr()),
            dbuf_y.size()*sizeof(double));
        if(!inp) Abort("Failed to read coords from txye file");

        m_params.t_coords.resize(dbuf_t.size());
        m_params.h_x_coords.resize(dbuf_x.size());
        m_params.h_y_coords.resize(dbuf_y.size());

        if (!std::is_sorted(dbuf_t.begin(), dbuf_t.end()) ||
            !std::is_sorted(dbuf_x.begin(), dbuf_x.end()) ||
            !std::is_sorted(dbuf_y.begin(), dbuf_y.end()))
            Abort("Coordinates are not sorted  in txye file");

        // Convert from double to amrex::Real
        std::transform(dbuf_t.begin(), dbuf_t.end(), m_params.t_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
        std::transform(dbuf_x.begin(), dbuf_x.end(), m_params.h_x_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
        std::transform(dbuf_y.begin(), dbuf_y.end(), m_params.h_y_coords.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
    }

    //Broadcast grid uniformity
    char is_grid_uniform = m_params.is_grid_uniform;
    ParallelDescriptor::Bcast(&is_grid_uniform, 1,
        ParallelDescriptor::IOProcessorNumber());
    m_params.is_grid_uniform = is_grid_uniform;

    //Broadcast grid size and coordinate sizes
    //When a non-uniform grid is used, nt, nx and ny are identical
    //to t_coords.size(), x_coords.size() and y_coords.size().
    //When a uniform grid is used, nt,nx and ny store the number of points
    //used for the mesh, while t_coords, x_coords and y_coords store the
    //extrems in each direaction. Thus t_coords and x_coords in this case
    //have size 2 and y_coords has size 1 in 2D and size 2 in 3D.
    int t_sizes[6] = {m_params.nt, m_params.nx, m_params.ny,
        static_cast<int>(m_params.t_coords.size()),
        static_cast<int>(m_params.h_x_coords.size()),
        static_cast<int>(m_params.h_y_coords.size())};
    ParallelDescriptor::Bcast(t_sizes, 6,
        ParallelDescriptor::IOProcessorNumber());
    m_params.nt = t_sizes[0]; m_params.nx = t_sizes[1]; m_params.ny = t_sizes[2];

    //Broadcast coordinates
    if(!ParallelDescriptor::IOProcessor()){
        m_params.t_coords.resize(t_sizes[3]);
        m_params.h_x_coords.resize(t_sizes[4]);
        m_params.h_y_coords.resize(t_sizes[5]);
    }
    ParallelDescriptor::Bcast(m_params.t_coords.dataPtr(),
        m_params.t_coords.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(m_params.h_x_coords.dataPtr(),
        m_params.h_x_coords.size(), ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(m_params.h_y_coords.dataPtr(),
        m_params.h_y_coords.size(), ParallelDescriptor::IOProcessorNumber());

    m_params.d_x_coords.resize(m_params.h_x_coords.size());
    m_params.d_y_coords.resize(m_params.h_y_coords.size());
    Gpu::copyAsync(Gpu::hostToDevice,
                   m_params.h_x_coords.begin(), m_params.h_x_coords.end(),
                   m_params.d_x_coords.begin());
    Gpu::copyAsync(Gpu::hostToDevice,
                   m_params.h_y_coords.begin(), m_params.h_y_coords.end(),
                   m_params.d_y_coords.begin());
    Gpu::synchronize();
}

std::pair<int,int>
WarpXLaserProfiles::FromTXYEFileLaserProfile::find_left_right_time_indices(amrex::Real t) const
{
    int idx_t_right;
    if(m_params.is_grid_uniform){
        const auto t_min = m_params.t_coords.front();
        const auto t_max = m_params.t_coords.back();
        const auto temp_idx_t_right = static_cast<int>(
            std::ceil( (m_params.nt-1)*(t-t_min)/(t_max-t_min)));
        idx_t_right = max(min(temp_idx_t_right, m_params.nt-1),1);
    }
    else{
        idx_t_right = std::distance(m_params.t_coords.begin(),
        std::upper_bound(m_params.t_coords.begin(),
            m_params.t_coords.end(), t));
    }
    return std::make_pair(idx_t_right-1, idx_t_right);
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::read_data_t_chuck(int t_begin, int t_end)
{
    amrex::Print() << Utils::TextMsg::Info(
        "Reading [" + std::to_string(t_begin) + ", " + std::to_string(t_end) +
        ") data chunk from " + m_params.txye_file_name);

    //Indices of the first and last timestep to read
    auto i_first = max(0, t_begin);
    auto i_last = min(t_end-1, m_params.nt-1);
    if(i_last-i_first+1 > static_cast<int>(m_params.E_data.size()))
        Abort("Data chunk to read from file is too large");

    Vector<Real> h_E_data(m_params.E_data.size());

    if(ParallelDescriptor::IOProcessor()){
        //Read data chunk
        std::ifstream inp(m_params.txye_file_name, std::ios::binary);
        if(!inp) Abort("Failed to open txye file");
        inp.exceptions(std::ios_base::failbit | std::ios_base::badbit);

        auto skip_amount = 1 +
            3*sizeof(uint32_t) +
            m_params.t_coords.size()*sizeof(double) +
            m_params.h_x_coords.size()*sizeof(double) +
            m_params.h_y_coords.size()*sizeof(double) +
            sizeof(double)*t_begin*m_params.nx*m_params.ny;
        inp.seekg(skip_amount);
        if(!inp) Abort("Failed to read field data from txye file");
        const int read_size = (i_last - i_first + 1)*
            m_params.nx*m_params.ny;
        Vector<double> buf_e(read_size);
        inp.read(reinterpret_cast<char*>(buf_e.dataPtr()), read_size*sizeof(double));
        if(!inp) Abort("Failed to read field data from txye file");
        std::transform(buf_e.begin(), buf_e.end(), h_E_data.begin(),
            [](auto x) {return static_cast<amrex::Real>(x);} );
    }

    //Broadcast E_data
    ParallelDescriptor::Bcast(h_E_data.dataPtr(),
        h_E_data.size(), ParallelDescriptor::IOProcessorNumber());

    Gpu::copyAsync(Gpu::hostToDevice,h_E_data.begin(),h_E_data.end(),m_params.E_data.begin());
    Gpu::synchronize();

    //Update first and last indices
    m_params.first_time_index = i_first;
    m_params.last_time_index = i_last;
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
    const auto tmp_e_max = m_common_params.e_max;
    const auto tmp_x_min = m_params.h_x_coords.front();
    const auto tmp_x_max = m_params.h_x_coords.back();
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    const auto tmp_y_min = m_params.h_y_coords.front();
    const auto tmp_y_max = m_params.h_y_coords.back();
#endif
    const auto tmp_nx = m_params.nx;
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    const auto tmp_ny = m_params.ny;
#endif
    const auto p_E_data = m_params.E_data.dataPtr();
    const auto tmp_idx_first_time = m_params.first_time_index;
    const int idx_t_right = idx_t_left+1;
    const auto t_left = idx_t_left*
        (m_params.t_coords.back()-m_params.t_coords.front())/(m_params.nt-1) +
        m_params.t_coords.front();
    const auto t_right = idx_t_right*
        (m_params.t_coords.back()-m_params.t_coords.front())/(m_params.nt-1) +
        m_params.t_coords.front();

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
        amplitude[i] = WarpXUtilAlgo::trilinear_interp(
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
            t, Xp[i], Yp[i])*tmp_e_max;

#elif defined(WARPX_DIM_XZ)
        //Interpolate amplitude
        const auto idx = [=](int i_interp, int j_interp){
            return (i_interp-tmp_idx_first_time) * tmp_nx + j_interp;
        };
        amplitude[i] = WarpXUtilAlgo::bilinear_interp(
            t_left, t_right,
            x_0, x_1,
            p_E_data[idx(idx_t_left, idx_x_left)],
            p_E_data[idx(idx_t_left, idx_x_right)],
            p_E_data[idx(idx_t_right, idx_x_left)],
            p_E_data[idx(idx_t_right, idx_x_right)],
            t, Xp[i])*tmp_e_max;
        amrex::ignore_unused(Yp);
#else
        // TODO: implement WARPX_DIM_1D_Z
        amrex::ignore_unused(x_0, x_1, tmp_e_max, p_E_data, tmp_idx_first_time,
                             t_left, t_right, Xp, Yp, t, idx_x_left);
#endif
        }
    );
}

void
WarpXLaserProfiles::FromTXYEFileLaserProfile::internal_fill_amplitude_nonuniform(
    const int idx_t_left,
    const int np,
    Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude) const
{
    // Copy member variables to tmp copies
    // and get pointers to underlying data for GPU.
    const auto tmp_e_max = m_common_params.e_max;
    const auto tmp_x_min = m_params.h_x_coords.front();
    const auto tmp_x_max = m_params.h_x_coords.back();
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    const auto tmp_y_min = m_params.h_y_coords.front();
    const auto tmp_y_max = m_params.h_y_coords.back();
#endif
    const auto p_x_coords = m_params.d_x_coords.dataPtr();
    const int tmp_x_coords_size = static_cast<int>(m_params.d_x_coords.size());
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
    const auto p_y_coords = m_params.d_y_coords.dataPtr();
    const int tmp_y_coords_size = static_cast<int>(m_params.d_y_coords.size());
#endif
    const auto p_E_data = m_params.E_data.dataPtr();
    const auto tmp_idx_first_time = m_params.first_time_index;
    const int idx_t_right = idx_t_left+1;
    const auto t_left = m_params.t_coords[idx_t_left];
    const auto t_right = m_params.t_coords[idx_t_right];

#if (defined WARPX_DIM_1D_Z)
    amrex::Abort(Utils::TextMsg::Err(
        "WarpXLaserProfiles::FromTXYEFileLaserProfile Not implemented for 1D"));
#endif

    // Loop through the macroparticle to calculate the proper amplitude
    amrex::ParallelFor(
    np,
    [=] AMREX_GPU_DEVICE (int ip) {
        //Amplitude is zero if we are out of bounds
        if (Xp[ip] <= tmp_x_min || Xp[ip] >= tmp_x_max){
            amplitude[ip] = 0.0_rt;
            return;
        }
#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
        if (Yp[ip] <= tmp_y_min || Yp[ip] >= tmp_y_max){
            amplitude[ip] = 0.0_rt;
            return;
        }
#else
        amrex::ignore_unused(Yp);
#endif

        //Find indices along x
        auto const p_x_right = WarpXUtilAlgo::upper_bound(
                p_x_coords, p_x_coords+tmp_x_coords_size, Xp[ip]);
        const int idx_x_right = p_x_right - p_x_coords;
        const int idx_x_left = idx_x_right - 1;

#if (defined(WARPX_DIM_3D) || (defined WARPX_DIM_RZ))
        //Find indices along y
        auto const p_y_right = WarpXUtilAlgo::upper_bound(
            p_y_coords, p_y_coords+tmp_y_coords_size, Yp[ip]);
        const int idx_y_right = p_y_right - p_y_coords;
        const int idx_y_left = idx_y_right - 1;

        //Interpolate amplitude
        const auto idx = [=](int i, int j, int k){
            return
                (i-tmp_idx_first_time)*tmp_x_coords_size*tmp_y_coords_size+
                j*tmp_y_coords_size + k;
        };
        amplitude[ip] = WarpXUtilAlgo::trilinear_interp(
            t_left, t_right,
            p_x_coords[idx_x_left], p_x_coords[idx_x_right],
            p_y_coords[idx_y_left], p_y_coords[idx_y_right],
            p_E_data[idx(idx_t_left, idx_x_left, idx_y_left)],
            p_E_data[idx(idx_t_left, idx_x_left, idx_y_right)],
            p_E_data[idx(idx_t_left, idx_x_right, idx_y_left)],
            p_E_data[idx(idx_t_left, idx_x_right, idx_y_right)],
            p_E_data[idx(idx_t_right, idx_x_left, idx_y_left)],
            p_E_data[idx(idx_t_right, idx_x_left, idx_y_right)],
            p_E_data[idx(idx_t_right, idx_x_right, idx_y_left)],
            p_E_data[idx(idx_t_right, idx_x_right, idx_y_right)],
            t, Xp[ip], Yp[ip])*tmp_e_max;

#elif (defined WARPX_DIM_XZ)
        //Interpolate amplitude
        const auto idx = [=](int i, int j){
            return (i-tmp_idx_first_time) * tmp_x_coords_size + j;
        };
        amplitude[ip] = WarpXUtilAlgo::bilinear_interp(
            t_left, t_right,
            p_x_coords[idx_x_left], p_x_coords[idx_x_right],
            p_E_data[idx(idx_t_left, idx_x_left)],
            p_E_data[idx(idx_t_left, idx_x_right)],
            p_E_data[idx(idx_t_right, idx_x_left)],
            p_E_data[idx(idx_t_right, idx_x_right)],
            t, Xp[ip])*tmp_e_max;
#else
        // TODO: implement WARPX_DIM_1D_Z
        amrex::ignore_unused(idx_x_left, idx_t_left, idx_t_right, tmp_e_max,
                             p_E_data, tmp_idx_first_time, t_left, t_right, t);
#endif
        }
    );
}
