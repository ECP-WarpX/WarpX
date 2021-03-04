/* Copyright 2019-2021 Yinjian Zhao, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleHistogram.H"
#include "WarpX.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Utils/WarpXUtil.H"
#include "Particles/Filter/FilterFunctors.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_ParticleReduce.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <limits>
#include <memory>


using namespace amrex;

struct NormalizationType {
    enum {
        no_normalization = 0,
        unity_particle_weight,
        max_to_unity,
        area_to_unity
    };
};

// constructor
ParticleHistogram::ParticleHistogram (std::string rd_name)
: ReducedDiags{rd_name}
{

    ParmParse pp(rd_name);

    // read species
    std::string selected_species_name;
    pp.get("species",selected_species_name);

    // read bin parameters
    pp.get("bin_number",m_bin_num);
    getWithParser(pp, "bin_max",   m_bin_max);
    getWithParser(pp, "bin_min",   m_bin_min);
    m_bin_size = (m_bin_max - m_bin_min) / m_bin_num;

    // read histogram function
    std::string function_string = "";
    Store_parserString(pp,"histogram_function(t,x,y,z,ux,uy,uz)",
                       function_string);
    m_parser = std::make_unique<ParserWrapper<m_nvars>>(
        makeParser(function_string,{"t","x","y","z","ux","uy","uz"}));

    // read normalization type
    std::string norm_string = "default";
    pp.query("normalization",norm_string);

    // set normalization type
    if ( norm_string == "default" ) {
        m_norm = NormalizationType::no_normalization;
    } else if ( norm_string == "unity_particle_weight" ) {
        m_norm = NormalizationType::unity_particle_weight;
    } else if ( norm_string == "max_to_unity" ) {
        m_norm = NormalizationType::max_to_unity;
    } else if ( norm_string == "area_to_unity" ) {
        m_norm = NormalizationType::area_to_unity;
    } else {
        Abort("Unknown ParticleHistogram normalization type.");
    }

    // get MultiParticleContainer class object
    const auto & mypc = WarpX::GetInstance().GetPartContainer();
    // get species names (std::vector<std::string>)
    auto const species_names = mypc.GetSpeciesNames();
    // select species
    for ( int i = 0; i < mypc.nSpecies(); ++i )
    {
        if ( selected_species_name == species_names[i] ){
            m_selected_species_id = i;
        }
    }
    // if m_selected_species_id is not modified
    if ( m_selected_species_id == -1 ){
        Abort("Unknown species for ParticleHistogram reduced diagnostic.");
    }

    // Read optional filter
    std::string buf;
    m_do_parser_filter = pp.query("filter_function(t,x,y,z,ux,uy,uz)", buf);
    if (m_do_parser_filter) {
        std::string filter_string = "";
        Store_parserString(pp,"filter_function(t,x,y,z,ux,uy,uz)", filter_string);
        m_parser_filter = std::make_unique<ParserWrapper<m_nvars>>(
                                     makeParser(filter_string,{"t","x","y","z","ux","uy","uz"}));
    }

    // resize data array
    m_data.resize(m_bin_num,0.0_rt);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs{m_path + m_rd_name + "." + m_extension, std::ofstream::out};
            // write header row
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";
            for (int i = 0; i < m_bin_num; ++i)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(3+i) + "]";
                Real b = m_bin_min + m_bin_size*(Real(i)+0.5_rt);
                ofs << "bin" + std::to_string(1+i)
                             + "=" + std::to_string(b) + "()";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that computes the histogram
void ParticleHistogram::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) return;

    // get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();

    // get time at level 0
    auto const t = warpx.gett_new(0);

    // get MultiParticleContainer class object
    const auto & mypc = warpx.GetPartContainer();

    // get WarpXParticleContainer class object
    auto & myspc = mypc.GetParticleContainer(m_selected_species_id);

    // get parser
    HostDeviceParser<m_nvars> fun_partparser = getParser(m_parser);

    // get filter parser
    HostDeviceParser<m_nvars> fun_filterparser = getParser(m_parser_filter);

    // declare local variables
    auto const num_bins = m_bin_num;
    Real const bin_min  = m_bin_min;
    Real const bin_size = m_bin_size;
    const bool is_unity_particle_weight =
        (m_norm == NormalizationType::unity_particle_weight) ? true : false;

    bool const do_parser_filter = m_do_parser_filter;

    // zero-out old data on the host
    std::fill(m_data.begin(), m_data.end(), amrex::Real(0.0));
    amrex::Gpu::DeviceVector< amrex::Real > d_data( m_data.size(), 0.0 );
    amrex::Real* const AMREX_RESTRICT dptr_data = d_data.dataPtr();

    int const nlevs = std::max(0, myspc.finestLevel()+1);
    for (int lev = 0; lev < nlevs; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
            for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti)
            {
                auto const GetPosition = GetParticlePosition(pti);

                auto & attribs = pti.GetAttribs();
                Real* const AMREX_RESTRICT d_w = attribs[PIdx::w].dataPtr();
                Real* const AMREX_RESTRICT d_ux = attribs[PIdx::ux].dataPtr();
                Real* const AMREX_RESTRICT d_uy = attribs[PIdx::uy].dataPtr();
                Real* const AMREX_RESTRICT d_uz = attribs[PIdx::uz].dataPtr();

                long const np = pti.numParticles();

                //Flag particles that need to be copied if they cross the z_slice
                amrex::ParallelFor(np,
                   [=] AMREX_GPU_DEVICE(int i)
                {
                    amrex::ParticleReal x, y, z;
                    GetPosition(i, x, y, z);
                    auto const w  = d_w[i];
                    auto const ux = d_ux[i] / PhysConst::c;
                    auto const uy = d_uy[i] / PhysConst::c;
                    auto const uz = d_uz[i] / PhysConst::c;

                    // don't count a particle if it is filtered out
                    if (do_parser_filter)
                        if (!fun_filterparser(t, x, y, z, ux, uy, uz))
                            return;
                    // continue function if particle is not filtered out
                    auto const f = fun_partparser(t, x, y, z, ux, uy, uz);

                    // determine particle bin
                    int const bin = int(Math::floor((f-bin_min)/bin_size));
                    if ( bin<0 || bin>=num_bins ) return; // discard if out-of-range

                    // add particle to histogram bin
                    //! @todo performance: on CPU, we are probably faster by
                    //        letting each thread compute its own histogram and
                    //        then we reduce the histograms after the loop
                    if ( is_unity_particle_weight ) {
                        amrex::HostDevice::Atomic::Add(&dptr_data[bin], 1.0_rt);
                    } else {
                        amrex::HostDevice::Atomic::Add(&dptr_data[bin], w);
                    }
                });
            }
        }
    }

    // blocking copy from device to host
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
        d_data.begin(), d_data.end(), m_data.begin());

    // reduced sum over mpi ranks
    ParallelDescriptor::ReduceRealSum
        (m_data.data(), m_data.size(), ParallelDescriptor::IOProcessorNumber());

    // normalize the maximum value to be one
    if ( m_norm == NormalizationType::max_to_unity )
    {
        Real f_max = 0.0_rt;
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( m_data[i] > f_max ) f_max = m_data[i];
        }
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( f_max > std::numeric_limits<Real>::min() ) m_data[i] /= f_max;
        }
        return;
    }

    // normalize the area (integral) to be one
    if ( m_norm == NormalizationType::area_to_unity )
    {
        Real f_area = 0.0_rt;
        for ( int i = 0; i < m_bin_num; ++i )
        {
            f_area += m_data[i] * m_bin_size;
        }
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( f_area > std::numeric_limits<Real>::min() ) m_data[i] /= f_area;
        }
        return;
    }

}
// end void ParticleHistogram::ComputeDiags
