/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ParticleHistogram.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include "WarpXUtil.H"
#include "AMReX_REAL.H"
#include "AMReX_ParticleReduce.H"
#include <limits>

using namespace amrex;

// constructor
ParticleHistogram::ParticleHistogram (std::string rd_name)
: ReducedDiags{rd_name}
{

    // read parameters
    ParmParse pp(rd_name);
    pp.get("species",m_species_name);
    pp.get("bin_number",m_bin_num);
    pp.get("bin_max",   m_bin_max);
    pp.get("bin_min",   m_bin_min);
    pp.query("normalization",m_norm);
    m_bin_size = (m_bin_max - m_bin_min) / m_bin_num;
    std::string function_string = "";
    Store_parserString(pp,"histogram_function(t,x,y,z,ux,uy,uz)",
                       function_string);
    m_parser.reset(new ParserWrapper<7>(
        makeParser(function_string,{"t","x","y","z","ux","uy","uz"})));

    // resize data array
    m_data.resize(m_bin_num,0.0);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs;
            ofs.open(m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app);
            // write header row
            ofs << "#";
            ofs << "histogram_function(t,x,y,z,ux,uy,uz) = ";
            ofs << function_string << std::endl;
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";
            for (int i = 0; i < m_bin_num; ++i)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(3+i) + "]";
                ofs << "bin" + std::to_string(1+i) + "()";
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
    if ( (step+1) % m_freq != 0 ) { return; }

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();

    // get number of species (int)
    auto const nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    auto const species_names = mypc.GetSpeciesNames();

    // get time at level 0
    auto const t = warpx.gett_new(0);

    // select species
    int i_s;
    for ( int i = 0; i < nSpecies; ++i )
    {
        if ( m_species_name == species_names[i] ) i_s = i;
    }

    // get WarpXParticleContainer class object
    auto & myspc = mypc.GetParticleContainer(i_s);

    using PType = typename WarpXParticleContainer::SuperParticleType;

    // get parser
    ParserWrapper<7> *fun_partparser = m_parser.get();

    for ( int i = 0; i < m_bin_num; ++i )
    {
        // compute the histogram
        m_data[i] = ReduceSum( myspc,
        [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
        {
            auto const w  = p.rdata(PIdx::w);
            auto const x  = p.pos(0);
            auto const y  = p.pos(1);
            auto const z  = p.pos(2);
            auto const ux = p.rdata(PIdx::ux);
            auto const uy = p.rdata(PIdx::uy);
            auto const uz = p.rdata(PIdx::uz);
            auto const f = (*fun_partparser)(t,x,y,z,ux,uy,uz);
            auto const f1 = m_bin_min + m_bin_size*i;
            auto const f2 = m_bin_min + m_bin_size*(i+1);
            if ( f > f1 && f < f2 ) {
                if ( m_norm == "unity_particle_weight" )
                    return 1.0;
                else
                    return w;
            } else
                return 0.0;
        });
        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum
            (m_data[i], ParallelDescriptor::IOProcessorNumber());
    }

    if ( m_norm == "max_to_unity" )
    {
        Real f_max = 0.0;
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( m_data[i] > f_max )
                f_max = m_data[i];
        }
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( f_max > std::numeric_limits<Real>::min() )
                m_data[i] /= f_max;
        }
    }
    else if ( m_norm == "area_to_unity" )
    {
        Real f_area = 0.0;
        for ( int i = 0; i < m_bin_num; ++i )
        {
            f_area += m_data[i] * m_bin_size;
        }
        for ( int i = 0; i < m_bin_num; ++i )
        {
            if ( f_area > std::numeric_limits<Real>::min() )
                m_data[i] /= f_area;
        }
    }
    else if ( m_norm == "default" )
    { /* do nothing */ }
    else if ( m_norm == "unity_particle_weight" )
    { /* do nothing */ }
    else { Abort("Unknown ParticleHistogram normalization type."); }

}
// end void ParticleHistogram::ComputeDiags
