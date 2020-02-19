/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Histogram.H"
#include "WarpX.H"
#include "WarpXConst.H"
#include <WarpXUtil.H>
#include "AMReX_REAL.H"
#include "AMReX_ParticleReduce.H"

using namespace amrex;

// constructor
Histogram::Histogram (std::string rd_name)
: ReducedDiags{rd_name}
{

    // read parameters
    ParmParse pp(rd_name);
    pp.get("species",m_species_name);
    pp.get("bin_number",m_bin_num);
    pp.get("bin_max",   m_bin_max);
    pp.get("bin_min",   m_bin_min);
    m_bin_size = (m_bin_max - m_bin_min) / m_bin_num;
    std::string function_string = "";
    Store_parserString(pp,"histogram_function(t,x,y,z,ux,uy,uz)",
                       function_string);
    m_parser.reset(new ParserWrapper<7>(
        makeParser(function_string,{"t","x","y","z","ux","uy","uz"})));

    // resize data array
    m_data.resize(m_bin_num+2,0.0);

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
void Histogram::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get MultiParticleContainer class object
    auto & mypc = warpx.GetPartContainer();

    // get number of species (int)
    auto nSpecies = mypc.nSpecies();

    // get species names (std::vector<std::string>)
    auto species_names = mypc.GetSpeciesNames();

    // get time at level 0
    auto t = warpx.gett_new(0);

    // select species
    int i_s;
    for ( int i = 0; i < nSpecies; ++i )
    {
        if ( m_species_name == species_names[i] ) { i_s = i; }
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
            auto w  = p.rdata(PIdx::w);
            auto x  = p.pos(0);
            auto y  = p.pos(1);
            auto z  = p.pos(2);
            auto ux = p.rdata(PIdx::ux);
            auto uy = p.rdata(PIdx::uy);
            auto uz = p.rdata(PIdx::uz);
            auto f = (*fun_partparser)(t,x,y,z,ux,uy,uz);
            auto f1 = m_bin_min + m_bin_size*i;
            auto f2 = m_bin_min + m_bin_size*(i+1);
            if ( f > f1 && f < f2 )
            { return w; }
            else
            { return 0.0; }
        });
        // reduced sum over mpi ranks
        ParallelDescriptor::ReduceRealSum
            (m_data[i], ParallelDescriptor::IOProcessorNumber());
    }

}
// end void Histogram::ComputeDiags
