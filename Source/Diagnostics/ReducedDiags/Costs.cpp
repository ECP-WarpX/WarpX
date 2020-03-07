/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Costs.H"
#include "WarpX.H"
#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
#include <AMReX_ParticleReduce.H>

#include <iostream>
#include <cmath>


using namespace amrex;

// constructor
Costs::Costs (std::string rd_name)
: ReducedDiags{rd_name}
{

    // RZ coordinate is not working
    #if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "Costs reduced diagnostics does not work for RZ coordinate.");
    #endif

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // read number of levels
    int nLevels = 0;
    ParmParse pp("amr");
    pp.query("max_level", nLevels);
    nLevels += 1;
    
    // get number of boxes over all levels
    int nBoxes = 0;
    for (int lev = 0; lev < nLevels; ++lev)
    {
        //const amrex::Vector<amrex::Real>* cost_heuristic = warpx.getCostsHeuristic(lev);
        nBoxes += 0; //??
    }

    // resize data array (for each box on each level we will save 6 data fields:
    // [cost, proc, lev, i, j, k])
    int nDataFields = 6;
    amrex::Print() << "nboxes: " << nBoxes << "\n";
    amrex::Print() << "nlevels: " << nLevels << "\n";
    m_data.resize(nDataFields*nBoxes,0.0);

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
            for (int boxNumber = 0; boxNumber < nBoxes; ++boxNumber)
            {
                ofs << m_sep;
                ofs << "[" + std::to_string(3 + nDataFields*boxNumber) + "]";
                ofs << "cost_box"+std::to_string(boxNumber);
                ofs << m_sep;
                ofs << "[" + std::to_string(4 + nDataFields*boxNumber) + "]";
                ofs << "proc_box"+std::to_string(boxNumber);
                ofs << m_sep;
                ofs << "[" + std::to_string(5 + nDataFields*boxNumber) + "]";
                ofs << "lev_box"+std::to_string(boxNumber);
                ofs << m_sep;
                ofs << "[" + std::to_string(6 + nDataFields*boxNumber) + "]";
                ofs << "i_box"+std::to_string(boxNumber);
                ofs << m_sep;
                ofs << "[" + std::to_string(7 + nDataFields*boxNumber) + "]";
                ofs << "j_box"+std::to_string(boxNumber);
                ofs << m_sep;
                ofs << "[" + std::to_string(8 + nDataFields*boxNumber) + "]";
                ofs << "k_box"+std::to_string(boxNumber);
            }
            ofs << std::endl;
            
            // close file
            ofs.close();
        }
    }
}
// end constructor

// function that computes costs
void Costs::ComputeDiags (int step)
{
    
    // Judge if the diags should be done
    if ( ((step+1) % m_freq != 0) ) { return; }

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // compute the costs if it's not a load balance step
    if ((step+1) % warpx.get_load_balance_int() != 0)
    {
        if (warpx.load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Heuristic)
        {
            warpx.ComputeCostsHeuristic();

        } else if (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            // costs are update via timers, nothing to do (for now)
        }
    }
    
    // keeps track of correct index in array over all boxes on all levels
    int shift = 0;

    // save data
    auto nLevels = warpx.finestLevel() + 1;
    for (int lev = 0; lev < nLevels; ++lev)
    {
        const amrex::DistributionMapping& dm = warpx.DistributionMap(lev);
        const amrex::Vector<amrex::Real>* cost_heuristic = warpx.getCostsHeuristic(lev);
        const MultiFab & Ex = warpx.getEfield(lev,0);
        for (MFIter mfi(Ex, false); mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            //m_data[mfi.index()] = 0.;
            // m_data[shift + mfi.index()*nDataFields + 0] = (*cost_heuristic)[mfi.index()];
            // m_data[shift + mfi.index()*nDataFields + 1] = dm[mfi.index()];
            // m_data[shift + mfi.index()*nDataFields + 2] = lev;
            // m_data[shift + mfi.index()*nDataFields + 3] = tbx.loVect()[0];
            // m_data[shift + mfi.index()*nDataFields + 4] = tbx.loVect()[1];
            // m_data[shift + mfi.index()*nDataFields + 5] = tbx.loVect()[2];                
        }
        
        // we looped through all the boxes on level lev, update the shift index
        shift += nDataFields*(cost_heuristic->size());
    }
        
    // parallel reduce to IO proc and get data over all procs
    if (ParallelDescriptor::IOProcessor())
    {
        //amrex::Vector<Real>::iterator it_m_data = m_data.begin();
        //amrex::Real* addr_it_m_data = &(*it_m_data);
        //ParallelAllReduce::Sum(addr_it_m_data, m_data.size(), ParallelContext::CommunicatorSub());
    }
    /* m_data now contains up-to-date values for:
     *  [cost, proc, lev, i, j, k] of box 0 at level 0,
     *   cost, proc, lev, i, j, k] of box 1 at level 0,
     *   cost, proc, lev, i, j, k] of box 2 at level 0,
     *  ...
     *   cost, proc, lev, i, j, k] of box 0 at level 1,
     *   cost, proc, lev, i, j, k] of box 1 at level 1,
     *   cost, proc, lev, i, j, k] of box 2 at level 1,
     *   ......] */

}
// end void Costs::ComputeDiags
