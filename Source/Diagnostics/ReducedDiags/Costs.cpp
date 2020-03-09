/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Costs.H"
#include "WarpX.H"
//#include "Utils/WarpXConst.H"

#include <AMReX_REAL.H>
//#include <AMReX_ParticleReduce.H>

#include <iostream>
//#include <cmath>


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

    // Does it make sense to get nBoxes at t=0?  What if it changes?
    // resize data array (for each box on each level we will save 6 data fields:
    // [cost, proc, lev, i, j, k])

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
            for (int boxNumber = 0; boxNumber < 1; ++boxNumber)
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

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();
    
    // Judge if the diags should be done
    // costs is initialized only if we're doing load balance
    if ( ((step+1) % m_freq != 0)
         || step <= 0
         || warpx.get_load_balance_int() < 1 ) { return; }

    // get number of boxes over all levels
    auto nLevels = warpx.finestLevel() + 1;
    int nBoxes = 0;
    for (int lev = 0; lev < nLevels; ++lev)
    {
        const amrex::Vector<amrex::Real>* cost_heuristic = warpx.getCostsHeuristic(lev);
        const MultiFab & Ex = warpx.getEfield(lev,0);
        // we need to compute the number of boxes
        nBoxes += Ex.size();
        //amrex::Print() << "nBoxes: " << nBoxes << "step: " << step << "\n";
    }
    //amrex::AllPrint() << "proc: " << ParallelDescriptor::MyProc() << " nDataFields*nBoxes" << nDataFields*nBoxes << "\n";    
    m_data.resize(nDataFields*nBoxes,0.0);
    // you just got to clear it
    m_data.assign(nDataFields*nBoxes,0.0);
    //amrex::AllPrint() << "proc: " << ParallelDescriptor::MyProc() << " m_data.size()" << m_data.size() << "\n";
    
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

    // amrex::Print()  << " m_data before";
    // for (int iter=0; iter<m_data.size(); ++iter) {
    //     amrex::Print() << " " << m_data[iter];
    // }
    // amrex::Print()  << "\n";
    
    // save data
    for (int lev = 0; lev < nLevels; ++lev)
    {
        const amrex::DistributionMapping& dm = warpx.DistributionMap(lev);
        const amrex::Vector<amrex::Real>* cost_heuristic = warpx.getCostsHeuristic(lev);
        const MultiFab & Ex = warpx.getEfield(lev,0);
        for (MFIter mfi(Ex, false); mfi.isValid(); ++mfi)
        {
            const Box& tbx = mfi.tilebox();
            m_data[shift + mfi.index()*nDataFields + 0] = (*cost_heuristic)[mfi.index()];
            m_data[shift + mfi.index()*nDataFields + 1] = dm[mfi.index()];
            m_data[shift + mfi.index()*nDataFields + 2] = lev;
            m_data[shift + mfi.index()*nDataFields + 3] = tbx.loVect()[0];
            m_data[shift + mfi.index()*nDataFields + 4] = tbx.loVect()[1];
            m_data[shift + mfi.index()*nDataFields + 5] = tbx.loVect()[2];
        }
        
        // we looped through all the boxes on level lev, update the shift index
        shift += nDataFields*(Ex.size());
    }

    // parallel reduce to IO proc and get data over all procs
    //if (ParallelDescriptor::IOProcessor())
    //{
    // amrex::Print()  << " m_data after";
    // for (int iter=0; iter<m_data.size(); ++iter) {
    //     amrex::Print() << " " << m_data[iter];
    // }
    // amrex::Print()  << "\n";
        
    amrex::Vector<amrex::Real>::iterator it_m_data = m_data.begin();
    amrex::Real* addr_it_m_data = &(*it_m_data);
    // Be sure, this does what you want
    //if (ParallelDescriptor::IOProcessor())
    //{
    // How I can do this only for IOproc?
    ParallelAllReduce::Sum(addr_it_m_data,
                           m_data.size(),
                           ParallelContext::CommunicatorSub());
        //}
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
