/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "Costs.H"
#include "WarpX.H"


using namespace amrex;

// constructor
Costs::Costs (std::string rd_name)
: ReducedDiags{rd_name}
{

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs;
            ofs.open(m_path + m_rd_name + "." + m_extension,
                     std::ofstream::out | std::ofstream::app);

            // write header row
            // for each box on each level we will save 6 data fields: [cost, proc, lev, i_low, j_low, k_low])
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";
            ofs << m_sep;
            for (int boxNumber=0; boxNumber<1; boxNumber++)
            {
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
                ofs << std::endl;
            }

            // close file
            ofs.close();
        }
    }
}
// end constructor

// function that gathers costs
void Costs::ComputeDiags (int step)
{

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();
    
    // judge if the diags should be done
    // costs is initialized only if we're doing load balance
    if ( ((step+1) % m_freq != 0) || warpx.get_load_balance_int() < 1 ) { return; }

    // get number of boxes over all levels
    auto nLevels = warpx.finestLevel() + 1;
    int nBoxes = 0;
    for (int lev = 0; lev < nLevels; ++lev)
    {
        const amrex::Vector<amrex::Real>* cost_heuristic = warpx.getCostsHeuristic(lev);
        nBoxes += cost_heuristic->size();
    }
    
    // resize and clear data array
    m_data.resize(nDataFields*nBoxes, 0.0);
    m_data.assign(nDataFields*nBoxes, 0.0);
    
    // if not a load balance step, we must compute costs
    if ((step+1) % warpx.get_load_balance_int() != 0)
    {
        if (warpx.load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Heuristic)
        {
            warpx.ComputeCostsHeuristic();

        } else if (WarpX::load_balance_costs_update_algo == LoadBalanceCostsUpdateAlgo::Timers)
        {
            // costs are update via timers, nothing to do
        }
    }
    
    // keeps track of correct index in array over all boxes on all levels
    int shift = 0;

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
        shift += nDataFields*(cost_heuristic->size());
    }

    // parallel reduce to IO proc and get data over all procs
    ParallelDescriptor::ReduceRealSum(m_data, ParallelDescriptor::IOProcessorNumber());
    
    //amrex::Vector<amrex::Real>::iterator it_m_data = m_data.begin();
    //amrex::Real* addr_it_m_data = &(*it_m_data);
    // ParallelAllReduce::Sum(addr_it_m_data,
    //                        m_data.size(),
    //                        ParallelContext::CommunicatorSub());

    /* m_data now contains up-to-date values for:
     *  [[cost, proc, lev, i, j, k] of box 0 at level 0,
     *   [cost, proc, lev, i, j, k] of box 1 at level 0,
     *   [cost, proc, lev, i, j, k] of box 2 at level 0,
     *   ...
     *   [cost, proc, lev, i, j, k] of box 0 at level 1,
     *   [cost, proc, lev, i, j, k] of box 1 at level 1,
     *   [cost, proc, lev, i, j, k] of box 2 at level 1,
     *   ......] */

}
// end void Costs::ComputeDiags
