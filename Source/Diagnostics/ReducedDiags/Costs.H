/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_COSTS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_COSTS_H_

#include "WarpX.H"
#include "ReducedDiags.H"

/**
 *  This class mainly contains a function that
 *  computes the costs used for loadbalance.
 */
class Costs : public ReducedDiags
{
public:

    /** number of data fields we save for each box 
     *  (cost, processor, level, i_low, j_low, k_low) */
    const int nDataFields = 6;

    /** used to keep track of max number of boxes over all timesteps; this allows
     *  to compute the number of NaNs required to fill jagged array into a 
     *  rectangular one */
    int nBoxesMax = -1;
    
    /** constructor
     *  @param[in] rd_name reduced diags names */
    Costs(std::string rd_name);

    /** This funciton computes the costs.
     *  Heuristic costs are updated at each load balance step.  If costs reduced 
     *  diagnostic is desired at a non-load balance step, the heuristic cost is
     *  updated.  Timer-based costs are recorded inline, so in this case there 
     *  is nothing to compute/update.  */
    virtual void ComputeDiags(int step) override final;

    /** write to file function for costs;  this differs from the base class 
     *  `ReducedDiags` in that it will fill in blank entries with NaN at the 
     *  final timestep, ensuring that the data array is not jagged
     *  @param[in] step time step */
    virtual void WriteToFile(int step) const override final;
    
};

#endif
