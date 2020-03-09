/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_COSTS_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_COSTS_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the costs used for loadbalance.
 */
class Costs : public ReducedDiags
{
public:

    static const int nDataFields = 6;
    
    /** constructor
     *  @param[in] rd_name reduced diags names */
    Costs(std::string rd_name);

    /** This funciton computes the costs.
     *  Description of the computation.  */
    virtual void ComputeDiags(int step) override final;

};

#endif
