/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLENUMBER_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLENUMBER_H_

#include "ReducedDiags.H"

/**
 *  This class mainly contains a function that computes the total number of macroparticles of each
 *  species.
 */
class ParticleNumber : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ParticleNumber(std::string rd_name);

    /** This function computes the total number of macroparticles of each species.
     *  @param [in] step current time step
     */
    virtual void ComputeDiags(int step) override final;

};

#endif // WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLENUMBER_H_
