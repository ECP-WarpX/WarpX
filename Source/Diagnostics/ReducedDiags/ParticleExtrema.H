/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEEXTREMA_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEEXTREMA_H_

#include "ReducedDiags.H"
#include <fstream>

/**
 *  This class mainly contains a function that
 *  computes the particle extrema of each species.
 */
class ParticleExtrema : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ParticleExtrema(std::string rd_name);

    /// name of species
    std::string m_species_name;

    /** This funciton computes the particle extrema
     */
    void ComputeDiags(int step) override final;

};

#endif
