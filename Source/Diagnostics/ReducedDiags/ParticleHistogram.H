/* Copyright 2019-2020 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEHISTOGRAM_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_PARTICLEHISTOGRAM_H_

#include "ReducedDiags.H"
#include <WarpXParserWrapper.H>
#include <fstream>

/**
 *  This class contains a function that
 *  computes a histogram of user defined quantity.
 */
class ParticleHistogram : public ReducedDiags
{
public:

    /** constructor
     *  @param[in] rd_name reduced diags names */
    ParticleHistogram(std::string rd_name);

    /// name of species
    std::string m_species_name;

    /// particle weight type
    std::string m_norm = "default";

    /// number of bins
    int m_bin_num;

    /// max and min bin values
    amrex::Real m_bin_max;
    amrex::Real m_bin_min;

    /// bin size
    amrex::Real m_bin_size;

    /// parser
    std::unique_ptr<ParserWrapper<7>> m_parser;

    /** This funciton computes a histogram of user defined quantity.
     *  \param [in] step current time step.
     */
    virtual void ComputeDiags(int step) override final;

};

#endif
