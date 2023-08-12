/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Yinjian Zhao, Arianna Formenti
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_DIAGNOSTICS_REDUCEDDIAGS_COLLIDERRELEVANT_H_
#define WARPX_DIAGNOSTICS_REDUCEDDIAGS_COLLIDERRELEVANT_H_

#include "ReducedDiags.H"

#include <map>
#include <string>
#include <vector>

/**
 *  This class mainly contains a function that
 *  computes the particle extrema of each species.
 */
class ColliderRelevant : public ReducedDiags
{
public:

    /**
     * constructor
     * @param[in] rd_name reduced diags names
     */
    ColliderRelevant(std::string rd_name);

    /// name of species
    std::vector<std::string> m_beam_name;

    /**
     * This funciton computes the particle extrema
     *
     * @param[in] step current time step
     */
    void ComputeDiags(int step) override final;

private:
    /// auxiliary structure to store headers and indices of the reduced diagnostics
    struct aux_header_index
    {
        std::string header;
        int idx;
    };

    /// map to store header texts and indices of the reduced diagnostics
    std::map<std::string, aux_header_index> m_headers_indices;
};

#endif
