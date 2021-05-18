/* Copyright 2019-2020 Maxence Thevenet, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "LoadBalanceCosts.H"
#include "LoadBalanceEfficiency.H"
#include "ParticleHistogram.H"
#include "BeamRelevant.H"
#include "ParticleEnergy.H"
#include "ParticleExtrema.H"
#include "FieldEnergy.H"
#include "FieldMaximum.H"
#include "RhoMaximum.H"
#include "ParticleNumber.H"
#include "MultiReducedDiags.H"

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include <fstream>
#include <map>
#include <algorithm>
#include <functional>
#include <iterator>

using namespace amrex;

// constructor
MultiReducedDiags::MultiReducedDiags ()
{

    // read reduced diags names
    ParmParse pp_warpx("warpx");
    m_plot_rd = pp_warpx.queryarr("reduced_diags_names", m_rd_names);

    // if names are not given, reduced diags will not be done
    if ( m_plot_rd == 0 ) { return; }

    using CS = const std::string& ;
    const auto reduced_diags_dictionary =
        std::map<std::string, std::function<std::unique_ptr<ReducedDiags>(CS)>>{
            {"ParticleEnergy",        [](CS s){return std::make_unique<ParticleEnergy>(s);}},
            {"FieldEnergy",           [](CS s){return std::make_unique<FieldEnergy>(s);}},
            {"FieldMaximum",          [](CS s){return std::make_unique<FieldMaximum>(s);}},
            {"RhoMaximum",            [](CS s){return std::make_unique<RhoMaximum>(s);}},
            {"BeamRelevant",          [](CS s){return std::make_unique<BeamRelevant>(s);}},
            {"LoadBalanceCosts",      [](CS s){return std::make_unique<LoadBalanceCosts>(s);}},
            {"LoadBalanceEfficiency", [](CS s){return std::make_unique<LoadBalanceEfficiency>(s);}},
            {"ParticleHistogram",     [](CS s){return std::make_unique<ParticleHistogram>(s);}},
            {"ParticleNumber",        [](CS s){return std::make_unique<ParticleNumber>(s);}},
            {"ParticleExtrema",       [](CS s){return std::make_unique<ParticleExtrema>(s);}}
        };
    // loop over all reduced diags and fill m_multi_rd with requested reduced diags
    std::transform(m_rd_names.begin(), m_rd_names.end(), std::back_inserter(m_multi_rd),
        [&](const auto& rd_name){
            ParmParse pp_rd_name(rd_name);

            // read reduced diags type
            std::string rd_type;
            pp_rd_name.query("type", rd_type);

            if(reduced_diags_dictionary.count(rd_type) == 0)
                Abort("No matching reduced diagnostics type found.");

            return reduced_diags_dictionary.at(rd_type)(rd_name);
        });
    // end loop over all reduced diags

}
// end constructor

// call functions to compute diags
void MultiReducedDiags::ComputeDiags (int step)
{
    // loop over all reduced diags
    for (int i_rd = 0; i_rd < static_cast<int>(m_rd_names.size()); ++i_rd)
    {
        m_multi_rd[i_rd] -> ComputeDiags(step);
    }
    // end loop over all reduced diags
}
// end void MultiReducedDiags::ComputeDiags

// function to write data
void MultiReducedDiags::WriteToFile (int step)
{

    // Only the I/O rank does
    if ( !ParallelDescriptor::IOProcessor() ) { return; }

    // loop over all reduced diags
    for (int i_rd = 0; i_rd < static_cast<int>(m_rd_names.size()); ++i_rd)
    {
        // Judge if the diags should be done
        if (!m_multi_rd[i_rd]->m_intervals.contains(step+1)) { continue; }

        // call the write to file function
        m_multi_rd[i_rd]->WriteToFile(step);

    }
    // end loop over all reduced diags
}
// end void MultiReducedDiags::WriteToFile
