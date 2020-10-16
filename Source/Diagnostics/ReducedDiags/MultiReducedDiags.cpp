/* Copyright 2019-2020 Maxence Thevenet, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "LoadBalanceCosts.H"
#include "ParticleHistogram.H"
#include "BeamRelevant.H"
#include "ParticleEnergy.H"
#include "FieldEnergy.H"
#include "FieldMaximum.H"
#include "ParticleNumber.H"
#include "MultiReducedDiags.H"

#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include <fstream>

using namespace amrex;

// constructor
MultiReducedDiags::MultiReducedDiags ()
{

    // read reduced diags names
    ParmParse pp("warpx");
    m_plot_rd = pp.queryarr("reduced_diags_names", m_rd_names);

    // if names are not given, reduced diags will not be done
    if ( m_plot_rd == 0 ) { return; }

    // resize
    m_multi_rd.resize(m_rd_names.size());

    // loop over all reduced diags
    for (int i_rd = 0; i_rd < static_cast<int>(m_rd_names.size()); ++i_rd)
    {

        ParmParse pp_rd(m_rd_names[i_rd]);

        // read reduced diags type
        std::string rd_type;
        pp_rd.query("type", rd_type);

        // match diags
        if (rd_type.compare("ParticleEnergy") == 0)
        {
            m_multi_rd[i_rd] =
                std::make_unique<ParticleEnergy>(m_rd_names[i_rd]);
        }
        else if (rd_type.compare("FieldEnergy") == 0)
        {
            m_multi_rd[i_rd] =
                std::make_unique<FieldEnergy>(m_rd_names[i_rd]);
        }
        else if (rd_type.compare("FieldMaximum") == 0)
        {
            m_multi_rd[i_rd] =
                std::make_unique<FieldMaximum>(m_rd_names[i_rd]);
        }
        else if (rd_type.compare("BeamRelevant") == 0)
        {
            m_multi_rd[i_rd] =
                std::make_unique<BeamRelevant>(m_rd_names[i_rd]);
        }
        else if (rd_type.compare("LoadBalanceCosts") == 0)
        {
            m_multi_rd[i_rd] =
                std::make_unique<LoadBalanceCosts>(m_rd_names[i_rd]);
        }
        else if (rd_type.compare("ParticleHistogram") == 0)
        {
            m_multi_rd[i_rd] =
                std::make_unique<ParticleHistogram>(m_rd_names[i_rd]);
        }
        else if (rd_type.compare("ParticleNumber") == 0)
        {
            m_multi_rd[i_rd]=
                std::make_unique<ParticleNumber>(m_rd_names[i_rd]);
        }
        else
        { Abort("No matching reduced diagnostics type found."); }
        // end if match diags

    }
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
