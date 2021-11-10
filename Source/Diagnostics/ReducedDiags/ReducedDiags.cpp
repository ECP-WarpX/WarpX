/* Copyright 2019-2020 Maxence Thevenet, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "ReducedDiags.H"

#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>

#include <fstream>
#include <iomanip>

using namespace amrex;

// constructor
ReducedDiags::ReducedDiags (std::string rd_name)
{
    m_rd_name = rd_name;

    BackwardCompatibility();

    ParmParse pp_rd_name(m_rd_name);

    // read path
    pp_rd_name.query("path", m_path);

    // read extension
    pp_rd_name.query("extension", m_extension);

    // check if it is a restart run
    std::string restart_chkfile = "";
    ParmParse pp_amr("amr");
    pp_amr.query("restart", restart_chkfile);
    m_IsNotRestart = restart_chkfile.empty();

    if (ParallelDescriptor::IOProcessor())
    {
        // create folder
        if (!UtilCreateDirectory(m_path, 0755))
        { CreateDirectoryFailed(m_path); }

        // replace / create output file
        if ( m_IsNotRestart ) // not a restart
        {
            std::ofstream ofs{m_path+m_rd_name+"."+m_extension, std::ios::trunc};
            ofs.close();
        }
    }

    // read reduced diags intervals
    std::vector<std::string> intervals_string_vec = {"1"};
    pp_rd_name.getarr("intervals", intervals_string_vec);
    m_intervals = IntervalsParser(intervals_string_vec);

    // read separator
    pp_rd_name.query("separator", m_sep);
}
// end constructor

void ReducedDiags::InitData ()
{
    // Defines an empty function InitData() to be overwritten if needed.
    // Function used to initialize data of the diagnostics after the WarpX
    // data structures are all set up
}

void ReducedDiags::LoadBalance ()
{
    // Defines an empty function LoadBalance() to be overwritten if needed.
    // Function used to redistribute parallel data of the diagnostics in
    // load balancing operations
}

void ReducedDiags::BackwardCompatibility ()
{
    amrex::ParmParse pp_rd_name(m_rd_name);
    std::vector<std::string> backward_strings;
    if (pp_rd_name.queryarr("frequency", backward_strings)){
        amrex::Abort("<reduced_diag_name>.frequency is no longer a valid option. "
                     "Please use the renamed option <reduced_diag_name>.intervals instead.");
    }
}

// write to file function
void ReducedDiags::WriteToFile (int step) const
{
    // open file
    std::ofstream ofs{m_path + m_rd_name + "." + m_extension,
        std::ofstream::out | std::ofstream::app};

    // write step
    ofs << step+1;

    ofs << m_sep;

    // set precision
    ofs << std::fixed << std::setprecision(14) << std::scientific;

    // write time
    ofs << WarpX::GetInstance().gett_new(0);

    // loop over data size and write
    for (int i = 0; i < static_cast<int>(m_data.size()); ++i)
    {
        ofs << m_sep;
        ofs << m_data[i];
    }
    // end loop over data size

    // end line
    ofs << std::endl;

    // close file
    ofs.close();
}
// end ReducedDiags::WriteToFile
