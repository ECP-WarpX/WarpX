#include "FlushFormat.H"

#include "WarpX.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_buildInfo.H>


void
FlushFormat::WriteJobInfo(const std::string& dir) const
{

    auto & warpx = WarpX::GetInstance();

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        // job_info file with details about the run
        std::ofstream jobInfoFile;
        std::string FullPathJobInfoFile = dir;

        const std::string PrettyLine = std::string(78, '=') + "\n";
//        std::string OtherLine = std::string(78, '-') + "\n";
//        std::string SkipSpace = std::string(8, ' ') + "\n";

        FullPathJobInfoFile += "/warpx_job_info";
        jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

        // job information
        jobInfoFile << PrettyLine;
        jobInfoFile << " WarpX Job Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "number of MPI processes: " << amrex::ParallelDescriptor::NProcs() << "\n";
#ifdef AMREX_USE_OMP
        jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

        jobInfoFile << "\n\n";

        // build information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Build Information\n";
        jobInfoFile << PrettyLine;

        jobInfoFile << "build date:    " << amrex::buildInfoGetBuildDate() << "\n";
        jobInfoFile << "build machine: " << amrex::buildInfoGetBuildMachine() << "\n";
        jobInfoFile << "build dir:     " << amrex::buildInfoGetBuildDir() << "\n";
        jobInfoFile << "AMReX dir:     " << amrex::buildInfoGetAMReXDir() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "COMP:          " << amrex::buildInfoGetComp() << "\n";
        jobInfoFile << "COMP version:  " << amrex::buildInfoGetCompVersion() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "C++ compiler:  " << amrex::buildInfoGetCXXName() << "\n";
        jobInfoFile << "C++ flags:     " << amrex::buildInfoGetCXXFlags() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "Fortran comp:  " << amrex::buildInfoGetFName() << "\n";
        jobInfoFile << "Fortran flags: " << amrex::buildInfoGetFFlags() << "\n";

        jobInfoFile << "\n";

        jobInfoFile << "Link flags:    " << amrex::buildInfoGetLinkFlags() << "\n";
        jobInfoFile << "Libraries:     " << amrex::buildInfoGetLibraries() << "\n";

        jobInfoFile << "\n";

        const char* githash1 = amrex::buildInfoGetGitHash(1);
        const char* githash2 = amrex::buildInfoGetGitHash(2);
        const char* githash3 = amrex::buildInfoGetGitHash(3);
        if (strlen(githash1) > 0) {
            jobInfoFile << "WarpX  git describe: " << githash1 << "\n";
        }
        if (strlen(githash2) > 0) {
            jobInfoFile << "AMReX  git describe: " << githash2 << "\n";
        }
        if (strlen(githash3) > 0) {
            jobInfoFile << "PICSAR git describe: " << githash3 << "\n";
        }

        jobInfoFile << "\n\n";

        // grid information
        jobInfoFile << PrettyLine;
        jobInfoFile << " Grid Information\n";
        jobInfoFile << PrettyLine;

        for (int i = 0; i <= warpx.finestLevel(); i++)
        {
            jobInfoFile << " level: " << i << "\n";
            jobInfoFile << "   number of boxes = " << warpx.boxArray(i).size() << "\n";
            jobInfoFile << "   maximum zones   = ";
            for (int n = 0; n < AMREX_SPACEDIM; n++)
            {
                jobInfoFile << warpx.Geom(i).Domain().length(n) << " ";
            }
            jobInfoFile << "\n\n";
        }

        jobInfoFile << " Boundary conditions\n";

        jobInfoFile << "   -x: " << "interior" << "\n";
        jobInfoFile << "   +x: " << "interior" << "\n";
        if (AMREX_SPACEDIM >= 2) {
            jobInfoFile << "   -y: " << "interior" << "\n";
            jobInfoFile << "   +y: " << "interior" << "\n";
        }
#if defined(WARPX_DIM_3D)
            jobInfoFile << "   -z: " << "interior" << "\n";
            jobInfoFile << "   +z: " << "interior" << "\n";
#endif

        jobInfoFile << "\n\n";


        // runtime parameters
        jobInfoFile << PrettyLine;
        jobInfoFile << " Inputs File Parameters\n";
        jobInfoFile << PrettyLine;

        amrex::ParmParse::dumpTable(jobInfoFile, true);

        jobInfoFile.close();
    }
}
