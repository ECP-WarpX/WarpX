/* Copyright 2019-2020 Axel Huebl, Ligia Diana Amorim, Maxence Thevenet
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "InjectorDensity.H"

#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/Algorithms/LinearInterpolation.H"

#include <AMReX_ParmParse.H>
#include <AMReX.H>

#ifdef AMREX_USE_SENSEI_INSITU
#   include <AMReX_AmrMeshInSituBridge.H>
#endif
#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

#include <algorithm>
#include <cctype>
#include <vector>

#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
namespace io = openPMD;
#endif

using namespace amrex;
namespace
{
}


void InjectorDensity::clear ()
{
    switch (type)
    {
    case Type::parser:
    {
        break;
    }
    case Type::predefined:
    {
        object.predefined.clear();
        break;
    }
    case Type::fromfile:
    {
        object.fromfile.clear();
        break;
    }
    default:
        return;
    }
}

InjectorDensityPredefined::InjectorDensityPredefined (
    std::string const& a_species_name) noexcept
{
    const ParmParse pp_species_name(a_species_name);

    std::vector<amrex::Real> v;
    // Read parameters for the predefined plasma profile.
    utils::parser::getArrWithParser(
        pp_species_name, "predefined_profile_params", v);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v.size() <= 6,
                                     "Too many parameters for InjectorDensityPredefined");
    for (int i = 0; i < static_cast<int>(v.size()); ++i) {
        p[i] = v[i];
    }

    // Parse predefined profile name, and update member variable profile.
    std::string which_profile_s;
    pp_species_name.query("predefined_profile_name", which_profile_s);
    std::transform(which_profile_s.begin(), which_profile_s.end(),
                   which_profile_s.begin(), ::tolower);
    if (which_profile_s == "parabolic_channel"){
        profile = Profile::parabolic_channel;
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(v.size() > 5,
            "InjectorDensityPredefined::parabolic_channel: not enough parameters");
    }
}

InjectorDensityFromFile::InjectorDensityFromFile (std::string const & a_species_name)
{

    //get the file path for reading the external density and the name of the field to be read
    const ParmParse pp (a_species_name);
    const ParmParse pp_warpx ("warpx");
    const ParmParse pp_amr ("amr");

    std::string read_density_from_path= "./";
    pp.query("read_density_from_path", external_density_path);
    pp.query("density_name", density);
    pp_amr.query("max_grid_size", max_grid_size);

    std::vector<int> cells(3);
    utils::parser::getArrWithParser(pp_amr, "n_cell", cells);

    for (int i = 0; i < cells.size(); i++) {
        std::cout << i;
        n_cell.push_back(cells[i]);
    }

    std::cout << "[1] File path is: " << external_density_path << '\n' ;
    std::cout << "[2] Density name is: " << density << '\n' ;
    std::cout << "[3] Max grid size is: " << max_grid_size << '\n' ;
    std::cout << "[4] N_cell for each dimension: " << n_cell[0] << " " << n_cell[1] << " " << n_cell[2] << '\n';

    // Get WarpX domain info
    WarpX& warpx = WarpX::GetInstance();
    amrex::Geometry const& geom0 = warpx.Geom(0);
    real_box = geom0.ProbDomain();
    auto const dx = geom0.CellSizeArray();
    auto geom_lo = geom0.ProbLoArray();
    auto geom_hi = geom0.ProbHiArray();
    lo0 = geom_lo[0];
    lo1 = geom_lo[1];
    lo2 = geom_lo[2];
    hi0 = geom_hi[0];
    hi1 = geom_hi[1];
    hi2 = geom_hi[2];

    // creating the mulitfab array
    IntVect rho_nodal_flag;
    rho_nodal_flag = IntVect( AMREX_D_DECL(1,1,1) );
    MultiFab *mf = warpx.get_pointer_rho_fp(0);
    const amrex::IntVect nodal_flag = mf->ixType().toIntVect();

    std::cout << rho_nodal_flag << '\n';

    // Read external field openPMD data
    io::Series series = io::Series(external_density_path, io::Access::READ_ONLY);
    auto iseries = series.iterations.begin()->second;

    io::Mesh P = iseries.meshes[density];
    io::MeshRecordComponent p_scalar0 = P[io::RecordComponent::SCALAR];
    //auto all_data0 = p_scalar0.loadChunk<double>();

    auto axisLabels = P.getAttribute("axisLabels").get<std::vector<std::string>>();
    auto fileGeom = P.getAttribute("geometry").get<std::string>();

    io::ParticleSpecies particle = iseries.particles[a_species_name];
    charge = particle["charge"][openPMD::RecordComponent::SCALAR].loadChunk<double>();
    std::cout << "[5] electron charge is " << charge << '\n';

    offset = P.gridGlobalOffset();
    spacing = P.gridSpacing< long double >();
    extent = p_scalar0.getExtent();
    const int extent0 = static_cast<int>(extent[0]);
    const int extent1 = static_cast<int>(extent[1]);
    const int extent2 = static_cast<int>(extent[2]);

    std::cout << "[7] Extent for each dimension: " << extent[0] << " " << extent[1] << " " << extent[2] << '\n';

    const long double file_dx = static_cast<amrex::Real>(spacing[0]);
    const long double file_dy = static_cast<amrex::Real>(spacing[1]);
    const long double file_dz = static_cast<amrex::Real>(spacing[2]);

    const io::Offset chunk_offset =  {0, 0, 0};  //{0,0,0};
    const io::Extent chunk_extent = {extent[0], extent[1], extent[2]};

    std::shared_ptr<double> P_chunk_data = p_scalar0.loadChunk<double>(chunk_offset,chunk_extent);
    series.flush();
    P_data_host = P_chunk_data.get();

    // Load data to GPU
    const size_t total_extent = size_t(extent[0]) * extent[1] * extent[2];
    amrex::Gpu::DeviceVector<double> P_data_gpu(total_extent);
    double* P_data = P_data_gpu.data();
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, P_data_host, P_data_host + total_extent, P_data);

    std::cout << "Offset is : (" << offset[0] << ", " << offset[1] <<  ", " << offset[2] << ")\n";

    // iterate over the external field and find the index that corresponds to the grid points
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        std::cout << "inside MFIter\n";
        const amrex::Box box = mfi.growntilebox();
        const amrex::Box tb = mfi.tilebox(rho_nodal_flag, mf->nGrowVect());
        std::cout << "before mffab definition\n";
        mffab = mf->array(mfi);
        std::cout << "after mffab definition\n";

        amrex::ParallelFor (tb,
            [box, this, dx,file_dx, file_dy, file_dz, extent0, extent1, extent2, P_data] AMREX_GPU_DEVICE (int i, int j, int k) {
                // i,j,k denote x,y,z indices in 3D xyz.
                // i,j denote r,z indices in 2D rz; k is just 0

                //for 3D only!!
                const int ii = i;
                amrex::Real x0, x1;
                std::cout << "inside parallel for";

                if ( box.type(0)==amrex::IndexType::CellIndex::NODE )
                    {x0 = static_cast<amrex::Real>(real_box.lo(0)) + ii*dx[0]; }
                else {
                    x0 = static_cast<amrex::Real>(real_box.lo(0)) + ii*dx[0] + 0.5_rt*dx[0]; }
                if (box.type(1)==amrex::IndexType::CellIndex::NODE )
                    {x1 = real_box.lo(1) + j*dx[1]; }
                else {
                    x1 = real_box.lo(1) + j*dx[1] + 0.5_rt*dx[1]; }

                amrex::Real x2;
                if ( box.type(2)==amrex::IndexType::CellIndex::NODE )
                    { x2 = real_box.lo(2) + k*dx[2]; }
                else { x2 = real_box.lo(2) + k*dx[2] + 0.5_rt*dx[2]; }

                // Get index of the external field array

                int const ix = std::floor( (x0-offset[0])/file_dx );
                int const iy = std::floor( (x1-offset[1])/file_dy );
                int const iz = std::floor( (x2-offset[2])/file_dz );

                // Get coordinates of external grid point
                auto const xx0 = offset[0] + ix * file_dx;
                auto const xx1 = offset[1] + iy * file_dy;
                auto const xx2 = offset[2] + iz * file_dz;

                const amrex::Array4<double> fc_array(P_data, {0,0,0}, {extent2, extent1, extent0}, 1);
                const double
                    f000 = fc_array(iz  , iy  , ix  ),
                    f001 = fc_array(iz+1, iy  , ix  ),
                    f010 = fc_array(iz  , iy+1, ix  ),
                    f011 = fc_array(iz+1, iy+1, ix  ),
                    f100 = fc_array(iz  , iy  , ix+1),
                    f101 = fc_array(iz+1, iy  , ix+1),
                    f110 = fc_array(iz  , iy+1, ix+1),
                    f111 = fc_array(iz+1, iy+1, ix+1);
                mffab(i,j,k) = static_cast<amrex::Real> (utils::algorithms::trilinear_interp<double>
                    (xx0, xx0+file_dx, xx1, xx1+file_dy, xx2, xx2+file_dz,
                     f000, f001, f010, f011, f100, f101, f110, f111,
                     x0, x1, x2) );
                //mffab(i,j,k) = static_cast<amrex::Real> (f000);
//#endif

        }
        ); //end ParallelFor
    }
    std::cout << "made it to the end\n" ;

}

// Note that we are not allowed to have non-trivial destructor.
// So we rely on clear() to free memory if needed.
void InjectorDensityPredefined::clear ()
{
}
void InjectorDensityFromFile::clear ()
{
}
