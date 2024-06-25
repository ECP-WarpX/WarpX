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
//#include <AMReX_MultiFabUtill.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

//#include <algorithm>
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
    std::cout << "start injector density from file\n" ;
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

    for (int i = 0; i < static_cast<int>(cells.size()); i++) {
        std::cout << i;
        n_cell.push_back(cells[i]);
    }

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

    std::cout << "lo0, lo1, lo2: " << lo0 << " " << lo1 << " " << lo2 << '\n';

    // creating the mulitfab array
    std::cout << "calling multifab\n";
    m_rho = warpx.get_pointer_rho_fp(0);
    m_rho->setVal(0);
    const amrex::IntVect nodal_flag = m_rho->ixType().toIntVect();

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_rho != nullptr, "Inside Assert");

    // Read external field openPMD data
    io::Series series = io::Series(external_density_path, io::Access::READ_ONLY);
    auto iseries = series.iterations.begin()->second;

    io::Mesh P = iseries.meshes[density];
    io::MeshRecordComponent p_scalar0 = P[io::RecordComponent::SCALAR];

    auto axisLabels = P.getAttribute("axisLabels").get<std::vector<std::string>>();
    auto fileGeom = P.getAttribute("geometry").get<std::string>();

    offset = P.gridGlobalOffset();
    spacing = P.gridSpacing< long double >();
    extent = p_scalar0.getExtent();
    const int extent0 = static_cast<int>(extent[0]);
    const int extent1 = static_cast<int>(extent[1]);
    const int extent2 = static_cast<int>(extent[2]);

    std::cout << "Extent for each dimension: " << extent[0] << " " << extent[1] << " " << extent[2] << '\n';

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
    amrex::AllPrint() << "Before MFIter Loop...\n";
    for (MFIter mfi(*m_rho, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box tb = mfi.tilebox(nodal_flag, m_rho->nGrowVect());
        amrex::Array4<amrex::Real> rho_array4 = m_rho->array(mfi);

        int start[3] = {rho_array4.begin.x, rho_array4.begin.y, rho_array4.begin.z};
        int end[3] = {rho_array4.end.x, rho_array4.end.y, rho_array4.end.z};
        amrex::IntVect start_vect(start);
        amrex::IntVect end_vect(end);
        amrex::Box mf_bx(start_vect, end_vect);

        m_cell_size = geom0.CellSizeArray();

        amrex::AllPrint() << "MultiFab meta data\n";
        amrex::AllPrint() << "  tb lo * dx: " << tb.smallEnd(0) * dx[0] << "  " << tb.smallEnd(1) * dx[1] << "  " << tb.smallEnd(2) * dx[2] << '\n';
        amrex::AllPrint() << "  cell size: " << m_cell_size[0] << "  " << m_cell_size[1] << "  " << m_cell_size[2] << '\n';
        amrex::AllPrint() << "  dx: " << dx[0] << "  " << dx[1] << "  " << dx[2] << '\n';

        amrex::AllPrint() << "Source offsets\n";
        amrex::AllPrint() << "  offsets=" << offset[0] << " " << offset[1] << " " << offset[2] << "\n";
        amrex::AllPrint() << "  file_dx=" << file_dx << " " << file_dy << " " << file_dz << "\n";


        amrex::ParallelFor (tb,
                            [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
                                // i,j,k denote x,y,z indices in 3D xyz
                                // i,j denote r,z indices in 2D rz; k is just 0

                                // TODO: so far for 3D only!!
                                amrex::Real x0, x1, x2;

                                // real position in the multifab that we interpolate to
                                x0 = real_box.lo(0) + amrex::Real(i)*m_cell_size[0];
                                if ( tb.type(0)==amrex::IndexType::CellIndex::CELL ) {
                                    x0 += 0.5_rt*m_cell_size[0];
                                }
                                x1 = real_box.lo(1) + amrex::Real(j)*m_cell_size[1];
                                if (tb.type(1)==amrex::IndexType::CellIndex::CELL ) {
                                    x1 += 0.5_rt*m_cell_size[1];
                                }
                                x2 = real_box.lo(2) + amrex::Real(k)*m_cell_size[2];
                                if ( tb.type(2)==amrex::IndexType::CellIndex::CELL ){
                                    x2 += 0.5_rt*m_cell_size[2];
                                }

                                // the next lower index in the data source at the current real position
                                int const ix = std::floor( (x0-offset[0])/file_dx );
                                int const iy = std::floor( (x1-offset[1])/file_dy );
                                int const iz = std::floor( (x2-offset[2])/file_dz );

                                // Get coordinates of data source grid point (real position)
                                auto const xx0 = ix * file_dx + offset[0];
                                auto const xx1 = iy * file_dy + offset[1];
                                auto const xx2 = iz * file_dz + offset[2];

                                // source data access
                                const amrex::Array4<double> src_data(P_data, {0, 0, 0}, {extent2, extent1, extent0}, 1);

                                if (iz < 0 | iz > extent2-2 | iy < 0 | iy > extent1-2 | ix < 0 | ix > extent0-2 ){
                                    rho_array4(i,j,k) = static_cast<amrex::Real> (0);
                                }else{
                                    // interpolate from closest source data points to rho MultiFab resolution
                                    const double
                                        f000 = src_data(iz    , iy    , ix    ),
                                        f001 = src_data(iz + 1, iy    , ix    ),
                                        f010 = src_data(iz    , iy + 1, ix    ),
                                        f011 = src_data(iz + 1, iy + 1, ix    ),
                                        f100 = src_data(iz    , iy    , ix + 1),
                                        f101 = src_data(iz + 1, iy    , ix + 1),
                                        f110 = src_data(iz    , iy + 1, ix + 1),
                                        f111 = src_data(iz + 1, iy + 1, ix + 1);
                                    rho_array4(i,j,k) = static_cast<amrex::Real> (utils::algorithms::trilinear_interp<double>
                                            (xx0, xx0+file_dx, xx1, xx1+file_dy, xx2, xx2+file_dz,
                                            f000, f001, f010, f011, f100, f101, f110, f111,
                                            x0, x1, x2) );
                                }
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
