//
// Created by juliette on 13/04/23.
//

#include "ParticleHistogram2D.H"

#include "Diagnostics/ReducedDiags/ReducedDiags.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_Config.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_TableData.H>
#include <AMReX_FabConv.H>

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <ostream>
#include <vector>
#include <string>
#include <openPMD/openPMD.hpp>


namespace io = openPMD;

using namespace amrex;

struct NormalizationType {
    enum {
        no_normalization = 0,
        unity_particle_weight,
        max_to_unity,
        area_to_unity
    };
};

// constructor
ParticleHistogram2D::ParticleHistogram2D (std::string rd_name)
        : ReducedDiags{rd_name}
{
    ParmParse pp_rd_name(rd_name);

    // read species
    std::string selected_species_name;
    pp_rd_name.get("species",selected_species_name);

    // read bin parameters
    utils::parser::getWithParser(pp_rd_name, "bin_number_abs",m_bin_num_abs);
    utils::parser::getWithParser(pp_rd_name, "bin_max_abs",   m_bin_max_abs);
    utils::parser::getWithParser(pp_rd_name, "bin_min_abs",   m_bin_min_abs);
    m_bin_size_abs = (m_bin_max_abs - m_bin_min_abs) / m_bin_num_abs;
    utils::parser::getWithParser(pp_rd_name, "bin_number_ord",m_bin_num_ord);
    utils::parser::getWithParser(pp_rd_name, "bin_max_ord",   m_bin_max_ord);
    utils::parser::getWithParser(pp_rd_name, "bin_min_ord",   m_bin_min_ord);
    m_bin_size_ord = (m_bin_max_ord - m_bin_min_ord) / m_bin_num_ord;

    // read histogram function
    utils::parser::Store_parserString(pp_rd_name,"histogram_function_abs(t,x,y,z,ux,uy,uz)",
                                      function_string_abs);
    m_parser_abs = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(function_string_abs,{"t","x","y","z","ux","uy","uz"}));
    utils::parser::Store_parserString(pp_rd_name,"histogram_function_ord(t,x,y,z,ux,uy,uz)",
                                      function_string_ord);
    m_parser_ord = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(function_string_ord,{"t","x","y","z","ux","uy","uz"}));

    // read normalization type
    std::string norm_string = "default";
    pp_rd_name.query("normalization",norm_string);

    // set normalization type
    if ( norm_string == "default" ) {
        m_norm = NormalizationType::no_normalization;
    } else if ( norm_string == "unity_particle_weight" ) {
        m_norm = NormalizationType::unity_particle_weight;
    } else if ( norm_string == "max_to_unity" ) {
        m_norm = NormalizationType::max_to_unity;
    } else if ( norm_string == "area_to_unity" ) {
        m_norm = NormalizationType::area_to_unity;
    } else {
        Abort(Utils::TextMsg::Err(
                "Unknown ParticleHistogram2D normalization type."));
    }

    // get MultiParticleContainer class object
    const auto & mypc = WarpX::GetInstance().GetPartContainer();
    // get species names (std::vector<std::string>)
    auto const species_names = mypc.GetSpeciesNames();
    // select species
    for ( int i = 0; i < mypc.nSpecies(); ++i )
    {
        if ( selected_species_name == species_names[i] ){
            m_selected_species_id = i;
        }
    }
    // if m_selected_species_id is not modified
    if ( m_selected_species_id == -1 ){
        Abort(Utils::TextMsg::Err(
                "Unknown species for ParticleHistogram2D reduced diagnostic."));
    }

    // Read optional filter
    std::string buf;
    m_do_parser_filter = pp_rd_name.query("filter_function(t,x,y,z,ux,uy,uz)", buf);
    if (m_do_parser_filter) {
        utils::parser::Store_parserString(
                pp_rd_name,"filter_function(t,x,y,z,ux,uy,uz)", filter_string);
        m_parser_filter = std::make_unique<amrex::Parser>(
                utils::parser::makeParser(filter_string,{"t","x","y","z","ux","uy","uz"}));
    }
/*
    // Read optional value function
    m_do_parser_value = pp_rd_name.query("value_function(t,x,y,z,ux,uy,uz,w)", buf);
    if (m_do_parser_value) {
        utils::parser::Store_parserString(
                pp_rd_name,"value_function(t,x,y,z,ux,uy,uz,w)", value_string);
        m_w = std::make_unique<amrex::Parser>(
                utils::parser::makeParser(value_string,{"t","x","y","z","ux","uy","uz","w"}));
    }*/
}
// end constructor

// function that computes the histogram
void ParticleHistogram2D::ComputeDiags (int step)
{
    // Judge if the diags should be done at this step
    if (!m_intervals.contains(step+1)) return;

    // resize data array
    Array<int,2> tlo{0,0}; // lower bounds
    Array<int,2> thi{m_bin_num_abs, m_bin_num_ord}; // upper bounds
    amrex::TableData<amrex::Real,2> d_data_2D(tlo, thi);
    m_h_data_2D = amrex::TableData<amrex::Real,2> (tlo, thi, The_Pinned_Arena());
    auto const& h_table_data = m_h_data_2D.table();

    // Initialize data on the host
    for (int i = tlo[0]; i <= thi[0]; ++i) {
        for (int j = tlo[1]; j <= thi[1]; ++j) {
            h_table_data(i,j) = 0.0_rt;
        }
    }

    d_data_2D.copy(m_h_data_2D);
    auto d_table = d_data_2D.table();

    // get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();

    // get time at level 0
    auto const t = warpx.gett_new(0);

    // get MultiParticleContainer class object
    const auto & mypc = warpx.GetPartContainer();

    // get WarpXParticleContainer class object
    auto & myspc = mypc.GetParticleContainer(m_selected_species_id);

    // get parser
    auto fun_partparser_abs =
            utils::parser::compileParser<m_nvars>(m_parser_abs.get());
    auto fun_partparser_ord =
            utils::parser::compileParser<m_nvars>(m_parser_ord.get());

    // get filter parser
    auto fun_filterparser =
            utils::parser::compileParser<m_nvars>(m_parser_filter.get());

    // declare local variables
    auto const num_bins_abs = m_bin_num_abs;
    auto const num_bins_ord = m_bin_num_ord;
    Real const bin_min_abs  = m_bin_min_abs;
    Real const bin_size_abs = m_bin_size_abs;
    Real const bin_min_ord  = m_bin_min_ord;
    Real const bin_size_ord = m_bin_size_ord;
    const bool is_unity_particle_weight =
            (m_norm == NormalizationType::unity_particle_weight) ? true : false;

    bool const do_parser_filter = m_do_parser_filter;

    int const nlevs = std::max(0, myspc.finestLevel()+1);
    for (int lev = 0; lev < nlevs; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        {
            for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti)
            {
                auto const GetPosition = GetParticlePosition(pti);

                auto & attribs = pti.GetAttribs();
                ParticleReal* const AMREX_RESTRICT d_w = attribs[PIdx::w].dataPtr();
                ParticleReal* const AMREX_RESTRICT d_ux = attribs[PIdx::ux].dataPtr();
                ParticleReal* const AMREX_RESTRICT d_uy = attribs[PIdx::uy].dataPtr();
                ParticleReal* const AMREX_RESTRICT d_uz = attribs[PIdx::uz].dataPtr();

                long const np = pti.numParticles();

                //Flag particles that need to be copied if they cross the z_slice
                amrex::ParallelFor(np,
                                   [=] AMREX_GPU_DEVICE(int i)
                                   {
                                       amrex::ParticleReal x, y, z;
                                       GetPosition(i, x, y, z);
                                       auto const w  = (amrex::Real)d_w[i];
                                       auto const ux = d_ux[i] / PhysConst::c;
                                       auto const uy = d_uy[i] / PhysConst::c;
                                       auto const uz = d_uz[i] / PhysConst::c;

                                       // don't count a particle if it is filtered out
                                       if (do_parser_filter)
                                           if (!fun_filterparser(t, x, y, z, ux, uy, uz))
                                               return;

                                       // continue function if particle is not filtered out
                                       auto const f_abs = fun_partparser_abs(t, x, y, z, ux, uy, uz);
                                       auto const f_ord = fun_partparser_ord(t, x, y, z, ux, uy, uz);

                                       // determine particle bin
                                       int const bin_abs = int(Math::floor((f_abs-bin_min_abs)/bin_size_abs));
                                       if ( bin_abs<0 || bin_abs>=num_bins_abs ) return; // discard if out-of-range

                                       int const bin_ord = int(Math::floor((f_ord-bin_min_ord)/bin_size_ord));
                                       if ( bin_ord<0 || bin_ord>=num_bins_ord ) return; // discard if out-of-range

                                       amrex::Real &data = d_table(bin_abs, bin_ord);
                                       if ( is_unity_particle_weight ) {
                                           amrex::HostDevice::Atomic::Add(&data, 1.0_rt);
                                       } else {
                                           amrex::HostDevice::Atomic::Add(&data, w);
                                       }
                                   });
            }
        }
    }

    // Copy data from GPU memory
    m_h_data_2D.copy(d_data_2D);

    // reduced sum over mpi ranks
    int size = static_cast<int> (d_data_2D.size());
    //int size = static_cast<int> (m_h_data_2D.size());
    ParallelDescriptor::ReduceRealSum
            (h_table_data.p, size, ParallelDescriptor::IOProcessorNumber());

    // Return for all that are not IO processor
    if ( !ParallelDescriptor::IOProcessor() ) { return; }
/*
    // normalize the maximum value to be one
    if ( m_norm == NormalizationType::max_to_unity )
    {
        Real f_max = 0.0_rt;
        for (int i = tlo[0]; i <= thi[0]; ++i) {
            for (int j = tlo[1]; j <= thi[1]; ++j) {
                if ( h_table_data(i,j) > f_max ) f_max = h_table_data(i,j);
            }}
        for (int i = tlo[0]; i <= thi[0]; ++i) {
            for (int j = tlo[1]; j <= thi[1]; ++j) {
                if ( f_max > std::numeric_limits<Real>::min() ) h_table_data(i,j) /= f_max;
            }}
        return;
    }

    // normalize the area (integral) to be one
    if ( m_norm == NormalizationType::area_to_unity )
    {
        Real f_area = 0.0_rt;
        for (int i = tlo[0]; i <= thi[0]; ++i) {
            for (int j = tlo[1]; j <= thi[1]; ++j) {
                f_area += h_table_data(i,j) * m_bin_size_abs * m_bin_size_ord;
            }}
        for (int i = tlo[0]; i <= thi[0]; ++i) {
            for (int j = tlo[1]; j <= thi[1]; ++j) {
                if ( f_area > std::numeric_limits<Real>::min() ) h_table_data(i,j) /= f_area;
            }}
    }
    */
}
// end void ParticleHistogram2D::ComputeDiags

void ParticleHistogram2D::WriteToFile (int step) const
{
    // only IO processor writes
    if ( !ParallelDescriptor::IOProcessor() ) { return; }

    // Create the OpenPMD series
    auto series = io::Series(
            m_path + "hist2D/openpmd_%06T.bp",
            io::Access::APPEND);
    auto i = series.iterations[step + 1];
    // record
    auto f_mesh = i.meshes["data"];
    f_mesh.setAttribute("function_abscissa", function_string_abs);
    f_mesh.setAttribute("function_ordinate", function_string_ord);
    f_mesh.setAttribute("filter", filter_string);

    // record components
    auto data = f_mesh[io::RecordComponent::SCALAR];

    // meta data
    f_mesh.setAxisLabels({function_string_ord, function_string_abs}); // ord, abs
    //f_mesh.setAxisLabels({"y", "x"}); // ord, abs
    std::vector< double > const& gridGlobalOffset = {0,0};
    f_mesh.setGridGlobalOffset(gridGlobalOffset);
    f_mesh.setGridSpacing<amrex::Real>({m_bin_size_ord, m_bin_size_abs});

    data.setPosition<amrex::Real>({0, 0});

    auto dataset = io::Dataset(
            io::determineDatatype<double>(),
            {static_cast<unsigned long>(m_bin_num_ord)+1, static_cast<unsigned long>(m_bin_num_abs)+1});
    data.resetDataset(dataset);

    // UNITS ?!
    auto const& h_table_data = m_h_data_2D.table();
    data.storeChunkRaw(
            h_table_data.p,
            {0, 0},
            {static_cast<unsigned long>(m_bin_num_ord)+1, static_cast<unsigned long>(m_bin_num_abs)+1});

    series.flush();
}
