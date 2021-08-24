/* Copyright 2019-2021 Axel Huebl, Junmin Gu
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpXOpenPMD.H"

#include "Diagnostics/ParticleDiag/ParticleDiag.H"
#include "FieldIO.H"
#include "Particles/Filter/FilterFunctors.H"
#include "Utils/RelativeCellPosition.H"
#include "Utils/WarpXAlgorithmSelection.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_ArrayOfStructs.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Config.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Particle.H>
#include <AMReX_Particles.H>
#include <AMReX_Periodicity.H>
#include <AMReX_StructOfArrays.H>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>

namespace detail
{
#ifdef WARPX_USE_OPENPMD
    /** \brief Convert a snake_case string to a camelCase one.
     *
     *  WarpX uses snake_case internally for some component
     *  names, but OpenPMD assumes "_" indicates vector or
     *  tensor fields.
     *
     * @return camelCase version of input
     */
    inline std::string
    snakeToCamel (const std::string& snake_string)
    {
        std::string camelString = snake_string;
        int n = camelString.length();
        for (int x = 0; x < n; x++)
        {
            if (x == 0)
            {
                std::transform(camelString.begin(), camelString.begin()+1, camelString.begin(),
                               [](unsigned char c){ return std::tolower(c); });
            }
            if (camelString[x] == '_')
            {
                std::string tempString = camelString.substr(x + 1, 1);
                std::transform(tempString.begin(), tempString.end(), tempString.begin(),
                               [](unsigned char c){ return std::toupper(c); });
                camelString.erase(x, 2);
                camelString.insert(x, tempString);
            }
        }

        return camelString;
    }

    /** Create the option string
     *
     * @return JSON option string for openPMD::Series
     */
    inline std::string
    getSeriesOptions (std::string const & operator_type,
                      std::map< std::string, std::string > const & operator_parameters)
    {
        std::string options;

        std::string op_parameters;
        for (const auto& kv : operator_parameters) {
            if (!op_parameters.empty()) op_parameters.append(",\n");
            op_parameters.append(std::string(12, ' '))         /* just pretty alignment */
                    .append("\"").append(kv.first).append("\": ")    /* key */
                    .append("\"").append(kv.second).append("\""); /* value (as string) */
        }
        if (!operator_type.empty()) {
            options = R"END(
{
  "adios2": {
    "dataset": {
      "operators": [
        {
          "type": ")END";
            options += operator_type + "\"";
        }
        if (!operator_type.empty() && !op_parameters.empty()) {
            options += R"END(
         ,"parameters": {
)END";
            options += op_parameters + "}";
        }
        if (!operator_type.empty())
            options += R"END(
        }
      ]
    }
  }
}
)END";
        if (options.empty()) options = "{}";
        return options;
    }

    /** Unclutter a real_names to openPMD record
     *
     * @param fullName name as in real_names variable
     * @return pair of openPMD record and component name
     */
    inline std::pair< std::string, std::string >
    name2openPMD ( std::string const& fullName )
    {
        std::string record_name = fullName;
        std::string component_name = openPMD::RecordComponent::SCALAR;

        // we use "_" as separator in names to group vector records
        std::size_t startComp = fullName.find_last_of("_");
        if( startComp != std::string::npos ) {  // non-scalar
            record_name = fullName.substr(0, startComp);
            component_name = fullName.substr(startComp + 1u);
        }
        return make_pair(record_name, component_name);
    }

    /** Return the component labels for particle positions
     */
    inline std::vector< std::string >
    getParticlePositionComponentLabels ()
    {
        using vs = std::vector< std::string >;
#if defined(WARPX_DIM_XZ)
        vs const positionComponents{"x", "z"};
#elif defined(WARPX_DIM_RZ)
        // note: although we internally store particle positions
        //       for AMReX in r,z and a theta attribute, we
        //       actually need them for algorithms (e.g. push)
        //       and I/O in Cartesian.
        //       Other attributes like momentum are consequently
        //       stored in x,y,z internally.
        vs const positionComponents{"x", "y", "z"};
#elif (AMREX_SPACEDIM==3)
        vs const positionComponents{"x", "y", "z"};
#else
#   error Unknown WarpX dimensionality.
#endif
        return positionComponents;
    }

    /** Return the axis (index) names of a mesh
     *
     * This will be returned in C order. This is inverse of the Fortran order
     * of the index labels for the AMReX FArrayBox.
     */
    inline std::vector< std::string >
    getFieldAxisLabels ()
    {
        using vs = std::vector< std::string >;

        // Fortran order of the index labels for the AMReX FArrayBox
#if defined(WARPX_DIM_XZ)
        vs const axisLabels{"x", "z"};
#elif defined(WARPX_DIM_RZ)
        // if we are start to write individual modes
        //vs const axisLabels{"r", "z"};
        // if we just write reconstructed 2D fields at theta=0
        vs const axisLabels{"x", "z"};
#elif (AMREX_SPACEDIM==3)
        vs const axisLabels{"x", "y", "z"};
#else
#   error Unknown WarpX dimensionality.
#endif

        // revert to C order (fastest varying index last)
        return {axisLabels.rbegin(), axisLabels.rend()};
    }

    /** Return the component names of a mesh
     */
    inline std::vector< std::string >
    getFieldComponentLabels ()
    {
        using vs = std::vector< std::string >;
#if defined(WARPX_DIM_RZ)
        // if we are start to write individual modes
        //vs const fieldComponents{"r", "z"};
        // if we just write reconstructed fields at theta=0
        vs const fieldComponents{"x", "y", "z"};
#else
        // note: 1D3V and 2D3V simulations still have 3 components for the fields
        vs const fieldComponents{"x", "y", "z"};
#endif
        return fieldComponents;
    }

    /** Get the openPMD physical dimensionality of a record
     *
     * @param record_name name of the openPMD record
     * @return map with base quantities and power scaling
     */
    inline std::map< openPMD::UnitDimension, double >
    getUnitDimension ( std::string const & record_name )
    {

        if( record_name == "position" ) return {
            {openPMD::UnitDimension::L,  1.}
        };
        else if( record_name == "positionOffset" ) return {
            {openPMD::UnitDimension::L,  1.}
        };
        else if( record_name == "momentum" ) return {
            {openPMD::UnitDimension::L,  1.},
            {openPMD::UnitDimension::M,  1.},
            {openPMD::UnitDimension::T, -1.}
        };
        else if( record_name == "charge" ) return {
            {openPMD::UnitDimension::T,  1.},
            {openPMD::UnitDimension::I,  1.}
        };
        else if( record_name == "mass" ) return {
            {openPMD::UnitDimension::M,  1.}
        };
        else if( record_name == "E" ) return {
            {openPMD::UnitDimension::L,  1.},
            {openPMD::UnitDimension::M,  1.},
            {openPMD::UnitDimension::T, -3.},
            {openPMD::UnitDimension::I, -1.},
        };
        else if( record_name == "B" ) return {
            {openPMD::UnitDimension::M,  1.},
            {openPMD::UnitDimension::I, -1.},
            {openPMD::UnitDimension::T, -2.}
        };
        else return {};
    }

    /** \brief For a given field that is to be written to an openPMD file,
     * set the metadata that indicates the physical unit.
     */
    inline void
    setOpenPMDUnit ( openPMD::Mesh mesh, const std::string field_name )
    {
        if (field_name[0] == 'E'){  // Electric field
            mesh.setUnitDimension({
                                          {openPMD::UnitDimension::L,  1},
                                          {openPMD::UnitDimension::M,  1},
                                          {openPMD::UnitDimension::T, -3},
                                          {openPMD::UnitDimension::I, -1},
                                  });
        } else if (field_name[0] == 'B'){ // Magnetic field
            mesh.setUnitDimension({
                                          {openPMD::UnitDimension::M,  1},
                                          {openPMD::UnitDimension::I, -1},
                                          {openPMD::UnitDimension::T, -2}
                                  });
        } else if (field_name[0] == 'j'){ // current
            mesh.setUnitDimension({
                                          {openPMD::UnitDimension::L, -2},
                                          {openPMD::UnitDimension::I,  1},
                                  });
        } else if (field_name.substr(0,3) == "rho"){ // charge density
            mesh.setUnitDimension({
                                          {openPMD::UnitDimension::L, -3},
                                          {openPMD::UnitDimension::I,  1},
                                          {openPMD::UnitDimension::T,  1},
                                  });
        }
    }
#endif // WARPX_USE_OPENPMD
}

#ifdef WARPX_USE_OPENPMD
WarpXOpenPMDPlot::WarpXOpenPMDPlot (
    openPMD::IterationEncoding ie,
    std::string openPMDFileType,
    std::string operator_type,
    std::map< std::string, std::string > operator_parameters,
    std::vector<bool> fieldPMLdirections)
  :m_Series(nullptr),
   m_Encoding(ie),
   m_OpenPMDFileType(std::move(openPMDFileType)),
   m_fieldPMLdirections(std::move(fieldPMLdirections))
{
  // pick first available backend if default is chosen
  if( m_OpenPMDFileType == "default" )
#if openPMD_HAVE_ADIOS2==1
    m_OpenPMDFileType = "bp";
#elif openPMD_HAVE_ADIOS1==1
    m_OpenPMDFileType = "bp";
#elif openPMD_HAVE_HDF5==1
    m_OpenPMDFileType = "h5";
#else
    m_OpenPMDFileType = "json";
#endif

    m_OpenPMDoptions = detail::getSeriesOptions(operator_type, operator_parameters);
}

WarpXOpenPMDPlot::~WarpXOpenPMDPlot ()
{
  if( m_Series )
  {
    m_Series->flush();
    m_Series.reset( nullptr );
  }
}

std::string
WarpXOpenPMDPlot::GetFileName (std::string& filepath)
{
  filepath.append("/");
  std::string filename = "openpmd";
  //
  // OpenPMD supports timestepped names
  //
  if (m_Encoding == openPMD::IterationEncoding::fileBased) {
      std::string fileSuffix = std::string("_%0") + std::to_string(m_file_min_digits) + std::string("T");
      filename = filename.append(fileSuffix);
  }
  filename.append(".").append(m_OpenPMDFileType);
  filepath.append(filename);
  return filename;
}

void WarpXOpenPMDPlot::SetStep (int ts, const std::string& dirPrefix, int file_min_digits,
                                bool isBTD)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ts >= 0 , "openPMD iterations are unsigned");

    m_dirPrefix = dirPrefix;
    m_file_min_digits = file_min_digits;

    if( ! isBTD ) {
        if (m_CurrentStep >= ts) {
            // note m_Series is reset in Init(), so using m_Series->iterations.contains(ts) is only able to check the
            // last written step in m_Series's life time, but not other earlier written steps by other m_Series
            std::string warnMsg =
                    " Warning from openPMD writer: Already written iteration:" + std::to_string(ts);
            std::cout << warnMsg << std::endl;
            amrex::Warning(warnMsg);
        }
    }

    m_CurrentStep = ts;
    Init(openPMD::Access::CREATE, isBTD);
}

void WarpXOpenPMDPlot::CloseStep (bool isBTD, bool isLastBTDFlush)
{
    // default close is true
    bool callClose = true;
    // close BTD file only when isLastBTDFlush is true
    if (isBTD and !isLastBTDFlush) callClose = false;
    if (callClose) {
        if (m_Series) {
            GetIteration(m_CurrentStep, isBTD).close();
        }

        // create a little helper file for ParaView 5.9+
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            // see Init()
            std::string filepath = m_dirPrefix;
            std::string const filename = GetFileName(filepath);

            std::ofstream pv_helper_file(m_dirPrefix + "/paraview.pmd");
            pv_helper_file << filename << std::endl;
            pv_helper_file.close();
        }
    }
}

void
WarpXOpenPMDPlot::Init (openPMD::Access access, bool isBTD)
{
    if( isBTD && m_Series != nullptr )
        return; // already open for this snapshot (aka timestep in lab frame)

    // either for the next ts file,
    // or init a single file for all ts
    std::string filepath = m_dirPrefix;
    GetFileName(filepath);

    // close a previously open series before creating a new one
    // see ADIOS1 limitation: https://github.com/openPMD/openPMD-api/pull/686
    if ( m_Encoding == openPMD::IterationEncoding::fileBased )
        m_Series = nullptr;
    else if ( m_Series != nullptr )
        return;

    if (amrex::ParallelDescriptor::NProcs() > 1) {
#if defined(AMREX_USE_MPI)
        m_Series = std::make_unique<openPMD::Series>(
                filepath, access,
                amrex::ParallelDescriptor::Communicator(),
                m_OpenPMDoptions
        );
        m_MPISize = amrex::ParallelDescriptor::NProcs();
        m_MPIRank = amrex::ParallelDescriptor::MyProc();
#else
        amrex::Abort("openPMD-api not built with MPI support!");
#endif
    } else {
        m_Series = std::make_unique<openPMD::Series>(filepath, access, m_OpenPMDoptions);
        m_MPISize = 1;
        m_MPIRank = 1;
    }

    m_Series->setIterationEncoding( m_Encoding );

    // input file / simulation setup author
    if( !WarpX::authors.empty())
        m_Series->setAuthor( WarpX::authors );
    // more natural naming for PIC
    m_Series->setMeshesPath( "fields" );
    // conform to ED-PIC extension of openPMD
    uint32_t const openPMD_ED_PIC = 1u;
    m_Series->setOpenPMDextension( openPMD_ED_PIC );
    // meta info
    m_Series->setSoftware( "WarpX", WarpX::Version() );
}

void
WarpXOpenPMDPlot::WriteOpenPMDParticles (const amrex::Vector<ParticleDiag>& particle_diags,
                                         const bool isBTD)
{
  WARPX_PROFILE("WarpXOpenPMDPlot::WriteOpenPMDParticles()");

  for (unsigned i = 0, n = particle_diags.size(); i < n; ++i) {
    WarpXParticleContainer* pc = particle_diags[i].getParticleContainer();
    auto tmp = ParticleBuffer::getTmpPC<amrex::PinnedArenaAllocator>(pc);
    // names of amrex::Real and int particle attributes in SoA data
    amrex::Vector<std::string> real_names;
    amrex::Vector<std::string> int_names;
    amrex::Vector<int> int_flags;
    amrex::Vector<int> real_flags;

    // see openPMD ED-PIC extension for namings
    // note: an underscore separates the record name from its component
    //       for non-scalar records
    real_names.push_back("weighting");

    real_names.push_back("momentum_x");
    real_names.push_back("momentum_y");
    real_names.push_back("momentum_z");

#ifdef WARPX_DIM_RZ
    real_names.push_back("theta");
#endif

    // get the names of the real comps
    real_names.resize(pc->NumRealComps());
    auto runtime_rnames = pc->getParticleRuntimeComps();
    for (auto const& x : runtime_rnames)
    {
        real_names[x.second+PIdx::nattribs] = detail::snakeToCamel(x.first);
    }

    // plot any "extra" fields by default
    real_flags = particle_diags[i].plot_flags;
    real_flags.resize(pc->NumRealComps(), 1);

    // and the names
    int_names.resize(pc->NumIntComps());
    auto runtime_inames = pc->getParticleRuntimeiComps();
    for (auto const& x : runtime_inames)
    {
        int_names[x.second+0] = detail::snakeToCamel(x.first);
    }

    // plot by default
    int_flags.resize(pc->NumIntComps(), 1);

      pc->ConvertUnits(ConvertDirection::WarpX_to_SI);

      RandomFilter const random_filter(particle_diags[i].m_do_random_filter,
                                       particle_diags[i].m_random_fraction);
      UniformFilter const uniform_filter(particle_diags[i].m_do_uniform_filter,
                                         particle_diags[i].m_uniform_stride);
      ParserFilter parser_filter(particle_diags[i].m_do_parser_filter,
                                 compileParser<ParticleDiag::m_nvars>
                                     (particle_diags[i].m_particle_filter_parser.get()),
                                 pc->getMass());
      parser_filter.m_units = InputUnits::SI;
      GeometryFilter const geometry_filter(particle_diags[i].m_do_geom_filter,
                                           particle_diags[i].m_diag_domain);

      using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
      tmp.copyParticles(*pc,
                        [=] AMREX_GPU_HOST_DEVICE (const SrcData& src, int ip, const amrex::RandomEngine& engine)
      {
          const SuperParticleType& p = src.getSuperParticle(ip);
          return random_filter(p, engine) * uniform_filter(p, engine)
                 * parser_filter(p, engine) * geometry_filter(p, engine);
      }, true);

    // real_names contains a list of all real particle attributes.
    // real_flags is 1 or 0, whether quantity is dumped or not.

    {
      DumpToFile(&tmp,
         particle_diags[i].getSpeciesName(),
         m_CurrentStep,
         real_flags,
         int_flags,
         real_names, int_names,
         pc->getCharge(), pc->getMass(),
         isBTD
      );
    }

    // Convert momentum back to WarpX units
    pc->ConvertUnits(ConvertDirection::SI_to_WarpX);
  }
}

void
WarpXOpenPMDPlot::DumpToFile (ParticleContainer* pc,
                    const std::string& name,
                    int iteration,
                    const amrex::Vector<int>& write_real_comp,
                    const amrex::Vector<int>& write_int_comp,
                    const amrex::Vector<std::string>& real_comp_names,
                    const amrex::Vector<std::string>&  int_comp_names,
                    amrex::ParticleReal const charge,
                    amrex::ParticleReal const mass,
                    const bool isBTD) const
{
  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Series != nullptr, "openPMD: series must be initialized");

  WarpXParticleCounter counter(pc);
  openPMD::Iteration currIteration = GetIteration(iteration, isBTD);

  openPMD::ParticleSpecies currSpecies = currIteration.particles[name];
  // meta data for ED-PIC extension
  currSpecies.setAttribute( "particleShape", double( WarpX::noz ) );
  // TODO allow this per direction in the openPMD standard, ED-PIC extension?
  currSpecies.setAttribute( "particleShapes", [](){
      return std::vector< double >{
          double(WarpX::nox),
#if AMREX_SPACEDIM==3
          double(WarpX::noy),
#endif
          double(WarpX::noz)
      };
  }() );
  currSpecies.setAttribute( "particlePush", [](){
      switch( WarpX::particle_pusher_algo ) {
          case ParticlePusherAlgo::Boris : return "Boris";
          case ParticlePusherAlgo::Vay : return "Vay";
          case ParticlePusherAlgo::HigueraCary : return "HigueraCary";
          default: return "other";
      }
  }() );
  currSpecies.setAttribute( "particleInterpolation", [](){
      switch( WarpX::field_gathering_algo ) {
          case GatheringAlgo::EnergyConserving : return "energyConserving";
          case GatheringAlgo::MomentumConserving : return "momentumConserving";
          default: return "other";
      }
  }() );
  currSpecies.setAttribute( "particleSmoothing", "none" );
  currSpecies.setAttribute( "currentDeposition", [](){
      switch( WarpX::current_deposition_algo ) {
          case CurrentDepositionAlgo::Esirkepov : return "Esirkepov";
          case CurrentDepositionAlgo::Vay : return "Vay";
          default: return "directMorseNielson";
      }
  }() );

  //
  // define positions & offsets
  //
  SetupPos(currSpecies, counter.GetTotalNumParticles(), charge, mass);
  SetupRealProperties(currSpecies, write_real_comp, real_comp_names, write_int_comp, int_comp_names, counter.GetTotalNumParticles());

  // open files from all processors, in case some will not contribute below
  m_Series->flush();

  for (auto currentLevel = 0; currentLevel <= pc->finestLevel(); currentLevel++)
    {
      uint64_t offset = static_cast<uint64_t>( counter.m_ParticleOffsetAtRank[currentLevel] );

      for (ParticleIter pti(*pc, currentLevel); pti.isValid(); ++pti) {
         auto const numParticleOnTile = pti.numParticles();
         uint64_t const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );

         // get position and particle ID from aos
         // note: this implementation iterates the AoS 4x...
         // if we flush late as we do now, we can also copy out the data in one go
         const auto& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
         {
           // Save positions
           auto const positionComponents = detail::getParticlePositionComponentLabels();
#if defined(WARPX_DIM_RZ)
           {
              std::shared_ptr<amrex::ParticleReal> z(
                      new amrex::ParticleReal[numParticleOnTile],
                      [](amrex::ParticleReal const *p) { delete[] p; }
              );
              for (auto i = 0; i < numParticleOnTile; i++)
                  z.get()[i] = aos[i].pos(1);  // {0: "r", 1: "z"}
              std::string const positionComponent = "z";
              currSpecies["position"]["z"].storeChunk(z, {offset}, {numParticleOnTile64});
           }

           //   reconstruct x and y from polar coordinates r, theta
           auto const& soa = pti.GetStructOfArrays();
           amrex::ParticleReal const* theta = soa.GetRealData(PIdx::theta).dataPtr();
           AMREX_ALWAYS_ASSERT_WITH_MESSAGE(theta != nullptr, "openPMD: invalid theta pointer.");
           AMREX_ALWAYS_ASSERT_WITH_MESSAGE(int(soa.GetRealData(PIdx::theta).size()) == numParticleOnTile,
                                            "openPMD: theta and tile size do not match");
           {
               std::shared_ptr< amrex::ParticleReal > x(
                       new amrex::ParticleReal[numParticleOnTile],
                       [](amrex::ParticleReal const *p){ delete[] p; }
               );
               std::shared_ptr< amrex::ParticleReal > y(
                       new amrex::ParticleReal[numParticleOnTile],
                       [](amrex::ParticleReal const *p){ delete[] p; }
               );
               for (auto i=0; i<numParticleOnTile; i++) {
                   auto const r = aos[i].pos(0);  // {0: "r", 1: "z"}
                   x.get()[i] = r * std::cos(theta[i]);
                   y.get()[i] = r * std::sin(theta[i]);
               }
               currSpecies["position"]["x"].storeChunk(x, {offset}, {numParticleOnTile64});
               currSpecies["position"]["y"].storeChunk(y, {offset}, {numParticleOnTile64});
           }
#else
           for (auto currDim = 0; currDim < AMREX_SPACEDIM; currDim++) {
                std::shared_ptr< amrex::ParticleReal > curr(
                    new amrex::ParticleReal[numParticleOnTile],
                    [](amrex::ParticleReal const *p){ delete[] p; }
                );
                for (auto i=0; i<numParticleOnTile; i++) {
                     curr.get()[i] = aos[i].pos(currDim);
                }
                std::string const positionComponent = positionComponents[currDim];
                currSpecies["position"][positionComponent].storeChunk(curr, {offset}, {numParticleOnTile64});
           }
#endif

           // save particle ID after converting it to a globally unique ID
           std::shared_ptr< uint64_t > ids(
               new uint64_t[numParticleOnTile],
               [](uint64_t const *p){ delete[] p; }
           );
           for (auto i=0; i<numParticleOnTile; i++) {
               ids.get()[i] = WarpXUtilIO::localIDtoGlobal( aos[i].id(), aos[i].cpu() );
           }
           auto const scalar = openPMD::RecordComponent::SCALAR;
           currSpecies["id"][scalar].storeChunk(ids, {offset}, {numParticleOnTile64});
         }
         //  save "extra" particle properties in AoS and SoA
         SaveRealProperty(pti,
             currSpecies,
             offset,
             write_real_comp, real_comp_names,
             write_int_comp, int_comp_names);

         offset += numParticleOnTile64;
      }
    }
    m_Series->flush();
}

void
WarpXOpenPMDPlot::SetupRealProperties (openPMD::ParticleSpecies& currSpecies,
                      const amrex::Vector<int>& write_real_comp,
                      const amrex::Vector<std::string>& real_comp_names,
                      const amrex::Vector<int>& write_int_comp,
                      const amrex::Vector<std::string>& int_comp_names,
                      unsigned long long np) const
{
    auto dtype_real = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np});
    auto dtype_int  = openPMD::Dataset(openPMD::determineDatatype<int>(), {np});

    //
    // the beam/input3d showed write_real_comp.size() = 16 while only 10 real comp names
    // so using the min to be safe.
    //
    auto const getComponentRecord = [&currSpecies](std::string const comp_name) {
        // handle scalar and non-scalar records by name
        std::string record_name, component_name;
        std::tie(record_name, component_name) = detail::name2openPMD(comp_name);
        return currSpecies[record_name][component_name];
    };
    auto const real_counter = std::min(write_real_comp.size(), real_comp_names.size());
    for (int i = 0; i < real_counter; ++i) {
      if (write_real_comp[i]) {
          getComponentRecord(real_comp_names[i]).resetDataset(dtype_real);
      }
    }
    auto const int_counter = std::min(write_int_comp.size(), int_comp_names.size());
    for (int i = 0; i < int_counter; ++i) {
        if (write_int_comp[i]) {
            getComponentRecord(int_comp_names[i]).resetDataset(dtype_int);
        }
    }

    std::set< std::string > addedRecords; // add meta-data per record only once
    for (auto idx=0; idx<m_NumSoARealAttributes; idx++) {
        auto ii = m_NumAoSRealAttributes + idx; // jump over AoS names
        if (write_real_comp[ii]) {
            // handle scalar and non-scalar records by name
            std::string record_name, component_name;
            std::tie(record_name, component_name) = detail::name2openPMD(real_comp_names[ii]);
            auto currRecord = currSpecies[record_name];

            // meta data for ED-PIC extension
            bool newRecord = false;
            std::tie(std::ignore, newRecord) = addedRecords.insert(record_name);
            if( newRecord ) {
                currRecord.setUnitDimension( detail::getUnitDimension(record_name) );
                if( record_name == "weighting" )
                    currRecord.setAttribute( "macroWeighted", 1u );
                else
                    currRecord.setAttribute( "macroWeighted", 0u );
                if( record_name == "momentum" || record_name == "weighting" )
                    currRecord.setAttribute( "weightingPower", 1.0 );
                else
                    currRecord.setAttribute( "weightingPower", 0.0 );
            }
        }
    }
    for (auto idx=0; idx<int_counter; idx++) {
        auto ii = m_NumAoSIntAttributes + idx; // jump over AoS names
        if (write_int_comp[ii]) {
            // handle scalar and non-scalar records by name
            std::string record_name, component_name;
            std::tie(record_name, component_name) = detail::name2openPMD(int_comp_names[ii]);
            auto currRecord = currSpecies[record_name];

            // meta data for ED-PIC extension
            bool newRecord = false;
            std::tie(std::ignore, newRecord) = addedRecords.insert(record_name);
            if( newRecord ) {
                currRecord.setUnitDimension( detail::getUnitDimension(record_name) );
                currRecord.setAttribute( "macroWeighted", 0u );
                if( record_name == "momentum" || record_name == "weighting" )
                    currRecord.setAttribute( "weightingPower", 1.0 );
                else
                    currRecord.setAttribute( "weightingPower", 0.0 );
            }
        }
    }
}

void
WarpXOpenPMDPlot::SaveRealProperty (ParticleIter& pti,
                       openPMD::ParticleSpecies& currSpecies,
                       unsigned long long const offset,
                       amrex::Vector<int> const& write_real_comp,
                       amrex::Vector<std::string> const& real_comp_names,
                       amrex::Vector<int> const& write_int_comp,
                       amrex::Vector<std::string> const& int_comp_names) const

{
  int numOutputReal = 0;
  int const totalRealAttrs = m_NumAoSRealAttributes + m_NumSoARealAttributes;

  for( int i = 0; i < totalRealAttrs; ++i )
    if( write_real_comp[i] )
      ++numOutputReal;

  auto const numParticleOnTile = pti.numParticles();
  uint64_t const numParticleOnTile64 = static_cast<uint64_t>( numParticleOnTile );
  auto const& aos = pti.GetArrayOfStructs();  // size =  numParticlesOnTile
  auto const& soa = pti.GetStructOfArrays();

  // first we concatinate the AoS into contiguous arrays
  {
    for( auto idx=0; idx<m_NumAoSRealAttributes; idx++ ) {
      if( write_real_comp[idx] ) {
          // handle scalar and non-scalar records by name
          std::string record_name, component_name;
          std::tie(record_name, component_name) = detail::name2openPMD(real_comp_names[idx]);
          auto currRecord = currSpecies[record_name];
          auto currRecordComp = currRecord[component_name];

          std::shared_ptr< amrex::ParticleReal > d(
              new amrex::ParticleReal[numParticleOnTile],
              [](amrex::ParticleReal const *p){ delete[] p; }
          );

          for( auto kk=0; kk<numParticleOnTile; kk++ )
               d.get()[kk] = aos[kk].rdata(idx);

          currRecordComp.storeChunk(d,
               {offset}, {numParticleOnTile64});
      }
    }
  }

  auto const getComponentRecord = [&currSpecies](std::string const comp_name) {
    // handle scalar and non-scalar records by name
    std::string record_name, component_name;
    std::tie(record_name, component_name) = detail::name2openPMD(comp_name);
    return currSpecies[record_name][component_name];
  };

  // here we the save the SoA properties (real)
  {
    auto const real_counter = std::min(write_real_comp.size(), real_comp_names.size());
    for (auto idx=0; idx<real_counter; idx++) {
      auto ii = m_NumAoSRealAttributes + idx;
      if (write_real_comp[ii]) {
        getComponentRecord(real_comp_names[ii]).storeChunk(openPMD::shareRaw(soa.GetRealData(idx)),
          {offset}, {numParticleOnTile64});
      }
    }
  }
  // and now SoA int properties
  {
    auto const int_counter = std::min(write_int_comp.size(), int_comp_names.size());
    for (auto idx=0; idx<int_counter; idx++) {
      auto ii = m_NumAoSIntAttributes + idx; // jump over AoS names
      if (write_int_comp[ii]) {
        getComponentRecord(int_comp_names[ii]).storeChunk(openPMD::shareRaw(soa.GetIntData(idx)),
          {offset}, {numParticleOnTile64});
      }
    }
  }
}


void
WarpXOpenPMDPlot::SetupPos (
    openPMD::ParticleSpecies& currSpecies,
    const unsigned long long& np,
    amrex::ParticleReal const charge,
    amrex::ParticleReal const mass) const
{
  auto const realType = openPMD::Dataset(openPMD::determineDatatype<amrex::ParticleReal>(), {np});
  auto const idType = openPMD::Dataset(openPMD::determineDatatype< uint64_t >(), {np});

  auto const positionComponents = detail::getParticlePositionComponentLabels();
  for( auto const& comp : positionComponents ) {
      currSpecies["positionOffset"][comp].resetDataset( realType );
      currSpecies["positionOffset"][comp].makeConstant( 0. );
      currSpecies["position"][comp].resetDataset( realType );
  }

  auto const scalar = openPMD::RecordComponent::SCALAR;
  currSpecies["id"][scalar].resetDataset( idType );
  currSpecies["charge"][scalar].resetDataset( realType );
  currSpecies["charge"][scalar].makeConstant( charge );
  currSpecies["mass"][scalar].resetDataset( realType );
  currSpecies["mass"][scalar].makeConstant( mass );

  // meta data
  currSpecies["position"].setUnitDimension( detail::getUnitDimension("position") );
  currSpecies["positionOffset"].setUnitDimension( detail::getUnitDimension("positionOffset") );
  currSpecies["charge"].setUnitDimension( detail::getUnitDimension("charge") );
  currSpecies["mass"].setUnitDimension( detail::getUnitDimension("mass") );

  // meta data for ED-PIC extension
  currSpecies["position"].setAttribute( "macroWeighted", 0u );
  currSpecies["position"].setAttribute( "weightingPower", 0.0 );
  currSpecies["positionOffset"].setAttribute( "macroWeighted", 0u );
  currSpecies["positionOffset"].setAttribute( "weightingPower", 0.0 );
  currSpecies["id"].setAttribute( "macroWeighted", 0u );
  currSpecies["id"].setAttribute( "weightingPower", 0.0 );
  currSpecies["charge"].setAttribute( "macroWeighted", 0u );
  currSpecies["charge"].setAttribute( "weightingPower", 1.0 );
  currSpecies["mass"].setAttribute( "macroWeighted", 0u );
  currSpecies["mass"].setAttribute( "weightingPower", 1.0 );
}


/*
 * Set up parameter for mesh container using the geometry (from level 0)
 *
 * @param [IN] meshes: openPMD-api mesh container
 * @param [IN] full_geom: field geometry
 *
 */
void
WarpXOpenPMDPlot::SetupFields ( openPMD::Container< openPMD::Mesh >& meshes,
                                amrex::Geometry& full_geom ) const
{
      // meta data for ED-PIC extension
      auto const period = full_geom.periodicity(); // TODO double-check: is this the proper global bound or of some level?
      std::vector<std::string> fieldBoundary(6, "reflecting");
      std::vector<std::string> particleBoundary(6, "absorbing");
#if AMREX_SPACEDIM != 3
      fieldBoundary.resize(4);
      particleBoundary.resize(4);
#endif

      for (auto i = 0u; i < fieldBoundary.size() / 2u; ++i)
          if (m_fieldPMLdirections.at(i))
              fieldBoundary.at(i) = "open";

      for (auto i = 0u; i < fieldBoundary.size() / 2u; ++i)
          if (period.isPeriodic(i)) {
              fieldBoundary.at(2u * i) = "periodic";
              fieldBoundary.at(2u * i + 1u) = "periodic";
              particleBoundary.at(2u * i) = "periodic";
              particleBoundary.at(2u * i + 1u) = "periodic";
          }

      meshes.setAttribute("fieldSolver", []() {
          switch (WarpX::maxwell_solver_id) {
              case MaxwellSolverAlgo::Yee :
                  return "Yee";
              case MaxwellSolverAlgo::CKC :
                  return "CK";
              case MaxwellSolverAlgo::PSATD :
                  return "PSATD";
              default:
                  return "other";
          }
      }());
      meshes.setAttribute("fieldBoundary", fieldBoundary);
      meshes.setAttribute("particleBoundary", particleBoundary);
      meshes.setAttribute("currentSmoothing", []() {
          if (WarpX::use_filter) return "Binomial";
          else return "none";
      }());
      if (WarpX::use_filter)
          meshes.setAttribute("currentSmoothingParameters", []() {
              std::stringstream ss;
              ss << "period=1;compensator=false";
              ss << ";numPasses_x=" << WarpX::filter_npass_each_dir[0];
#if (AMREX_SPACEDIM == 3)
              ss << ";numPasses_y=" << WarpX::filter_npass_each_dir[1];
              ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[2];
#else
              ss << ";numPasses_z=" << WarpX::filter_npass_each_dir[1];
#endif
              std::string currentSmoothingParameters = ss.str();
              return currentSmoothingParameters;
          }());
      meshes.setAttribute("chargeCorrection", []() {
          if (WarpX::do_dive_cleaning) return "hyperbolic"; // TODO or "spectral" or something? double-check
          else return "none";
      }());
      if (WarpX::do_dive_cleaning)
          meshes.setAttribute("chargeCorrectionParameters", "period=1");
}


/*
 * Setup component properties  for a field mesh
 * @param [IN]: mesh          a mesh field
 * @param [IN]: full_geom     geometry for the mesh
 * @param [IN]: mesh_comp     a component for the mesh
 */
void
WarpXOpenPMDPlot::SetupMeshComp (openPMD::Mesh& mesh,
                                 amrex::Geometry& full_geom,
                                 openPMD::MeshRecordComponent& mesh_comp) const
{
       amrex::Box const & global_box = full_geom.Domain();
       auto const global_size = getReversedVec(global_box.size());
       // - Grid spacing
       std::vector<double> const grid_spacing = getReversedVec(full_geom.CellSize());
       // - Global offset
       std::vector<double> const global_offset = getReversedVec(full_geom.ProbLo());
       // - AxisLabels
       std::vector<std::string> axis_labels = detail::getFieldAxisLabels();

       // Prepare the type of dataset that will be written
       openPMD::Datatype const datatype = openPMD::determineDatatype<amrex::Real>();
       auto const dataset = openPMD::Dataset(datatype, global_size);

       mesh.setDataOrder(openPMD::Mesh::DataOrder::C);
       mesh.setAxisLabels(axis_labels);
       mesh.setGridSpacing(grid_spacing);
       mesh.setGridGlobalOffset(global_offset);
       mesh.setAttribute("fieldSmoothing", "none");
       mesh_comp.resetDataset(dataset);

}

/*
 * Get component names of a field for openPMD-api book-keeping
 * Level is reflected as _lvl<meshLevel>
 *
 * @param meshLevel [IN]:    level of mesh
 * @param varname [IN]:      name from WarpX
 * @param field_name [OUT]:  field name for openPMD-api output
 * @param comp_name [OUT]:   comp name for openPMD-api output
 */
void
WarpXOpenPMDPlot::GetMeshCompNames (int meshLevel,
                                    const std::string& varname,
                                    std::string& field_name,
                                    std::string& comp_name) const
{
    if (varname.size() >= 2u ) {
        std::string const varname_1st = varname.substr(0u, 1u); // 1st character
        std::string const varname_2nd = varname.substr(1u, 1u); // 2nd character

        // Check if this field is a vector. If so, then extract the field name
        std::vector< std::string > const vector_fields = {"E", "B", "j"};
        std::vector< std::string > const field_components = detail::getFieldComponentLabels();
        for( std::string const& field : vector_fields ) {
            for( std::string const& component : field_components ) {
                if( field.compare( varname_1st ) == 0 &&
                    component.compare( varname_2nd ) == 0 )
                {
                    field_name = varname_1st + varname.substr(2); // Strip component
                    comp_name = varname_2nd;
                }
            }
        }
    }

    if ( 0 == meshLevel )
      return;

    field_name += std::string("_lvl").append(std::to_string(meshLevel));
}
/*
 * Write Field with all mesh levels
 *
 */
void
WarpXOpenPMDPlot::WriteOpenPMDFieldsAll ( //const std::string& filename,
                      const std::vector<std::string>& varnames,
                      const amrex::Vector<amrex::MultiFab>& mf,
                      amrex::Vector<amrex::Geometry>& geom,
                      const int iteration,
                      const double time, bool isBTD,
                      const amrex::Geometry& full_BTD_snapshot ) const
{
  //This is AMReX's tiny profiler. Possibly will apply it later
  WARPX_PROFILE("WarpXOpenPMDPlot::WriteOpenPMDFields()");

  AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_Series != nullptr, "openPMD series must be initialized");

  // is this either a regular write (true) or the first write in a
  // backtransformed diagnostic (BTD):
  bool const first_write_to_iteration = ! m_Series->iterations.contains( iteration );

  // meta data
  openPMD::Iteration series_iteration = GetIteration(m_CurrentStep, isBTD);

  // collective open
  series_iteration.open();

  auto meshes = series_iteration.meshes;
  if (first_write_to_iteration) {
    // lets see whether full_geom varies from geom[0]   xgeom[1]
    series_iteration.setTime( time );
  }

  for (int i=0; i<geom.size(); i++) {
    amrex::Geometry full_geom = geom[i];
    if( isBTD )
      full_geom = full_BTD_snapshot;

    // setup is called once. So it uses property "period" from first
    // geometry for <all> field levels.
    if ( (0 == i) && first_write_to_iteration )
      SetupFields(meshes, full_geom);

    amrex::Box const & global_box = full_geom.Domain();
    auto const global_size = getReversedVec(global_box.size());

    int const ncomp = mf[i].nComp();
    for ( int icomp=0; icomp<ncomp; icomp++ ) {
          std::string const & varname = varnames[icomp];

          // assume fields are scalar unless they match the following match of known vector fields
          std::string field_name = varname;
          std::string comp_name = openPMD::MeshRecordComponent::SCALAR;
          GetMeshCompNames( i, varname, field_name, comp_name );

          auto mesh = meshes[field_name];
          auto mesh_comp = mesh[comp_name];
          if ( first_write_to_iteration )
          {
             SetupMeshComp( mesh, full_geom, mesh_comp );
             detail::setOpenPMDUnit( mesh, field_name );

             auto relative_cell_pos = utils::getRelativeCellPosition(mf[i]);     // AMReX Fortran index order
             std::reverse( relative_cell_pos.begin(), relative_cell_pos.end() ); // now in C order
             mesh_comp.setPosition( relative_cell_pos );
           }

           // Loop through the multifab, and store each box as a chunk,
           // in the openPMD file.
           for( amrex::MFIter mfi(mf[i]); mfi.isValid(); ++mfi )
           {
                amrex::FArrayBox const& fab = mf[i][mfi];
                amrex::Box const& local_box = fab.box();

                // Determine the offset and size of this chunk
                amrex::IntVect const box_offset = local_box.smallEnd() - global_box.smallEnd();
                auto const chunk_offset = getReversedVec( box_offset );
                auto const chunk_size = getReversedVec( local_box.size() );

                amrex::Real const * local_data = fab.dataPtr( icomp );
                mesh_comp.storeChunk( openPMD::shareRaw(local_data),
                                      chunk_offset, chunk_size );
            }
    } // icomp loop
    // Flush data to disk after looping over all components
    m_Series->flush();
  } // levels loop (i)
}
#endif // WARPX_USE_OPENPMD



//
//
//
WarpXParticleCounter::WarpXParticleCounter (ParticleContainer* pc)
{
  m_MPISize = amrex::ParallelDescriptor::NProcs();
  m_MPIRank = amrex::ParallelDescriptor::MyProc();

  m_ParticleCounterByLevel.resize(pc->finestLevel()+1);
  m_ParticleOffsetAtRank.resize(pc->finestLevel()+1);
  m_ParticleSizeAtRank.resize(pc->finestLevel()+1);

  for (auto currentLevel = 0; currentLevel <= pc->finestLevel(); currentLevel++)
    {
      long numParticles = 0; // numParticles in this processor

      for (ParticleIter pti(*pc, currentLevel); pti.isValid(); ++pti) {
          auto numParticleOnTile = pti.numParticles();
          numParticles += numParticleOnTile;
      }

      unsigned long long offset=0; // offset of this level
      unsigned long long sum=0; // numParticles in this level (sum from all processors)

      GetParticleOffsetOfProcessor(numParticles, offset,  sum);

      m_ParticleCounterByLevel[currentLevel] = sum;
      m_ParticleOffsetAtRank[currentLevel] = offset;
      m_ParticleSizeAtRank[currentLevel] = numParticles;

      // adjust offset, it should be numbered after particles from previous levels
      for (auto lv=0; lv<currentLevel; lv++)
    m_ParticleOffsetAtRank[currentLevel] += m_ParticleCounterByLevel[lv];

      m_Total += sum;
    }
}


// get the offset in the overall particle id collection
//
// note: this is a MPI-collective operation
//
// input: num of particles  of from each   processor
//
// output:
//     offset within <all> the particles in the comm
//     sum of all particles in the comm
//
void
WarpXParticleCounter::GetParticleOffsetOfProcessor (
    const long& numParticles,
    unsigned long long& offset,
    unsigned long long& sum
) const
{
    offset = 0;
#if defined(AMREX_USE_MPI)
    std::vector<long> result(m_MPISize, 0);
    amrex::ParallelGather::Gather (numParticles, result.data(), -1, amrex::ParallelDescriptor::Communicator());

    sum = 0;
    int const num_results = result.size();
    for (int i=0; i<num_results; i++) {
        sum += result[i];
        if (i<m_MPIRank)
            offset += result[i];
    }
#else
    sum = numParticles;
#endif
}
