#ifndef WARPX_FLUSHFORMATPLOTPLUS_H_
#define WARPX_FLUSHFORMATPLOTPLUS_H_

#include "FlushFormat.H"

#include "Diagnostics/ParticleDiag/ParticleDiag_fwd.H"
#include "Particles/Filter/FilterFunctors.H"
#include "Particles/ParticleIO.H"
#include "Utils/Parser/ParserUtils.H"

#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#include <iosfwd>
#include <string>


#include <AMReX_PlotFileUtil.H> // for AMReXWithOpenPMD
#include <AMReX_MultiFab.H> // for StoreMesh's multi  fab

#ifdef AMREX_USE_OPENPMD_API
class AMReXWithOpenPMD
{
public:
  AMReXWithOpenPMD(const std::string& prefix);
  ~AMReXWithOpenPMD();

  void SetWriter(amrex::openpmd_api::AMReX_openPMDWriter* w);

  bool InitLocalHandler(const std::string& prefix);

  void StoreMesh (const amrex::Vector<const amrex::MultiFab*> &mf,
          const amrex::Vector<std::string> &varnames,
          const amrex::Vector<amrex::Geometry> &geom,
          amrex::Real time);

  std::unique_ptr<amrex::openpmd_api::AMReX_openPMDHandler> m_UserHandler;

private:
  std::string m_Prefix = "";
};



/**
 * \brief This class aims at dumping diags data to disk via the AMReX-openPMD I/O module
 * In particular, function WriteToFile takes fields and particles as input arguments,
 * and writes data to file.
 */
class FlushFormatPlotPlus : public FlushFormat
{
public:
    /** Flush fields and particles to plotfile */
    virtual void WriteToFile (
        const amrex::Vector<std::string>& varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
        const std::string prefix, int file_min_digits,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        const bool use_pinned_pc = false,
        bool isBTD = false, int snapshotID = -1,
        int bufferID = 1, int numBuffers = 1,
        const amrex::Geometry& full_BTD_snapshot = amrex::Geometry(),
        bool isLastBTDFlush = false,
        const amrex::Vector<int>& totalParticlesFlushedAlready = amrex::Vector<int>() ) const override;

  /*
    void WriteJobInfo(const std::string& dir) const;
    void WriteWarpXHeader(const std::string& name, amrex::Vector<amrex::Geometry>& geom) const;
    void WriteAllRawFields (const bool plot_raw_fields, const int nlevels,
                            const std::string& plotfilename,
                            const bool plot_raw_fields_guards) const;
  */
    /** \brief Write particles data to file.
     * \param[in] filename name of output directory
     * \param[in] particle_diags Each element of this vector handles output of 1 species.
     * \param[in] time the simulation time on the coarsest level
     * \param[in] isBTD whether this is a back-transformed diagnostic
     */
    void WriteParticles(const std::string& filename,
                        const amrex::Vector<ParticleDiag>& particle_diags,
                        const amrex::Real time,
                        bool isBTD = false) const;

    ~FlushFormatPlotPlus() {}

  std::unique_ptr< AMReXWithOpenPMD > m_Writer = std::make_unique< AMReXWithOpenPMD >("");

  void GetNames(const ParticleDiag& part_diag,
        amrex::Vector<std::string>& real_names,
        amrex::Vector<std::string>& int_names,
        amrex::Vector<int>& int_flags,
        amrex::Vector<int>& real_flags
        ) const;

  void DefaultWriter(const std::string& prefix,
             int iteration,
             const amrex::Vector<std::string> varnames,
             const amrex::Vector<amrex::MultiFab>& mf,
             amrex::Vector<amrex::Geometry>& geom,
             const double time,
             const amrex::Vector<ParticleDiag>& particle_diags,
             const bool use_pinned_pc
             ) const;
  void BTDWriter(const std::string& prefix,
         int iteration,
         const amrex::Vector<std::string> varnames,
         const amrex::Vector<amrex::MultiFab>& mf,
         amrex::Vector<amrex::Geometry>& geom,
         const double time,
         const amrex::Vector<ParticleDiag>& particle_diags,
         const bool use_pinned_pc,
         bool isLastBTDFlush,
         const amrex::Vector<int>& totalParticlesFlushedAlready
         ) const;

//
// copies ptls
// general purpose for both BTD and default use cases
//
template< typename T_ParticleContainer >
  void CopyPtls(PinnedMemoryParticleContainer& tmp,
                const double time,
                WarpXParticleContainer* pc,
                T_ParticleContainer* ptlSrc,
                const ParticleDiag& part_diag) const
  {
      const auto mass = pc->AmIA<PhysicalSpecies::photon>() ? PhysConst::m_e : pc->getMass();
      RandomFilter const random_filter(part_diag.m_do_random_filter,
                                       part_diag.m_random_fraction);
      UniformFilter const uniform_filter(part_diag.m_do_uniform_filter,
                                         part_diag.m_uniform_stride);

      ParserFilter parser_filter(part_diag.m_do_parser_filter,
                                 utils::parser::compileParser<ParticleDiag::m_nvars>
                                   (part_diag.m_particle_filter_parser.get()),
                                 pc->getMass(), time);
      parser_filter.m_units = InputUnits::SI;
      GeometryFilter const geometry_filter(part_diag.m_do_geom_filter,
                                           part_diag.m_diag_domain);

      particlesConvertUnits(ConvertDirection::WarpX_to_SI, ptlSrc, mass);
      using SrcData = WarpXParticleContainer::ParticleTileType::ConstParticleTileDataType;
      tmp.copyParticles(*ptlSrc,
            [random_filter,uniform_filter,parser_filter,geometry_filter]
            AMREX_GPU_HOST_DEVICE
            (const SrcData& src, int ip, const amrex::RandomEngine& engine)
            {
                const SuperParticleType& p = src.getSuperParticle(ip);
                return random_filter(p, engine) * uniform_filter(p, engine)
                        * parser_filter(p, engine) * geometry_filter(p, engine);
            }, true);
        particlesConvertUnits(ConvertDirection::SI_to_WarpX, ptlSrc, mass);
  }
};
#endif // AMREX_USE_OPENPMD_API
#endif // WARPX_FLUSHFORMATPLOTPLUS_H_
