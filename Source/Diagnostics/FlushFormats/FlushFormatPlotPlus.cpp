#include "FlushFormatPlotPlus.H"

#include "Particles/ParticleIO.H"
#include "Diagnostics/ParticleDiag/ParticleDiag.H"
#include "Particles/Filter/FilterFunctors.H"
#include "Particles/WarpXParticleContainer.H"
#include "Particles/PinnedMemoryParticleContainer.H"
#include "Utils/Interpolate.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "WarpX.H"
#include "Diagnostics/MultiDiagnostics.H"

#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Config.H>
#include <AMReX_GpuAllocators.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_IntVect.H>
#include <AMReX_MakeType.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleIO.H>
#include <AMReX_Particles.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>
#include <AMReX_buildInfo.H>

#ifdef AMREX_USE_OMP
#   include <omp.h>
#endif

#include <algorithm>
#include <array>
#include <cstring>
#include <fstream>
#include <map>
#include <memory>
#include <utility>
#include <vector>


#ifdef AMREX_USE_OPENPMD_API


#include "warpxWriter.H"
#include "warpxBTD.H"

class StepMgr
{
public:
  StepMgr(int step, AMReXWithOpenPMD* owner)
    :m_Step(step),
     m_Owner(owner)
  {
    m_Owner = owner;
    m_Owner->m_UserHandler->m_Writer->SetStep(m_Step);
  }
  ~StepMgr()
  {
    m_Owner->m_UserHandler->m_Writer->CloseStep(m_Step);
  }

private:
  int m_Step;
  AMReXWithOpenPMD* m_Owner;
};


using namespace amrex;


void
FlushFormatPlotPlus::WriteToFile (
    const amrex::Vector<std::string>& varnames,
    const amrex::Vector<amrex::MultiFab>& mf,
    amrex::Vector<amrex::Geometry>& geom,
    const amrex::Vector<int> iteration, const double time,
    const amrex::Vector<ParticleDiag>& particle_diags,
    [[maybe_unused]] int nlev,
    const std::string prefix,
    int file_min_digits,
    [[maybe_unused]] bool plot_raw_fields,
    [[maybe_unused]] bool plot_raw_fields_guards,
    const bool use_pinned_pc,
    bool isBTD,
    int snapshotID,
    [[maybe_unused]] int bufferID,
    [[maybe_unused]] int numBuffers,
    const amrex::Geometry& full_BTD_snapshot,
    bool isLastBTDFlush, const amrex::Vector<int>& totalParticlesFlushedAlready) const
{
    WARPX_PROFILE("FlushFormatPlotPlus::WriteToFile()");

    //const std::string& filename = amrex::Concatenate(prefix, iteration[0], file_min_digits);
    // move the name of plotfile here for classic plotfiles
    std::string plot_name = amrex::Concatenate(prefix, snapshotID, file_min_digits);
    plot_name = plot_name+"/buffer";
    const std::string& filename = amrex::Concatenate(plot_name, iteration[0], file_min_digits);

    { // amrex-bp version
      WARPX_PROFILE("FlushFormatPlotPlus_OpenPMDPlotFiles()");

      std::string opmdPrefix = prefix;
      if (isBTD)
      {
          amrex::Vector<amrex::Geometry> btd_geoms(geom.size());
          for (auto i=0; i<geom.size(); i++)
              btd_geoms[i] = full_BTD_snapshot;
          BTDWriter(opmdPrefix, snapshotID, varnames, mf, btd_geoms, time, particle_diags, use_pinned_pc, isLastBTDFlush, totalParticlesFlushedAlready);
      }
      else
        DefaultWriter(opmdPrefix, iteration[0], varnames, mf, geom, time, particle_diags, use_pinned_pc);

    } // end of amrex-bp version
}

void FlushFormatPlotPlus::BTDWriter(const std::string& prefix,
                                    int output_iteration,
                                    const amrex::Vector<std::string> varnames,
                                    const amrex::Vector<amrex::MultiFab>& mf,
                                    amrex::Vector<amrex::Geometry>& geom,
                                    const double time,
                                    const amrex::Vector<ParticleDiag>& particle_diags,
                                    const bool /*use_pinned_pc */,
                                    bool isLastBTDFlush,
                                    const amrex::Vector<int>& totalParticlesFlushedAlready
                                    ) const
{
    auto & warpx = WarpX::GetInstance();

    if ( m_Writer->InitLocalHandler(prefix) )
      {
         AMReX_warpxBTDWriter* testWriter = new AMReX_warpxBTDWriter(warpx.getPMLdirections(),
                                                                     openPMD::IterationEncoding::groupBased
                                                                    );
         m_Writer->SetWriter(testWriter);
         testWriter->EnforceEncoding();
      }

    AMReX_warpxBTDWriter* btdWriter =  (AMReX_warpxBTDWriter*) (m_Writer->m_UserHandler->m_Writer.get());
    if (isLastBTDFlush) {
      btdWriter  ->SetLastFlush();
    }

    StepMgr sm(output_iteration, m_Writer.get());

    // write mesh
    if (!mf.empty()) {
        m_Writer->StoreMesh(//parms.nlevs,
            amrex::GetVecOfConstPtrs(mf),
            varnames,
            geom,
            static_cast<Real>(time)
            );
    }

    // write particles
    for (unsigned whichDiag = 0, n = particle_diags.size(); whichDiag < n; ++whichDiag)
    {
      const ParticleDiag& part_diag = particle_diags[whichDiag];
      WarpXParticleContainer* pc = part_diag.getParticleContainer();

      PinnedMemoryParticleContainer* pinned_pc = part_diag.getPinnedParticleContainer();
      PinnedMemoryParticleContainer tmp;
      tmp = pinned_pc->make_alike<amrex::PinnedArenaAllocator>();

      CopyPtls(tmp, time, pc, pinned_pc, part_diag);

      tmp.CountParticles();

      btdWriter -> AssignPtlOffset(totalParticlesFlushedAlready[whichDiag]);

      Vector<std::string> real_names;
      Vector<std::string> int_names;
      Vector<int> int_flags;
      Vector<int> real_flags;
      GetNames(part_diag, real_names, int_names, int_flags, real_flags);

      m_Writer->m_UserHandler->m_Writer->DumpParticles(tmp,
                          part_diag.getSpeciesName(),
                          real_flags,
                          int_flags,
                          real_names,
                          int_names,
                    [=] ([[maybe_unused]] auto& ppc, openPMD::ParticleSpecies& currSpecies, unsigned long long localTotal)
                    {
                      btdWriter->SetConstantMassCharge(currSpecies, localTotal, pc->getCharge(),  pc->getMass());
                    },
                    [=] (auto& pti, openPMD::ParticleSpecies& currSpecies, unsigned long long offset)
                    {
                      btdWriter->SavePosId_RZ(pti, currSpecies, offset);
                    });

    }
}

void FlushFormatPlotPlus::GetNames(const ParticleDiag& part_diag,
                   Vector<std::string>& real_names,
                   Vector<std::string>& int_names,
                   Vector<int>& int_flags,
                   Vector<int>& real_flags
                   ) const
{
  WarpXParticleContainer* pc = part_diag.getParticleContainer();
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
  for (auto const& x : runtime_rnames) { real_names[x.second+PIdx::nattribs] = x.first; }

  // plot any "extra" fields by default
  real_flags = part_diag.m_plot_flags;
  real_flags.resize(pc->NumRealComps(), 1);

  // and the names
  int_names.resize(pc->NumIntComps());
  auto runtime_inames = pc->getParticleRuntimeiComps();
  for (auto const& x : runtime_inames) { int_names[x.second+0] = x.first; }

  // plot by default
  int_flags.resize(pc->NumIntComps(), 1);
}

void FlushFormatPlotPlus::DefaultWriter(const std::string& prefix,
                    int iteration,
                    const amrex::Vector<std::string> varnames,
                    const amrex::Vector<amrex::MultiFab>& mf,
                    amrex::Vector<amrex::Geometry>& geom,
                    const double time,
                    const amrex::Vector<ParticleDiag>& particle_diags,
                    const bool use_pinned_pc
                    ) const
{
    auto & warpx = WarpX::GetInstance();

    if (m_Writer->InitLocalHandler(prefix)) {
      AMReX_warpxWriter* testWriter = new AMReX_warpxWriter(warpx.getPMLdirections());
      m_Writer->SetWriter(testWriter);
    }

    AMReX_warpxWriter* warpxWriter =  (AMReX_warpxWriter*) (m_Writer->m_UserHandler->m_Writer.get());

    int output_iteration = iteration;
    StepMgr sm(output_iteration, m_Writer.get());

    // write mesh
    m_Writer->StoreMesh(//parms.nlevs,
            amrex::GetVecOfConstPtrs(mf),
            varnames,
            geom,
            static_cast<Real>(time)
            );
    // write particles
    for (auto& part_diag : particle_diags) {
      WarpXParticleContainer* pc = part_diag.getParticleContainer();

      PinnedMemoryParticleContainer* pinned_pc = part_diag.getPinnedParticleContainer();
      PinnedMemoryParticleContainer tmp;
      if (use_pinned_pc) {
          tmp = pinned_pc->make_alike<amrex::PinnedArenaAllocator>();
          CopyPtls(tmp, time, pc, pinned_pc, part_diag);
      } else {
          tmp = pc->make_alike<amrex::PinnedArenaAllocator>();
          CopyPtls(tmp, time, pc, pc, part_diag);
      }

    tmp.CountParticles();

    Vector<std::string> real_names;
    Vector<std::string> int_names;
    Vector<int> int_flags;
    Vector<int> real_flags;
    GetNames(part_diag, real_names, int_names, int_flags, real_flags);
    m_Writer->m_UserHandler->m_Writer->DumpParticles(tmp,
                          part_diag.getSpeciesName(),
                          real_flags,
                          int_flags,
                          real_names,
                          int_names,

                    [=] ([[maybe_unused]] auto& ppc, openPMD::ParticleSpecies& currSpecies, unsigned long long localTotal)
                    {
                      warpxWriter->SetConstantMassCharge(currSpecies, localTotal, pc->getCharge(),  pc->getMass());
                    },
                    [=] (auto& pti, openPMD::ParticleSpecies& currSpecies, unsigned long long offset)
                    {
                      warpxWriter->SavePosId_RZ(pti, currSpecies, offset);
                    });
    }
}




AMReXWithOpenPMD::AMReXWithOpenPMD(const std::string& prefix)
  :m_Prefix(prefix)
{
  // warpx has multiple diags, each should maintain its own handler
  m_UserHandler = openpmd_api::InitUserHandler(prefix);
}

void AMReXWithOpenPMD::SetWriter(amrex::openpmd_api::AMReX_openPMDWriter* w)
{
  BL_ASSERT ( m_UserHandler != nullptr );
  BL_ASSERT ( w != nullptr );

  m_UserHandler->SetWriter(w);
}

AMReXWithOpenPMD::~AMReXWithOpenPMD()
{
  openpmd_api::CloseUserHandler(m_UserHandler);
}

bool AMReXWithOpenPMD::InitLocalHandler(const std::string& prefix)
{
  if (m_Prefix.compare(prefix) == 0)
    return false;

  m_Prefix = prefix;
  m_UserHandler = openpmd_api::InitUserHandler(prefix);
  return true;
}

void AMReXWithOpenPMD::StoreMesh (const Vector<const MultiFab*> &mf,
                  const Vector<std::string> &varnames,
                  const Vector<Geometry> &geom,
                  Real time)
{
  if ((m_UserHandler == nullptr) || (m_UserHandler->m_Writer == nullptr))
    return;
  BL_ASSERT ( geom.size() == mf.size() );
  BL_ASSERT ( mf[0]->nComp() <= varnames.size() );

  m_UserHandler->m_Writer->WriteMesh(varnames,
                     mf,
                     geom,
                     time);
}

#endif //#ifdef AMREX_USE_OPENPMD_API
