#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"
#include <AMReX_AmrParticles.H>

#include "warpxWriter.H"
using namespace amrex;


////////////////////////////
//class AMReX_warpxBTDWriter final : public amrex::openpmd_api::AMReX_openPMDWriter
class AMReX_warpxBTDWriter final : public AMReX_warpxWriter
{
 public:
  AMReX_warpxBTDWriter(std::vector<bool> fieldPMLdirections,
               openPMD::IterationEncoding ie =  openPMD::IterationEncoding::groupBased,
               std::string filetype = "default",
               std::string openpmdOptions = "{}");

  //~AMReX_warpxBTDWriter() override = default;
  ~AMReX_warpxBTDWriter() override;

  // BTD does not support streams
  openPMD::Iteration GetIteration (int const iteration) const override
    {
      return m_Series->iterations[iteration];
    }

  unsigned long long GetGrandOffset() const override
  {
    return m_AssignedPtlOffset;
  }

  /*
  */
  bool AllocatePtlProperties(openPMD::ParticleSpecies& currSpecies,
                 const amrex::Vector<int>& write_real_comp,
                 const amrex::Vector<std::string>& real_comp_names,
                 const amrex::Vector<int>& write_int_comp,
                 const amrex::Vector<std::string>& int_comp_names,
                 const unsigned long long np) const override;


  void SetupConstant(openPMD::ParticleSpecies& currSpecies,
             const unsigned long long& np) const override;

  void CloseStep(int ts) override;

  void Init(openPMD::Access access) override;

  //
  // in warpx, btd ignores "geom" and takes in a single geometry "full_BTD_snapshot"
  // and loops over all levels in the vector "geom"
  // might as well supply equivalent vector of full_BTD_snapshot instead of both "geom" and "full_BTD_snapshot"
  //
  //   TODO  -- need to find out whether BTD is only for levels if all levels are using full_BTD_snapshot in warpX
  //
  void WriteMesh(const std::vector<std::string>& varnames,
         const amrex::Vector<const amrex::MultiFab*>& mf,
         const amrex::Vector<amrex::Geometry>& geom,
         //const int iteration,
         const double time ) const override;
  // BTD specific:
  void SetLastFlush() {m_LastFlush = true;}
  bool IsLastFlush() {return m_LastFlush;}
  void AssignPtlOffset(unsigned long long m) { m_AssignedPtlOffset = m; }

  void SetConstantMassCharge(openPMD::ParticleSpecies& currSpecies,
                 const unsigned long long& np,
                 amrex::ParticleReal const charge,
                 amrex::ParticleReal const mass) const override;

 private:
  bool m_LastFlush = false;
  unsigned long long m_AssignedPtlOffset = 0;
  bool m_FirstSightOfStep = false;
}; // class warpxBTDwriter




AMReX_warpxBTDWriter::AMReX_warpxBTDWriter(
                       std::vector<bool> fieldPMLdirections,
                       openPMD::IterationEncoding ie,
                       std::string filetype,
                       std::string openpmdOptions)
  :AMReX_warpxWriter(fieldPMLdirections, ie, filetype, openpmdOptions)
{
  m_openPMDDatasetOptions="{ \"resizable\": true }";

  if (ie == openPMD::IterationEncoding::variableBased)
    {
      std::string warnMsg = " Unable to support BTD with streaming. Using GroupBased ";
      m_openPMDEncoding = openPMD::IterationEncoding::groupBased;
    }
}

AMReX_warpxBTDWriter::~AMReX_warpxBTDWriter()
{}


void AMReX_warpxBTDWriter::CloseStep(int ts)
{
  if  (!m_LastFlush)
    return;

  AMReX_openPMDWriter::CloseStep(ts);
  m_LastFlush = false; // this step is over
}

void AMReX_warpxBTDWriter::Init(openPMD::Access access)
{
  if (m_Series != nullptr) {
    m_FirstSightOfStep = ! m_Series->iterations.contains( m_CurrentStep );
    return;
  }
  m_FirstSightOfStep = true;
  AMReX_openPMDWriter::Init(access);
}

bool AMReX_warpxBTDWriter::AllocatePtlProperties(openPMD::ParticleSpecies& currSpecies,
                          const amrex::Vector<int>& write_real_comp,
                          const amrex::Vector<std::string>& real_comp_names,
                          const amrex::Vector<int>& write_int_comp,
                          const amrex::Vector<std::string>& int_comp_names,
                          const unsigned long long np) const
{
  if (np > 0)
    {
      return AMReX_openPMDWriter::AllocatePtlProperties(currSpecies, write_real_comp, real_comp_names,
                         write_int_comp, int_comp_names, np + GetGrandOffset());
    }
  //  np == 0
  if ( m_LastFlush && ( m_AssignedPtlOffset == 0 ) )
    {
      // properties were never allocated
      //AMReX_openPMDWriter::SetupPos(currSpecies, np);
      return AMReX_openPMDWriter::AllocatePtlProperties(currSpecies, write_real_comp, real_comp_names, write_int_comp, int_comp_names, np);
    }

  if (m_LastFlush) {
      return true; // no need to do more allocation
  }

  return false;
}

/*
*/
void AMReX_warpxBTDWriter::SetupConstant(openPMD::ParticleSpecies& currSpecies,
                     const unsigned long long& np) const
{
  if (m_LastFlush)
    {
      AMReX_openPMDWriter::SetupConstant(currSpecies, np + GetGrandOffset());
    }
}

void AMReX_warpxBTDWriter::SetConstantMassCharge(openPMD::ParticleSpecies& currSpecies,
                         const unsigned long long& np,
                         amrex::ParticleReal const charge,
                         amrex::ParticleReal const mass) const
{
  if (m_LastFlush)
    {
      BL_PROFILE("SetChargeMass(BTD)");
      AMReX_warpxWriter::SetConstantMassCharge(currSpecies,  np + GetGrandOffset(),
                        charge, mass);
    }
}



void AMReX_warpxBTDWriter::WriteMesh(const std::vector<std::string>& varnames,
                     const amrex::Vector<const amrex::MultiFab*>& mf,
                     const amrex::Vector<amrex::Geometry>& geom,
                     //const int iteration,
                     const double time ) const
{
  BL_PROFILE("WriteMesh(BTD)");
  if (m_FirstSightOfStep)  {
    AMReX_openPMDWriter::WriteMesh(varnames, mf, geom, /*iteration,*/ time);
    return;
  }

  // field and meshcomp were already setup for this iteration
  // only need to
  openPMD::Iteration series_iteration = GetIteration(m_CurrentStep);
  series_iteration.open();

  auto meshes = series_iteration.meshes;
  for (int lev=0; lev < geom.size(); lev++)
    {
      amrex::Geometry full_geom = geom[lev];
      //CompSetup(lev, meshes, full_geom, varnames, mf[lev]);
      CompStorage(lev, meshes, full_geom, varnames, mf[lev]);
#ifdef AMREX_USE_GPU
      amrex::Gpu::streamSynchronize();
#endif
      m_Series->flush();
    } // for lev loop
}
