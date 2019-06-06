#ifndef WARPX_PhotonParticleContainer_H_
#define WARPX_PhotonParticleContainer_H_

#include <PhysicalParticleContainer.H>
#include <AMReX_Vector.H>
#include <memory>

#ifdef WARPX_QED
    #include "breit_wheeler_engine_wrapper.h"
#endif

class PhotonParticleContainer
    : public PhysicalParticleContainer
{
public:

#ifdef WARPX_QED
        //A copy of the BW engine should be passed to the constructor
    PhotonParticleContainer (amrex::AmrCore* amr_core,
                        int ispecies,
                        const std::string& name,
                        std::shared_ptr<warpx_breit_wheeler_engine> bw_engine);
#else
    PhotonParticleContainer (amrex::AmrCore* amr_core,
                                    int ispecies,
                                    const std::string& name);
#endif

    virtual ~PhotonParticleContainer () {}

    virtual void InitData() override;

    virtual void Evolve (int lev,
                         const amrex::MultiFab& Ex,
                         const amrex::MultiFab& Ey,
                         const amrex::MultiFab& Ez,
                         const amrex::MultiFab& Bx,
                         const amrex::MultiFab& By,
                         const amrex::MultiFab& Bz,
                         amrex::MultiFab& jx,
                         amrex::MultiFab& jy,
                         amrex::MultiFab& jz,
                         amrex::MultiFab* cjx,
                         amrex::MultiFab* cjy,
                         amrex::MultiFab* cjz,
                         amrex::MultiFab* rho,
                         amrex::MultiFab* crho,
                         const amrex::MultiFab* cEx,
                         const amrex::MultiFab* cEy,
                         const amrex::MultiFab* cEz,
                         const amrex::MultiFab* cBx,
                         const amrex::MultiFab* cBy,
                         const amrex::MultiFab* cBz,
                         amrex::Real t,
                         amrex::Real dt) override;

    virtual void PushPX(WarpXParIter& pti,
                        amrex::Cuda::ManagedDeviceVector<amrex::Real>& xp,
                        amrex::Cuda::ManagedDeviceVector<amrex::Real>& yp,
                        amrex::Cuda::ManagedDeviceVector<amrex::Real>& zp,
                        amrex::Cuda::ManagedDeviceVector<amrex::Real>& giv,
                        amrex::Real dt) override;

  virtual void DepositCurrent(WarpXParIter& pti,
                              RealVector& wp,
                              RealVector& uxp,
                              RealVector& uyp,
                              RealVector& uzp,
                              amrex::MultiFab& jx,
                              amrex::MultiFab& jy,
                              amrex::MultiFab& jz,
                              amrex::MultiFab* cjx,
                              amrex::MultiFab* cjy,
                              amrex::MultiFab* cjz,
                              const long np_current,
                              const long np,
                              int thread_num,
                              int lev,
                              amrex::Real dt ) override {};

//This function initialises the optical depth of the photons
#ifdef WARPX_QED
    virtual void InitOpticalDepth(
                                    WarpXParIter& pti,
                                    warpx_breit_wheeler_engine& engine);
#endif

private:

//Needs a copy of the BW engine
#ifdef WARPX_QED
    std::shared_ptr<warpx_breit_wheeler_engine> bw_engine;
#endif

    amrex::Real size_in_inches = 1.;

};

#endif // #ifndef WARPX_PhotonParticleContainer_H_
