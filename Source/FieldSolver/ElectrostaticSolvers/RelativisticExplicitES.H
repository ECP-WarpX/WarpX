
#ifndef WARPX_RELATIVISTICEXPLICITSOLVERES_H_
#define WARPX_RELATIVISTICEXPLICITSOLVERES_H_

#include "ElectrostaticBase.H"
#include "Particles/WarpXParticleContainer.H"


class RelativisticExplicitES final : public ElectrostaticBase
{
public:

    RelativisticExplicitES (int nlevs_max) : ElectrostaticBase (nlevs_max) { ReadParameters();}

    void ComputeSpaceChargeField (
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_fp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_cp,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& charge_buf,
        amrex::Vector< std::unique_ptr<amrex::MultiFab> >& phi_fp,
        MultiParticleContainer& mpc,
        MultiFluidContainer* mfl,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Efield_fp,
        amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Bfield_fp
    );

    void AddSpaceChargeField (
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& charge_buf,
        WarpXParticleContainer& pc,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield
    );

    void AddBoundaryField (
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield,
        amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield
    );
};

#endif // WARPX_RELATIVISTICEXPLICITSOLVERES_H_
