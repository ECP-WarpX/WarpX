
#ifndef WARPX_RELATIVISTICEXPLICITSOLVERES_H_
#define WARPX_RELATIVISTICEXPLICITSOLVERES_H_

#include "ElectrostaticSolver.H"
#include "Particles/WarpXParticleContainer.H"


class RelativisticExplicitES final : public ElectrostaticSolver
{
public:

    RelativisticExplicitES (int nlevs_max) : ElectrostaticSolver (nlevs_max) {
        ReadParameters();
    }

    void InitData () override;

    void ComputeSpaceChargeField (
        [[maybe_unused]] amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_fp,
        [[maybe_unused]] amrex::Vector< std::unique_ptr<amrex::MultiFab> >& rho_cp,
        [[maybe_unused]] amrex::Vector< std::unique_ptr<amrex::MultiFab> >& charge_buf,
        [[maybe_unused]] amrex::Vector< std::unique_ptr<amrex::MultiFab> >& phi_fp,
        [[maybe_unused]] MultiParticleContainer& mpc,
        [[maybe_unused]] MultiFluidContainer* mfl,
        [[maybe_unused]] amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Efield_fp,
        [[maybe_unused]] amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> >& Bfield_fp
    ) override;

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
