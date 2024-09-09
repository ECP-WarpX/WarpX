
#ifndef Relativistic_Explicit_Solver_ES_H_
#define Relativistic_Explicit_Solver_ES_H_

#include "ElectrostaticSolver.H"
#include "Particles/WarpXParticleContainer.H"


class RelativisticExplicitES : public ElectrostaticBase 
{
public: 

    void
    RelativisticExplicitES::ComputeSpaceChargeField (amrex::Vector<std::unique_ptr<amrex::MultiFab> > charge_buf,
                                                    WarpXParticleContainer& pc, 
                                                    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield, 
                                                    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield);

    void
    RelativisticExplicitES::AddSpaceChargeField (amrex::Vector<std::unique_ptr<amrex::MultiFab> > charge_buf,
                                                WarpXParticleContainer& pc, 
                                                amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield, 
                                                amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield);
    void
    RelativisticExplicitES::AddBoundaryField (amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Efield, 
                                            amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3>>& Bfield);
}

#endif