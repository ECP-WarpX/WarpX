#include "ElectrostaticSolver.H"

class LabFrameExplicitES : public ElectrostaticBase
{
    public:
   
    void ComputeSpaceChargeField (amrex::Vector< std::array< std::unique_ptr<amrex::MultiFab>, 3> > Efield_fp, MultiParticleContainer& mpc, MultiFluidContainer mfl );
    void computePhiTriDiagonal ( const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
    amrex::Vector<std::unique_ptr<amrex::MultiFab> >& phi);
    // Multifab for electrostatic potential
    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_phi_fp;
    amrex::Vector< std::unique_ptr<amrex::MultiFab> > m_rho_fp;
}
