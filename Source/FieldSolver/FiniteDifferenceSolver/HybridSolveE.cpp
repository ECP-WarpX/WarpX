#include "FiniteDifferenceSolver.H"

using namespace amrex;

void FiniteDifferenceSolver::HybridSolveE (
    std::array< std::unique_ptr<amrex::MultiFab>, 3>& Efield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3> const& Bfield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& Jfield,
    std::unique_ptr<amrex::MultiFab> const& rhofield,
    std::array< std::unique_ptr<amrex::MultiFab>, 3 > const& edge_lengths )
{
#ifndef AMREX_USE_EB
    amrex::ignore_unused(edge_lengths);
#endif

    Print() << "kom tog hier in..." << std::endl;
}