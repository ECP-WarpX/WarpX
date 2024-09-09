#include "ElectrostaticSolver.H"
#include "PoissonBoundaryHandler.H"

ElectrostaticSolver::ElectrostaticSolver ()
{
    m_boundaryhandler = std::make_unique<PoissonBoundaryHandler>();
}
