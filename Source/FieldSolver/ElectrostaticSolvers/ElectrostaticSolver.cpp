#include "ElectrostaticSolver.H"
#include "PoissonBoundaryHandler.H"

ElectrostaticSolver::ElectrostaticSolver ()
{
    m_poisson_boundary_handler = std::make_unique<PoissonBoundaryHandler>();
}
