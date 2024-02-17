#ifndef _NONLINEAR_SOLVER_LIBRARY_H_
#define _NONLINEAR_SOLVER_LIBRARY_H_

#include "PicardSolver.H"
#include "NewtonSolver.H"

/**
  * \brief struct to select the nonlinear solver for implicit schemes
  */
enum NonlinearSolverType {
    Picard = 0,
    Newton = 1
};

#endif
