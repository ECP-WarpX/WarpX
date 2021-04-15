/* Copyright 2019 David Grote, Luca Fedeli, Remi Lehe
 * Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef UTILS_WARPXALGORITHMSELECTION_H_
#define UTILS_WARPXALGORITHMSELECTION_H_

#include <AMReX_ParmParse.H>
#include <string>

/**
  * \brief struct to determine the computational medium, i.e., vacuum or material/macroscopic
           default is vacuum.
  */
struct MediumForEM {
    enum {
        Vacuum = 0,
        Macroscopic = 1
    };
};

/**
  * \brief struct to select algorithm for macroscopic Maxwell solver
           LaxWendroff (semi-implicit) represents sigma*E = sigma*0.5*(E^(n) + E^(n+1))
           Backward Euler (fully-implicit) represents sigma*E = sigma*E^(n+1)
           default is Backward Euler as it is more robust.
  */
struct MacroscopicSolverAlgo {
    enum {
        BackwardEuler = 0,
        LaxWendroff = 1
    };
};

struct MaxwellSolverAlgo {
    enum {
        Yee = 0,
        CKC = 1,
        PSATD = 2
    };
};

struct ElectrostaticSolverAlgo {
    enum {
        None = 0,
        Relativistic = 1,
        LabFrame = 2
    };
};

struct ParticlePusherAlgo {
    enum {
        Boris = 0,
        Vay = 1,
        HigueraCary = 2
    };
};

struct CurrentDepositionAlgo {
    enum {
         Esirkepov = 0,
         Direct = 1,
         Vay = 2
    };
};

struct ChargeDepositionAlgo {
    // Only the Standard algorithm is implemented
    enum {
         Standard = 0
    };
};

struct GatheringAlgo {
    enum {
         EnergyConserving = 0,
         MomentumConserving
    };
};

/** Strategy to compute weights for use in load balance.
 */
struct LoadBalanceCostsUpdateAlgo {
    enum {
        Timers    = 0, //!< load balance according to in-code timer-based weights (i.e., with  `costs`)
        Heuristic = 1, /**< load balance according to weights computed from number of cells
                             and number of particles per box (i.e., with `costs_heuristic`)*/
        GpuClock  = 2
    };
};

/** Field boundary conditions at the domain boundary
 */
struct FieldBoundaryType {
    enum {
        PML = 0,
        Periodic = 1,
        PEC = 2,     //!< perfect electric conductor (PEC) with E_tangential=0
        PMC = 3      //!< perfect magnetic conductor (PMC) with B_tangential=0
    };
};

/** Particle boundary conditions at the domain boundary
 */
struct ParticleBoundaryType {
    enum {
        Absorbing = 0,     //!< particles crossing domain boundary are removed
        Open = 1,          //!< particles cross domain boundary leave with damped j
        Reflecting = 2,     //!< particles are reflected
        Periodic = 3
    };
};

int
GetAlgorithmInteger( amrex::ParmParse& pp, const char* pp_search_key );

/** Select BC Type for fields, if field=true
 *  else select BCType for particles.
 */
int
GetBCTypeInteger( std::string BCType, bool field );

#endif // UTILS_WARPXALGORITHMSELECTION_H_
