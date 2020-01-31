#ifndef UTILS_WARPXALGORITHMSELECTION_H_
#define UTILS_WARPXALGORITHMSELECTION_H_

#include <AMReX_ParmParse.H>
#include <string>

struct MaxwellSolverAlgo {
    enum {
        Yee = 0,
        CKC = 1
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
         Direct = 1
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

// Implementation of equation (13) of (Vay et al, JCP 243, 2013):
// - use rho via Gauss law and discrete continuity equation (last two terms)
// - use rho via Gauss law (second last term) and J for all remaining terms
struct PsatdPushAlgo {
    enum {
	 Standard = 0,
	 Mixed = 1
    };
};

int
GetAlgorithmInteger( amrex::ParmParse& pp, const char* pp_search_key );

#endif // UTILS_WARPXALGORITHMSELECTION_H_
