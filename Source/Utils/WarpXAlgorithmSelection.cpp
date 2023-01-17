/* Copyright 2019-2020 Axel Huebl, David Grote, Luca Fedeli
 * Remi Lehe, Weiqun Zhang, Yinjian Zhao
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "WarpXAlgorithmSelection.H"

#include <AMReX.H>

#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cstring>
#include <cctype>
#include <map>
#include <utility>

// Define dictionary with correspondance between user-input strings,
// and corresponding integer for use inside the code

const std::map<std::string, int> electromagnetic_solver_algo_to_int = {
    {"none",    ElectromagneticSolverAlgo::None },
    {"yee",     ElectromagneticSolverAlgo::Yee },
    {"ckc",     ElectromagneticSolverAlgo::CKC },
    {"psatd",   ElectromagneticSolverAlgo::PSATD },
    {"ect",     ElectromagneticSolverAlgo::ECT },
    {"default", ElectromagneticSolverAlgo::Yee }
};

const std::map<std::string, int> electrostatic_solver_algo_to_int = {
    {"none", ElectrostaticSolverAlgo::None },
    {"relativistic", ElectrostaticSolverAlgo::Relativistic},
    {"labframe-electromagnetostatic", ElectrostaticSolverAlgo::LabFrameElectroMagnetostatic},
    {"labframe", ElectrostaticSolverAlgo::LabFrame},
    {"default", ElectrostaticSolverAlgo::None }
};

const std::map<std::string, int> particle_pusher_algo_to_int = {
    {"boris",   ParticlePusherAlgo::Boris },
    {"vay",     ParticlePusherAlgo::Vay },
    {"higuera", ParticlePusherAlgo::HigueraCary },
    {"default", ParticlePusherAlgo::Boris }
};

const std::map<std::string, int> current_deposition_algo_to_int = {
    {"esirkepov", CurrentDepositionAlgo::Esirkepov },
    {"direct",    CurrentDepositionAlgo::Direct },
    {"vay",       CurrentDepositionAlgo::Vay },
    {"default",   CurrentDepositionAlgo::Esirkepov } // NOTE: overwritten for PSATD below
};

const std::map<std::string, int> charge_deposition_algo_to_int = {
    {"standard",   ChargeDepositionAlgo::Standard },
    {"default",    ChargeDepositionAlgo::Standard }
};

const std::map<std::string, int> gathering_algo_to_int = {
    {"energy-conserving",   GatheringAlgo::EnergyConserving },
    {"momentum-conserving", GatheringAlgo::MomentumConserving },
    {"default",             GatheringAlgo::EnergyConserving }
};

const std::map<std::string, int> psatd_solution_type_to_int = {
    {"first-order", PSATDSolutionType::FirstOrder},
    {"second-order", PSATDSolutionType::SecondOrder},
    {"default", PSATDSolutionType::SecondOrder}
};

const std::map<std::string, int> J_in_time_to_int = {
    {"constant", JInTime::Constant},
    {"linear", JInTime::Linear},
    {"default", JInTime::Constant}
};

const std::map<std::string, int> rho_in_time_to_int = {
    {"constant", RhoInTime::Constant},
    {"linear", RhoInTime::Linear},
    {"default", RhoInTime::Linear}
};

const std::map<std::string, int> load_balance_costs_update_algo_to_int = {
    {"timers",    LoadBalanceCostsUpdateAlgo::Timers },
    {"gpuclock",  LoadBalanceCostsUpdateAlgo::GpuClock },
    {"heuristic", LoadBalanceCostsUpdateAlgo::Heuristic },
    {"default",   LoadBalanceCostsUpdateAlgo::Timers }
};

const std::map<std::string, int> MaxwellSolver_medium_algo_to_int = {
    {"vacuum", MediumForEM::Vacuum},
    {"macroscopic", MediumForEM::Macroscopic},
    {"default", MediumForEM::Vacuum}
};

const std::map<std::string, int> MacroscopicSolver_algo_to_int = {
    {"backwardeuler", MacroscopicSolverAlgo::BackwardEuler},
    {"laxwendroff", MacroscopicSolverAlgo::LaxWendroff},
    {"default", MacroscopicSolverAlgo::BackwardEuler}
};

const std::map<std::string, int> FieldBCType_algo_to_int = {
    {"pml",      FieldBoundaryType::PML},
    {"periodic", FieldBoundaryType::Periodic},
    {"pec",      FieldBoundaryType::PEC},
    {"pmc",      FieldBoundaryType::PMC},
    {"damped",   FieldBoundaryType::Damped},
    {"absorbing_silver_mueller", FieldBoundaryType::Absorbing_SilverMueller},
    {"neumann",  FieldBoundaryType::Neumann},
    {"none",     FieldBoundaryType::None},
    {"default",  FieldBoundaryType::PML}
};

const std::map<std::string, ParticleBoundaryType> ParticleBCType_algo_to_enum = {
    {"absorbing",  ParticleBoundaryType::Absorbing},
    {"open",       ParticleBoundaryType::Open},
    {"reflecting", ParticleBoundaryType::Reflecting},
    {"periodic",   ParticleBoundaryType::Periodic},
    {"default",    ParticleBoundaryType::Absorbing}
};

const std::map<std::string, int> ReductionType_algo_to_int = {
    {"maximum",  ReductionType::Maximum},
    {"minimum",  ReductionType::Minimum},
    {"integral", ReductionType::Sum}
};

int
GetAlgorithmInteger( amrex::ParmParse& pp, const char* pp_search_key ){

    // Read user input ; use "default" if it is not found
    std::string algo = "default";
    pp.query( pp_search_key, algo );
    // Convert to lower case
    std::transform(algo.begin(), algo.end(), algo.begin(), ::tolower);

    // Pick the right dictionary
    std::map<std::string, int> algo_to_int;
    if (0 == std::strcmp(pp_search_key, "maxwell_solver")) {
        algo_to_int = electromagnetic_solver_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "do_electrostatic")) {
        algo_to_int = electrostatic_solver_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "particle_pusher")) {
        algo_to_int = particle_pusher_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "current_deposition")) {
        algo_to_int = current_deposition_algo_to_int;
        if (WarpX::electromagnetic_solver_id == ElectromagneticSolverAlgo::PSATD)
            algo_to_int["default"] = CurrentDepositionAlgo::Direct;
    } else if (0 == std::strcmp(pp_search_key, "charge_deposition")) {
        algo_to_int = charge_deposition_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "field_gathering")) {
        algo_to_int = gathering_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "solution_type")) {
        algo_to_int = psatd_solution_type_to_int;
    } else if (0 == std::strcmp(pp_search_key, "J_in_time")) {
        algo_to_int = J_in_time_to_int;
    } else if (0 == std::strcmp(pp_search_key, "rho_in_time")) {
        algo_to_int = rho_in_time_to_int;
    } else if (0 == std::strcmp(pp_search_key, "load_balance_costs_update")) {
        algo_to_int = load_balance_costs_update_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "em_solver_medium")) {
        algo_to_int = MaxwellSolver_medium_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "macroscopic_sigma_method")) {
        algo_to_int = MacroscopicSolver_algo_to_int;
    } else if (0 == std::strcmp(pp_search_key, "reduction_type")) {
        algo_to_int = ReductionType_algo_to_int;
    } else {
        std::string pp_search_string = pp_search_key;
        amrex::Abort("Unknown algorithm type: " + pp_search_string);
    }

    // Check if the user-input is a valid key for the dictionary
    if (algo_to_int.count(algo) == 0) {
        // Not a valid key ; print error message
        std::string pp_search_string = pp_search_key;
        std::string error_message = "Invalid string for algo." + pp_search_string
            + ": " + algo + ".\nThe valid values are:\n";
        for ( const auto &valid_pair : algo_to_int ) {
            if (valid_pair.first != "default"){
                error_message += " - " + valid_pair.first + "\n";
            }
        }
        amrex::Abort(error_message);
    }

    // If the input is a valid key, return the value
    return algo_to_int[algo];
}

int
GetFieldBCTypeInteger( std::string BCType ){
    std::transform(BCType.begin(), BCType.end(), BCType.begin(), ::tolower);

    if (FieldBCType_algo_to_int.count(BCType) == 0) {
        std::string error_message = "Invalid string for field/particle BC. : " + BCType                         + "\nThe valid values are : \n";
        for (const auto &valid_pair : FieldBCType_algo_to_int) {
            if (valid_pair.first != "default"){
                error_message += " - " + valid_pair.first + "\n";
            }
        }
        amrex::Abort(error_message);
    }
    // return FieldBCType_algo_to_int[BCType]; // This operator cannot be used for a const map
    return FieldBCType_algo_to_int.at(BCType);
}

ParticleBoundaryType
GetParticleBCTypeInteger( std::string BCType ){
    std::transform(BCType.begin(), BCType.end(), BCType.begin(), ::tolower);

    if (ParticleBCType_algo_to_enum.count(BCType) == 0) {
        std::string error_message = "Invalid string for particle BC. : " + BCType + "\nThe valid values are : \n";
        for (const auto &valid_pair : ParticleBCType_algo_to_enum) {
            if (valid_pair.first != "default"){
                error_message += " - " + valid_pair.first + "\n";
            }
        }
        amrex::Abort(error_message);
    }
    // return ParticleBCType_algo_to_enum[BCType]; // This operator cannot be used for a const map
    return ParticleBCType_algo_to_enum.at(BCType);
}
