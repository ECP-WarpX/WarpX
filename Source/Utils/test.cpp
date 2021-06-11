#include "MsgLogger.H"

#include <mpi.h>


#include <iostream>
#include <sstream>
#include <array>

const auto arr_1 = std::array<bool,4>{false, false, false, false};
const auto arr_2 = std::array<bool,4>{true, false, false, false};
const auto arr_3 = std::array<bool,4>{false, true, false, false};
const auto arr_4 = std::array<bool,4>{true, true, true, true};

const int master_rank = 0;

int main(int argc, char** argv)
{
    MPI_Init(NULL, NULL);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MsgLogger::Logger logger;

    logger.record_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::low,
        "Particles", "species_1 has 0 particles at time 0");

    logger.record_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::medium,
        "Laser", "laser_1 has e_max=0 and will be disabled");

    logger.record_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::high,
        "Parallelization", "Parallelization is extremely inefficient!");

    logger.record_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::medium,
        "Laser", "laser_1 is out of the simulation box and will be disabled");

    logger.record_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::medium,
        "Algorithms", "algorithm AAA is experimental!");

    logger.record_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::medium,
        "MultiLineComment", "A\nmulti-line \ncomment!");

    logger.record_collective_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::low,
        "W1", "Collective warning 1. Affected tasks:",
        (world_rank < 4) ? arr_1[world_rank] : false,
        MPI_COMM_WORLD);

    logger.record_collective_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::low,
        "W2", "Collective warning 2.\nAffected tasks:",
        (world_rank < 4) ? arr_2[world_rank] : false,
        MPI_COMM_WORLD);

    logger.record_collective_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::low,
        "W3", "Collective warning 3.\nAffected tasks:",
        (world_rank < 4) ? arr_3[world_rank] : false,
        MPI_COMM_WORLD);

    logger.record_collective_entry(
        MsgLogger::Type::warning, MsgLogger::Importance::low,
        "W4", "Collective warning 4.\nAffected tasks:",
        (world_rank < 4) ? arr_4[world_rank] : false,
        MPI_COMM_WORLD);

    if (world_rank == master_rank){
        std::stringstream ss;
        logger.print_warnings(ss);
        std::cout << ss.str();
    }

    MPI_Finalize();
}