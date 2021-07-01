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

    MsgLogger::Logger logger (world_rank == master_rank);

    logger.record_entry(
        MsgLogger::Importance::low,
        "Particles", "species_1 has 0 particles at time 0");

    logger.record_entry(
        MsgLogger::Importance::medium,
        "Laser", "laser_1 has e_max=0 and will be disabled");

    logger.collective_gather_entries("ONE");

    logger.record_entry(
        MsgLogger::Importance::high,
        "Parallelization", "Parallelization is extremely inefficient!");

    logger.record_entry(
        MsgLogger::Importance::medium,
        "Laser", "laser_1 is out of the simulation box and will be disabled");

    logger.record_entry(
        MsgLogger::Importance::medium,
        "Algorithms", "algorithm AAA is experimental!");

    logger.record_entry(
        MsgLogger::Importance::medium,
        "MultiLineComment", "A\nmulti-line \ncomment!");

    logger.collective_gather_entries("TWO");

    if ( (world_rank < 4) ? arr_1[world_rank] : false ){
        logger.record_entry(MsgLogger::Importance::low,
            "W1", "Collective warning 1");
    }

    if ( (world_rank < 4) ? arr_2[world_rank] : false ){
        logger.record_entry(MsgLogger::Importance::low,
            "W1", "Collective warning 2");
    }

    logger.collective_gather_entries("THREE");


    if ( (world_rank < 4) ? arr_3[world_rank] : false ){
        logger.record_entry(MsgLogger::Importance::low,
            "W2", "Collective warning 2");
    }

    if ( (world_rank < 4) ? arr_4[world_rank] : false ){
        logger.record_entry(MsgLogger::Importance::low,
            "W2", "Collective warning 4");
    }

    for (int i = 0; i < world_size; ++i){
        if (world_rank == i){
            std::stringstream ss;
            logger.debug_raw_print(ss);
            std::cout << "RANK " << i << std::endl;
            std::cout << ss.str() << std::endl;
            std::cout << "======= \n \n" << std::endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
}