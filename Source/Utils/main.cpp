#include "MsgLogger.H"

#include <iostream>
#include <sstream>

int main()
{
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
        "Algorithms", "algorithm xxx is experimental!");

    std::stringstream ss;
    logger.print_warnings(ss);

    std::cout << ss.str();
}