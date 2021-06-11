/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MSG_LOGGER_H_
#define WARPX_MSG_LOGGER_H_

#include <mpi.h> //for tests

#include <vector>
#include <map>
#include <string>
#include <sstream>

namespace MsgLogger{

enum class Type{
    warning
};

enum class Importance{
    low,
    medium,
    high
};

struct WarnStyle{
    std::string header;
    std::string low_prefix;
    std::string medium_prefix;
    std::string high_prefix;
    std::string new_line_prefix;
    std::string topic_left;
    std::string topic_right;
    std::string footer;
    std::string no_warning_msg;
};

const auto default_warn_style =
    WarnStyle{
        "******************* WARNINGS! *********************",
        "* [!  ]", "* [!! ]", "* [!!!]",
        "*       ",
        " [ ", " ] ",
        "***************************************************",
        "* No warnings have been raised!"};

class Logger
{

public:
    Logger();

    void record_entry(
        Type type,
        Importance importance,
        std::string topic,
        std::string text);

    void record_collective_entry(
        Type type,
        Importance importance,
        std::string topic,
        std::string text,
        bool affects_me,
        MPI_Comm comm,
        int ranks_per_line = 15);

    void print_warnings(std::stringstream& ss,
        const WarnStyle& ws = default_warn_style);

private:

    void
    aux_print_entries(
        const std::string& prefix,
        const std::string& topic,
        const std::vector<std::string>& entries,
        const WarnStyle& ws,
        std::stringstream& ss) const;

    std::string
    aux_msg_formatter(
        const std::string& msg,
        const int first_line_offset,
        const WarnStyle& ws) const;

    std::map<Type,
        std::map<Importance,
            std::map<std::string,
                std::vector<std::string>>>> m_entries;
};

}

#endif //WARPX_MSG_LOGGER_H_
