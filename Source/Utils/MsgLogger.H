/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MSG_LOGGER_H_
#define WARPX_MSG_LOGGER_H_

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

const auto new_line_skip = std::string{"*      "};

class Logger
{

public:
    Logger(){}

    void record_entry(
        Type type,
        Importance importance,
        std::string topic,
        std::string text);

    void print_warnings(std::stringstream& ss);

private:

    void
    aux_print_entries(
        const std::string& prefix,
        const std::string& topic,
        const std::vector<std::string>& entries,
        std::stringstream& ss);

    std::map<Type,
        std::map<Importance,
            std::map<std::string,
                std::vector<std::string>>>> m_entries;
};

}

#endif //WARPX_MSG_LOGGER_H_
