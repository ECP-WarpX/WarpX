/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_MSG_LOGGER_FWD_H
#define WARPX_MSG_LOGGER_FWD_H

#include <AMReX.H>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Utils{
namespace MsgLogger{

    enum class Priority;

    inline constexpr const char* PriorityToString(const Priority& priority);

    inline Priority StringToPriority(const std::string& priority);

    struct Msg;

    struct MsgWithCounter;

    struct MsgWithCounterAndRanks;

    class Logger;
}
}

#endif //WARPX_MSG_LOGGER_FWD_H
