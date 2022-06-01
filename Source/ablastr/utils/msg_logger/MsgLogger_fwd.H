/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef ABLASTR_MSG_LOGGER_FWD_H
#define ABLASTR_MSG_LOGGER_FWD_H

namespace ablastr::utils::msg_logger
{
    enum class Priority;

    struct Msg;
    struct MsgWithCounter;
    struct MsgWithCounterAndRanks;

    class Logger;
}

#endif //ABLASTR_MSG_LOGGER_FWD_H
