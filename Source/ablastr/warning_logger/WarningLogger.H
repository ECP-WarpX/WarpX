/* Copyright 2022 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef ABLASTR_WARNING_LOGGER_H_
#define ABLASTR_WARNING_LOGGER_H_

#include <string>

namespace ablastr::warning_logger
{
    /** Priority is recorded together with warnings. It influences
    * the display order and the appearance of a warning message.
    */
    enum class WarnPriority
    {
        /** Low priority warning message */
        low,
        /** Medium priority warning message */
        medium,
        /** High priority warning message */
        high
    };

    struct WarnParams{
        bool record_warnings = true;
        bool always_warn_immediately = false;
        bool abort_on_warning = false;
        WarnPriority abort_on_warning_threshold = WarnPriority::high;
    };

    void InitWarningLogger(WarnParams warn_params = WarnParams{});

    bool IsWarningLoggerInit();

    WarnParams GetWarningLoggerParams();

    void SetWarningLoggerParams(WarnParams warn_params);

    /**
    * \brief This function records a warning message.
    *
    * @param[in] topic a string to identify the topic of the warning (e.g., "parallelization", "pbc", "particles"...)
    * @param[in] text the text of the warning message
    * @param[in] priority priority of the warning message ("medium" by default)
    */
    void RecordWarning(
        std::string topic,
        std::string text,
        WarnPriority priority);

    /**
    * \brief This function prints all the warning messages collected on the present MPI rank
    * (i.e., this is not a collective call). This function is mainly intended for debug purposes.
    *
    * @param[in] when a string to mark when the warnings are printed out (it appears in the warning list)
    * @return a string containing the "local" warning list
    */
    std::string PrintLocalWarnings(
            const std::string& when);

    /**
    * \brief This function prints all the warning messages collected by all the MPI ranks
    * (i.e., this is a collective call). Only the I/O rank prints the message.
    *
    * @param[in] when a string to mark when the warnings are printed out (it appears in the warning list)
    * @return a string containing the "global" warning list
    */
    std::string PrintGlobalWarnings(
        const std::string& when);

    /**
    * \brief This function reads warning messages from the inputfile. It is intended for
    * debug&testing purposes
    *
    * @param[in, out] params the inputfile parser
    */
    void DebugReadWarningsFromInput(amrex::ParmParse& params);
}

#endif //ABLASTR_WARNING_LOGGER_H_
