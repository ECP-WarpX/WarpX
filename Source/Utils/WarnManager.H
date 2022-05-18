/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_WARN_MANAGER_H_
#define WARPX_WARN_MANAGER_H_

#include "WarnManager_fwd.H"

#include "MsgLogger/MsgLogger_fwd.H"

#include <AMReX_ParmParse.H>

#include <memory>
#include <string>
#include <vector>

namespace Utils
{
    /**
    * The class WarnManager manages warning messages in WarpX,
    * providing methods to record warnings, and print warning
    * lists.
    */
    class WarnManager
    {
        public:

        /**
        * The constructor.
        */
        WarnManager();

        /**
        * \brief This function records a warning message.
        *
        * @param[in] topic a string to identify the topic of the warning (e.g., "parallelization", "pbc", "particles"...)
        * @param[in] text the text of the warning message
        * @param[in] priority priority of the warning message ("medium" by default)
         */
        void record_warning(
            std::string topic,
            std::string text,
            MsgLogger::Priority priority);

        /**
        * \brief This function prints all the warning messages collected on the present MPI rank
        * (i.e., this is not a collective call). This function is mainly intended for debug purposes.
        *
        * @param[in] when a string to mark when the warnings are printed out (it appears in the warning list)
        * @return a string containing the "local" warning list
        */
        std::string print_local_warnings(
            const std::string& when) const;

        /**
        * \brief This function prints all the warning messages collected by all the MPI ranks
        * (i.e., this is a collective call). Only the I/O rank prints the message.
        *
        * @param[in] when a string to mark when the warnings are printed out (it appears in the warning list)
        * @return a string containing the "global" warning list
        */
        std::string print_global_warnings(
            const std::string& when) const;

        /**
        * \brief This function reads warning messages from the inputfile. It is intended for
        * debug&testing purposes
        *
        * @param[in, out] params the inputfile parser
        */
        void debug_read_warnings_from_input(amrex::ParmParse& params);

        static const int warn_line_size = 80 /*! Maximum line length to be used in formatting warning list*/;
        static const int warn_tab_size = 5 /*! Tabulation size to be used in formatting warning list*/;

        private:

        /**
        * \brief This function generates a string for a single entry of the warning list
        * for a MessageWithCounter struct (i.e., a warning message paired with a counter storing
        * how many times the warning has been raised)
        *
        * @param[in] msg_with_counter a MessageWithCounter
        * @return a string containing the warning message
        */
        std::string print_warn_msg(
            const MsgLogger::MsgWithCounter& msg_with_counter) const;

        /**
        * \brief This function generates a string for a single entry of the warning list
        * for a MsgWithCounterAndRanks struct (i.e., a warning message paired with a counter storing
        * how many times the warning has been raised and info on which ranks have raised the warning)
        *
        * @param[in] msg_with_counter_and_ranks a MsgWithCounterAndRanks
        * @return a string containing the warning message
        */
        std::string print_warn_msg(
            const MsgLogger::MsgWithCounterAndRanks& msg_with_counter_and_ranks) const;

        /**
        * \brief This function generates the header of the warning messages list
        *
        * @param[in] when a string to mark when the warnings are printed out (it appears in the warning list)
        * @param[in] line_size maximum line length to be used in formatting warning list
        * @param[in] is_global flag: true if the header is for a global warning list, false otherwise
        * @return a string containing the header of the warning list
        */
        static std::string get_header(
            const std::string& when,
            const int line_size,
            const bool is_global);

        /**
        * \brief This function formats each line of a warning message text
        *
        * @param[in] msg the warning message text
        * @param[in] line_size maximum line length to be used in formatting warning list
        * @param[in] tab_size tabulation size to be used in formatting warning list
        * @return a string containing the formatted warning message text
        */
        static std::string msg_formatter(
            const std::string& msg,
            const int line_size,
            const int tab_size);

        int m_rank = 0 /*! MPI rank (appears in the warning list)*/;
        std::unique_ptr<MsgLogger::Logger> m_p_logger /*! The Logger stores all the warning messages*/;
    };
}

#endif //WARPX_WARN_MANAGER_H_
