#ifndef WARPX_WARN_MANAGER_H_
#define WARPX_WARN_MANAGER_H_

#include "WarnManager_fwd.H"

#include "MsgLogger/MsgLogger.H"

#include <AMReX_ParmParse.H>

#include <string>
#include <vector>

namespace Utils
{
    class WarnManager
    {

        public:

        WarnManager();

        void record_warning(
            std::string topic,
            std::string text,
            MsgLogger::Priority priority);

        std::string print_local_warnings(
            const std::string& when) const;

        std::string print_global_warnings(
            const std::string& when) const;

        void debug_read_warnings_from_input(amrex::ParmParse& params);

        static const int warn_line_size = 80;
        static const int warn_tab_size = 5;

        private:

        std::vector<MsgLogger::MsgWithCounter>
        aux_sort_messages(
            std::vector<MsgLogger::MsgWithCounter>&& all_msg_with_counter) const;

        std::vector<MsgLogger::MsgWithCounterAndRanks>
        aux_sort_messages(
            std::vector<MsgLogger::MsgWithCounterAndRanks>&& all_msg_with_counter_and_ranks) const;

        std::string aux_get_header(
            const std::string& when,
            int line_size,
            bool is_global) const;

        std::string aux_print_warn_msg(
            const MsgLogger::MsgWithCounter& msg_with_counter) const;

        std::string aux_print_warn_msg(
            const MsgLogger::MsgWithCounterAndRanks& msg_with_counter_and_ranks) const;

        std::string aux_msg_formatter(
            const std::string& msg,
            const int line_size,
            const int tab_size) const;

        int m_rank = 0;
        MsgLogger::Logger m_logger;
    };
}

#endif //WARPX_WARN_MANAGER_H_
