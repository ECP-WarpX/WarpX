#ifndef WARPX_WARN_MANAGER_H_
#define WARPX_WARN_MANAGER_H_

#include "WarnManager_fwd.H"

#include "MsgLogger/MsgLogger.H"

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
            const std::string& when);

        static const int warn_line_size = 80;
        static const int warn_tab_size = 5;

        private:

        std::vector<MsgLogger::MsgWithCounter>
        aux_sort_messages(
            std::vector<MsgLogger::MsgWithCounter>&& all_msg_with_counter);

        std::string aux_get_header(
            const std::string& when,
            int line_size,
            bool is_global);

        std::string aux_print_warn_msg(
            const MsgLogger::MsgWithCounter& msg_with_counter);

        std::string aux_print_line(
            const std::string& str, int line_size, int tab_size);

        MsgLogger::Logger m_logger;
    };
}

#endif //WARPX_WARN_MANAGER_H_
