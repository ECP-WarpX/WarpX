#ifndef WARPX_MSG_LOGGER_H_
#define WARPX_MSG_LOGGER_H_

#include <AMReX.H>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Utils{
namespace MsgLogger{

    enum class Priority
    {
        high,
        medium,
        low
    };

    inline constexpr const char* PriorityToString(const Priority& priority)
    {
        if(priority == Priority::high)
            return "high";
        else if (priority == Priority::medium)
            return "medium";
        else
            return "low";
    }

    inline Priority StringToPriority(const std::string& priority)
    {
        if(priority == "high")
            return Priority::high;
        else if (priority == "medium")
            return Priority::medium;
        else if (priority == "low")
            return Priority::low;
        else
            amrex::Abort("Priority string not recognized");
        return Priority::low;
    }

    struct Msg
    {
        std::string topic;
        std::string text;
        Priority priority;

        std::vector<char> serialize() const;
        static Msg deserialize(std::vector<char>::const_iterator& it);
    };

    struct MsgWithCounter
    {
        Msg msg;
        int counter;

        std::vector<char> serialize() const;
        static MsgWithCounter deserialize(std::vector<char>::const_iterator& it);
    };

    struct MsgWithCounterAndRanks
    {
        MsgWithCounter msg_with_counter;
        bool all_ranks;
        std::vector<int> ranks;

        std::vector<char> serialize() const;
        static MsgWithCounterAndRanks deserialize(std::vector<char>::const_iterator& it);
    };

    constexpr bool operator<(const Msg& l, const Msg& r)
    {
        return (l.priority < r.priority) ||
        ((l.priority == r.priority) && (l.topic < r.topic)) ||
        ((l.priority == r.priority) && (l.topic == r.topic) && (l.text < r.text));
    }

    class Logger
    {

        public:

        Logger();

        void record_msg(Msg msg);

        std::vector<Msg> get_msg_list() const;

        std::vector<MsgWithCounter> get_msg_with_counter_list() const;

        std::vector<MsgWithCounterAndRanks>
        collective_gather_msg_with_counter_and_ranks() const;

        private:

        std::vector<char> serialize_msg_list(
            const std::vector<Msg>& msg_list) const;
        static std::vector<Msg> deserialize_msg_list(
            const std::vector<char>& serialized);

        std::pair<int,int> aux_find_gather_rank_and_items(
            int how_many_items) const;

        std::vector<MsgWithCounterAndRanks>
        aux_get_msg_with_counter_and_ranks() const;

        std::vector<char>
        collective_gather_get_serialized_gather_rank_msg_list(
            const std::vector<Msg>& my_list,
            int gather_rank) const;

        std::vector<char>
        aux_create_data_package(
            const std::vector<int>& gather_rank_msg_counters,
            const std::map<Msg, int>& to_send
        ) const;

        void
        update_list_with_packaged_data(
            std::vector<MsgWithCounterAndRanks>& list_with_ranks,
            const std::vector<char>& all_data,
            const std::vector<int>& disps) const;

        int m_rank = 0;
        int m_num_procs = 0;
        int m_io_rank = 0;
        bool m_am_i_io = false;
        std::map<Msg, int> m_messages;
    };
}
}

#endif //WARPX_MSG_LOGGER_H_
