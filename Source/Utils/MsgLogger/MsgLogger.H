#ifndef WARPX_MSG_LOGGER_H_
#define WARPX_MSG_LOGGER_H_

#include <map>
#include <string>
#include <vector>

namespace Utils{
namespace MsgLogger{

    enum class Priority
    {
        high,
        medium,
        low
    };

    constexpr const char* PriorityToString(const Priority& priority)
    {
        if(priority == Priority::high)
            return "high";
        else if (priority == Priority::medium)
            return "medium";
        else
            return "low";
    }

    struct Msg
    {
        std::string topic;
        std::string text;
        Priority priority;
    };

    struct MsgWithCounter
    {
        Msg msg;
        int counter;

        std::vector<char> serialize() const;
        static MsgWithCounter deserialize(const std::vector<char>& serialized);
    };

    struct MsgWithCounterAndRanks
    {
        MsgWithCounter msg_with_counter;
        std::vector<int> ranks;
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

        std::vector<MsgWithCounter> get_msg_list() const;

        std::vector<MsgWithCounterAndRanks> collective_gather_msg_lists() const;

        private:

        std::vector<char> serialize_msg_list(
            const std::vector<MsgWithCounter>& msg_list) const;

        aux_broadcast_io_rank_list();

        std::map<Msg, int> m_messages;
    };
}
}

#endif //WARPX_MSG_LOGGER_H_
