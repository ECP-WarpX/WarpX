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
        bool all_affected = true;
        std::vector<int> affected_ranks{};
    };

    /*std::vector<char> serialize_msg(
        const Msg& msg);

    Msg deserialize_msg(
        const std::vector<char>& msg);

    std::vector<char> serialize_msg_with_counter(
        const MsgWithCounter& msg_with_counter);

    MsgWithCounter deserialize_msg_with_counter(
        std::vector<char>& serialized_msg_with_counter);*/

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

        void parallel_merge_msg_lists();

        private:

        std::map<Msg, int> m_messages;
    };
}
}

#endif //WARPX_MSG_LOGGER_H_
