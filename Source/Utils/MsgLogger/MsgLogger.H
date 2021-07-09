#ifndef WARPX_MSG_LOGGER_H_
#define WARPX_MSG_LOGGER_H_

#include <map>
#include <string>

namespace Utils{
namespace MsgLogger{

    enum class Priority
    {
        low,
        medium,
        high
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

    constexpr bool operator<(const Msg& l, const Msg& r)
    {
        return  (l.priority < r.priority) ||
        (l.topic < r.topic) ||
        (l.text < r.text);
    }

    class Logger
    {

        public:

        Logger();

        void record_msg(Msg msg);

        private:

        std::map<Msg, int> m_messages;
    };
}
}

#endif //WARPX_MSG_LOGGER_H_
