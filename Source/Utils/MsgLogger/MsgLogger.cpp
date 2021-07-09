
#include "MsgLogger.H"

#include <iostream> //DEBUG

using namespace Utils::MsgLogger;

Logger::Logger(){}

void Logger::record_msg(Msg msg)
{
    m_messages[msg]++;
    std::cout << "DEBUG :" << "\n" <<
        "    topic   : " << msg.topic << "\n" <<
        "    text    : " << msg.text << "\n" <<
        "    priority: " << PriorityToString(msg.priority) << "\n" <<
        "    count   : " << m_messages[msg] << "\n" <<
        "____________________________________\n";
}
