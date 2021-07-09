
#include "MsgLogger.H"

#include <iostream> //DEBUG

using namespace Utils::MsgLogger;

Logger::Logger(){}

void Logger::record_msg(Msg msg)
{
    m_messages[msg]++;
    int priority_int = -1;
    if (msg.priority == Priority::high) priority_int = 0;
    if (msg.priority == Priority::medium) priority_int = 1;
    if (msg.priority == Priority::low) priority_int = 2;
    std::cout << "DEBUG :" <<
        "    topic   : " << msg.topic << "\n" <<
        "    text    : " << msg.text << "\n" <<
        "    priority: " << priority_int << "\n" <<
        "    count   : " << m_messages[msg] << "\n";
}
