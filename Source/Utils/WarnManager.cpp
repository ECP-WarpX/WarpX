
#include "WarnManager.H"

using namespace Utils;
using namespace Utils::MsgLogger;

WarnManager::WarnManager(){}

void WarnManager::record_warning(
            std::string topic,
            std::string text,
            Priority priority)
{
    m_logger.record_msg(Msg{topic, text, priority});
}
