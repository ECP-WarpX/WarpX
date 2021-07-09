#ifndef WARPX_WARN_MANAGER_H_
#define WARPX_WARN_MANAGER_H_

#include "WarnManager_fwd.H"

#include "MsgLogger/MsgLogger.H"

#include <string>


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

        private:

        MsgLogger::Logger m_logger;
    };
}

#endif //WARPX_WARN_MANAGER_H_
