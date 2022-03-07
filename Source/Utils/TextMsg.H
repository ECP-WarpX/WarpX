/* Copyright 2022 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_TEXT_MSG_H_
#define WARPX_TEXT_MSG_H_

#include <string>

#define WARPX_ALWAYS_ASSERT_WITH_MESSAGE(EX,MSG) (EX)?((void)0):Utils::TextMsg::Assert( # EX , __FILE__, __LINE__ , MSG)

namespace Utils
{
namespace TextMsg
{
    /** \brief This function formats a text message as an error message,
    *   adding the '### ERROR: ' prefix and (by default) performing text wrapping
    *
    * @param[in] msg the string to be formatterd
    * @param[in] do_text_wrapping if true, the message is automatically
    *
    * @return the formatted message
    */
    std::string Err(const std::string& msg, const bool do_text_wrapping = true);

    /** \brief This function formats a text message as an info message,
    *   adding the '### INFO: ' prefix and (by default) performing text wrapping
    *
    * @param[in] msg the string to be formatterd
    * @param[in] do_text_wrapping if true, the message is automatically
    *
    * @return the formatted message
    */
    std::string Info(const std::string& msg, const bool do_text_wrapping = true);

    /** \brief This function formats a text message as a warning message,
    *   adding the '### WARN: ' prefix and (by default) performing text wrapping.
    *   Warning: this format is not used by the WarningLogger, which has an internal,
    *   dedicated, formatter.
    *
    * @param[in] msg the string to be formatterd
    * @param[in] do_text_wrapping if true, the message is automatically
    *
    * @return the formatted message
    */
    std::string Warn(const std::string& msg, const bool do_text_wrapping = true);

    void Assert(const char* ex, const char* file, const int line, const std::string& msg);
}
}

#endif //WARPX_TEXT_MSG_H_
