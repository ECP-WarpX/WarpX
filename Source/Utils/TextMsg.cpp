/* Copyright 2022 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "TextMsg.H"

#include "WarpXUtil.H"

#include <AMReX_Print.H>

#include <sstream>

using namespace Utils::TextMsg;

namespace
{
    constexpr auto err_prefix       = "### ERROR   : ";
    constexpr auto warn_prefix      = "!!! WARNING : ";
    constexpr auto info_prefix      = "--- INFO    : ";
    constexpr auto err_line_prefix  = "#            ";
    constexpr auto line_prefix      = "             ";
    constexpr auto line_length   = 66;

    std::string Msg(
        const std::string& msg,
        const std::string& msg_prefix,
        const std::string& msg_line_prefix,
        const int msg_line_length,
        const bool do_text_wrapping)
    {
        if(!do_text_wrapping){
            return msg_prefix + msg + "\n";
        }

        const auto wrapped_text = WarpXUtilStr::automatic_text_wrap(
            msg, msg_line_length);

        std::stringstream ss_out;

        std::for_each(std::begin(wrapped_text), std::end(wrapped_text),
            [&,ln=0](const auto& line) mutable {
                ss_out << ((ln++ == 0) ? msg_prefix : msg_line_prefix);
                ss_out << line << "\n";
            });

        return ss_out.str();
    }
}

std::string Err(const std::string& msg, const bool do_text_wrapping = true)
{
    return ::Msg(
        msg, ::err_prefix, ::err_line_prefix, ::line_length, do_text_wrapping);
}

std::string Info(const std::string& msg, const bool do_text_wrapping = true)
{
    return ::Msg(
        msg, ::info_prefix, ::line_prefix, ::line_length, do_text_wrapping);
}

std::string Warn(const std::string& msg, const bool do_text_wrapping = true)
{
    return ::Msg(
        msg, ::warn_prefix, ::line_prefix, ::line_length, do_text_wrapping);
}
