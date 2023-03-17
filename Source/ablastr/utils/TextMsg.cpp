/* Copyright 2022 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "TextMsg.H"

#include <AMReX.H>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>


namespace
{
    constexpr auto err_prefix       = "### ERROR   : ";
    constexpr auto warn_prefix      = "!!! WARNING : ";
    constexpr auto info_prefix      = "--- INFO    : ";
    constexpr auto err_line_prefix  = "#            ";
    constexpr auto line_prefix      = "             ";
    constexpr auto line_length   = 66;

    std::string Msg (
        const std::string& msg,
        const std::string& msg_prefix,
        const std::string& msg_line_prefix,
        const int msg_line_length,
        const bool do_text_wrapping)
    {
        if(!do_text_wrapping){
            return msg_prefix + msg + "\n";
        }

        const auto wrapped_text = ablastr::utils::automatic_text_wrap(
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

std::string
ablastr::utils::TextMsg::Err (const std::string& msg, const bool do_text_wrapping)
{
    return ::Msg(
        msg, ::err_prefix, ::err_line_prefix, ::line_length, do_text_wrapping);
}

std::string
ablastr::utils::TextMsg::Info (const std::string& msg, const bool do_text_wrapping)
{
    return ::Msg(
        msg, ::info_prefix, ::line_prefix, ::line_length, do_text_wrapping);
}

std::string
ablastr::utils::TextMsg::Warn (const std::string& msg, const bool do_text_wrapping)
{
    return ::Msg(
        msg, ::warn_prefix, ::line_prefix, ::line_length, do_text_wrapping);
}

void
ablastr::utils::TextMsg::Assert (const char* ex, const char* file, const int line, const std::string& msg)
{
    const auto n_msg = "\n" + Err(msg);
    amrex::Assert(ex , file, line , n_msg.c_str());
}

std::vector< std::string >
ablastr::utils::automatic_text_wrap (
    const std::string& text, const int max_line_length)
{

    auto ss_text = std::stringstream{text};
    auto wrapped_text_lines = std::vector< std::string >{};

    std::string line;
    while(std::getline(ss_text, line,'\n')){

        auto ss_line = std::stringstream{line};
        int counter = 0;
        std::stringstream ss_line_out;
        std::string word;

        while (ss_line >> word){
            const auto wlen = static_cast<int>(word.length());

            if(counter == 0){
                ss_line_out << word;
                counter += wlen;
            }
            else{
                if (counter + wlen < max_line_length){
                    ss_line_out << " " << word;
                    counter += (wlen+1);
                }
                else{
                    wrapped_text_lines.push_back(ss_line_out.str());
                    ss_line_out.str("");
                    ss_line_out << word;
                    counter = wlen;
                }
            }
        }

        wrapped_text_lines.push_back(ss_line_out.str());
    }

    return wrapped_text_lines;
}
