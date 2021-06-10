/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "MsgLogger.H"

#include <AMReX_Config.H>

#include <algorithm>

using namespace MsgLogger;

void Logger::record_entry(
    Type type,
    Importance importance,
    std::string topic,
    std::string text)
{
    #ifdef AMREX_USE_OMP
    #pragma omp critical
    #endif
    {
        m_entries[type][importance][topic].push_back(text);
    }
}

void Logger::print_warnings(std::stringstream& ss)
{
    auto& all_warnings = m_entries[Type::warning];

    if (all_warnings.empty()){
        ss << "\n * No warnings have been raised!\n\n";
        return;
    }

    const auto low_prefix           = "* [!  ]";
    const auto medium_prefix        = "* [!! ]";
    const auto high_prefix          = "* [!!!]";
    const auto new_line_skip        = "*       ";
    const auto max_line_length      = 60;

    ss << "\n";
    ss << "******************* WARNINGS! ********************* \n";

    for (auto& by_topic : all_warnings[Importance::high])
        aux_print_entries(high_prefix, by_topic.first, by_topic.second, max_line_length, ss);

    for (auto& by_topic : all_warnings[Importance::medium])
        aux_print_entries(medium_prefix, by_topic.first, by_topic.second, max_line_length, ss);

    for (auto& by_topic : all_warnings[Importance::low])
        aux_print_entries(low_prefix, by_topic.first, by_topic.second, max_line_length, ss);

    ss << "*************************************************** \n";
    ss << "\n\n";

}

void
Logger::aux_print_entries(
    const std::string& prefix,
    const std::string& topic,
    const std::vector<std::string>& entries,
    const int max_line_length,
    std::stringstream& ss) const
{
    const int first_line_offset;
    for(const auto& msg : entries){
        ss << prefix
            << " [ " << topic << " ] "
            << aux_print_formatter(
                msg,
                prefix,
                max_line_length,
                first_line_offset) << "\n";
    }
}

std::string
Logger::aux_msg_formatter(
        std::string msg,
        std::string& new_line_skip,
        const int max_line_length) const
{
    std::stringstream ss_out;

    std::stringstream ss_in{msg};
    std::string line;
    while(std::getline(msg,line,'\n')){
        std::stringstream ss_line{line};
        std::string word;
        while (std::getline(ss_line, word, " "))
        {

        }


        ss_out << '\n';
    }

    return ss_out.str();
}