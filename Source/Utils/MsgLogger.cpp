/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "MsgLogger.H"

#include <algorithm>

using namespace MsgLogger;

void Logger::record_entry(
    Type type,
    Importance importance,
    std::string topic,
    std::string text)
{
    m_entries[type][importance][topic].push_back(text);
}

void Logger::print_warnings(std::stringstream& ss)
{

    auto& all_warnings = m_entries[Type::warning];

    const auto& low_prefix = "* [!  ]";
    const auto& medium_prefix = "* [!! ]";
    const auto& high_prefix = "* [!!!]";

    ss << "******************* WARNINGS! ********************* \n";

    for (auto& by_topic : all_warnings[Importance::high])
        aux_print_entries(high_prefix, by_topic.first, by_topic.second, ss);

    for (auto& by_topic : all_warnings[Importance::medium])
        aux_print_entries(medium_prefix, by_topic.first, by_topic.second, ss);

    for (auto& by_topic : all_warnings[Importance::low])
        aux_print_entries(low_prefix, by_topic.first, by_topic.second, ss);

    ss << "*************************************************** \n";
}

void
Logger::aux_print_entries(
    const std::string& prefix,
    const std::string& topic,
    const std::vector<std::string>& entries,
    std::stringstream& ss)
{
    for(const auto& msg : entries){
        ss << prefix
            << " [ " << topic << " ] "
            << msg << "\n";
    }
}