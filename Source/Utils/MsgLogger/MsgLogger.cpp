/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "MsgLogger.H"

#include <algorithm>

using namespace MsgLogger;

namespace MsgLogger{
std::string importance_to_string(Importance imp)
{
    switch(imp){
        case Importance::low :
            return std::string{"low"};
        case Importance::medium :
            return std::string{"medium"};
        case Importance::high:
            return std::string{"high"};
        default:
            return std::string{"unknown!"};
    }
}
}

Logger::Logger(bool am_I_IO_proc):
    m_am_I_IO_proc{am_I_IO_proc} {};

void Logger::record_entry(
    Importance importance,
    std::string topic,
    std::string text)
{
    m_entries[importance][topic].push_back(
        Entry{text,1});
}

void Logger::collective_gather_entries(
        const std::string& tag_name,
        bool flush_non_IO_procs)
{
    m_entries_archive.push_back(
        TaggedEntryMap{tag_name, m_entries}
    );

    m_entries.clear();
}

void Logger::debug_raw_print(std::stringstream& ss)
{
    for (auto& tagged_map : m_entries_archive){
        ss << "--> TAG: " << tagged_map.tag << "\n";
        debug_raw_print_entry_map(tagged_map.entry_map, ss);
        ss << "________________________________\n";
    }

    if(!m_entries.empty()){
        ss << "--> UNMERGED: \n";
        debug_raw_print_entry_map(m_entries, ss);
        ss << "________________________________\n";
    }
}

void Logger::debug_raw_print_entry_map(
    EntryMap& entry_map, std::stringstream& ss)
{
    for (auto& by_imp : entry_map){
        ss << "   IMPORTANCE: " << importance_to_string(by_imp.first) << "\n";
        for (auto& by_topic : by_imp.second){
            ss << "      TOPIC: " << by_topic.first << "\n";
            for (auto& msg: by_topic.second){
                ss << "         MSG: " << msg.msg << " [" << msg.counter << "]\n";
            }
        }
    }
}


