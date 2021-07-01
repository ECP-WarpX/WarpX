/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "MsgLogger.H"

//#include <AMReX_Config.H>
//#include <AMReX_ParallelDescriptor.H>

#include <algorithm>

using namespace MsgLogger;

Logger::Logger()
{
    bool m_am_I_IO_proc = true;
};

void Logger::record_entry(
    Type type,
    Importance importance,
    std::string topic,
    std::string text)
{
    m_entries[type][importance][topic].push_back(
        Entry{text,1});
}

void Logger::print_all_warnings(
    std::stringstream& ss,
    const WarnStyle& ws)
{

    ss << "\n";
    ss << ws.header << "\n";

    for (auto& tagged_entry_map : m_entries_archive){
        ss << ws.tag_prefix << tagged_entry_map.tag << "\n";

        ss << 
    }

    auto& all_warnings = m_entries[Type::warning];

    if (all_warnings.empty()){
        ss << ws.no_warning_msg << "\n";
        return;
    }



    for (auto& by_topic : all_warnings[Importance::high])
        aux_print_entries(ws.high_prefix, by_topic.first, by_topic.second, ws, ss);

    for (auto& by_topic : all_warnings[Importance::medium])
        aux_print_entries(ws.medium_prefix, by_topic.first, by_topic.second, ws, ss);

    for (auto& by_topic : all_warnings[Importance::low])
        aux_print_entries(ws.low_prefix, by_topic.first, by_topic.second, ws, ss);

    ss << ws.footer << "\n";
}

void
Logger::aux_print_entries(
    const std::string& prefix,
    const std::string& topic,
    const std::vector<Entry>& entries,
    const WarnStyle& ws,
    std::stringstream& ss) const
{
    const std::string first_line_prefix =
        prefix + ws.topic_left + topic + ws.topic_right;
    const auto first_line_offset = first_line_prefix.length();
    for(const auto& entry : entries){
        ss << first_line_prefix <<
            aux_msg_formatter(
                entry.msg,
                first_line_offset,
                ws);
    }
}

std::string
Logger::aux_msg_formatter(
        const std::string& msg,
        const int first_line_offset,
        const WarnStyle& ws) const
{
    std::stringstream ss_out;

    std::stringstream ss_in{msg};
    std::string line;
    bool is_first = true;
    while(std::getline(ss_in, line,'\n')){
        if(!is_first) ss_out << ws.new_line_prefix;
        ss_out << line << '\n';
        is_first = false;
    }
    return ss_out.str();
}


/*
void Logger::record_collective_entry(
        Type type,
        Importance importance,
        std::string topic,
        std::string text,
        bool affects_me,
        MPI_Comm comm,
        int ranks_per_line)
{
    bool is_someone_affected = false;

    int world_rank;
    MPI_Comm_rank(comm, &world_rank);
    int world_size;
    MPI_Comm_size(comm, &world_size);

    auto flag_vector = std::vector<char>(world_size);

    const char c_affects_me = affects_me;
    MPI_Allgather(&c_affects_me, 1, MPI_CHAR, flag_vector.data(), 1, MPI_CHAR, comm);

    for(auto& el : flag_vector){
        if(el){
            is_someone_affected = true;
            break;
        }
    }

    if(!is_someone_affected) return;

    std::stringstream ss_collective_text;

    ss_collective_text << text << "\n";

    int counter = 0;

    for (int rr = 0; rr < flag_vector.size(); ++rr){
        if (flag_vector[rr]){
            if(counter != 0) ss_collective_text << " ";
            ss_collective_text << rr;
            counter++;
        }
        if(counter >= ranks_per_line){
            counter = 0;
            ss_collective_text << "\n";
        }
    }

    if (counter != 0) ss_collective_text << "\n";

    m_entries[type][importance][topic].push_back(ss_collective_text.str());
}
*/