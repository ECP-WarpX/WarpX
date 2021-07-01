/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MSG_LOGGER_H_
#define WARPX_MSG_LOGGER_H_

#include <mpi.h>

#include <vector>
#include <map>
#include <string>
#include <sstream>

namespace MsgLogger{

enum class Importance
{
    low,
    medium,
    high
};

std::string importance_to_string(Importance imp);

struct Entry
{
    std::string msg;
    int counter = 0;
};

typedef
    std::map<Importance,
        std::map<std::string,
            std::vector<Entry>>> EntryMap;

struct TaggedEntryMap
{
    std::string tag;
    EntryMap entry_map;
};

class Logger
{

public:
    Logger(bool am_I_IO_proc);

    void record_entry(
        Importance importance,
        std::string topic,
        std::string text);

    void collective_gather_entries(
        const std::string& tag_name,
        bool flush_non_IO_procs = true);

    void debug_raw_print(std::stringstream& ss);

private:

    void debug_raw_print_entry_map(EntryMap& entry_map, std::stringstream& ss);

    bool m_am_I_IO_proc = false;

    EntryMap m_entries;
    std::vector<TaggedEntryMap> m_entries_archive;
};

}

#endif //WARPX_MSG_LOGGER_H_
