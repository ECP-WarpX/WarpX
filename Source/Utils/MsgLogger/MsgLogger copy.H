/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_MSG_LOGGER_H_
#define WARPX_MSG_LOGGER_H_

#include <vector>
#include <map>
#include <string>
#include <sstream>

namespace MsgLogger{

enum class Type
{
    warning
};

enum class Importance
{
    low,
    medium,
    high
};

struct WarnStyle
{
    std::string header;
    std::string tag_prefix;
    std::string tag_delimiter;
    std::string low_prefix;
    std::string medium_prefix;
    std::string high_prefix;
    std::string new_line_prefix;
    std::string topic_left;
    std::string topic_right;
    std::string footer;
    std::string no_warning_msg;
};

const auto default_warn_style =
    WarnStyle{
        "******************* WARNINGS **********************",
        "*---> ",
        "*__________________________________________________",

        "* [!  ]", "* [!! ]", "* [!!!]",
        "*       ",
        " [ ", " ] ",
        "***************************************************",
        "* No warnings have been raised!"};

struct Entry
{
    std::string msg;
    int counter = 0;
}

typedef
    std::map<Type,
        std::map<Importance,
            std::map<std::string,
                std::vector<Entry>>>> EntryMap;

struct TaggedEntryMap
{
    std::string tag;
    EntryMap entry_map;
}

class Logger
{

public:
    Logger();

    void record_entry(
        Type type,
        Importance importance,
        std::string topic,
        std::string text);

    void collective_gather_entries(
        const std::string& tag_name,
        bool flush_non_IO_procs = true);

    void print_all_warnings(std::stringstream& ss,
        const WarnStyle& ws = default_warn_style);

private:

    void

    aux_print_warn_epoch

    aux_print_entries(
        const std::string& prefix,
        const std::string& topic,
        const std::vector<Entry>& entries,
        const WarnStyle& ws,
        std::stringstream& ss) const;

    std::string
    aux_msg_formatter(
        const std::string& msg,
        const int first_line_offset,
        const WarnStyle& ws) const;

    EntryMap m_entries;
    std::vector<TaggedEntryMap> m_entries_archive;

    bool m_am_I_IO_proc = false;
};

}

#endif //WARPX_MSG_LOGGER_H_
