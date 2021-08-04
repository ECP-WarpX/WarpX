
#include "WarnManager.H"

#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <sstream>

using namespace Utils;
using namespace Utils::MsgLogger;

WarnManager::WarnManager(){
    m_rank = amrex::ParallelDescriptor::MyProc();
}

void WarnManager::record_warning(
            std::string topic,
            std::string text,
            Priority priority)
{
    m_logger.record_msg(Msg{topic, text, priority});
}

std::string WarnManager::print_local_warnings(const std::string& when) const
{
    const auto all_warnings = aux_sort_messages(m_logger.get_msg_list());

    std::stringstream ss;

    ss << "\n" << aux_get_header(when, warn_line_size, false);

    if(all_warnings.size() == 0){
        ss << "No recorded warnings.\n";
    }
    else{
        for(const auto& warn_msg : all_warnings){
            ss << aux_print_warn_msg(warn_msg);
            ss << "*\n";
        }
    }

    ss << std::string(warn_line_size, '*') << "\n\n" ;

    return ss.str();
}

std::string WarnManager::print_global_warnings(const std::string& when) const
{
    auto all_warnings = m_logger.collective_gather_msg_lists();

    if(m_rank != amrex::ParallelDescriptor::IOProcessorNumber())
        return "[see I/O rank message]";

    all_warnings = aux_sort_messages(std::move(all_warnings));

    std::stringstream ss;

    ss << "\n" << aux_get_header(when, warn_line_size, true);

    if(all_warnings.size() == 0){
        ss << "No recorded warnings.\n";
    }
    else{
        for(const auto warn_msg : all_warnings){
            ss << aux_print_warn_msg(warn_msg);
            ss << "*\n";
        }
    }

    ss << std::string(warn_line_size, '*') << "\n\n" ;

    return ss.str();
}

void WarnManager::debug_read_warnings_from_input(amrex::ParmParse& params)
{
    std::vector<std::string> warnings;
    params.queryarr("test_warnings", warnings);

    for (const auto& warn : warnings){
        amrex::ParmParse pp_warn(warn);

        std::string topic;
        pp_warn.query("topic", topic);

        std::string msg;
        pp_warn.query("msg", msg);

        std::string spriority;
        pp_warn.query("priority", spriority);
        Priority priority = StringToPriority(spriority);

        int all_involved = 0;
        pp_warn.query("all_involved", all_involved);
        if(all_involved){
            this->record_warning(topic, msg, priority);
        }
        else{
            std::vector<int> who_involved;
            pp_warn.queryarr("who_involved", who_involved);
            if(std::find (who_involved.begin(), who_involved.end(), m_rank)
                != who_involved.end()){
                this->record_warning(topic, msg, priority);
            }
        }
    }

}

std::vector<MsgLogger::MsgWithCounter>
WarnManager::aux_sort_messages(
        std::vector<MsgLogger::MsgWithCounter>&& all_msg_with_counter) const
{
    std::sort(all_msg_with_counter.begin(), all_msg_with_counter.end(),
        [](const MsgWithCounter& a, const MsgWithCounter& b){
            return a.msg < b.msg;});
    return all_msg_with_counter;
}

std::vector<MsgLogger::MsgWithCounterAndRanks>
WarnManager::aux_sort_messages(
        std::vector<MsgLogger::MsgWithCounterAndRanks>&&
        all_msg_with_counter_and_ranks) const
{
    std::sort(all_msg_with_counter_and_ranks.begin(),
        all_msg_with_counter_and_ranks.end(),
        [](const MsgWithCounterAndRanks& a, const MsgWithCounterAndRanks& b){
            return a.msg_with_counter.msg < b.msg_with_counter.msg;});
    return all_msg_with_counter_and_ranks;
}

std::string WarnManager::aux_get_header(
    const std::string& when,
    int line_size,
    bool is_global) const
{
    const std::string warn_header{"**** WARNINGS "};

    std::stringstream ss;

    ss << warn_header <<
        std::string(line_size - static_cast<int>(warn_header.length()), '*') << "\n" ;

    if(is_global){
        ss << "* GLOBAL warning list  after " << " [ " <<  when << " ]\n*\n";
    }
    else{
        auto const mpi_rank = amrex::ParallelDescriptor::MyProc();
        ss << "* LOCAL" << " ( rank # " << mpi_rank << " ) "
            << " warning list  after " <<  when << "\n*\n";
    }

    return ss.str();
}

std::string
WarnManager::aux_print_warn_msg(
    const MsgLogger::MsgWithCounter& msg_with_counter) const
{
    std::stringstream ss;
    ss << "* --> ";
    if (msg_with_counter.msg.priority == MsgLogger::Priority::high)
        ss << "[!!!]";
    else if (msg_with_counter.msg.priority == MsgLogger::Priority::medium)
        ss << "[!! ]";
    else if (msg_with_counter.msg.priority == MsgLogger::Priority::low)
        ss << "[!  ]";
    else
        ss << "[???]";

    ss << " [" + msg_with_counter.msg.topic << "] ";

    if(msg_with_counter.counter == 2){
        ss << "[raised twice]\n";
    }
    else if(msg_with_counter.counter == 1){
        ss << "[raised once]\n";
    }
    else{
        ss << "[raised " << msg_with_counter.counter << " times]\n";
    }

    ss << aux_msg_formatter(msg_with_counter.msg.text, warn_line_size, warn_tab_size);

    return ss.str();
}

std::string
WarnManager::aux_print_warn_msg(
    const MsgLogger::MsgWithCounterAndRanks& msg_with_counter_and_ranks) const
{
    std::stringstream ss;
    ss << this->aux_print_warn_msg(msg_with_counter_and_ranks.msg_with_counter);

    return ss.str();
}

std::string
WarnManager::aux_msg_formatter(
        const std::string& msg,
        const int line_size,
        const int tab_size) const
{
    const auto prefix = "*" + std::string(tab_size, ' ');
    const auto prefix_length = static_cast<int>(prefix.length());

    std::stringstream ss_out;
    std::stringstream ss_msg{msg};

    std::string line;
    std::string word;

    while(std::getline(ss_msg, line,'\n')){
        ss_out << prefix;

        std::stringstream ss_line{line};
        int counter = prefix_length;

        while (ss_line >> word){
            const auto wlen = static_cast<int>(word.length());

            if(counter == prefix_length){
                ss_out << word;
                counter += wlen;
            }
            else{
                if (counter + wlen < line_size){
                    ss_out << " " << word;
                    counter += (wlen+1);
                }
                else{
                    ss_out << "\n" << prefix << word;
                    counter = prefix_length + wlen;
                }
            }
        }

        ss_out << '\n';
    }


    return ss_out.str();
}
