
#include "WarnManager.H"

#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <sstream>

using namespace Utils;
using namespace Utils::MsgLogger;

WarnManager::WarnManager(){}

void WarnManager::record_warning(
            std::string topic,
            std::string text,
            Priority priority)
{
    m_logger.record_msg(Msg{topic, text, priority});
}

std::string WarnManager::print_local_warnings(const std::string& when)
{
    auto all_warnings = aux_sort_messages(m_logger.get_msg_list());

    std::stringstream ss;

    ss << "\n" << aux_get_header(when, warn_line_size, false);

    if(all_warnings.size() == 0){
        ss << "No recorded warnings.\n";
    }
    else{
        for(const auto warn_msg : all_warnings){
            ss << aux_print_warn_msg(warn_msg);
        }
    }

    ss << std::string(warn_line_size, '*') << "\n\n" ;

    return ss.str();
}

std::vector<MsgLogger::MsgWithCounter>
WarnManager::aux_sort_messages(
        std::vector<MsgLogger::MsgWithCounter>&& all_msg_with_counter)
{
    std::sort(all_msg_with_counter.begin(), all_msg_with_counter.end(),
        [](const MsgWithCounter& a, const MsgWithCounter& b){
            return a.msg < b.msg;});
    return all_msg_with_counter;
}

std::string WarnManager::aux_get_header(
    const std::string& when,
    int line_size,
    bool is_global)
{
    const std::string warn_header{"**** WARNINGS "};

    std::stringstream ss;

    ss << warn_header <<
        std::string(line_size - static_cast<int>(warn_header.length()), '*') << "\n" ;

    if(is_global){
        ss << "* GLOBAL warning list  after " << " [ " <<  when << " ]\n";
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
    const MsgLogger::MsgWithCounter& msg_with_counter)
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

    ss << aux_print_line(msg_with_counter.msg.text, warn_line_size, warn_tab_size);

    ss << "*\n";

    return ss.str();
}

std::string
WarnManager::aux_print_line(
    const std::string& str, int line_size, int tab_size)
{
    /* TODO */
    return "*" + std::string(tab_size, ' ') + str + "\n";
}
