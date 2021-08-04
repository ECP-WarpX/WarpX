
#include "MsgLogger.H"

#include "MsgLoggerSerialization.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <iostream>
#include <sstream>

using namespace Utils::MsgLogger;

std::vector<char> MsgWithCounter::serialize() const
{
    std::vector<char> serialized_msg;

    put_in(this->msg.topic, serialized_msg);

    put_in(this->msg.text, serialized_msg);

    const int int_priority = static_cast<int>(this->msg.priority);
    put_in(int_priority, serialized_msg);

    put_in(this->counter, serialized_msg);

    return serialized_msg;
}

MsgWithCounter MsgWithCounter::deserialize (const std::vector<char>& serialized)
{
    MsgWithCounter msg_with_counter;

    auto it = serialized.begin();

    msg_with_counter.msg.topic = get_out<std::string> (it);
    msg_with_counter.msg.text = get_out<std::string> (it);
    msg_with_counter.msg.priority = static_cast<Priority> (get_out<int> (it));
    msg_with_counter.counter = get_out<typeof(counter)> (it);

    return msg_with_counter;
}

Logger::Logger(){
    m_rank = amrex::ParallelDescriptor::MyProc();
    m_num_procs = amrex::ParallelDescriptor::NProcs();
    m_io_rank = amrex::ParallelDescriptor::IOProcessorNumber();
    m_am_i_io = (m_rank == m_io_rank);
}

void Logger::record_msg(Msg msg)
{
    m_messages[msg]++;
}

std::vector<MsgWithCounter> Logger::get_msg_list() const
{
    auto res = std::vector<MsgWithCounter>{};

    for (auto msg : m_messages){
        res.emplace_back(MsgWithCounter{msg.first, msg.second});
    }

    return res;
}

std::vector<MsgWithCounterAndRanks> Logger::collective_gather_msg_lists() const
{

    const auto my_list =
        get_msg_list();
    const auto how_many_items = my_list.size();

    int gather_rank = 0;
    int gather_rank_how_many = 0;
    std::tie(gather_rank, gather_rank_how_many) =
        aux_find_gather_rank_and_items(how_many_items);

    if(gather_rank_how_many == 0)
        return std::vector<MsgWithCounterAndRanks>{};

    auto list_with_ranks = std::vector<MsgWithCounterAndRanks>{};

    for (const auto& el : my_list){
        list_with_ranks.emplace_back(MsgWithCounterAndRanks{
            el,
            std::vector<int>{m_io_rank},
        });
    }

    return list_with_ranks;
}


std::vector<char> Logger::serialize_msg_list(
    const std::vector<MsgWithCounter>& msg_list) const
{
    std::vector<char> serialized;

    const auto how_many = static_cast<int> (msg_list.size());
    put_in (how_many, serialized);

    for (auto msg_w_counter : msg_list){
        const auto serialized_msg_w_counter = msg_w_counter.serialize();
        serialized.insert(serialized.end(),
            serialized_msg_w_counter.begin(),
            serialized_msg_w_counter.end());
    }
    return serialized;
}

std::pair<int,int> Logger::aux_find_gather_rank_and_items(int how_many_items) const
{
    int max_items = 0;
    int max_rank = 0;

    auto num_msg =
        amrex::ParallelDescriptor::Gather(how_many_items, m_io_rank);

    if (m_am_i_io){
        const auto it_max = std::max_element(num_msg.begin(), num_msg.end());
        max_items = *it_max;
        max_rank = (max_items == how_many_items) ?
            m_io_rank : it_max - num_msg.begin();
    }

    auto pack = std::array<int,2>{max_rank, max_items};
    amrex::ParallelDescriptor::Bcast(pack.data(), 2, m_io_rank);

    return std::make_pair(max_rank, max_items);
}

