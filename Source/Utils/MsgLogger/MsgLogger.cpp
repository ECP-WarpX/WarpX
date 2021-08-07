
#include "MsgLogger.H"

#include "MsgLoggerSerialization.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <iostream>
#include <sstream>

using namespace Utils::MsgLogger;

std::vector<char> Msg::serialize() const
{
    std::vector<char> serialized_msg;

    put_in(this->topic, serialized_msg);
    put_in(this->text, serialized_msg);
    const int int_priority = static_cast<int>(this->priority);
    put_in(int_priority, serialized_msg);

    return serialized_msg;
}

Msg Msg::deserialize (std::vector<char>::const_iterator& it)
{
    Msg msg;

    msg.topic = get_out<std::string> (it);
    msg.text = get_out<std::string> (it);
    msg.priority = static_cast<Priority> (get_out<int> (it));

    return msg;
}

std::vector<char> MsgWithCounter::serialize() const
{
    std::vector<char> serialized_msg_with_counter;

    put_in(this->counter, serialized_msg_with_counter);
    const auto serialized_msg =  msg.serialize();
    serialized_msg_with_counter.insert(serialized_msg_with_counter.end(),
       serialized_msg.begin(), serialized_msg.end());

    return serialized_msg_with_counter;
}

MsgWithCounter MsgWithCounter::deserialize (std::vector<char>::const_iterator& it)
{
    MsgWithCounter msg_with_counter;

    msg_with_counter.counter = get_out<typeof(counter)> (it);
    msg_with_counter.msg = Msg::deserialize(it);

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

std::vector<Msg> Logger::get_msg_list() const
{
    auto res = std::vector<Msg>{};

    for (auto msg_w_counter : m_messages){
        res.push_back(msg_w_counter.first);
    }

    return res;
}

std::vector<MsgWithCounter> Logger::get_msg_with_counter_list() const
{
    auto res = std::vector<MsgWithCounter>{};

    for (auto msg : m_messages){
        res.emplace_back(MsgWithCounter{msg.first, msg.second});
    }

    return res;
}



std::vector<MsgWithCounterAndRanks>
Logger::collective_gather_msg_with_counter_and_ranks() const
{
    if (m_num_procs == 1)
        return aux_get_msg_with_counter_and_ranks();

    const auto my_list =
        get_msg_list();
    const auto how_many_items = my_list.size();

    int gather_rank = 0;
    int gather_rank_how_many = 0;
    std::tie(gather_rank, gather_rank_how_many) =
        aux_find_gather_rank_and_items(how_many_items);

    if(gather_rank_how_many == 0)
        return std::vector<MsgWithCounterAndRanks>{};

    const bool is_gather_rank = (m_rank == gather_rank);

    auto serialized_gather_rank_msg_list = std::vector<char>{};
    int size_serialized_gather_rank_msg_list = 0;
    if (is_gather_rank){
        serialized_gather_rank_msg_list = serialize_msg_list(my_list);
        size_serialized_gather_rank_msg_list =
            serialized_gather_rank_msg_list.size();
    }
    amrex::ParallelDescriptor::Bcast(
        &size_serialized_gather_rank_msg_list, 1, gather_rank);
    if (!is_gather_rank)
        serialized_gather_rank_msg_list.resize(
            size_serialized_gather_rank_msg_list);
    amrex::ParallelDescriptor::Bcast(
        serialized_gather_rank_msg_list.data(),
        size_serialized_gather_rank_msg_list, gather_rank);

    if (!is_gather_rank){
        const auto gather_rank_msg_list =
            Logger::deserialize_msg_list(serialized_gather_rank_msg_list);

        for (auto ee : gather_rank_msg_list){
            std::cout << "@@@@@@@@@@@@@@@@@ " << ee.text << std::endl;
        }
    }

    auto list_with_ranks = std::vector<MsgWithCounterAndRanks>{};

    for (const auto& el : my_list){
        list_with_ranks.emplace_back(MsgWithCounterAndRanks{
            el,
            1,
            std::vector<int>{m_io_rank}
        });
    }

    return list_with_ranks;
}


std::vector<char> Logger::serialize_msg_list(
    const std::vector<Msg>& msg_list) const
{
    std::vector<char> serialized;

    const auto how_many = static_cast<int> (msg_list.size());
    put_in (how_many, serialized);

    for (auto msg : msg_list){
        const auto serialized_msg = msg.serialize();
        serialized.insert(serialized.end(),
            serialized_msg.begin(),
            serialized_msg.end());
    }
    return serialized;
}

std::vector<Msg> Logger::deserialize_msg_list(
    const std::vector<char>& serialized)
{

    auto it = serialized.begin();

    const auto how_many = get_out<int>(it);
    auto msg_list = std::vector<Msg>{};
    msg_list.reserve(how_many);

    for (int i = 0; i < how_many; ++i){
        msg_list.emplace_back(Msg::deserialize(it));
    }

    return msg_list;
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

std::vector<MsgWithCounterAndRanks>
Logger::aux_get_msg_with_counter_and_ranks() const
{
    std::vector<MsgWithCounterAndRanks> res;
    for (const auto& el : m_messages)
    {
        res.emplace_back(
            MsgWithCounterAndRanks{
                MsgWithCounter{el.first, el.second},
                    std::vector<int>{m_rank}});
    }
    return res;
}
