#include "MsgLogger.H"

#include "MsgLoggerSerialization.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <iostream>
#include <sstream>
#include <numeric>

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

    put_in_vec(msg.serialize(), serialized_msg_with_counter);
    put_in(this->counter, serialized_msg_with_counter);

    return serialized_msg_with_counter;
}

MsgWithCounter MsgWithCounter::deserialize (std::vector<char>::const_iterator& it)
{
    MsgWithCounter msg_with_counter;

    msg_with_counter.msg = Msg::deserialize(it);
    msg_with_counter.counter = get_out<int> (it);

    return msg_with_counter;
}

std::vector<char> MsgWithCounterAndRanks::serialize() const
{
    std::vector<char> serialized_msg_with_counter_and_ranks;

    put_in_vec(this->msg_with_counter.serialize(), serialized_msg_with_counter_and_ranks);
    put_in(this->all_ranks, serialized_msg_with_counter_and_ranks);
    put_in_vec(this->ranks, serialized_msg_with_counter_and_ranks);

    return serialized_msg_with_counter_and_ranks;
}

MsgWithCounterAndRanks
MsgWithCounterAndRanks::deserialize (std::vector<char>::const_iterator& it)
{
    MsgWithCounterAndRanks msg_with_counter_and_ranks;

    msg_with_counter_and_ranks.msg_with_counter = MsgWithCounter::deserialize(it);
    msg_with_counter_and_ranks.all_ranks = get_out<bool>(it);
    msg_with_counter_and_ranks.ranks = get_out_vec<int>(it);

    return msg_with_counter_and_ranks;
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

    const auto serialized_gather_rank_msg_list =
        collective_gather_get_serialized_gather_rank_msg_list(
            my_list, gather_rank);

    std::vector<int> gather_rank_msg_counters(gather_rank_how_many);

    auto to_send = std::map<Msg, int>{m_messages};

    auto package_lengths = std::vector<int>{};
    auto all_data = std::vector<char>{};
    if (!is_gather_rank){
        const auto gather_rank_msg_list =
            Logger::deserialize_msg_list(serialized_gather_rank_msg_list);

        int counter = 0;
        for (const auto& msg : gather_rank_msg_list){
            const auto pp = to_send.find(msg);
            if (pp != to_send.end()){
                gather_rank_msg_counters[counter] += pp->second;
                to_send.erase(msg);
            }
            counter++;
        }

        const auto package = aux_create_data_package(
            gather_rank_msg_counters,
            to_send
        );

        amrex::ParallelDescriptor::Gather(static_cast<int>(package.size()), gather_rank);
        const auto dummy_recc = std::vector<int>{};
        const auto dummy_recv = static_cast<char*>(nullptr);
        const auto dummy_disps = std::vector<int>{};
        amrex::ParallelDescriptor::Gatherv(
            package.data(), package.size(), dummy_recv,  dummy_recc,
            dummy_disps, gather_rank);
    }
    else{
        const int zero_size = 0;
        package_lengths =
            amrex::ParallelDescriptor::Gather(zero_size, gather_rank);

        auto disps = std::vector<int>(package_lengths.size());
        std::partial_sum(package_lengths.begin(), package_lengths.end(), disps.begin());
        std::rotate(disps.rbegin(), disps.rbegin()+1, disps.rend());
        disps[0] = 0;

        all_data.resize(std::accumulate(package_lengths.begin(),
            package_lengths.end(),0));
        if(all_data.size() == 0) all_data.resize(1);
        const auto dummy_send = static_cast<char*>(nullptr);
        amrex::ParallelDescriptor::Gatherv(
            dummy_send, 0, all_data.data(), package_lengths,
            disps, gather_rank);
    }

    auto list_with_ranks = std::vector<MsgWithCounterAndRanks>{};

    if(is_gather_rank){
        //aux_update_list_with_packaged_data();

        for (const auto& msg_with_counter : get_msg_with_counter_list()){
            MsgWithCounterAndRanks msg_with_counter_and_ranks;
            msg_with_counter_and_ranks.msg_with_counter = msg_with_counter;
            msg_with_counter_and_ranks.all_ranks = true;
            msg_with_counter_and_ranks.ranks = std::vector<int>{};

            list_with_ranks.emplace_back(msg_with_counter_and_ranks);
        }
    }


    if (gather_rank != m_io_rank){
        if(is_gather_rank){
            std::vector<char> package;
            for (const auto& el: list_with_ranks)
                put_in_vec<char>(el.serialize(), package);

            auto package_size = static_cast<int>(package.size());
            amrex::ParallelDescriptor::Send(&package_size, 1, m_io_rank, 0);
            amrex::ParallelDescriptor::Send(package, m_io_rank, 1);
            int list_size = static_cast<int>(my_list.size());
            amrex::ParallelDescriptor::Send(&list_size, 1, m_io_rank, 2);
        }
        else if (m_rank == m_io_rank){
            int vec_size = 0;
            amrex::ParallelDescriptor::Recv(&vec_size, 1, gather_rank, 0);
            std::vector<char> package(vec_size);
            amrex::ParallelDescriptor::Recv(package, gather_rank, 1);
            int list_size = 0;
            amrex::ParallelDescriptor::Recv(&list_size, 1, gather_rank, 2);
            const auto vv = package;
            auto it = vv.begin();
            for (int i = 0; i < list_size; ++i){
                list_with_ranks.emplace_back(
                    MsgWithCounterAndRanks::deserialize(it)
                );
            }
        }
    }

    return list_with_ranks;
}

std::vector<char>
Logger::collective_gather_get_serialized_gather_rank_msg_list(
    const std::vector<Msg>& my_list,
    int gather_rank) const
{
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

    return serialized_gather_rank_msg_list;
}


std::vector<char> Logger::serialize_msg_list(
    const std::vector<Msg>& msg_list) const
{
    std::vector<char> serialized;

    const auto how_many = static_cast<int> (msg_list.size());
    put_in (how_many, serialized);

    int i = 0;
    for (auto msg : msg_list){
        put_in_vec(msg.serialize(), serialized);
        std::cout << i++ << "/" << how_many << " :::::: "  << "|" <<
        msg.text << "|" <<
        msg.topic << "|" << std::endl;
    }

    std::cout << "ççççççççççççççç " << serialized.size() << std::endl;

    std::cout << "xx " << (int)serialized[0] << (int) serialized[1]
        << (int)serialized[2] << (int) serialized[3]
        << (int)serialized[4] << (int) serialized[5]
        << (int)serialized[6] << (int) serialized[7] << std::endl;

    return serialized;
}

std::vector<Msg> Logger::deserialize_msg_list(
    const std::vector<char>& serialized)
{
    std::cout << "yy " << (int)serialized[0] << (int) serialized[1]
        << (int)serialized[2] << (int) serialized[3]
        << (int)serialized[4] << (int) serialized[5]
        << (int)serialized[6] << (int) serialized[7] << std::endl;
    auto it = serialized.begin();

    const auto how_many = get_out<int>(it);
    auto msg_list = std::vector<Msg>{};

    msg_list.reserve(how_many);

    std::cout << "ççççççççççççççç " << serialized.size() << " " << it - serialized.begin() << std::endl;

    for (int i = 0; i < how_many; ++i){
        msg_list.emplace_back(Msg::deserialize(it));
        std::cout << i << "/" << how_many << " :::::: "  << "|" <<
        msg_list.back().text << "|" <<
        msg_list.back().topic << "|" << std::endl;
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

    return std::make_pair(pack[0], pack[1]);
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
                true,
                std::vector<int>{m_rank}});
    }
    return res;
}

std::vector<char>
Logger::aux_create_data_package(
    const std::vector<int>& gather_rank_msg_counters,
    const std::map<Msg, int>& to_send
    ) const
{
    std::vector<char> package;
    put_in(gather_rank_msg_counters.size(), package);
    for (const auto& el : gather_rank_msg_counters)
        put_in(el, package);
    put_in(to_send.size(), package);
    for (const auto& el : to_send)
        put_in_vec<char>(MsgWithCounter{el.first, el.second}.serialize(), package);

    return package;
}

