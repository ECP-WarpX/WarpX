
#include "MsgLogger.H"

#include "MsgLoggerSerialization.H"

#include <iostream> //DEBUG
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

Logger::Logger(){}

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

    const auto my_serialized_list =
        serialize_msg_list(this->get_msg_list());
    const auto how_many_items = my_serialized_list.size();

    const auto gather_rank = aux_find_gather_rank();

    return std::vector<MsgWithCounterAndRanks>{};
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

int Logger::aux_find_gather_rank() const
{
    return 0;
}

