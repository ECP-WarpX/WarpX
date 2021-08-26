/* Copyright 2021 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_MSG_LOGGER_H_
#define WARPX_MSG_LOGGER_H_

#include <AMReX.H>

#include <map>
#include <string>
#include <utility>
#include <vector>

namespace Utils{
namespace MsgLogger{

    /** Priority is recorded together with messages. It influences
    * the display order and the appearance of a message.
    */
    enum class Priority
    {
        /** Low priority message */
        low,
        /** Medium priority message */
        medium,
        /** High priority message */
        high
    };

    /**
    * \brief This function converts a Priority into the corresponding
    * string (e.g, Priority::low --> "low")
    *
    * @param[in] priority the priority
    * @return the corresponding string
    */
    std::string PriorityToString(const Priority& priority);

    /**
    * \brief This function converts a string into the corresponding
    * priority (e.g, "low" --> Priority::low)
    *
    * @param[in] priority_string the priority string
    * @return the corresponding priority
    */
    Priority StringToPriority(const std::string& priority_string);

    /**
    * This struct represents a message, which is composed by
    * a topic, a text and a priority. It also provides methods for
    * serialization and deserialization.
    */
    struct Msg
    {
        std::string topic  /*! The message topic*/;
        std::string text   /*! The message text*/;
        Priority priority  /*! The priority of the message*/;

        /**
        * \brief This function returns a byte representation of the struct
        *
        * @return a byte vector
        */
        std::vector<char> serialize() const;

        /**
        * \brief This function generates a Msg struct from a byte vector
        *
        * @param[in] it iterator of a byte array
        * @return a Msg struct
        */
        static Msg deserialize(std::vector<char>::const_iterator& it);
    };

    /**
    * This struct represents a message with counter, which is composed
    * by a message and a counter. The latter is intended to store the
    * number of times a message is recorded. The struct also provides
    * methods for serialization and deserialization.
    */
    struct MsgWithCounter
    {
        Msg msg      /*! A message*/;
        int counter  /*! The counter*/;

        /**
        * \brief This function returns a byte representation of the struct
        *
        * @return a byte vector
        */
        std::vector<char> serialize() const;

        /**
        * \brief This function generates a MsgWithCounter struct from a byte vector
        *
        * @param[in] it iterator of a byte array
        * @return a MsgWithCounter struct
        */
        static MsgWithCounter deserialize(std::vector<char>::const_iterator& it);
    };

    /**
    * This struct represents a message with counter and ranks, which is
    * composed by a message with counter, a bool flag and a std::vector<int>.
    * The bool flag is used to store if a message is emitted by all the ranks.
    * The std::vector<int> is used to store the affected ranks
    * (note: when we switch to C++17, should we consider variants?).
    * The struct also provides methods for serialization and deserialization.
    */
    struct MsgWithCounterAndRanks
    {
        MsgWithCounter msg_with_counter /*! A message with counter*/;
        bool all_ranks                  /*! Flag to store if message is emitted by all ranks*/;
        std::vector<int> ranks          /*! Affected ranks*/;

        /**
        * \brief This function returns a byte representation of the struct
        *
        * @return a byte vector
        */
        std::vector<char> serialize() const;

        /**
        * \brief This function generates a MsgWithCounterAndRanks struct from a byte vector
        *
        * @param[in] it iterator of a byte array
        * @return a MsgWithCounterAndRanks struct
        */
        static MsgWithCounterAndRanks deserialize(std::vector<char>::const_iterator& it);
    };

    /**
    * \brief This implemets the < operator for Msg.
    * Warning messages are first ordered by priority (warning: high < medium < low
    * to give precedence to higher priorities), then by topic (alphabetically),
    * and finally by text (alphabetically).
    *
    * @param[in] l a Msg
    * @param[in] r a Msg
    * @return true if l<r, false otherwise
    */
    constexpr bool operator<(const Msg& l, const Msg& r)
    {
        return
            (l.priority > r.priority) ||
            ((l.priority == r.priority) && (l.topic < r.topic)) ||
            ((l.priority == r.priority) && (l.topic == r.topic) && (l.text < r.text));
    }

    /**
    * This class is responsible for storing messages and merging messages
    * collected by different processes.
    */
    class Logger
    {
        public:

        /**
        * \brief The constructor.
        */
        Logger();

        /**
        * \brief This function records a message
        *
        * @param[in] msg a Msg struct
        */
        void record_msg(Msg msg);

        /**
        * \brief This function returns a list of the recorded messages
        *
        * @return the list of messages
        */
        std::vector<Msg> get_msg_list() const;

        /**
        * \brief This function returns a list of the recorded messages
        * with the corresponding counter
        *
        * @return the list of messages with counters
        */
        std::vector<MsgWithCounter> get_msg_with_counter_list() const;

        /**
        * \brief This (collective) function generates a list of messages
        * with counters and emitting ranks by gathering data from
        * all the ranks
        *
        * @return the list of messages with counters and ranks if I/O rank, an empty vector otherwise
        */
        std::vector<MsgWithCounterAndRanks>
        collective_gather_msg_with_counter_and_ranks() const;

        private:

        std::vector<char> serialize_msg_list(
            const std::vector<Msg>& msg_list) const;
        static std::vector<Msg> deserialize_msg_list(
            const std::vector<char>& serialized);

        std::pair<int,int> aux_find_gather_rank_and_items(
            int how_many_items) const;

        std::vector<MsgWithCounterAndRanks>
        aux_get_msg_with_counter_and_ranks() const;

        std::vector<char>
        collective_gather_get_serialized_gather_rank_msg_list(
            const std::vector<Msg>& my_list,
            int gather_rank) const;

        std::vector<char>
        aux_create_data_package(
            const std::vector<int>& gather_rank_msg_counters,
            const std::map<Msg, int>& to_send
        ) const;

        void
        update_list_with_packaged_data(
            std::vector<MsgWithCounterAndRanks>& list_with_ranks,
            const std::vector<char>& all_data,
            const std::vector<int>& disps) const;

        int m_rank = 0         /*! MPI rank of the current process*/;
        int m_num_procs = 0    /*! Number of MPI ranks*/;
        int m_io_rank = 0      /*! Rank of the I/O process*/;
        bool m_am_i_io = false /*! Flag to store if the process is responsible for I/O*/;

        std::map<Msg, int> m_messages /*! This stores a map to associate warning messages with the corresponding counters*/;
    };
}
}

#endif //WARPX_MSG_LOGGER_H_
