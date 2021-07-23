#ifndef WARPX_MSG_LOGGER_SERIALIZATION_H_
#define WARPX_MSG_LOGGER_SERIALIZATION_H_

#include <array>
#include <cstring>
#include <string>
#include <vector>

namespace Utils{
namespace MsgLogger{

    /**
    * This function transform a variable of type T into a vector of chars holding its
    * byte representation and it appends this vector at the end of an
    * existing vector of chars. T must be either a trivially copyable type or an std::string
    * (see specialization below)
    *
    * @tparam T the variable type (must be trivially copyable)
    * @param[in] val a variable of type T to be serialized
    * @param[in, out] vec a reference to the vector to which the byte representation of val is appended
    */
    template <typename T>
    void put_in(const T& val, std::vector<char>& vec)
    {
        static_assert(std::is_trivially_copyable<T>(), "Cannot serialize \
            non-trivally copyable types.");

        const auto* ptr_val = reinterpret_cast<const char*>(&val);
        vec.insert(vec.end(), ptr_val, ptr_val+sizeof(T));
    }

    template <>
    void put_in<std::string> (const std::string& val, std::vector<char>& vec)
    {
        const char* c_str = val.c_str();
        const int length = static_cast<int>(val.length());

        vec.insert(vec.end(), c_str, c_str+length);
    }


    /**
    * This function extract a variable of type T from a byte vector, at the position
    * given by a std::vector<char> iterator. The iterator is then advanced according to
    * the number of bytes read from the byte vector. T must be either a trivially copyable type
    * or an std::string (see specialization below).
    *
    * @tparam T the variable type (must be trivially copyable)
    * @param[in, out] it the iterator to a byte vector
    * @return the variable extracted from the byte array
    */
    template<typename T>
    T get_out(std::vector<char>::const_iterator& it)
    {
        static_assert(std::is_trivially_copyable<T>(), "Cannot extract \
            non-trivally copyable types from char vectors.");

        auto temp = std::array<char, sizeof(T)>{};
        std::copy(it, it + sizeof(T), temp.begin());
        it += sizeof(T);
        T res;
        std::memcpy(&res, temp.data(), sizeof(T));
        return res;
    }

    template<>
    std::string get_out<std::string> (std::vector<char>::const_iterator& it)
    {
        const int length = get_out<int> (it);
        const auto str = std::string{it, it+length};
        it += length;
        return str;
    }
}
}

#endif //WARPX_MSG_LOGGER_SERIALIZATION_H_
