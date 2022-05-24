#ifndef QED_TABLES_ARG_PARSER_
#define QED_TABLES_ARG_PARSER_

#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace ArgParser
{
    enum class ArgType { NoArg, Integer, Double, String };

    using ArgVal = std::optional<std::variant<int, double, std::string>>;
    using Key = std::tuple<std::string, ArgType, std::string>;
    using ParsedArgs = std::unordered_map<std::string, ArgVal>;

    template <typename T>
    T GetVal(const ArgVal& arg_val)
    {
        return std::get<T>(*arg_val);
    }

    ParsedArgs parse_args (const std::vector<Key>& keys, const int argc, char const* const* argv);
};


#endif //QED_TABLES_ARG_PARSER_
