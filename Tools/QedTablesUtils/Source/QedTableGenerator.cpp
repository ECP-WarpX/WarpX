#include "QedTablesArgParser.H"

#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

using namespace ArgParser;
using namespace std;

const auto line_commands = vector<ArgParser::Key>{
    {"-h"                  , ArgType::NoArg  , "Prints all command line arguments"},
    {"--table"             , ArgType::String , "Either BW (Breit-Wheeler) or QS (Quantum Synchrotron)"},
    {"--mode"              , ArgType::String,  "Precision of the calculations: either DP (double) or SP (single)"},
    {"--dndt_chi_min"      , ArgType::Double , "minimum chi parameter for the dNdt table"},
    {"--dndt_chi_max"      , ArgType::Double , "maximum chi parameter for the dNdt table"},
    {"--dndt_how_many"     , ArgType::Integer, "number of points in the dNdt table"},
    {"--pair_chi_min"      , ArgType::Double , "minimum chi for the pair production table (BW only)"},
    {"--pair_chi_max"      , ArgType::Double , "maximum chi for the pair production table (BW only)"},
    {"--pair_chi_how_many" , ArgType::Integer, "number of chi points in the pair production table (BW only)"},
    {"--pair_frac_how_many", ArgType::Integer, "number of frac points in the pair production table (BW only)"},
    {"--em_chi_min"        , ArgType::Double , "minimum chi for the photon emission table (QS only)"},
    {"--em_chi_max"        , ArgType::Double , "maximum chi for the photon emission production table (QS only)"},
    {"--em_chi_how_many"   , ArgType::Integer, "number of chi points in the photon emission table (QS only)"},
    {"--em_frac_how_many"  , ArgType::Integer, "number of frac points in the photon emission table (QS only)"},
    {"-o"                  , ArgType::String , "filename to save the lookup table"}
};

void GenerateTable (const ParsedArgs& args);
void GenerateTableBW (const ParsedArgs& args);
void GenerateTableQS (const ParsedArgs& args);

bool IsDoublePrecision(const string& mode);

void PrintHelp (const vector<ArgParser::Key>& ll);

template <typename ContainerType, typename ElementType>
bool Contains (const ContainerType& container, const ElementType& el);

void AbortWithMessage(const string& msg);

void SuccessExit();

int main (int argc, char** argv)
{
    cout << "### QED Table Generator ###" << endl;
    const auto args_map = parse_args(line_commands, argc, argv);

    if (args_map.empty() || Contains(args_map, "-h")){
        PrintHelp(line_commands);
    }
    else{
        GenerateTable(args_map);
    }

    SuccessExit();
}


void GenerateTable (const ParsedArgs& args)
{
    if (!Contains(args, "--table"))
        AbortWithMessage("'--table' argument must be provided");

    if (GetVal<string>(args.at("--table")) == "BW"s)
        GenerateTableBW(args);
    else if (GetVal<string>(args.at("--table")) == "QS"s)
        GenerateTableQS(args);
    else
        AbortWithMessage("'--table' must be eiter 'QS' or 'BW'");
}

void GenerateTableBW (const ParsedArgs& args)
{
    const auto is_double = IsDoublePrecision(GetVal<string>(args.at("--mode")));
}

void GenerateTableQS (const ParsedArgs& args)
{
    const auto is_double = IsDoublePrecision(GetVal<string>(args.at("--mode")));
}

bool IsDoublePrecision(const string& mode)
{
    if (mode == "DP"s)
        return true;
    else if (mode == "SP")
        return false;
    else
        AbortWithMessage("'--mode' must be eiter 'DP' or 'SP'");

    return true;
}

void PrintHelp (const vector<ArgParser::Key>& cmd_list)
{
    cout << "Command line options: " << endl;

    for (const auto& el : cmd_list){
        const auto type = get<1>(el);
        string stype = "[??????]";
        if (type == ArgType::NoArg)
            stype = "[NO ARG]";
        else if (type == ArgType::String)
            stype = "[STRING]";
        else if (type == ArgType::Double)
            stype = "[DOUBLE]";
        else if (type == ArgType::Integer)
            stype = "[INTEGR]";

        cout << get<0>(el) <<
        " " << stype << " " <<
        get<2>(el) << endl;
    }

}

template <typename ContainerType, typename ElementType>
bool Contains (const ContainerType& container, const ElementType& el)
{
    return container.find(el) != std::end(container);
}

void AbortWithMessage(const std::string& msg)
{
    cout << "### ABORT : " << msg << std::endl;
    cout << "___________________________" << endl;
    exit(1);
}

void SuccessExit()
{
    cout << "___________________________" << endl;
    exit(0);
}