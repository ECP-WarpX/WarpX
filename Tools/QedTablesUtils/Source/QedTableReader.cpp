#include "QedTableCommons.H"
#include "ArgParser/QedTablesArgParser.H"

#include <picsar_qed/physics/breit_wheeler/breit_wheeler_engine_tables.hpp>
#include <picsar_qed/physics/quantum_sync/quantum_sync_engine_tables.hpp>
#include <picsar_qed/utils/serialization.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace ArgParser;
using namespace std;

namespace pxr_sr = picsar::multi_physics::utils::serialization;
namespace pxr_bw = picsar::multi_physics::phys::breit_wheeler;
namespace pxr_qs = picsar::multi_physics::phys::quantum_sync;

const auto line_commands = vector<ArgParser::Key>{
    {"-h"                  , ArgType::NoArg  , "Prints all command line arguments"},
    {"-i"                  , ArgType::String , "Name of the file to open"},
    {"--table"             , ArgType::String , "Either BW (Breit-Wheeler) or QS (Quantum Synchrotron)"},
    {"--mode"              , ArgType::String,  "Precision of the calculations: either DP (double) or SP (single)"},
    {"-o"                  , ArgType::String , "filename to save the lookup table in human-readable format"}
};

void ReadTable (const ParsedArgs& args);
template <typename RealType>
void ReadTableBW (const string& input_file, const string& outfile_name);
template <typename RealType>
void ReadTableQS (const string& input_file, const string& outfile_name);

int main (int argc, char** argv)
{
    cout << "### QED Table Reader ###" << endl;
    const auto args_map = ParseArgs(line_commands, argc, argv);

    if (args_map.empty() || Contains(args_map, "-h")){
        PrintHelp(line_commands);
    }
    else{
        ReadTable(args_map);
    }

    SuccessExit();
}

void ReadTable (const ParsedArgs& args)
{
    if (!Contains(args, "--table"))
        AbortWithMessage("'--table' argument must be provided");

    if (!Contains(args, "--mode"))
        AbortWithMessage("'--mode' argument must be provided");

    if (!Contains(args, "-o"))
        AbortWithMessage("'-o' argument must be provided");

    if (!Contains(args, "-i"))
        AbortWithMessage("'-i' argument must be provided");

    const auto which_table = GetVal<string>(args.at("--table"));
    const auto mode = GetVal<string>(args.at("--mode"));
    const auto input_file = GetVal<string>(args.at("-i"));
    const auto outfile_name = GetVal<string>(args.at("-o"));

    bool use_double = false;

    if (mode == "DP"s)
        use_double = true;
    else if (mode == "SP"s)
        use_double = false;
    else
        AbortWithMessage("'--mode' must be eiter 'DP' or 'SP'");

    if (which_table == "BW"s){
        if (use_double)
            ReadTableBW<double>(input_file, outfile_name);
        else
            ReadTableBW<float>(input_file, outfile_name);
    }
    else if (which_table == "QS"s)
        if (use_double)
            ReadTableQS<double>(input_file, outfile_name);
        else
            ReadTableQS<float>(input_file, outfile_name);
    else
        AbortWithMessage("'--table' must be eiter 'QS' or 'BW'");
}

template<typename RealType>
void ReadTableBW (const string& input_file, const string& outfile_name)
{
    cout << "    Reading BW table " <<
        (is_same<RealType, double>::value ? "(double "s : "(single "s) << " precision)\n"s;

    auto ifs = ifstream(input_file, std::ios::binary);
    auto raw_data = vector<char>(istreambuf_iterator<char>(ifs), {});
    ifs.close();

    auto raw_iter = raw_data.begin();
    const auto size_first = pxr_sr::get_out<uint64_t>(raw_iter);
    if(size_first <= 0 || size_first >= raw_data.size() )
        AbortWithMessage("Something is wrong with " + input_file);

    const auto raw_dndt_table = vector<char>{
        raw_iter, raw_iter+size_first};

    const auto raw_pair_prod_table = vector<char>{
        raw_iter+size_first, raw_data.end()};

    const auto dndt_table =
        pxr_bw::dndt_lookup_table<RealType, vector<RealType>>{raw_dndt_table};
    const auto pair_prod_table =
        pxr_bw::pair_prod_lookup_table<RealType, vector<RealType>>{raw_pair_prod_table};

    if (!dndt_table.is_init() || !pair_prod_table.is_init())
        AbortWithMessage("Something went wrong with lookup table initialization");

    auto of_dndt = ofstream{outfile_name + "_dndt"};
    of_dndt << std::setprecision(std::numeric_limits<RealType>::digits10 + 1);
    const auto coord_dndt = dndt_table.get_all_coordinates();
    for (const auto& cc : coord_dndt )
        of_dndt << cc << " " << dndt_table.interp(cc) << "\n";
    of_dndt.close();

    auto of_pair = ofstream{outfile_name + "_pair"};
    of_pair << std::setprecision(std::numeric_limits<RealType>::digits10 + 1);
    const auto coord_pair = pair_prod_table.get_all_coordinates();
    for (const auto& cc : coord_pair ){
            of_pair << cc[0] << " " << cc[1]
                << " " << pair_prod_table.interp(cc[0], cc[1]) << "\n";
    }
    of_pair.close();
}

template<typename RealType>
void ReadTableQS (const string& input_file, const string& outfile_name)
{
    cout << "    Reading QS table " <<
        (is_same<RealType, double>::value ? "(double "s : "(single "s) << " precision)\n"s;

    auto ifs = ifstream(input_file, std::ios::binary);
    auto raw_data = vector<char>(istreambuf_iterator<char>(ifs), {});
    ifs.close();

    auto raw_iter = raw_data.begin();
    const auto size_first = pxr_sr::get_out<uint64_t>(raw_iter);
    if(size_first <= 0 || size_first >= raw_data.size() )
        AbortWithMessage("Something is wrong with " + input_file);

    const auto raw_dndt_table = vector<char>{
        raw_iter, raw_iter+size_first};

    const auto raw_phot_em_table = vector<char>{
        raw_iter+size_first, raw_data.end()};

    const auto dndt_table =
        pxr_qs::dndt_lookup_table<RealType, vector<RealType>>{raw_dndt_table};
    const auto phot_em_table =
        pxr_qs::photon_emission_lookup_table<RealType, vector<RealType>>{raw_phot_em_table};

    if (!dndt_table.is_init() || !phot_em_table.is_init())
        AbortWithMessage("Something went wrong with lookup table initialization");

    auto of_dndt = ofstream{outfile_name + "_dndt"};
    of_dndt << std::setprecision(std::numeric_limits<RealType>::digits10 + 1);
    const auto coord_dndt = dndt_table.get_all_coordinates();
    for (const auto& cc : coord_dndt )
        of_dndt << cc << " " << dndt_table.interp(cc) << "\n";
    of_dndt.close();

    auto of_phot_em = ofstream{outfile_name + "_phot_em"};
    of_phot_em << std::setprecision(std::numeric_limits<RealType>::digits10 + 1);
    const auto coord_phot_em = phot_em_table.get_all_coordinates();
    for (const auto& cc : coord_phot_em ){
            of_phot_em << cc[0] << " " << cc[1]
                << " " << phot_em_table.interp(cc[0], cc[1]) << "\n";
    }
    of_phot_em.close();

    return;
}
