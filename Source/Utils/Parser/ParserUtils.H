/* Copyright 2019-2022 Andrew Myers, Burlen Loring, Luca Fedeli
 * Maxence Thevenet, Remi Lehe, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_UTILS_PARSER_H_
#define WARPX_UTILS_PARSER_H_

namespace utils::parser
{

}

/**
* \brief Do a safe cast of a real to an int
* This ensures that the float value is within the range of ints and if not,
* raises an exception.
*
* \param x Real value to cast
* \param real_name String, the name of the variable being casted to use in the error message
*/
int
safeCastToInt(amrex::Real x, const std::string& real_name);

/**
* \brief Initialize an amrex::Parser object from a string containing a math expression
*
* \param parse_function String to read to initialize the parser.
* \param varnames A list of predefined independent variables
*/
amrex::Parser makeParser (std::string const& parse_function, amrex::Vector<std::string> const& varnames);

/** Parse a string and return as a real
 *
 * In case the string cannot be interpreted as Real,
 * this function ... <throws an exception? aborts with error message?>
 *
 * \param str The string to be parsed
 * \return representation as real
 */
double parseStringtoReal(std::string str);

/** Parse a string and return an int
 *
 * In case the string cannot be interpreted as Real,
 * this function ... <throws an exception? aborts with error message?>
  *
 * \param str The string to be parsed
 * \param name For integers, the name, to be used in error messages
 * \return rounded closest integer
 */
int parseStringtoInt(std::string str, std::string name);

template <int N>
amrex::ParserExecutor<N> compileParser (amrex::Parser const* parser)
{
    if (parser) {
        return parser->compile<N>();
    } else {
        return amrex::ParserExecutor<N>{};
    }
}

/** Similar to amrex::ParmParse::query, but also supports math expressions for the value.
 *
 * amrex::ParmParse::query reads a name and a value from the input file. This function does the
 * same, and applies the amrex::Parser to the value, so the user has the choice to specify a value or
 * a math expression (including user-defined constants).
 * Works for amrex::Real numbers and integers.
 *
 * \param[in] a_pp amrex::ParmParse object
 * \param[in] str name of the parameter to read
 * \param[out] val where the value queried and parsed is stored, either a scalar or vector
 */
template <typename T>
int queryWithParser (const amrex::ParmParse& a_pp, char const * const str, T& val)
{
    // call amrex::ParmParse::query, check if the user specified str.
    std::string tmp_str;
    int is_specified = a_pp.query(str, tmp_str);
    if (is_specified)
    {
        // If so, create a parser object and apply it to the value provided by the user.
        std::string str_val;
        Store_parserString(a_pp, str, str_val);

        auto parser = makeParser(str_val, {});

        if (std::is_same<T, int>::value) {

            val = safeCastToInt(std::round(parser.compileHost<0>()()), str);
        }
        else {
            val = static_cast<T>(parser.compileHost<0>()());
        }
    }
    // return the same output as amrex::ParmParse::query
    return is_specified;
}

template <typename T>
int queryArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<T>& val)
{
    // call amrex::ParmParse::query, check if the user specified str.
    std::vector<std::string> tmp_str_arr;
    int is_specified = a_pp.queryarr(str, tmp_str_arr);
    if (is_specified)
    {
        // If so, create parser objects and apply them to the values provided by the user.
        int const n = static_cast<int>(tmp_str_arr.size());
        val.resize(n);
        for (int i=0 ; i < n ; i++) {
            auto parser = makeParser(tmp_str_arr[i], {});
            if (std::is_same<T, int>::value) {
                val[i] = safeCastToInt(std::round(parser.compileHost<0>()()), str);
            }
            else {
                val[i] = static_cast<T>(parser.compileHost<0>()());
            }
        }
    }
    // return the same output as amrex::ParmParse::query
    return is_specified;
}

/** Similar to amrex::ParmParse::query, but also supports math expressions for the value.
 *
 * amrex::ParmParse::query reads a name and a value from the input file. This function does the
 * same, and applies the amrex::Parser to the value, so the user has the choice to specify a value or
 * a math expression (including user-defined constants).
 * Works for amrex::Real numbers and integers.
 *
 * \param[in] a_pp amrex::ParmParse object
 * \param[in] str name of the parameter to read
 * \param[out] val where the value queried and parsed is stored, either a scalar or vector
 * \param[in] start_ix start index in the list of inputs values (optional with arrays, default is
 * amrex::ParmParse::FIRST for starting with the first input value)
 * \param[in] num_val number of input values to use (optional with arrays, default is
 * amrex::ParmParse::LAST for reading until the last input value)
 */
template <typename T>
int queryArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<T>& val,
                        const int start_ix, const int num_val)
{
    // call amrex::ParmParse::query, check if the user specified str.
    std::vector<std::string> tmp_str_arr;
    int is_specified = a_pp.queryarr(str, tmp_str_arr, start_ix, num_val);
    if (is_specified)
    {
        // If so, create parser objects and apply them to the values provided by the user.
        int const n = static_cast<int>(tmp_str_arr.size());
        val.resize(n);
        for (int i=0 ; i < n ; i++) {
            auto parser = makeParser(tmp_str_arr[i], {});
            if (std::is_same<T, int>::value) {
                val[i] = safeCastToInt(std::round(parser.compileHost<0>()()), str);
            }
            else {
                val[i] = static_cast<T>(parser.compileHost<0>()());
            }
        }
    }
    // return the same output as amrex::ParmParse::query
    return is_specified;
}

/** Similar to amrex::ParmParse::get, but also supports math expressions for the value.
 *
 * amrex::ParmParse::get reads a name and a value from the input file. This function does the
 * same, and applies the Parser to the value, so the user has the choice to specify a value or
 * a math expression (including user-defined constants).
 * Works for amrex::Real numbers and integers.
 *
 * \param[in] a_pp amrex::ParmParse object
 * \param[in] str name of the parameter to read
 * \param[out] val where the value queried and parsed is stored
 */
template <typename T>
void getWithParser (const amrex::ParmParse& a_pp, char const * const str, T& val)
{
    // If so, create a parser object and apply it to the value provided by the user.
    std::string str_val;
    Store_parserString(a_pp, str, str_val);

    auto parser = makeParser(str_val, {});
    if (std::is_same<T, int>::value) {
        val = safeCastToInt(std::round(parser.compileHost<0>()()), str);
    }
    else {
        val = static_cast<T>(parser.compileHost<0>()());
    }
}

template <typename T>
void getArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<T>& val)
{
    // Create parser objects and apply them to the values provided by the user.
    std::vector<std::string> tmp_str_arr;
    a_pp.getarr(str, tmp_str_arr);

    int const n = static_cast<int>(tmp_str_arr.size());
    val.resize(n);
    for (int i=0 ; i < n ; i++) {
        auto parser = makeParser(tmp_str_arr[i], {});
        if (std::is_same<T, int>::value) {
            val[i] = safeCastToInt(std::round(parser.compileHost<0>()()), str);
        }
        else {
            val[i] = static_cast<T>(parser.compileHost<0>()());
        }
    }
}

/** Similar to amrex::ParmParse::get, but also supports math expressions for the value.
 *
 * amrex::ParmParse::get reads a name and a value from the input file. This function does the
 * same, and applies the Parser to the value, so the user has the choice to specify a value or
 * a math expression (including user-defined constants).
 * Works for amrex::Real numbers and integers.
 *
 * \param[in] a_pp amrex::ParmParse object
 * \param[in] str name of the parameter to read
 * \param[out] val where the value queried and parsed is stored
 * \param[in] start_ix start index in the list of inputs values (optional with arrays, default is
 * amrex::ParmParse::FIRST for starting with the first input value)
 * \param[in] num_val number of input values to use (optional with arrays, default is
 * amrex::ParmParse::LAST for reading until the last input value)
 */
template <typename T>
void getArrWithParser (const amrex::ParmParse& a_pp, char const * const str, std::vector<T>& val,
                       const int start_ix, const int num_val)
{
    // Create parser objects and apply them to the values provided by the user.
    std::vector<std::string> tmp_str_arr;
    a_pp.getarr(str, tmp_str_arr, start_ix, num_val);

    int const n = static_cast<int>(tmp_str_arr.size());
    val.resize(n);
    for (int i=0 ; i < n ; i++) {
        auto parser = makeParser(tmp_str_arr[i], {});
        if (std::is_same<T, int>::value) {
            val[i] = safeCastToInt(std::round(parser.compileHost<0>()()), str);
        }
        else {
            val[i] = static_cast<T>(parser.compileHost<0>()());
        }
    }
}