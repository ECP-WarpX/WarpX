/* Copyright 2022 Andrew Myers, Burlen Loring, Luca Fedeli
 * Maxence Thevenet, Remi Lehe, Revathi Jambunathan
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "IntervalsParser.H"

#include "ParserUtils.H"
#include "Utils/Strings/StringUtils.H"
#include "Utils/TextMsg.H"

#include <AMReX_Utility.H>

#include <algorithm>

utils::parser::SliceParser::SliceParser (const std::string& instr)
{
    namespace utils_str = utils::strings;

    // split string and trim whitespaces
    auto insplit = utils_str::split<std::vector<std::string>>(instr, m_separator, true);

    if(insplit.size() == 1){ // no colon in input string. The input is the period.
        m_period = parseStringtoInt(insplit[0], "interval period");}
    else if(insplit.size() == 2) // 1 colon in input string. The input is start:stop
    {
        if (!insplit[0].empty()){
            m_start = parseStringtoInt(insplit[0], "interval start");}
        if (!insplit[1].empty()){
            m_stop = parseStringtoInt(insplit[1], "interval stop");}
    }
    else // 2 colons in input string. The input is start:stop:period
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            insplit.size() == 3,
            instr + "' is not a valid syntax for a slice.");
        if (!insplit[0].empty()){
            m_start = parseStringtoInt(insplit[0], "interval start");}
        if (!insplit[1].empty()){
            m_stop = parseStringtoInt(insplit[1], "interval stop");}
        if (!insplit[2].empty()){
            m_period = parseStringtoInt(insplit[2], "interval period");}
    }
}


bool utils::parser::SliceParser::contains (const int n) const
{
    if (m_period <= 0) {return false;}
    return (n - m_start) % m_period == 0 && n >= m_start && n <= m_stop;
}


int utils::parser::SliceParser::nextContains (const int n) const
{
    if (m_period <= 0) {return std::numeric_limits<int>::max();}
    int next = m_start;
    if (n >= m_start) {next = ((n-m_start)/m_period + 1)*m_period+m_start;}
    if (next > m_stop) {next = std::numeric_limits<int>::max();}
    return next;
}


int utils::parser::SliceParser::previousContains (const int n) const
{
    if (m_period <= 0) {return false;}
    int previous = ((std::min(n-1,m_stop)-m_start)/m_period)*m_period+m_start;
    if ((n < m_start) || (previous < 0)) {previous = 0;}
    return previous;
}


int utils::parser::SliceParser::getPeriod () const {return m_period;}


int utils::parser::SliceParser::getStart () const {return m_start;}


int utils::parser::SliceParser::getStop () const {return m_stop;}


utils::parser::IntervalsParser::IntervalsParser (
    const std::vector<std::string>& instr_vec)
{
    namespace utils_str = utils::strings;

    std::string inconcatenated;
    for (const auto& instr_element : instr_vec) inconcatenated +=instr_element;

    auto insplit = utils_str::split<std::vector<std::string>>(inconcatenated, m_separator);

    for(const auto& inslc : insplit)
    {
        SliceParser temp_slice(inslc);
        m_slices.push_back(temp_slice);
        if ((temp_slice.getPeriod() > 0) &&
               (temp_slice.getStop() >= temp_slice.getStart())) m_activated = true;
    }
}


bool utils::parser::IntervalsParser::contains (const int n) const
{
    return std::any_of(m_slices.begin(), m_slices.end(),
        [&](const auto& slice){return slice.contains(n);});
}


int utils::parser::IntervalsParser::nextContains (const int n) const
{
    int next = std::numeric_limits<int>::max();
    for(const auto& slice: m_slices){
        next = std::min(slice.nextContains(n),next);
    }
    return next;
}


int utils::parser::IntervalsParser::previousContains (const int n) const
{
    int previous = 0;
    for(const auto& slice: m_slices){
        previous = std::max(slice.previousContains(n),previous);
    }
    return previous;
}


int utils::parser::IntervalsParser::previousContainsInclusive (
    const int n) const
{
    if (contains(n)){return n;}
    else {return previousContains(n);}
}


int utils::parser::IntervalsParser::localPeriod (const int n) const
{
    return nextContains(n) - previousContainsInclusive(n);
}


bool utils::parser::IntervalsParser::isActivated () const
{
    return m_activated;
}
