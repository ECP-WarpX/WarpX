/* Copyright 2022 The ABLASTR Community
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: Neil Zaim, Luca Fedeli, Axel Huebl
 */
#include "ablastr/utils/IntervalsParser.H"

#include "ablastr/utils/text/StringUtils.H"
//include "ParserUtils.H"
#include "Utils/TextMsg.H"

#include <AMReX_Utility.H>

#include <algorithm>


ablastr::utils::parser::SliceParser::SliceParser (
    const std::string& instr,
    const bool isBTD
):
    m_isBTD{isBTD}
{
    // split string and trim whitespaces
    auto insplit = ablastr::utils::text::split_string<std::vector<std::string>>(
        instr, m_separator, true);

    if(insplit.size() == 1){ // no colon in input string. The input is the period.
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_isBTD, "must specify interval stop for BTD");
        m_period = parseStringtoInt(insplit[0], "interval period");}
    else if(insplit.size() == 2) // 1 colon in input string. The input is start:stop
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_isBTD || !insplit[1].empty(), "must specify interval stop for BTD");
        if (!insplit[0].empty()){
            m_start = parseStringtoInt(insplit[0], "interval start");}
        if (!insplit[1].empty()){
            m_stop = parseStringtoInt(insplit[1], "interval stop");}
    }
    else // 2 colons in input string. The input is start:stop:period
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(!m_isBTD || !insplit[1].empty(), "must specify interval stop for BTD");
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


bool ablastr::utils::parser::SliceParser::contains (const int n) const
{
    if (m_period <= 0) {return false;}
    return (n - m_start) % m_period == 0 && n >= m_start && n <= m_stop;
}


int ablastr::utils::parser::SliceParser::nextContains (const int n) const
{
    if (m_period <= 0) {return std::numeric_limits<int>::max();}
    int next = m_start;
    if (n >= m_start) {next = ((n-m_start)/m_period + 1)*m_period+m_start;}
    if (next > m_stop) {next = std::numeric_limits<int>::max();}
    return next;
}


int ablastr::utils::parser::SliceParser::previousContains (const int n) const
{
    if (m_period <= 0) {return false;}
    int previous = ((std::min(n-1,m_stop)-m_start)/m_period)*m_period+m_start;
    if ((n < m_start) || (previous < 0)) {previous = 0;}
    return previous;
}


int ablastr::utils::parser::SliceParser::getPeriod () const
{
    return m_period;
}


int ablastr::utils::parser::SliceParser::getStart () const
{
    return m_start;
}


int ablastr::utils::parser::SliceParser::getStop () const
{
    return m_stop;
}


int
ablastr::utils::parser::SliceParser::numContained () const
{
    return (m_stop - m_start) / m_period + 1;
}

ablastr::utils::parser::IntervalsParser::IntervalsParser (
    const std::vector<std::string>& instr_vec
)
{
    std::string inconcatenated;
    for (const auto& instr_element : instr_vec) { inconcatenated +=instr_element; }

    auto insplit = ablastr::utils::text::split_string<std::vector<std::string>>(
        inconcatenated, m_separator);

    for(const auto& inslc : insplit)
    {
        const SliceParser temp_slice(inslc);
        m_slices.push_back(temp_slice);
        if ((temp_slice.getPeriod() > 0) &&
               (temp_slice.getStop() >= temp_slice.getStart())) { m_activated = true; }
    }
}


bool
ablastr::utils::parser::IntervalsParser::contains (const int n) const
{
    return std::any_of(m_slices.begin(), m_slices.end(),
        [&](const auto& slice){return slice.contains(n);});
}


int
ablastr::utils::parser::IntervalsParser::nextContains (const int n) const
{
    int next = std::numeric_limits<int>::max();
    for(const auto& slice: m_slices){
        next = std::min(slice.nextContains(n),next);
    }
    return next;
}


int
ablastr::utils::parser::IntervalsParser::previousContains (const int n) const
{
    int previous = 0;
    for(const auto& slice: m_slices){
        previous = std::max(slice.previousContains(n),previous);
    }
    return previous;
}


int
ablastr::utils::parser::IntervalsParser::previousContainsInclusive (
    const int n) const
{
    if (contains(n)){return n;}
    else {return previousContains(n);}
}


int
ablastr::utils::parser::IntervalsParser::localPeriod (const int n) const
{
    return nextContains(n) - previousContainsInclusive(n);
}


bool
ablastr::utils::parser::IntervalsParser::isActivated () const
{
    return m_activated;
}
