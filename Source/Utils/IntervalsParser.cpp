#include "IntervalsParser.H"
#include "TextMsg.H"
#include "WarpXUtil.H"

#include <AMReX_Utility.H>

#include <algorithm>
#include <memory>

SliceParser::SliceParser (const std::string& instr, const bool isBTD)
{
    m_isBTD = isBTD;
    // split string and trim whitespaces
    auto insplit = WarpXUtilStr::split<std::vector<std::string>>(instr, m_separator, true);

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

bool SliceParser::contains (const int n) const
{
    if (m_period <= 0) {return false;}
    return (n - m_start) % m_period == 0 && n >= m_start && n <= m_stop;
}

int SliceParser::nextContains (const int n) const
{
    if (m_period <= 0) {return std::numeric_limits<int>::max();}
    int next = m_start;
    if (n >= m_start) {next = ((n-m_start)/m_period + 1)*m_period+m_start;}
    if (next > m_stop) {next = std::numeric_limits<int>::max();}
    return next;
}

int SliceParser::previousContains (const int n) const
{
    if (m_period <= 0) {return false;}
    int previous = ((std::min(n-1,m_stop)-m_start)/m_period)*m_period+m_start;
    if ((n < m_start) || (previous < 0)) {previous = 0;}
    return previous;
}

int SliceParser::getPeriod () const {return m_period;}

int SliceParser::getStart () const {return m_start;}

int SliceParser::getStop () const {return m_stop;}

int SliceParser::numContained () const {return (m_stop - m_start) / m_period + 1;}

IntervalsParser::IntervalsParser (const std::vector<std::string>& instr_vec)
{
    std::string inconcatenated;
    for (const auto& instr_element : instr_vec) inconcatenated +=instr_element;

    auto insplit = WarpXUtilStr::split<std::vector<std::string>>(inconcatenated, m_separator);

    for(const auto& inslc : insplit)
    {
        SliceParser temp_slice(inslc);
        m_slices.push_back(temp_slice);
        if ((temp_slice.getPeriod() > 0) &&
               (temp_slice.getStop() >= temp_slice.getStart())) m_activated = true;
    }
}

bool IntervalsParser::contains (const int n) const
{
    return std::any_of(m_slices.begin(), m_slices.end(),
        [&](const auto& slice){return slice.contains(n);});
}

int IntervalsParser::nextContains (const int n) const
{
    int next = std::numeric_limits<int>::max();
    for(const auto& slice: m_slices){
        next = std::min(slice.nextContains(n),next);
    }
    return next;
}

int IntervalsParser::previousContains (const int n) const
{
    int previous = 0;
    for(const auto& slice: m_slices){
        previous = std::max(slice.previousContains(n),previous);
    }
    return previous;
}

int IntervalsParser::previousContainsInclusive (const int n) const
{
    if (contains(n)){return n;}
    else {return previousContains(n);}
}

int IntervalsParser::localPeriod (const int n) const
{
    return nextContains(n) - previousContainsInclusive(n);
}

bool IntervalsParser::isActivated () const {return m_activated;}

BTDIntervalsParser::BTDIntervalsParser (const std::vector<std::string>& instr_vec)
{
    std::string inconcatenated;
    for (const auto& instr_element : instr_vec) inconcatenated +=instr_element;

    auto insplit = WarpXUtilStr::split<std::vector<std::string>>(inconcatenated, m_separator);

    std::cout << "BTD intervals:\n";
    for(const auto& inslc : insplit)
    {
        bool isBTD = true;
        SliceParser temp_slice(inslc, isBTD);
        ///
        std::cout << "slice start: " << temp_slice.getStart() << "\n";
        std::cout << "slice stop: " << temp_slice.getStop() << "\n";
        std::cout << "slice period: " << temp_slice.getPeriod() << "\n";
        std::cout << "slice size: " << temp_slice.numContained() << "\n";
        ///
        m_slices.push_back(temp_slice);
        if ((temp_slice.getPeriod() > 0) &&
               (temp_slice.getStop() >= temp_slice.getStart())) m_activated = true;
    }
    int n_slices = m_slices.size();
    std::cout << "n slices = " << n_slices << "\n";
    m_slice_starting_i_buffer = std::vector<int> (n_slices);
    for (int ii = 0; ii < n_slices-1; ii++)
    {
        m_slice_starting_i_buffer[ii+1] = m_slice_starting_i_buffer[ii] + m_slices[ii].numContained();
    }
    for (int ii = 0; ii < n_slices; ii++)
    {
        std::cout << "slice ii contains i_buffers starting with " << m_slice_starting_i_buffer[ii] << "\n";
    }
    m_n_snapshots = 0;
    for (auto& slice : m_slices)
    {
        m_n_snapshots += slice.numContained();
    }
}

int BTDIntervalsParser::NumSnapshots () { return m_n_snapshots; }

int BTDIntervalsParser::GetSliceIndex(int i_buffer)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(0<=i_buffer < m_n_snapshots,
        "invalid i_buffer=" + std::to_string(i_buffer) + " submitted.  Require 0<= i_buffer < " + std::to_string(m_n_snapshots));
    int slice_index = 0;
    while (i_buffer >= m_slice_starting_i_buffer[slice_index+1]
            && static_cast<unsigned long>(slice_index) < m_slices.size()-1)
    {
        slice_index++;
    }
    return slice_index;
}

int BTDIntervalsParser::GetBTDIteration(int i_buffer)
{
    const auto slice_ind = GetSliceIndex(i_buffer);
    auto slice_start_ind = m_slice_starting_i_buffer[slice_ind];
    auto slice = m_slices[slice_ind];
    auto slice_start = slice.getStart();
    auto slice_period = slice.getPeriod();
    int BTDiteration = slice_start + (i_buffer - slice_start_ind) * slice_period;
    return BTDiteration;
}
