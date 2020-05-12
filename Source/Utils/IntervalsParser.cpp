#include "WarpXUtil.H"
#include "IntervalsParser.H"

SliceParser::SliceParser (const std::string& instr)
{
    const std::string assert_msg = "ERROR: '" + instr + "' is not a valid syntax for a slice.";

    // split string and trim whitespaces
    auto insplit = WarpXUtilStr::split<std::vector<std::string>>(instr, m_separator, true);

    if(insplit.size() == 1){ // no colon in input string. The input is the period.
        WarpXUtilMsg::AlwaysAssert(amrex::is_it<int>(insplit[0], m_period),assert_msg);}
    else if(insplit.size() == 2) // 1 colon in input string. The input is start:stop
    {
        if (!insplit[0].empty()){
            WarpXUtilMsg::AlwaysAssert(amrex::is_it<int>(insplit[0], m_start),assert_msg);}
        if (!insplit[1].empty()){
            WarpXUtilMsg::AlwaysAssert(amrex::is_it<int>(insplit[1], m_stop),assert_msg);}
    }
    else // 2 colons in input string. The input is start:stop:period
    {
        WarpXUtilMsg::AlwaysAssert(insplit.size() == 3,assert_msg);
        if (!insplit[0].empty()){
            WarpXUtilMsg::AlwaysAssert(amrex::is_it<int>(insplit[0], m_start),assert_msg);}
        if (!insplit[1].empty()){
            WarpXUtilMsg::AlwaysAssert(amrex::is_it<int>(insplit[1], m_stop),assert_msg);}
        if (!insplit[2].empty()){
            WarpXUtilMsg::AlwaysAssert(amrex::is_it<int>(insplit[2], m_period),assert_msg);}
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

IntervalsParser::IntervalsParser (const std::string& instr)
{
    auto insplit = WarpXUtilStr::split<std::vector<std::string>>(instr, m_separator);

    for(int i=0; i<insplit.size(); i++)
    {
        SliceParser temp_slice(insplit[i]);
        m_slices.push_back(temp_slice);
        if ((temp_slice.getPeriod() > 0) &&
               (temp_slice.getStop() >= temp_slice.getStart())) m_activated = true;
    }
}

bool IntervalsParser::contains (const int n) const
{
    for(int i=0; i<m_slices.size(); i++){
        if (m_slices[i].contains(n)) return true;
    }
    return false;
}

int IntervalsParser::nextContains (const int n) const
{
    int next = std::numeric_limits<int>::max();
    for(int i=0; i<m_slices.size(); i++){
        next = std::min(m_slices[i].nextContains(n),next);
    }
    return next;
}

int IntervalsParser::previousContains (const int n) const
{
    int previous = 0;
    for(int i=0; i<m_slices.size(); i++){
        previous = std::max(m_slices[i].previousContains(n),previous);
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
