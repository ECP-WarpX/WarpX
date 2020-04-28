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
    else{
        return (n - m_start) % m_period == 0 && n >= m_start && n <= m_stop;}
}

int SliceParser::getPeriod () const {return m_period;}

IntervalsParser::IntervalsParser (const std::string& instr)
{
    auto insplit = WarpXUtilStr::split<std::vector<std::string>>(instr, m_separator);

    for(int i=0; i<insplit.size(); i++)
    {
        SliceParser temp_slice(insplit[i]);
        m_slices.push_back(temp_slice);
        if (temp_slice.getPeriod() > 0) m_activated = true;
    }
}

bool IntervalsParser::contains (const int n) const
{
    for(int i=0; i<m_slices.size(); i++){
        if (m_slices[i].contains(n)) return true;
    }
    return false;
}

bool IntervalsParser::isActivated () const {return m_activated;}