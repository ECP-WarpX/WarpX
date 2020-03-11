#ifndef WARPX_MULTIDIAGNOSTICS_H_
#define WARPX_MULTIDIAGNOSTICS_H_

#include "FullDiagnostics.H"

class MultiDiagnostics
{
public:
    MultiDiagnostics ();
    void ReadParameters ();
    void InitData ();
    void FilterComputePackFlush ();
private:
    amrex::Vector<std::unique_ptr<Diagnostics> > alldiags;
    int ndiags = 0;
    std::vector<std::string> diags_names;
    enum struct DiagTypes {Full};
    std::vector<DiagTypes> diags_types;
};

#endif // WARPX_MULTIDIAGNOSTICS_H_
