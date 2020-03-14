#ifndef WARPX_DIAGNOSTICS_H_
#define WARPX_DIAGNOSTICS_H_

#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

class Diagnostics
{
public:
    Diagnostics (int i, std::string name);
    void FilterAndPack ();
    void Flush ();
    virtual void FlushRaw ();
    void InitData();
private:
    std::string diag_name;
    int diag_index;
    int const ncomp = 9; // E, B, J
    amrex::Vector< std::string > varnames
    {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz"};
    // levels    < Ex Ey Ez Bx etc.                 >
    amrex::Vector< amrex::Vector <amrex::MultiFab*> > allfields;
    amrex::Vector< amrex::MultiFab > mf_avg;
    amrex::Vector< amrex::MultiFab* > mf_out;
    amrex::Vector< amrex::Geometry > output_geom;
    int nlev;
};

#endif // WARPX_DIAGNOSTICS_H_
