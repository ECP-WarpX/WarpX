#ifndef WARPX_DIAGNOSTICS_H_
#define WARPX_DIAGNOSTICS_H_

#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

class Diagnostics
{
public:
    Diagnostics (int i, std::string name);
    void Filter ();
    virtual void Compute () = 0;
    virtual void PackFields ();
    void Flush ();
    virtual void FlushRaw ();
    void InitData();
private:
    std::string diag_name;
    int diag_index;
    int const ncomp = 9; // E, B, J
    amrex::Vector< std::string > varnames {"Ex", "Ey", "Ez", "Bx", "By", "Bz", "jx", "jy", "jz"};
    // levels    < Ex Ey Ez Bx etc.                 >
    amrex::Vector< amrex::Vector <amrex::MultiFab*> > allfields;
    amrex::Vector< amrex::MultiFab > mf_avg;
    amrex::Vector< amrex::Geometry > output_geom;
    amrex::Vector<const amrex::MultiFab*>& output_mf;
    int finest_level;
};

#endif // WARPX_DIAGNOSTICS_H_
