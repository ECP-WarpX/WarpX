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
private:
    std::string diag_name;
    int diag_index;
    // levels    < Ex Ey Ez Bx etc.                 >
    amrex::Vector< amrex::Vector <amrex::MultiFab*> > allfields;
    amrex::Vector< amrex::MultiFab > mf_avg;
};

#endif // WARPX_DIAGNOSTICS_H_
