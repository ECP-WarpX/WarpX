#include <AMReX_RealVect.H>

class RealMatrix : public std::array<amrex::RealVect, AMREX_SPACEDIM>
{
public:

    RealMatrix
    Inverse ();

    amrex::RealVect
    Dot (amrex::RealVect x);

    RealMatrix
    Dot (RealMatrix B);

};