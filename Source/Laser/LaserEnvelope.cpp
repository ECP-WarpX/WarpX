#include "LaserEnvelope.H"

using namespace amrex;

void WarpX::AllocateLaserEnvelope (int lev, const BoxArray& ba, const DistributionMapping& dm, const IntVect& ngA, const bool aux_is_nodal)
{
    A_nodal_flag = IntVect::TheNodeVector();
    A_ncomps = ptr.2
    AllocInitMultiFab(rho_fp[lev], amrex::convert(ba, A_nodal_flag), dm, A_ncomps, ngA, tag("A_fp"), 0.0_rt);
}