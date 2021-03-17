#include "WarpX.H"

#include <AMReX_Config.H>
#ifdef AMREX_USE_EB
#   include <AMReX_EB2.H>
#   include <AMReX_ParmParse.H>
#endif


void
WarpX::InitEB ()
{
#ifdef AMREX_USE_EB
    BL_PROFILE("InitEB");

    amrex::ParmParse pp_eb2("eb2");
    if (!pp_eb2.contains("geom_type")) {
        pp_eb2.add("geom_type", "all_regular"); // use all_regular by default
    }
    amrex::EB2::Build(Geom(maxLevel()), maxLevel(), maxLevel());
#endif
}
