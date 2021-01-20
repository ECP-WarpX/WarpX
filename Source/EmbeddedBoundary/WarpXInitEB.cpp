#include "WarpX.H"

void
WarpX::InitEB ()
{
#ifdef AMREX_USE_EB
    BL_PROFILE("InitEB");

    amrex::ParmParse pp("eb2");
    if (!pp.contains("geom_type")) {
        pp.add("geom_type", "all_regular"); // use all_regrlar by default
    }
    amrex::EB2::Build(Geom(maxLevel()), maxLevel(), maxLevel());
#endif
}
