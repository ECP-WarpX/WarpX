
#include "Diagnostics.H"
#include "WarpX.H"

using namespace amrex;

Diagnostics::Diagnostics (int i, std::string name)
    : diag_index(i), diag_name(name)
{
    allfields.resize(1);
    allfields[0].resize(1);
    // auto & warpx = WarpX::GetInstance();
    // allfields[0] = warpx.Efield_fp[0][0].get();
    // allfields[0][0] = warpx.getEfieldp(0,0);
    // MultiFab* foo = warpx.getEfieldp(0,0);
    // my_Ex = warpx.getEfield(0,0);
    // allfields[0][0].reset();
    // allfields[0][0] = foo;
}

void
Diagnostics::Filter () {}

void
Diagnostics::PackFields () {}

void
Diagnostics::Flush () {}

void
Diagnostics::FlushRaw () {}
