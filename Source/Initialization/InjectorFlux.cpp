/* Copyright 2023 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "InjectorFlux.H"

using namespace amrex;

void InjectorFlux::clear ()
{
    switch (type)
    {
    case Type::constant:
    case Type::parser:
    {
        break;
    }
    }
}
