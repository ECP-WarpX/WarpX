/* Copyright 2019 David Grote, Maxence Thevenet, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX_py.H"

std::map< std::string, WARPX_CALLBACK_PY_FUNC_0 > warpx_callback_py_map;

// Execute Python callbacks of the type given by the input string
void ExecutePythonCallback ( std::string name )
{
    if (warpx_callback_py_map.count(name) == 1u) {
        WARPX_PROFILE("warpx_py_"+name);
        warpx_callback_py_map[name]();
    }
}
