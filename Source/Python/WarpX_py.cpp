/* Copyright 2019-2022 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: David Grote, Maxence Thevenet, Weiqun Zhang, Roelof Groenewald
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX_py.H"

std::map< std::string, std::function<void()> > warpx_callback_py_map;

void InstallPythonCallback ( std::string name, std::function<void()> callback )
{
    warpx_callback_py_map[name] = callback;
}

bool IsPythonCallbackInstalled ( std::string name )
{
    return (warpx_callback_py_map.count(name) == 1u);
}

// Execute Python callbacks of the type given by the input string
void ExecutePythonCallback ( std::string name )
{
    if ( IsPythonCallbackInstalled(name) ) {
        WARPX_PROFILE("warpx_py_"+name);
        warpx_callback_py_map[name]();
    }
}

void ClearPythonCallback ( std::string name )
{
    warpx_callback_py_map.erase(name);
}
