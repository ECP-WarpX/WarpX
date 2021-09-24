/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef TEMPERATURE_PROPERTIES_H_
#define TEMPERATURE_PROPERTIES_H_

//#include <Utils/WarpXUtil.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>


enum TemperatureInitType {ConstantValue, ParserFunction};

/*
 * Struct to store temperature properties, for use in momentum initialization.
 */
struct TemperatureProperties
{
    TemperatureProperties (amrex::ParmParse& pp);

    TemperatureInitType m_type;

    amrex::Real m_temperature;
    std::unique_ptr<amrex::Parser> m_ptr_temperature_parser;
};

#endif
