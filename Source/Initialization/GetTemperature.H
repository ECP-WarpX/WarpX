/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef GET_TEMPERATURE_H_
#define GET_TEMPERATURE_H_

#include <Initialization/TemperatureProperties.H>

/*
 * Functor to return temperature at a point on the grid, either constant (m_temperature)
 * or a spatially varying value computed using the parser function (m_temperature_parser)
 */
struct GetTemperature
{
    TemperatureInitType m_type;

    amrex::Real m_temperature;
    amrex::ParserExecutor<3> m_temperature_parser;

    GetTemperature (TemperatureProperties const& temp) noexcept;

    AMREX_GPU_HOST_DEVICE
    amrex::Real operator() (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept;
};
#endif
