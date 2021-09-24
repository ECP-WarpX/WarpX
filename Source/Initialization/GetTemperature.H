/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef GET_TEMPERATURE_H_
#define GET_TEMPERATURE_H_

#include <TemperatureProperties.h>

/*
 * Functor to return temperature at a point on the grid, either constant (m_temperature)
 * or a spatially varying value computed using the parser function (m_temperature_parser)
 */
struct GetTemperature
{
    TemperatureInitType m_type;

    amrex::Real m_temperature;
    amrex::ParserExecutor<3> m_temperature_parser;

    GetTemperature (InitTemperature const& temp) noexcept
    {
        m_type = temp.m_type;
        if (m_type == ConstantValue) {
            m_temperature = temp.m_temperature;
        }
        else if (m_type == ParserFunction) {
            m_temperature_parser = temp.m_ptr_temperature_parser->compile<3>();
        }
    }

    AMREX_GPU_HOST_DEVICE
    amrex::Real operator() (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (m_type)
        {
            case (ConstantValue):
            {
                return m_temperature;
            }
            case (ParserFunction):
            {
                return m_temperature_parser(x,y,z);
            }
            default:
            {
                amrex::Abort("Get initial temperature: unknown type");
                return 0.0;
            }
        }

    }
};
#endif
