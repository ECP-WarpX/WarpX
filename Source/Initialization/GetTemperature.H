/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

enum TemperatureInitType { Constant, Parser };

/*
 * Class to store initial temperature for Maxwell-Boltzmann and Maxwell-Juttner
 * particle momentum distributions.
 */
struct GetTemperature
{
    TemperatureInitType m_type;

    amrex::Real m_temperature;
    amrex::ParserExecutor<3> m_temp_parser;

    amrex::Real operator () (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        if (m_type == Constant)
        {
            return m_temperature;
        }
        else if (m_type == Parser) {
            return m_temp_parser(x, y, z);
        }
        else
        {
            amrex::Abort("TemperatureInitType not known! \n");
        }
    }
};
