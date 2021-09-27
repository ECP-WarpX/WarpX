/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "GetTemperature.H"

GetTemperature::GetTemperature (TemperatureProperties const& temp) noexcept {
    m_type = temp.m_type;
    if (m_type == ConstantValue) {
        m_temperature = temp.m_temperature;
    }
    else if (m_type == ParserFunction) {
        m_temperature_parser = temp.m_ptr_temperature_parser->compile<3>();
    }
}

amrex::Real GetTemperature::operator() (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
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
