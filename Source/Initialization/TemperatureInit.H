/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef TEMPERATURE_INIT_H_
#define TEMPERATURE_INIT_H_

#include <AMReX_REAL.H>
#include <AMReX_Parser.H>

/* Class that stores and returns a constant temperature */
struct TemperatureInitConstant
{
    TemperatureInitConstant (amrex::Real a_temperature) noexcept
        : m_temperature(a_temperature) {}

    amrex::Real getTemperature(amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return m_temperature;
    }
private:
    amrex::Real m_temperature;
};

struct TemperatureInitParser
{
    TemperatureInitParser (amrex::ParserExecutor<3> const& a_temperature_parser) noexcept
        : m_temperature_parser(a_temperature_parser) {}

    amrex::Real getTemperature(amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return m_temperature_parser(x, y, z);
    }
private:
    amrex::ParserExecutor<3> m_temperature_parser;
};

/*
 * Base struct for temperature initialization. TemperatureInit contains a union (called Object) 
 * that holds any one instance of:
 * - TemperatureInitConstant      : generates constant temperature
 * - TemperatureInitParser        : generates temperature from parser
 */
struct TemperatureInit
{
    TemperatureInit (TemperatureInitConstant* t, amrex::Real a_temperature) : 
        type(Type::constant), object(t, a_temperature) {};

    TemperatureInit (TemperatureInitParser* t,
            amrex::ParserExecutor<3> const& a_temperature_parser) :
        type(Type::parser), object(t, a_temperature_parser) {};

    amrex::Real getTemperature (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (type)
        {
            case Type::constant:
            {
                return object.constant.getTemperature(x,y,z);
            }
            case Type::parser:
            {
                return object.parser.getTemperature(x,y,z);
            }
            default:
            {
                amrex::Abort("Get initial temperature: unknown type");
                return 0.0;
            }
        }

    }

    void clear ();

private:
    enum struct Type {constant, parser};
    Type type;

    union Object {
        Object (TemperatureInitConstant*, amrex::Real a_temperature) noexcept 
            : constant(a_temperature) {}
        Object (TemperatureInitParser*, amrex::ParserExecutor<3> const& a_temperature) noexcept
            : parser(a_temperature) {}
        TemperatureInitConstant constant;
        TemperatureInitParser parser;
    };
    Object object;
};

// In order for GetTemperature to be trivially copyable, its destructor
// must be trivial.  So we have to rely on a custom deleter for unique_ptr.
struct TemperatureInitDeleter {
    void operator () (TemperatureInit* p) const {
        if (p) {
            p->clear();
            delete p;
        }
    }
};
#endif
