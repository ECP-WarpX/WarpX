/* Copyright 2021 Hannah Klion
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef TEMPERATURE_PROPERTIES_H_
#define TEMPERATURE_PROPERTIES_H_

#include <Utils/WarpXUtil.H>

#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

/* Type of temperature initialization. Used by TemperatureProperties and GetTemperature. */
enum TemperatureInitType {TempConstantValue, TempParserFunction};

/**
 * \brief Struct to store temperature properties, for use in momentum initialization.
 *
 * Reads in and stores temperature used to initialize the Maxwell-Boltzmann and Maxwell-Juttner
 * momentum distributions in InjectorMomentum. The information is read from the parameters of
 * the species being initialized, and will be accessed by GetTemperature.
 */
struct TemperatureProperties
{
    /**
     * \brief Read runtime parameters to populate constant or spatially-varying temperature
     * information
     *
     * \param[in] pp: Reference to the parameter parser object for the species being initialized
     */
    TemperatureProperties (amrex::ParmParse& pp);

    /* Type of temperature initialization */
    TemperatureInitType m_type;

    /* Constant temperature value, if m_type == TempConstantValue */
    amrex::Real m_temperature;
    /* Storage of the parser function, if m_type == TempParserFunction */
    std::unique_ptr<amrex::Parser> m_ptr_temperature_parser;
};

#endif
