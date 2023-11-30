/* Copyright 2023 Luca Fedeli
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef EXTERNAL_FIELD_H_
#define EXTERNAL_FIELD_H_

#include "ExternalField_fwd.H"

#include <AMReX_Array.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

#include <memory>
#include <string>

enum class ExternalFieldType
{
    default_zero,
    constant,
    parse_ext_grid_function,
    read_from_file
};

/**
 * \brief Struct to store data related to external electromagnetic fields
 * (flags, field values, and field parsers)
 */
struct ExternalFieldParams
{
    //! Initial electric field on the grid
    amrex::GpuArray<amrex::Real,3> E_external_grid = {0,0,0};
    //! Initial magnetic field on the grid
    amrex::GpuArray<amrex::Real,3> B_external_grid = {0,0,0};

    //! Initialization type for external magnetic field on the grid
    ExternalFieldType B_ext_grid_type = ExternalFieldType::default_zero;
    //! Initialization type for external electric field on the grid
    ExternalFieldType E_ext_grid_type = ExternalFieldType::default_zero;

    //! User-defined parser to initialize x-component of the magnetic field on the grid
    std::unique_ptr<amrex::Parser> Bxfield_parser;
    //! User-defined parser to initialize y-component of the magnetic field on the grid
    std::unique_ptr<amrex::Parser> Byfield_parser;
    //! User-defined parser to initialize z-component of the magnetic field on the grid
    std::unique_ptr<amrex::Parser> Bzfield_parser;
    //! User-defined parser to initialize x-component of the electric field on the grid
    std::unique_ptr<amrex::Parser> Exfield_parser;
    //! User-defined parser to initialize y-component of the electric field on the grid
    std::unique_ptr<amrex::Parser> Eyfield_parser;
    //! User-defined parser to initialize z-component of the electric field on the grid
    std::unique_ptr<amrex::Parser> Ezfield_parser;

     //! Path of the file where external fields are stored
    std::string external_fields_path;
};

/**
 * \brief This function reads the ExternalField data, initializes an ExternalFieldParams struct and returns
 * a unique pointer containing the initialized struct. "pp_warpx" must point at
 * the "warpx" parameter group in the inputfile.
 *
 */
[[nodiscard]] std::unique_ptr<ExternalFieldParams>
ReadExternalFieldParams(const amrex::ParmParse& pp_warpx);

#endif //EXTERNAL_FIELD_H_
