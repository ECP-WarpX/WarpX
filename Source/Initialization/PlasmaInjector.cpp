/* Copyright 2019-2020 Andrew Myers, Axel Huebl, Cameron Yang
 * David Grote, Luca Fedeli, Maxence Thevenet
 * Remi Lehe, Revathi Jambunathan, Weiqun Zhang
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "PlasmaInjector.H"

#include "Initialization/GetTemperature.H"
#include "Initialization/GetVelocity.H"
#include "Initialization/InjectorDensity.H"
#include "Initialization/InjectorMomentum.H"
#include "Initialization/InjectorPosition.H"
#include "Particles/SpeciesPhysicalProperties.H"
#include "Utils/Parser/ParserUtils.H"
#include "Utils/TextMsg.H"
#include "Utils/WarpXConst.H"
#include "WarpX.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_Config.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_Print.H>
#include <AMReX_RandomEngine.H>
#include <AMReX_REAL.H>

#include <algorithm>
#include <cctype>
#include <map>
#include <memory>
#include <sstream>
#include <utility>
#include <vector>

namespace {
    void StringParseAbortMessage(const std::string& var,
                                 const std::string& name) {
        std::stringstream stringstream;
        std::string string;
        stringstream << var << " string '" << name << "' not recognized.";
        string = stringstream.str();
        amrex::Abort(Utils::TextMsg::Err(string.c_str()));
    }
}

PlasmaInjector::PlasmaInjector (int ispecies, const std::string& name)
    : species_id(ispecies), species_name(name)
{
    using namespace amrex::literals;

    amrex::ParmParse pp_species_name(species_name);

#ifdef AMREX_USE_GPU
    static_assert(std::is_trivially_copyable<InjectorPosition>::value,
                  "InjectorPosition must be trivially copyable");
    static_assert(std::is_trivially_copyable<InjectorDensity>::value,
                  "InjectorDensity must be trivially copyable");
    static_assert(std::is_trivially_copyable<InjectorMomentum>::value,
                  "InjectorMomentum must be trivially copyable");
#endif

    pp_species_name.query("radially_weighted", radially_weighted);
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(radially_weighted, "ERROR: Only radially_weighted=true is supported");

    // Unlimited boundaries
    xmin = std::numeric_limits<amrex::Real>::lowest();
    ymin = std::numeric_limits<amrex::Real>::lowest();
    zmin = std::numeric_limits<amrex::Real>::lowest();

    xmax = std::numeric_limits<amrex::Real>::max();
    ymax = std::numeric_limits<amrex::Real>::max();
    zmax = std::numeric_limits<amrex::Real>::max();

    // NOTE: When periodic boundaries are used, default injection range is set to mother grid dimensions.
    const amrex::Geometry& geom = WarpX::GetInstance().Geom(0);
    if( geom.isPeriodic(0) ) {
#       ifndef WARPX_DIM_1D_Z
        xmin = geom.ProbLo(0);
        xmax = geom.ProbHi(0);
#       else
        zmin = geom.ProbLo(0);
        zmax = geom.ProbHi(0);
#       endif
    }

#   ifndef WARPX_DIM_1D_Z
    if( geom.isPeriodic(1) ) {
#       ifndef WARPX_DIM_3D
        zmin = geom.ProbLo(1);
        zmax = geom.ProbHi(1);
#       else
        ymin = geom.ProbLo(1);
        ymax = geom.ProbHi(1);
#       endif
    }
#       endif

#   ifdef WARPX_DIM_3D
    if( geom.isPeriodic(2) ) {
        zmin = geom.ProbLo(2);
        zmax = geom.ProbHi(2);
    }
#   endif

    utils::parser::queryWithParser(pp_species_name, "xmin", xmin);
    utils::parser::queryWithParser(pp_species_name, "ymin", ymin);
    utils::parser::queryWithParser(pp_species_name, "zmin", zmin);
    utils::parser::queryWithParser(pp_species_name, "xmax", xmax);
    utils::parser::queryWithParser(pp_species_name, "ymax", ymax);
    utils::parser::queryWithParser(pp_species_name, "zmax", zmax);

    utils::parser::queryWithParser(pp_species_name, "density_min", density_min);
    utils::parser::queryWithParser(pp_species_name, "density_max", density_max);

    std::string physical_species_s;
    bool species_is_specified = pp_species_name.query("species_type", physical_species_s);
    if (species_is_specified){
        const auto physical_species_from_string = species::from_string( physical_species_s );
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(physical_species_from_string,
            physical_species_s + " does not exist!");
        physical_species = physical_species_from_string.value();
        charge = species::get_charge( physical_species );
        mass = species::get_mass( physical_species );
    }

    // Parse injection style
    std::string injection_style = "none";
    pp_species_name.query("injection_style", injection_style);
    std::transform(injection_style.begin(),
                   injection_style.end(),
                   injection_style.begin(),
                   ::tolower);

    // parse charge and mass
    const bool charge_is_specified =
        utils::parser::queryWithParser(pp_species_name, "charge", charge);
    const bool mass_is_specified =
        utils::parser::queryWithParser(pp_species_name, "mass", mass);

    if ( charge_is_specified && species_is_specified ){
        ablastr::warn_manager::WMRecordWarning("Species",
            "Both '" + species_name +  ".charge' and " +
                species_name + ".species_type' are specified.\n" +
                species_name + ".charge' will take precedence.\n");

    }
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        charge_is_specified ||
        species_is_specified ||
        (injection_style == "external_file"),
        "Need to specify at least one of species_type or charge for species '" +
        species_name + "'."
    );

    if ( mass_is_specified && species_is_specified ){
        ablastr::warn_manager::WMRecordWarning("Species",
            "Both '" + species_name +  ".mass' and " +
                species_name + ".species_type' are specified.\n" +
                species_name + ".mass' will take precedence.\n");
    }

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        mass_is_specified ||
        species_is_specified ||
        (injection_style == "external_file"),
        "Need to specify at least one of species_type or mass for species '" +
        species_name + "'."
    );

    num_particles_per_cell_each_dim.assign(3, 0);

    if (injection_style == "none") {
        return;
    } else if (injection_style == "singleparticle") {
        utils::parser::getArrWithParser(
            pp_species_name, "single_particle_pos", single_particle_pos, 0, 3);
        utils::parser::getArrWithParser(
            pp_species_name, "single_particle_vel", single_particle_vel, 0, 3);
        for (auto& x : single_particle_vel) {
            x *= PhysConst::c;
        }
        utils::parser::getWithParser(
            pp_species_name, "single_particle_weight", single_particle_weight);
        add_single_particle = true;
        return;
    } else if (injection_style == "multipleparticles") {
        utils::parser::getArrWithParser(
            pp_species_name, "multiple_particles_pos_x", multiple_particles_pos_x);
        utils::parser::getArrWithParser(
            pp_species_name, "multiple_particles_pos_y", multiple_particles_pos_y);
        utils::parser::getArrWithParser(
            pp_species_name, "multiple_particles_pos_z", multiple_particles_pos_z);
        utils::parser::getArrWithParser(
            pp_species_name, "multiple_particles_vel_x", multiple_particles_vel_x);
        utils::parser::getArrWithParser(
            pp_species_name, "multiple_particles_vel_y", multiple_particles_vel_y);
        utils::parser::getArrWithParser(
            pp_species_name, "multiple_particles_vel_z", multiple_particles_vel_z);
        utils::parser::getArrWithParser(
            pp_species_name, "multiple_particles_weight", multiple_particles_weight);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            ((multiple_particles_pos_x.size() == multiple_particles_pos_y.size()) &&
             (multiple_particles_pos_x.size() == multiple_particles_pos_z.size()) &&
             (multiple_particles_pos_x.size() == multiple_particles_vel_x.size()) &&
             (multiple_particles_pos_x.size() == multiple_particles_vel_y.size()) &&
             (multiple_particles_pos_x.size() == multiple_particles_vel_z.size()) &&
             (multiple_particles_pos_x.size() == multiple_particles_weight.size())),
            "Error: The multiple particles source quantities must all have the same number of elements");
        for (auto& vx : multiple_particles_vel_x) { vx *= PhysConst::c; }
        for (auto& vy : multiple_particles_vel_y) { vy *= PhysConst::c; }
        for (auto& vz : multiple_particles_vel_z) { vz *= PhysConst::c; }
        add_multiple_particles = true;
        return;
    } else if (injection_style == "gaussian_beam") {

        utils::parser::getWithParser(pp_species_name, "x_m", x_m);
        utils::parser::getWithParser(pp_species_name, "y_m", y_m);
        utils::parser::getWithParser(pp_species_name, "z_m", z_m);
        utils::parser::getWithParser(pp_species_name, "x_rms", x_rms);
        utils::parser::getWithParser(pp_species_name, "y_rms", y_rms);
        utils::parser::getWithParser(pp_species_name, "z_rms", z_rms);
        utils::parser::queryWithParser(pp_species_name, "x_cut", x_cut);
        utils::parser::queryWithParser(pp_species_name, "y_cut", y_cut);
        utils::parser::queryWithParser(pp_species_name, "z_cut", z_cut);
        utils::parser::getWithParser(pp_species_name, "q_tot", q_tot);
        utils::parser::getWithParser(pp_species_name, "npart", npart);
        pp_species_name.query("do_symmetrize", do_symmetrize);
        gaussian_beam = true;
        parseMomentum(pp_species_name);
#if defined(WARPX_DIM_XZ)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( y_rms > 0._rt,
        "Error: Gaussian beam y_rms must be strictly greater than 0 in 2D "
        "(it is used when computing the particles' weights from the total beam charge)");
#elif defined(WARPX_DIM_1D_Z)
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( x_rms > 0._rt,
        "Error: Gaussian beam x_rms must be strictly greater than 0 in 1D "
        "(it is used when computing the particles' weights from the total beam charge)");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE( y_rms > 0._rt,
        "Error: Gaussian beam y_rms must be strictly greater than 0 in 1D "
        "(it is used when computing the particles' weights from the total beam charge)");
#endif

    }
    // Depending on injection type at runtime, initialize inj_pos
    // so that inj_pos->getPositionUnitBox calls
    // InjectorPosition[Random or Regular].getPositionUnitBox.
    else if (injection_style == "nrandompercell") {
        utils::parser::getWithParser(
            pp_species_name, "num_particles_per_cell", num_particles_per_cell);
#if WARPX_DIM_RZ
        if (WarpX::n_rz_azimuthal_modes > 1) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            num_particles_per_cell>=2*WarpX::n_rz_azimuthal_modes,
            "Error: For accurate use of WarpX cylindrical gemoetry the number "
            "of particles should be at least two times n_rz_azimuthal_modes "
            "(Please visit PR#765 for more information.)");
        }
#endif
        // Construct InjectorPosition with InjectorPositionRandom.
        h_inj_pos = std::make_unique<InjectorPosition>(
            (InjectorPositionRandom*)nullptr,
            xmin, xmax, ymin, ymax, zmin, zmax);
        parseDensity(pp_species_name);
        parseMomentum(pp_species_name);
    } else if (injection_style == "nfluxpercell") {
        surface_flux = true;
        utils::parser::getWithParser(
            pp_species_name, "num_particles_per_cell", num_particles_per_cell_real);
#ifdef WARPX_DIM_RZ
        if (WarpX::n_rz_azimuthal_modes > 1) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            num_particles_per_cell_real>=2*WarpX::n_rz_azimuthal_modes,
            "Error: For accurate use of WarpX cylindrical geometry the number "
            "of particles should be at least two times n_rz_azimuthal_modes "
            "(Please visit PR#765 for more information.)");
        }
#endif
        utils::parser::getWithParser(
            pp_species_name, "surface_flux_pos", surface_flux_pos);
        utils::parser::queryWithParser(
            pp_species_name, "flux_tmin", flux_tmin);
        utils::parser::queryWithParser(
            pp_species_name, "flux_tmax", flux_tmax);
        std::string flux_normal_axis_string;
        pp_species_name.get("flux_normal_axis", flux_normal_axis_string);
        flux_normal_axis = -1;
#ifdef WARPX_DIM_RZ
        if      (flux_normal_axis_string == "r" || flux_normal_axis_string == "R") {
            flux_normal_axis = 0;
        }
        if      (flux_normal_axis_string == "t" || flux_normal_axis_string == "T") {
            flux_normal_axis = 1;
        }
#else
#    ifndef WARPX_DIM_1D_Z
        if      (flux_normal_axis_string == "x" || flux_normal_axis_string == "X") {
            flux_normal_axis = 0;
        }
#    endif
#endif
#ifdef WARPX_DIM_3D
        if (flux_normal_axis_string == "y" || flux_normal_axis_string == "Y") {
            flux_normal_axis = 1;
        }
#endif
        if (flux_normal_axis_string == "z" || flux_normal_axis_string == "Z") {
            flux_normal_axis = 2;
        }
#ifdef WARPX_DIM_3D
        std::string flux_normal_axis_help = "'x', 'y', or 'z'.";
#else
#    ifdef WARPX_DIM_RZ
        std::string flux_normal_axis_help = "'r' or 'z'.";
#    elif WARPX_DIM_XZ
        std::string flux_normal_axis_help = "'x' or 'z'.";
#    else
        std::string flux_normal_axis_help = "'z'.";
#    endif
#endif
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(flux_normal_axis >= 0,
            "Error: Invalid value for flux_normal_axis. It must be " + flux_normal_axis_help);
        pp_species_name.get("flux_direction", flux_direction);
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(flux_direction == +1 || flux_direction == -1,
            "Error: flux_direction must be -1 or +1.");
        // Construct InjectorPosition with InjectorPositionRandom.
        h_inj_pos = std::make_unique<InjectorPosition>(
            (InjectorPositionRandomPlane*)nullptr,
            xmin, xmax, ymin, ymax, zmin, zmax,
            flux_normal_axis);
        parseDensity(pp_species_name);
        parseMomentum(pp_species_name);
    } else if (injection_style == "nuniformpercell") {
        // Note that for RZ, three numbers are expected, r, theta, and z.
        // For 2D, only two are expected. The third is overwritten with 1.
        // For 1D, only one is expected. The second and third are overwritten with 1.
#if defined(WARPX_DIM_1D_Z)
        constexpr int num_required_ppc_each_dim = 1;
#elif defined(WARPX_DIM_XZ)
        constexpr int num_required_ppc_each_dim = 2;
#else
        constexpr int num_required_ppc_each_dim = 3;
#endif
        utils::parser::getArrWithParser(
            pp_species_name, "num_particles_per_cell_each_dim",
            num_particles_per_cell_each_dim, 0, num_required_ppc_each_dim);
#if WARPX_DIM_XZ
        num_particles_per_cell_each_dim.push_back(1);
#endif
#if WARPX_DIM_1D_Z
        num_particles_per_cell_each_dim.push_back(1); // overwrite 2nd number with 1
        num_particles_per_cell_each_dim.push_back(1); // overwrite 3rd number with 1
#endif
#if WARPX_DIM_RZ
        if (WarpX::n_rz_azimuthal_modes > 1) {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            num_particles_per_cell_each_dim[1]>=2*WarpX::n_rz_azimuthal_modes,
            "Error: For accurate use of WarpX cylindrical geometry the number "
            "of particles in the theta direction should be at least two times "
            "n_rz_azimuthal_modes (Please visit PR#765 for more information.)");
        }
#endif
        // Construct InjectorPosition from InjectorPositionRegular.
        h_inj_pos = std::make_unique<InjectorPosition>(
            (InjectorPositionRegular*)nullptr,
            xmin, xmax, ymin, ymax, zmin, zmax,
            amrex::Dim3{num_particles_per_cell_each_dim[0],
                num_particles_per_cell_each_dim[1],
                num_particles_per_cell_each_dim[2]});
        num_particles_per_cell = num_particles_per_cell_each_dim[0] *
                                 num_particles_per_cell_each_dim[1] *
                                 num_particles_per_cell_each_dim[2];
        parseDensity(pp_species_name);
        parseMomentum(pp_species_name);
    } else if (injection_style == "external_file") {
#ifndef WARPX_USE_OPENPMD
        amrex::Abort(Utils::TextMsg::Err(
            "WarpX has to be compiled with USE_OPENPMD=TRUE to be able"
            " to read the external openPMD file with species data"));
#endif
        external_file = true;
        std::string str_injection_file;
        pp_species_name.get("injection_file", str_injection_file);
        // optional parameters
        utils::parser::queryWithParser(pp_species_name, "q_tot", q_tot);
        utils::parser::queryWithParser(pp_species_name, "z_shift",z_shift);

#ifdef WARPX_USE_OPENPMD
        if (amrex::ParallelDescriptor::IOProcessor()) {
            m_openpmd_input_series = std::make_unique<openPMD::Series>(
                str_injection_file, openPMD::Access::READ_ONLY);

            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_openpmd_input_series->iterations.size() == 1u,
                "External file should contain only 1 iteration\n");
            openPMD::Iteration it = m_openpmd_input_series->iterations.begin()->second;
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                it.particles.size() == 1u,
                "External file should contain only 1 species\n");
            std::string const ps_name = it.particles.begin()->first;
            openPMD::ParticleSpecies ps = it.particles.begin()->second;

            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                ps.contains("charge") || charge_is_specified || species_is_specified,
                std::string("'") + ps_name +
                ".injection_file' does not contain a 'charge' species record. "
                "Please specify '" + ps_name + ".charge' or "
                "'" + ps_name + ".species_type' in your input file!\n");
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                ps.contains("mass") || mass_is_specified || species_is_specified,
                std::string("'") + ps_name +
                ".injection_file' does not contain a 'mass' species record. "
                "Please specify '" + ps_name + ".mass' or "
                "'" + ps_name + ".species_type' in your input file!\n");

            if (charge_is_specified) {
                ablastr::warn_manager::WMRecordWarning("Species",
                    "Both '" + ps_name + ".charge' and '" +
                        ps_name + ".injection_file' specify a charge.\n'" +
                        ps_name + ".charge' will take precedence.\n");
            }
            else if (species_is_specified) {
                ablastr::warn_manager::WMRecordWarning("Species",
                    "Both '" + ps_name + ".species_type' and '" +
                        ps_name + ".injection_file' specify a charge.\n'" +
                        ps_name + ".species_type' will take precedence.\n");
            }
            else {
                // TODO: Add ASSERT_WITH_MESSAGE to test if charge is a constant record
                amrex::ParticleReal const p_q =
                    ps["charge"][openPMD::RecordComponent::SCALAR].loadChunk<amrex::ParticleReal>().get()[0];
                double const charge_unit = ps["charge"][openPMD::RecordComponent::SCALAR].unitSI();
                charge = p_q * charge_unit;
            }
            if (mass_is_specified) {
                ablastr::warn_manager::WMRecordWarning("Species",
                    "Both '" + ps_name + ".mass' and '" +
                        ps_name + ".injection_file' specify a charge.\n'" +
                        ps_name + ".mass' will take precedence.\n");
            }
            else if (species_is_specified) {
                ablastr::warn_manager::WMRecordWarning("Species",
                    "Both '" + ps_name + ".species_type' and '" +
                        ps_name + ".injection_file' specify a mass.\n'" +
                        ps_name + ".species_type' will take precedence.\n");
            }
            else {
                // TODO: Add ASSERT_WITH_MESSAGE to test if mass is a constant record
                amrex::ParticleReal const p_m =
                    ps["mass"][openPMD::RecordComponent::SCALAR].loadChunk<amrex::ParticleReal>().get()[0];
                double const mass_unit = ps["mass"][openPMD::RecordComponent::SCALAR].unitSI();
                mass = p_m * mass_unit;
            }
        } // IOProcessor

        // Broadcast charge and mass to non-IO processors
        if (!charge_is_specified && !species_is_specified)
            amrex::ParallelDescriptor::Bcast(&charge, 1,
                amrex::ParallelDescriptor::IOProcessorNumber());
        if (!mass_is_specified && !species_is_specified)
            amrex::ParallelDescriptor::Bcast(&mass, 1,
                amrex::ParallelDescriptor::IOProcessorNumber());
#else
        amrex::Abort(Utils::TextMsg::Err(
            "Plasma injection via external_file requires openPMD support: "
            "Add USE_OPENPMD=TRUE when compiling WarpX."));
#endif  // WARPX_USE_OPENPMD

    } else {
        StringParseAbortMessage("Injection style", injection_style);
    }

    if (h_inj_pos) {
#ifdef AMREX_USE_GPU
        d_inj_pos = static_cast<InjectorPosition*>
            (amrex::The_Arena()->alloc(sizeof(InjectorPosition)));
        amrex::Gpu::htod_memcpy_async(d_inj_pos, h_inj_pos.get(), sizeof(InjectorPosition));
#else
        d_inj_pos = h_inj_pos.get();
#endif
    }

    if (h_inj_rho) {
#ifdef AMREX_USE_GPU
        d_inj_rho = static_cast<InjectorDensity*>
            (amrex::The_Arena()->alloc(sizeof(InjectorDensity)));
        amrex::Gpu::htod_memcpy_async(d_inj_rho, h_inj_rho.get(), sizeof(InjectorDensity));
#else
        d_inj_rho = h_inj_rho.get();
#endif
    }

    if (h_inj_mom) {
#ifdef AMREX_USE_GPU
        d_inj_mom = static_cast<InjectorMomentum*>
            (amrex::The_Arena()->alloc(sizeof(InjectorMomentum)));
        amrex::Gpu::htod_memcpy_async(d_inj_mom, h_inj_mom.get(), sizeof(InjectorMomentum));
#else
        d_inj_mom = h_inj_mom.get();
#endif
    }

    amrex::Gpu::synchronize();
}

PlasmaInjector::~PlasmaInjector ()
{
#ifdef AMREX_USE_GPU
    if (d_inj_pos) {
        amrex::The_Arena()->free(d_inj_pos);
    }
    if (d_inj_rho) {
        amrex::The_Arena()->free(d_inj_rho);
    }
    if (d_inj_mom) {
        amrex::The_Arena()->free(d_inj_mom);
    }
#endif
}

// Depending on injection type at runtime, initialize inj_rho
// so that inj_rho->getDensity calls
// InjectorPosition[Constant or Custom or etc.].getDensity.
void PlasmaInjector::parseDensity (amrex::ParmParse& pp)
{
    // parse density information
    std::string rho_prof_s;
    pp.get("profile", rho_prof_s);
    std::transform(rho_prof_s.begin(), rho_prof_s.end(),
                   rho_prof_s.begin(), ::tolower);
    if (rho_prof_s == "constant") {
        utils::parser::getWithParser(pp, "density", density);
        // Construct InjectorDensity with InjectorDensityConstant.
        h_inj_rho.reset(new InjectorDensity((InjectorDensityConstant*)nullptr, density));
    } else if (rho_prof_s == "custom") {
        // Construct InjectorDensity with InjectorDensityCustom.
        h_inj_rho.reset(new InjectorDensity((InjectorDensityCustom*)nullptr, species_name));
    } else if (rho_prof_s == "predefined") {
        // Construct InjectorDensity with InjectorDensityPredefined.
        h_inj_rho.reset(new InjectorDensity((InjectorDensityPredefined*)nullptr,species_name));
    } else if (rho_prof_s == "parse_density_function") {
        utils::parser::Store_parserString(
            pp, "density_function(x,y,z)", str_density_function);
        // Construct InjectorDensity with InjectorDensityParser.
        density_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(str_density_function,{"x","y","z"}));
        h_inj_rho.reset(new InjectorDensity((InjectorDensityParser*)nullptr,
            density_parser->compile<3>()));
    } else {
        //No need for profile definition if external file is used
        std::string injection_style = "none";
        pp.query("injection_style", injection_style);
        if (injection_style != "external_file") {
            StringParseAbortMessage("Density profile type", rho_prof_s);
        }
    }
}

// Depending on injection type at runtime, initialize inj_mom
// so that inj_mom->getMomentum calls
// InjectorMomentum[Constant or Custom or etc.].getMomentum.
void PlasmaInjector::parseMomentum (amrex::ParmParse& pp)
{
    using namespace amrex::literals;

    // parse momentum information
    std::string mom_dist_s;
    pp.get("momentum_distribution_type", mom_dist_s);
    std::transform(mom_dist_s.begin(),
                   mom_dist_s.end(),
                   mom_dist_s.begin(),
                   ::tolower);
    if (mom_dist_s == "at_rest") {
        constexpr amrex::Real ux = 0._rt;
        constexpr amrex::Real uy = 0._rt;
        constexpr amrex::Real uz = 0._rt;
        // Construct InjectorMomentum with InjectorMomentumConstant.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumConstant*)nullptr, ux, uy, uz));
    } else if (mom_dist_s == "constant") {
        amrex::Real ux = 0._rt;
        amrex::Real uy = 0._rt;
        amrex::Real uz = 0._rt;
        utils::parser::queryWithParser(pp, "ux", ux);
        utils::parser::queryWithParser(pp, "uy", uy);
        utils::parser::queryWithParser(pp, "uz", uz);
        // Construct InjectorMomentum with InjectorMomentumConstant.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumConstant*)nullptr, ux, uy, uz));
    } else if (mom_dist_s == "custom") {
        // Construct InjectorMomentum with InjectorMomentumCustom.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumCustom*)nullptr, species_name));
    } else if (mom_dist_s == "gaussian") {
        amrex::Real ux_m = 0._rt;
        amrex::Real uy_m = 0._rt;
        amrex::Real uz_m = 0._rt;
        amrex::Real ux_th = 0._rt;
        amrex::Real uy_th = 0._rt;
        amrex::Real uz_th = 0._rt;
        utils::parser::queryWithParser(pp, "ux_m", ux_m);
        utils::parser::queryWithParser(pp, "uy_m", uy_m);
        utils::parser::queryWithParser(pp, "uz_m", uz_m);
        utils::parser::queryWithParser(pp, "ux_th", ux_th);
        utils::parser::queryWithParser(pp, "uy_th", uy_th);
        utils::parser::queryWithParser(pp, "uz_th", uz_th);
        // Construct InjectorMomentum with InjectorMomentumGaussian.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumGaussian*)nullptr,
                                             ux_m, uy_m, uz_m, ux_th, uy_th, uz_th));
    } else if (mom_dist_s == "gaussianflux") {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(surface_flux,
            "Error: gaussianflux can only be used with injection_style = NFluxPerCell");
        amrex::Real ux_m = 0._rt;
        amrex::Real uy_m = 0._rt;
        amrex::Real uz_m = 0._rt;
        amrex::Real ux_th = 0._rt;
        amrex::Real uy_th = 0._rt;
        amrex::Real uz_th = 0._rt;
        utils::parser::queryWithParser(pp, "ux_m", ux_m);
        utils::parser::queryWithParser(pp, "uy_m", uy_m);
        utils::parser::queryWithParser(pp, "uz_m", uz_m);
        utils::parser::queryWithParser(pp, "ux_th", ux_th);
        utils::parser::queryWithParser(pp, "uy_th", uy_th);
        utils::parser::queryWithParser(pp, "uz_th", uz_th);
        // Construct InjectorMomentum with InjectorMomentumGaussianFlux.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumGaussianFlux*)nullptr,
                                             ux_m, uy_m, uz_m, ux_th, uy_th, uz_th,
                                             flux_normal_axis, flux_direction));
    } else if (mom_dist_s == "maxwell_boltzmann"){
        h_mom_temp = std::make_unique<TemperatureProperties>(pp);
        GetTemperature getTemp(*h_mom_temp.get());
        h_mom_vel = std::make_unique<VelocityProperties>(pp);
        GetVelocity getVel(*h_mom_vel.get());
        // Construct InjectorMomentum with InjectorMomentumBoltzmann.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumBoltzmann*)nullptr, getTemp, getVel));
    } else if (mom_dist_s == "maxwell_juttner"){
        h_mom_temp = std::make_unique<TemperatureProperties>(pp);
        GetTemperature getTemp(*h_mom_temp.get());
        h_mom_vel = std::make_unique<VelocityProperties>(pp);
        GetVelocity getVel(*h_mom_vel.get());
        // Construct InjectorMomentum with InjectorMomentumJuttner.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumJuttner*)nullptr, getTemp, getVel));
    } else if (mom_dist_s == "radial_expansion") {
        amrex::Real u_over_r = 0._rt;
        utils::parser::queryWithParser(pp, "u_over_r", u_over_r);
        // Construct InjectorMomentum with InjectorMomentumRadialExpansion.
        h_inj_mom.reset(new InjectorMomentum
                        ((InjectorMomentumRadialExpansion*)nullptr, u_over_r));
    } else if (mom_dist_s == "parse_momentum_function") {
        utils::parser::Store_parserString(pp, "momentum_function_ux(x,y,z)",
            str_momentum_function_ux);
        utils::parser::Store_parserString(pp, "momentum_function_uy(x,y,z)",
            str_momentum_function_uy);
        utils::parser::Store_parserString(pp, "momentum_function_uz(x,y,z)",
            str_momentum_function_uz);
        // Construct InjectorMomentum with InjectorMomentumParser.
        ux_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(str_momentum_function_ux, {"x","y","z"}));
        uy_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(str_momentum_function_uy, {"x","y","z"}));
        uz_parser = std::make_unique<amrex::Parser>(
            utils::parser::makeParser(str_momentum_function_uz, {"x","y","z"}));
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumParser*)nullptr,
                                             ux_parser->compile<3>(),
                                             uy_parser->compile<3>(),
                                             uz_parser->compile<3>()));
    } else {
        //No need for momentum definition if external file is used
        std::string injection_style = "none";
        pp.query("injection_style", injection_style);
        if (injection_style != "external_file") {
            StringParseAbortMessage("Momentum distribution type", mom_dist_s);
        }
    }
}

amrex::XDim3 PlasmaInjector::getMomentum (amrex::Real x,
                                          amrex::Real y,
                                          amrex::Real z) const noexcept
{
    return h_inj_mom->getMomentum(x, y, z, amrex::RandomEngine{}); // gamma*beta
}

bool PlasmaInjector::insideBounds (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
{
    return (x < xmax and x >= xmin and
            y < ymax and y >= ymin and
            z < zmax and z >= zmin);
}

bool PlasmaInjector::overlapsWith (const amrex::XDim3& lo,
                                   const amrex::XDim3& hi) const noexcept
{
    return ! (   (xmin > hi.x) || (xmax < lo.x)
              || (ymin > hi.y) || (ymax < lo.y)
              || (zmin > hi.z) || (zmax < lo.z) );
}

InjectorPosition*
PlasmaInjector::getInjectorPosition ()
{
    return d_inj_pos;
}

InjectorDensity*
PlasmaInjector::getInjectorDensity ()
{
    return d_inj_rho;
}

InjectorMomentum*
PlasmaInjector::getInjectorMomentum ()
{
    return d_inj_mom;
}
