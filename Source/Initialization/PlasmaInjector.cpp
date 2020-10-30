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
#include "Particles/SpeciesPhysicalProperties.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <functional>
#include <sstream>
#include <string>
#include <memory>


using namespace amrex;

namespace {
    void StringParseAbortMessage(const std::string& var,
                                 const std::string& name) {
        std::stringstream stringstream;
        std::string string;
        stringstream << var << " string '" << name << "' not recognized.";
        string = stringstream.str();
        amrex::Abort(string.c_str());
    }

    Real parseChargeName(const ParmParse& pp, const std::string& name) {
        Real result;
        if (name == "q_e") {
            return PhysConst::q_e;
        } else if (pp.query("charge", result)) {
            return result;
        } else {
            StringParseAbortMessage("Charge", name);
            return 0.0;
        }
    }

    Real parseChargeString(const ParmParse& pp, const std::string& name) {
        if(name.substr(0, 1) == "-")
            return -1.0 * parseChargeName(pp, name.substr(1, name.size() - 1));
        return parseChargeName(pp, name);
    }

    Real parseMassString(const ParmParse& pp, const std::string& name) {
        Real result;
        if (name == "m_e") {
            return PhysConst::m_e;
        } else if (name == "m_p"){
            return PhysConst::m_p;
        } else if (name == "inf"){
            return std::numeric_limits<double>::infinity();
        } else if (pp.query("mass", result)) {
            return result;
        } else {
            StringParseAbortMessage("Mass", name);
            return 0.0;
        }
    }
}

PlasmaInjector::PlasmaInjector () {}

PlasmaInjector::PlasmaInjector (int ispecies, const std::string& name)
    : species_id(ispecies), species_name(name)
{
    ParmParse pp(species_name);

    static_assert(std::is_trivially_copyable<InjectorPosition>::value,
                  "InjectorPosition must be trivially copyable");
    static_assert(std::is_trivially_copyable<InjectorDensity>::value,
                  "InjectorDensity must be trivially copyable");
    static_assert(std::is_trivially_copyable<InjectorMomentum>::value,
                  "InjectorMomentum must be trivially copyable");

    pp.query("radially_weighted", radially_weighted);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(radially_weighted, "ERROR: Only radially_weighted=true is supported");

    // Unlimited boundaries
    xmin = std::numeric_limits<amrex::Real>::lowest();
    ymin = std::numeric_limits<amrex::Real>::lowest();
    zmin = std::numeric_limits<amrex::Real>::lowest();

    xmax = std::numeric_limits<amrex::Real>::max();
    ymax = std::numeric_limits<amrex::Real>::max();
    zmax = std::numeric_limits<amrex::Real>::max();

    // NOTE: When periodic boundaries are used, default injection range is set to mother grid dimensions.
    const Geometry& geom = WarpX::GetInstance().Geom(0);
    if( geom.isPeriodic(0) ) {
        xmin = geom.ProbLo(0);
        xmax = geom.ProbHi(0);
    }

    if( geom.isPeriodic(1) ) {
#       ifndef WARPX_DIM_3D
        zmin = geom.ProbLo(1);
        zmax = geom.ProbHi(1);
#       else
        ymin = geom.ProbLo(1);
        ymax = geom.ProbHi(1);
#       endif
    }

#   ifdef WARPX_DIM_3D
    if( geom.isPeriodic(2) ) {
        zmin = geom.ProbLo(2);
        zmax = geom.ProbHi(2);
    }
#   endif

    pp.query("xmin", xmin);
    pp.query("ymin", ymin);
    pp.query("zmin", zmin);
    pp.query("xmax", xmax);
    pp.query("ymax", ymax);
    pp.query("zmax", zmax);

    pp.query("density_min", density_min);
    pp.query("density_max", density_max);

    std::string physical_species_s;
    bool species_is_specified = pp.query("species_type", physical_species_s);
    if (species_is_specified){
        physical_species = species::from_string( physical_species_s );
        // charge = SpeciesCharge[physical_species];
        charge = species::get_charge( physical_species );
        // mass = SpeciesMass[physical_species];
        mass = species::get_mass( physical_species );
    }

    std::string s_inj_style;
    pp.query("injection_style", s_inj_style);

    // parse charge and mass
    std::string charge_s;
    std::string mass_s;
    bool charge_is_specified = pp.query("charge", charge_s);
    bool mass_is_specified = pp.query("mass", mass_s);

    if (charge_is_specified){
        std::transform(charge_s.begin(),
                       charge_s.end(),
                       charge_s.begin(),
                       ::tolower);
        charge = parseChargeString(pp, charge_s);
    }
    if ( charge_is_specified && species_is_specified ){
        Print() << "WARNING: Both '" << species_name << ".charge' and "
                << species_name << ".species_type' are specified\n'"
                << species_name << ".charge' will take precedence.\n";
    }
    if (!charge_is_specified && !species_is_specified && s_inj_style != "external_file"){
        // external file will throw own assertions below if charge cannot be found
        amrex::Abort("Need to specify at least one of species_type or charge");
    }

    if (mass_is_specified){
        std::transform(mass_s.begin(),
                       mass_s.end(),
                       mass_s.begin(),
                       ::tolower);
        mass = parseMassString(pp, mass_s);
    }
    if ( mass_is_specified && species_is_specified ){
        Print() << "WARNING: Both '" << species_name << ".mass' and "
                << species_name << ".species_type' are specified\n'"
                << species_name << ".mass' will take precedence.\n";
    }
    if (!mass_is_specified && !species_is_specified && s_inj_style != "external_file"){
        // external file will throw own assertions below if mass cannot be found
        amrex::Abort("Need to specify at least one of species_type or mass");
    }

    // parse injection style
    std::string part_pos_s;
    pp.get("injection_style", part_pos_s);
    std::transform(part_pos_s.begin(),
                   part_pos_s.end(),
                   part_pos_s.begin(),
                   ::tolower);
    num_particles_per_cell_each_dim.assign(3, 0);
    if (part_pos_s == "python") {
        return;
    } else if (part_pos_s == "singleparticle") {
        pp.getarr("single_particle_pos", single_particle_pos, 0, 3);
        pp.getarr("single_particle_vel", single_particle_vel, 0, 3);
        for (auto& x : single_particle_vel) {
            x *= PhysConst::c;
        }
        pp.get("single_particle_weight", single_particle_weight);
        add_single_particle = true;
        return;
    } else if (part_pos_s == "gaussian_beam") {
        pp.get("x_m", x_m);
        pp.get("y_m", y_m);
        pp.get("z_m", z_m);
        pp.get("x_rms", x_rms);
        pp.get("y_rms", y_rms);
        pp.get("z_rms", z_rms);
        pp.query("x_cut", x_cut);
        pp.query("y_cut", y_cut);
        pp.query("z_cut", z_cut);
        pp.get("q_tot", q_tot);
        pp.get("npart", npart);
        pp.query("do_symmetrize", do_symmetrize);
        gaussian_beam = true;
        parseMomentum(pp);
    }
    // Depending on injection type at runtime, initialize inj_pos
    // so that inj_pos->getPositionUnitBox calls
    // InjectorPosition[Random or Regular].getPositionUnitBox.
    else if (part_pos_s == "nrandompercell") {
        pp.query("num_particles_per_cell", num_particles_per_cell);
#if WARPX_DIM_RZ
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            num_particles_per_cell>=2*WarpX::n_rz_azimuthal_modes,
            "Error: For accurate use of WarpX cylindrical gemoetry the number "
            "of particles should be at least two times n_rz_azimuthal_modes "
            "(Please visit PR#765 for more information.)");
#endif
        // Construct InjectorPosition with InjectorPositionRandom.
        h_inj_pos = std::make_unique<InjectorPosition>(
            (InjectorPositionRandom*)nullptr,
            xmin, xmax, ymin, ymax, zmin, zmax);
        parseDensity(pp);
        parseMomentum(pp);
    } else if (part_pos_s == "nuniformpercell") {
        // Note that for RZ, three numbers are expected, r, theta, and z.
        // For 2D, only two are expected. The third is overwritten with 1.
        num_particles_per_cell_each_dim.assign(3, 1);
        pp.getarr("num_particles_per_cell_each_dim", num_particles_per_cell_each_dim);
#if WARPX_DIM_XZ
        num_particles_per_cell_each_dim[2] = 1;
#endif
#if WARPX_DIM_RZ
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            num_particles_per_cell_each_dim[1]>=2*WarpX::n_rz_azimuthal_modes,
            "Error: For accurate use of WarpX cylindrical gemoetry the number "
            "of particles in the theta direction should be at least two times "
            "n_rz_azimuthal_modes (Please visit PR#765 for more information.)");
#endif
        // Construct InjectorPosition from InjectorPositionRegular.
        h_inj_pos = std::make_unique<InjectorPosition>(
            (InjectorPositionRegular*)nullptr,
            xmin, xmax, ymin, ymax, zmin, zmax,
            Dim3{num_particles_per_cell_each_dim[0],
                num_particles_per_cell_each_dim[1],
                num_particles_per_cell_each_dim[2]});
        num_particles_per_cell = num_particles_per_cell_each_dim[0] *
                                 num_particles_per_cell_each_dim[1] *
                                 num_particles_per_cell_each_dim[2];
        parseDensity(pp);
        parseMomentum(pp);
    } else if (part_pos_s == "external_file") {
#ifndef WARPX_USE_OPENPMD
        amrex::Abort("WarpX has to be compiled with USE_OPENPMD=TRUE to be able"
                     " to read the external openPMD file with species data");
#endif
        external_file = true;
        std::string str_injection_file;
        pp.get("injection_file", str_injection_file);
        // optional parameters
        pp.query("q_tot", q_tot);
        pp.query("z_shift",z_shift);

#ifdef WARPX_USE_OPENPMD
        if (ParallelDescriptor::IOProcessor()) {
            m_openpmd_input_series = std::make_unique<openPMD::Series>(
                str_injection_file, openPMD::Access::READ_ONLY);

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                m_openpmd_input_series->iterations.size() == 1u,
                "External file should contain only 1 iteration\n");
            openPMD::Iteration it = m_openpmd_input_series->iterations.begin()->second;
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                it.particles.size() == 1u,
                "External file should contain only 1 species\n");
            std::string const ps_name = it.particles.begin()->first;
            openPMD::ParticleSpecies ps = it.particles.begin()->second;

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ps.contains("charge") || charge_is_specified || species_is_specified,
                std::string("'") + ps_name +
                ".injection_file' does not contain a 'charge' species record. "
                "Please specify '" + ps_name + ".charge' or "
                "'" + ps_name + ".species_type' in your input file!\n");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ps.contains("mass") || mass_is_specified || species_is_specified,
                std::string("'") + ps_name +
                ".injection_file' does not contain a 'mass' species record. "
                "Please specify '" + ps_name + ".mass' or "
                "'" + ps_name + ".species_type' in your input file!\n");

            if (charge_is_specified) {
                Print() << "WARNING: Both '" << ps_name << ".charge' and '"
                        << ps_name << ".injection_file' specify a charge.\n'"
                        << ps_name << ".charge' will take precedence.\n";
            }
            else if (species_is_specified) {
                Print() << "WARNING: Both '" << ps_name << ".species_type' and '"
                        << ps_name << ".injection_file' specify a charge.\n'"
                        << ps_name << ".species_type' will take precedence.\n";
            }
            else {
                // TODO: Add ASSERT_WITH_MESSAGE to test if charge is a constant record
                ParticleReal const p_q = ps["charge"][openPMD::RecordComponent::SCALAR].loadChunk<ParticleReal>().get()[0];
                double const charge_unit = ps["charge"][openPMD::RecordComponent::SCALAR].unitSI();
                charge = p_q * charge_unit;
            }
            if (mass_is_specified) {
                Print() << "WARNING: Both '" << ps_name << ".mass' and '"
                        << ps_name << ".injection_file' specify a mass.\n'"
                        << ps_name << ".mass' will take precedence.\n";
            }
            else if (species_is_specified) {
                Print() << "WARNING: Both '" << ps_name << ".species_type' and '"
                        << ps_name << ".injection_file' specify a mass.\n'"
                        << ps_name << ".species_type' will take precedence.\n";
            }
            else {
                // TODO: Add ASSERT_WITH_MESSAGE to test if mass is a constant record
                ParticleReal const p_m = ps["mass"][openPMD::RecordComponent::SCALAR].loadChunk<ParticleReal>().get()[0];
                double const mass_unit = ps["mass"][openPMD::RecordComponent::SCALAR].unitSI();
                mass = p_m * mass_unit;
            }
        } // IOProcessor

        // Broadcast charge and mass to non-IO processors
        if (!charge_is_specified && !species_is_specified)
            ParallelDescriptor::Bcast(&charge, 1,
                ParallelDescriptor::IOProcessorNumber());
        if (!mass_is_specified && !species_is_specified)
            ParallelDescriptor::Bcast(&mass, 1,
                ParallelDescriptor::IOProcessorNumber());
#else
        Abort("Plasma injection via external_file requires openPMD support: "
                     "Add USE_OPENPMD=TRUE when compiling WarpX.\n");
#endif  // WARPX_USE_OPENPMD

    } else {
        StringParseAbortMessage("Injection style", part_pos_s);
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
void PlasmaInjector::parseDensity (ParmParse& pp)
{
    // parse density information
    std::string rho_prof_s;
    pp.get("profile", rho_prof_s);
    std::transform(rho_prof_s.begin(), rho_prof_s.end(),
                   rho_prof_s.begin(), ::tolower);
    if (rho_prof_s == "constant") {
        pp.get("density", density);
        // Construct InjectorDensity with InjectorDensityConstant.
        h_inj_rho.reset(new InjectorDensity((InjectorDensityConstant*)nullptr, density));
    } else if (rho_prof_s == "custom") {
        // Construct InjectorDensity with InjectorDensityCustom.
        h_inj_rho.reset(new InjectorDensity((InjectorDensityCustom*)nullptr, species_name));
    } else if (rho_prof_s == "predefined") {
        // Construct InjectorDensity with InjectorDensityPredefined.
        h_inj_rho.reset(new InjectorDensity((InjectorDensityPredefined*)nullptr,species_name));
    } else if (rho_prof_s == "parse_density_function") {
        Store_parserString(pp, "density_function(x,y,z)", str_density_function);
        // Construct InjectorDensity with InjectorDensityParser.
        h_inj_rho.reset(new InjectorDensity((InjectorDensityParser*)nullptr,
                                            makeParser(str_density_function,{"x","y","z"})));
    } else {
        //No need for profile definition if external file is used
        std::string s_inj_style;
        pp.query("injection_style", s_inj_style);
        if (s_inj_style != "external_file") {
            StringParseAbortMessage("Density profile type", rho_prof_s);
        }
    }
}

// Depending on injection type at runtime, initialize inj_mom
// so that inj_mom->getMomentum calls
// InjectorMomentum[Constant or Custom or etc.].getMomentum.
void PlasmaInjector::parseMomentum (ParmParse& pp)
{
    // parse momentum information
    std::string mom_dist_s;
    pp.get("momentum_distribution_type", mom_dist_s);
    std::transform(mom_dist_s.begin(),
                   mom_dist_s.end(),
                   mom_dist_s.begin(),
                   ::tolower);
    if (mom_dist_s == "constant") {
        Real ux = 0.;
        Real uy = 0.;
        Real uz = 0.;
        pp.query("ux", ux);
        pp.query("uy", uy);
        pp.query("uz", uz);
        // Construct InjectorMomentum with InjectorMomentumConstant.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumConstant*)nullptr, ux,uy, uz));
    } else if (mom_dist_s == "custom") {
        // Construct InjectorMomentum with InjectorMomentumCustom.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumCustom*)nullptr, species_name));
    } else if (mom_dist_s == "gaussian") {
        Real ux_m = 0.;
        Real uy_m = 0.;
        Real uz_m = 0.;
        Real ux_th = 0.;
        Real uy_th = 0.;
        Real uz_th = 0.;
        pp.query("ux_m", ux_m);
        pp.query("uy_m", uy_m);
        pp.query("uz_m", uz_m);
        pp.query("ux_th", ux_th);
        pp.query("uy_th", uy_th);
        pp.query("uz_th", uz_th);
        // Construct InjectorMomentum with InjectorMomentumGaussian.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumGaussian*)nullptr,
                                             ux_m, uy_m, uz_m, ux_th, uy_th, uz_th));
    } else if (mom_dist_s == "maxwell_boltzmann"){
        Real beta = 0.;
        Real theta = 10.;
        int dir = 0;
        std::string direction = "x";
        pp.query("beta", beta);
        if(beta < 0){
            amrex::Abort("Please enter a positive beta value. Drift direction is set with <s_name>.bulk_vel_dir = 'x' or '+x', '-x', 'y' or '+y', etc.");
        }
        pp.query("theta", theta);
        pp.query("bulk_vel_dir", direction);
        if(direction[0] == '-'){
            beta = -beta;
        }
        if((direction == "x" || direction[1] == 'x') ||
           (direction == "X" || direction[1] == 'X')){
            dir = 0;
        } else if ((direction == "y" || direction[1] == 'y') ||
                   (direction == "Y" || direction[1] == 'Y')){
            dir = 1;
        } else if ((direction == "z" || direction[1] == 'z') ||
                   (direction == "Z" || direction[1] == 'Z')){
            dir = 2;
        } else{
            std::stringstream stringstream;
            stringstream << "Cannot interpret <s_name>.bulk_vel_dir input '" << direction << "'. Please enter +/- x, y, or z with no whitespace between the sign and other character.";
            direction = stringstream.str();
            amrex::Abort(direction.c_str());
        }
        // Construct InjectorMomentum with InjectorMomentumBoltzmann.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumBoltzmann*)nullptr, theta, beta, dir));
    } else if (mom_dist_s == "maxwell_juttner"){
        Real beta = 0.;
        Real theta = 10.;
        int dir = 0;
        std::string direction = "x";
        pp.query("beta", beta);
        if(beta < 0){
            amrex::Abort("Please enter a positive beta value. Drift direction is set with <s_name>.bulk_vel_dir = 'x' or '+x', '-x', 'y' or '+y', etc.");
        }
        pp.query("theta", theta);
        pp.query("bulk_vel_dir", direction);
        if(direction[0] == '-'){
            beta = -beta;
        }
        if((direction == "x" || direction[1] == 'x') ||
           (direction == "X" || direction[1] == 'X')){
            dir = 0;
        } else if ((direction == "y" || direction[1] == 'y') ||
                   (direction == "Y" || direction[1] == 'Y')){
            dir = 1;
        } else if ((direction == "z" || direction[1] == 'z') ||
                   (direction == "Z" || direction[1] == 'Z')){
            dir = 2;
        } else{
            std::stringstream stringstream;
            stringstream << "Cannot interpret <s_name>.bulk_vel_dir input '" << direction << "'. Please enter +/- x, y, or z with no whitespace between the sign and other character.";
            direction = stringstream.str();
            amrex::Abort(direction.c_str());
        }
        // Construct InjectorMomentum with InjectorMomentumJuttner.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumJuttner*)nullptr, theta, beta, dir));
    } else if (mom_dist_s == "radial_expansion") {
        Real u_over_r = 0.;
        pp.query("u_over_r", u_over_r);
        // Construct InjectorMomentum with InjectorMomentumRadialExpansion.
        h_inj_mom.reset(new InjectorMomentum
                        ((InjectorMomentumRadialExpansion*)nullptr, u_over_r));
    } else if (mom_dist_s == "parse_momentum_function") {
        Store_parserString(pp, "momentum_function_ux(x,y,z)",
                                               str_momentum_function_ux);
        Store_parserString(pp, "momentum_function_uy(x,y,z)",
                                               str_momentum_function_uy);
        Store_parserString(pp, "momentum_function_uz(x,y,z)",
                                               str_momentum_function_uz);
        // Construct InjectorMomentum with InjectorMomentumParser.
        h_inj_mom.reset(new InjectorMomentum((InjectorMomentumParser*)nullptr,
                                             makeParser(str_momentum_function_ux,{"x","y","z"}),
                                             makeParser(str_momentum_function_uy,{"x","y","z"}),
                                             makeParser(str_momentum_function_uz,{"x","y","z"})));
    } else {
        //No need for momentum definition if external file is used
        std::string s_inj_style;
        pp.query("injection_style", s_inj_style);
        if (s_inj_style != "external_file") {
            StringParseAbortMessage("Momentum distribution type", mom_dist_s);
        }
    }
}

XDim3 PlasmaInjector::getMomentum (Real x, Real y, Real z) const noexcept
{
    return h_inj_mom->getMomentum(x, y, z, amrex::RandomEngine{}); // gamma*beta
}

bool PlasmaInjector::insideBounds (Real x, Real y, Real z) const noexcept
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
