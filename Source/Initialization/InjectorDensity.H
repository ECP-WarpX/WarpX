/* Copyright 2019 Axel Huebl, Maxence Thevenet, Weiqun Zhang
 * Michael Rowan
 *
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef INJECTOR_DENSITY_H_
#define INJECTOR_DENSITY_H_

#include "Utils/WarpXConst.H"


#include <AMReX.H>
//include <WarpX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_Config.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_INT.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <string>

#ifdef WARPX_USE_OPENPMD
#   include <openPMD/openPMD.hpp>
namespace io = openPMD;
#endif

#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>
#include <memory>
#include <string>
#include <utility>

//using namespace amrex;
//namespace

// struct whose getDensity returns constant density.
struct InjectorDensityConstant
{
    InjectorDensityConstant (amrex::Real a_rho) noexcept : m_rho(a_rho) {}

    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return m_rho;
    }

private:
    amrex::Real m_rho;
};

// struct whose getDensity returns local density computed from parser.
struct InjectorDensityParser
{
    InjectorDensityParser (amrex::ParserExecutor<3> const& a_parser) noexcept
        : m_parser(a_parser) {}

    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return m_parser(x,y,z);
    }

    amrex::ParserExecutor<3> m_parser;
};

// struct whose getDensity returns local density computed from predefined profile.
struct InjectorDensityPredefined
{
    InjectorDensityPredefined (std::string const& a_species_name) noexcept;

    void clear ();

    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        // Choices for profile are:
        // - parabolic_channel
        switch (profile)
        {
        case Profile::parabolic_channel:
        {
            // These are cast as double to ensure sufficient precision in the
            // initialized profile density profile; without these, single and
            // double precision versions of the executable can show disagreement
            // From testing, it seems that (for at least some setups), it is only
            // necessary to case n0 as double to get good agreement between single
            // and double precision, but just in case, all are cast as double.
            const double z_start   = p[0];
            const double ramp_up   = p[1];
            const double plateau   = p[2];
            const double ramp_down = p[3];
            const double rc        = p[4];
            const double n0        = p[5];
            const double kp = PhysConst::q_e/PhysConst::c
                *std::sqrt( n0/(PhysConst::m_e*PhysConst::ep0) );
            double n;

            // Longitudinal profile, normalized to 1
            if        ((z-z_start)>=0               and
                       (z-z_start)<ramp_up ) {
                n = 0.5*(1.-std::cos(MathConst::pi*(z-z_start)/ramp_up));
            } else if ((z-z_start)>=ramp_up         and
                       (z-z_start)< ramp_up+plateau ) {
                n = 1.;
            } else if ((z-z_start)>=ramp_up+plateau and
                       (z-z_start)< ramp_up+plateau+ramp_down) {
                n = 0.5*(1.+std::cos(MathConst::pi*((z-z_start)-ramp_up-plateau)/ramp_down));
            } else {
                n = 0.;
            }
            // Multiply by transverse profile, and physical density
            n *= n0*(1.+4.*(x*x+y*y)/(kp*kp*rc*rc*rc*rc));
            return static_cast<amrex::Real>(n);
        }
        default:
            amrex::Abort("InjectorDensityPredefined: how did we get here?");
            return amrex::Real(0.0);
        }
    }

private:
    enum struct Profile { null, parabolic_channel };
    Profile profile{Profile::null};
    amrex::GpuArray<amrex::Real,6> p;
};

struct InjectorDensityFromFile
{
    InjectorDensityFromFile (std::string const& a_species_name);

    //~InjectorDensityFromFile() = default;

    void clear();


    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {

//
//        auto series = io::Series(file_path, io::Access::READ_ONLY);
//        auto iseries  = series.iterations.begin()->second;
//        auto P = iseries.meshes[dens];
//        auto p_scalar = P[io::RecordComponent::SCALAR];
//        auto particle = iseries.particles[species];
//
//        std::shared_ptr<double> charge = particle["charge"][openPMD::RecordComponent::SCALAR].loadChunk<double>();
//        series.flush();
//
//        auto fileGeom = P.getAttribute("geometry").get< std::string >();
//        const auto offset = P.gridGlobalOffset();
//        const auto spacing = P.gridSpacing< long double >();
//        const auto extent = p_scalar.getExtent();
//
//        const io::Offset chunk_offset = {0,0,0};
//        const io::Extent chunk_extent = {extent[0], extent[1], extent[2]};
//
//        auto P_chunk_data = p_scalar.loadChunk< double > (chunk_offset, chunk_extent);
//        series.flush();
//
//          WarpX& warpx = WarpX::GetInstance();
//          amrex::Geometry const& geom0 = warpx.Geom(0);
//          const amrex::RealBox& real_box = geom0.ProbDomain();
//        const auto dx = geom0.CellSizeArray();
//        //const amrex::IntVect nodal_flag = mf->ixType().toIntVect();

        //convert x,y,z to corresponding i,j,k to call the correct index
//        auto series = io::Series(external_density_path, io::Access::READ_ONLY);
//        auto iseries  = series.iterations.begin()->second;
//        auto P = iseries.meshes[density];
//        auto p_scalar = P[io::RecordComponent::SCALAR];
//        auto particle = iseries.particles[species];

        int i = std::floor((x-real_box.lo(0))/spacing[0]);
        int j = std::floor((y-real_box.lo(1))/spacing[1]);
        int k = std::floor((z-real_box.lo(2))/spacing[2]);

        double val = mffab(i, j, k);
//        double val = 3.3e26;

        return static_cast<amrex::Real>(val);
    }
//private:
    std::string external_density_path;
    std::string density;
    amrex::Array4<double> mffab;
    amrex::RealBox real_box;
    std::vector<long double> spacing;
    std::vector<int> n_cell;
    int max_grid_size;
};





// Base struct for density injector.
// InjectorDensity contains a union (called Object) that holds any one
// instance of:
// - InjectorDensityConstant  : to generate constant density;
// - InjectorDensityParser    : to generate density from parser;
// - InjectorDensityPredefined: to generate density from predefined profile;
// - InjectorDensityFromFile  : to read density from an external file;
// The choice is made at runtime, depending in the constructor called.
// This mimics virtual functions.
struct InjectorDensity
{
    // This constructor stores a InjectorDensityConstant in union object.
    InjectorDensity (InjectorDensityConstant* t, amrex::Real a_rho)
        : type(Type::constant),
          object(t,a_rho)
    { }

    // This constructor stores a InjectorDensityParser in union object.
    InjectorDensity (InjectorDensityParser* t, amrex::ParserExecutor<3> const& a_parser)
        : type(Type::parser),
          object(t,a_parser)
    { }

    // This constructor stores a InjectorDensityPredefined in union object.
    InjectorDensity (InjectorDensityPredefined* t, std::string const& a_species_name)
        : type(Type::predefined),
          object(t,a_species_name)
    { }
    //This constructor stores a InjectorDensityFromFile in union object
    InjectorDensity (InjectorDensityFromFile* t, std::string const& a_species_name)
            : type(Type::fromfile),
              object(t, a_species_name)
    { }

    // Explicitly prevent the compiler from generating copy constructors
    // and copy assignment operators.
    InjectorDensity (InjectorDensity const&) = delete;
    InjectorDensity (InjectorDensity&&) = delete;
    void operator= (InjectorDensity const&) = delete;
    void operator= (InjectorDensity &&) = delete;

    // Default destructor
//    ~InjectorDensity ()
//    {
//        // TODO: clean up fromfile
//    }

    ~InjectorDensity () = default;

    void clear ();

    // call getDensity from the object stored in the union
    // (the union is called Object, and the instance is called object).
    [[nodiscard]]
    AMREX_GPU_HOST_DEVICE
    amrex::Real
    getDensity (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (type)
        {
        case Type::parser:
        {
            return object.parser.getDensity(x,y,z);
        }
        case Type::constant:
        {
            return object.constant.getDensity(x,y,z);
        }
        case Type::predefined:
        {
            return object.predefined.getDensity(x,y,z);
        }
        case Type::fromfile:
        {
            return object.fromfile.getDensity(x,y,z);
        }
        default:
        {
            amrex::Abort("InjectorDensity: unknown type");
            return 0.0;
        }
        }
    }

private:
    enum struct Type { constant, predefined, parser, fromfile };
    Type type;

    // An instance of union Object constructs and stores any one of
    // the objects declared (constant or parser or predefined or fromfile).
    union Object {
        Object (InjectorDensityConstant*, amrex::Real a_rho) noexcept
            : constant(a_rho) {}
        Object (InjectorDensityParser*, amrex::ParserExecutor<3> const& a_parser) noexcept
            : parser(a_parser) {}
        Object (InjectorDensityPredefined*, std::string const& a_species_name) noexcept
            : predefined(a_species_name) {}
        Object (InjectorDensityFromFile*, std::string const& a_species_name) noexcept
            : fromfile(a_species_name) {}
        ~Object() { }
        InjectorDensityConstant   constant;
        InjectorDensityParser     parser;
        InjectorDensityPredefined predefined;
        InjectorDensityFromFile   fromfile;
    };
    Object object;
};

// In order for InjectorDensity to be trivially copyable, its destructor
// must be trivial.  So we have to rely on a custom deleter for unique_ptr.
struct InjectorDensityDeleter {
    void operator () (InjectorDensity* p) const {
        if (p) {
            p->clear();
            delete p;
        }
//        if (external_density_path) {
//            external_density_path->clear();
//            density->clear();
//            mffab->clear();
//            real_box->clear();
//            spacing->clear();
//
//            delete external_density_path;
//            delete density;
//            delete mffab;
//            delete real_box;
//            delete spacing;
//        }
    }
};

#endif
