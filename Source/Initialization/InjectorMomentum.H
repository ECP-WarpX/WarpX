/* Copyright 2019-2020 Andrew Myers, Axel Huebl, Cameron Yang,
 * Maxence Thevenet, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef INJECTOR_MOMENTUM_H_
#define INJECTOR_MOMENTUM_H_

#include "CustomMomentumProb.H"
#include "GetTemperature.H"
#include "GetVelocity.H"
#include "TemperatureProperties.H"
#include "VelocityProperties.H"
#include "Utils/WarpXConst.H"

#include <AMReX.H>
#include <AMReX_Config.H>
#include <AMReX_Dim3.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Parser.H>
#include <AMReX_Random.H>

#include <AMReX_BaseFwd.H>

#include <cmath>
#include <string>

// struct whose getMomentum returns constant momentum.
struct InjectorMomentumConstant
{
    InjectorMomentumConstant (amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz) noexcept
        : m_ux(a_ux), m_uy(a_uy), m_uz(a_uz) {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real, amrex::Real, amrex::Real,
                 amrex::RandomEngine const&) const noexcept
    {
        return amrex::XDim3{m_ux,m_uy,m_uz};
    }

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real, amrex::Real, amrex::Real) const noexcept
    {
        return amrex::XDim3{m_ux,m_uy,m_uz};
    }

private:
    amrex::Real m_ux, m_uy, m_uz;
};

// struct whose getMomentum returns momentum for 1 particle, from random
// gaussian distribution.
struct InjectorMomentumGaussian
{
    InjectorMomentumGaussian (amrex::Real a_ux_m, amrex::Real a_uy_m,
                              amrex::Real a_uz_m, amrex::Real a_ux_th,
                              amrex::Real a_uy_th, amrex::Real a_uz_th) noexcept
        : m_ux_m(a_ux_m), m_uy_m(a_uy_m), m_uz_m(a_uz_m),
          m_ux_th(a_ux_th), m_uy_th(a_uy_th), m_uz_th(a_uz_th)
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real /*x*/, amrex::Real /*y*/, amrex::Real /*z*/,
                 amrex::RandomEngine const& engine) const noexcept
    {
        return amrex::XDim3{amrex::RandomNormal(m_ux_m, m_ux_th, engine),
                            amrex::RandomNormal(m_uy_m, m_uy_th, engine),
                            amrex::RandomNormal(m_uz_m, m_uz_th, engine)};
    }

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real /*x*/, amrex::Real /*y*/, amrex::Real /*z*/) const noexcept
    {
        return amrex::XDim3{m_ux_m, m_uy_m, m_uz_m};
    }

private:
    amrex::Real m_ux_m, m_uy_m, m_uz_m;
    amrex::Real m_ux_th, m_uy_th, m_uz_th;
};


// struct whose getMomentum returns momentum for 1 particle, from random
// gaussian flux distribution in the specified direction.
// Along the normal axis, the distribution is v*Gaussian,
// with the sign set by flux_direction.
struct InjectorMomentumGaussianFlux
{
    InjectorMomentumGaussianFlux (amrex::Real a_ux_m, amrex::Real a_uy_m,
                                  amrex::Real a_uz_m, amrex::Real a_ux_th,
                                  amrex::Real a_uy_th, amrex::Real a_uz_th,
                                  int a_flux_normal_axis, int a_flux_direction) noexcept
        : m_ux_m(a_ux_m), m_uy_m(a_uy_m), m_uz_m(a_uz_m),
          m_ux_th(a_ux_th), m_uy_th(a_uy_th), m_uz_th(a_uz_th),
          m_flux_normal_axis(a_flux_normal_axis),
          m_flux_direction(a_flux_direction)
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real /*x*/, amrex::Real /*y*/, amrex::Real /*z*/,
                 amrex::RandomEngine const& engine) const noexcept
    {
        using namespace amrex::literals;
        // Generate the v*Gaussian distribution.
        amrex::Real const urand = 1._rt - amrex::Random(engine);
        amrex::Real ur = std::sqrt(2._rt*std::log(1._rt/urand));
        if (m_flux_direction < 0) ur = -ur;

        // Note: Here, in RZ geometry, the variables `ux` and `uy` actually
        // correspond to the radial and azimuthal component of the momentum
        // (and e.g.`m_flux_normal_axis==1` corresponds to v*Gaussian along theta)
        amrex::Real const ux = (m_flux_normal_axis == 0 ? m_ux_th*ur : amrex::RandomNormal(m_ux_m, m_ux_th, engine));
        amrex::Real const uy = (m_flux_normal_axis == 1 ? m_uy_th*ur : amrex::RandomNormal(m_uy_m, m_uy_th, engine));
        amrex::Real const uz = (m_flux_normal_axis == 2 ? m_uz_th*ur : amrex::RandomNormal(m_uz_m, m_uz_th, engine));
        return amrex::XDim3{ux, uy, uz};
    }

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real /*x*/, amrex::Real /*y*/, amrex::Real /*z*/) const noexcept
    {
        return amrex::XDim3{m_ux_m, m_uy_m, m_uz_m};
    }

private:
    amrex::Real m_ux_m, m_uy_m, m_uz_m;
    amrex::Real m_ux_th, m_uy_th, m_uz_th;
    int m_flux_normal_axis;
    int m_flux_direction;
};

// struct whose getMomentum returns momentum for 1 particle with relativistic
// drift velocity beta, from the Maxwell-Boltzmann distribution.
struct InjectorMomentumBoltzmann
{
    // Constructor whose inputs are:
    // a reference to the initial temperature container t,
    // a reference to the initial velocity container b
    InjectorMomentumBoltzmann(GetTemperature const& t, GetVelocity const& b) noexcept
        : velocity(b), temperature(t)
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real const x, amrex::Real const y, amrex::Real const z,
                 amrex::RandomEngine const& engine) const noexcept
    {
        // Calculate the local temperature and check if it's too
        // high for Boltzmann or less than zero
        amrex::Real const theta = temperature(x,y,z);
        if (theta < 0) {
            amrex::Abort("Negative temperature parameter theta encountered, which is not allowed");
        }
        // Calculate local velocity and abort if |beta|>=1
        amrex::Real const beta = velocity(x,y,z);
        if (beta <= -1 || beta >= 1) {
            amrex::Abort("beta = v/c magnitude greater than or equal to 1");
        }
        // Calculate the value of vave from the local temperature
        amrex::Real const vave = std::sqrt(amrex::Real(2) * theta);
        int const dir = velocity.direction();
        amrex::Real x1, x2, gamma;
        amrex::Real u[3];
        x1 = amrex::Random(engine);
        x2 = amrex::Random(engine);
        // Each value of sqrt(-log(x1))*sin(2*pi*x2) is a sample from a Gaussian
        // distribution with sigma = average velocity / c
        // using the Box-Mueller Method.
        u[(dir+1)%3] = vave*std::sqrt(-std::log(x1)) *std::sin(2*MathConst::pi*x2);
        u[(dir+2)%3] = vave*std::sqrt(-std::log(x1)) *std::cos(2*MathConst::pi*x2);
        u[dir] = vave*std::sqrt(-std::log(amrex::Random(engine)))*
            std::sin(2*MathConst::pi*amrex::Random(engine));
        gamma = u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
        gamma = std::sqrt(1+gamma);
        // The following condition is equtaion 32 in Zenitani 2015
        // (Phys. Plasmas 22, 042116) , called the flipping method. It
        // transforms the intergral: d3x' -> d3x  where d3x' is the volume
        // element for positions in the boosted frame. The particle positions
        // and densities can be initialized in the simulation frame.
        // The flipping method can transform any symmetric distribution from one
        // reference frame to another moving at a relative velocity of beta.
        // An equivalent alternative to this method native to WarpX would be to
        // initialize the particle positions and densities in the frame moving
        // at speed beta, and then perform a Lorentz transform on the positions
        // and MB sampled velocities to the simulation frame.
        x1 = amrex::Random(engine);
        if(-beta*u[dir]/gamma > x1)
        {
          u[dir] = -u[dir];
        }
        // This Lorentz transform is equation 17 in Zenitani.
        // It transforms the integral d3u' -> d3u
        // where d3u' is the volume element for momentum in the boosted frame.
        u[dir] = 1/std::sqrt(1-beta*beta)*(u[dir]+gamma*beta);
        // Note that if beta = 0 then the flipping method and Lorentz transform
        // have no effect on the u[dir] direction.
        return amrex::XDim3 {u[0],u[1],u[2]};
    }

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real const x, amrex::Real const y, amrex::Real const z) const noexcept
    {
        using namespace amrex;
        Real u[3];
        for (int idim = 0; idim < 3; ++idim) u[idim] = 0.0_rt;
        const Real beta = velocity(x,y,z);
        int const dir = velocity.direction();
        const Real gamma = static_cast<amrex::Real>(1./sqrt(1+beta*beta));
        u[dir] = gamma*beta;
        return XDim3 {u[0],u[1],u[2]};
    }

private:
    GetVelocity velocity;
    GetTemperature temperature;
};

// struct whose getMomentum returns momentum for 1 particle with relativistc
// drift velocity beta, from the Maxwell-Juttner distribution. Method is from
// Zenitani 2015 (Phys. Plasmas 22, 042116).
struct InjectorMomentumJuttner
{
    // Constructor whose inputs are:
    // a reference to the initial temperature container t,
    // a reference to the initial velocity container b
    InjectorMomentumJuttner(GetTemperature const& t, GetVelocity const& b) noexcept
        : velocity(b), temperature(t)
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real const x, amrex::Real const y, amrex::Real const z,
                 amrex::RandomEngine const& engine) const noexcept
    {
        // Sobol method for sampling MJ Speeds,
        // from Zenitani 2015 (Phys. Plasmas 22, 042116).
        amrex::Real x1, x2, gamma;
        amrex::Real u [3];
        amrex::Real const theta = temperature(x,y,z);
        // Check if temperature is too low to do sampling method. Abort for now,
        // in future should implement an alternate method e.g. inverse transform
        if (theta < 0.1) {
            amrex::Abort("Temeprature parameter theta is less than minimum 0.1 allowed for Maxwell-Juttner");
        }
        // Calculate local velocity and abort if |beta|>=1
        amrex::Real const beta = velocity(x,y,z);
        if (beta <= -1 || beta >= 1) {
            amrex::Abort("beta = v/c magnitude greater than or equal to 1");
        }
        int const dir = velocity.direction();
        x1 = static_cast<amrex::Real>(0.);
        gamma = static_cast<amrex::Real>(0.);
        u[dir] = static_cast<amrex::Real>(0.);
        // This condition is equation 10 in Zenitani,
        // though x1 is defined differently.
        while(u[dir]-gamma <= x1)
        {
            u[dir] = -theta*
                std::log(amrex::Random(engine)*amrex::Random(engine)*amrex::Random(engine));
            gamma = std::sqrt(1+u[dir]*u[dir]);
            x1 = theta*std::log(amrex::Random(engine));
        }
        // The following code samples a random unit vector
        // and multiplies the result by speed u[dir].
        x1 = amrex::Random(engine);
        x2 = amrex::Random(engine);
        // Direction dir is an input parameter that sets the boost direction:
        // 'x' -> d = 0, 'y' -> d = 1, 'z' -> d = 2.
        u[(dir+1)%3] = 2*u[dir]*std::sqrt(x1*(1-x1))*std::sin(2*MathConst::pi*x2);
        u[(dir+2)%3] = 2*u[dir]*std::sqrt(x1*(1-x1))*std::cos(2*MathConst::pi*x2);
        // The value of dir is the boost direction to be transformed.
        u[dir] = u[dir]*(2*x1-1);
        x1 = amrex::Random(engine);
        // The following condition is equtaion 32 in Zenitani, called
        // The flipping method. It transforms the intergral: d3x' -> d3x
        // where d3x' is the volume element for positions in the boosted frame.
        // The particle positions and densities can be initialized in the
        // simulation frame with this method.
        // The flipping method can similarly transform any
        // symmetric distribution from one reference frame to another moving at
        // a relative velocity of beta.
        // An equivalent alternative to this method native to WarpX
        // would be to initialize the particle positions and densities in the
        // frame moving at speed beta, and then perform a Lorentz transform
        // on their positions and MJ sampled velocities to the simulation frame.
        if(-beta*u[dir]/gamma>x1)
        {
            u[dir] = -u[dir];
        }
        // This Lorentz transform is equation 17 in Zenitani.
        // It transforms the integral d3u' -> d3u
        // where d3u' is the volume element for momentum in the boosted frame.
        u[dir] = 1/std::sqrt(1-beta*beta)*(u[dir]+gamma*beta);
        // Note that if beta = 0 then the flipping method and Lorentz transform
        // have no effect on the u[dir] direction.
        return amrex::XDim3 {u[0],u[1],u[2]};
    }

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real const x, amrex::Real const y, amrex::Real const z) const noexcept
    {
        using namespace amrex;
        Real u[3];
        for (int idim = 0; idim < 3; ++idim) u[idim] = 0.0_rt;
        Real const beta = velocity(x,y,z);
        int const dir = velocity.direction();
        const Real gamma = static_cast<Real>(1./sqrt(1.+beta*beta));
        u[dir] = gamma*beta;
        return XDim3 {u[0],u[1],u[2]};
    }

private:
    GetVelocity velocity;
    GetTemperature temperature;
};

/**
 * \brief struct whose getMomentum returns momentum for 1 particle, for
 * radial expansion.
 *
 * Note - u_over_r is expected to be the normalized momentum gamma*beta
 * divided by the physical position in SI units.
**/
struct InjectorMomentumRadialExpansion
{
    InjectorMomentumRadialExpansion (amrex::Real a_u_over_r) noexcept
        : u_over_r(a_u_over_r)
        {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z,
                 amrex::RandomEngine const&) const noexcept
    {
        return {x*u_over_r, y*u_over_r, z*u_over_r};
    }

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return {x*u_over_r, y*u_over_r, z*u_over_r};
    }

private:
    amrex::Real u_over_r;
};

// struct whose getMomentumm returns local momentum computed from parser.
struct InjectorMomentumParser
{
    InjectorMomentumParser (amrex::ParserExecutor<3> const& a_ux_parser,
                            amrex::ParserExecutor<3> const& a_uy_parser,
                            amrex::ParserExecutor<3> const& a_uz_parser) noexcept
        : m_ux_parser(a_ux_parser), m_uy_parser(a_uy_parser),
          m_uz_parser(a_uz_parser) {}

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z,
                 amrex::RandomEngine const&) const noexcept
    {
        return amrex::XDim3{m_ux_parser(x,y,z),m_uy_parser(x,y,z),m_uz_parser(x,y,z)};
    }

    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        return amrex::XDim3{m_ux_parser(x,y,z),m_uy_parser(x,y,z),m_uz_parser(x,y,z)};
    }

    amrex::ParserExecutor<3> m_ux_parser, m_uy_parser, m_uz_parser;
};

// Base struct for momentum injector.
// InjectorMomentum contains a union (called Object) that holds any one
// instance of:
// - InjectorMomentumConstant       : to generate constant density;
// - InjectorMomentumGaussian       : to generate gaussian distribution;
// - InjectorMomentumGaussianFlux   : to generate v*gaussian distribution;
// - InjectorMomentumRadialExpansion: to generate radial expansion;
// - InjectorMomentumParser         : to generate momentum from parser;
// The choice is made at runtime, depending in the constructor called.
// This mimics virtual functions.
struct InjectorMomentum
{
    // This constructor stores a InjectorMomentumConstant in union object.
    InjectorMomentum (InjectorMomentumConstant* t,
                      amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz)
        : type(Type::constant),
          object(t, a_ux, a_uy, a_uz)
    { }

    // This constructor stores a InjectorMomentumParser in union object.
    InjectorMomentum (InjectorMomentumParser* t,
                      amrex::ParserExecutor<3> const& a_ux_parser,
                      amrex::ParserExecutor<3> const& a_uy_parser,
                      amrex::ParserExecutor<3> const& a_uz_parser)
        : type(Type::parser),
          object(t, a_ux_parser, a_uy_parser, a_uz_parser)
    { }

    // This constructor stores a InjectorMomentumGaussian in union object.
    InjectorMomentum (InjectorMomentumGaussian* t,
                      amrex::Real a_ux_m, amrex::Real a_uy_m, amrex::Real a_uz_m,
                      amrex::Real a_ux_th, amrex::Real a_uy_th, amrex::Real a_uz_th)
        : type(Type::gaussian),
          object(t,a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th)
    { }

    // This constructor stores a InjectorMomentumGaussianFlux in union object.
    InjectorMomentum (InjectorMomentumGaussianFlux* t,
                      amrex::Real a_ux_m, amrex::Real a_uy_m, amrex::Real a_uz_m,
                      amrex::Real a_ux_th, amrex::Real a_uy_th, amrex::Real a_uz_th,
                      int a_flux_normal_axis, int a_flux_direction)
        : type(Type::gaussianflux),
          object(t,a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th,a_flux_normal_axis,a_flux_direction)
    { }

    InjectorMomentum (InjectorMomentumBoltzmann* t,
                       GetTemperature const& temperature, GetVelocity const& velocity)
         : type(Type::boltzmann),
           object(t, temperature, velocity)
    { }

     // This constructor stores a InjectorMomentumJuttner in union object.
     InjectorMomentum (InjectorMomentumJuttner* t,
                       GetTemperature const& temperature, GetVelocity const& velocity)
         : type(Type::juttner),
           object(t, temperature, velocity)
    { }

    // This constructor stores a InjectorMomentumCustom in union object.
    InjectorMomentum (InjectorMomentumCustom* t,
                      std::string const& a_species_name)
        : type(Type::custom),
          object(t, a_species_name)
    { }

    // This constructor stores a InjectorMomentumRadialExpansion in union object.
    InjectorMomentum (InjectorMomentumRadialExpansion* t,
                      amrex::Real u_over_r)
        : type(Type::radial_expansion),
          object(t, u_over_r)
    { }

    // Explicitly prevent the compiler from generating copy constructors
    // and copy assignment operators.
    InjectorMomentum (InjectorMomentum const&) = delete;
    InjectorMomentum (InjectorMomentum&&) = delete;
    void operator= (InjectorMomentum const&) = delete;
    void operator= (InjectorMomentum &&) = delete;

    void clear ();

    // call getMomentum from the object stored in the union
    // (the union is called Object, and the instance is called object).
    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getMomentum (amrex::Real x, amrex::Real y, amrex::Real z,
                 amrex::RandomEngine const& engine) const noexcept
    {
        switch (type)
        {
        case Type::parser:
        {
            return object.parser.getMomentum(x,y,z,engine);
        }
        case Type::gaussian:
        {
            return object.gaussian.getMomentum(x,y,z,engine);
        }
        case Type::gaussianflux:
        {
            return object.gaussianflux.getMomentum(x,y,z,engine);
        }
        case Type::boltzmann:
        {
            return object.boltzmann.getMomentum(x,y,z,engine);
        }
        case Type::juttner:
        {
            return object.juttner.getMomentum(x,y,z,engine);
        }
        case Type::constant:
        {
            return object.constant.getMomentum(x,y,z,engine);
        }
        case Type::radial_expansion:
        {
            return object.radial_expansion.getMomentum(x,y,z,engine);
        }
        case Type::custom:
        {
            return object.custom.getMomentum(x,y,z,engine);
        }
        default:
        {
            amrex::Abort("InjectorMomentum: unknown type");
            return {0.0,0.0,0.0};
        }
        }
    }

    // call getBulkMomentum from the object stored in the union
    // (the union is called Object, and the instance is called object).
    AMREX_GPU_HOST_DEVICE
    amrex::XDim3
    getBulkMomentum (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
        switch (type)
        {
        case Type::parser:
        {
            return object.parser.getBulkMomentum(x,y,z);
        }
        case Type::gaussian:
        {
            return object.gaussian.getBulkMomentum(x,y,z);
        }
        case Type::gaussianflux:
        {
            return object.gaussianflux.getBulkMomentum(x,y,z);
        }
        case Type::boltzmann:
        {
            return object.boltzmann.getBulkMomentum(x,y,z);
        }
        case Type::juttner:
        {
            return object.juttner.getBulkMomentum(x,y,z);
        }
        case Type::constant:
        {
            return object.constant.getBulkMomentum(x,y,z);
        }
        case Type::radial_expansion:
        {
            return object.radial_expansion.getBulkMomentum(x,y,z);
        }
        case Type::custom:
        {
            return object.custom.getBulkMomentum(x,y,z);
        }
        default:
        {
            amrex::Abort("InjectorMomentum: unknown type");
            return {0.0,0.0,0.0};
        }
        }
    }

    enum struct Type { constant, custom, gaussian, gaussianflux, boltzmann, juttner, radial_expansion, parser};
    Type type;

private:

    // An instance of union Object constructs and stores any one of
    // the objects declared (constant or custom or gaussian or
    // radial_expansion or parser).
    union Object {
        Object (InjectorMomentumConstant*,
                amrex::Real a_ux, amrex::Real a_uy, amrex::Real a_uz) noexcept
            : constant(a_ux,a_uy,a_uz) {}
        Object (InjectorMomentumCustom*,
                std::string const& a_species_name) noexcept
            : custom(a_species_name) {}
        Object (InjectorMomentumGaussian*,
                amrex::Real a_ux_m, amrex::Real a_uy_m,
                amrex::Real a_uz_m, amrex::Real a_ux_th,
                amrex::Real a_uy_th, amrex::Real a_uz_th) noexcept
            : gaussian(a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th) {}
        Object (InjectorMomentumGaussianFlux*,
                amrex::Real a_ux_m, amrex::Real a_uy_m,
                amrex::Real a_uz_m, amrex::Real a_ux_th,
                amrex::Real a_uy_th, amrex::Real a_uz_th,
                int a_flux_normal_axis, int a_flux_direction) noexcept
            : gaussianflux(a_ux_m,a_uy_m,a_uz_m,a_ux_th,a_uy_th,a_uz_th,a_flux_normal_axis,a_flux_direction) {}
        Object (InjectorMomentumBoltzmann*,
                GetTemperature const& t, GetVelocity const& b) noexcept
            : boltzmann(t,b) {}
        Object (InjectorMomentumJuttner*,
                GetTemperature const& t, GetVelocity const& b) noexcept
            : juttner(t,b) {}
        Object (InjectorMomentumRadialExpansion*,
                amrex::Real u_over_r) noexcept
            : radial_expansion(u_over_r) {}
        Object (InjectorMomentumParser*,
                amrex::ParserExecutor<3> const& a_ux_parser,
                amrex::ParserExecutor<3> const& a_uy_parser,
                amrex::ParserExecutor<3> const& a_uz_parser) noexcept
            : parser(a_ux_parser, a_uy_parser, a_uz_parser) {}
        InjectorMomentumConstant constant;
        InjectorMomentumCustom   custom;
        InjectorMomentumGaussian gaussian;
        InjectorMomentumGaussianFlux gaussianflux;
        InjectorMomentumBoltzmann boltzmann;
        InjectorMomentumJuttner juttner;
        InjectorMomentumRadialExpansion radial_expansion;
        InjectorMomentumParser   parser;
    };
    Object object;
};

// In order for InjectorMomentum to be trivially copyable, its destructor
// must be trivial.  So we have to rely on a custom deleter for unique_ptr.
struct InjectorMomentumDeleter {
    void operator () (InjectorMomentum* p) const {
        if (p) {
            p->clear();
            delete p;
        }
    }
};

#endif
