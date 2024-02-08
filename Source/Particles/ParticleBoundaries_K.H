/* Copyright 2021 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBOUNDARIES_K_H_
#define PARTICLEBOUNDARIES_K_H_

#include "ParticleBoundaries.H"

#include <AMReX_AmrCore.H>

namespace ApplyParticleBoundaries {

    /* \brief Applies the boundary condition on a specific axis
     *        This is called by apply_boundaries.
     */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    apply_boundary (amrex::ParticleReal& x, amrex::Real xmin, amrex::Real xmax,
                    bool& change_sign_ux, bool& particle_lost,
                    ParticleBoundaryType xmin_bc, ParticleBoundaryType xmax_bc,
                    amrex::Real refl_probability_xmin, amrex::Real refl_probability_xmax,
                    amrex::RandomEngine const& engine )
    {
        if (x < xmin) {
            if (xmin_bc == ParticleBoundaryType::Open) {
                particle_lost = true;
            }
            else if (xmin_bc == ParticleBoundaryType::Absorbing) {
                if (refl_probability_xmin == 0 || amrex::Random(engine) > refl_probability_xmin) {
                    particle_lost = true;
                }
                else
                {
                    x = 2*xmin - x;
                    change_sign_ux = true;
                }
            }
            else if (xmin_bc == ParticleBoundaryType::Reflecting) {
                x = 2*xmin - x;
                change_sign_ux = true;
            }
        }
        else if (x > xmax) {
            if (xmax_bc == ParticleBoundaryType::Open) {
                particle_lost = true;
            }
            else if (xmax_bc == ParticleBoundaryType::Absorbing) {
                if (refl_probability_xmax == 0 || amrex::Random(engine) > refl_probability_xmax) {
                    particle_lost = true;
                }
                else
                {
                    x = 2*xmax - x;
                    change_sign_ux = true;
                }
            }
            else if (xmax_bc == ParticleBoundaryType::Reflecting) {
                x = 2*xmax - x;
                change_sign_ux = true;
            }
        }
    }


namespace {
    /** Return u sampled according to the probability distribution:
      * p(u) \propto u \exp(-(u-u_m)^2/2u_th^2)
      *
      * @param u_m Central momentum
      * @param u_th Momentum spread
      * @param engine Object used to generate random numbers
      */
    [[nodiscard]]
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    amrex::Real
    generateGaussianFluxDist( amrex::Real u_m, amrex::Real u_th, amrex::RandomEngine const& engine ) {

        using namespace amrex::literals;

        // Momentum to be returned at the end of this function
        amrex::Real u = 0._rt;

        const amrex::Real abs_u_m = std::abs(u_m);

        if (u_th == 0._rt) {
            u = u_m; // Trivial case ; avoids division by 0 in the rest of the code below
        } else if (abs_u_m < 0.6*u_th) {
            // Mean velocity magnitude is less than thermal velocity
            // Use the distribution u*exp(-u**2*(1-abs(u_m)/u_th)/(2*u_th**2)) as an approximation
            // and then use the rejection method to correct it
            // ( stop rejecting with probability exp(-abs(u_m)/(2*u_th**3)*(u-sign(u_m)*u_th)**2) )
            // Note that this is the method that is used in the common case u_m=0
            const amrex::Real umsign = std::copysign(1._rt, u_m);
            const amrex::Real approx_u_th = u_th/std::sqrt( 1._rt - abs_u_m/u_th );
            const amrex::Real reject_prefactor = (abs_u_m/u_th)/(2._rt*u_th*u_th); // To save computation
            bool reject = true;
            while (reject) {
                // Generates u according to u*exp(-u**2/(2*approx_u_th**2)),
                // using the method of the inverse cumulative function
                amrex::Real xrand = 1._rt - amrex::Random(engine); // ensures urand > 0
                u = approx_u_th * std::sqrt(2._rt*std::log(1._rt/xrand));
                // Rejection method
                xrand = amrex::Random(engine);
                if (xrand < std::exp(-reject_prefactor*(u - umsign*u_th)*(u - umsign*u_th))) { reject = false; }
            }
        } else {
            // Mean velocity magnitude is greater than thermal velocity
            // Use the distribution exp(-(u-u_m-u_th**2/abs(u_m))**2/(2*u_th**2)) as an approximation
            // and then use the rejection method to correct it
            // ( stop rejecting with probability (u/abs(u_m))*exp(1-(u/abs(u_m))) ; note
            // that this number is always between 0 and 1 )
            // Note that in the common case `u_m = 0`, this rejection method
            // is not used, and the above rejection method is used instead.
            bool reject = true;
            const amrex::Real approx_u_m = u_m + u_th*u_th/abs_u_m;
            const amrex::Real inv_um = 1._rt/abs_u_m; // To save computation
            while (reject) {
                // Approximate distribution: normal distribution, where we only retain positive u
                u = -1._rt;
                while (u < 0) {
                    u = amrex::RandomNormal(approx_u_m, u_th, engine);
                }
                // Rejection method
                const amrex::Real xrand = amrex::Random(engine);
                if (xrand < u*inv_um* std::exp(1._rt - u*inv_um)) { reject = false; }
            }
        }

        return u;
    }
}

    /* \brief Applies absorbing or reflecting boundary condition to the input particles, along all axis.
     *        For reflecting boundaries, the position of the particle is changed appropriately and
     *        the sign of the velocity is changed (depending on the reflect_all_velocities flag).
     *        For absorbing, a flag is set whether the particle has been lost (it is up to the calling
     *        code to take appropriate action to remove any lost particles). Absorbing boundaries can
     *        be given a reflection coefficient for stochastic reflection of particles, this
     *        coefficient is zero by default.
     *        Note that periodic boundaries are handled in AMReX code.
     *
     * \param x, xmin, xmax: particle x position, location of x boundary
     * \param y, ymin, ymax: particle y position, location of y boundary (3D only)
     * \param z, zmin, zmax: particle z position, location of z boundary
     * \param ux, uy, uz: particle momenta
     * \param particle_lost: output, flags whether the particle was lost
     * \param boundaries: object with boundary condition settings
    */
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void
    apply_boundaries (
#ifndef WARPX_DIM_1D_Z
              amrex::ParticleReal& x, amrex::Real xmin, amrex::Real xmax,
#endif
#if (defined WARPX_DIM_3D) || (defined WARPX_DIM_RZ)
                      amrex::ParticleReal& y,
#endif
#if (defined WARPX_DIM_3D)
                      amrex::Real ymin, amrex::Real ymax,
#endif
                      amrex::ParticleReal& z, amrex::Real zmin, amrex::Real zmax,
                      amrex::ParticleReal& ux, amrex::ParticleReal& uy, amrex::ParticleReal& uz,
                      bool& particle_lost,
                      ParticleBoundaries::ParticleBoundariesData const& boundaries,
                      amrex::RandomEngine const& engine)
    {
        bool change_sign_ux = false;
        bool change_sign_uy = false;
        bool change_sign_uz = false;

#ifndef WARPX_DIM_1D_Z
        apply_boundary(x, xmin, xmax, change_sign_ux, particle_lost,
                       boundaries.xmin_bc, boundaries.xmax_bc,
                       boundaries.reflection_model_xlo(-ux), boundaries.reflection_model_xhi(ux),
                       engine);
#endif
#ifdef WARPX_DIM_3D
        apply_boundary(y, ymin, ymax, change_sign_uy, particle_lost,
                       boundaries.ymin_bc, boundaries.ymax_bc,
                       boundaries.reflection_model_ylo(-uy), boundaries.reflection_model_yhi(uy),
                       engine);
#endif
        apply_boundary(z, zmin, zmax, change_sign_uz, particle_lost,
                       boundaries.zmin_bc, boundaries.zmax_bc,
                       boundaries.reflection_model_zlo(-uz), boundaries.reflection_model_zhi(uz),
                       engine);

        if (boundaries.reflect_all_velocities && (change_sign_ux | change_sign_uy | change_sign_uz)) {
            change_sign_ux = true;
            change_sign_uy = true;
            change_sign_uz = true;
        }
#ifdef WARPX_DIM_RZ
        // Note that the reflection of the position does "r = 2*rmax - r", but this is only approximate.
        // The exact calculation requires the position at the start of the step.
        if (change_sign_ux && change_sign_uy) {
            ux = -ux;
            uy = -uy;
        } else if (change_sign_ux) {
            // Reflect only ur
            // Note that y is theta
            amrex::Real ur = ux*std::cos(y) + uy*std::sin(y);
            const amrex::Real ut = -ux*std::sin(y) + uy*std::cos(y);
            ur = -ur;
            ux = ur*std::cos(y) - ut*std::sin(y);
            uy = ur*std::sin(y) + ut*std::cos(y);
        }
#else
        ux = PhysConst::c*amrex::RandomNormal(0, boundaries.u_th, engine);
        uy = PhysConst::c*amrex::RandomNormal(0, boundaries.u_th, engine);
        // Generate particles with a Gaussian flux distribution in the direction normal to the boundary
        if (change_sign_ux) { ux = -PhysConst::c*generateGaussianFluxDist(0, boundaries.u_th, engine ) * ( (ux>0)-(ux<0) ); }
        if (change_sign_uy) { uy = -PhysConst::c*generateGaussianFluxDist(0, boundaries.u_th, engine ) * ( (uy>0)-(uy<0) ); }
#endif
        uz = PhysConst::c*amrex::RandomNormal(0, boundaries.u_th, engine);
        if (change_sign_uz) { uz = -PhysConst::c*generateGaussianFluxDist(0, boundaries.u_th, engine ) * ( (uz>0)-(uz<0) ); }
    }

}
#endif
