
/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_PARTICLES_RADIATION_H
#define WARPX_PARTICLES_RADIATION_H

#include "RadiationHandler_fwd.H"

#include "Particles/WarpXParticleContainer_fwd.H"

#include <ablastr/math/Complex.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <memory>
#include <string>

/* \brief CollisionHandler is a class that contains the
 * calculation of radiation for particles at each frequencies and angles.
 */
class RadiationHandler
{
public:
    RadiationHandler (const amrex::Array<amrex::Real,3>& center);

    /* Perform the calculation of the radiation */
    //void doRadiation (amrex::Real dt, MultiParticleContainer* mypc);
    void add_radiation_contribution(const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc, amrex::Real current_time);

    void gather_and_write_radiation(const std::string& filename);

private:

    // Frequency range
    amrex::Array<amrex::Real,2> m_omega_range;
    int m_omega_points;
    amrex::Gpu::DeviceVector<amrex::Real> m_omegas;

    // Dimensions of the detector
    amrex::Array<amrex::Real,2> m_theta_range;

    amrex::Real m_det_distance;
    amrex::Array<int,2> m_det_pts;
    amrex::Array<amrex::Real,3> m_det_direction;
    amrex::Array<amrex::Real,3> m_det_orientation;

    amrex::Gpu::DeviceVector<amrex::Real> det_pos_x;
    amrex::Gpu::DeviceVector<amrex::Real> det_pos_y;
    amrex::Gpu::DeviceVector<amrex::Real> det_pos_z;

    amrex::Gpu::DeviceVector<ablastr::math::Complex> m_radiation_data;

};
#endif // WARPX_PARTICLES_RADIATION_H
