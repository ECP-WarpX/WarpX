
/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_PARTICLES_RADIATION_H
#define WARPX_PARTICLES_RADIATION_H

#include "RadiationHandler_fwd.H"

#include "Particles/Pusher/GetAndSetPosition.H"
#include "Particles/Sorting/SortingUtils.H"
#include "WarpX.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>


#include <memory>
#include <string>

/* \brief CollisionHandler is a class that contains the
 * calculation of radiation for particles at each frequencies and angles.
 */
class RadiationHandler
{
public:
    RadiationHandler ();

    /* Perform the calculation of the radiation */
    //void doRadiation (amrex::Real dt, MultiParticleContainer* mypc);
    void add_radiation_contribution(const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc, amrex::Real current_time);
    void add_detector();
    void gather_and_write_radiation(const std::string& filename);
    void Integral_overtime();


private:

    // Frequency range
    amrex::Vector<amrex::Real> m_omega_range;

    // Dimensions of the detector
    amrex::Vector<amrex::Real> m_theta_range;

    int m_omega_points;

    //put a detector
    int m_get_a_detector;
    // Resolution of the detector
    amrex::Vector<amrex::Real>  m_d_theta;
    double m_d_omega;
    //The vector with the resolution of the detector
    amrex::Vector<amrex::Real> m_d_det;
    amrex::Real m_det_distance;

    //Define the Fab with the datas
    int ncomp;

    amrex::Vector<int> m_det_pts;
    amrex::Vector<int>  m_det_direction;
    amrex::Vector<amrex::Real> m_pos_det;
    amrex::Vector<amrex::Real> n;

    //resolution of the grid of the detectir
    amrex::Vector<amrex::Real> m_d_d;

    //
    amrex::Gpu::DeviceVector<amrex::Real> m_radiation_data;
    amrex::Gpu::DeviceVector<amrex::Real> m_radiation_calculation;

    //
    amrex::Real dephas;

    //
    amrex::Vector<amrex::Vector<amrex::Real>> det_bornes;
    amrex::Vector<amrex::Vector<amrex::Real>> det_pos;
    amrex::Vector<amrex::Real> omega_calc;

};
#endif // WARPX_PARTICLES_RADIATION_H
