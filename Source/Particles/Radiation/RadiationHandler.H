
/* Copyright 2023 Thomas Clark, Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef WARPX_PARTICLES_RADIATION_H_
#define WARPX_PARTICLES_RADIATION_H_
#include "Utils/Math/Complex.H"
#include "Particles/WarpXParticleContainer.H"
#include "Particles/Pusher/GetAndSetPosition.H"


#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Box.H>
#include <AMReX_Vector.H>
#include <AMReX_FArrayBox.H>

#include <memory>
#include <string>

/* \brief CollisionHandler is a light weight class that contains the
 * calculation of radiation for particles at each frequencies and angles.
 */
class RadiationHandler
{
public:
    RadiationHandler ();

    /* Perform the calculation of the radiation */
    //void doRadiation (amrex::Real dt, MultiParticleContainer* mypc);    
    void add_radiation_contribution(const amrex::Real dt, std::unique_ptr<WarpXParticleContainer>& pc);
    void add_detector(amrex::IntVect m_det_pts, amrex::IntVect m_omega_rg, amrex::Vector<amrex::Real> m_d_det, int m_omega_pts, double m_d_distance, amrex::IntVect m_d_direction, amrex::Vector<amrex::Real> m_th_rg);
    void keepoldmomentum(std::unique_ptr<WarpXParticleContainer>& pc);

private:

    // Frequency range
    amrex::IntVect m_omega_range;
 
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
    int m_det_distance;

    //Define the Fab with the datas
    int ncomp;

    amrex::IntVect m_det_pts;
    amrex::IntVect m_det_direction;
    amrex::Vector<amrex::Real> m_pos_det;
};
#endif // WARPX_PARTICLES_RADIATION_H_

