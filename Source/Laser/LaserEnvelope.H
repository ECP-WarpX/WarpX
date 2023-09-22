/* Copyright 2023 The WarpX Community
 *
 * This file is part of WarpX.
 *
 * Authors: Adam Bensoubaya, Edoardo Zoni, Remi Lehe
 * License: BSD-3-Clause-LBNL
 */
#include <complex>
#include <iostream>
#include <complex>
#include <vector>

#include "LaserEnvelope_fwd.H"
#include "WarpX.H"

#include <AMReX.H>
#include <iostream>
#include <vector>
#include <AMReX_REAL.H>
#include <Utils/WarpX_Complex.H>
#include <AMReX_AmrCore.H>

class LaserEnvelope
{
public:

    /**
     * \brief Constructor
     */
    LaserEnvelope (const int nlevs_max);

    /**
     * TODO
     */
    void ReadParameters ();

    /**
     * TODO
     */
    void AllocateMFs (const int nlevs_max);

    /**
     * TODO
     */
    void AllocateLevelMFs (
        const int lev,
        const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm,
        const int ncomps,
        const amrex::IntVect& ngA,
        const amrex::IntVect& A_nodal_flag);

    /**
     * TODO
     */
    void ClearLevel (const int lev);

    /**
     * TODO
     */
    Complex FillAmplitude (const amrex::Real x, const amrex::Real y, const amrex::Real z);

    /**
     * TODO
     */
    void InitData (const int finestLevel);

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> A_laser_envelope;

    amrex::Vector<amrex::Real> m_direction;
    amrex::Real m_e_max;
    amrex::Real m_wavelength;
    amrex::Real m_polarization;
    amrex::Real m_position;
    amrex::Real m_profile_duration;
    amrex::Real m_profile;
    amrex::Real m_waist;
    amrex::Real m_duration;
    amrex::Real m_z_peak;
    amrex::Real m_theta_stc;
    amrex::Real m_focal_distance;
    amrex::Real m_zeta = 0._rt;
    amrex::Real m_beta = 0._rt;
    amrex::Real m_phi2 = 0._rt;
    amrex::Real m_phi0 = 0._rt;
};
