#include <complex>
#include <iostream>
#include <complex>
#include <vector>

#include <AMReX.H>
#include <iostream>
#include <vector>
#include <AMReX_REAL.H>
#include <Utils/WarpX_Complex.H>
#include <AMReX_AmrCore.H>

class LaserEnvelope {
public:
    LaserEnvelope(amrex::AmrCore* amr_core, int ispecies, const std::string& name);
    amrex::Vector< amrex::FabArray< amrex::BaseFab <Complex> > > A_laser_envelope;
    void AllocateLaserEnvelope (int lev, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::IntVect& ngA, const bool aux_is_nodal);
    amrex::Real m_direction;
    amrex::Real m_e_max;
    amrex::Real m_polarization;
    amrex::Real m_position;
    amrex::Real m_profile_duration;
    amrex::Real m_profile_focal_distance;
    amrex::Real m_profile;
    amrex::Real m_profile_t_peak;
    amrex::Real m_profile_waist;
    amrex::Real m_wavelength;
    const std::string& m_laser_name;
    };
