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
    LaserEnvelope(const int nlevs_max);

    /**
     * TODO
     */
    void ReadParameters();

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
    void InitData ();

    amrex::Vector<std::unique_ptr<amrex::MultiFab>> A_laser_envelope;

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
};
