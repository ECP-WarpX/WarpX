#include <complex>
#include <iostream>
#include <complex>
#include <vector>


#include <iostream>
#include <vector>

class LaserEnvelope {
private:
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > A_potential_vector;
};

/**
 * Gaussian laser profile
 */
class GaussianLaserProfile : public ILaserProfile
{

public:
    void
    init (
        const amrex::ParmParse& ppl,
        CommonLaserParameters params) override final;

    //No update needed
    void
    update (amrex::Real /*t */) override final {}

    void
    fill_amplitude (
        const int np,
        amrex::Real const * AMREX_RESTRICT const Xp,
        amrex::Real const * AMREX_RESTRICT const Yp,
        amrex::Real t,
        amrex::Real * AMREX_RESTRICT const amplitude) const override final;

private:
    struct {
        amrex::Real waist          = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real duration       = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real t_peak         = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real focal_distance = std::numeric_limits<amrex::Real>::quiet_NaN();
        amrex::Real zeta           = 0;
        amrex::Real beta           = 0;
        amrex::Real phi2           = 0;
        amrex::Real phi0           = 0;

        amrex::Vector<amrex::Real> stc_direction; //! Direction of the spatio-temporal couplings
        amrex::Real theta_stc; //! Angle between polarization (p_X) and direction of spatiotemporal coupling (stc_direction)
    } m_params;

    CommonLaserParameters m_common_params;
};