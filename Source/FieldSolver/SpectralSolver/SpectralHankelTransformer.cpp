#include <WarpXConst.H>
#include <SpectralHankelTransformer.H>

SpectralHankelTransformer::SpectralHankelTransformer (int const nr_nodes,
                                                      int const modes,
                                                      amrex::Real const rmax)
: nr(nr_nodes), n_rz_azimuthal_modes(modes)
{

    dht0.resize(n_rz_azimuthal_modes);
    dhtp.resize(n_rz_azimuthal_modes);
    dhtm.resize(n_rz_azimuthal_modes);
    kr.resize(nr*n_rz_azimuthal_modes);

    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {
        dht0[mode].reset( new HankelTransform(mode  , mode, nr, rmax) );
        dhtp[mode].reset( new HankelTransform(mode+1, mode, nr, rmax) );
        dhtm[mode].reset( new HankelTransform(mode-1, mode, nr, rmax) );

        // Save all of the kr's in one place to allow easy access later
        amrex::Real *kr_array = kr.dataPtr();
        auto const & nu = dht0[mode]->getSpectralFrequencies();
        auto const & nu_array = nu.dataPtr();
        amrex::ParallelFor(nr,
        [=] AMREX_GPU_DEVICE (int ir)
        {
            int const ii = ir + mode*nr;
            kr_array[ii] = 2*MathConst::pi*nu_array[ir];
        });
    }

}

// Converts a scalar field from the physical to the spectral space for all modes
void
SpectralHankelTransformer::PhysicalToSpectral_Scalar (amrex::Box const & box,
                                                      amrex::FArrayBox const & F_physical,
                                                      amrex::FArrayBox const & G_spectral)
{
    // The Hankel transform is purely real, so the real and imaginary parts of
    // F can be transformed separately, so a simple loop over components
    // can be done.
    const int nz = box.length(1);
    for (int icomp=0 ; icomp < 2*n_rz_azimuthal_modes-1 ; icomp++) {
        amrex::FArrayBox F_physical_i(F_physical, amrex::make_alias, icomp, 1);
        amrex::FArrayBox G_spectral_i(G_spectral, amrex::make_alias, icomp, 1);
        amrex::Array4<amrex::Real> const & F_physical_i_array = F_physical_i.array();
        amrex::Array4<amrex::Real> const & G_spectral_i_array = G_spectral_i.array();
        int const mode = (icomp + 1)/2;
        dht0[mode]->HankelForwardTransform(nz, F_physical_i_array, G_spectral_i_array);
    }
}

// Converts a vector field from the physical to the spectral space for all modes
void
SpectralHankelTransformer::PhysicalToSpectral_Vector (amrex::Box const & box,
                                                      amrex::FArrayBox const & F_r_physical,
                                                      amrex::FArrayBox const & F_t_physical,
                                                      amrex::FArrayBox const & G_p_spectral,
                                                      amrex::FArrayBox const & G_m_spectral)
{
    const int nz = box.length(1);

    amrex::Array4<const amrex::Real> const & F_r_physical_array = F_r_physical.array();
    amrex::Array4<const amrex::Real> const & F_t_physical_array = F_t_physical.array();

    amrex::FArrayBox temp_p_r(box);
    amrex::FArrayBox temp_p_c(box);
    amrex::FArrayBox temp_m_r(box);
    amrex::FArrayBox temp_m_c(box);

    amrex::Array4<amrex::Real> const & temp_p_r_array = temp_p_r.array();
    amrex::Array4<amrex::Real> const & temp_p_c_array = temp_p_c.array();
    amrex::Array4<amrex::Real> const & temp_m_r_array = temp_m_r.array();
    amrex::Array4<amrex::Real> const & temp_m_c_array = temp_m_c.array();

    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {
        int icomp;
        if (mode == 0) {
            icomp = 0;
        } else {
            icomp = 2*mode - 1;
        }

        if (mode == 0) {
            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::Real const r_real = F_r_physical_array(i,j,0,icomp);
                amrex::Real const r_imag = 0.;
                amrex::Real const t_real = F_t_physical_array(i,j,0,icomp);
                amrex::Real const t_imag = 0.;
                // Combine the values
                temp_p_r_array(i,j,0) = 0.5*(r_real + t_imag);
                temp_p_c_array(i,j,0) = 0.5*(r_imag - t_real);
                temp_m_r_array(i,j,0) = 0.5*(r_real - t_imag);
                temp_m_c_array(i,j,0) = 0.5*(r_imag + t_real);
            });
        } else {
            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                amrex::Real const r_real = F_r_physical_array(i,j,0,icomp);
                amrex::Real const r_imag = F_r_physical_array(i,j,0,icomp+1);
                amrex::Real const t_real = F_t_physical_array(i,j,0,icomp);
                amrex::Real const t_imag = F_t_physical_array(i,j,0,icomp+1);
                // Combine the values
                temp_p_r_array(i,j,0) = 0.5*(r_real + t_imag);
                temp_p_c_array(i,j,0) = 0.5*(r_imag - t_real);
                temp_m_r_array(i,j,0) = 0.5*(r_real - t_imag);
                temp_m_c_array(i,j,0) = 0.5*(r_imag + t_real);
            });
        }

        amrex::FArrayBox G_p_spectral_r(G_p_spectral, amrex::make_alias, icomp  , 1);
        amrex::FArrayBox G_m_spectral_r(G_m_spectral, amrex::make_alias, icomp  , 1);
        amrex::FArrayBox G_p_spectral_c(G_p_spectral, amrex::make_alias, icomp+1, 1);
        amrex::FArrayBox G_m_spectral_c(G_m_spectral, amrex::make_alias, icomp+1, 1);

        auto const & G_p_spectral_r_array = G_p_spectral_r.array();
        auto const & G_p_spectral_c_array = G_p_spectral_c.array();
        auto const & G_m_spectral_r_array = G_m_spectral_r.array();
        auto const & G_m_spectral_c_array = G_m_spectral_c.array();

        dhtp[mode]->HankelForwardTransform(nz, temp_p_r_array, G_p_spectral_r_array);
        dhtp[mode]->HankelForwardTransform(nz, temp_p_c_array, G_p_spectral_c_array);
        dhtm[mode]->HankelForwardTransform(nz, temp_m_r_array, G_m_spectral_r_array);
        dhtm[mode]->HankelForwardTransform(nz, temp_m_c_array, G_m_spectral_c_array);

    }
}

// Converts a scalar field from the spectral to the physical space for all modes
void
SpectralHankelTransformer::SpectralToPhysical_Scalar (amrex::Box const & box,
                                                      amrex::FArrayBox const & G_spectral,
                                                      amrex::FArrayBox const & F_physical)
{
    // The Hankel inverse transform is purely real, so the real and imaginary parts of
    // F can be transformed separately, so a simple loop over components
    // can be done.
    const int nz = box.length(1);
    for (int icomp=0 ; icomp < 2*n_rz_azimuthal_modes-1 ; icomp++) {
        amrex::FArrayBox G_spectral_i(G_spectral, amrex::make_alias, icomp, 1);
        amrex::FArrayBox F_physical_i(F_physical, amrex::make_alias, icomp, 1);
        amrex::Array4<amrex::Real> const & G_spectral_i_array = G_spectral_i.array();
        amrex::Array4<amrex::Real> const & F_physical_i_array = F_physical_i.array();
        int const mode = (icomp + 1)/2;
        dht0[mode]->HankelInverseTransform(nz, G_spectral_i_array, F_physical_i_array);
    }
}

// Converts a vector field from the spectral to the physical space for all modes
void
SpectralHankelTransformer::SpectralToPhysical_Vector (amrex::Box const & box,
                                                      amrex::FArrayBox const& G_p_spectral,
                                                      amrex::FArrayBox const& G_m_spectral,
                                                      amrex::FArrayBox const& F_r_physical,
                                                      amrex::FArrayBox const& F_t_physical)
{
    const int nz = box.length(1);

    amrex::Array4<const amrex::Real> const & G_p_spectral_array = G_p_spectral.array();
    amrex::Array4<const amrex::Real> const & G_m_spectral_array = G_m_spectral.array();

    amrex::FArrayBox temp_r_r(box);
    amrex::FArrayBox temp_t_r(box);
    amrex::FArrayBox temp_r_c(box);
    amrex::FArrayBox temp_t_c(box);

    amrex::Array4<amrex::Real> const & temp_r_r_array = temp_r_r.array();
    amrex::Array4<amrex::Real> const & temp_t_r_array = temp_t_r.array();
    amrex::Array4<amrex::Real> const & temp_r_c_array = temp_r_c.array();
    amrex::Array4<amrex::Real> const & temp_t_c_array = temp_t_c.array();

    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {
        int icomp;
        if (mode == 0) {
            icomp = 0;
        } else {
            icomp = 2*mode - 1;
        }

        // Note that a litte time could be saved by skipping the complex part for mode 0
        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            amrex::Real const p_real = G_p_spectral_array(i,j,0,icomp);
            amrex::Real const m_real = G_m_spectral_array(i,j,0,icomp);
            amrex::Real const p_imag = G_p_spectral_array(i,j,0,icomp+1);
            amrex::Real const m_imag = G_m_spectral_array(i,j,0,icomp+1);
            // Combine the values
            temp_r_r_array(i,j,0) = p_real + m_real;
            temp_t_r_array(i,j,0) = m_imag - p_imag;
            temp_r_c_array(i,j,0) = p_imag + m_imag;
            temp_t_c_array(i,j,0) = p_real - m_real;
        });

        amrex::FArrayBox F_r_physical_r(F_r_physical, amrex::make_alias, icomp  , 1);
        amrex::FArrayBox F_t_physical_r(F_t_physical, amrex::make_alias, icomp  , 1);
        auto const & F_r_physical_r_array = F_r_physical_r.array();
        auto const & F_t_physical_r_array = F_t_physical_r.array();
        dhtp[mode]->HankelInverseTransform(nz, temp_r_r_array, F_r_physical_r_array);
        dhtm[mode]->HankelInverseTransform(nz, temp_t_r_array, F_t_physical_r_array);

        if (mode > 0) {
            amrex::FArrayBox F_r_physical_c(F_r_physical, amrex::make_alias, icomp+1, 1);
            amrex::FArrayBox F_t_physical_c(F_t_physical, amrex::make_alias, icomp+1, 1);
            auto const & F_r_physical_c_array = F_r_physical_c.array();
            auto const & F_t_physical_c_array = F_t_physical_c.array();
            dhtp[mode]->HankelInverseTransform(nz, temp_r_c_array, F_r_physical_c_array);
            dhtm[mode]->HankelInverseTransform(nz, temp_t_c_array, F_t_physical_c_array);
        }

    }
}
