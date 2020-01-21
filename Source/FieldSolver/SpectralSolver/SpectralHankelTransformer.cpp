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

    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {
        dht0[mode].reset( new HankelTransform(mode  , mode, nr, rmax) );
        dhtp[mode].reset( new HankelTransform(mode+1, mode, nr, rmax) );
        dhtm[mode].reset( new HankelTransform(mode-1, mode, nr, rmax) );
    }

    ExtractKrArray();

}

// This needs to be separate since the ParallelFor cannot be in the constructor.
void
SpectralHankelTransformer::ExtractKrArray ()
{
    kr.resize(nr*n_rz_azimuthal_modes);

    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {

        // Save all of the kr's in one place to allow easy access later
        amrex::Real *kr_array = kr.dataPtr();
        auto const & kr_m = dht0[mode]->getSpectralWavenumbers();
        auto const & kr_m_array = kr_m.dataPtr();
        int const nr_temp = nr;
        amrex::ParallelFor(nr,
        [=] AMREX_GPU_DEVICE (int ir)
        {
            int const ii = ir + mode*nr_temp;
            kr_array[ii] = kr_m_array[ir];
        });
    }
}
// Converts a scalar field from the physical to the spectral space for all modes
void
SpectralHankelTransformer::PhysicalToSpectral_Scalar (amrex::Box const & box,
                                                      amrex::FArrayBox const & F_physical,
                                                      amrex::FArrayBox       & G_spectral)
{
    // The Hankel transform is purely real, so the real and imaginary parts of
    // F can be transformed separately, so a simple loop over components
    // can be done.
    // Note that F_physical does not include the imaginary part of mode 0,
    // but G_spectral does.
    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {
        int const mode_r = 2*mode;
        int const mode_i = 2*mode + 1;
        if (mode == 0) {
            int const icomp = 0;
            dht0[mode]->HankelForwardTransform(F_physical, icomp, G_spectral, mode_r);
            G_spectral.setVal(0., mode_i);
        } else {
            int const icomp = 2*mode - 1;
            dht0[mode]->HankelForwardTransform(F_physical, icomp  , G_spectral, mode_r);
            dht0[mode]->HankelForwardTransform(F_physical, icomp+1, G_spectral, mode_i);
        }
    }
}

// Converts a vector field from the physical to the spectral space for all modes
void
SpectralHankelTransformer::PhysicalToSpectral_Vector (amrex::Box const & box,
                                                      amrex::FArrayBox & F_r_physical,
                                                      amrex::FArrayBox & F_t_physical,
                                                      amrex::FArrayBox & G_p_spectral,
                                                      amrex::FArrayBox & G_m_spectral)
{
    // Note that F and G include the imaginary part of mode 0
    // F will be overwritten

    amrex::Array4<amrex::Real> const & F_r_physical_array = F_r_physical.array();
    amrex::Array4<amrex::Real> const & F_t_physical_array = F_t_physical.array();

    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {

        int const mode_r = 2*mode;
        int const mode_i = 2*mode + 1;

        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            amrex::Real const r_real = F_r_physical_array(i,j,k,mode_r);
            amrex::Real const r_imag = F_r_physical_array(i,j,k,mode_i);
            amrex::Real const t_real = F_t_physical_array(i,j,k,mode_r);
            amrex::Real const t_imag = F_t_physical_array(i,j,k,mode_i);
            // Combine the values
            // temp_p = (F_r - I*F_t)/2
            // temp_m = (F_r + I*F_t)/2
            F_r_physical_array(i,j,k,mode_r) = 0.5*(r_real + t_imag);
            F_r_physical_array(i,j,k,mode_i) = 0.5*(r_imag - t_real);
            F_t_physical_array(i,j,k,mode_r) = 0.5*(r_real - t_imag);
            F_t_physical_array(i,j,k,mode_i) = 0.5*(r_imag + t_real);
        });

        amrex::Gpu::streamSynchronize();

        dhtp[mode]->HankelForwardTransform(F_r_physical, mode_r, G_p_spectral, mode_r);
        dhtp[mode]->HankelForwardTransform(F_r_physical, mode_i, G_p_spectral, mode_i);
        dhtm[mode]->HankelForwardTransform(F_t_physical, mode_r, G_m_spectral, mode_r);
        dhtm[mode]->HankelForwardTransform(F_t_physical, mode_i, G_m_spectral, mode_i);

    }
}

// Converts a scalar field from the spectral to the physical space for all modes
void
SpectralHankelTransformer::SpectralToPhysical_Scalar (amrex::Box const & box,
                                                      amrex::FArrayBox const & G_spectral,
                                                      amrex::FArrayBox       & F_physical)
{
    // The Hankel inverse transform is purely real, so the real and imaginary parts of
    // F can be transformed separately, so a simple loop over components
    // can be done.
    // Note that F_physical does not include the imaginary part of mode 0,
    // but G_spectral does.

    amrex::Gpu::streamSynchronize();

    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {
        int const mode_r = 2*mode;
        int const mode_i = 2*mode + 1;
        if (mode == 0) {
            int const icomp = 0;
            dht0[mode]->HankelInverseTransform(G_spectral, mode_r, F_physical, icomp);
        } else {
            int const icomp = 2*mode - 1;
            dht0[mode]->HankelInverseTransform(G_spectral, mode_r, F_physical, icomp);
            dht0[mode]->HankelInverseTransform(G_spectral, mode_i, F_physical, icomp+1);
        }
    }
}

// Converts a vector field from the spectral to the physical space for all modes
void
SpectralHankelTransformer::SpectralToPhysical_Vector (amrex::Box const & box,
                                                      amrex::FArrayBox const& G_p_spectral,
                                                      amrex::FArrayBox const& G_m_spectral,
                                                      amrex::FArrayBox      & F_r_physical,
                                                      amrex::FArrayBox      & F_t_physical)
{
    // Note that F and G include the imaginary part of mode 0

    amrex::Array4<amrex::Real> const & F_r_physical_array = F_r_physical.array();
    amrex::Array4<amrex::Real> const & F_t_physical_array = F_t_physical.array();

    for (int mode=0 ; mode < n_rz_azimuthal_modes ; mode++) {

        int const mode_r = 2*mode;
        int const mode_i = 2*mode + 1;

        amrex::Gpu::streamSynchronize();

        dhtp[mode]->HankelInverseTransform(G_p_spectral, mode_r, F_r_physical, mode_r);
        dhtp[mode]->HankelInverseTransform(G_p_spectral, mode_i, F_r_physical, mode_i);
        dhtm[mode]->HankelInverseTransform(G_m_spectral, mode_r, F_t_physical, mode_r);
        dhtm[mode]->HankelInverseTransform(G_m_spectral, mode_i, F_t_physical, mode_i);

        amrex::Gpu::streamSynchronize();

        amrex::ParallelFor(box,
        [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            amrex::Real const p_real = F_r_physical_array(i,j,k,mode_r);
            amrex::Real const p_imag = F_r_physical_array(i,j,k,mode_i);
            amrex::Real const m_real = F_t_physical_array(i,j,k,mode_r);
            amrex::Real const m_imag = F_t_physical_array(i,j,k,mode_i);
            // Combine the values
            // F_r =    G_p + G_m
            // F_t = I*(G_p - G_m)
            F_r_physical_array(i,j,k,mode_r) =  p_real + m_real;
            F_r_physical_array(i,j,k,mode_i) =  p_imag + m_imag;
            F_t_physical_array(i,j,k,mode_r) = -p_imag + m_imag;
            F_t_physical_array(i,j,k,mode_i) =  p_real - m_real;
        });

    }
}
