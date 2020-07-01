#include "BackTransformFunctor.H"
#include "WarpX.H"
using namespace amrex;

BackTransformFunctor::BackTransformFunctor (amrex::MultiFab const * mf_src, int lev,
                                            const int ncomp, const int num_buffers,
                                            const amrex::IntVect crse_ratio)
    : ComputeDiagFunctor(ncomp, crse_ratio), m_mf_src(mf_src), m_lev(lev), m_num_buffers(num_buffers)
{
    InitData();
}

void
BackTransformFunctor::operator ()(amrex::MultiFab& mf_dst, int dcomp, const int i_buffer) const
{
    amrex::Print() << " we are in BTFunctor operator \n";
    auto& warpx = WarpX::GetInstance();
    amrex::Real gamma_boost = warpx.gamma_boost;
    int moving_window_dir = warpx.moving_window_dir;
    amrex::Real beta_boost = std::sqrt( 1._rt - 1._rt/( gamma_boost * gamma_boost) );
    // Get the right slice of each field in the CC MultiFab, transform it and store it in the output.
    // 1. ncomp, scomp = 0, dcomp = dcomp, bool interpolate = true
    // 2. slice = amrex::get_slice_data
    amrex::MultiFab slice;
    LorentzTransformZ( slice, gamma_boost, beta_boost);
    
}

void
BackTransformFunctor::PrepareFunctorData (int i_buffer,
                          bool ZSliceInDomain, amrex::Real current_z_boost, 
                          amrex::Box buffer_box )
{
    m_buffer_box[i_buffer] = buffer_box;
    m_current_z_boost[i_buffer] = current_z_boost;
    m_perform_backtransform[i_buffer] = 0;
    if (ZSliceInDomain) m_perform_backtransform[i_buffer] = 1;
}

void
BackTransformFunctor::InitData ()
{

     m_buffer_box.resize( m_num_buffers );
     m_current_z_boost.resize( m_num_buffers );
     m_perform_backtransform.resize( m_num_buffers );     
}

void
BackTransformFunctor::LorentzTransformZ (amrex::MultiFab& data, amrex::Real gamma_boost,
                                         amrex::Real beta_boost) const
{
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(data, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& tbx = mfi.tilebox();
        amrex::Array4< amrex::Real > arr = data[mfi].array();
        amrex::Real clight = PhysConst::c;
        amrex::Real inv_clight = 1.0_rt/clight;
        // arr(x,y,z,comp) has ten-components namely,
        // Ex Ey Ez Bx By Bz jx jy jz rho in that order.
        amrex::ParallelFor( tbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                // Back-transform the transverse electric and magnetic fields.
                // Note that the z-components, Ez, Bz, are not changed by the transform.
                amrex::Real e_lab, b_lab, j_lab, rho_lab;
                // Transform Ex_boost (ncomp=0) & By_boost (ncomp=4) to lab-frame
                e_lab = gamma_boost * ( arr(i, j, k, 0)
                                        + beta_boost * clight * arr(i, j, k, 4) );
                b_lab = gamma_boost * ( arr(i, j, k, 4)
                                        + beta_boost * inv_clight * arr(i, j, k, 0) );
                // Store lab-frame data in-place
                arr(i, j, k, 0) = e_lab;
                arr(i, j, k, 4) = b_lab;
                
                // Transform Ey_boost (ncomp=1) & Bx_boost (ncomp=3) to lab-frame
                e_lab = gamma_boost * ( arr(i, j, k, 1)
                                        - beta_boost * clight * arr(i, j, k, 3) );
                b_lab = gamma_boost * ( arr(i, j, k, 3)
                                        - beta_boost * inv_clight * arr(i, j, k, 1) );
                // Store lab-frame data in-place
                arr(i, j, k, 1) = e_lab;
                arr(i, j, k, 3) = b_lab;
           
                // Transform charge density (ncomp=9)
                // and z-component of current density (ncomp=8)
                j_lab = gamma_boost * ( arr(i, j, k, 8)
                                        + beta_boost * clight * arr(i, j, k, 9) );
                rho_lab = gamma_boost * ( arr(i, j, k, 9)
                                          + beta_boost * inv_clight * arr(i, j, k, 8) );
                // Store lab-frame jz and rho in-place 
                arr(i, j, k, 8) = j_lab;
                arr(i, j, k, 9) = rho_lab;


                
                
            }      
        );
    }    

}
