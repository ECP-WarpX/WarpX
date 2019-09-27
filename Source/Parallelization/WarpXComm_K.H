#ifndef WARPX_COMM_K_H_
#define WARPX_COMM_K_H_

#include <AMReX_FArrayBox.H>

template <bool is_nodal>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int getInterpolationShift( int i_fine, int i_coarse ){
    return 0;
}
template <bool is_nodal>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real getInterpolationWeight( int i_fine, int i_coarse ){
    return 0;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int getInterpolationShift <true> ( int i_fine, int i_coarse ){
    return +1;
}
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real getInterpolationWeight <true> ( int i_fine, int i_coarse ){
    amrex::Real w = (i_fine == i_coarse) ? 0.0 : 0.5;
    return w;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
int getInterpolationShift <false> ( int i_fine, int i_coarse ){
    int i0 = (i_fine == 2*i_coarse) ? -1 : 1;
    return i0;
}
template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real getInterpolationWeight <false> ( int i_fine, int i_coarse ){
    return amrex::Real(0.75);
}


template< bool is_nodal_x, bool is_nodal_y, bool is_nodal_z >
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void CalculateAuxiliaryField( int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Fa,
                            amrex::Array4<amrex::Real const> const& Ff,
                            amrex::Array4<amrex::Real const> const& Fc )
{
    using namespace amrex;

    // j, k, l are the indices on the fine grid
    // Get the corresponding indices on the coarse grid
    // For cell-centered grid: this is the nearest coarse point
    // For nodal grid: this is either the same point, or the point to the left
    int jg = amrex::coarsen(j,2);
    int kg = amrex::coarsen(k,2);
    int lg = amrex::coarsen(l,2);

    // Find the index of the other nearest coarse point,
    // and the corresponding weights, in each direction
    int const j0 = getInterpolationShift<is_nodal_x>( j, jg );
    Real const wx = getInterpolationWeight<is_nodal_x>( j, jg );
    Real const w0x = Real(1.) - wx;

    int const k0 = getInterpolationShift<is_nodal_x>( k, kg );
    Real const wy = getInterpolationWeight<is_nodal_x>( k, kg );
    Real const w0y = Real(1.) - wy;

#if (AMREX_SPACEDIM == 3)
    int const l0 = getInterpolationShift<is_nodal_x>( l, lg );
    Real const wz = getInterpolationWeight<is_nodal_x>( l, lg );
    Real const w0z = Real(1.) - wz;
#endif

    Fa(j,k,l) = Ff(j,k,l) +
    // Compute weighted interpolation of the coarse field,
    // by averaging between the two closest coarse points in each direction
#if (AMREX_SPACEDIM == 3)
               wz*(
                        wy*(   wx*Fc(jg   ,kg   ,lg  )
                            + w0x*Fc(jg+j0,kg   ,lg  ) )
                     + w0y*(   wx*Fc(jg   ,kg+k0,lg  )
                            + w0x*Fc(jg+j0,kg+k0,lg  ) )
                   )
             + w0z*(
                        wy*(   wx*Fc(jg   ,kg   ,lg+l0)
                            + w0x*Fc(jg+j0,kg   ,lg+l0) )
                     + w0y*(   wx*Fc(jg   ,kg+k0,lg+l0)
                            + w0x*Fc(jg+j0,kg+k0,lg+l0) )
                   );
#else
                        wy*(   wx*Fc(jg   ,kg   ,lg  )
                            + w0x*Fc(jg+j0,kg   ,lg  ) )
                     + w0y*(   wx*Fc(jg   ,kg+k0,lg  )
                            + w0x*Fc(jg+j0,kg+k0,lg  ) );
#endif
}

#endif //WARPX_COMM_K_H_
