#ifndef WARPX_COMM_K_H_
#define WARPX_COMM_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_x (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Bxa,
                            amrex::Array4<amrex::Real const> const& Bxf,
                            amrex::Array4<amrex::Real const> const& Bxc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;
    Bxa(j,k,l) = owx * Bxc(jg,kg,lg) + wx * Bxc(jg+1,kg,lg) + Bxf(j,k,l);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_y (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Bya,
                            amrex::Array4<amrex::Real const> const& Byf,
                            amrex::Array4<amrex::Real const> const& Byc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    // Note that for 2d, l=0, because the amrex convention is used here.

#if (AMREX_SPACEDIM == 3)
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;
    Bya(j,k,l) = owy * Byc(jg,kg,lg) + wy * Byc(jg,kg+1,lg) + Byf(j,k,l);
#else
    Bya(j,k,l) = Byc(jg,kg,lg) + Byf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_z (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Bza,
                            amrex::Array4<amrex::Real const> const& Bzf,
                            amrex::Array4<amrex::Real const> const& Bzc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    // Note that for 2d, l=0, because the amrex convention is used here.

#if (AMREX_SPACEDIM == 3)
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;
    Bza(j,k,l) = owz * Bzc(jg,kg,lg) + owz * Bzc(jg,kg,lg+1) + Bzf(j,k,l);
#else
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;
    Bza(j,k,l) = owy * Bzc(jg,kg,lg) + owy * Bzc(jg,kg+1,lg) + Bzf(j,k,l);
#endif
}

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




AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_x (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Exa,
                            amrex::Array4<amrex::Real const> const& Exf,
                            amrex::Array4<amrex::Real const> const& Exc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    // Cell-centered in x: determine whether we will average with
    // the coarse point to the right or left of the nearest coarse point
    int j0 = (j == 2*jg) ? -1 : 1;
    // Nodal in y: determine if the coarse point is at the same position
    // as a fine point
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 3)
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;
    Exa(j,k,l) =  0.75*( owy * owz * Exc(jg  ,kg  ,lg  )
                       +  wy * owz * Exc(jg  ,kg+1,lg  )
                       + owy *  wz * Exc(jg  ,kg  ,lg+1)
                       +  wy *  wz * Exc(jg  ,kg+1,lg+1) )
                + 0.25*( owy * owz * Exc(jg+j0,kg  ,lg  )
                       +  wy * owz * Exc(jg+j0,kg+1,lg  )
                       + owy *  wz * Exc(jg+j0,kg  ,lg+1)
                       +  wy *  wz * Exc(jg+j0,kg+1,lg+1) )
        + Exf(j,k,l);
#else
    Exa(j,k,l) =  0.75*( owy * Exc(jg   ,kg  ,lg  )
                       +  wy * Exc(jg   ,kg+1,lg  ) )
                + 0.25*( owy * Exc(jg+j0,kg  ,lg  )
                       +  wy * Exc(jg+j0,kg+1,lg  ) )
        + Exf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_y (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Eya,
                            amrex::Array4<amrex::Real const> const& Eyf,
                            amrex::Array4<amrex::Real const> const& Eyc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

#if (AMREX_SPACEDIM == 3)
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;
    Eya(j,k,l) = owx * owz * Eyc(jg  ,kg  ,lg  )
        +         wx * owz * Eyc(jg+1,kg  ,lg  )
        +        owx *  wz * Eyc(jg  ,kg  ,lg+1)
        +         wx *  wz * Eyc(jg+1,kg  ,lg+1)
        + Eyf(j,k,l);
#else
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;
    Eya(j,k,l) = owx * owy * Eyc(jg  ,kg  ,lg)
        +         wx * owy * Eyc(jg+1,kg  ,lg)
        +        owx *  wy * Eyc(jg  ,kg+1,lg)
        +         wx *  wy * Eyc(jg+1,kg+1,lg)
        + Eyf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_z (int j, int k, int l,
                            amrex::Array4<amrex::Real> const& Eza,
                            amrex::Array4<amrex::Real const> const& Ezf,
                            amrex::Array4<amrex::Real const> const& Ezc)
{
    using namespace amrex;

    int lg = amrex::coarsen(l,2);
    int kg = amrex::coarsen(k,2);
    int jg = amrex::coarsen(j,2);

    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

#if (AMREX_SPACEDIM == 3)
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;
    Eza(j,k,l) = owx * owy * Ezc(jg  ,kg  ,lg  )
        +         wx * owy * Ezc(jg+1,kg  ,lg  )
        +        owx *  wy * Ezc(jg  ,kg+1,lg  )
        +         wx *  wy * Ezc(jg+1,kg+1,lg  )
        + Ezf(j,k,l);
#else
    Eza(j,k,l) = owx * Ezc(jg,kg,lg) + wx * Ezc(jg+1,kg,lg) + Ezf(j,k,l);
#endif
}

#endif
