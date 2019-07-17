#ifndef FIELDGATHER_H_
#define FIELDGATHER_H_

#include "ShapeFactors.H"

using namespace amrex;

/* \brief Current Deposition for thread thread_num
 * /param xp, yp, zp   : Pointer to arrays of particle positions.
 * \param wp           : Pointer to array of particle weights.
 * \param uxp uyp uzp  : Pointer to arrays of particle momentum.
 * \param jx_arr       : Array4 of current density, either full array or tile.
 * \param jy_arr       : Array4 of current density, either full array or tile.
 * \param jz_arr       : Array4 of current density, either full array or tile.
 * \param offset       : Index of first particle for which current is deposited
 * \param np_to_gather : Number of particles for which current is deposited.
                        Particles [offset,offset+np_tp_depose] deposit current.
 * \param dx           : 3D cell size
 * \param xyzmin       : Physical lower bounds of domain.
 * \param lo           : Index lower bounds of domain.
 * \param stagger_shift: 0 if nodal, 0.5 if staggered.
 * /param q            : species charge.
 */
template <int depos_order, int lower_in_v>
void doGatherShapeN(const Real * const xp, 
                    const Real * const yp, 
                    const Real * const zp,
                    Real * const Exp, Real * const Eyp, Real * const Ezp,
                    Real * const Bxp, Real * const Byp, Real * const Bzp,
                    const amrex::Array4<const amrex::Real>& ex_arr, 
                    const amrex::Array4<const amrex::Real>& ey_arr, 
                    const amrex::Array4<const amrex::Real>& ez_arr,
                    const amrex::Array4<const amrex::Real>& bx_arr, 
                    const amrex::Array4<const amrex::Real>& by_arr, 
                    const amrex::Array4<const amrex::Real>& bz_arr,
                    const long offset, const long np_to_gather, 
                    const std::array<amrex::Real,3>& dx,
                    const std::array<Real, 3> xyzmin,
                    const Dim3 lo,
                    const amrex::Real stagger_shift)
{
    const Real dxi = 1.0/dx[0];
    const Real dzi = 1.0/dx[2];
#if (AMREX_SPACEDIM == 2)
    const Real invvol = dxi*dzi;
#else // (AMREX_SPACEDIM == 3)
    const Real dyi = 1.0/dx[1];
    const Real invvol = dxi*dyi*dzi;
#endif

    const Real xmin = xyzmin[0];
    const Real ymin = xyzmin[1];
    const Real zmin = xyzmin[2];

    // Loop over particles and deposit into jx_arr, jy_arr and jz_arr
    Print()<<"in my function\n";
    ParallelFor( np_to_gather,
                 [=] AMREX_GPU_DEVICE (long ip) {

                     // --- Compute shape factors
                     // x direction
                     // Get particle position after 1/2 push back in position
                     const Real xmid = (xp[ip]-xmin)*dxi;
                     // Compute shape factors for node-centered quantities
                     Real AMREX_RESTRICT sx [depos_order + 1];
                     // j: leftmost grid point (node-centered) that the particle touches
                     const int j  = compute_shape_factor<depos_order>(sx, xmid);
                     // Compute shape factors for cell-centered quantities
                     Real AMREX_RESTRICT sx0[depos_order + 1 - lower_in_v];
                     // j0: leftmost grid point (cell-centered) that the particle touches
                     const int j0 = compute_shape_factor<depos_order - lower_in_v>(sx0, xmid-stagger_shift);
                     
#if (AMREX_SPACEDIM == 3)
                     // y direction
                     const Real ymid= (yp[ip]-ymin)*dyi;
                     Real AMREX_RESTRICT sy [depos_order + 1];
		     const int k  = compute_shape_factor<depos_order>(sy,  ymid);
                     Real AMREX_RESTRICT sy0[depos_order + 1 - lower_in_v];
		     const int k0 = compute_shape_factor<depos_order-lower_in_v>(sy0, ymid-stagger_shift);
#endif
                     // z direction
                     const Real zmid= (zp[ip]-zmin)*dzi;
                     Real AMREX_RESTRICT sz [depos_order + 1];
		     const int l  = compute_shape_factor<depos_order>(sz,  zmid);
                     Real AMREX_RESTRICT sz0[depos_order + 1 - lower_in_v];
		     const int l0 = compute_shape_factor<depos_order - lower_in_v>(sz0, zmid-stagger_shift);

                     // Deposit current into jx_arr, jy_arr and jz_arr
             Exp[ip] = 0;
             Eyp[ip] = 0;
             Ezp[ip] = 0;
             Bxp[ip] = 0;
             Byp[ip] = 0;
             Bzp[ip] = 0;
#if (AMREX_SPACEDIM == 2)
/*
 -- ex(ip) = ex(ip) + sx0(jj)*sz (ll)*exg(j0+jj, l +ll)
 -- ey(ip) = ey(ip) + sx (jj)*sz (ll)*eyg(j +jj, l +ll)
 -- ez(ip) = ez(ip) + sx (jj)*sz0(ll)*ezg(j +jj, l0+ll)
 -- bx(ip) = bx(ip) + sx (jj)*sz0(ll)*bxg(j +jj, l0+ll)
by(ip) = by(ip) + sx0(jj)*sz0(ll)*byg(j0+jj, l0+ll)
 -- bz(ip) = bz(ip) + sx0(jj)*sz (ll)*bzg(j0+jj, l +ll)
 */
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int ix=0; ix<=depos_order; ix++){
                             Eyp[ip] += sx[ix]*sz[iz]*
                                 ey_arr(lo.x+j+ix, lo.y+l+iz, 0);
                         }
                     }
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                             Exp[ip] += sx0[ix]*sz[iz]*
                                 ex_arr(lo.x+j0+ix, lo.y+l +iz, 0);
                             Bzp[ip] += sx0[ix]*sz[iz]*
                                 bz_arr(lo.x+j0+ix, lo.y+l +iz, 0);
                         }
                     }
                     for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                         for (int ix=0; ix<=depos_order; ix++){
                             Ezp[ip] += sx[ix]*sz0[iz]*
                                 ez_arr(lo.x+j+ix, lo.y+l0 +iz, 0);
                             Bxp[ip] += sx[ix]*sz0[iz]*
                                 bx_arr(lo.x+j+ix, lo.y+l0 +iz, 0);
                         }
                     }
                     for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                         for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                             Byp[ip] += sx0[ix]*sz0[iz]*
                                 by_arr(lo.x+j0+ix, lo.y+l0+iz, 0);
                         }
                     }
                     */
#else // (AMREX_SPACEDIM == 3)
/*
ex(ip) = ex(ip) + sx0(jj)*sy(kk)*sz(ll)*exg(j0+jj, k+kk, l+ll)
ey(ip) = ey(ip) + sx(jj)*sy0(kk)*sz(ll)*eyg(j+jj, k0+kk, l+ll)
ez(ip) = ez(ip) + sx(jj)*sy(kk)*sz0(ll)*ezg(j+jj, k+kk, l0+ll)
bx(ip) = bx(ip) + sx(jj)*sy0(kk)*sz0(ll)*bxg(j+jj, k0+kk, l0+ll)
by(ip) = by(ip) + sx0(jj)*sy(kk)*sz0(ll)*byg(j0+jj, k+kk, l0+ll)
bz(ip) = bz(ip) + sx0(jj)*sy0(kk)*sz(ll)*bzg(j0+jj, k0+kk, l+ll)
 */
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int iy=0; iy<=depos_order; iy++){
                             for (int ix=0; ix<=depos_order-lower_in_v; ix++){
                                 
                             }
                         }
                     }
                                 /*
                                 amrex::Gpu::Atomic::Add(
                                     &jx_arr(lo.x+j0+ix, lo.y+k +iy, lo.z+l +iz),
                                     sx0[ix]*sy [iy]*sz [iz]*wqx);
                             }
                         }
                     }
                     for (int iz=0; iz<=depos_order; iz++){
                         for (int iy=0; iy<=depos_order-lower_in_v; iy++){
                             for (int ix=0; ix<=depos_order; ix++){
                                 amrex::Gpu::Atomic::Add(
                                     &jy_arr(lo.x+j +ix, lo.y+k0+iy, lo.z+l +iz), 
                                     sx [ix]*sy0[iy]*sz [iz]*wqy);
                             }
                         }
                     }
                     for (int iz=0; iz<=depos_order-lower_in_v; iz++){
                         for (int iy=0; iy<=depos_order; iy++){
                             for (int ix=0; ix<=depos_order; ix++){
                                 amrex::Gpu::Atomic::Add(
                                     &jz_arr(lo.x+j +ix, lo.y+k +iy, lo.z+l0+iz),
                                     sx [ix]*sy [iy]*sz0[iz]*wqz);
                                 */
                             }
                         }
                     }
#endif
                 }
        );
}

#endif // FIELDGATHER_H_
