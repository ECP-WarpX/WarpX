#ifndef WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_
#define WARPX_PARTICLES_GATHER_GETEXTERNALFIELDS_H_

#include "Particles/Pusher/GetAndSetPosition.H"

#include "Particles/WarpXParticleContainer_fwd.H"
#include "Utils/Algorithms/LinearInterpolation.H"
#include "Utils/WarpXConst.H"

#include "AcceleratorLattice/LatticeElementFinder.H"

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Parser.H>
#include <AMReX_REAL.H>

#include <optional>

#include <fstream>


struct ExternalFieldFromFile3DCyl 
{
    const amrex::RealVect dx;
    const amrex::RealVect origin;
    const amrex::Array4<amrex::Real const> arrayr;
    const amrex::Array4<amrex::Real const> arrayz;
    const amrex::Gpu::DeviceVector<amrex::Real> datar;
    const amrex::Gpu::DeviceVector<amrex::Real> dataz;
    public:
    ExternalFieldFromFile3DCyl(
        const amrex::RealVect dx_,
        const amrex::RealVect origin_,
        const amrex::Dim3 shape_,
        amrex::Gpu::DeviceVector<amrex::Real>& datar_,
        amrex::Gpu::DeviceVector<amrex::Real>& dataz_
        ) :
        dx(dx_), 
        origin(origin_), 
        arrayr(datar_.data(), {0,0,0}, shape_, 1),
        arrayz(dataz_.data(), {0,0,0}, shape_, 1),
        datar(std::move(datar_)), 
        dataz(std::move(dataz_))
    {
        assert(arrayr.begin.x == arrayz.begin.x);
        assert(arrayr.begin.y == arrayz.begin.y);
        assert(arrayr.begin.z == arrayz.begin.z);
        assert(arrayr.end.x == arrayz.end.x);
        assert(arrayr.end.y == arrayz.end.y);
        assert(arrayr.end.z == arrayz.end.z);
    }

    amrex::Real interp_(const amrex::Array4<amrex::Real const> array, int ir, int iz, amrex::Real r, amrex::Real z) const noexcept {
        const amrex::Real
                    f00 = array(0, iz  , ir  ),
                    f01 = array(0, iz  , ir+1),
                    f10 = array(0, iz+1, ir  ),
                    f11 = array(0, iz+1, ir+1);
        
        return static_cast<amrex::Real>(utils::algorithms::bilinear_interp<double>
                    (ir * dx[0], (ir+1) * dx[0], iz * dx[1], (iz+1) * dx[1],
                     f00, f10, f01, f11,
                     r, z));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::RealVect interpolate(amrex::Real x, amrex::Real y, amrex::Real z) const noexcept {
        auto r = std::sqrt(x*x + y*y);
        amrex::Real rel_r = r - origin[0];
        amrex::Real rel_z = z - origin[1];
        // Get index of the external field array
        int ir = floor( rel_r / dx[0] );
        int iz = floor( rel_z / dx[1] );

        if (ir < arrayr.begin.z || iz < arrayr.begin.y || (ir + 1) >= arrayr.end.z || (iz + 1) >= arrayr.end.y)
            return amrex::RealVect();

        auto Br = interp_(arrayr, ir, iz, rel_r, rel_z);
        auto Bz = interp_(arrayz, ir, iz, rel_r, rel_z);

        amrex::RealVect retvect = {
            Br * x / r,
            Br * y / r, 
            Bz
        };

#ifdef DEBUG_EXTERNAL_FIELD_INTERPOLATION
        { 
            std::ofstream myfile;
            myfile.open ("/tmp/ptcls-interpolate-rz.csv", std::ios_base::app);
            myfile << ir << " " << iz << " " << r << " " << z << " " << Br << " " << Bz << std::endl;
            myfile.close();   
        }
#endif

        return retvect;
    }
};


/** \brief Functor class that assigns external
 *         field values (E and B) to particles.
*/
struct GetExternalEBField
{
    enum ExternalFieldInitType { None, Parser, RepeatedPlasmaLens, ReadFromFile, Unknown };

    GetExternalEBField () = default;

    GetExternalEBField (const WarpXParIter& a_pti, long a_offset = 0) noexcept;

    ExternalFieldInitType m_Etype;
    ExternalFieldInitType m_Btype;

    amrex::ParticleReal m_gamma_boost;
    amrex::ParticleReal m_uz_boost;

    amrex::ParserExecutor<4> m_Exfield_partparser;
    amrex::ParserExecutor<4> m_Eyfield_partparser;
    amrex::ParserExecutor<4> m_Ezfield_partparser;
    amrex::ParserExecutor<4> m_Bxfield_partparser;
    amrex::ParserExecutor<4> m_Byfield_partparser;
    amrex::ParserExecutor<4> m_Bzfield_partparser;

    GetParticlePosition<PIdx> m_get_position;
    amrex::Real m_time;

    amrex::ParticleReal m_repeated_plasma_lens_period;
    const amrex::ParticleReal* AMREX_RESTRICT m_repeated_plasma_lens_starts = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_repeated_plasma_lens_lengths = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_repeated_plasma_lens_strengths_E = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_repeated_plasma_lens_strengths_B = nullptr;
    int m_n_lenses;
    amrex::Real m_dt;
    const amrex::ParticleReal* AMREX_RESTRICT m_ux = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uy = nullptr;
    const amrex::ParticleReal* AMREX_RESTRICT m_uz = nullptr;

    const ExternalFieldFromFile3DCyl* Bfield_file_external_particle_cyl = nullptr;

    std::optional<LatticeElementFinderDevice> d_lattice_element_finder;

    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isNoOp () const { return (m_Etype == None && m_Btype == None && !d_lattice_element_finder.has_value()); }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator () (long i,
                      amrex::ParticleReal& field_Ex,
                      amrex::ParticleReal& field_Ey,
                      amrex::ParticleReal& field_Ez,
                      amrex::ParticleReal& field_Bx,
                      amrex::ParticleReal& field_By,
                      amrex::ParticleReal& field_Bz) const noexcept
    {
        using namespace amrex::literals;

        if (d_lattice_element_finder) {
            // Note that the "*" is needed since d_lattice_element_finder is optional
            (*d_lattice_element_finder)(i, field_Ex, field_Ey, field_Ez,
                                           field_Bx, field_By, field_Bz);
        }

        if (m_Etype == None && m_Btype == None) { return; }

        amrex::ParticleReal Ex = 0._prt;
        amrex::ParticleReal Ey = 0._prt;
        amrex::ParticleReal Ez = 0._prt;
        amrex::ParticleReal Bx = 0._prt;
        amrex::ParticleReal By = 0._prt;
        amrex::ParticleReal Bz = 0._prt;

        constexpr amrex::ParticleReal inv_c2 = 1._prt/(PhysConst::c*PhysConst::c);

        if (m_Etype == ExternalFieldInitType::Parser)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);
            amrex::Real lab_time = m_time;
            if (m_gamma_boost > 1._prt) {
                lab_time = m_gamma_boost*m_time + m_uz_boost*z*inv_c2;
                z = m_gamma_boost*z + m_uz_boost*m_time;
            }
            Ex = m_Exfield_partparser((amrex::ParticleReal) x, (amrex::ParticleReal) y, (amrex::ParticleReal) z, lab_time);
            Ey = m_Eyfield_partparser((amrex::ParticleReal) x, (amrex::ParticleReal) y, (amrex::ParticleReal) z, lab_time);
            Ez = m_Ezfield_partparser((amrex::ParticleReal) x, (amrex::ParticleReal) y, (amrex::ParticleReal) z, lab_time);
        }

        if (m_Btype == ExternalFieldInitType::Parser)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);
            amrex::Real lab_time = m_time;
            if (m_gamma_boost > 1._prt) {
                lab_time = m_gamma_boost*m_time + m_uz_boost*z*inv_c2;
                z = m_gamma_boost*z + m_uz_boost*m_time;
            }
            Bx = m_Bxfield_partparser(x, y, z, lab_time);
            By = m_Byfield_partparser(x, y, z, lab_time);
            Bz = m_Bzfield_partparser(x, y, z, lab_time);
        }
        else if (m_Btype == ExternalFieldInitType::ReadFromFile)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);
            if (m_gamma_boost > 1._prt) {
                z = m_gamma_boost*z + m_uz_boost*m_time;
            }

            auto B = Bfield_file_external_particle_cyl->interpolate(x, y, z);

            Bx = B[0];
            By = B[1];
            Bz = B[2];
        }

        if (m_Etype == RepeatedPlasmaLens ||
            m_Btype == RepeatedPlasmaLens)
        {
            amrex::ParticleReal x, y, z;
            m_get_position(i, x, y, z);

            const amrex::ParticleReal uxp = m_ux[i];
            const amrex::ParticleReal uyp = m_uy[i];
            const amrex::ParticleReal uzp = m_uz[i];

            const amrex::ParticleReal gamma = std::sqrt(1._prt + (uxp*uxp + uyp*uyp + uzp*uzp)*inv_c2);
            const amrex::ParticleReal vzp = uzp/gamma;

            // the current slice in z between now and the next time step
            amrex::ParticleReal zl = z;
            amrex::ParticleReal zr = z + vzp*m_dt;

            if (m_gamma_boost > 1._prt) {
                zl = m_gamma_boost*zl + m_uz_boost*m_time;
                zr = m_gamma_boost*zr + m_uz_boost*(m_time + m_dt);
            }

            // the plasma lens periods do not start before z=0
            if (zl > 0) {
                // find which is the next lens
                int i_lens = static_cast<int>(std::floor(zl/m_repeated_plasma_lens_period));
                if (i_lens < m_n_lenses) {
                    amrex::ParticleReal const lens_start = m_repeated_plasma_lens_starts[i_lens] + i_lens*m_repeated_plasma_lens_period;
                    amrex::ParticleReal const lens_end = lens_start + m_repeated_plasma_lens_lengths[i_lens];

                    // Calculate the residence correction
                    // frac will be 1 if the step is completely inside the lens, between 0 and 1
                    // when entering or leaving the lens, and otherwise 0.
                    // This accounts for the case when particles step over the element without landing in it.
                    // This assumes that vzp > 0.
                    amrex::ParticleReal const zl_bounded = std::min(std::max(zl, lens_start), lens_end);
                    amrex::ParticleReal const zr_bounded = std::min(std::max(zr, lens_start), lens_end);
                    amrex::ParticleReal const frac = ((zr - zl) == 0._rt ? 1._rt : (zr_bounded - zl_bounded)/(zr - zl));

                    // Note that "+=" is used since the fields may have been set above
                    // if a different E or Btype was specified.
                    Ex += x*frac*m_repeated_plasma_lens_strengths_E[i_lens];
                    Ey += y*frac*m_repeated_plasma_lens_strengths_E[i_lens];
                    Bx += +y*frac*m_repeated_plasma_lens_strengths_B[i_lens];
                    By += -x*frac*m_repeated_plasma_lens_strengths_B[i_lens];
                }
            }

        }

        if (m_gamma_boost > 1._prt) {
            // Transform the fields to the boosted frame
            const amrex::ParticleReal Ex_boost = m_gamma_boost*Ex - m_uz_boost*By;
            const amrex::ParticleReal Ey_boost = m_gamma_boost*Ey + m_uz_boost*Bx;
            const amrex::ParticleReal Bx_boost = m_gamma_boost*Bx + m_uz_boost*Ey*inv_c2;
            const amrex::ParticleReal By_boost = m_gamma_boost*By - m_uz_boost*Ex*inv_c2;
            Ex = Ex_boost;
            Ey = Ey_boost;
            Bx = Bx_boost;
            By = By_boost;
        }

        field_Ex += Ex;
        field_Ey += Ey;
        field_Ez += Ez;
        field_Bx += Bx;
        field_By += By;
        field_Bz += Bz;

    }
};

#endif
