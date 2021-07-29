#ifndef WARPX_FLUSHFORMAT_H_
#define WARPX_FLUSHFORMAT_H_

#include "FlushFormat_fwd.H"

#include "Diagnostics/ParticleDiag/ParticleDiag.H"
#include "Particles/MultiParticleContainer.H"
#include "WarpX.H"

#include <AMReX_AmrParticles.H>

namespace ParticleBuffer {
    template <template<class> class Allocator>
    using BufferType = amrex::AmrParticleContainer<0, 0, PIdx::nattribs, 0, Allocator>;

    template <template<class> class Allocator>
    BufferType<Allocator> getTmpPC (const WarpXParticleContainer* pc)
    {
        BufferType<Allocator> tmp(&WarpX::GetInstance());
        // add runtime real comps to tmp
        for (int ic = 0; ic < pc->NumRuntimeRealComps(); ++ic) { tmp.AddRealComp(false); }

        // add runtime int comps to tmp
        for (int ic = 0; ic < pc->NumRuntimeIntComps(); ++ic) { tmp.AddIntComp(false); }

        return tmp;
    }

}

=======
>>>>>>> 1e198d47f7d9c6a8a991a2ddc89159677653eb1f
class FlushFormat
{
public:
    /** Flush fields and particles to file */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
        const std::string prefix, const int file_min_digits,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        bool plot_raw_rho, bool plot_raw_F,
        bool isBTD = false, int snapshotID = -1,
        const amrex::Geometry& full_BTD_snapshot = amrex::Geometry(),
        bool isLastBTDFlush = false) const = 0;

     virtual ~FlushFormat() {}
};

#endif // WARPX_FLUSHFORMAT_H_
