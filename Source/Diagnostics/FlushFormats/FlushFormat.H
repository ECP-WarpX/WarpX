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
        if ( pc->DoFieldIonization()) {
            tmp.AddIntComp(false);
        }
#ifdef WARPX_QED
        if ( pc->has_breit_wheeler() ) {
            tmp.AddRealComp(false);
        }
        if ( pc->has_quantum_sync() ) {
            tmp.AddRealComp(false);
        }
        return tmp;
#endif
    }

}

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
