/* Copyright 2021 Andrew Myers
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef PARTICLEBUFFER_H_
#define PARTICLEBUFFER_H_

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

#endif /*PARTICLEBUFFER_H_*/
