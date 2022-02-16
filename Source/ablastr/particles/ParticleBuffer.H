/* Copyright 2021-2022 Andrew Myers, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef ABLASTR_PARTICLEBUFFER_H_
#define ABLASTR_PARTICLEBUFFER_H_

#include <AMReX_AmrCore.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_ParGDB.H>

#include <type_traits>


namespace ablastr {
namespace particles {

/** Create an amrex::AmrParticleContainer with different Allocator
 *
 * This creates a new AMReX particle container type with same attributes
 * as an existing particle container. But at the same time, it can change
 * its allocator. This is helpful when creating temporary particle buffers
 * for filter operations and device-to-host copies.
 *
 * @tparam Allocator AMReX allocator, e.g., amrex::PinnedArenaAllocator
 */
template <template<class> class Allocator=amrex::DefaultAllocator>
struct ParticleBuffer {
    /** type trait to translate one particle container to another, with changed allocator */
    template <typename T_PC>
    using ContainerLike = amrex::ParticleContainer<T_PC::NStructReal, T_PC::NStructInt, T_PC::NArrayReal, T_PC::NArrayInt, Allocator>;

    /** create an empty particle container
     *
     * Has the same compile-time and runtime components as an existing particle
     * container.
     *
     * @param pc the particle container to task as blueprint
     * @return an empty particle container
     * */
    template <typename T_PC>
    static ContainerLike<std::remove_cv_t<T_PC> >
    make_container (T_PC * pc)
    {
        auto gdb = const_cast<amrex::ParGDBBase*>(pc->GetParGDB());
        ContainerLike<std::remove_cv_t<T_PC> > tmp(gdb);

        // add runtime real comps to tmp
        for (int ic = 0; ic < pc->NumRuntimeRealComps(); ++ic) { tmp.AddRealComp(false); }

        // add runtime int comps to tmp
        for (int ic = 0; ic < pc->NumRuntimeIntComps(); ++ic) { tmp.AddIntComp(false); }

        return tmp;
    }
};

} // namespace particles
} // namespace ablastr

#endif // ABLASTR_PARTICLEBUFFER_H_
