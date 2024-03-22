/* Copyright 2024 The WarpX Community
 *
 * Authors: Axel Huebl
 * License: BSD-3-Clause-LBNL
 */

#include "Python/pyWarpX.H"

#include <Particles/NamedComponentParticleContainer.H>

/** Register particle Container class that names components.
 *
 * @tparam T_Allocator the AMReX memory arena (e.g., pinned or default)
 */
template <template<class> class T_Allocator=amrex::DefaultAllocator>
void init_NamedComponentParticleContainer_alloc (py::module& m, std::string type_name)
{
    py::class_<
        NamedComponentParticleContainer<T_Allocator>,
        amrex::ParticleContainerPureSoA<PIdx::nattribs, 0, T_Allocator>
    > npc (m, type_name.c_str());
    npc
        .def("add_real_comp",
            [](NamedComponentParticleContainer<T_Allocator>& pc, const std::string& name, bool comm) {
                pc.AddRealComp(name, comm);
            },
            py::arg("name"), py::arg("comm")
        )
        .def("add_int_comp",
             [](NamedComponentParticleContainer<T_Allocator>& pc, const std::string& name, bool comm) {
                 pc.AddIntComp(name, comm);
             },
             py::arg("name"), py::arg("comm")
        )

        .def_property_readonly("real_comp_names",
            [](NamedComponentParticleContainer<T_Allocator>& pc)
            {
                return pc.getParticleComps();
            }
        )
        .def_property_readonly("int_comp_names",
            [](NamedComponentParticleContainer<T_Allocator>& pc)
            {
                return pc.getParticleiComps();
            }
        )

        .def("get_comp_index",
            [](NamedComponentParticleContainer<T_Allocator>& pc, std::string comp_name)
            {
                auto particle_comps = pc.getParticleComps();
                return particle_comps.at(comp_name);
            },
            py::arg("comp_name")
        )
        .def("get_icomp_index",
            [](NamedComponentParticleContainer<T_Allocator>& pc, std::string comp_name)
            {
                auto particle_comps = pc.getParticleiComps();
                return particle_comps.at(comp_name);
            },
            py::arg("comp_name")
        )
    ;
}

void init_NamedComponentParticleContainer (py::module& m)
{
    init_NamedComponentParticleContainer_alloc<amrex::DefaultAllocator>(m, "NamedComponentParticleContainer");
    init_NamedComponentParticleContainer_alloc<amrex::PinnedArenaAllocator>(m, "PinnedMemoryParticleContainer");
}
