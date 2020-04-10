/* Copyright 2020 Luca Fedeli, Neil Zaim
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef FILTER_CREATE_TRANSFORM_FROM_FAB_H_
#define FILTER_CREATE_TRANSFORM_FROM_FAB_H_

/**
 * \brief Apply a filter on a list of FABs, then create and apply a transform
 * operation to the particles depending on the output of the filter.
 *
 * \tparam N number of particles created in the dst(s) in each cell
 * \tparam DstTile the dst particle tile type
 * \tparam FAB the src FAB type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam CreateFunc1 the create function type for dst1
 * \tparam CreateFunc2 the create function type for dst2
 * \tparam TransFunc the transform function type
 *
 * \param[in,out] dst1 the first destination tile
 * \param[in,out] dst2 the second destination tile
 * \param[in] box the box where the particles are created
 * \param[in] src_FAB A FAB defined on box that is used in the transform function
 * \param[in] mask pointer to the mask - 1 means create, 0 means don't create
 * \param[in] dst1_index the location at which to starting writing the result to dst1
 * \param[in] dst2_index the location at which to starting writing the result to dst2
 * \param[in] create1 callable that defines what will be done for the create step for dst1.
 * \param[in] create2 callable that defines what will be done for the create step for dst2.
 * \param[in] transform callable that defines the transformation to apply on dst1 and dst2.
 *
 * \return num_added the number of particles that were written to dst 1 and dst 2.
 */
template <int N, typename DstTile, typename FAB, typename Index,
          typename CreateFunc1, typename CreateFunc2, typename TransFunc,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
Index filterCreateTransformFromFAB (DstTile& dst1, DstTile& dst2, const amrex::Box box,
                                const FAB *src_FAB, const Index* mask,
                                const Index dst1_index, const Index dst2_index,
                                CreateFunc1&& create1, CreateFunc2&& create2,
                                TransFunc&& transform) noexcept
{
    using namespace amrex;

    const auto ncells = box.volume();
    if (ncells == 0) return 0;

    const auto box_sizex = box.length(0);
    const auto box_sizey = box.length(1);

    auto & warpx = WarpX::GetInstance();
    const int level_0 = 0;
    Geometry const & geom = warpx.Geom(level_0);

    Real xlo_global = geom.ProbLo(0);
    Real ylo_global = geom.ProbLo(1);
    Real dx         = geom.CellSize(0);
    Real dy         = geom.CellSize(1);
#if (AMREX_SPACEDIM == 3)
    Real zlo_global = geom.ProbLo(2);
    Real dz         = geom.CellSize(2);
#endif

    auto arrNumPartCreation = src_FAB->array();
    Gpu::DeviceVector<Index> offsets(ncells);
    Gpu::exclusive_scan(mask, mask+ncells, offsets.begin());

    Index last_mask, last_offset;
    Gpu::copyAsync(Gpu::deviceToHost, mask+ncells-1, mask + ncells, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+ncells-1,
                                      offsets.data()+ncells, &last_offset);

    Gpu::streamSynchronize();

    const Index num_added = N*(last_mask + last_offset);
    dst1.resize(std::max(dst1_index + num_added, dst1.numParticles()));
    dst2.resize(std::max(dst2_index + num_added, dst2.numParticles()));

    auto p_offsets = offsets.dataPtr();

    const auto dst1_data = dst1.getParticleTileData();
    const auto dst2_data = dst2.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( ncells, i,
    {
        if (mask[i])
        {
            const IntVect iv = box.atOffset(i);
            const int j = iv[0];
            const int k = iv[1];
    // Currently all particles are created on nodes. This makes it useless
    // to use N>1 (for now).
            const Real x = xlo_global + j*dx;
#if (AMREX_SPACEDIM == 3)
            const int l = iv[2];
            const Real y = ylo_global + k*dy;
            const Real z = zlo_global + l*dz;
#else
            const int l = 0;
            const int y = 0;
            const Real z = ylo_global + k*dy;
#endif
            for (int n = 0; n < N; ++n)
            {
                create1(dst1_data, N*p_offsets[i] + dst1_index + n, x, y, z);
                create2(dst2_data, N*p_offsets[i] + dst2_index + n, x, y, z);
            }
            transform(dst1_data, dst2_data, N*p_offsets[i] + dst1_index,
                    N*p_offsets[i] + dst2_index, N, arrNumPartCreation(j,k,l));
        }
    });

    Gpu::streamSynchronize();
    return num_added;
}

/**
 * \brief Apply a filter on a list of FABs, then create and apply a transform
 * operation to the particles depending on the output of the filter.
 *
 * \tparam N number of particles created in the dst(s) in each cell
 * \tparam DstTile the dst particle tile type
 * \tparam FABs the src array of Array4 type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam FilterFunc the filter function type
 * \tparam CreateFunc1 the create function type for dst1
 * \tparam CreateFunc2 the create function type for dst2
 * \tparam TransFunc the transform function type
 *
 * \param[in,out] dst1 the first destination tile
 * \param[in,out] dst2 the second destination tile
 * \param[in] box the box where the particles are created
 * \param[in] src_FABs An Array of Array4 (e.g. EM fields) defined on box on which
 *            the filter operation is applied/
 * \param[in] dst1_index the location at which to starting writing the result to dst1
 * \param[in] dst2_index the location at which to starting writing the result to dst2
 * \param[in] filter a callable returning a value > 0 if particles are to be created
 *            in the considered cell.
 * \param[in] create1 callable that defines what will be done for the create step for dst1.
 * \param[in] create2 callable that defines what will be done for the create step for dst2.
 * \param[in] transform callable that defines the transformation to apply on dst1 and dst2.
 *
 * \return num_added the number of particles that were written to dst 1 and dst 2.
 */
template <int N, typename DstTile, typename FABs, typename Index,
          typename FilterFunc, typename CreateFunc1, typename CreateFunc2,
           typename TransFunc>
Index filterCreateTransformFromFAB (DstTile& dst1, DstTile& dst2, const amrex::Box box,
                                const FABs& src_FABs, const Index dst1_index,
                                const Index dst2_index, FilterFunc&& filter,
                                CreateFunc1&& create1, CreateFunc2&& create2,
                                TransFunc && transform) noexcept
{
    using namespace amrex;

    FArrayBox NumPartCreation(box, 1);
    Elixir tmp_eli = NumPartCreation.elixir(); // Do we actually need this elixir?
                            // The tests on GPU seem to be working well without it.
    auto arrNumPartCreation = NumPartCreation.array();

    const Dim3 lo = lbound(box);
    const Dim3 hi = ubound(box);

    const auto box_sizex = box.length(0);
    const auto box_sizey = box.length(1);

    const auto ncells = box.volume();
    if (ncells == 0) return 0;

    Gpu::DeviceVector<Index> mask(ncells);

    auto p_mask = mask.dataPtr();

    amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i, int j, int k){
        arrNumPartCreation(i,j,k) = filter(src_FABs,i,j,k);
        const IntVect iv(AMREX_D_DECL(i,j,k));
        const auto mask_position = box.index(iv);
        p_mask[mask_position] = (arrNumPartCreation(i,j,k) > 0);
    });

    return filterCreateTransformFromFAB<N>(dst1, dst2, box, &NumPartCreation,
                                        mask.dataPtr(), dst1_index, dst2_index,
                                        std::forward<CreateFunc1>(create1),
                                        std::forward<CreateFunc2>(create2),
                                        std::forward<TransFunc>(transform));
}

#endif // FILTER_CREATE_TRANSFORM_FROM_FAB_H_
