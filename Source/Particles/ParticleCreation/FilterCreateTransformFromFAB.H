#ifndef FILTER_CREATE_TRANSFORM_FROM_FAB_H_
#define FILTER_CREATE_TRANSFORM_FROM_FAB_H_

template <int N, typename DstTile, typename FAB, typename Index,
          typename CreateFunc, typename TransFunc,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
Index filterCreateTransformFromFAB (DstTile& dst1, DstTile& dst2, amrex::Box box,
                                     FAB *src_FAB, Index* mask,
                                    Index dst1_index, Index dst2_index,
                                    CreateFunc&& create, TransFunc&& transform) noexcept
{
    using namespace amrex;

    auto ncells = box.volume();
    if (ncells == 0) return 0;

    const auto box_sizex = box.length(0);
    const auto box_sizey = box.length(1);

    auto & warpx = WarpX::GetInstance();
    const int level_0 = 0;
    Geometry const & geom = warpx.Geom(level_0);

    Real xlo_global = geom.ProbLo(0);
    Real ylo_global = geom.ProbLo(1);
    Real zlo_global = geom.ProbLo(2);
    Real dx         = geom.CellSize(0);
    Real dy         = geom.CellSize(1);
    Real dz         = geom.CellSize(2);

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
            IntVect iv = box.atOffset(i);
            int j = iv[0];
            int k = iv[1];
            int l = iv[2];
            Real x = xlo_global + j*dx; // I think this should correspond to node positions right ?
#if (AMREX_SPACEDIM == 3)
            Real y = ylo_global + k*dy;
            Real z = zlo_global + l*dz;
#else
            Real y = 0.;
            Real z = ylo_global + k*dy;
#endif
            for (int n = 0; n < N; ++j)
            {
                create(dst1_data, N*p_offsets[i] + dst1_index + n, x, y, z);
                create(dst2_data, N*p_offsets[i] + dst2_index + n, x, y, z);
            }
            transform(dst1_data, N*p_offsets[i] + dst1_index, N,
                       arrNumPartCreation(j,k,l));
            transform(dst2_data, N*p_offsets[i] + dst2_index, N,
                       arrNumPartCreation(j,k,l));
        }
    });

    Gpu::streamSynchronize();
    return num_added;
}

template <int N, typename DstTile, typename FABs, typename Index,
          typename PredFunc, typename CreateFunc, typename TransFunc>
Index filterCreateTransformFromFAB (DstTile& dst1, DstTile& dst2, amrex::Box box,
                                    const FABs *src_FABs, Index dst1_index,
                                    Index dst2_index,  PredFunc&& filter,
                                    CreateFunc&& create,
                                    TransFunc && transform) noexcept
{
    using namespace amrex;

    // FAB representing the number of particles to create
    FArrayBox NumPartCreation(box, 1);
    auto arrNumPartCreation = NumPartCreation.array();

    const Dim3 lo = lbound(box);
    const Dim3 hi = ubound(box);

    const auto box_sizex = box.length(0);
    const auto box_sizey = box.length(1);

    auto ncells = box.volume();
    if (ncells == 0) return 0;

    Gpu::DeviceVector<Index> mask(ncells);

    auto p_mask = mask.dataPtr();

    amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i, int j, int k){
        arrNumPartCreation(i,j,k) = filter(src_FABs,i,j,k);
        IntVect iv(AMREX_D_DECL(i,j,k));
        const auto mask_position = box.index(iv);
        p_mask[mask_position] = (arrNumPartCreation(i,j,k) > 0);
    });

    return filterCreateTransformFromFAB<N>(dst1, dst2, box, &NumPartCreation,
                                        mask.dataPtr(), dst1_index, dst2_index,
                                        std::forward<CreateFunc>(create),
                                        std::forward<TransFunc>(transform));
}

#endif
