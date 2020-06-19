/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef QED_PAIR_GENERATION_H_
#define QED_PAIR_GENERATION_H_

#include "Utils/WarpXConst.H"
#include "Particles/WarpXParticleContainer.H"
#include "Particles/Gather/GetExternalFields.H"
#include "Particles/Gather/FieldGather.H"
#include "Particles/Pusher/GetAndSetPosition.H"
#include "QEDInternals/BreitWheelerEngineWrapper.H"

/** @file
 *
 * This file contains the implementation of the elementary process
 * functors needed for Breit-Wheeler pair generation (one photon generates
 * and electron-positron pair).
 */

/**
 * \brief Filter functor for the Breit Wheeler process
 */
class PairGenerationFilterFunc
{
public:

    /**
    * \brief Constructor of the PairGenerationFilterFunc functor.
    *
    * @param[in] opt_depth_runtime_comp index of the optical depth component
    */
    PairGenerationFilterFunc(int const opt_depth_runtime_comp)
        : m_opt_depth_runtime_comp(opt_depth_runtime_comp)
        {}

    /**
    * \brief Functor call. This method determines if a given (photon) particle
    * should undergo pair generation.
    *
    * @param[in] ptd particle tile data
    * @param[in] i particle index
    * @return true if a pair has to be generated, false otherwise
    */
    template <typename PData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool operator() (const PData& ptd, int const i) const noexcept
    {
        using namespace amrex;

        const amrex::ParticleReal opt_depth =
            ptd.m_runtime_rdata[m_opt_depth_runtime_comp][i];
        return (opt_depth < 0.0_rt);
    }

private:
    int m_opt_depth_runtime_comp = 0; /*!< Index of the optical depth component of the species.*/
};

/**
 * \brief Transform functor for the Breit-Wheeler process
 */
class PairGenerationTransformFunc
{
public:

    /**
    * \brief Constructor of the PairGenerationTransformFunc functor.
    *
    * A BreitWheelerGeneratePairs functor is passed by value. However, it contains
    * only few integer and real parameters and few pointers to the raw data of the
    * lookup tables. Therefore, it should be rather lightweight to copy.
    *
    * @param[in] generate_functor functor to be called to determine the properties of the generated pairs
    */
    PairGenerationTransformFunc(BreitWheelerGeneratePairs const generate_functor,
                                const WarpXParIter& a_pti, int lev, int ngE,
                                amrex::FArrayBox const& exfab,
                                amrex::FArrayBox const& eyfab,
                                amrex::FArrayBox const& ezfab,
                                amrex::FArrayBox const& bxfab,
                                amrex::FArrayBox const& byfab,
                                amrex::FArrayBox const& bzfab,
                                amrex::Array<amrex::Real,3> v_galilean,
                                int a_offset = 0);

    /**
    * \brief Functor call. It determines the properties of the generated pair
    * and it sets to -1 the id of the source photon
    *
    * @param[in,out] dst1 target species 1 (either electrons or positrons)
    * @param[in,out] dst2 target species 2 (either electrons or positrons)
    * @param[in] src source species (photons)
    * @param[in] i_src particle index of the source species
    * @param[in] i_dst1 particle index of target species 1
    * @param[in] i_dst2 particle index of target species 2
    */
    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst1, DstData& dst2, SrcData& src,
        int const i_src, int const i_dst1, int const i_dst2) const noexcept
    {
        using namespace amrex;

        constexpr ParticleReal me = PhysConst::m_e;
        constexpr ParticleReal one_over_me = 1._prt/me;

        const ParticleReal w = src.m_rdata[PIdx::w][i_src];
        const ParticleReal ux = src.m_rdata[PIdx::ux][i_src];
        const ParticleReal uy = src.m_rdata[PIdx::uy][i_src];
        const ParticleReal uz = src.m_rdata[PIdx::uz][i_src];

        // gather E and B
        amrex::ParticleReal xp, yp, zp;
        m_get_position(i_src, xp, yp, zp);

        amrex::ParticleReal ex = 0._rt, ey = 0._rt, ez = 0._rt;
        m_get_externalE(i_src, ex, ey, ez);

        amrex::ParticleReal bx = 0._rt, by = 0._rt, bz = 0._rt;
        m_get_externalB(i_src, bx, by, bz);

        doGatherShapeN(xp, yp, zp, ex, ey, ez, bx, by, bz,
                       m_ex_arr, m_ey_arr, m_ez_arr, m_bx_arr, m_by_arr, m_bz_arr,
                       m_ex_type, m_ey_type, m_ez_type, m_bx_type, m_by_type, m_bz_type,
                       m_dx_arr, m_xyzmin_arr, m_lo, m_n_rz_azimuthal_modes,
                       m_nox, m_l_lower_order_in_v);

        const auto px = ux*me;
        const auto py = uy*me;
        const auto pz = uz*me;

        auto e_w = 0.0_rt;
        auto p_w = 0.0_rt;
        auto e_px = 0.0_rt;
        auto e_py = 0.0_rt;
        auto e_pz = 0.0_rt;
        auto p_px = 0.0_rt;
        auto p_py = 0.0_rt;
        auto p_pz = 0.0_rt;

        //Despite the names of the variables, positrons and electrons
        //can be exchanged, since the physical process is completely
        //symmetric with respect to this exchange.
        m_generate_functor.operator()<1>(
            px, py, pz,
            ex, ey, ez,
            bx, by, bz,
            w,
            &e_px, &e_py, &e_pz,
            &p_px, &p_py, &p_pz,
            &e_w, &p_w);

        dst1.m_rdata[PIdx::w][i_dst1] = e_w;
        dst1.m_rdata[PIdx::ux][i_dst1] = e_px*one_over_me;
        dst1.m_rdata[PIdx::uy][i_dst1] = e_py*one_over_me;
        dst1.m_rdata[PIdx::uz][i_dst1] = e_pz*one_over_me;

        dst2.m_rdata[PIdx::w][i_dst2] = p_w;
        dst2.m_rdata[PIdx::ux][i_dst2] = p_px*one_over_me;
        dst2.m_rdata[PIdx::uy][i_dst2] = p_py*one_over_me;
        dst2.m_rdata[PIdx::uz][i_dst2] = p_pz*one_over_me;

        src.m_aos[i_src].id() = -1; //destroy photon after pair generation

    }

private:

    const BreitWheelerGeneratePairs
    m_generate_functor; /*!< A copy of the functor to generate pairs. It contains only pointers to the lookup tables.*/

    GetParticlePosition m_get_position;
    GetExternalEField m_get_externalE;
    GetExternalBField m_get_externalB;

    amrex::Array4<const amrex::Real> m_ex_arr;
    amrex::Array4<const amrex::Real> m_ey_arr;
    amrex::Array4<const amrex::Real> m_ez_arr;
    amrex::Array4<const amrex::Real> m_bx_arr;
    amrex::Array4<const amrex::Real> m_by_arr;
    amrex::Array4<const amrex::Real> m_bz_arr;

    amrex::IndexType m_ex_type;
    amrex::IndexType m_ey_type;
    amrex::IndexType m_ez_type;
    amrex::IndexType m_bx_type;
    amrex::IndexType m_by_type;
    amrex::IndexType m_bz_type;

    amrex::GpuArray<amrex::Real, 3> m_dx_arr;
    amrex::GpuArray<amrex::Real, 3> m_xyzmin_arr;

    int m_l_lower_order_in_v;
    int m_nox;
    int m_n_rz_azimuthal_modes;

    amrex::Dim3 m_lo;
};

#endif //QED_PAIR_GENERATION_H_
