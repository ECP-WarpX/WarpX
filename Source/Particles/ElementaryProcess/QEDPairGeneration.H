/* Copyright 2019 Luca Fedeli
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#ifndef QED_PAIR_GENERATION_H_
#define QED_PAIR_GENERATION_H_

#include "WarpXConst.H"
#include "WarpXParticleContainer.H"

#include <BreitWheelerEngineWrapper.H>
#include <QuantumSyncEngineWrapper.H>

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
    PairGenerationFilterFunc(int const opt_depth_runtime_comp):
        m_opt_depth_runtime_comp{opt_depth_runtime_comp}
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
        const auto opt_depth =
            ptd.m_runtime_rdata[m_opt_depth_runtime_comp][i];
        return (opt_depth < 0.0);
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
    PairGenerationTransformFunc(BreitWheelerGeneratePairs const generate_functor):
        m_generate_functor{generate_functor}
    {}

    /**
    * \brief method to set which target species have QED enabled and thus need to initialize
    * their optical depth component.
    *
    * @param[in] enable1 true if target species 1 has QED effects enabled
    * @param[in] enable2 true if target species 2 has QED effects enabled
    * @param[in] target1_opt_depth_runtime_comp index of the optical depth component of target species 1
    * @param[in] target2_opt_depth_runtime_comp index of the optical depth component of target species 2
    * @param[in] opt_depth_functor functor to initialize the optical depth of the target species
    */
    void enable_opt_depth_for_targets(const bool enable1, const bool enable2,
     const int target1_opt_depth_runtime_comp, const int target2_opt_depth_runtime_comp,
     const QuantumSynchrotronGetOpticalDepth&& opt_depth_functor)
    {
        if(enable1){
            m_target1_has_quantum_sync = true;
            m_target1_opt_depth_runtime_comp = target1_opt_depth_runtime_comp;
        }

        if(enable2){
            m_target2_has_quantum_sync = true;
            m_target2_opt_depth_runtime_comp = target2_opt_depth_runtime_comp;
        }

        if(enable1 || enable2){
            m_target_opt_depth_functor = opt_depth_functor;
        }

    }

    /**
    * \brief Functor call. It determines the properties of the generated pair.
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
    void operator() (DstData& dst1, DstData& dst2, const SrcData& src,
        int const i_src, int const i_dst1, int const i_dst2) const noexcept
    {

        constexpr amrex::ParticleReal me = PhysConst::m_e;
        const amrex::ParticleReal one_over_me = 1./me;

        const amrex::ParticleReal w = src.m_rdata[PIdx::w][i_src];
        const amrex::ParticleReal ux = src.m_rdata[PIdx::ux][i_src];
        const amrex::ParticleReal uy = src.m_rdata[PIdx::uy][i_src];
        const amrex::ParticleReal uz = src.m_rdata[PIdx::uz][i_src];
        const amrex::ParticleReal ex = src.m_rdata[PIdx::Ex][i_src];
        const amrex::ParticleReal ey = src.m_rdata[PIdx::Ey][i_src];
        const amrex::ParticleReal ez = src.m_rdata[PIdx::Ez][i_src];
        const amrex::ParticleReal bx = src.m_rdata[PIdx::Bx][i_src];
        const amrex::ParticleReal by = src.m_rdata[PIdx::By][i_src];
        const amrex::ParticleReal bz = src.m_rdata[PIdx::Bz][i_src];

        const auto px = ux*me;
        const auto py = uy*me;
        const auto pz = uz*me;

        auto e_w = amrex::ParticleReal{0.0};
        auto p_w = amrex::ParticleReal{0.0};
        auto e_px = amrex::ParticleReal{0.0};
        auto e_py = amrex::ParticleReal{0.0};
        auto e_pz = amrex::ParticleReal{0.0};
        auto p_px = amrex::ParticleReal{0.0};
        auto p_py = amrex::ParticleReal{0.0};
        auto p_pz = amrex::ParticleReal{0.0};

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

        src.m_rdata[PIdx::ux][i_src] = px*one_over_me;
        src.m_rdata[PIdx::uy][i_src] = py*one_over_me;
        src.m_rdata[PIdx::uz][i_src] = pz*one_over_me;

        dst1.m_rdata[PIdx::w][i_dst1] = e_w;
        dst1.m_rdata[PIdx::ux][i_dst1] = e_px*one_over_me;
        dst1.m_rdata[PIdx::uy][i_dst1] = e_py*one_over_me;
        dst1.m_rdata[PIdx::uz][i_dst1] = e_pz*one_over_me;

        dst2.m_rdata[PIdx::w][i_dst2] = p_w;
        dst2.m_rdata[PIdx::ux][i_dst2] = p_px*one_over_me;
        dst2.m_rdata[PIdx::uy][i_dst2] = p_py*one_over_me;
        dst2.m_rdata[PIdx::uz][i_dst2] = p_pz*one_over_me;

        //Initialize the optical depth component of the target species if needed.
        if(m_target1_has_quantum_sync){
            dst1.m_runtime_rdata[m_target1_opt_depth_runtime_comp][i_dst1] =
                m_target_opt_depth_functor();
        }

        if(m_target2_has_quantum_sync){
            dst2.m_runtime_rdata[m_target2_opt_depth_runtime_comp][i_dst2] =
                m_target_opt_depth_functor();
        }

    }

private:

    const BreitWheelerGeneratePairs
        m_generate_functor; /*!< A copy of the functor to generate pairs. It contains only pointers to the lookup tables.*/

    bool m_target1_has_quantum_sync = false; /*!< Flag to tell if target species 1 has QED effects enabled*/
    int m_target1_opt_depth_runtime_comp = 0; /*!< Index of the optical depth component of target species 1*/
    bool m_target2_has_quantum_sync = false; /*!< Flag to tell if target species 2 has QED effects enabled*/
    int m_target2_opt_depth_runtime_comp = 0; /*!< Index of the optical depth component of target species 2*/
    QuantumSynchrotronGetOpticalDepth
        m_target_opt_depth_functor; /*!< A copy of the functor to initialize the optical depth of the target species.*/


};

#endif //QED_PAIR_GENERATION_H_
