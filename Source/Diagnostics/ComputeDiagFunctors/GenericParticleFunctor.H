/* Copyright 2023 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_GENERICPARTICLEFUNCTOR_H_
#define WARPX_GENERICPARTICLEFUNCTOR_H_

#include "ComputeParticleDiagFunctor.H"
#include "Particles/Pusher/GetAndSetPosition.H"

#include <AMReX.H>
#include <AMReX_AmrParticles.H>




/**
 * \brief Filter to select particles that correspond to a z-slice of the corresponding lab-frame.
 */
struct SelectParticles
{
    using TmpParticles = WarpXParticleContainer::TmpParticles;

    /**
     * \brief Constructor of SelectParticles functor.
     *
     * @param[in] a_pti              WarpX particle iterator
     * @param[in] tmp_particle_data  temporary particle data
     * @param[in] current_z_boost    current z-position of the slice in boosted frame
     * @param[in] old_z_boost        previous z-position of the slice in boosted frame
     * @param[in] a_offset           index offset for particles to be selected
     */
    SelectParticles( const WarpXParIter& a_pti, TmpParticles& tmp_particle_data,
                     amrex::Real current_z_boost, amrex::Real old_z_boost,
                     int a_offset = 0);

    /**
     * \brief Functor call. This method determines if a given particle should be selected
     * for Lorentz transformation in obtaining the lab-frame data. The particles that
     * with positions that correspond to the specific z-slice in boosted frame are selected.
     *
     * @tparam SrcData type of source data
     * @param[in] src SrcData particle tile data
     * @param[in] i       particle index
     * @return 1 if particles is selected for transformation, else 0
     */
    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int i) const noexcept
    {
        amrex::ignore_unused(src);
        amrex::ParticleReal xp, yp, zp;
        m_get_position(i, xp, yp, zp);
        int Flag = 0;
        if ( ( (zp >= m_current_z_boost) && (zpold[i] <= m_old_z_boost) ) ||
             ( (zp <= m_current_z_boost) && (zpold[i] >= m_old_z_boost) ))
        {    Flag = 1;
        }
        return Flag;
    }

    /** Object to extract the positions of the macroparticles inside a ParallelFor kernel */
    GetParticlePosition<PIdx> m_get_position;
    /** Z coordinate in boosted frame that corresponds to a give snapshot*/
    amrex::Real m_current_z_boost;
    /** Previous Z coordinate in boosted frame that corresponds to a give snapshot*/
    amrex::Real m_old_z_boost;
    /** Particle z coordinate in boosted frame*/
    amrex::ParticleReal* AMREX_RESTRICT zpold = nullptr;
};

/**
 * \brief Functor to select particles and store them in particle buffers
 */
class
GenericParticleFunctor final : public ComputeParticleDiagFunctor
{
public:
    GenericParticleFunctor(WarpXParticleContainer *pc_src, std::string species_name, int num_buffers);
    /** Computes the Lorentz transform of source particles to obtain lab-frame data in pc_dst*/
    void operator () (PinnedMemoryParticleContainer& pc_dst, int &TotalParticleCounter, int i_buffer) const override;
    void InitData() override;
    /** \brief Prepare data required to back-transform particle attribtutes for lab-frame snapshot, i_buffer
     *
     * \param[in] i_buffer           index of the snapshot
     * \param[in] z_slice_in_domain  if the z-slice at current_z_boost is within the bounds of
     *            the boosted-frame and lab-frame domain. The particles are transformed
     *            only if this value is true.
     * \param[in] old_z_boost        previous z-position of the slice in boosted frame
     * \param[in] current_z_boost    z co-ordinate of the slice selected in boosted-frame.
     * \param[in] t_lab              current time in lab-frame for snapshot, i_buffer.
     * \param[in] snapshot_full      if the current snapshot, with index, i_buffer, is
     *            full (1) or not (0). If it is full, then Lorentz-transform is not performed
     *            by setting m_perform_backtransform to 0 for the corresponding ith snapshot.
     */
    void PrepareFunctorData ( int i_buffer, bool z_slice_in_domain, amrex::Real old_z_boost,
                              amrex::Real current_z_boost, amrex::Real t_lab,
                              int snapshot_full) override;
private:
    /** Source particle container */
    WarpXParticleContainer* m_pc_src = nullptr;
    /** String containing species name of particles being written */
    std::string m_species_name;
};

#endif //
