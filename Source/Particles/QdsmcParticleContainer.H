/* Copyright 2024 Marco Acciarri (Helion eEnergy Inc.)
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_QdsmcParticleContainer_H_
#define WARPX_QdsmcParticleContainer_H_

#include "Evolve/WarpXDtType.H"
#include "Evolve/WarpXPushType.H"

#include "Particles/Gather/ScaleFields.H"
#include "Particles/Resampling/Resampling.H"
#include "WarpXParticleContainer.H"

#include <AMReX_GpuContainers.H>
#include <AMReX_Particles.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <AMReX_BaseFwd.H>
#include <AMReX_AmrCoreFwd.H>

#include <memory>
#include <string>


/**
 * QdsmcParticleContainer is the ParticleContainer class containing 
 * fictitious particles needed for QDSMC algorithms used to push the
 * electron temperature used in the HybridPICModel and for neutral species
 *
 * QdsmcParticleContainer inherits from WarpXParticleContainer.
 */

class QdsmcParticleContainer
    : public WarpXParticleContainer
{
public:

    QdsmcParticleContainer (amrex::AmrCore* amr_core,
                               int ispecies,
                               const std::string& name);

    QdsmcParticleContainer (amrex::AmrCore* amr_core);

    /** This function queries deprecated input parameters and abort
     *  the run if one of them is specified. */
    void BackwardCompatibility ();

    ~QdsmcParticleContainer () override = default;

    QdsmcParticleContainer (QdsmcParticleContainer const &)              = delete;
    QdsmcParticleContainer& operator= (QdsmcParticleContainer const & )  = delete;
    QdsmcParticleContainer(QdsmcParticleContainer&& )                    = default;
    QdsmcParticleContainer& operator=(QdsmcParticleContainer&& )         = default;

    void InitData () override;

    void ReadHeader (std::istream& is) override;

    void WriteHeader (std::ostream& os) const override;

    void InitIonizationModule () override;

    /*
     * \brief Returns a pointer to the i'th plasma injector.
     */
    PlasmaInjector* GetPlasmaInjector (int i) override;

    /**
     * \brief Evolve is the central function PhysicalParticleContainer that
     * advances plasma particles for a time dt (typically one timestep).
     *
     * \param fields the WarpX field register
     * \param lev level on which particles are living
     * \param current_fp_string current coarse or fine patch identifier in fields
     * \param t current physical time
     * \param dt time step by which particles are advanced
     * \param a_dt_type type of time step (used for sub-cycling)
     * \param skip_deposition Skip the charge and current deposition.
     * \param push_type Type of particle push, explicit or implicit. Defaults to explicit
     *
     * Evolve iterates over particle iterator (each box) and performs filtering,
     * field gather, particle push and current deposition for all particles
     * in the box.
     */
    void Evolve (ablastr::fields::MultiFabRegister& fields,
                 int lev,
                 const std::string& current_fp_string,
                 amrex::Real t,
                 amrex::Real dt,
                 DtType a_dt_type=DtType::Full,
                 bool skip_deposition=false,
                 PushType push_type=PushType::Explicit) override;

    virtual void PushPX (WarpXParIter& pti,
                         amrex::FArrayBox const * exfab,
                         amrex::FArrayBox const * eyfab,
                         amrex::FArrayBox const * ezfab,
                         amrex::FArrayBox const * bxfab,
                         amrex::FArrayBox const * byfab,
                         amrex::FArrayBox const * bzfab,
                         amrex::IntVect ngEB, int /*e_is_nodal*/,
                         long offset,
                         long np_to_push,
                         int lev, int gather_lev,
                         amrex::Real dt, ScaleFields scaleFields,
                         DtType a_dt_type=DtType::Full);

    void PostRestart () final {}

    // Inject particles in Box 'part_box'
    virtual void AddParticles (int lev);

    /**
     * \brief Default initialize runtime attributes in a tile. This routine does not initialize the
     * first n_external_attr_real real attributes and the first n_external_attr_int integer
     * attributes, which have been in principle externally set elsewhere.
     *
     * @param[inout] pinned_tile the tile in which attributes are initialized
     * @param[in] n_external_attr_real The number of real attributes that have been externally set.
     * These are NOT initialized by this function.
     * @param[in] n_external_attr_int The number of integer attributes that have been externally set.
     * These are NOT initialized by this function.
     * @param[in] engine the random engine, used in initialization of QED optical depths
     */
    void DefaultInitializeRuntimeAttributes (
        typename ContainerLike<amrex::PinnedArenaAllocator>::ParticleTileType& pinned_tile,
        int n_external_attr_real,
        int n_external_attr_int) final;


    std::vector<std::string> getUserIntAttribs () const override {
        return m_user_int_attribs;
    }

    std::vector<std::string> getUserRealAttribs () const override {
        return m_user_real_attribs;
    }

    amrex::Vector< amrex::Parser* > getUserIntAttribParser () const override {
        return GetVecOfPtrs(m_user_int_attrib_parser);
    }

    amrex::Vector< amrex::Parser* > getUserRealAttribParser () const override {
        return GetVecOfPtrs(m_user_real_attrib_parser);
    }

protected:

    std::string species_name;

    /* Vector of user-defined integer attributes for species, species_name */
    std::vector<std::string> m_user_int_attribs;
    /* Vector of user-defined real attributes for species, species_name */
    std::vector<std::string> m_user_real_attribs;
    /* Vector of user-defined parser for initializing user-defined integer attributes */
    amrex::Vector< std::unique_ptr<amrex::Parser> > m_user_int_attrib_parser;
    /* Vector of user-defined parser for initializing user-defined real attributes */
    amrex::Vector< std::unique_ptr<amrex::Parser> > m_user_real_attrib_parser;

};

#endif
