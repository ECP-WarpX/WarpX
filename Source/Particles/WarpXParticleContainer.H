/* Copyright 2019-2020 Andrew Myers, Axel Huebl, David Grote
 * Jean-Luc Vay, Junmin Gu, Luca Fedeli
 * Maxence Thevenet, Remi Lehe, Revathi Jambunathan
 * Weiqun Zhang, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_WarpXParticleContainer_H_
#define WARPX_WarpXParticleContainer_H_

#include "WarpXParticleContainer_fwd.H"

#include "Evolve/WarpXDtType.H"
#include "Particles/ParticleBoundaries.H"
#include "SpeciesPhysicalProperties.H"
#include "Utils/TextMsg.H"

#ifdef WARPX_QED
#    include "ElementaryProcess/QEDInternals/BreitWheelerEngineWrapper_fwd.H"
#    include "ElementaryProcess/QEDInternals/QuantumSyncEngineWrapper_fwd.H"
#endif
#include "MultiParticleContainer_fwd.H"

#include <AMReX_Array.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_INT.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>
#include <AMReX_REAL.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>
#include <AMReX_AmrCoreFwd.H>

#include <array>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <utility>

namespace ParticleStringNames
{
    const std::map<std::string, int> to_index = {
        {"w",     PIdx::w    },
        {"ux",    PIdx::ux   },
        {"uy",    PIdx::uy   },
        {"uz",    PIdx::uz   },
#ifdef WARPX_DIM_RZ
        {"theta", PIdx::theta}
#endif
    };
}

class WarpXParIter
    : public amrex::ParIter<0,0,PIdx::nattribs>
{
public:
    using amrex::ParIter<0,0,PIdx::nattribs>::ParIter;

    WarpXParIter (ContainerType& pc, int level);

    WarpXParIter (ContainerType& pc, int level, amrex::MFItInfo& info);

    const std::array<RealVector, PIdx::nattribs>& GetAttribs () const {
        return GetStructOfArrays().GetRealData();
    }

    std::array<RealVector, PIdx::nattribs>& GetAttribs () {
        return GetStructOfArrays().GetRealData();
    }

    const RealVector& GetAttribs (int comp) const {
        return GetStructOfArrays().GetRealData(comp);
    }

    RealVector& GetAttribs (int comp) {
        return GetStructOfArrays().GetRealData(comp);
    }

    IntVector& GetiAttribs (int comp) {
        return GetStructOfArrays().GetIntData(comp);
    }
};

/**
 * WarpXParticleContainer is the base polymorphic class from which all concrete
 * particle container classes (that store a collection of particles) derive. Derived
 * classes can be used for plasma particles, photon particles, or non-physical
 * particles (e.g., for the laser antenna).
 * It derives from amrex::ParticleContainer<0,0,PIdx::nattribs>, where the
 * template arguments stand for the number of int and amrex::Real SoA and AoS
 * data in amrex::Particle.
 *  - AoS amrex::Real: x, y, z (default), 0 additional (first template
 *    parameter)
 *  - AoS int: id, cpu (default), 0 additional (second template parameter)
 *  - SoA amrex::Real: PIdx::nattribs (third template parameter), see PIdx for
 * the list.
 *
 * WarpXParticleContainer contains the main functions for initialization,
 * interaction with the grid (field gather and current deposition) and particle
 * push.
 *
 * Note: many functions are pure virtual (meaning they MUST be defined in
 * derived classes, e.g., Evolve) or actual functions (e.g. CurrentDeposition).
 */
class WarpXParticleContainer
    : public amrex::ParticleContainer<0,0,PIdx::nattribs>
{
public:
    friend MultiParticleContainer;

    // amrex::StructOfArrays with DiagIdx::nattribs amrex::ParticleReal components
    // and 0 int components for the particle data.
    using DiagnosticParticleData = amrex::StructOfArrays<DiagIdx::nattribs, 0>;
    // DiagnosticParticles is a vector, with one element per MR level.
    // DiagnosticParticles[lev] is typically a key-value pair where the key is
    // a pair [grid_index, tile_index], and the value is the corresponding
    // DiagnosticParticleData (see above) on this tile.
    using DiagnosticParticles = amrex::Vector<std::map<std::pair<int, int>, DiagnosticParticleData> >;

    WarpXParticleContainer (amrex::AmrCore* amr_core, int ispecies);
    virtual ~WarpXParticleContainer() {}

    virtual void InitData () = 0;

    /**
     * Evolve is the central WarpXParticleContainer function that advances
     * particles for a time dt (typically one timestep). It is a pure virtual
     * function for flexibility.
     */
    virtual void Evolve (int lev,
                         const amrex::MultiFab& Ex, const amrex::MultiFab& Ey, const amrex::MultiFab& Ez,
                         const amrex::MultiFab& Bx, const amrex::MultiFab& By, const amrex::MultiFab& Bz,
                         amrex::MultiFab& jx, amrex::MultiFab& jy, amrex::MultiFab& jz,
                         amrex::MultiFab* cjx, amrex::MultiFab* cjy, amrex::MultiFab* cjz,
                         amrex::MultiFab* rho, amrex::MultiFab* crho,
                         const amrex::MultiFab* cEx, const amrex::MultiFab* cEy, const amrex::MultiFab* cEz,
                         const amrex::MultiFab* cBx, const amrex::MultiFab* cBy, const amrex::MultiFab* cBz,
                         amrex::Real t, amrex::Real dt, DtType a_dt_type=DtType::Full, bool skip_deposition=false) = 0;

    virtual void PostRestart () = 0;

    virtual void GetParticleSlice(const int /*direction*/, const amrex::Real /*z_old*/,
                                  const amrex::Real /*z_new*/, const amrex::Real /*t_boost*/,
                                  const amrex::Real /*t_lab*/, const amrex::Real /*dt*/,
                                  DiagnosticParticles& /*diagnostic_particles*/) {}

    void AllocData ();

    ///
    /// This pushes the particle positions by one half time step.
    /// It is used to desynchronize the particles after initializaton
    /// or when restarting from a checkpoint.
    ///
    void PushX (         amrex::Real dt);
    void PushX (int lev, amrex::Real dt);

    ///
    /// This pushes the particle momenta by dt.
    ///
    virtual void PushP (int lev, amrex::Real dt,
                        const amrex::MultiFab& Ex,
                        const amrex::MultiFab& Ey,
                        const amrex::MultiFab& Ez,
                        const amrex::MultiFab& Bx,
                        const amrex::MultiFab& By,
                        const amrex::MultiFab& Bz) = 0;

    /**
     * \brief Deposit current density.
     *
     * \param[in,out] J vector of current densities (one three-dimensional array of pointers
     *                to MultiFabs per mesh refinement level)
     * \param[in] dt: Time step for particle level (is used when temporarily
     *                modifying the particle positions, either within the
     *                Esirkepov or when `relative_t` is different than 0
     * \param[in] relative_t: Time at which to deposit J, relative to the time
     *                of the current positions of the particles (expressed as
     *                a fraction of dt). When different than 0, the particle
     *                position will be temporarily modified to match the
     *                time of the deposition.
     */
    void DepositCurrent (amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& J,
                         const amrex::Real dt, const amrex::Real relative_t);

    /**
     * \brief Deposit charge density.
     *
     * \param[in,out] rho vector of charge densities (one pointer to MultiFab per mesh refinement level)
     * \param[in] local if false, exchange the data in the guard cells after the deposition
     * \param[in] reset if true, reset all values of rho to zero
     * \param[in] do_rz_volume_scaling whether to scale the final density by some volume norm in RZ geometry
     * \param[in] interpolate_across_levels whether to average down from the fine patch to the coarse patch
     * \param[in] icomp component of the MultiFab where rho is deposited (old, new)
     */
    void DepositCharge (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rho,
                        const bool local = false, const bool reset = false,
                        const bool do_rz_volume_scaling = false,
                        const bool interpolate_across_levels = true,
                        const int icomp = 0);

    std::unique_ptr<amrex::MultiFab> GetChargeDensity(int lev, bool local = false);

    virtual void DepositCharge (WarpXParIter& pti,
                               RealVector const & wp,
                               const int * const ion_lev,
                               amrex::MultiFab* rho,
                               int icomp,
                               const long offset,
                               const long np_to_depose,
                               int thread_num,
                               int lev,
                               int depos_lev);

    virtual void DepositCurrent (WarpXParIter& pti,
                                RealVector const & wp,
                                RealVector const & uxp,
                                RealVector const & uyp,
                                RealVector const & uzp,
                                int const * const ion_lev,
                                amrex::MultiFab* const jx,
                                amrex::MultiFab* const jy,
                                amrex::MultiFab* const jz,
                                long const offset,
                                long const np_to_depose,
                                int const thread_num,
                                int const lev,
                                int const depos_lev,
                                amrex::Real const dt,
                                amrex::Real const relative_time);

    // If particles start outside of the domain, ContinuousInjection
    // makes sure that they are initialized when they enter the domain, and
    // NOT before. Virtual function, overriden by derived classes.
    // Current status:
    // PhysicalParticleContainer: implemented.
    // LaserParticleContainer: implemented.
    // RigidInjectedParticleContainer: not implemented.
    virtual void ContinuousInjection(const amrex::RealBox& /*injection_box*/) {}
    // Update optional sub-class-specific injection location.
    virtual void UpdateContinuousInjectionPosition(amrex::Real /*dt*/) {}

    // Inject a continuous flux of particles from a defined plane
    virtual void ContinuousFluxInjection(amrex::Real /*t*/, amrex::Real /*dt*/) {}

    ///
    /// This returns the total charge for all the particles in this ParticleContainer.
    /// This is needed when solving Poisson's equation with periodic boundary conditions.
    ///
    amrex::Real sumParticleCharge(bool local = false);

    std::array<amrex::Real, 3> meanParticleVelocity(bool local = false);

    amrex::Real maxParticleVelocity(bool local = false);

    void AddNParticles (int lev,
                        int n, const amrex::ParticleReal* x, const amrex::ParticleReal* y, const amrex::ParticleReal* z,
                        const amrex::ParticleReal* vx, const amrex::ParticleReal* vy, const amrex::ParticleReal* vz,
                        int nattr, const amrex::ParticleReal* attr, int uniqueparticles, amrex::Long id=-1);

    virtual void ReadHeader (std::istream& is) = 0;

    virtual void WriteHeader (std::ostream& os) const = 0;

    virtual void ConvertUnits (ConvertDirection /*convert_dir*/){}

    static void ReadParameters ();

    static void BackwardCompatibility ();

    /** \brief Apply particle BC.
     *
     * \param[in] boundary_conditions Type of boundary conditions. For now, only absorbing or none
     * are supported
     */
    void ApplyBoundaryConditions ();

    bool do_splitting = false;
    bool initialize_self_fields = false;
    amrex::Real self_fields_required_precision = amrex::Real(1.e-11);
    amrex::Real self_fields_absolute_tolerance = amrex::Real(0.0);
    int self_fields_max_iters = 200;
    int self_fields_verbosity = 2;

    // split along diagonals (0) or axes (1)
    int split_type = 0;

    using amrex::ParticleContainer<0, 0, PIdx::nattribs>::AddRealComp;
    using amrex::ParticleContainer<0, 0, PIdx::nattribs>::AddIntComp;

    void AddRealComp (const std::string& name, bool comm=true)
    {
        auto search = particle_comps.find(name);
        if (search == particle_comps.end()) {
            particle_comps[name] = NumRealComps();
            particle_runtime_comps[name] = NumRealComps() - PIdx::nattribs;
            AddRealComp(comm);
        } else {
            amrex::Print() << Utils::TextMsg::Info(
                name + " already exists in particle_comps, not adding.");
        }
    }

    void AddIntComp (const std::string& name, bool comm=true)
    {
        auto search = particle_icomps.find(name);
        if (search == particle_icomps.end()) {
            particle_icomps[name] = NumIntComps();
            particle_runtime_icomps[name] = NumIntComps() - 0;
            AddIntComp(comm);
        } else {
            amrex::Print() << Utils::TextMsg::Info(
                name + " already exists in particle_icomps, not adding.");
        }
    }

    int doBackTransformedDiagnostics () const { return do_back_transformed_diagnostics; }
    /** Whether back-transformed diagnostics need to be performed for a particular species.
     *
     * \param[in] do_back_transformed_particles The parameter to set if back-transformed particles are set to true/false
     */
    void SetDoBackTransformedParticles(const bool do_back_transformed_particles) {
        m_do_back_transformed_particles = do_back_transformed_particles;
    }

    std::map<std::string, int> getParticleComps () const noexcept { return particle_comps;}
    std::map<std::string, int> getParticleiComps () const noexcept { return particle_icomps;}
    std::map<std::string, int> getParticleRuntimeComps () const noexcept { return particle_runtime_comps;}
    std::map<std::string, int> getParticleRuntimeiComps () const noexcept { return particle_runtime_icomps;}

    //amrex::Real getCharge () {return charge;}
    amrex::ParticleReal getCharge () const {return charge;}
    //amrex::Real getMass () {return mass;}
    amrex::ParticleReal getMass () const {return mass;}

    int DoFieldIonization() const { return do_field_ionization; }

#ifdef WARPX_QED
    //Species for which QED effects are relevant should override these methods
    virtual bool has_quantum_sync() const {return false;}
    virtual bool has_breit_wheeler() const {return false;}

    int DoQED() const { return has_quantum_sync() || has_breit_wheeler(); }
#else
    int DoQED() const { return false; }
#endif

    /* \brief This function tests if the current species
    *  is of a given PhysicalSpecies (specified as a template parameter).
    * @tparam PhysSpec the PhysicalSpecies to test against
    * @return the result of the test
    */
    template<PhysicalSpecies PhysSpec>
    bool AmIA () const noexcept {return (physical_species == PhysSpec);}

    /**
    * \brief This function returns a string containing the name of the species type
    */
    std::string getSpeciesTypeName () const {return species::get_name(physical_species);}

    amrex::Vector<amrex::Real> get_v_galilean () {return m_v_galilean;}

    /**
     * \brief Virtual method to resample the species. Overriden by PhysicalParticleContainer only.
     * Empty body is here because making the method purely virtual would mean that we need to
     * override the method for every derived class. Note that in practice this function is never
     * called because resample() is only called for PhysicalParticleContainers.
     */
    virtual void resample (const int /*timestep*/) {}

    /**
     * When using runtime components, AMReX requires to touch all tiles
     * in serial and create particles tiles with runtime components if
     * they do not exist (or if they were defined by default, i.e.,
     * without runtime component).
     */
     void defineAllParticleTiles () noexcept;

protected:
    amrex::Vector<amrex::Real> m_v_galilean = amrex::Vector<amrex::Real>(3, amrex::Real(0.));
    std::map<std::string, int> particle_comps;
    std::map<std::string, int> particle_icomps;
    std::map<std::string, int> particle_runtime_comps;
    std::map<std::string, int> particle_runtime_icomps;

    int species_id;

    amrex::Real charge;
    amrex::Real mass;
    PhysicalSpecies physical_species;

    // Controls boundaries for particles exiting the domain
    ParticleBoundaries m_boundary_conditions;

    //! instead of depositing (current, charge) on the finest patch level, deposit to the coarsest grid
    bool m_deposit_on_main_grid = false;

    //! instead of gathering fields from the finest patch level, gather from the coarsest
    bool m_gather_from_main_grid = false;

    int do_not_push = 0;
    int do_not_deposit = 0;
    int do_not_gather = 0;

    // Whether to allow particles outside of the simulation domain to be
    // initialized when they enter the domain.
    // This is currently required because continuous injection does not
    // support all features allowed by direct injection.
    int do_continuous_injection = 0;

    int do_field_ionization = 0;
    int ionization_product;
    std::string ionization_product_name;
    int ion_atomic_number;
    int ionization_initial_level = 0;
    amrex::Gpu::DeviceVector<amrex::Real> ionization_energies;
    amrex::Gpu::DeviceVector<amrex::Real> adk_power;
    amrex::Gpu::DeviceVector<amrex::Real> adk_prefactor;
    amrex::Gpu::DeviceVector<amrex::Real> adk_exp_prefactor;
    std::string physical_element;

    int do_resampling = 0;

    int do_back_transformed_diagnostics = 1;
    /** Whether back-transformed diagnostics is turned on for the corresponding species.*/
    bool m_do_back_transformed_particles = false;

#ifdef WARPX_QED
    //Species can receive a shared pointer to a QED engine (species for
    //which this is relevant should override these functions)
    virtual void
    set_breit_wheeler_engine_ptr(std::shared_ptr<BreitWheelerEngine>){}
    virtual void
    set_quantum_sync_engine_ptr(std::shared_ptr<QuantumSynchrotronEngine>){}

    int m_qed_breit_wheeler_ele_product;
    std::string m_qed_breit_wheeler_ele_product_name;
    int m_qed_breit_wheeler_pos_product;
    std::string m_qed_breit_wheeler_pos_product_name;
    int m_qed_quantum_sync_phot_product;
    std::string m_qed_quantum_sync_phot_product_name;

#endif
    amrex::Vector<amrex::FArrayBox> local_rho;
    amrex::Vector<amrex::FArrayBox> local_jx;
    amrex::Vector<amrex::FArrayBox> local_jy;
    amrex::Vector<amrex::FArrayBox> local_jz;

public:
    using PairIndex = std::pair<int, int>;
    using TmpParticleTile = std::array<amrex::Gpu::DeviceVector<amrex::ParticleReal>,
                                       TmpIdx::nattribs>;
    using TmpParticles = amrex::Vector<std::map<PairIndex, TmpParticleTile> >;

    TmpParticles getTmpParticleData () const noexcept {return tmp_particle_data;}
protected:
    TmpParticles tmp_particle_data;

private:
    virtual void particlePostLocate(ParticleType& p, const amrex::ParticleLocData& pld,
                                    const int lev) override;

};

#endif
