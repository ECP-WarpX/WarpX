#ifndef WARPX_FLUSHFORMATCATALYST_H_
#define WARPX_FLUSHFORMATCATALYST_H_

#include "FlushFormatInSitu.H"

#include "Diagnostics/ParticleDiag/ParticleDiag_fwd.H"

#ifdef AMREX_USE_CONDUIT
#   include <AMReX_Conduit_Blueprint.H>
#endif
#include <AMReX_Geometry.H>
#include <AMReX_Vector.H>

#include <AMReX_BaseFwd.H>

#ifdef AMREX_USE_CATALYST
#   include <catalyst.hpp>
#endif

#include <string>

/**
 * \brief This class aims at dumping performing in-situ analysis and visualization
 * with Catalyst. Catalyst initialize and finalize are called in the constructor 
 * and destructor respectivelyThe function WriteToFile takes in the particles, 
 * meshes, and fields, and formats them to be used by given pipeline scripts. All 
 * exports are defined and executed externally through the given catalyst pipeline 
 * and implementation. 
 */
class FlushFormatCatalyst : public FlushFormatInSitu
{
public:
    // Initailize
    FlushFormatCatalyst();

    /** Send particle, mesh, and field information through the catalyst pipeline */
    virtual void WriteToFile (
        const amrex::Vector<std::string> varnames,
        const amrex::Vector<amrex::MultiFab>& mf,
        amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<int> iteration, const double time,
        const amrex::Vector<ParticleDiag>& particle_diags, int nlev,
        const std::string prefix, int file_min_digits,
        bool plot_raw_fields,
        bool plot_raw_fields_guards,
        const bool use_pinned_pc = false,
        bool isBTD = false, int snapshotID = -1,
        int bufferID = 1, int numBuffers = 1,
        const amrex::Geometry& full_BTD_snapshot = amrex::Geometry(),
        bool isLastBTDFlush = false,
        const amrex::Vector<int>& totalParticlesFlushedAlready = amrex::Vector<int>() ) const override;

    // Finalize
    ~FlushFormatCatalyst();
};

#endif