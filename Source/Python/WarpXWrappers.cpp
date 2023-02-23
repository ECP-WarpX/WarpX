/* Copyright 2019 Andrew Myers, Axel Huebl, David Grote
 * Luca Fedeli, Maxence Thevenet, Remi Lehe
 * Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BoundaryConditions/PML.H"
#include "Initialization/WarpXAMReXInit.H"
#include "Particles/MultiParticleContainer.H"
#include "Particles/ParticleBoundaryBuffer.H"
#include "Particles/WarpXParticleContainer.H"
#include "Utils/WarpXProfilerWrapper.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"
#include "WarpXWrappers.H"
#include "WarpX_py.H"

#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_ArrayOfStructs.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuControl.H>
#include <AMReX_IndexType.H>
#include <AMReX_IntVect.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PODVector.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>
#include <AMReX_StructOfArrays.H>

#include <array>
#include <cstdlib>

namespace
{
    amrex::Real** getMultiFabPointers (amrex::MultiFab& mf, int *num_boxes, int *ncomps, int **ngrowvect, int **shapes)
    {
        *ncomps = mf.nComp();
        *num_boxes = mf.local_size();
        int shapesize = AMREX_SPACEDIM;
        *ngrowvect = static_cast<int*>(malloc(sizeof(int)*shapesize));
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            (*ngrowvect)[j] = mf.nGrow(j);
        }
        if (mf.nComp() > 1) shapesize += 1;
        *shapes = static_cast<int*>(malloc(sizeof(int)*shapesize * (*num_boxes)));
        auto data =
            static_cast<amrex::Real**>(malloc((*num_boxes) * sizeof(amrex::Real*)));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi ) {
            int i = mfi.LocalIndex();
            data[i] = mf[mfi].dataPtr();
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                (*shapes)[shapesize*i+j] = mf[mfi].box().length(j);
            }
            if (mf.nComp() > 1) (*shapes)[shapesize*i+AMREX_SPACEDIM] = mf.nComp();
        }
        return data;
    }
    int* getMultiFabLoVects (const amrex::MultiFab& mf, int *num_boxes, int **ngrowvect)
    {
        int shapesize = AMREX_SPACEDIM;
        *ngrowvect = static_cast<int*>(malloc(sizeof(int)*shapesize));
        for (int j = 0; j < AMREX_SPACEDIM; ++j) {
            (*ngrowvect)[j] = mf.nGrow(j);
        }
        *num_boxes = mf.local_size();
        auto loVects = static_cast<int*>(malloc((*num_boxes)*AMREX_SPACEDIM * sizeof(int)));

        int i = 0;
        for ( amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi, ++i ) {
            const int* loVect = mf[mfi].loVect();
            for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                loVects[AMREX_SPACEDIM*i+j] = loVect[j];
            }
        }
        return loVects;
    }
    // Copy the nodal flag data and return the copy:
    // the nodal flag data should not be modifiable from Python.
    int* getFieldNodalFlagData ( const amrex::MultiFab* mf )
    {
        if (mf == nullptr) return nullptr;
        const amrex::IntVect nodal_flag( mf->ixType().toIntVect() );
        auto *nodal_flag_data = static_cast<int*>(malloc(AMREX_SPACEDIM * sizeof(int)));

        constexpr int NODE = amrex::IndexType::NODE;

        for (int i=0 ; i < AMREX_SPACEDIM ; i++) {
            nodal_flag_data[i] = (nodal_flag[i] == NODE ? 1 : 0);
        }
        return nodal_flag_data;
    }
}

    int warpx_Real_size()
    {
        return (int)sizeof(amrex::Real);
    }

    int warpx_ParticleReal_size()
    {
        return (int)sizeof(amrex::ParticleReal);
    }

    int warpx_nSpecies()
    {
        const auto & mypc = WarpX::GetInstance().GetPartContainer();
        return mypc.nSpecies();
    }

    bool warpx_use_fdtd_nci_corr()
    {
        return WarpX::use_fdtd_nci_corr;
    }

    int warpx_galerkin_interpolation()
    {
        return WarpX::galerkin_interpolation;
    }

    int warpx_nComps()
    {
        return PIdx::nattribs;
    }

    int warpx_nCompsSpecies(const char* char_species_name)
    {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);
        return myspc.NumRealComps();
    }

    int warpx_SpaceDim()
    {
        return AMREX_SPACEDIM;
    }

    void amrex_init (int argc, char* argv[])
    {
        warpx_amrex_init(argc, argv);
    }

    void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm mpicomm)
    {
        warpx_amrex_init(argc, argv, true, mpicomm);
    }

    void amrex_finalize (int /*finalize_mpi*/)
    {
        amrex::Finalize();
    }

    void warpx_init ()
    {
        WarpX& warpx = WarpX::GetInstance();
        warpx.InitData();
        ExecutePythonCallback("afterinit");
        ExecutePythonCallback("particleloader");
    }

    void warpx_finalize ()
    {
        WarpX::ResetInstance();
    }

    void warpx_set_callback_py (
        const char* char_callback_name, WARPX_CALLBACK_PY_FUNC_0 callback)
    {
        const std::string callback_name(char_callback_name);
        warpx_callback_py_map[callback_name] = callback;
    }

    void warpx_clear_callback_py (const char* char_callback_name)
    {
        const std::string callback_name(char_callback_name);
        warpx_callback_py_map.erase(callback_name);
    }

    void warpx_evolve (int numsteps)
    {
        WarpX& warpx = WarpX::GetInstance();
        warpx.Evolve(numsteps);
    }

    void warpx_addNParticles(
        const char* char_species_name, int lenx, amrex::ParticleReal const * x,
        amrex::ParticleReal const * y, amrex::ParticleReal const * z,
        amrex::ParticleReal const * vx, amrex::ParticleReal const * vy,
        amrex::ParticleReal const * vz, const int nattr_real,
        amrex::ParticleReal const * attr_real, const int nattr_int,
        int const * attr_int, int uniqueparticles)
    {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);
        const int lev = 0;
        myspc.AddNParticles(lev, lenx, x, y, z, vx, vy, vz, nattr_real, attr_real,
                            nattr_int, attr_int, uniqueparticles);
    }

    void warpx_ConvertLabParamsToBoost()
    {
      ConvertLabParamsToBoost();
    }

    void warpx_ReadBCParams()
    {
      ReadBCParams();
    }

    void warpx_CheckGriddingForRZSpectral()
    {
        CheckGriddingForRZSpectral();
    }

    amrex::Real warpx_getProbLo(int dir)
    {
      WarpX& warpx = WarpX::GetInstance();
      const amrex::Geometry& geom = warpx.Geom(0);
      return geom.ProbLo(dir);
    }

    amrex::Real warpx_getProbHi(int dir)
    {
      WarpX& warpx = WarpX::GetInstance();
      const amrex::Geometry& geom = warpx.Geom(0);
      return geom.ProbHi(dir);
    }

    amrex::Real warpx_getCellSize(int dir, int lev) {
        const std::array<amrex::Real,3>& dx = WarpX::CellSize(lev);
        return dx[dir];
    }

    long warpx_getNumParticles(const char* char_species_name, const bool local) {
        const auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);
        // the first argument below is to only count valid particles
        return myspc.TotalNumberOfParticles(true, local);
    }

#define WARPX_GET_FIELD(FIELD, GETTER) \
    amrex::Real** FIELD(int lev, int direction, \
                        int *return_size, int *ncomps, int **ngrowvect, int **shapes) { \
        auto * mf = GETTER(lev, direction); \
        if (mf != nullptr) { \
            return getMultiFabPointers(*mf, return_size, ncomps, ngrowvect, shapes); \
        } else { \
            return nullptr; \
        } \
    }

#define WARPX_GET_LOVECTS(FIELD, GETTER) \
    int* FIELD(int lev, int direction, \
               int *return_size, int **ngrowvect) { \
        auto * mf = GETTER(lev, direction); \
        if (mf != nullptr) { \
            return getMultiFabLoVects(*mf, return_size, ngrowvect); \
        } else { \
            return nullptr; \
        } \
    }

    WARPX_GET_FIELD(warpx_getEfield, WarpX::GetInstance().get_pointer_Efield_aux)
    WARPX_GET_FIELD(warpx_getEfieldCP, WarpX::GetInstance().get_pointer_Efield_cp)
    WARPX_GET_FIELD(warpx_getEfieldFP, WarpX::GetInstance().get_pointer_Efield_fp)

    WARPX_GET_FIELD(warpx_getBfield, WarpX::GetInstance().get_pointer_Bfield_aux)
    WARPX_GET_FIELD(warpx_getBfieldCP, WarpX::GetInstance().get_pointer_Bfield_cp)
    WARPX_GET_FIELD(warpx_getBfieldFP, WarpX::GetInstance().get_pointer_Bfield_fp)

    WARPX_GET_FIELD(warpx_getEdgeLengths, WarpX::GetInstance().get_pointer_edge_lengths)
    WARPX_GET_FIELD(warpx_getFaceAreas, WarpX::GetInstance().get_pointer_face_areas)

    WARPX_GET_FIELD(warpx_getCurrentDensity, WarpX::GetInstance().get_pointer_current_fp)
    WARPX_GET_FIELD(warpx_getCurrentDensityCP, WarpX::GetInstance().get_pointer_current_cp)
    WARPX_GET_FIELD(warpx_getCurrentDensityFP, WarpX::GetInstance().get_pointer_current_fp)

    WARPX_GET_FIELD(warpx_getVectorPotentialFP, WarpX::GetInstance().get_pointer_vector_potential_fp)

    WARPX_GET_LOVECTS(warpx_getEfieldLoVects, WarpX::GetInstance().get_pointer_Efield_aux)
    WARPX_GET_LOVECTS(warpx_getEfieldCPLoVects, WarpX::GetInstance().get_pointer_Efield_cp)
    WARPX_GET_LOVECTS(warpx_getEfieldFPLoVects, WarpX::GetInstance().get_pointer_Efield_fp)

    WARPX_GET_LOVECTS(warpx_getBfieldLoVects, WarpX::GetInstance().get_pointer_Bfield_aux)
    WARPX_GET_LOVECTS(warpx_getBfieldCPLoVects, WarpX::GetInstance().get_pointer_Bfield_cp)
    WARPX_GET_LOVECTS(warpx_getBfieldFPLoVects, WarpX::GetInstance().get_pointer_Bfield_fp)

    WARPX_GET_LOVECTS(warpx_getCurrentDensityLoVects, WarpX::GetInstance().get_pointer_current_fp)
    WARPX_GET_LOVECTS(warpx_getCurrentDensityCPLoVects, WarpX::GetInstance().get_pointer_current_cp)
    WARPX_GET_LOVECTS(warpx_getCurrentDensityFPLoVects, WarpX::GetInstance().get_pointer_current_fp)
    WARPX_GET_LOVECTS(warpx_getVectorPotentialFPLoVects, WarpX::GetInstance().get_pointer_vector_potential_fp)

    WARPX_GET_LOVECTS(warpx_getEdgeLengthsLoVects, WarpX::GetInstance().get_pointer_edge_lengths)
    WARPX_GET_LOVECTS(warpx_getFaceAreasLoVects, WarpX::GetInstance().get_pointer_face_areas)

    int* warpx_getEx_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_Efield_aux(0,0) );}
    int* warpx_getEy_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_Efield_aux(0,1) );}
    int* warpx_getEz_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_Efield_aux(0,2) );}
    int* warpx_getBx_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_Bfield_aux(0,0) );}
    int* warpx_getBy_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_Bfield_aux(0,1) );}
    int* warpx_getBz_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_Bfield_aux(0,2) );}
    int* warpx_getJx_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_current_fp(0,0) );}
    int* warpx_getJy_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_current_fp(0,1) );}
    int* warpx_getJz_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_current_fp(0,2) );}
    int* warpx_getAx_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_vector_potential_fp(0,0) );}
    int* warpx_getAy_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_vector_potential_fp(0,1) );}
    int* warpx_getAz_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_vector_potential_fp(0,2) );}
    int* warpx_getRho_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_rho_fp(0) );}
    int* warpx_getPhi_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_phi_fp(0) );}
    int* warpx_getF_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_F_fp(0) );}
    int* warpx_getG_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_G_fp(0) );}
    int* warpx_get_edge_lengths_x_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_edge_lengths(0, 0) );}
    int* warpx_get_edge_lengths_y_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_edge_lengths(0, 1) );}
    int* warpx_get_edge_lengths_z_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_edge_lengths(0, 2) );}
    int* warpx_get_face_areas_x_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_face_areas(0, 0) );}
    int* warpx_get_face_areas_y_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_face_areas(0, 1) );}
    int* warpx_get_face_areas_z_nodal_flag() {return getFieldNodalFlagData( WarpX::GetInstance().get_pointer_face_areas(0, 2) );}

#define WARPX_GET_SCALAR(SCALAR, GETTER) \
    amrex::Real** SCALAR(int lev, \
                         int *return_size, int *ncomps, int **ngrowvect, int **shapes) { \
        auto * mf = GETTER(lev); \
        if (mf != nullptr) { \
            return getMultiFabPointers(*mf, return_size, ncomps, ngrowvect, shapes); \
        } else { \
            return nullptr; \
        } \
    }

#define WARPX_GET_LOVECTS_SCALAR(SCALAR, GETTER) \
    int* SCALAR(int lev, \
                int *return_size, int **ngrowvect) { \
        auto * mf = GETTER(lev); \
        if (mf != nullptr) { \
            return getMultiFabLoVects(*mf, return_size, ngrowvect); \
        } else { \
            return nullptr; \
        } \
    }

    WARPX_GET_SCALAR(warpx_getChargeDensityCP, WarpX::GetInstance().get_pointer_rho_cp)
    WARPX_GET_SCALAR(warpx_getChargeDensityFP, WarpX::GetInstance().get_pointer_rho_fp)

    WARPX_GET_LOVECTS_SCALAR(warpx_getChargeDensityCPLoVects, WarpX::GetInstance().get_pointer_rho_cp)
    WARPX_GET_LOVECTS_SCALAR(warpx_getChargeDensityFPLoVects, WarpX::GetInstance().get_pointer_rho_fp)

    WARPX_GET_SCALAR(warpx_getPhiFP, WarpX::GetInstance().get_pointer_phi_fp)

    WARPX_GET_LOVECTS_SCALAR(warpx_getPhiFPLoVects, WarpX::GetInstance().get_pointer_phi_fp)

    // F and G
    WARPX_GET_SCALAR(warpx_getFfieldCP, WarpX::GetInstance().get_pointer_F_cp)
    WARPX_GET_SCALAR(warpx_getFfieldFP, WarpX::GetInstance().get_pointer_F_fp)
    WARPX_GET_LOVECTS_SCALAR(warpx_getFfieldCPLoVects, WarpX::GetInstance().get_pointer_F_cp)
    WARPX_GET_LOVECTS_SCALAR(warpx_getFfieldFPLoVects, WarpX::GetInstance().get_pointer_F_fp)
    WARPX_GET_SCALAR(warpx_getGfieldCP, WarpX::GetInstance().get_pointer_G_cp)
    WARPX_GET_SCALAR(warpx_getGfieldFP, WarpX::GetInstance().get_pointer_G_fp)
    WARPX_GET_LOVECTS_SCALAR(warpx_getGfieldCPLoVects, WarpX::GetInstance().get_pointer_G_cp)
    WARPX_GET_LOVECTS_SCALAR(warpx_getGfieldFPLoVects, WarpX::GetInstance().get_pointer_G_fp)

#define WARPX_GET_FIELD_PML(FIELD, GETTER) \
    amrex::Real** FIELD(int lev, int direction, \
                        int *return_size, int *ncomps, int **ngrowvect, int **shapes) { \
        auto * pml = WarpX::GetInstance().GetPML(lev); \
        if (!pml) return nullptr; \
        auto * mf = (pml->GETTER()[direction]); \
        if (!mf) return nullptr; \
        return getMultiFabPointers(*mf, return_size, ncomps, ngrowvect, shapes); \
    }

#define WARPX_GET_LOVECTS_PML(FIELD, GETTER) \
    int* FIELD(int lev, int direction, \
               int *return_size, int **ngrowvect) { \
        auto * pml = WarpX::GetInstance().GetPML(lev); \
        if (!pml) return nullptr; \
        auto * mf = (pml->GETTER()[direction]); \
        if (!mf) return nullptr; \
        return getMultiFabLoVects(*mf, return_size, ngrowvect); \
    }

    WARPX_GET_FIELD_PML(warpx_getEfieldCP_PML, GetE_cp)
    WARPX_GET_FIELD_PML(warpx_getEfieldFP_PML, GetE_fp)
    WARPX_GET_FIELD_PML(warpx_getBfieldCP_PML, GetB_cp)
    WARPX_GET_FIELD_PML(warpx_getBfieldFP_PML, GetB_fp)
    WARPX_GET_FIELD_PML(warpx_getCurrentDensityCP_PML, Getj_cp)
    WARPX_GET_FIELD_PML(warpx_getCurrentDensityFP_PML, Getj_fp)
    WARPX_GET_LOVECTS_PML(warpx_getEfieldCPLoVects_PML, GetE_cp)
    WARPX_GET_LOVECTS_PML(warpx_getEfieldFPLoVects_PML, GetE_fp)
    WARPX_GET_LOVECTS_PML(warpx_getBfieldCPLoVects_PML, GetB_cp)
    WARPX_GET_LOVECTS_PML(warpx_getBfieldFPLoVects_PML, GetB_fp)
    WARPX_GET_LOVECTS_PML(warpx_getCurrentDensityCPLoVects_PML, Getj_cp)
    WARPX_GET_LOVECTS_PML(warpx_getCurrentDensityFPLoVects_PML, Getj_fp)

#define WARPX_GET_SCALAR_PML(SCALAR, GETTER) \
    amrex::Real** SCALAR(int lev, \
                        int *return_size, int *ncomps, int **ngrowvect, int **shapes) { \
        auto * pml = WarpX::GetInstance().GetPML(lev); \
        if (!pml) return nullptr; \
        auto * mf = pml->GETTER(); \
        if (!mf) return nullptr; \
        return getMultiFabPointers(*mf, return_size, ncomps, ngrowvect, shapes); \
    }

#define WARPX_GET_LOVECTS_PML_SCALAR(SCALAR, GETTER) \
    int* SCALAR(int lev, \
               int *return_size, int **ngrowvect) { \
        auto * pml = WarpX::GetInstance().GetPML(lev); \
        if (!pml) return nullptr; \
        auto * mf = pml->GETTER(); \
        if (!mf) return nullptr; \
        return getMultiFabLoVects(*mf, return_size, ngrowvect); \
    }

    // F and G
    WARPX_GET_SCALAR_PML(warpx_getFfieldCP_PML, GetF_cp)
    WARPX_GET_SCALAR_PML(warpx_getFfieldFP_PML, GetF_fp)
    WARPX_GET_LOVECTS_PML_SCALAR(warpx_getFfieldCPLoVects_PML, GetF_cp)
    WARPX_GET_LOVECTS_PML_SCALAR(warpx_getFfieldFPLoVects_PML, GetF_fp)
    WARPX_GET_SCALAR_PML(warpx_getGfieldCP_PML, GetG_cp)
    WARPX_GET_SCALAR_PML(warpx_getGfieldFP_PML, GetG_fp)
    WARPX_GET_LOVECTS_PML_SCALAR(warpx_getGfieldCPLoVects_PML, GetG_cp)
    WARPX_GET_LOVECTS_PML_SCALAR(warpx_getGfieldFPLoVects_PML, GetG_fp)

    int* warpx_getF_pml_nodal_flag ()
    {
        auto * pml = WarpX::GetInstance().GetPML(0);
        if (!pml) return nullptr;
        return getFieldNodalFlagData(pml->GetF_fp());
    }

    int* warpx_getG_pml_nodal_flag ()
    {
        auto * pml = WarpX::GetInstance().GetPML(0);
        if (!pml) return nullptr;
        return getFieldNodalFlagData(pml->GetG_fp());
    }

    amrex::ParticleReal** warpx_getParticleStructs(
            const char* char_species_name, int lev,
            int* num_tiles, int** particles_per_tile) {
        const auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);

        *num_tiles = myspc.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<amrex::ParticleReal**>(malloc(*num_tiles*sizeof(typename WarpXParticleContainer::ParticleType*)));
        int i = 0;
        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti, ++i) {
            auto& aos = pti.GetArrayOfStructs();
            data[i] = (amrex::ParticleReal*) aos.data();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    amrex::ParticleReal** warpx_getParticleArrays (
            const char* char_species_name, const char* char_comp_name,
            int lev, int* num_tiles, int** particles_per_tile ) {

        const auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);

        int comp = warpx_getParticleCompIndex(char_species_name, char_comp_name);

        *num_tiles = myspc.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<amrex::ParticleReal**>(malloc(*num_tiles*sizeof(amrex::ParticleReal*)));
        int i = 0;
        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti, ++i) {
            auto& soa = pti.GetStructOfArrays();
            data[i] = (amrex::ParticleReal*) soa.GetRealData(comp).dataPtr();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    void warpx_convert_id_to_long (amrex::Long* ids, const WarpXParticleContainer::ParticleType* pstructs, int size)
    {
        amrex::Long* d_ptr = nullptr;
#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<amrex::Long> d_ids(size);
        d_ptr = d_ids.data();
#else
        d_ptr = ids;
#endif
        amrex::ParallelFor(size, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            d_ptr[i] = pstructs[i].id();
        });
#ifdef AMREX_USE_GPU
        amrex::Gpu::dtoh_memcpy(ids, d_ptr, size*sizeof(amrex::Long));
#endif
    }

    void warpx_convert_cpu_to_int (int* cpus, const WarpXParticleContainer::ParticleType* pstructs, int size)
    {
        int* d_ptr = nullptr;
#ifdef AMREX_USE_GPU
        amrex::Gpu::DeviceVector<int> d_cpus(size);
        d_ptr = d_cpus.data();
#else
        d_ptr = cpus;
#endif
        amrex::ParallelFor(size, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            d_ptr[i] = pstructs[i].cpu();
        });
#ifdef AMREX_USE_GPU
        amrex::Gpu::dtoh_memcpy(cpus, d_ptr, size*sizeof(int));
#endif
    }

    int warpx_getParticleCompIndex (
         const char* char_species_name, const char* char_comp_name )
    {
        const auto & mypc = WarpX::GetInstance().GetPartContainer();

        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);

        const std::string comp_name(char_comp_name);
        auto particle_comps = myspc.getParticleComps();

        return particle_comps.at(comp_name);
    }

    void warpx_addRealComp(const char* char_species_name,
        const char* char_comp_name, bool comm=true)
    {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);

        const std::string comp_name(char_comp_name);
        myspc.AddRealComp(comp_name, comm);

        mypc.defineAllParticleTiles();
    }

    amrex::Real warpx_sumParticleCharge(const char* char_species_name, const bool local)
    {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);
        return myspc.sumParticleCharge(local);
    }

    int warpx_getParticleBoundaryBufferSize(const char* species_name, int boundary, bool local)
    {
        const std::string name(species_name);
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        return particle_buffers.getNumParticlesInContainer(species_name, boundary, local);
    }

    int** warpx_getParticleBoundaryBufferScrapedSteps(const char* species_name, int boundary, int lev,
                     int* num_tiles, int** particles_per_tile)
    {
        const std::string name(species_name);
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        auto& particle_buffer = particle_buffers.getParticleBuffer(species_name, boundary);

        const int comp = particle_buffer.NumIntComps() - 1;

        *num_tiles = particle_buffer.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<int**>(malloc(*num_tiles*sizeof(int*)));
        int i = 0;
        for (amrex::ParIter<0,0,PIdx::nattribs, 0, amrex::PinnedArenaAllocator> pti(particle_buffer, lev); pti.isValid(); ++pti, ++i) {
            auto& soa = pti.GetStructOfArrays();
            data[i] = (int*) soa.GetIntData(comp).dataPtr();
            (*particles_per_tile)[i] = pti.numParticles();
        }

        return data;
    }

    amrex::ParticleReal** warpx_getParticleBoundaryBuffer(const char* species_name, int boundary, int lev,
                     int* num_tiles, int** particles_per_tile, const char* comp_name)
    {
        const std::string name(species_name);
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        auto& particle_buffer = particle_buffers.getParticleBuffer(species_name, boundary);

        const int comp = warpx_getParticleCompIndex(species_name, comp_name);

        *num_tiles = particle_buffer.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<amrex::ParticleReal**>(malloc(*num_tiles*sizeof(amrex::ParticleReal*)));
        int i = 0;
        for (amrex::ParIter<0,0,PIdx::nattribs, 0, amrex::PinnedArenaAllocator> pti(particle_buffer, lev); pti.isValid(); ++pti, ++i) {
            auto& soa = pti.GetStructOfArrays();
            data[i] = (amrex::ParticleReal*) soa.GetRealData(comp).dataPtr();
            (*particles_per_tile)[i] = pti.numParticles();
        }

        return data;
    }

    amrex::ParticleReal** warpx_getParticleBoundaryBufferStructs(const char* species_name, int boundary, int lev,
                     int* num_tiles, int** particles_per_tile)
    {
        const std::string name(species_name);
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        auto& particle_buffer = particle_buffers.getParticleBuffer(species_name, boundary);

        *num_tiles = particle_buffer.numLocalTilesAtLevel(lev);
        *particles_per_tile = static_cast<int*>(malloc(*num_tiles*sizeof(int)));
        memset(*particles_per_tile, 0, *num_tiles*sizeof(int));

        auto data = static_cast<amrex::ParticleReal**>(malloc(*num_tiles*sizeof(typename WarpXParticleContainer::ParticleType*)));
        int i = 0;
        for (amrex::ParIter<0,0,PIdx::nattribs, 0, amrex::PinnedArenaAllocator> pti(particle_buffer, lev); pti.isValid(); ++pti, ++i) {
            auto& aos = pti.GetArrayOfStructs();
            data[i] = (amrex::ParticleReal*) aos.data();
            (*particles_per_tile)[i] = pti.numParticles();
        }
        return data;
    }

    void warpx_clearParticleBoundaryBuffer () {
        auto& particle_buffers = WarpX::GetInstance().GetParticleBoundaryBuffer();
        particle_buffers.clearParticles();
    }

    void warpx_depositChargeDensity (const char* char_species_name, int lev) {
        // this function is used to deposit a given species' charge density
        // in the rho_fp multifab which can then be accessed from python via
        // pywarpx.fields.RhoFPWrapper()
        WarpX& warpx = WarpX::GetInstance();
        const auto & mypc = warpx.GetPartContainer();
        const std::string species_name(char_species_name);
        auto & myspc = mypc.GetParticleContainerFromName(species_name);
        auto * rho_fp = warpx.get_pointer_rho_fp(lev);

        if (rho_fp == nullptr) {
            ablastr::warn_manager::WMRecordWarning(
                "WarpXWrappers", "rho_fp is not allocated",
                ablastr::warn_manager::WarnPriority::low
            );
            return;
        }

        for (WarpXParIter pti(myspc, lev); pti.isValid(); ++pti)
        {
            const long np = pti.numParticles();
            auto& wp = pti.GetAttribs(PIdx::w);
            // Do this unconditionally, ignoring myspc.do_not_deposit, to support diagnostic uses
            myspc.DepositCharge(pti, wp, nullptr, rho_fp, 0, 0, np, 0, lev, lev);
        }
#ifdef WARPX_DIM_RZ
        warpx.ApplyInverseVolumeScalingToChargeDensity(rho_fp, lev);
#endif
    }

    void warpx_ComputeDt () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.ComputeDt();
    }
    void warpx_MoveWindow (int step,bool move_j) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.MoveWindow(step, move_j);
    }

    void warpx_EvolveE (amrex::Real dt) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveE(dt);
    }
    void warpx_EvolveB (amrex::Real dt, DtType a_dt_type) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.EvolveB(dt, a_dt_type);
    }
    void warpx_FillBoundaryE () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.FillBoundaryE(warpx.getngEB());
    }
    void warpx_FillBoundaryB () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.FillBoundaryB(warpx.getngEB());
    }
    void warpx_SyncRho () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.SyncRho();
    }
    void warpx_SyncCurrent (
        const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_fp,
        const amrex::Vector<std::array<std::unique_ptr<amrex::MultiFab>,3>>& J_cp) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.SyncCurrent(J_fp, J_cp);
    }
    void warpx_UpdateAuxilaryData () {
        WarpX& warpx = WarpX::GetInstance();
        warpx.UpdateAuxilaryData();
    }
    void warpx_PushParticlesandDepose (amrex::Real cur_time) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.PushParticlesandDepose(cur_time);
    }

    int warpx_getistep (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.getistep(lev);
    }
    void warpx_setistep (int lev, int ii) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.setistep(lev, ii);
    }
    amrex::Real warpx_gett_new (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.gett_new(lev);
    }
    void warpx_sett_new (int lev, amrex::Real time) {
        WarpX& warpx = WarpX::GetInstance();
        warpx.sett_new(lev, time);
    }
    amrex::Real warpx_getdt (int lev) {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.getdt(lev);
    }

    int warpx_maxStep () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.maxStep();
    }
    amrex::Real warpx_stopTime () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.stopTime();
    }

    int warpx_finestLevel () {
        WarpX& warpx = WarpX::GetInstance();
        return warpx.finestLevel();
    }

    int warpx_getMyProc () {
        return amrex::ParallelDescriptor::MyProc();
    }

    int warpx_getNProcs () {
        return amrex::ParallelDescriptor::NProcs();
    }

    void warpx_setPotentialEB (const char * char_potential) {
        WarpX& warpx = WarpX::GetInstance();
        const std::string potential(char_potential);
        warpx.m_poisson_boundary_handler.setPotentialEB(potential);
    }

    void mypc_Redistribute () {
        auto & mypc = WarpX::GetInstance().GetPartContainer();
        mypc.Redistribute();
    }
