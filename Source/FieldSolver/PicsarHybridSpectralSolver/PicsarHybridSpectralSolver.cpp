/* Copyright 2019 Axel Huebl, Maxence Thevenet, Remi Lehe
 * Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"
#include "WarpX_f.H"
#include <AMReX_iMultiFab.H>
#include <fftw3-mpi.h>

using namespace amrex;

void
WarpX::AllocLevelDataFFT (int lev)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 0, "PSATD doesn't work with mesh refinement yet");

    InitFFTComm(lev);

    BoxArray ba_fp_fft;
    DistributionMapping dm_fp_fft;
    FFTDomainDecomposition(lev, ba_fp_fft, dm_fp_fft, ba_valid_fp_fft[lev], domain_fp_fft[lev],
                           geom[lev].Domain());

    // rho2 has one extra ghost cell, so that it's safe to deposit charge density after
    // pushing particle.

    Efield_fp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_fp_fft,Ex_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Efield_fp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_fp_fft,Ey_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Efield_fp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_fp_fft,Ez_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Bfield_fp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_fp_fft,Bx_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Bfield_fp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_fp_fft,By_nodal_flag),
                                             dm_fp_fft, 1, 0));
    Bfield_fp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_fp_fft,Bz_nodal_flag),
                                             dm_fp_fft, 1, 0));
    current_fp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_fp_fft,jx_nodal_flag),
                                              dm_fp_fft, 1, 0));
    current_fp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_fp_fft,jy_nodal_flag),
                                              dm_fp_fft, 1, 0));
    current_fp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_fp_fft,jz_nodal_flag),
                                              dm_fp_fft, 1, 0));
    rho_fp_fft[lev].reset(new MultiFab(amrex::convert(ba_fp_fft,IntVect::TheNodeVector()),
                                       dm_fp_fft, 2, 0));

    if (lev > 0)
    {
        BoxArray ba_cp_fft;
        DistributionMapping dm_cp_fft;
        FFTDomainDecomposition(lev, ba_cp_fft, dm_cp_fft, ba_valid_cp_fft[lev], domain_cp_fft[lev],
                               amrex::coarsen(geom[lev].Domain(),2));

        Efield_cp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_cp_fft,Ex_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Efield_cp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_cp_fft,Ey_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Efield_cp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_cp_fft,Ez_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Bfield_cp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_cp_fft,Bx_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Bfield_cp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_cp_fft,By_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        Bfield_cp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_cp_fft,Bz_nodal_flag),
                                                 dm_cp_fft, 1, 0));
        current_cp_fft[lev][0].reset(new MultiFab(amrex::convert(ba_cp_fft,jx_nodal_flag),
                                                  dm_cp_fft, 1, 0));
        current_cp_fft[lev][1].reset(new MultiFab(amrex::convert(ba_cp_fft,jy_nodal_flag),
                                                  dm_cp_fft, 1, 0));
        current_cp_fft[lev][2].reset(new MultiFab(amrex::convert(ba_cp_fft,jz_nodal_flag),
                                                  dm_cp_fft, 1, 0));
        rho_cp_fft[lev].reset(new MultiFab(amrex::convert(ba_cp_fft,IntVect::TheNodeVector()),
                                           dm_cp_fft, 2, 0));
    }

}

/** \brief Create MPI sub-communicators for each FFT group,
 *         and put them in PICSAR module
 *
 * These communicators are passed to the parallel FFTW library, in order
 * to perform a global FFT within each FFT group.
 */
void
WarpX::InitFFTComm (int lev)
{
    int nprocs = ParallelDescriptor::NProcs();
    ngroups_fft = std::min(ngroups_fft, nprocs);

    // # of processes in the subcommunicator
    int np_fft = nprocs / ngroups_fft;
// TODO:
//    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(np_fft*ngroups_fft == nprocs,
//        "Number of processes must be divisible by number of FFT groups");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngroups_fft == 1,
        "Number of FFT groups should be 1 at this point.");

    int myproc = ParallelDescriptor::MyProc();
    // my color in ngroups_fft subcommunicators.  0 <= color_fft < ngroups_fft
    color_fft[lev] = myproc / np_fft;
    MPI_Comm_split(ParallelDescriptor::Communicator(), color_fft[lev], myproc, &comm_fft[lev]);
}

/** \brief Perform domain decomposition for the FFTW
 *
 *  Attribute one (unique) box to each proc, in such a way that:
 *    - The global domain is divided among FFT groups,
 *      with additional guard cells around each FFT group
 *    - The domain associated to an FFT group (with its guard cells)
 *      is further divided in sub-subdomains along z, so as to distribute
 *      it among the procs within an FFT group
 *
 *  The attribution is done by setting (within this function):
 *  - ba_fft: the BoxArray representing the final set of sub-domains for the FFT
 *            (includes/covers the guard cells of the FFT groups)
 *  - dm_fft: the mapping between these sub-domains and the corresponding proc
 *              (imposes one unique box for each proc)
 *  - ba_valid: the BoxArray that contains valid part of the sub-domains of ba_fft
 *            (i.e. does not include/cover the guard cells of the FFT groups)
 *  - domain_fft: a Box that represent the domain of the FFT group for the current proc
 */
void
WarpX::FFTDomainDecomposition (int lev, BoxArray& ba_fft, DistributionMapping& dm_fft,
                               BoxArray& ba_valid, Box& domain_fft, const Box& domain)
{

    IntVect nguards_fft(AMREX_D_DECL(nox_fft/2,noy_fft/2,noz_fft/2));

    int nprocs = ParallelDescriptor::NProcs();

    BoxList bl(domain, ngroups_fft);  // This does a multi-D domain decomposition for groups
    AMREX_ALWAYS_ASSERT(bl.size() == ngroups_fft);
    const Vector<Box>& bldata = bl.data();

    // This is the domain for the FFT sub-group (including guard cells)
    domain_fft = amrex::grow(bldata[color_fft[lev]], nguards_fft);
    // Ask FFTW to chop the current FFT sub-group domain in the z-direction
    // and give a chunk to each MPI rank in the current sub-group.
    int nz_fft, z0_fft;
    const ptrdiff_t nx_global = domain_fft.length(0)+1;
    const ptrdiff_t ny_global = domain_fft.length(1)+1;
    const ptrdiff_t nz_global = domain_fft.length(2)+1;
    ptrdiff_t local_n0, local_0_start;
    auto alloc_local = fftw_mpi_local_size_3d(
        nz_global, ny_global, nx_global/2+1, comm_fft[lev],
        &local_n0, &local_0_start);
    nz_fft = local_n0;
    z0_fft = local_0_start;

    // Each MPI rank adds a box with its chunk of the FFT grid
    // (given by the above decomposition) to the list `bx_fft`,
    // then list is shared among all MPI ranks via AllGather
    Vector<Box> bx_fft;
    if (nz_fft > 0) {
        Box b = domain_fft;
        b.setRange(AMREX_SPACEDIM-1, z0_fft+domain_fft.smallEnd(AMREX_SPACEDIM-1), nz_fft);
        bx_fft.push_back(b);
    } else {
        // Add empty box for the AllGather call
        bx_fft.push_back(Box());
    }
    amrex::AllGatherBoxes(bx_fft);
    AMREX_ASSERT(bx_fft.size() == ParallelDescriptor::NProcs());
    // Build pmap and bx_fft without the empty boxes
    Vector<int> pmap;
    for (int i = 0; i < bx_fft.size(); ++i) {
        if (bx_fft[i].ok()) {
            pmap.push_back(i);
        }
    }
    bx_fft.erase(std::remove_if(bx_fft.begin(),bx_fft.end(),
                                [](Box const& b) { return b.isEmpty(); }),
                 bx_fft.end());
    AMREX_ASSERT(bx_fft.size() == pmap.size());

    // Define the AMReX objects for the FFT grid: BoxArray and DistributionMapping
    ba_fft.define(BoxList(std::move(bx_fft)));
    dm_fft.define(std::move(pmap));

    // For communication between WarpX normal domain and FFT domain, we need to create a
    // special BoxArray ba_valid
    const Box foobox(-nguards_fft-2, -nguards_fft-2);

    BoxList bl_valid; // List of boxes: will be filled by the valid part of the subdomains of ba_fft
    bl_valid.reserve(ba_fft.size());
    int np_fft = nprocs / ngroups_fft;
    for (int i = 0; i < ba_fft.size(); ++i)
    {
        int igroup = dm_fft[i] / np_fft; // This should be consistent with InitFFTComm
        const Box& bx = ba_fft[i] & bldata[igroup]; // Intersection with the domain of
                                                    // the FFT group *without* guard cells
        if (bx.ok())
        {
            bl_valid.push_back(bx);
        }
        else
        {
            bl_valid.push_back(foobox);
        }
    }

    ba_valid.define(std::move(bl_valid));
}


void
WarpX::FreeFFT (int lev)
{
    fftw_mpi_cleanup();
    if (comm_fft[lev] != MPI_COMM_NULL) {
        MPI_Comm_free(&comm_fft[lev]);
    }
    comm_fft[lev] = MPI_COMM_NULL;
}
