
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>

#include <WarpX.H>

using namespace amrex;

void
WarpX::InitSpaceChargeField ()
{

    // Get parameters of the grid
    const int lev = 0;
    const Geometry& gm = Geom(lev);
    const DistributionMapping& dm = DistributionMap(lev);
    BoxArray nba = boxArray(lev);
    nba.surroundingNodes();
    const Box& gbox = gm.Domain();

    // Get particle average momentum
    // TODO: Do a global reduction
    Real uz_avg = 1.;
    Real gamma = std::sqrt( 1 + uz_avg*uz_avg );
    Real beta_z = uz_avg/gamma;

    // Create a new boosted geometry object, that will be used to
    // solve the relativistic Poisson equation
    RealBox domain = gm.ProbDomain();
    const int dir = AMREX_SPACEDIM-1; // boost along the last axis (z direction)
    domain.setHi( dir, gamma*domain.hi(dir) );
    domain.setLo( dir, gamma*domain.lo(dir) );
    auto boosted_geom = Geometry( gm.Domain(), &domain );

    // Deposit particle charge density (source of Poisson solver)
    bool local = false;
    const std::unique_ptr<MultiFab>& rho = mypc->GetChargeDensity(lev, local);

    // Allocate the field for the potential
    MultiFab phi(nba, dm, 1, 0);
    phi.setVal(0.);

#ifndef USE_OPENBC_POISSON

    // Call amrex's multigrid solver

    // Define the boundary conditions
    Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
    for (int idim=0; idim<AMREX_SPACEDIM; idim++){
        if ( gm.isPeriodic(idim) ) {
            lobc[idim] = LinOpBCType::Periodic;
            hibc[idim] = LinOpBCType::Periodic;
        } else {
            // Use Dirichlet boundary condition by default.
            // Ideally, we would often want open boundary conditions here.
            lobc[idim] = LinOpBCType::Dirichlet;
            hibc[idim] = LinOpBCType::Dirichlet;
        }
    }

    // Define the linear operator (Poisson operator)
    MLNodeLaplacian linop( {boosted_geom}, {nba}, {dm} );
    linop.setDomainBC( lobc, hibc );
    BoxArray cba = nba;
    cba.enclosedCells();
    MultiFab sigma(cba, dm, 1, 0);
    sigma.setVal(PhysConst::ep0);
    linop.setSigma(0, sigma);

    // Solve the Poisson equation
    MLMG mlmg(linop);
    const Real reltol = 1.e-3;
    Vector<MultiFab*> phi_vec;
    phi_vec.resize(1);
    phi_vec[0] = &phi;
    Vector<MultiFab const*> rho_vec;
    rho_vec.resize(1);
    rho_vec[0] = &(*rho);
    mlmg.solve( phi_vec, rho_vec, reltol, 0.0);

#else
    // Call the openBC Poisson solver
    int lohi[6];
    warpx_openbc_decompose(gbox.loVect(), gbox.hiVect(), lohi, lohi+3);
    int nprocs = ParallelDescriptor::NProcs();
    int myproc = ParallelDescriptor::MyProc();

    Vector<int> alllohi(6*nprocs,100000);
    MPI_Allgather(lohi, 6, MPI_INT, alllohi.data(), 6, MPI_INT,
                    ParallelDescriptor::Communicator());
    BoxList bl{IndexType::TheNodeType()};
    for (int i = 0; i < nprocs; ++i)
    {
        bl.push_back(Box(IntVect(alllohi[6*i  ],alllohi[6*i+1],alllohi[6*i+2]),
                         IntVect(alllohi[6*i+3],alllohi[6*i+4],alllohi[6*i+5]),
                         IndexType::TheNodeType()));
    }
    BoxArray ba{bl};

    Vector<int> iprocmap(nprocs+1);
    std::iota(iprocmap.begin(), iprocmap.end(), 0);
    iprocmap.back() = myproc;

    DistributionMapping dm{iprocmap};

    MultiFab rho_openbc(ba, dm, 1, 0);
    MultiFab phi_openbc(ba, dm, 1, 0);

    rho_openbc.setVal(0.0);
    rho_openbc.copy(*rho, 0, 0, 1, rho->nGrow(), 0, gm.periodicity(), FabArrayBase::ADD);

    static_assert(AMREX_SPACEDIM == 3, "Openbc is 3D only");
    BL_ASSERT(finestLevel() == 0);
    const Real* dx = boosted_geom.CellSize();
    warpx_openbc_potential(rho_openbc[myproc].dataPtr(), phi_openbc[myproc].dataPtr(), dx);
    phi.copy(phi_openbc, gm.periodicity());
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(phi, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Real* dx = gm.CellSize();
        const Real inv_dx = 1./dx[0];
#if (AMREX_SPACEDIM == 3)
        const Real inv_dy = 1./dx[1];
        const Real inv_dz = 1./dx[2];
#else
        const Real inv_dz = 1./dx[1];
#endif
        const Real inv_gamma2 = 1./(gamma*gamma);

        const Box& tbx  = mfi.tilebox(Ex_nodal_flag);
        const Box& tby  = mfi.tilebox(Ey_nodal_flag);
        const Box& tbz  = mfi.tilebox(Ez_nodal_flag);

        const auto& phi_arr = phi[mfi].array();
        const auto& Ex_arr = (*Efield_fp[lev][0])[mfi].array();
        const auto& Ey_arr = (*Efield_fp[lev][1])[mfi].array();
        const auto& Ez_arr = (*Efield_fp[lev][2])[mfi].array();

#if (AMREX_SPACEDIM == 3)
        amrex::ParallelFor( tbx, tby, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Ex_arr(i,j,k) += inv_dx*( phi_arr(i+1,j,k) - phi_arr(i,j,k) );
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Ey_arr(i,j,k) += inv_dy*( phi_arr(i,j+1,k) - phi_arr(i,j,k) );
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Ez_arr(i,j,k) += inv_gamma2*inv_dz*( phi_arr(i,j,k+1) - phi_arr(i,j,k) );
            }
        );
#else
        amrex::ParallelFor( tbx, tbz,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Ex_arr(i,j,k) += inv_dx*( phi_arr(i+1,j,k) - phi_arr(i,j,k) );
            },
            [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                Ez_arr(i,j,k) += inv_gamma2*inv_dz*( phi_arr(i,j+1,k) - phi_arr(i,j,k) );
            }
        );
#endif
    }
}
