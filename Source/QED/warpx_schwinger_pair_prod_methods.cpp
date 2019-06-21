//This file impements methods of WarpX to provide Schwinger Pair Production process

#include "WarpX.H"
#include "schwinger_pair_prod_wrapper.h"

void WarpX::DoSchwingerPairProduction(amrex::Real a_dt)
{
    //Takes the fields of the coarsest grid
    const amrex::MultiFab&  Exmf = *Efield_fp[0][0].get();
    const amrex::MultiFab&  Eymf = *Efield_fp[0][1].get();
    const amrex::MultiFab&  Ezmf = *Efield_fp[0][2].get();
    const amrex::MultiFab&  Bxmf = *Bfield_fp[0][0].get();
    const amrex::MultiFab&  Bymf = *Bfield_fp[0][1].get();
    const amrex::MultiFab&  Bzmf = *Bfield_fp[0][2].get();

    //Cell dimensions
    const amrex::Real* dx = geom[0].CellSize();

    amrex::iMultiFab nparticles(amrex::convert(Exmf.boxArray(),amrex::IntVect::TheCellVector()),
                    Exmf.DistributionMap(), 1, 0);

    nparticles.setVal(0);

    for (amrex::MFIter mfi(Exmf); mfi.isValid(); ++mfi){
        const auto& bx = amrex::enclosedCells(mfi.validbox());
        const auto& Ex = Exmf.array(mfi);
        const auto& Ey = Eymf.array(mfi);
        const auto& Ez = Ezmf.array(mfi);
        const auto& Bx = Bxmf.array(mfi);
        const auto& By = Bymf.array(mfi);
        const auto& Bz = Bzmf.array(mfi);
        const auto& np = nparticles.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Averaging fields to cell centers
            amrex::Real Bxc = 0.5*(Bx(i,j,k) + Bx(i+1,j  ,k  ));
            amrex::Real Byc = 0.5*(By(i,j,k) + By(i  ,j+1,k  ));
            amrex::Real Bzc = 0.5*(Bz(i,j,k) + Bz(i  ,j  ,k+1));
            amrex::Real Exc = 0.25*(Ex(i,j,k)+Ex(i  ,j+1,k)+Ex(i,j  ,k+1)+Ex(i  ,j+1,k+1));
            amrex::Real Eyc = 0.25*(Ey(i,j,k)+Ey(i+1,j  ,k)+Ey(i,j  ,k+1)+Ey(i+1,j  ,k+1));
            amrex::Real Ezc = 0.25*(Ez(i,j,k)+Ez(i+1,j  ,k)+Ez(i,j+1,k  )+Ez(i+1,j+1,k  ));


            size_t num_particles_to_add;
            amrex::Real weight; //Weight will be 1/volume
            warpx_schwinger_pair_engine::internal_generate_pairs_multiple(
                Exc, Eyc, Ezc, Bxc, Byc, Bzc,
                dx[0], dx[1], dx[2], a_dt,
                &num_particles_to_add, &weight,
                1.0, //Lambda is irrelevant. Since SI units are used, it can be set to any value
                amrex::Random()
                );

                np(i,j,k) += num_particles_to_add;
        });
    }
}
