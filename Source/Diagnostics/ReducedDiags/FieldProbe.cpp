/* Copyright 2021 Lorenzo Giacomel
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldProbe.H"
#include "Particles/Gather/FieldGather.H"

#include "Utils/IntervalsParser.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_Array.H>
#include <AMReX_Config.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <AMReX_Reduce.H>
#include <AMReX_Geometry.H>

#include <cmath>
#include <ostream>

using namespace amrex;

// constructor
FieldProbe::FieldProbe (std::string rd_name)
: ReducedDiags{rd_name}
{
    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "FieldProbe reduced diagnostics does not work for RZ coordinate.");
#endif

    // read number of levels
    int nLevel = 0;
    amrex::ParmParse pp_amr("amr");
    pp_amr.query("max_level", nLevel);
    nLevel += 1;

    amrex::ParmParse pp_rd_name(rd_name);
    getWithParser(pp_rd_name, "x_probe", x_probe);
#if (AMREX_SPACEDIM == 3)
    getWithParser(pp_rd_name, "y_probe", y_probe);
#endif
    getWithParser(pp_rd_name, "z_probe", z_probe);

    pp_rd_name.query("raw_fields", raw_fields);

    pp_rd_name.query("interp_order", interp_order);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(interp_order <= WarpX::nox ,
                                     "Field probe interp_order should be lower than algo.particle_shape");
    // resize data array
    m_data.resize(noutputs*nLevel, 0.0_rt);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs{m_path + m_rd_name + "." + m_extension, std::ofstream::out};
            // write header row
            int c = 0;
            ofs << "#";
            ofs << "[" << c++ << "]step()";
            ofs << m_sep;
            ofs << "[" << c++ << "]time(s)";
            for (int lev = 0; lev < nLevel; ++lev)
            {
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ex_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ey_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ez_lev" + std::to_string(lev) + " (V/m)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Bx_lev" + std::to_string(lev) + " (T)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_By_lev" + std::to_string(lev) + " (T)";
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Bz_lev" + std::to_string(lev) + " (T)";
            }
            ofs << std::endl;
            // close file
            ofs.close();
        }
    }
}
// end constructor

// function that computes field values at probe position
void FieldProbe::ComputeDiags (int step)
{
    // Judge if the diags should be done
    if (!m_intervals.contains(step+1)) { return; }

    // get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();

    // get number of level
    const auto nLevel = warpx.finestLevel() + 1;

    // vector where we store the field values
    amrex::Vector<amrex::Real> fp_values(noutputs, 0);

    // loop over refinement levels
    for (int lev = 0; lev < nLevel; ++lev)
    {

        const amrex::Geometry& gm = warpx.Geom(lev);
        const auto prob_lo = gm.ProbLo();
        const auto prob_hi = gm.ProbHi();

#if (AMREX_SPACEDIM == 2)
        m_probe_in_domain = x_probe >= prob_lo[0] and x_probe < prob_hi[0] and
                            z_probe >= prob_lo[1] and z_probe < prob_hi[1];
#else
        m_probe_in_domain = x_probe >= prob_lo[0] and x_probe < prob_hi[0] and
                            y_probe >= prob_lo[1] and y_probe < prob_hi[1] and
                            z_probe >= prob_lo[2] and z_probe < prob_hi[2];
#endif
        if(lev == 0) m_probe_in_domain_lev_0 = m_probe_in_domain;

        if( !m_probe_in_domain ) break;

        const auto cell_size = gm.CellSizeArray();

        const int i_probe = static_cast<int>(amrex::Math::floor((x_probe - prob_lo[0]) / cell_size[0]));
#if (AMREX_SPACEDIM == 2)
        const int j_probe = static_cast<int>(amrex::Math::floor((z_probe - prob_lo[1]) / cell_size[1]));
        const int k_probe = 0;
#elif(AMREX_SPACEDIM == 3)
        const int j_probe = static_cast<int>(amrex::Math::floor((y_probe - prob_lo[1]) / cell_size[1]));
        const int k_probe = static_cast<int>(amrex::Math::floor((z_probe - prob_lo[2]) / cell_size[2]));
#endif
        // get MultiFab data at lev
        const amrex::MultiFab &Ex = warpx.getEfield(lev, 0);
        const amrex::MultiFab &Ey = warpx.getEfield(lev, 1);
        const amrex::MultiFab &Ez = warpx.getEfield(lev, 2);
        const amrex::MultiFab &Bx = warpx.getBfield(lev, 0);
        const amrex::MultiFab &By = warpx.getBfield(lev, 1);
        const amrex::MultiFab &Bz = warpx.getBfield(lev, 2);

        // Prepare interpolation of field components to probe_position
        // The arrays below store the index type (staggering) of each MultiFab.
        amrex::IndexType const Extype = Ex.ixType();
        amrex::IndexType const Eytype = Ey.ixType();
        amrex::IndexType const Eztype = Ez.ixType();
        amrex::IndexType const Bxtype = Bx.ixType();
        amrex::IndexType const Bytype = By.ixType();
        amrex::IndexType const Bztype = Bz.ixType();

        int probe_proc = -1;
        // MFIter loop to interpolate fields to probe position
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(Ex, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box &box = amrex::enclosedCells(mfi.nodaltilebox());

            if (box.contains(i_probe, j_probe, k_probe)) {
                // Make the box cell centered in preparation for the interpolation (and to avoid
                // including ghost cells in the calculation)
                const auto &arrEx = Ex[mfi].array();
                const auto &arrEy = Ey[mfi].array();
                const auto &arrEz = Ez[mfi].array();
                const auto &arrBx = Bx[mfi].array();
                const auto &arrBy = By[mfi].array();
                const auto &arrBz = Bz[mfi].array();

                amrex::Real Ex_interp = 0._rt, Ey_interp = 0._rt, Ez_interp = 0._rt;
                amrex::Real Bx_interp = 0._rt, By_interp = 0._rt, Bz_interp = 0._rt;

                amrex::Vector<amrex::Real> v_galilean{amrex::Vector<amrex::Real>(3, amrex::Real(0.))};
                const auto &xyzmin = WarpX::GetInstance().LowerCornerWithGalilean(box, v_galilean,
                                                                                  lev);
                const std::array<Real, 3> &dx = WarpX::CellSize(lev);
                const amrex::GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
                const amrex::GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};

                //Interpolating to the probe point
#if (AMREX_SPACEDIM == 2)
                constexpr amrex::Real y_probe = 0._rt;
#endif
                doGatherShapeN(x_probe, y_probe, z_probe, Ex_interp, Ey_interp, Ez_interp,
                               Bx_interp, By_interp, Bz_interp, arrEx, arrEy, arrEz, arrBx,
                               arrBy, arrBz, Extype, Eytype, Eztype, Bxtype, Bytype, Bztype, dx_arr,
                               xyzmin_arr, amrex::lbound(box), WarpX::n_rz_azimuthal_modes,
                               interp_order, false);

                // Either save the interpolated fields or the raw fields depending on the raw_fields flag
                fp_values[index_Ex] = raw_fields ? arrEx(i_probe, j_probe, k_probe) : Ex_interp;
                fp_values[index_Ey] = raw_fields ? arrEy(i_probe, j_probe, k_probe) : Ey_interp;
                fp_values[index_Ez] = raw_fields ? arrEz(i_probe, j_probe, k_probe) : Ez_interp;
                fp_values[index_Bx] = raw_fields ? arrBx(i_probe, j_probe, k_probe) : Bx_interp;
                fp_values[index_By] = raw_fields ? arrBy(i_probe, j_probe, k_probe) : By_interp;
                fp_values[index_Bz] = raw_fields ? arrBz(i_probe, j_probe, k_probe) : Bz_interp;

                probe_proc = amrex::ParallelDescriptor::MyProc();
            }
        }

        //All the processors have probe_proc = -1 except the one that contains the point, which
        //has probe_proc equal to a number >=0. Therefore, ReduceIntMax communicates to all the
        //processors the rank of the processor which contains the point
        amrex::ParallelDescriptor::ReduceIntMax(probe_proc);

        if(probe_proc != amrex::ParallelDescriptor::IOProcessorNumber() and probe_proc != -1) {
            if (amrex::ParallelDescriptor::MyProc() == probe_proc) {
                amrex::ParallelDescriptor::Send(fp_values.dataPtr(),
                                                noutputs,
                                                amrex::ParallelDescriptor::IOProcessorNumber(),
                                                0);
            }
            if (amrex::ParallelDescriptor::MyProc()
            == amrex::ParallelDescriptor::IOProcessorNumber()) {
                amrex::ParallelDescriptor::Recv(fp_values.dataPtr(), noutputs, probe_proc, 0);
            }
        }

        // Fill output array
        m_data[lev * noutputs + index_Ex] = fp_values[index_Ex];
        m_data[lev * noutputs + index_Ey] = fp_values[index_Ey];
        m_data[lev * noutputs + index_Ez] = fp_values[index_Ez];
        m_data[lev * noutputs + index_Bx] = fp_values[index_Bx];
        m_data[lev * noutputs + index_By] = fp_values[index_By];
        m_data[lev * noutputs + index_Bz] = fp_values[index_Bz];
    }
    // end loop over refinement levels

    /* m_data now contains up-to-date values for:
     *  [probe(Ex),probe(Ey),probe(Ez),
     *   probe(Bx),probe(By),probe(Bz)] */
}
// end void FieldProbe::ComputeDiags

void FieldProbe::WriteToFile (int step) const
{
    if(m_probe_in_domain_lev_0){
        ReducedDiags::WriteToFile (step);
    }
}
