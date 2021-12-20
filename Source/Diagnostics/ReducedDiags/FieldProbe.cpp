/* Copyright 2021 Lorenzo Giacomel, Tiberius Rheaume, Axel Huebl
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldProbe.H"
#include "FieldProbeParticleContainer.H"
#include "Particles/Gather/FieldGather.H"

#include "Utils/IntervalsParser.H"
#include "Utils/WarpXConst.H"
#include "Utils/WarpXUtil.H"
#include "WarpX.H"

#include <AMReX_Array.H>
#include <AMReX_Config.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParIter.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <AMReX_Reduce.H>
#include <AMReX_Geometry.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Vector.H>

#include <cmath>
#include <ostream>
#include <string>
#include <unordered_map>

using namespace amrex;

// constructor

FieldProbe::FieldProbe (std::string rd_name)
: ReducedDiags{rd_name}, m_probe(&WarpX::GetInstance())
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

    /* Obtain input data from parsing inputs file.
     * For the case of a single particle:
     *     Define x, y, and z of particle
     *     Define whether or not to integrate fields
     * For the case of a line detector:
     *     Define x, y, and z of end of line point 1
     *     Define x, y, and z of end of line point 2
     *     Define resolution to determine number of particles
     *     Define whether ot not to integrate fields
     * For the case of a plane detector:
     *     Define a vector normal to the detector plane
     *     Define a vector in the "up" direction of the plane
     *     Define the size of the plane (width of half square)
     *     Define resolution to determine number of particles
     *     Define whether ot not to integrate fields
     */
    amrex::ParmParse pp_rd_name(rd_name);
    getWithParser(pp_rd_name, "x_probe", x_probe);
    getWithParser(pp_rd_name, "y_probe", y_probe);
    getWithParser(pp_rd_name, "z_probe", z_probe);
    pp_rd_name.query("integrate", m_field_probe_integrate);
    pp_rd_name.query("raw_fields", raw_fields);
    pp_rd_name.query("interp_order", interp_order);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(interp_order <= WarpX::nox ,
                                     "Field probe interp_order should be less than or equal to algo.particle_shape");
    // resize data array
    m_data.resize(noutputs, 0.0_rt);

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
            // maps FieldProbe observables to units
            std::unordered_map< int, std::string > u_map;

            if (m_field_probe_integrate)
            {
                u_map =
                {
                    {FieldProbePIdx::Ex , " (V*s/m) "},
                    {FieldProbePIdx::Ey , " (V*s/m) "},
                    {FieldProbePIdx::Ez , " (V*s/m) "},
                    {FieldProbePIdx::Bx , " (T*s) "},
                    {FieldProbePIdx::By , " (T*s) "},
                    {FieldProbePIdx::Bz , " (T*s) "},
                    {FieldProbePIdx::S , " (W*s/m^2) "}
                };
            }
            else
            {
                u_map =
                {
                    {FieldProbePIdx::Ex , " (V/m) "},
                    {FieldProbePIdx::Ey , " (V/m) "},
                    {FieldProbePIdx::Ez , " (V/m) "},
                    {FieldProbePIdx::Bx , " (T) "},
                    {FieldProbePIdx::By , " (T) "},
                    {FieldProbePIdx::Bz , " (T) "},
                    {FieldProbePIdx::S , " (W/m^2) "}
                };
            }
            for (int lev = 0; lev < nLevel; ++lev)
            {
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ex_lev" + std::to_string(lev) + u_map[FieldProbePIdx::Ex];
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ey_lev" + std::to_string(lev) + u_map[FieldProbePIdx::Ey];
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Ez_lev" + std::to_string(lev) + u_map[FieldProbePIdx::Ez];
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Bx_lev" + std::to_string(lev) + u_map[FieldProbePIdx::Bx];
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_By_lev" + std::to_string(lev) + u_map[FieldProbePIdx::By];
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_Bz_lev" + std::to_string(lev) + u_map[FieldProbePIdx::Bz];
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_S_lev" + std::to_string(lev) + u_map[FieldProbePIdx::S]; //update all units if integrating (might be energy / m^2)
            }
            ofs << std::endl;

            // close file
            ofs.close();
        }
    }
} // end constructor

void FieldProbe::InitData ()
{
    //create 1D array for X, Y, and Z of particles
    amrex::Vector<amrex::ParticleReal> xpos;
    amrex::Vector<amrex::ParticleReal> ypos;
    amrex::Vector<amrex::ParticleReal> zpos;

    // for now, only one MPI rank adds a probe particle
    if (ParallelDescriptor::IOProcessor())
    {
        xpos.push_back(x_probe);
        ypos.push_back(y_probe);
        zpos.push_back(z_probe);
    }

    // add np particles on lev 0 to m_probe
    m_probe.AddNParticles(0, xpos, ypos, zpos);
}

void FieldProbe::LoadBalance() {
    m_probe.Redistribute();
}

bool FieldProbe::ProbeInDomain () const
{
    // get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();
    int const lev = 0;
    const amrex::Geometry& gm = warpx.Geom(lev);
    const auto prob_lo = gm.ProbLo();
    const auto prob_hi = gm.ProbHi();

    /*
     * Determine if probe exists within simulation boundaries. During 2D simulations,
     * y values will be set to 0 making it unnecessary to check. Generally, the second
     * value in a position array will be the y value, but in the case of 2D, prob_lo[1]
     * and prob_hi[1] refer to z. This is a result of warpx.Geom(lev).
     */
#if defined(WARPX_DIM_1D_Z)
    return z_probe >= prob_lo[1] && z_probe < prob_hi[1];
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    return x_probe >= prob_lo[0] && x_probe < prob_hi[0] &&
           z_probe >= prob_lo[1] && z_probe < prob_hi[1];
#else
    return x_probe >= prob_lo[0] && x_probe < prob_hi[0] &&
           y_probe >= prob_lo[1] && y_probe < prob_hi[1] &&
           z_probe >= prob_lo[2] && z_probe < prob_hi[2];
#endif
}

void FieldProbe::ComputeDiags (int step)
{
    // Judge if the diags should be done
    if (!m_field_probe_integrate)
    {
        if (!m_intervals.contains(step+1)) { return; }
    }
    // get a reference to WarpX instance
    auto & warpx = WarpX::GetInstance();

    // get number of mesh-refinement levels
    const auto nLevel = warpx.finestLevel() + 1;

    // loop over refinement levels
    for (int lev = 0; lev < nLevel; ++lev)
    {
        const amrex::Geometry& gm = warpx.Geom(lev);
        const auto prob_lo = gm.ProbLo();

        // get MultiFab data at lev
        const amrex::MultiFab &Ex = warpx.getEfield(lev, 0);
        const amrex::MultiFab &Ey = warpx.getEfield(lev, 1);
        const amrex::MultiFab &Ez = warpx.getEfield(lev, 2);
        const amrex::MultiFab &Bx = warpx.getBfield(lev, 0);
        const amrex::MultiFab &By = warpx.getBfield(lev, 1);
        const amrex::MultiFab &Bz = warpx.getBfield(lev, 2);

        /*
         * Prepare interpolation of field components to probe_position
         * The arrays below store the index type (staggering) of each MultiFab.
         */
        amrex::IndexType const Extype = Ex.ixType();
        amrex::IndexType const Eytype = Ey.ixType();
        amrex::IndexType const Eztype = Ez.ixType();
        amrex::IndexType const Bxtype = Bx.ixType();
        amrex::IndexType const Bytype = By.ixType();
        amrex::IndexType const Bztype = Bz.ixType();

        //defined for use in determining which CPU contains the particles
        int probe_proc = -1;

        // loop over each particle
        // TODO: add OMP parallel as in PhysicalParticleContainer::Evolve
        using MyParIter = FieldProbeParticleContainer::iterator;
        for (MyParIter pti(m_probe, lev); pti.isValid(); ++pti)
        {
            const auto getPosition = GetParticlePosition(pti);

            if( ProbeInDomain() )
            {
                const auto cell_size = gm.CellSizeArray();
                const int i_probe = static_cast<int>(amrex::Math::floor((x_probe - prob_lo[0]) / cell_size[0]));
#if defined(WARPX_DIM_1D_Z)
                const int j_probe = 0;
                const int k_probe = 0;
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
                const int j_probe = static_cast<int>(amrex::Math::floor((z_probe - prob_lo[1]) / cell_size[1]));
                const int k_probe = 0;
#elif defined(WARPX_DIM_3D)
                const int j_probe = static_cast<int>(amrex::Math::floor((y_probe - prob_lo[1]) / cell_size[1]));
                const int k_probe = static_cast<int>(amrex::Math::floor((z_probe - prob_lo[2]) / cell_size[2]));
#endif
                const auto &arrEx = Ex[pti].array();
                const auto &arrEy = Ey[pti].array();
                const auto &arrEz = Ez[pti].array();
                const auto &arrBx = Bx[pti].array();
                const auto &arrBy = By[pti].array();
                const auto &arrBz = Bz[pti].array();

                /*
                 * Make the box cell centered in preparation for the interpolation (and to avoid
                 * including ghost cells in the calculation)
                 */
                amrex::Box box = pti.tilebox();
                box.grow(Ex.nGrowVect());

                //preparing to write data to particle
                auto& attribs = pti.GetStructOfArrays().GetRealData();
                ParticleReal* const AMREX_RESTRICT part_Ex = attribs[FieldProbePIdx::Ex].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_Ey = attribs[FieldProbePIdx::Ey].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_Ez = attribs[FieldProbePIdx::Ez].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_Bx = attribs[FieldProbePIdx::Bx].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_By = attribs[FieldProbePIdx::By].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_Bz = attribs[FieldProbePIdx::Bz].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_S = attribs[FieldProbePIdx::S].dataPtr();

                amrex::Vector<amrex::Real> v_galilean{amrex::Vector<amrex::Real>(3, amrex::Real(0.))};
                const auto &xyzmin = WarpX::GetInstance().LowerCornerWithGalilean(box, v_galilean, lev);
                const std::array<Real, 3> &dx = WarpX::CellSize(lev);

                const amrex::GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
                const amrex::GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};
                const Dim3 lo = lbound(box);

                // Interpolating to the probe positions for each particle
                // Temporarily defining modes and interp outside ParallelFor to avoid GPU compilation errors.
                const int temp_modes = WarpX::n_rz_azimuthal_modes;
                const int temp_interp_order = interp_order;
                const bool temp_raw_fields = raw_fields;
                const bool temp_field_probe_integrate = m_field_probe_integrate;
                amrex::Real const dt = WarpX::GetInstance().getdt(lev);

                long const np = pti.numParticles();
                amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
                {
                    amrex::ParticleReal xp, yp, zp;
                    getPosition(ip, xp, yp, zp);

                    amrex::ParticleReal Exp = 0._rt, Eyp = 0._rt, Ezp = 0._rt;
                    amrex::ParticleReal Bxp = 0._rt, Byp = 0._rt, Bzp = 0._rt;

                    // first gather E and B to the particle positions
                    if (temp_raw_fields)
                    {
                        Exp = arrEx(i_probe, j_probe, k_probe);
                        Eyp = arrEy(i_probe, j_probe, k_probe);
                        Ezp = arrEz(i_probe, j_probe, k_probe);
                        Bxp = arrBx(i_probe, j_probe, k_probe);
                        Byp = arrBy(i_probe, j_probe, k_probe);
                        Bzp = arrBz(i_probe, j_probe, k_probe);
                    }
                    else
                        doGatherShapeN(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                   arrEx, arrEy, arrEz, arrBx, arrBy, arrBz,
                                   Extype, Eytype, Eztype, Bxtype, Bytype, Bztype,
                                   dx_arr, xyzmin_arr, lo, temp_modes,
                                   temp_interp_order, false);

                    //Calculate the Poynting Vector S
                    amrex::Real const sraw[3]{
                        Exp * Bzp - Ezp * Byp,
                        Ezp * Bxp - Exp * Bzp,
                        Exp * Byp - Eyp * Bxp
                    };
                    amrex::ParticleReal const S = (1._rt / PhysConst::mu0)  * sqrt(sraw[0] * sraw[0] + sraw[1] * sraw[1] + sraw[2] * sraw[2]);

                    /*
                     * Determine whether or not to integrate field data.
                     * If not integrating, store instantaneous values.
                     */
                    if (temp_field_probe_integrate)
                    {
                        // store values on particles
                        part_Ex[ip] += Exp * dt; //remember to add lorentz transform
                        part_Ey[ip] += Eyp * dt; //remember to add lorentz transform
                        part_Ez[ip] += Ezp * dt; //remember to add lorentz transform
                        part_Bx[ip] += Bxp * dt; //remember to add lorentz transform
                        part_By[ip] += Byp * dt; //remember to add lorentz transform
                        part_Bz[ip] += Bzp * dt; //remember to add lorentz transform
                        part_S[ip] += S * dt; //remember to add lorentz transform
                    }
                    else
                    {
                        // Either save the interpolated fields or the raw fields depending on the raw_fields flag
                        part_Ex[ip] = Exp; //remember to add lorentz transform
                        part_Ey[ip] = Eyp; //remember to add lorentz transform
                        part_Ez[ip] = Ezp; //remember to add lorentz transform
                        part_Bx[ip] = Bxp; //remember to add lorentz transform
                        part_By[ip] = Byp; //remember to add lorentz transform
                        part_Bz[ip] = Bzp; //remember to add lorentz transform
                        part_S[ip] = S; //remember to add lorentz transform
                    }
                });// ParallelFor Close

                // this check is here because for m_field_probe_integrate == True, we always compute
                // but we only write when we truly are in an output interval step
                if (m_intervals.contains(step+1)) {
                    for (int ip = 0; ip < np; ip++) {
                        // Fill output array
                        m_data[ip * noutputs + FieldProbePIdx::Ex] = part_Ex[ip];
                        m_data[ip * noutputs + FieldProbePIdx::Ey] = part_Ey[ip];
                        m_data[ip * noutputs + FieldProbePIdx::Ez] = part_Ez[ip];
                        m_data[ip * noutputs + FieldProbePIdx::Bx] = part_Bx[ip];
                        m_data[ip * noutputs + FieldProbePIdx::By] = part_By[ip];
                        m_data[ip * noutputs + FieldProbePIdx::Bz] = part_Bz[ip];
                        m_data[ip * noutputs + FieldProbePIdx::S] = part_S[ip];
                    }
                    /* m_data now contains up-to-date values for:
                     *  [Ex, Ey, Ez, Bx, By, Bz, and S] */
                }

                // do we have the one and only probe particle on our MPI rank?
                if (np > 0)
                    probe_proc = amrex::ParallelDescriptor::MyProc();
            }

        } // end particle iterator loop

        // make sure data is in m_data
        Gpu::synchronize();

        // this check is here because for m_field_probe_integrate == True, we always compute
        // but we only write when we truly are in an output interval step
        if (m_intervals.contains(step+1)) {
            /*
             * All the processors have probe_proc = -1 except the one that contains the point, which
             * has probe_proc equal to a number >=0. Therefore, ReduceIntMax communicates to all the
             * processors the rank of the processor which contains the point
             */
            amrex::ParallelDescriptor::ReduceIntMax(probe_proc);

            if (probe_proc != amrex::ParallelDescriptor::IOProcessorNumber() and probe_proc != -1) {
                if (amrex::ParallelDescriptor::MyProc() == probe_proc) {
                    amrex::ParallelDescriptor::Send(m_data.data(), noutputs,
                                                    amrex::ParallelDescriptor::IOProcessorNumber(),
                                                    0);
                }
                if (amrex::ParallelDescriptor::MyProc()
                    == amrex::ParallelDescriptor::IOProcessorNumber()) {
                    amrex::ParallelDescriptor::Recv(m_data.data(), noutputs, probe_proc, 0);
                }
            }
        } // send to IO Processor
    }// end loop over refinement levels
} // end void FieldProbe::ComputeDiags

void FieldProbe::WriteToFile (int step) const
{
    if (ProbeInDomain() && amrex::ParallelDescriptor::IOProcessor())
    {
        ReducedDiags::WriteToFile (step);
    }
}
