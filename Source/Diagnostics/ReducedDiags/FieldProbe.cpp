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

    /*
     * Obtain input data from parsing inputs file.
     * For the case of a single particle:
     *     Define x, y, and z of particle
     *     Define whether or not to integrate fields
     * For the case of a line detector:
     *     Define x, y, and z of end of line point 1
     *     Define x, y, and z of end of line point 2
     *     Define resoliton to determine number of particles
     *     Define whether ot not to integrate fields
     * For the case of a plane detector:
     *     Define a vector normal to the detector plane
     *     Define a vector in the "up" direction of the plane
     *     Define the size of the plane (width of half square)
     *     Define resoliton to determine number of particles
     *     Define whether ot not to integrate fields
     */

    amrex::ParmParse pp_rd_name(rd_name);
    getWithParser(pp_rd_name, "x_probe", x_probe);
    getWithParser(pp_rd_name, "y_probe", y_probe);
    getWithParser(pp_rd_name, "z_probe", z_probe);
    pp_rd_name.query("integrate", field_probe_integrate);
    pp_rd_name.query("raw_fields", raw_fields);
    pp_rd_name.query("interp_order", interp_order);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(interp_order <= WarpX::nox ,
                                     "Field probe interp_order should be lower or equal than algo.particle_shape");
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
                ofs << m_sep;
                ofs << "[" << c++ << "]probe_S_lev" + std::to_string(lev) + " (W/m^2)"; //update all units if integrating (might be energy / m^2)

            }
            ofs << std::endl;

            // close file

            ofs.close();
        }
    }
}//end constructor

void FieldProbe::InitData()
{

    //create 1D array for X, Y, and Z of particles

    amrex::ParticleReal xpos[1]{x_probe};
    amrex::ParticleReal ypos[1]{y_probe};
    amrex::ParticleReal zpos[1]{z_probe};

    //add np partciles on lev 0 to m_probe

    m_probe.AddNParticles(0, np, xpos, ypos, zpos);
}
// function that computes field values at probe position

void FieldProbe::ComputeDiags (int step)
{
    // Judge if the diags should be done

    if (!m_intervals.contains(step+1)) { return; }

    // get a reference to WarpX instance

    auto & warpx = WarpX::GetInstance();

    // get number of level

    const auto nLevel = warpx.finestLevel() + 1;

    // vector to store the field values

    amrex::Vector<amrex::Real> fp_values(noutputs, 0);

   	// loop over refinement levels

    for (int lev = 0; lev < nLevel; ++lev)
    {
        const amrex::Geometry& gm = warpx.Geom(lev);
        const auto prob_lo = gm.ProbLo();
        const auto prob_hi = gm.ProbHi();

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

	    //loop over each particle

    	using MyParIter = amrex::ParIter<0,0,7,0>;
       	for (MyParIter pti(m_probe, lev); pti.isValid(); ++pti) 
        {
            const auto getPosition = GetParticlePosition(pti);

            /*
             * Determine if probe exists within simulation boundaries. During 2D simulations,
             * y values will be set to 0 making it unnecessary to check. Generally, the second
             * value in a position array will be the y value, but in the case of 2D, prob_lo[1]
             * and prob_hi[1] refer to z. This is a result of warpx.Geom(lev).
             */ 

#if (AMREX_SPACEDIM == 2)
            m_probe_in_domain = x_probe >= prob_lo[0] and x_probe < prob_hi[0] and
                                z_probe >= prob_lo[1] and z_probe < prob_hi[1];
#else
            m_probe_in_domain = x_probe >= prob_lo[0] and x_probe < prob_hi[0] and
                                y_probe >= prob_lo[1] and y_probe < prob_hi[1] and
                                z_probe >= prob_lo[2] and z_probe < prob_hi[2];
#endif

            if(lev == 0) m_probe_in_domain_lev_0 = m_probe_in_domain;

            if( m_probe_in_domain ) 
            {
	            /*
                 * Make the box cell centered in preparation for the interpolation (and to avoid
                 * including ghost cells in the calculation)
                 */

                const auto &arrEx = Ex[pti].array();
                const auto &arrEy = Ey[pti].array();
                const auto &arrEz = Ez[pti].array();
                const auto &arrBx = Bx[pti].array();
                const auto &arrBy = By[pti].array();
                const auto &arrBz = Bz[pti].array();

                amrex::Box box = pti.tilebox();
                box.grow(Ex.nGrowVect());

                //preparing to write data to particle

                auto& attribs = pti.GetStructOfArrays().GetRealData();
                ParticleReal* const AMREX_RESTRICT part_Ex = attribs[static_cast<int>(ParticleVal::Ex)].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_Ey = attribs[static_cast<int>(ParticleVal::Ey)].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_Ez = attribs[static_cast<int>(ParticleVal::Ez)].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_Bx = attribs[static_cast<int>(ParticleVal::Bx)].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_By = attribs[static_cast<int>(ParticleVal::By)].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_Bz = attribs[static_cast<int>(ParticleVal::Bz)].dataPtr();
                ParticleReal* const AMREX_RESTRICT part_S = attribs[static_cast<int>(ParticleVal::S)].dataPtr();

                amrex::Vector<amrex::Real> v_galilean{amrex::Vector<amrex::Real>(3, amrex::Real(0.))};
                const auto &xyzmin = WarpX::GetInstance().LowerCornerWithGalilean(box, v_galilean, lev);
                const std::array<Real, 3> &dx = WarpX::CellSize(lev);

                const amrex::GpuArray<amrex::Real, 3> dx_arr = {dx[0], dx[1], dx[2]};
                const amrex::GpuArray<amrex::Real, 3> xyzmin_arr = {xyzmin[0], xyzmin[1], xyzmin[2]};

                // Interpolating to the probe positions for each particle

                amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
                {
                    amrex::ParticleReal xp, yp, zp;
                    getPosition(ip, xp, yp, zp);

                    amrex::ParticleReal Exp = 0._rt, Eyp = 0._rt, Ezp = 0._rt;
                    amrex::ParticleReal Bxp = 0._rt, Byp = 0._rt, Bzp = 0._rt;
                    amrex::ParticleReal S = 0._rt;

                    // first gather E and B to the particle positions

                    doGatherShapeN(xp, yp, zp, Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                                   arrEx, arrEy, arrEz, arrBx, arrBy, arrBz,
                                   Extype, Eytype, Eztype, Bxtype, Bytype, Bztype,
                                   dx_arr, xyzmin_arr, amrex::lbound(box), WarpX::n_rz_azimuthal_modes,
                                   interp_order, false);

	                //Calculate S

                    amrex::Real sraw[3]{
		            Exp * Bzp - Ezp * Byp,
		            Ezp * Bxp - Exp * Bzp,
		            Exp * Byp - Eyp * Bxp};
    	            S = (1._rt / PhysConst::mu0)  * sqrt(sraw[0] * sraw[0] + sraw[1] * sraw[1] + sraw[2] * sraw[2]);

                    /*
                     * Determine whether or not to integrate field data.
                     * If not integrating, store instantaneous values.
                     */

                    if (field_probe_integrate != 1)
                    {

                        //Store Values on Particles

                        part_Ex[ip] = Exp; //remember to add lorentz transform
                        part_Ey[ip] = Eyp; //remember to add lorentz transform
                        part_Ez[ip] = Ezp; //remember to add lorentz transform
                        part_Bx[ip] = Bxp; //remember to add lorentz transform
                        part_By[ip] = Byp; //remember to add lorentz transform
                        part_Bz[ip] = Bzp; //remember to add lorentz transform
                        part_S[ip] = S; //remember to add lorentz transform
    	            }
                    else
                    {
                        part_Ex[ip] = Exp + part_Ex[ip]; //remember to add lorentz transform
                        part_Ey[ip] = Eyp + part_Ey[ip]; //remember to add lorentz transform
                        part_Ez[ip] = Ezp + part_Ez[ip]; //remember to add lorentz transform
                        part_Bx[ip] = Bxp + part_Bx[ip]; //remember to add lorentz transform
                        part_By[ip] = Byp + part_By[ip]; //remember to add lorentz transform
                        part_Bz[ip] = Bzp + part_Bz[ip]; //remember to add lorentz transform
                        part_S[ip] = S + part_S[ip]; //remember to add lorentz transform
                    }

                    // Fill output array

                    m_data[0 * noutputs + static_cast<int>(ParticleVal::Ex)] = part_Ex[ip];
                    m_data[0 * noutputs + static_cast<int>(ParticleVal::Ey)] = part_Ey[ip];
                    m_data[0 * noutputs + static_cast<int>(ParticleVal::Ez)] = part_Ez[ip];
                    m_data[0 * noutputs + static_cast<int>(ParticleVal::Bx)] = part_Bx[ip];
                    m_data[0 * noutputs + static_cast<int>(ParticleVal::By)] = part_By[ip];
                    m_data[0 * noutputs + static_cast<int>(ParticleVal::Bz)] = part_Bz[ip];
                    m_data[0 * noutputs + static_cast<int>(ParticleVal::S)] = part_S[ip];
                });// ParallelFor Close  

                probe_proc = amrex::ParallelDescriptor::MyProc();

            }

	    }//End MyParIter

        /*
         * All the processors have probe_proc = -1 except the one that contains the point, which
         * has probe_proc equal to a number >=0. Therefore, ReduceIntMax communicates to all the
         * processors the rank of the processor which contains the point
         */

        amrex::ParallelDescriptor::ReduceIntMax(probe_proc);

        if(probe_proc != amrex::ParallelDescriptor::IOProcessorNumber() and probe_proc != -1) 
        {
            if (amrex::ParallelDescriptor::MyProc() == probe_proc) 
            {
                amrex::ParallelDescriptor::Send(fp_values.dataPtr(), noutputs, 
                                amrex::ParallelDescriptor::IOProcessorNumber(),
                                0);
            }
            if (amrex::ParallelDescriptor::MyProc()
            == amrex::ParallelDescriptor::IOProcessorNumber()) 
            {
                amrex::ParallelDescriptor::Recv(fp_values.dataPtr(), noutputs, probe_proc, 0);
            }
        }
    }// end loop over refinement levels

}// end void FieldProbe::ComputeDiags

/* 
 * m_data now contains up-to-date values for:
 *  [probe(Ex),probe(Ey),probe(Ez),
 *   probe(Bx),probe(By),probe(Bz)]
 */

//Utilize ReducedDiags output generation funtion

void FieldProbe::WriteToFile (int step) const
{
    if(m_probe_in_domain_lev_0)
    {
        ReducedDiags::WriteToFile (step);
    }
}
