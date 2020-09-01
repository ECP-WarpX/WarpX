/* Copyright 2019-2020 Neil Zaim, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */

#include "FieldReduced.H"
#include "WarpX.H"

using namespace amrex;

// constructor
FieldReduced::FieldReduced (const std::string& rd_name, const std::vector<std::string>& field_type)
: ReducedDiags{rd_name}
{

    // RZ coordinate is not working
#if (defined WARPX_DIM_RZ)
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(false,
        "FieldReduced reduced diagnostics does not work for RZ coordinate.");
#endif

    for (int i = 0; i < field_type.size(); ++i)
    {
        std::string type = field_type[i];

        if (type.compare("FieldEnergy") == 0) {m_fieldEnergy = true;}
        else if (type.compare("MaxField") == 0) {m_maxField = true;}
        else { Abort("No matching type found for reduced field diagnostic."); }
    }

    // read number of levels
    int nLevel = 0;
    ParmParse pp("amr");
    pp.query("max_level", nLevel);
    nLevel += 1;

    // resize data array
    int data_size = 0;
    constexpr int noutputs_fieldEnergy = 3; // total energy, E-field energy and B-field energy
    constexpr int noutputs_maxField = 8; // max of Ex,Ey,Ez,|E|,Bx,By,Bz and |B|
    if (m_fieldEnergy) {data_size += noutputs_fieldEnergy*nLevel;}
    if (m_maxField)
    {
        m_offset_maxField = data_size;
        data_size += noutputs_maxField*nLevel;
    }
    m_data.resize(data_size,0.0);

    if (ParallelDescriptor::IOProcessor())
    {
        if ( m_IsNotRestart )
        {
            // open file
            std::ofstream ofs;
            ofs.open(m_path + m_rd_name + "." + m_extension,
                std::ofstream::out | std::ofstream::app);
            // write header row
            ofs << "#";
            ofs << "[1]step()";
            ofs << m_sep;
            ofs << "[2]time(s)";

            if (m_fieldEnergy)
            {
                constexpr int shift_total = 3;
                constexpr int shift_E = 4;
                constexpr int shift_B = 5;
                for (int lev = 0; lev < nLevel; ++lev)
                {
                    ofs << m_sep;
                    ofs << "[" + std::to_string(shift_total+noutputs_fieldEnergy*lev) + "]";
                    ofs << "total_lev"+std::to_string(lev)+"(J)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(shift_E+noutputs_fieldEnergy*lev) + "]";
                    ofs << "E_lev"+std::to_string(lev)+"(J)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(shift_B+noutputs_fieldEnergy*lev) + "]";
                    ofs << "B_lev"+std::to_string(lev)+"(J)";
                }
            }

            if (m_maxField)
            {
                constexpr int shift_Ex = 3;
                constexpr int shift_Ey = 4;
                constexpr int shift_Ez = 5;
                constexpr int shift_absE = 6;
                constexpr int shift_Bx = 7;
                constexpr int shift_By = 8;
                constexpr int shift_Bz = 9;
                constexpr int shift_absB = 10;
                for (int lev = 0; lev < nLevel; ++lev)
                {
                    ofs << m_sep;
                    ofs << "[" + std::to_string(m_offset_maxField+
                                                shift_Ex+noutputs_maxField*lev) + "]";
                    ofs << "max_Ex_lev"+std::to_string(lev)+" (V/m)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(m_offset_maxField+
                                                shift_Ey+noutputs_maxField*lev) + "]";
                    ofs << "max_Ey_lev"+std::to_string(lev)+" (V/m)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(m_offset_maxField+
                                                shift_Ez+noutputs_maxField*lev) + "]";
                    ofs << "max_Ez_lev"+std::to_string(lev)+" (V/m)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(m_offset_maxField+
                                                shift_absE+noutputs_maxField*lev) + "]";
                    ofs << "max_|E|_lev"+std::to_string(lev)+" (V/m)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(m_offset_maxField+
                                                shift_Bx+noutputs_maxField*lev) + "]";
                    ofs << "max_Bx_lev"+std::to_string(lev)+" (T)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(m_offset_maxField+
                                                shift_By+noutputs_maxField*lev) + "]";
                    ofs << "max_By_lev"+std::to_string(lev)+" (T)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(m_offset_maxField+
                                                shift_Bz+noutputs_maxField*lev) + "]";
                    ofs << "max_Bz_lev"+std::to_string(lev)+" (T)";
                    ofs << m_sep;
                    ofs << "[" + std::to_string(m_offset_maxField+
                                                shift_absB+noutputs_maxField*lev) + "]";
                    ofs << "max_|B|_lev"+std::to_string(lev)+" (T)";
                }
            }

            ofs << std::endl;
            // close file
            ofs.close();
        }
    }

}
// end constructor

// function that computes field energy and/or maximum field values
void FieldReduced::ComputeDiags (int step)
{

    // Judge if the diags should be done
    if ( (step+1) % m_freq != 0 ) { return; }

    // get WarpX class object
    auto & warpx = WarpX::GetInstance();

    // get number of level
    const auto nLevel = warpx.finestLevel() + 1;

    // loop over refinement levels
    for (int lev = 0; lev < nLevel; ++lev)
    {

        // get MultiFab data at lev
        const MultiFab & Ex = warpx.getEfield(lev,0);
        const MultiFab & Ey = warpx.getEfield(lev,1);
        const MultiFab & Ez = warpx.getEfield(lev,2);
        const MultiFab & Bx = warpx.getBfield(lev,0);
        const MultiFab & By = warpx.getBfield(lev,1);
        const MultiFab & Bz = warpx.getBfield(lev,2);

        if (m_fieldEnergy)
        {
            // get cell size
            Geometry const & geom = warpx.Geom(lev);
#if (AMREX_SPACEDIM == 2)
            auto dV = geom.CellSize(0) * geom.CellSize(1);
#elif (AMREX_SPACEDIM == 3)
            auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
#endif

            // compute E squared
            Real const tmpEx = Ex.norm2(0,geom.periodicity());
            Real const tmpEy = Ey.norm2(0,geom.periodicity());
            Real const tmpEz = Ez.norm2(0,geom.periodicity());
            Real const Es = tmpEx*tmpEx + tmpEy*tmpEy + tmpEz*tmpEz;

            // compute B squared
            Real const tmpBx = Bx.norm2(0,geom.periodicity());
            Real const tmpBy = By.norm2(0,geom.periodicity());
            Real const tmpBz = Bz.norm2(0,geom.periodicity());
            Real const Bs = tmpBx*tmpBx + tmpBy*tmpBy + tmpBz*tmpBz;

            constexpr int noutputs_fieldEnergy = 3; // total energy, E-field energy and B-field energy
            constexpr int index_total = 0;
            constexpr int index_E = 1;
            constexpr int index_B = 2;

            // save data
            m_data[lev*noutputs_fieldEnergy+index_E] = 0.5_rt * Es * PhysConst::ep0 * dV;
            m_data[lev*noutputs_fieldEnergy+index_B] = 0.5_rt * Bs / PhysConst::mu0 * dV;
            m_data[lev*noutputs_fieldEnergy+index_total] = m_data[lev*noutputs_fieldEnergy+index_E] +
                                                           m_data[lev*noutputs_fieldEnergy+index_B];
        }

        if (m_maxField)
        {
            constexpr int noutputs_maxField = 8; // max of Ex,Ey,Ez,|E|,Bx,By,Bz and |B|
            constexpr int index_Ex = 0;
            constexpr int index_Ey = 1;
            constexpr int index_Ez = 2;
            constexpr int index_absE = 3;
            constexpr int index_Bx = 4;
            constexpr int index_By = 5;
            constexpr int index_Bz = 6;
            constexpr int index_absB = 7;
            // get Maximums of E field components
            m_data[m_offset_maxField+lev*noutputs_maxField+index_Ex] = Ex.norm0();
            m_data[m_offset_maxField+lev*noutputs_maxField+index_Ey] = Ey.norm0();
            m_data[m_offset_maxField+lev*noutputs_maxField+index_Ez] = Ez.norm0();

            // get Maximums of B field components
            m_data[m_offset_maxField+lev*noutputs_maxField+index_Bx] = Bx.norm0();
            m_data[m_offset_maxField+lev*noutputs_maxField+index_By] = By.norm0();
            m_data[m_offset_maxField+lev*noutputs_maxField+index_Bz] = Bz.norm0();

            // Create temporary MultiFAB to be filled with |E| and |B| squared
            // Note that allocating a full MultiFAB is probably not optimal for memory usage
            const int ncomp = 1;
            const int ngrow = 0;  // no ghost cells for temporary MultiFAB
            MultiFab mftemp(amrex::convert(Ex.boxArray(), IntVect{AMREX_D_DECL(0,0,0)}),
                                                Ex.DistributionMap(), ncomp, ngrow);
            // Temporary MultiFab is cell-centered so that it can be filled for any staggering
            // (possibly with the sum of components which do not have the same staggering).

            // MFIter loop to fill temporary MultiFAB with |E| squared.
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
            for ( MFIter mfi(mftemp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                const Box& box = mfi.tilebox();

                const auto& arrEx = Ex[mfi].array();
                const auto& arrEy = Ey[mfi].array();
                const auto& arrEz = Ez[mfi].array();
                auto arrtemp = mftemp[mfi].array();

                amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i, int j, int k){
                arrtemp(i,j,k) = arrEx(i,j,k)*arrEx(i,j,k) + arrEy(i,j,k)*arrEy(i,j,k)
                                + arrEz(i,j,k)*arrEz(i,j,k);
            });
            }
            m_data[m_offset_maxField+lev*noutputs_maxField+index_absE] = std::sqrt(mftemp.max(0));

            // MFIter loop to fill temporary MultiFAB with |B| squared.
 #ifdef _OPENMP
 #pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
 #endif
            for ( MFIter mfi(mftemp, TilingIfNotGPU()); mfi.isValid(); ++mfi )
            {
                const Box& box = mfi.tilebox();

                const auto& arrBx = Bx[mfi].array();
                const auto& arrBy = By[mfi].array();
                const auto& arrBz = Bz[mfi].array();
                auto arrtemp = mftemp[mfi].array();

                amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i, int j, int k){
                arrtemp(i,j,k) = arrBx(i,j,k)*arrBx(i,j,k) + arrBy(i,j,k)*arrBy(i,j,k)
                                + arrBz(i,j,k)*arrBz(i,j,k);
            });
            }
            m_data[m_offset_maxField+lev*noutputs_maxField+index_absB] = std::sqrt(mftemp.max(0));
        }

    }
    // end loop over refinement levels

    /* if fieldEnergy is activated,
     *  m_data now contains up-to-date values for:
     *  [total field energy at level 0,
     *   electric field energy at level 0,
     *   magnetic field energy at level 0,
     *   total field energy at level 1,
     *   electric field energy at level 1,
     *   magnetic field energy at level 1,
     *   ......] */

    /* if maxField is activated,
     *  m_data now contains up-to-date values for:
     *  [.....,
     *  max(Ex),max(Ey),max(Ez),max(|E|),
     *   max(Bx),max(By),max(Bz),max(|B|)]
     * for each level */

}
// end void FieldReduced::ComputeDiags
