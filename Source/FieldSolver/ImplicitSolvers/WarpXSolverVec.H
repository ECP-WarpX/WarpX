/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef _WarpXSolverVec_H_
#define _WarpXSolverVec_H_

#include "Utils/TextMsg.H"

#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
#include <AMReX_Geometry.H>
#include <AMReX_IntVect.H>
#include <AMReX_LayoutData.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <memory>
#include <ostream>
#include <vector>

/*
 *   This is a wrapper class around a Vector of array of pointers to MultiFabs that contains
 *   basic math operators and functionality needed to interact with linear solvers in
 *   AMReX, such as GMRES.
 *
 *   The size of the Vector is the number of amr levels. Hardcoded for 1 right now. The size of
 *   the array is the number of MultiFabs. It is hardcoded for 3 right now as it is only
 *   used for the electric field in the electromagnetic implicit time solver. In the future,
 *   the array size can be made a template parameter so that this class can be used for other
 *   solver vectors, such as electrostatic (array size 1) or Darwin (array size 4).
 */

class WarpXSolverVec
{
public:

    WarpXSolverVec()
    {
        m_is_defined = false;
    }

    WarpXSolverVec(const WarpXSolverVec& a_vec)
    {
        m_is_defined = false;
        Define(a_vec);
    }

    ~WarpXSolverVec() = default;

    using value_type = amrex::Real;
    using RT = value_type;

    inline
    bool IsDefined () const { return m_is_defined; }

    inline
    void Define (const WarpXSolverVec& a_vec)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            a_vec.IsDefined(),
            "WarpXSolverVec::Define(a_vec) called with undefined a_vec");
        Define( a_vec.getVec() );
        return;
    }

    inline
    void Define ( const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& a_solver_vec )
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !IsDefined(),
            "WarpXSolverVec::Define() called on undefined WarpXSolverVec");
        m_field_vec.resize(m_num_amr_levels);
        const int lev = 0;
        for (int n=0; n<3; n++) {
            const amrex::MultiFab& mf_model = *a_solver_vec[lev][n];
            m_field_vec[lev][n] = std::make_unique<amrex::MultiFab>(mf_model.boxArray(), mf_model.DistributionMap(),
                                                                      mf_model.nComp(), amrex::IntVect::Zero);
        }
        m_is_defined = true;
        return;
    }

    void SetDotMask( const amrex::Vector<amrex::Geometry>&  a_Geom );
    RT dotProduct( const WarpXSolverVec&  a_X ) const;

    //
    // JRA: make functions below are derived over from amrex/Src/LinearSolvers/MLMG/AMReX_MLLinOp.H
    //
    inline
    void make ( amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& mf, amrex::IntVect const& ng) const
    {
        mf.clear();
        mf.resize(m_num_amr_levels);
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int idim = 0; idim < 3; ++idim) {
                mf[lev][idim] = make(lev, idim, ng);
            }
        }
    }

    [[nodiscard]] std::unique_ptr<amrex::MultiFab> make (int amrlev, int idim, amrex::IntVect const& ng) const
    {
        AMREX_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::make () called on undefined WarpXSolverVec");
        const amrex::BoxArray& grids = m_field_vec[amrlev][idim]->boxArray();
        const amrex::DistributionMapping& dmap = m_field_vec[amrlev][idim]->DistributionMap();
        const amrex::FabFactory<amrex::FArrayBox>& factory = m_field_vec[amrlev][idim]->Factory();
        return std::make_unique<amrex::MultiFab>( grids, dmap, m_ncomp, ng, amrex::MFInfo(), factory );
    }

    inline
    void Copy ( const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& a_solver_vec )
    {
        AMREX_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::Copy() called on undefined WarpXSolverVec");
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int n = 0; n < 3; ++n) {
                amrex::MultiFab::Copy(*m_field_vec[lev][n], *a_solver_vec[lev][n], 0, 0, m_ncomp, amrex::IntVect::Zero);
            }
        }
    }

    inline
    void Copy ( const WarpXSolverVec&  a_vec )
    {
        AMREX_ASSERT_WITH_MESSAGE(
            a_vec.IsDefined(),
            "WarpXSolverVec::Copy(a_vec) called with undefined a_vec");
        if (!IsDefined()) { Define(a_vec); }
        const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& field_vec = a_vec.getVec();
        Copy(field_vec);
    }

    inline
    void operator= ( const WarpXSolverVec&  a_vec )
    {
      Copy(a_vec);
    }

    inline
    void operator+= ( const WarpXSolverVec&  a_vec )
    {
        m_field_vec[0][0]->plus(*(a_vec.getVec()[0][0]), 0, 1, 0);
        m_field_vec[0][1]->plus(*(a_vec.getVec()[0][1]), 0, 1, 0);
        m_field_vec[0][2]->plus(*(a_vec.getVec()[0][2]), 0, 1, 0);
    }

    inline
    WarpXSolverVec operator+ ( const WarpXSolverVec&  a_vec ) const
    {
        WarpXSolverVec retval(*this);
        retval = *this;
        retval += a_vec;
        return retval;
    }

    inline
    void operator-= (const WarpXSolverVec& a_vec)
    {
        m_field_vec[0][0]->minus(*(a_vec.getVec()[0][0]), 0, 1, 0);
        m_field_vec[0][1]->minus(*(a_vec.getVec()[0][1]), 0, 1, 0);
        m_field_vec[0][2]->minus(*(a_vec.getVec()[0][2]), 0, 1, 0);
    }

    inline
    WarpXSolverVec operator- ( const WarpXSolverVec& a_vec ) const
    {
        WarpXSolverVec retval(*this);
        retval = *this;
        retval -= a_vec;
        return retval;
    }

    inline
    void linComb (const RT a, const WarpXSolverVec& X, const RT b, const WarpXSolverVec& Y)
    {
        for (int n=0; n<3; n++) {
            amrex::MultiFab::LinComb(*m_field_vec[0][n], a, *X.getVec()[0][n], 0,
                                                         b, *Y.getVec()[0][n], 0,
                                                         0, 1, 0);
        }
    }

    inline
    void increment (const WarpXSolverVec& a_X, const RT a_a)
    {
        for (int n=0; n<3; n++) {
            amrex::MultiFab::Saxpy( *m_field_vec[0][n], a_a, *a_X.getVec()[0][n],
                                    0, 0, 1, amrex::IntVect::Zero );
        }
    }

    inline
    void scale (RT a_a)
    {
        m_field_vec[0][0]->mult(a_a, 0, 1);
        m_field_vec[0][1]->mult(a_a, 0, 1);
        m_field_vec[0][2]->mult(a_a, 0, 1);
    }

    inline
    void operator*= (RT a_a)
    {
      scale(a_a);
    }

    inline
    void zero () { setVal(0.0); }

    inline
    void ones () { setVal(RT(1.0)); }

    inline
    void setVal ( const RT  a_val )
    {
        AMREX_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::ones() called on undefined WarpXSolverVec");
        const int lev = 0;
        m_field_vec[lev][0]->setVal(a_val);
        m_field_vec[lev][1]->setVal(a_val);
        m_field_vec[lev][2]->setVal(a_val);
    }

    inline
    RT norm2 () const
    {
        RT norm = dotProduct(*this);
        return std::sqrt(norm);
    }

    const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& getVec() const {return m_field_vec;}
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& getVec() {return m_field_vec;}

    // clearDotMask() must be called by the highest class that owns WarpXSovlerVec()
    // after it is done being used ( typically in the destructor ) to avoid the
    // following error message after the simulation finishes:
    // malloc_consolidate(): unaligned fastbin chunk detected
    static void clearDotMask() { m_dotMask.clear(); }

private:

    bool  m_is_defined;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_field_vec;
    static constexpr int m_ncomp = 1;
    static constexpr int m_num_amr_levels = 1;

    static bool m_dot_mask_defined;
    static amrex::Vector<std::array<std::unique_ptr<amrex::iMultiFab>,3>> m_dotMask;

};

template <typename T>
WarpXSolverVec operator* (T a_a, const WarpXSolverVec& a_vec)
{
  WarpXSolverVec retval(a_vec);
  retval = a_vec;
  retval.scale(a_a);
  return retval;
}

#endif
