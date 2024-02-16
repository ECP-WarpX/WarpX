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

    inline
    bool IsDefined() const { return m_is_defined; }

    inline
    void Define(const WarpXSolverVec& a_vec)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            a_vec.IsDefined(),
            "WarpXSolverVec::Define(a_vec) called with undefined a_vec");
        Define( a_vec.getVec() );
        m_is_defined = true;
        return;
    }

    inline
    void Define( const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& a_solver_vec )
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            !IsDefined(),
            "WarpXSolverVec::Define() called on undefined WarpXSolverVec");
        const int lev = 0;
        m_field_vec.resize(1); // size is number of levels
        for (int dir=0; dir<3; dir++) {
            const amrex::MultiFab& mf_model = *a_solver_vec[lev][dir];
            m_field_vec[lev][dir] = std::make_unique<amrex::MultiFab>(mf_model.boxArray(), mf_model.DistributionMap(),
                                                                      mf_model.nComp(), amrex::IntVect::Zero);
            //amrex::Print() << "mf_model.boxArray() = " << mf_model.boxArray() << std::endl;
            //amrex::Print() << "m_field_.boxArray() = " << m_field_vec[lev][dir]->boxArray() << std::endl;
        }
        m_is_defined = true;
        return;
    }
    
    inline 
    void Copy( const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& a_solver_vec )
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::Copy() called on undefined WarpXSolverVec");
        //CH_assert(a_vec.getVectorSize() == m_size);
        // NEED A WAY TO ASSSERT A_VEC IS SAME TYPE AS THIS VEC
        const int ncomps = 1;
        amrex::MultiFab::Copy(*m_field_vec[0][0], *a_solver_vec[0][0], 0, 0, ncomps, amrex::IntVect::Zero);
        amrex::MultiFab::Copy(*m_field_vec[0][1], *a_solver_vec[0][1], 0, 0, ncomps, amrex::IntVect::Zero);
        amrex::MultiFab::Copy(*m_field_vec[0][2], *a_solver_vec[0][2], 0, 0, ncomps, amrex::IntVect::Zero);
    }

    inline 
    void Copy(const WarpXSolverVec& a_vec)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            a_vec.IsDefined(),
            "WarpXSolverVec::Copy(a_vec) called with undefined a_vec");
        if (!IsDefined()) {
            Define(a_vec);
        } else {
            //CH_assert(a_vec.getVectorSize() == m_size);
            // NEED A WAY TO ASSSERT A_VEC IS SAME TYPE AS THIS VEC
        }
        const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& field_vec = a_vec.getVec();
        Copy( field_vec );
    }

    inline 
    void operator= (const WarpXSolverVec& a_vec)
    {
      Copy(a_vec);
    }

    inline
    void operator+= (const WarpXSolverVec& a_vec)
    {
        m_field_vec[0][0]->plus(*(a_vec.getVec()[0][0]), 0, 1, 0);
        m_field_vec[0][1]->plus(*(a_vec.getVec()[0][1]), 0, 1, 0);
        m_field_vec[0][2]->plus(*(a_vec.getVec()[0][2]), 0, 1, 0);
    }

    inline
    WarpXSolverVec operator+ ( const WarpXSolverVec& a_vec ) const
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
    void zero()
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::zero() called on undefined WarpXSolverVec");
        const int lev = 0;
        m_field_vec[lev][0]->setVal(0.0);
        m_field_vec[lev][1]->setVal(0.0);
        m_field_vec[lev][2]->setVal(0.0);
    }

    inline
    void ones()
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::ones() called on undefined WarpXSolverVec");
        const int lev = 0;
        m_field_vec[lev][0]->setVal(1.0);
        m_field_vec[lev][1]->setVal(1.0);
        m_field_vec[lev][2]->setVal(1.0);
    }

    inline
    amrex::Real norm() const
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::norm() called on undefined WarpXSolverVec");
        amrex::Real nF0 = m_field_vec[0][0]->norm0(0, 0);
        amrex::Real nF1 = m_field_vec[0][1]->norm0(0, 0);
        amrex::Real nF2 = m_field_vec[0][2]->norm0(0, 0);
        amrex::Real normF = std::sqrt(nF0*nF0 + nF1*nF1 + nF2*nF2);
        //amrex::Real normF = std::max(nF0, std::max(nF1,nF2));
        return normF;
    }

    const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& getVec() const {return m_field_vec;}
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& getVec() {return m_field_vec;}

  private:

    bool  m_is_defined;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_field_vec;
    
};

#endif
