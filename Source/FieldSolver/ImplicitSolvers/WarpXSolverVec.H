/* Copyright 2024 Justin Angus
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WarpXSolverVec_H_
#define WarpXSolverVec_H_

#include "Utils/TextMsg.H"

#include <ablastr/utils/SignalHandling.H>
#include <ablastr/warn_manager/WarnManager.H>

#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_BLassert.H>
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

#include "FieldSolver/Fields.H"

/**
 * \brief This is a wrapper class around a Vector of array of pointers to MultiFabs that
 *  contains basic math operators and functionality needed to interact with nonlinear
 *  solvers in WarpX and linear solvers in AMReX, such as GMRES.
 *
 *  The size of the Vector is the number of amr levels. Hardcoded for 1 right now. The
 *  size of the array is the number of MultiFabs. It is hardcoded for 3 right now as it
 *  is only used for the electric field in the implicit electromagnetic time solvers. In
 *  the future, the array size can be made a template parameter so that this class can
 *  be used for other solver vectors, such as electrostatic (array size 1) or Darwin (array size 4).
 */

class WarpX;
class WarpXSolverVec
{
public:

    WarpXSolverVec() = default;

    WarpXSolverVec(const WarpXSolverVec&) = delete;

    ~WarpXSolverVec() = default;

    using value_type = amrex::Real;
    using RT = value_type;

    [[nodiscard]] inline bool IsDefined () const { return m_is_defined; }

    void Define ( WarpX*                    a_WarpX,
            const warpx::fields::FieldType  a_solver_vec_type );

    inline
    void Define (const WarpXSolverVec& a_vec)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            a_vec.IsDefined(),
            "WarpXSolverVec::Define(a_vec) called with undefined a_vec");
        Define( a_vec.m_WarpX, a_vec.getVecType() );
    }

    [[nodiscard]] RT dotProduct( const WarpXSolverVec&  a_X ) const;

    inline
    void Copy ( const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& a_solver_vec,
                const warpx::fields::FieldType                                           a_solver_vec_type )
    {

        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::Copy() called on undefined WarpXSolverVec");
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            a_solver_vec_type==m_solver_vec_type,
            "WarpXSolverVec::Copy() called with vecs of different types");
        for (int lev=0; lev<m_num_amr_levels; ++lev) {
            for (int n=0; n<3; ++n) {
                amrex::MultiFab::Copy( *m_solver_vec[lev][n], *a_solver_vec[lev][n], 0, 0, m_ncomp,
                                       amrex::IntVect::TheZeroVector() );
            }
        }
    }

    inline
    void Copy ( const WarpXSolverVec&  a_vec )
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            a_vec.IsDefined(),
            "WarpXSolverVec::Copy(a_vec) called with undefined a_vec");
        if (IsDefined()) {
            WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
                a_vec.m_solver_vec_type==m_solver_vec_type,
                "WarpXSolverVec::Copy(a_vec) called with vecs of different types");
        }
        else { Define(a_vec); }
        const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& solver_vec = a_vec.getVec();
        Copy(solver_vec, a_vec.getVecType());
    }

    // Prohibit Copy assignment operator
    WarpXSolverVec& operator= ( const WarpXSolverVec&  a_vec ) = delete;

    // Move assignment operator
    WarpXSolverVec(WarpXSolverVec&&) noexcept = default;
    WarpXSolverVec& operator= ( WarpXSolverVec&&  a_vec ) noexcept
    {
        if (this != &a_vec) {
            m_solver_vec = std::move(a_vec.m_solver_vec);
            m_solver_vec_type = a_vec.m_solver_vec_type;
            m_is_defined = true;
        }
        return *this;
    }

    inline
    void operator+= ( const WarpXSolverVec&  a_vec )
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            a_vec.m_solver_vec_type==m_solver_vec_type,
            "WarpXSolverVec operator += a_vec called with solver vecs of different types");
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int n=0; n<3; n++) {
                m_solver_vec[lev][n]->plus(*(a_vec.getVec()[lev][n]), 0, 1, 0);
            }
        }
    }

    inline
    void operator-= (const WarpXSolverVec& a_vec)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            a_vec.m_solver_vec_type==m_solver_vec_type,
            "WarpXSolverVec operator -= a_vec called with solver vecs of different types");
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int n=0; n<3; n++) {
                m_solver_vec[lev][n]->minus(*(a_vec.getVec()[lev][n]), 0, 1, 0);
            }
        }
    }

    /**
     * \brief Y = a*X + b*Y
     */
    inline
    void linComb (const RT a, const WarpXSolverVec& X, const RT b, const WarpXSolverVec& Y)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            X.m_solver_vec_type==m_solver_vec_type &&
            Y.m_solver_vec_type==m_solver_vec_type,
            "WarpXSolverVec::linComb(a,X,b,Y) called with solver vecs of different types");
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int n=0; n<3; n++) {
                amrex::MultiFab::LinComb(*m_solver_vec[lev][n], a, *X.getVec()[lev][n], 0,
                                                                b, *Y.getVec()[lev][n], 0,
                                                                0, 1, 0);
            }
        }
    }

    /**
     * \brief Increment Y by a*X (Y += a*X)
     */
    void increment (const WarpXSolverVec& X, const RT a)
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            X.m_solver_vec_type==m_solver_vec_type,
            "WarpXSolverVec::increment(X,z) called with solver vecs of different types");
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int n=0; n<3; n++) {
                amrex::MultiFab::Saxpy( *m_solver_vec[lev][n], a, *X.getVec()[lev][n],
                                        0, 0, 1, amrex::IntVect::TheZeroVector() );
            }
        }
    }

    /**
     * \brief Scale Y by a (Y *= a)
     */
    inline
    void scale (RT a_a)
    {
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int n=0; n<3; n++) {
                m_solver_vec[lev][n]->mult(a_a, 0, 1);
            }
        }
    }

    inline
    void zero () { setVal(0.0); }

    inline
    void setVal ( const RT  a_val )
    {
        WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
            IsDefined(),
            "WarpXSolverVec::ones() called on undefined WarpXSolverVec");
        for (int lev = 0; lev < m_num_amr_levels; ++lev) {
            for (int n=0; n<3; n++) {
                m_solver_vec[lev][n]->setVal(a_val);
            }
        }
    }

    [[nodiscard]] inline RT norm2 () const
    {
        auto const norm = dotProduct(*this);
        return std::sqrt(norm);
    }

    [[nodiscard]] const amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& getVec() const {return m_solver_vec;}
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > >& getVec() {return m_solver_vec;}

    // solver vec types are for now limited to warpx::fields::FieldType
    warpx::fields::FieldType getVecType () const { return m_solver_vec_type; }

    // clearDotMask() must be called by the highest class that owns WarpXSolverVec()
    // after it is done being used ( typically in the destructor ) to avoid the
    // following error message after the simulation finishes:
    // malloc_consolidate(): unaligned fastbin chunk detected
    static void clearDotMask() { m_dotMask.clear(); }

private:

    bool m_is_defined = false;
    amrex::Vector<std::array< std::unique_ptr<amrex::MultiFab>, 3 > > m_solver_vec;
    warpx::fields::FieldType m_solver_vec_type;

    static constexpr int m_ncomp = 1;
    static constexpr int m_num_amr_levels = 1;

    inline static bool m_warpx_ptr_defined = false;
    inline static WarpX* m_WarpX = nullptr;
    void SetWarpXPointer( WarpX*  a_WarpX );

    inline static bool m_dot_mask_defined = false;
    inline static amrex::Vector<std::array<std::unique_ptr<amrex::iMultiFab>,3>> m_dotMask;
    void SetDotMask();

};

#endif
