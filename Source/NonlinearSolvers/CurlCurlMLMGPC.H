/* Copyright 2024 Debojyoti Ghosh
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef CURL_CURL_MLMG_PC_H_
#define CURL_CURL_MLMG_PC_H_

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLCurlCurl.H>
#include <AMReX_MLMG.H>
#include <AMReX_GMRES_MLMG.H>

#include "FieldSolver/Fields.H"
#include "Utils/WarpXConst.H"
#include "Preconditioner.H"

/**
 * \brief Curl-curl Preconditioner
 *
 *  Preconditioner that solves the curl-curl equation for the E-field, given
 *  a RHS. Uses AMReX's curl-curl linear operator and multigrid solver.
 *
 *  This class is templated on a solution-type class T and an operator class Ops.
 *
 *  The Ops class must have the following function:
 *      + Return number of AMR levels
 *      + Return the amrex::Geometry object given an AMR level
 *      + Return hi and lo linear operator boundaries
 *      + Return the time step factor (theta) for the time integration scheme
 *
 *  The T class must have the following functions:
 *      + Return underlying vector of amrex::MultiFab arrays
 */

template <class T, class Ops>
class CurlCurlMLMGPC : public Preconditioner<T,Ops>
{
    public:

        using RT = typename T::value_type;

        /**
         * \brief Default constructor
         */
        CurlCurlMLMGPC () = default;

        /**
         * \brief Default destructor
         */
        ~CurlCurlMLMGPC () override = default;

        /**
         * \brief Define the preconditioner
         */
        virtual void Define (const T&, Ops* const) override;

        /**
         * \brief Update the preconditioner
         */
        virtual void Update (const T&) override;

        /**
         * \brief Apply (solve) the preconditioner given a RHS
         */
        virtual void Apply (T&, const T&) override;

        /**
         * \brief Print parameters
         */
        virtual void printParameters() const override;

        /**
         * \brief Check if the nonlinear solver has been defined.
         */
        [[nodiscard]] inline bool IsDefined () const override { return m_is_defined; }

    protected:

        using MFArr = amrex::Array<amrex::MultiFab,3>;

        bool m_is_defined = false;

        bool m_verbose = true;
        bool m_bottom_verbose = false;
        bool m_agglomeration = true;
        bool m_consolidation = true;
        bool m_use_gmres = false;
        bool m_use_gmres_pc = true;

        int m_max_iter = 300;
        int m_max_coarsening_level = 30;

        RT m_atol = 1.0e-16;
        RT m_rtol = 1.0e-10;

        Ops* m_ops;

        int m_num_amr_levels = 0;
        amrex::Vector<amrex::Geometry> m_geom;
        amrex::Vector<amrex::BoxArray> m_grids;
        amrex::Vector<amrex::DistributionMapping> m_dmap;
        amrex::IntVect m_gv;

        amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> m_bc_lo;
        amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> m_bc_hi;

        std::unique_ptr<amrex::LPInfo> m_info;
        std::unique_ptr<amrex::MLCurlCurl> m_curl_curl;
        std::unique_ptr<amrex::MLMGT<MFArr>> m_solver;
        std::unique_ptr<amrex::GMRESMLMGT<MFArr>> m_gmres_solver;

        /**
         * \brief Read parameters
         */
        void readParameters();

        /**
         * \brief Define staggered-grid MultiFabs for AMReX's curl-curl operator given WarpX's fields
         */
         void defineMLCCAMFFromWarpXAMF ( amrex::Array<amrex::MultiFab,3>&,
                                          const amrex::Array<std::unique_ptr<amrex::MultiFab>,3>&,
                                          const amrex::IntVect& );

        /**
         * \brief Copy to AMReX's staggered-grid MultiFabs from WarpX's staggered-grid MultiFabs
         */
         void copyMLCCAMFFromWarpXAMF ( amrex::Array<amrex::MultiFab,3>&,
                                        const amrex::Array<std::unique_ptr<amrex::MultiFab>,3>&,
                                        const amrex::IntVect& );

        /**
         * \brief Copy to WarpX's staggered-grid MultiFabs from AMReX's staggered-grid MultiFabs
         */
         void copyWarpXAMFFromMLCCAMF ( amrex::Array<std::unique_ptr<amrex::MultiFab>,3>&,
                                        const amrex::Array<amrex::MultiFab,3>&,
                                        const amrex::IntVect& );

    private:

};

template <class T, class Ops>
void CurlCurlMLMGPC<T,Ops>::printParameters() const
{
    using namespace amrex;
    Print() << PreconditionerTypes::curl_curl_mlmg << " verbose:              " << (m_verbose?"true":"false") << "\n";
    Print() << PreconditionerTypes::curl_curl_mlmg << " bottom verbose:       " << (m_bottom_verbose?"true":"false") << "\n";
    Print() << PreconditionerTypes::curl_curl_mlmg << " max iter:             " << m_max_iter << "\n";
    Print() << PreconditionerTypes::curl_curl_mlmg << " agglomeration:        " << m_agglomeration << "\n";
    Print() << PreconditionerTypes::curl_curl_mlmg << " consolidation:        " << m_consolidation << "\n";
    Print() << PreconditionerTypes::curl_curl_mlmg << " max_coarsening_level: " << m_max_coarsening_level << "\n";
    Print() << PreconditionerTypes::curl_curl_mlmg << " absolute tolerance:   " << m_atol << "\n";
    Print() << PreconditionerTypes::curl_curl_mlmg << " relative tolerance:   " << m_rtol << "\n";
    Print() << PreconditionerTypes::curl_curl_mlmg << " use GMRES:            " << (m_use_gmres?"true":"false") << "\n";
    if (m_use_gmres) {
        Print() << PreconditionerTypes::curl_curl_mlmg
                << " use PC for GMRES:     "
                << (m_use_gmres_pc?"true":"false") << "\n";
    }
}

template <class T, class Ops>
void CurlCurlMLMGPC<T,Ops>::readParameters()
{
    amrex::ParmParse pp(PreconditionerTypes::curl_curl_mlmg);
    pp.query("verbose", m_verbose);
    pp.query("bottom_verbose", m_bottom_verbose);
    pp.query("max_iter", m_max_iter);
    pp.query("agglomeration", m_agglomeration);
    pp.query("consolidation", m_consolidation);
    pp.query("max_coarsening_level", m_max_coarsening_level);
    pp.query("absolute_tolerance",  m_atol);
    pp.query("relative_tolerance",  m_rtol);
    pp.query("use_gmres",  m_use_gmres);
    pp.query("use_gmres_pc",  m_use_gmres_pc);

    return;
}

template <class T, class Ops>
void CurlCurlMLMGPC<T,Ops>::Define ( const T& a_U,
                                     Ops* const a_ops )
{
    using namespace amrex;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        !IsDefined(),
        "CurlCurlMLMGPC::Define() called on defined object" );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        (a_ops != nullptr),
        "CurlCurlMLMGPC::Define(): a_ops is nullptr" );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_U.getArrayVecType()==warpx::fields::FieldType::Efield_fp,
        "CurlCurlMLMGPC::Define() must be called with Efield_fp type");

    m_ops = a_ops;
    // read preconditioner parameters
    readParameters();

    // create info object for curl-curl op
    m_info = std::make_unique<LPInfo>();
    m_info->setAgglomeration(m_agglomeration);
    m_info->setConsolidation(m_consolidation);
    m_info->setMaxCoarseningLevel(m_max_coarsening_level);

    // Get data vectors from a_U
    auto& u_mfarrvec = a_U.getArrayVec();

    // Set number of AMR levels and create geometry, grids, and
    // distribution mapping vectors.
    m_num_amr_levels = m_ops->numAMRLevels();
    m_geom.resize(m_num_amr_levels);
    m_grids.resize(m_num_amr_levels);
    m_dmap.resize(m_num_amr_levels);
    for (int n = 0; n < m_num_amr_levels; n++) {
        m_geom[n] = m_ops->GetGeometry(n);
        m_dmap[n] = u_mfarrvec[n][0]->DistributionMap();

        BoxArray ba = u_mfarrvec[n][0]->boxArray();
        m_grids[n] = ba.enclosedCells();
    }

    // Construct the curl-curl linear operator and set its BCs
    m_curl_curl = std::make_unique<MLCurlCurl>(m_geom, m_grids, m_dmap, *m_info);
    m_curl_curl->setDomainBC(m_ops->GetLinOpBCLo(), m_ops->GetLinOpBCHi());
    m_curl_curl->setScalars(1.0, 1.0);

    // Construct the MLMG solver
    m_solver = std::make_unique<MLMGT<MFArr>>(*m_curl_curl);
    m_solver->setMaxIter(m_max_iter);
    m_solver->setVerbose(static_cast<int>(m_verbose));
    m_solver->setBottomVerbose(static_cast<int>(m_bottom_verbose));

    // If using GMRES solver, construct it
    if (m_use_gmres) {
        m_gmres_solver = std::make_unique<GMRESMLMGT<MFArr>>(*m_solver);
        m_gmres_solver->usePrecond(m_use_gmres_pc);
        m_gmres_solver->setPrecondNumIters(m_max_iter);
        m_gmres_solver->setVerbose(static_cast<int>(m_verbose));
    }

    m_is_defined = true;
    return;
}

template <class T, class Ops>
void CurlCurlMLMGPC<T,Ops>::Update (const T& a_U)
{
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        IsDefined(),
        "CurlCurlMLMGPC::Update() called on undefined object" );

    // a_U is not needed for a linear operator
    amrex::ignore_unused(a_U);

    // set the coefficients alpha and beta for curl-curl op
    RT alpha = (m_ops->theta()*this->m_dt*PhysConst::c) * (m_ops->theta()*this->m_dt*PhysConst::c);
    RT beta = RT(1.0);
    m_curl_curl->setScalars(alpha, beta);
    if (m_verbose) {
        amrex::Print() << "Updating " << PreconditionerTypes::curl_curl_mlmg
                       << ": dt = " << this->m_dt << ", "
                       << " coefficients: "
                       << "alpha = " << alpha << ", "
                       << "beta = " << beta << "\n";
    }

    return;
}

template <class T, class Ops>
void CurlCurlMLMGPC<T,Ops>::defineMLCCAMFFromWarpXAMF ( amrex::Array<amrex::MultiFab,3>& a_mfarr_mlcc,
                                                        const amrex::Array<std::unique_ptr<amrex::MultiFab>,3>& a_mfarr_warpx,
                                                        const amrex::IntVect& ng )
{
#if defined(WARPX_DIM_1D_Z)
    WARPX_ABORT_WITH_MESSAGE("CurlCurlMLMGPC::defineMLCCAMFFromWarpXAMF not yet implemented for 1D");
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    // In WarpX, missing dimension is y
    // In AMReX, missing dimension is z
    {
        a_mfarr_mlcc[0].define( a_mfarr_warpx[0]->boxArray(),
                                a_mfarr_warpx[0]->DistributionMap(),
                                a_mfarr_warpx[0]->nComp(),
                                ng );
    }
    {
        a_mfarr_mlcc[1].define( a_mfarr_warpx[2]->boxArray(),
                                a_mfarr_warpx[2]->DistributionMap(),
                                a_mfarr_warpx[2]->nComp(),
                                ng );
    }
    {
        a_mfarr_mlcc[2].define( a_mfarr_warpx[1]->boxArray(),
                                a_mfarr_warpx[1]->DistributionMap(),
                                a_mfarr_warpx[1]->nComp(),
                                ng );
    }
#elif defined(WARPX_DIM_3D)
    for (int d = 0; d < 3; d++) {
        a_mfarr_mlcc[d].define( a_mfarr_warpx[d]->boxArray(),
                                a_mfarr_warpx[d]->DistributionMap(),
                                a_mfarr_warpx[d]->nComp(),
                                ng );
    }
#endif
}

template <class T, class Ops>
void CurlCurlMLMGPC<T,Ops>::copyMLCCAMFFromWarpXAMF ( amrex::Array<amrex::MultiFab,3>& a_dst,
                                                      const amrex::Array<std::unique_ptr<amrex::MultiFab>,3>& a_src,
                                                      const amrex::IntVect& ng )
{
#if defined(WARPX_DIM_1D_Z)
    WARPX_ABORT_WITH_MESSAGE("CurlCurlMLMGPC::copyMLCCAMFFromWarpXAMF not yet implemented for 1D");
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    // In WarpX, missing dimension is y
    // In AMReX, missing dimension is z
    amrex::MultiFab::Copy( a_dst[0], *(a_src[0]), 0, 0, a_dst[0].nComp(), ng);
    amrex::MultiFab::Copy( a_dst[1], *(a_src[2]), 0, 0, a_dst[0].nComp(), ng);
    amrex::MultiFab::Copy( a_dst[2], *(a_src[1]), 0, 0, a_dst[0].nComp(), ng);
#elif defined(WARPX_DIM_3D)
    for (int d = 0; d < 3; d++) {
        amrex::MultiFab::Copy( a_dst[d], *(a_src[d]), 0, 0, ncomp, ng);
    }
#endif
}

template <class T, class Ops>
void CurlCurlMLMGPC<T,Ops>::copyWarpXAMFFromMLCCAMF ( amrex::Array<std::unique_ptr<amrex::MultiFab>,3>& a_dst,
                                                      const amrex::Array<amrex::MultiFab,3>& a_src,
                                                      const amrex::IntVect& ng )
{
#if defined(WARPX_DIM_1D_Z)
    WARPX_ABORT_WITH_MESSAGE("CurlCurlMLMGPC::copyMLCCAMFFromWarpXAMF not yet implemented for 1D");
#elif defined(WARPX_DIM_XZ) || defined(WARPX_DIM_RZ)
    // In WarpX, missing dimension is y
    // In AMReX, missing dimension is z
    amrex::MultiFab::Copy( *(a_dst[0]), a_src[0], 0, 0, a_dst[0]->nComp(), ng);
    amrex::MultiFab::Copy( *(a_dst[1]), a_src[2], 0, 0, a_dst[0]->nComp(), ng);
    amrex::MultiFab::Copy( *(a_dst[2]), a_src[1], 0, 0, a_dst[0]->nComp(), ng);
#elif defined(WARPX_DIM_3D)
    for (int d = 0; d < 3; d++) {
        amrex::MultiFab::Copy( *(a_dst[d]), a_src[d], 0, 0, ncomp, ng);
    }
#endif
}

template <class T, class Ops>
void CurlCurlMLMGPC<T,Ops>::Apply (T& a_x, const T& a_b)
{
    using namespace amrex;

    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        IsDefined(),
        "CurlCurlMLMGPC::Apply() called on undefined object" );
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_x.getArrayVecType()==warpx::fields::FieldType::Efield_fp,
        "CurlCurlMLMGPC::Apply() - a_x must be Efield_fp type");
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        a_b.getArrayVecType()==warpx::fields::FieldType::Efield_fp,
        "CurlCurlMLMGPC::Apply() - a_b must be Efield_fp type");

    // Get the data vectors
    auto& b_mfarrvec = a_b.getArrayVec();
    auto& x_mfarrvec = a_x.getArrayVec();
    WARPX_ALWAYS_ASSERT_WITH_MESSAGE(
        ((b_mfarrvec.size() == m_num_amr_levels) && (x_mfarrvec.size() == m_num_amr_levels)),
        "Error in CurlCurlMLMGPC::Apply() - mismatch in number of levels." );

    for (int n = 0; n < m_num_amr_levels; n++) {

        // Copy initial guess to local object
        Array<MultiFab,3> solution, rhs;
        defineMLCCAMFFromWarpXAMF(solution, x_mfarrvec[n], IntVect::TheUnitVector());
        defineMLCCAMFFromWarpXAMF(rhs, b_mfarrvec[n], IntVect::TheZeroVector());
        copyMLCCAMFFromWarpXAMF(solution, x_mfarrvec[n], IntVect::TheZeroVector());
        copyMLCCAMFFromWarpXAMF(rhs, b_mfarrvec[n], IntVect::TheZeroVector());

        m_curl_curl->prepareRHS({&rhs});
        if (m_use_gmres) {
            m_gmres_solver->solve(solution, rhs, m_rtol, m_atol);
        } else {
            m_solver->solve({&solution}, {&rhs}, m_rtol, m_atol);
        }

        // Copy solution in local object to a_x
        copyWarpXAMFFromMLCCAMF(x_mfarrvec[n], solution, IntVect::TheZeroVector());

    }

    return;
}

#endif
